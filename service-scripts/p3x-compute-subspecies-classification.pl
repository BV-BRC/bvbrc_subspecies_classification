#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long::Descriptive;
use GenomeTypeObject;
use File::Temp;
use JSON::XS;
use IPC::Run qw(run);
use Time::HiRes qw(gettimeofday);

my ($opt, $usage) = describe_options("%c %o [< in] [> out]",
    [ "in|i=s", "Input GTO" ],
    [ "out|o=s", "Output GTO" ],
    [ "dry_run|n", "Dry run", { default => 0 } ],
    [ "help|h", "Print help" ],
);

print($usage->text), exit if $opt->help;
print($usage->text), exit 1 if @ARGV != 0;

chomp(my $hostname = `hostname -f`);

my @VIRUS_TAXON_RULES = (
    # [ taxon_id, virus_type ]
    # Matched by walking ncbi_lineage and checking tax_id
    [ 12637, 'DENGUE' ],  # Dengue virus
    [ 11234, 'MEASLES' ], # Measles morbillivirus
    [ 10244, 'MPOX' ],    # Monkeypox virus
);

# Influenza A handled separately
# Influenza A species tax_id (Alphainfluenzavirus influenzae)
# https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=2955291
my $INFLUENZA_A_TAXID = 2955291;
my @INFLUENZA_SUBTYPE_RULES = (
    [ qr/^H5/i, sub {'INFLUENZAH5'} ],
    [ qr/^H3N2/i, sub {'INFLUENZAH3N2'} ],
    [ qr/^H3/i, sub {'SWINEH3'} ],
    [ qr/^H1/i, sub {
        my $c = shift;
        $c eq 'usa' ? 'SWINEH1US' : 'SWINEH1'
    } ],
);

sub get_gto {
    my ($opt) = @_;
    my $gto =
        $opt->in ? GenomeTypeObject->new({ file => $opt->in })
            : GenomeTypeObject->new({ file => \*STDIN });

    $gto or die "Error reading GenomeTypeObject\n";
    return $gto;
}

sub lineage_has_name {
    my ($gto, $needle) = @_;
    for my $entry (@{$gto->{ncbi_lineage} // []}) {
        my ($name) = @$entry;
        return 1 if defined($name) && $name =~ /\Q$needle\E/i;
    }
    return 0;
}

sub parse_comment_metadata {
    my ($gto) = @_;
    my %meta;

    for my $contig (@{$gto->{contigs} // []}) {

        my $locus = $contig->{genbank_locus} // {};

        # capture GenBank metadata
        $meta{organism} //= $locus->{organism};
        $meta{source} //= $locus->{source};
        $meta{definition} //= $locus->{definition};

        # existing comment parsing
        my $comments = $locus->{comment} // [];
        for my $line (@$comments) {
            if ($line =~ /^\s*([^:]+?)\s*::\s*(.*?)\s*$/) {
                my ($k, $v) = ($1, $2);
                $k =~ s/\s+/_/g;
                $meta{lc($k)} = $v if defined($v) && $v ne '';
            }
        }
    }

    return \%meta;
}

sub subtype_from_gto {
    my ($gto, $meta) = @_;
    return $meta->{subtype} if $meta->{subtype};

    my $name = $gto->{scientific_name} // '';
    return uc($1) if $name =~ /\((H\d+N\d+)\)/i;
    return undef;
}

sub country_from_gto {
    my ($meta) = @_;
    return $meta->{country};
}

sub has_ha_segment {
    my ($gto, $meta) = @_;

    for my $contig (@{$gto->{contigs} // []}) {
        my $rep = $contig->{replicon_type} // '';
        my $def = $contig->{genbank_locus}->{definition} // '';

        return 1 if $rep =~ /\bHA\b/i;
        return 1 if $def =~ /\bsegment\s+4\b/i;
        return 1 if $def =~ /\bhemagglutinin\b/i;
    }

    return 1 if (($meta->{segment_name} // '') =~ /^HA$/i);
    return 0;
}

sub virus_text_for_detection {
    my ($gto, $meta) = @_;

    my @texts;

    push @texts, $gto->{scientific_name}
        if defined $gto->{scientific_name};

    push @texts, $meta->{organism} if defined $meta->{organism};
    push @texts, $meta->{source} if defined $meta->{source};
    push @texts, $meta->{definition} if defined $meta->{definition};

    # remove duplicates
    my %seen;
    @texts = grep {defined($_) && !$seen{$_}++} @texts;

    return @texts;
}

sub determine_virus_type {
    my ($gto) = @_;

    my $genome_id = $gto->{id} // 'unknown';
    my $meta = parse_comment_metadata($gto);
    my $country = lc(country_from_gto($meta) // '');
    my $lineage = $gto->{ncbi_lineage} // [];

    my $is_influenza_a = 0;
    my $subtype_name = '';

    for my $entry (@$lineage) {
        my ($name, $tax_id, $rank) = @$entry;

        # Detect Influenza A by genus tax_id
        if (defined $tax_id && $tax_id == $INFLUENZA_A_TAXID) {
            $is_influenza_a = 1;
            print STDERR "[$genome_id] Detected Influenza A via lineage\n";
        }

        # Capture subtype name from serotype rank, e.g. "H1N1 subtype"
        if (defined $rank && $rank eq 'serotype' && defined $name) {
            $subtype_name = $name;
            print STDERR "[$genome_id] Found serotype in lineage: '$subtype_name'\n";
        }
    }

    if ($is_influenza_a) {
        if (!has_ha_segment($gto, $meta)) {
            print STDERR "[$genome_id] Influenza A detected but no HA segment found; skipping classification\n";
            return undef;
        }
        print STDERR "[$genome_id] HA segment confirmed\n";

        if (!$subtype_name) {
            print STDERR "[$genome_id] Influenza A detected but no serotype rank found in lineage; skipping classification\n";
            return undef;
        }

        # Extract the H-type from serotype name, fall back to comment metadata
        my ($subtype) = ($subtype_name =~ /^(H\d+N?\d*)/i);
        if (!defined $subtype) {
            print STDERR "[$genome_id] Could not parse subtype from serotype '$subtype_name'; skipping classification\n";
            return undef;
        }
        print STDERR "[$genome_id] Parsed subtype: '$subtype' (country: '$country')\n";

        for my $rule (@INFLUENZA_SUBTYPE_RULES) {
            my ($regex, $cb) = @$rule;
            if ($subtype =~ $regex) {
                my $virus_type = $cb->($country);
                print STDERR "[$genome_id] Matched influenza subtype rule for '$subtype' -> virus_type='$virus_type'\n";
                return $virus_type;
            }
        }

        print STDERR "[$genome_id] Influenza A subtype '$subtype' did not match any known subtype rule; skipping classification\n";
        return undef; # Influenza A but unrecognized subtype
    }

    # Non-influenza: match by tax_id against taxon rules
    my %taxon_lookup = map {$_->[0] => $_->[1]} @VIRUS_TAXON_RULES;
    for my $entry (@$lineage) {
        my ($name, $tax_id) = @$entry;
        if (defined $tax_id && exists $taxon_lookup{$tax_id}) {
            my $virus_type = $taxon_lookup{$tax_id};
            print STDERR "[$genome_id] Matched non-influenza taxon rule: '$name' (tax_id=$tax_id) -> virus_type='$virus_type'\n";
            return $virus_type;
        }
    }

    print STDERR "[$genome_id] No matching virus type found in lineage; skipping classification\n";
    return undef;
}

sub virus_key_from_type {
    my ($virus_type) = @_;
    my %map = (
        SWINEH1       => "h1",
        SWINEH1US     => "h1us",
        SWINEH3       => "h3",
        INFLUENZAH3N2 => "h3",
        INFLUENZAH5   => "h5",
        DENGUE        => "dengue",
        MPOX          => "monkeypox",
        MEASLES       => "measles",
    );
    return $map{$virus_type};
}

sub clade_field_for_type {
    my ($virus_type) = @_;
    my %clades = (
        "dengue"    => "subtype",
        "h1"        => "h1_clade_global",
        "h1us"      => "h1_clade_us",
        "h3"        => "h3_clade",
        "h5"        => "h5_clade",
        "monkeypox" => "clade",
        "measles"   => "subclade",
    );
    my $key = virus_key_from_type($virus_type);
    return $clades{$key};
}

sub contig_is_ha {
    my ($contig) = @_;
    my $rep = $contig->{replicon_type} // '';
    my $def = $contig->{genbank_locus}->{definition} // '';
    return 1 if $rep =~ /\bHA\b/i;
    return 1 if $def =~ /\bsegment\s+4\b/i;
    return 1 if $def =~ /\bhemagglutinin\b/i;
    return 0;
}

sub primary_result_for_genome {
    my ($gto, $virus_type, $results) = @_;

    my @nonempty = grep {
        defined($_->{classification}) && $_->{classification} ne ''
    } @$results;

    return undef unless @nonempty;

    # Influenza: use HA-segment result
    if ($virus_type =~ /^(SWINEH1|SWINEH1US|SWINEH3|INFLUENZAH5)$/) {
        my %ha_contigs = map {
            my $c = $_;
            contig_is_ha($c) ? (($c->{id} // '') => 1) : ()
        } @{$gto->{contigs} // []};

        for my $row (@nonempty) {
            my $query = $row->{query} // '';
            my (undef, $contig_id) = split(/\|/, $query, 2);
            return $row if $contig_id && $ha_contigs{$contig_id};
        }
    }

    # Non-influenza or fallback: if all equal, use that value
    my %uniq = map {(($_->{classification} // '') => 1)} @nonempty;
    my @vals = grep {$_ ne ''} keys %uniq;
    return $nonempty[0] if @vals == 1;

    # Final fallback: first non-empty result
    return $nonempty[0];
}

sub fasta_from_gto {
    my ($gto) = @_;
    my $genome_id = $gto->{id} // "unknown";

    my $fasta = '';
    for my $contig (@{$gto->{contigs} // []}) {
        my $cid = $contig->{id} // "contig";
        my $dna = $contig->{dna} // '';
        next unless $dna ne '';
        $fasta .= ">$genome_id|$cid\n$dna\n";
    }

    die "No contig DNA found in GTO\n" unless $fasta ne '';
    return $fasta;
}

sub run_classifier {
    my ($opt, $gto, $virus_type, $fasta_string) = @_;

    my $tmpdir = File::Temp->newdir(CLEANUP => 0);
    my $job_file = "$tmpdir/job.json";

    my $job = {
        virus_type       => $virus_type,
        input_source     => "fasta_data",
        input_fasta_data => $fasta_string,
        output_file      => ($gto->{id} // "unknown"),
    };

    open(my $jfh, ">", $job_file) or die "Cannot write $job_file: $!";
    print $jfh encode_json($job);
    close($jfh);

    my @cmd = (
        "run_subspecies_classification",
        "-j", $job_file,
        "-o", "$tmpdir",
        "--result-only",
    );

    print STDERR "Invoke classifier: @cmd\n";

    if ($opt->dry_run) {
        return {
            results => [
                {
                    query          => $gto->{id} // "unknown",
                    classification => "DRY_RUN_CLASSIFICATION",
                }
            ],
            command => \@cmd,
        };
    }

    my ($stdout, $stderr);
    my $ok = run(\@cmd, ">", \$stdout, "2>", \$stderr);
    if (!$ok) {
        die "Classifier failed: @cmd\nstderr:\n$stderr\n";
    }

    my $decoded;
    eval {
        $decoded = decode_json($stdout);
    };
    die "Could not parse classifier JSON output.\nstdout:\n$stdout\nstderr:\n$stderr\n" if $@;

    $decoded->{command} = \@cmd;
    $decoded->{results} //= [];

    return $decoded;
}

sub attach_classification {
    my ($gto, $virus_type, $result, $hostname) = @_;

    my $results = $result->{results} // [];
    return unless @$results;

    my $primary = primary_result_for_genome($gto, $virus_type, $results);
    return unless $primary;

    my $classification = $primary->{classification} // "";

    # Ignore empty or unassigned results
    if (!$classification || $classification =~ /^unassigned$/i) {
        print STDERR "Classification result is empty or unassigned; skipping update.\n";
        return;
    }

    my $clade_field = clade_field_for_type($virus_type)
        or die "No clade field for virus type $virus_type\n";

    my $event = {
        tool_name      => "p3x-compute-subspecies-classification",
        parameters     => [ @{$result->{command} // []} ],
        execution_time => scalar gettimeofday,
        hostname       => $hostname,
    };

    print STDERR "Genome: ", ($gto->{id} // "unknown"), "\n";
    print STDERR "Virus type: $virus_type\n";
    print STDERR "Classification: $classification\n";
    print STDERR "Writing to GTO field: $clade_field\n";

    my $event_id = $gto->add_analysis_event($event);

    # Write genome-level field
    $gto->{$clade_field} = $classification;
}

my $gto = get_gto($opt);

my $virus_type = determine_virus_type($gto);
if (!$virus_type) {
    print STDERR "Genome not eligible for subspecies classification; skipping.\n";
}
else {
    my $fasta = fasta_from_gto($gto);
    my $result = run_classifier($opt, $gto, $virus_type, $fasta);
    attach_classification($gto, $virus_type, $result, $hostname);
}

if ($opt->out) {
    $gto->destroy_to_file($opt->out);
}
else {
    $gto->destroy_to_file(\*STDOUT);
}