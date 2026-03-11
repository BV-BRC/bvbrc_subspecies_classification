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
    ["in|i=s", "Input GTO"],
    ["out|o=s", "Output GTO"],
    ["dry_run|n", "Dry run", { default => 0 }],
    ["help|h", "Print help"],
);

print($usage->text), exit if $opt->help;
print($usage->text), exit 1 if @ARGV != 0;

chomp(my $hostname = `hostname -f`);

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
        my $comments = $contig->{genbank_locus}->{comment} // [];
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

sub determine_virus_type {
    my ($gto) = @_;

    my $meta = parse_comment_metadata($gto);
    my $name = $gto->{scientific_name} // '';
    my $subtype = subtype_from_gto($gto, $meta);
    my $country = lc(country_from_gto($meta) // '');

    if (lineage_has_name($gto, "Influenza A virus") || $name =~ /influenza a virus/i) {
        return undef unless has_ha_segment($gto, $meta);

        if (defined($subtype) && $subtype =~ /^H1/i) {
            return $country eq 'usa' ? 'SWINEH1US' : 'SWINEH1';
        }
        elsif (defined($subtype) && $subtype =~ /^H3/i) {
            return 'SWINEH3';
        }
        elsif (defined($subtype) && $subtype =~ /^H5/i) {
            return 'INFLUENZAH5';
        }
    }

    if ($name =~ /dengue/i || lineage_has_name($gto, "Dengue")) {
        return 'DENGUE';
    }

    if ($name =~ /(monkeypox|mpox)/i || lineage_has_name($gto, "Monkeypox")) {
        return 'MPOX';
    }

    return undef;
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
                    query => $gto->{id} // "unknown",
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

    my $event = {
        tool_name      => "p3x-compute-subspecies-classification",
        parameters     => [@{$result->{command} // []}],
        execution_time => scalar gettimeofday,
        hostname       => $hostname,
    };

    my $event_id = $gto->add_analysis_event($event);

    $gto->{classifications} //= [];

    for my $row (@$results) {
        my $classification = $row->{classification} // "";
        next unless $classification ne "";

        push @{$gto->{classifications}}, {
            name        => "SubspeciesClassification",
            version     => "result-only-v1",
            description => "Viral subtype/clade classification",
            comment     => "virus_type=$virus_type; query=" .
                ($row->{query} // "") .
                "; classification=$classification",
            event_id    => $event_id,
        };
    }

    $gto->{computed_subspecies_classification} = {
        virus_type => $virus_type,
        results    => $results,
        event_id   => $event_id,
    };
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
} else {
    $gto->destroy_to_file(\*STDOUT);
}