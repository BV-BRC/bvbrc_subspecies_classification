#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long::Descriptive;
use GenomeTypeObject;
use JSON::XS;
use File::Spec;
use IPC::Run qw(run);
use Time::HiRes qw(gettimeofday);
use Cwd;

my ($opt, $usage) = describe_options("%c %o",
    [ "in|i=s", "Input GTO" ],
    [ "out|o=s", "Output GTO" ],
    [ "work_dir|w=s", "Working directory" ],
    [ "help|h", "Print help" ]
);

print($usage->text), exit if $opt->help;

my $gto = GenomeTypeObject->new({ file => $opt->in })
    or die "Could not read GTO\n";

# Influenza A species tax_id (Alphainfluenzavirus influenzae)
# https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=2955291
my $INFLUENZA_A_TAXID = 2955291;

# H5N1 serotype tax_id
# https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=102793
my $H5N1_TAXID = 102793;

# -------------------------------------------------
# determine if genome is H5N1
# -------------------------------------------------

sub is_h5n1 {
    my ($gto) = @_;

    my $genome_id = $gto->{id} // 'unknown';
    my $lineage = $gto->{ncbi_lineage} // [];

    # Use lineage to confirm Influenza A
    my $is_influenza_a = 0;
    for my $entry (@$lineage) {
        my ($name, $tax_id, $rank) = @$entry;
        if (defined $tax_id && $tax_id == $INFLUENZA_A_TAXID) {
            $is_influenza_a = 1;
            print STDERR "[$genome_id] Detected Influenza A via lineage (tax_id=$INFLUENZA_A_TAXID)\n";
            last;
        }
    }

    unless ($is_influenza_a) {
        print STDERR "[$genome_id] Not Influenza A; skipping\n";
        return 0;
    }

    # Build source list: genome name first, then per-contig organism and definition
    my @sources;

    push @sources, [ 'genome name', $gto->{scientific_name} ]
        if defined $gto->{scientific_name};

    for my $contig (@{$gto->{contigs} // []}) {
        my $locus = $contig->{genbank_locus} // {};
        push @sources, [ 'organism name', $locus->{organism} ]
            if defined $locus->{organism};
        push @sources, [ 'genbank definition', $locus->{definition} ]
            if defined $locus->{definition};
    }

    for my $source (@sources) {
        my ($label, $text) = @$source;
        if (my ($subtype) = $text =~ /\((H\d+N\d+)\)/i) {
            print STDERR "[$genome_id] Subtype '$subtype' found in $label: '$text'\n";
            if (uc($subtype) eq 'H5N1') {
                print STDERR "[$genome_id] Confirmed H5N1 from $label\n";
                return 1;
            }
            else {
                print STDERR "[$genome_id] Subtype is '$subtype', not H5N1; skipping\n";
                return 0;
            }
        }
    }

    print STDERR "[$genome_id] Could not find subtype in any source; skipping\n";
    return 0;
}

# -------------------------------------------------
# map contigs to influenza segments
# -------------------------------------------------

sub segment_from_contig {
    my ($contig) = @_;

    my $rep = lc($contig->{replicon_type} // "");
    my $def = lc($contig->{genbank_locus}->{definition} // "");

    return "PB2" if $rep =~ /pb2/ || $def =~ /segment\s*1/;
    return "PB1" if $rep =~ /pb1/ || $def =~ /segment\s*2/;
    return "PA" if $rep =~ /pa/ || $def =~ /segment\s*3/;
    return "HA" if $rep =~ /ha/ || $def =~ /segment\s*4|hemagglutinin/;
    return "NP" if $rep =~ /np/ || $def =~ /segment\s*5/;
    return "NA" if $rep =~ /na/ || $def =~ /segment\s*6|neuraminidase/;
    return "MP" if $rep =~ /mp|m1|m2/ || $def =~ /segment\s*7/;
    return "NS" if $rep =~ /ns/ || $def =~ /segment\s*8/;

    return undef;
}

# canonical influenza segment order
my @SEGMENT_ORDER = qw(PB2 PB1 PA HA NP NA MP NS);

# -------------------------------------------------
# build FASTA file
# -------------------------------------------------

sub build_fasta {
    my ($gto, $work_dir) = @_;

    my $genome_id = $gto->{id} // "unknown";
    my %segments;

    for my $contig (@{$gto->{contigs} // []}) {

        my $seg = segment_from_contig($contig);
        next unless $seg;

        $segments{$seg} = $contig;
    }

    my @found = keys %segments;

    if (@found != 8) {
        die "Genome does not contain all 8 influenza segments\n";
    }

    my $strain = $gto->{strain} // $gto->{scientific_name} // "unknown";
    $strain =~ s/[^A-Za-z0-9]+/-/g;
    $strain =~ s/^-|-$//g;

    my $fasta = File::Spec->catfile($work_dir, "$strain.fasta");

    open(my $fh, ">", $fasta) or die "Cannot write $fasta\n";

    for my $seg (@SEGMENT_ORDER) {

        my $contig = $segments{$seg};
        my $dna = $contig->{dna} // next;
        my $contig_id = $contig->{id} // $genome_id;

        print $fh ">$contig_id $strain $seg\n";
        print $fh "$dna\n";
    }

    close($fh);

    return $fasta;
}

# -------------------------------------------------
# run genoflu
# -------------------------------------------------

sub run_genoflu {
    my ($fasta, $work_dir) = @_;

    my $original_dir = Cwd::cwd();
    chdir($work_dir) or die "Cannot chdir to $work_dir: $!";

    my @cmd = (
        "GenoFLU",
        "-f", $fasta
    );

    print STDERR "Running: @cmd\n";

    my ($stdout, $stderr);
    run(\@cmd, ">", \$stdout, "2>", \$stderr)
        or die "GenoFLU failed:\n$stderr\n";

    chdir($original_dir) or die "Cannot chdir back to $original_dir: $!";

    opendir(my $dh, $work_dir) or die $!;
    my ($stats) = grep {/_stats\.tsv$/} readdir($dh);
    closedir($dh);

    die "No stats file produced\n" unless $stats;

    return File::Spec->catfile($work_dir, $stats);
}

# -------------------------------------------------
# read genotype from TSV
# -------------------------------------------------

sub read_genotype {
    my ($file) = @_;

    open(my $fh, "<", $file) or die "Cannot open stats file $file: $!";

    my $header = <$fh>;
    chomp $header;
    my @cols = split(/\t/, $header);

    my %idx;
    for my $i (0 .. $#cols) {
        $idx{$cols[$i]} = $i;
    }

    my $line = <$fh>;
    chomp $line;

    my @vals = split(/\t/, $line);

    close($fh);

    return $vals[$idx{"Genotype"}];
}

# -------------------------------------------------
# attach genotype to GTO
# -------------------------------------------------

sub attach_genotype {
    my ($gto, $genotype) = @_;

    my $genome_id = $gto->{id} // 'unknown';

    if (!$genotype || $genotype =~ /^Not assigned/) {
        print STDERR "[$genome_id] Skipping GTO update — no valid genotype to attach\n";
        return;
    }

    my $event = {
        tool_name      => "p3x-compute-subclade",
        execution_time => scalar gettimeofday,
    };

    my $event_id = $gto->add_analysis_event($event);

    $gto->{subclade} = $genotype;
}

if (!is_h5n1($gto)) {
    print STDERR "Genome is not H5N1 — skipping\n";
}
else {
    my $work_dir = File::Temp->newdir(CLEANUP => 0);

    my $fasta = build_fasta($gto, $work_dir);

    my $stats = run_genoflu($fasta, $work_dir);

    my $genotype = read_genotype($stats);

    attach_genotype($gto, $genotype);
}

if ($opt->out) {
    $gto->destroy_to_file($opt->out);
}
else {
    $gto->destroy_to_file(\*STDOUT);
}