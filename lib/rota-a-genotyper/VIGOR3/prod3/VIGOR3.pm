#######################################################################################
#
# Copyright (c) 2009 - 2015 J. Craig Venter Institute.
#   This file is part of JCVI VIGOR
# 
#   JCVI VIGOR is free software: you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation, either version 3 of the License, or
#   (at your option) any later version.
#   
#   JCVI VIGOR is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#   
#   You should have received a copy of the GNU General Public License
#   along with JCVI VIGOR.  If not, see <http://www.gnu.org/licenses/>.
# 
# Contributors:
#     Shiliang Wang - Initial idea and implementation.
#     Jeff Hoover - Redesigning, refactoring, and expanding the scope.
#     Susmita Shrivastava and Neha Gupta - Creation and curation of sequence databases.
#     Paolo Amedeo and Danny Katzel - Maintenance and further improvements.
#
#######################################################################################

use strict;
use warnings;
our $myBin;
our $myData;
our $myConf;
use Cwd ('getcwd', 'realpath');
use Data::Dumper;

########################################
# globals
my $rootspace = "$myBin/vigorscratch";

use constant CODON_LN => 3;
use constant REF_BUFFER => 125; ## Half of the scanning window on reference sequence (used for finding ribosomal slippage site)
use constant SLIP_RANGE => 100; ## Scanning window on candidate sequence for finding slippage site.
use constant MAX_VARIATION_VAL => 5;
use constant MIN_GAP_LN => 20; ## Minimumlength of a stretch of Ns to be recognized as a sequence gap
use constant DEBUG => 0;

#my $rootspace = "/home/jhoover/VIRAL/VIGOR";

my %codons;
my %splice_pairs;
my %reference_seqs;
my %refmat;
my $refmatdb = "";
my %parameters;
my $configData;
our $flavorData;
my $flavorDB;
my ( $next_sequence_buffer, $next_sequence_eof );
my (
    @sequence_fragments, $sequence_fragcount,
    $frag_sequence_id,   @sequence_gaps
);
my (
    %gene_variations, %required_genes,
    %optional_genes,  %noncanonical_splicing_genes
);
my @roman_numerals =
  ( "i", "ii", "iii", "iv", "v", "vi", "vii", "viii", "ix", "x" );
my $roman_regexp = join( "|", @roman_numerals );

sub initialize_defaults {
    ########################################
    # use default codons/splice sites
    set_default_codons();
    set_default_splice_pairs();

    ########################################
    # set parameter defaults

    # enforce JCVI rules
    set_parameter( "jcvi_rules",     1 );    # 1=on, 0=off
    set_parameter( "use_locus_tags", 1 );    # 1=on, 0=off

    # ribosomal slippage
    set_parameter( "ribosomal_slippage_genes", "" )
      ;    # geneid=spliceform,geneid=spliceform...
    set_parameter( "slippage_motif", '[NT][NT][NT]AAAC' );    # reg. exp
    set_parameter( "slippage_frameshift", -1);
    set_parameter( "slippage_offset", 0);

    # rna editing
    set_parameter( "rna_editing_genes", "" )
      ; # default instruction, geneid=fs/old/new/note/~~geneid=fs/old/new/note/...

    # splicing
    set_parameter( "spliced_genes", "" )
      ;    # geneid=spliceform,geneid=spliceform...
    set_parameter( "min_intron_size",             30 );      # nucleotides
    set_parameter( "max_intron_size",             7000 );    # nucleotides
    set_parameter( "min_exon_size",               3 );       # nucleotides
    set_parameter( "noncanonical_splicing_genes", "" )
      ;    # comma separated, geneid:donor+acceptor;donor+acceptor,...
    set_parameter( "splicesite_extend_range", 0 );    # nucleotides

    # candidate regions
    set_parameter( "reference_db",    "$myData/virus_db" );
    set_parameter( "candidate_refdb", "" );
    set_parameter( "tile_size",       10000000 );             # nucleotides
    set_parameter( "tile_overlap",    1000000 );              # nucleotides
    set_parameter( "candidate_blastopts", "-p blastx -M BLOSUM45 -e 1E-10 -g F -F \"\" -b 2000 -v 0" );
    set_parameter( "candidate_evalue",            "" );
    set_parameter( "min_candidate_pctsimilarity", 60 );
    set_parameter( "min_candidate_sbjcoverage",   33 );
    set_parameter( "candidate_overlap_threshold", 0.8 );    # 80% of query range
    set_parameter( "selectivity",                 0.985 )
      ;    # 0-1, higher value, fewer permutations selected
    set_parameter( "variation", 0.975 )
      ;    # 0-1, higher value, less variation from reference allowed

    # frameshift detection
    set_parameter( "frameshift_sensitivity", 1 )
      ;    # 0=ignore frameshifts, 1=normal, 2=sensitive

    # bl2seq
    set_parameter( "bl2seq_extend_span", 50 );
    set_parameter( "bl2seq_extend_hsp",  25 );

    # genes
    set_parameter( "default_gene_variation", 1 );    # 0 (least) - 5 (most)
    set_parameter( "default_gene_required",  0 )
      ;    # 1 to default genes to required
    set_parameter( "default_gene_optional", 0 )
      ;    # 1 to default genes to optional
    set_parameter( "gene_variations", "" )
      ;    # comma separated, gene_name:value (0-5)
    set_parameter( "required_genes", "" );    # 1=default required
    set_parameter( "optional_genes", "" );    # 1=default optional

    # polyproteins
    set_parameter( "mature_pep_refdb",         "$myData/virus_mp" );
    set_parameter( "mature_pep_mincoverage",   50 );
    set_parameter( "mature_pep_minsimilarity", 40 );
    set_parameter( "mature_pep_minidentity",   25 );

    # pseudogene reporting
    set_parameter( "min_pseudogene_identity",   70 );
    set_parameter( "min_pseudogene_similarity", 80 );
    set_parameter( "min_pseudogene_coverage",   80 );
}

sub turnoff_refmatch_requirements {
    set_parameter( "min_candidate_pctsimilarity", 0 );
    set_parameter( "min_candidate_sbjcoverage",   0 );
    set_parameter( "mature_pep_mincoverage",      0 );
    set_parameter( "mature_pep_minsimilarity",    0 );
    set_parameter( "mature_pep_minidentity",      0 );
    set_parameter( "min_pseudogene_identity",     0 );
    set_parameter( "min_pseudogene_similarity",   0 );
    set_parameter( "min_pseudogene_coverage",     0 );
}

########################################
# manipulate parameter set

sub load_config_db {
    my $db  = "$myBin/VIGOR3.db";
    my $dbh = connectSQLite($db);
    $configData = querySQLArrayHash( $dbh, "select * from vigor_params" );
    our $flavorData = querySQLArrayHash( $dbh, "select * from flavor_db" );
    for my $row (@$flavorData) {
        my $flavor = lc( $$row{flavor} );
        my $db     = "$myData/$$row{db}";
        $$flavorDB{$flavor} = $db;
    }
    $dbh->disconnect;

}

sub set_parameter {
    my ( $param_name, $param_value ) = @_;

    if ( $param_name eq "gene_variations" ) {
        my %empty;
        %gene_variations = %empty;
        my @tmp = split /,/, $param_value;
        for my $variation (@tmp) {
            my ( $gene_name, $variation_value ) = split /:/, $variation;
            if    ( $variation_value < 0 ) { $variation_value = 0 }
            elsif ( $variation_value > MAX_VARIATION_VAL ) { $variation_value = MAX_VARIATION_VAL }
            $gene_variations{$gene_name} = $variation_value;
        }
    }
    elsif ( $param_name eq "required_genes" ) {
        my %empty;
        %required_genes = %empty;
        my @tmp = split /,/, $param_value;
        for my $gene (@tmp) {
            $required_genes{$gene} = 1;
        }
    }
    elsif ( $param_name eq "optional_genes" ) {
        my %empty;
        %optional_genes = %empty;
        my @tmp = split /,/, $param_value;
        for my $gene (@tmp) {
            $optional_genes{$gene} = 1;
        }
    }
    elsif ( $param_name eq "noncanonical_splicing_genes" ) {
        my %empty;
        %noncanonical_splicing_genes = %empty;
        my @tmp = split /,/, $param_value;
        for my $value (@tmp) {
            my ( $gene_name, $splice_values ) = split /:/, $value;
            for my $splice ( split /;/, $splice_values ) {
                my ( $donor, $acceptor ) = split /\+/, uc $splice;
                $noncanonical_splicing_genes{$gene_name}{$donor}{$acceptor} = 1;
            }
        }
    }
    elsif ( $param_name eq "frameshift_sensitivity" ) {
        if ( $param_value != 0 && $param_value != 1 && $param_value != 2 ) {
            die
"\nValue \"$param_value\" is out of range for parameter $param_name, expected 0, 1, or 2\n";
        }
    }

    $parameters{ lc $param_name } = $param_value;

    return;
}

sub get_parameter {
    my ($name) = @_;
    return $parameters{ lc $name };
}

sub show_parameters {
    my ($RPT) = @_;
    print "\nparameters\n";
    if ( defined $RPT ) { print $RPT "\nparameters\n" }
    for my $param ( sort keys %parameters ) {
        print "\t$param\t$parameters{$param}\n";
        if ( defined $RPT ) { print $RPT "\t$param\t$parameters{$param}\n" }
    }
    print "\n";
    if ( defined $RPT ) { print $RPT "\n" }

    return;
}

sub parse_parameters {
    my ($paramstring) = @_;

    #print "param string $paramstring\n";
    for my $param ( split /~~/, $paramstring ) {

        #print "param $param\n";
        if ( $param eq "" ) { next }
        my $brk = index( $param, "=" );
        if ( $brk <= 0 ) {
            die "\ninvalid parameter \"$param\" in \"$paramstring\"\n";
        }
        my $param_name  = substr( $param, 0, $brk );
        my $param_value = substr( $param, $brk + 1 );

        #print "parama name=$param_name  value=$param_value\n";
        set_parameter( $param_name, $param_value );
    }
}

########################################
# manipulate codon set
sub set_default_codons {

    $codons{'GC.'}         = 'A';    # Alanine
    $codons{'TG[TCY]'}     = 'C';    # Cysteine
    $codons{'GA[TCY]'}     = 'D';    # Aspartic Acid
    $codons{'GA[AGR]'}     = 'E';    # Glutamic Acid
    $codons{'TT[TCY]'}     = 'F';    # Phenylalanine
    $codons{'GG.'}         = 'G';    # Glycine
    $codons{'CA[TCY]'}     = 'H';    # Histidine
    $codons{'AT[TCAY]'}    = 'I';    # Isoleucine
    $codons{'AA[AGR]'}     = 'K';    # Lysine
    $codons{'TT[AGR]|CT.'} = 'L';    # Leucine
    $codons{'ATG'}         = 'M';    # Methionine
    $codons{'AA[TCY]'}     = 'N';    # Asparagine
    $codons{'CC.'}         = 'P';    # Proline
    $codons{'CA[AGR]'}     = 'Q';    # Glutamine
    $codons{'CG.|AG[AGR]'} = 'R';    # Arginine
    $codons{'TC.|AG[TCY]'} = 'S';    # Serine
    $codons{'AC.'}         = 'T';    # Threonine
    $codons{'GT.'}         = 'V';    # Valine
    $codons{'TGG'}         = 'W';    # Tryptophan
    $codons{'TA[TCY]'}     = 'Y';    # Tyrosine
    $codons{'TA[AGR]|TGA'} = '*';    # Stop
}

sub get_codons {
    return \%codons;
}

########################################
# manipulate splice donor/acceptor sets
sub set_default_splice_pairs {
    $splice_pairs{GT}{AG} = 3;
    $splice_pairs{GC}{AG} = 1.5;
    $splice_pairs{AT}{AC} = 1.0;
    $splice_pairs{GT}{CA} = 0.5;
    return;
}

sub get_splice_pairs {
    return \%splice_pairs;
}

sub get_splice_donors {
    return \%splice_pairs;
}

sub get_splice_acceptors {
    my %acceptors;
    for my $d ( keys %splice_pairs ) {
        my %acctmp = %{ $splice_pairs{$d} };
        for my $a ( keys %acctmp ) {
            $acceptors{$a}{$d} = $splice_pairs{$d}{$a};
        }
    }
    return \%acceptors;
}

########################################
# reference sequences
sub get_reference_seq {
    my ( $seq_id, $refs ) = @_;

    if ( defined $refs ) {
        return $$refs{$seq_id};
    }
    else {

        #        if ( ! defined $reference_seqs{$seq_id} ) {
        #            print "UNKNOWN REF=$seq_id\n";
        #        }
        return $reference_seqs{$seq_id};
    }
}

sub get_reference_seqs {
    return \%reference_seqs;
}

sub set_reference_seqs {
    my (%refs) = @_;
    %reference_seqs = %refs;
}

sub load_reference_db {

    my $vigorspace   = get_parameter("vigorspace");
    my $reference_db = get_parameter("reference_db");
    my $pepref_fasta = "$vigorspace/pepref.fasta";
    unlink $pepref_fasta;

    # convert gbk/fasta to peptide reference database
    if ( $reference_db =~ /\.gb[kf]$/i ) {
        &runCmd($vigorspace, "$myBin/perl $myBin/dbutils/gb2VigorFasta.pl $reference_db > $pepref_fasta");
        &runCmd($vigorspace, "$myBin/formatdb -i $pepref_fasta -p T -l $vigorspace/formatdb_pepref_fasta.log");

        # update references to point to blastdb
        set_parameter( "reference_db", $pepref_fasta );
        $reference_db = $pepref_fasta;
    }

    # load reference sequences
    else {
        my $cmd = "$myBin/fastacmd -d $reference_db -D 1 2> /dev/null | sed 's/\t/ /g' | sed 's/^>gnl|[^ ]* */>/' > $pepref_fasta";
        &runCmd($vigorspace, $cmd);
    }
    my %refs = loadFasta($pepref_fasta);
    set_reference_seqs(%refs);
    my $refcnt           = scalar keys %refs;
    my $reference_dbsize = get_db_size($reference_db);
    set_parameter( "reference_dbsize", $reference_dbsize );

    # convert gbk/fasta to peptide reference database
    my $mp_reference_db = get_parameter("mature_pep_refdb");
    if ( $mp_reference_db =~ /\.gb[kf]$/i ) {
        my $mpref_fasta = "$vigorspace/matpepref.fasta";
        my $cmd = "$myBin/perl $myBin/dbutils/gb2VigorFasta.pl -m $mp_reference_db > $mpref_fasta";
        &runCmd($vigorspace, $cmd);
        set_parameter( "mature_pep_refdb", $mpref_fasta );
    }
}

########################################
# miscellaneous utilities
sub minval {
    my (@vals) = @_;

    if ( !@vals ) { return undef }
    my $val = $vals[0];
    for my $i ( 1 .. @vals - 1 ) {
        if ( $vals[$i] < $val ) { $val = $vals[$i] }
    }

    return $val;
}

sub maxval {
    my (@vals) = @_;

    if ( !@vals ) { return undef }
    my $val = $vals[0];
    for my $i ( 1 .. @vals - 1 ) {
        if ( $vals[$i] > $val ) { $val = $vals[$i] }
    }

    return $val;
}

sub sign {
    my ($value) = @_;

    if ( $value < 0 ) {
        return -1;
    }
    elsif ( $value > 0 ) {
        return 1;
    }
    else {
        return 0;
    }
}

sub create_workspace {
    my $vigorspace;
    if ( ! -e $rootspace ) { mkdir $rootspace }

    $vigorspace = $rootspace . "/" . rand(1000000);

    #print "vigorspace=$vigorspace\n";
    while ( -e $vigorspace || $vigorspace !~ /\/[0-9]+\.[0-9]+$/) {
        $vigorspace = $rootspace . "/" . rand(1000000);

        #print "vigorspace=$vigorspace\n";
    }

    mkdir $vigorspace;
    set_parameter( "vigorspace", $vigorspace );

    return $vigorspace;
}

########################################
# various comparators used in sorts

sub compare_qry_positions {
    my ( $a, $b ) = @_;

    if ( $$a{orientation} > $$b{orientation} ) {
        return 1;
    }
    elsif ( $$a{orientation} < $$b{orientation} ) {
        return -1;
    }

    my $abegin = $$a{query_left};
    my $aend   = $$a{query_right};
    my $bbegin = $$b{query_left};
    my $bend   = $$b{query_right};

    if ( $$a{orientation} == -1 ) {
        $abegin = -$$a{query_right};
        $aend   = -$$a{query_left};
        $bbegin = -$$b{query_right};
        $bend   = -$$b{query_left};
    }

    if ( $abegin < $bbegin ) {
        return -1;
    }
    elsif ( $abegin > $bbegin ) {
        return 1;
    }
    elsif ( $aend < $bend ) {
        return -1;
    }
    elsif ( $aend > $bend ) {
        return -1;
    }
    else {
        return 0;
    }
}

sub order_left_to_right {
    my ( $a, $b ) = @_;
    if ( $$a{query_id} lt $$b{query_id} ) {
        return -1;
    }
    if ( $$a{query_left} < $$b{query_left} ) {
        return -1;
    }
    elsif ( $$a{query_left} > $$b{query_left} ) {
        return 1;
    }
    elsif ( $$a{query_right} < $$b{query_right} ) {
        return -1;
    }
    elsif ( $$a{query_right} > $$b{query_right} ) {
        return 1;
    }
    elsif ( $$a{orientation} > $$b{orientation} ) {
        return -1;
    }
    elsif ( $$a{orientation} < $$b{orientation} ) {
        return 1;
    }
    else {
        return 0;
    }
}

sub compare_gene_ids {
    my ( $a, $b ) = @_;
    return compare_ids( $$a{gene_id}, $$b{gene_id} );
}

sub compare_ids {
    my ( $aid, $bid ) = @_;

    my (@atmp) = split /\./, $aid;
    my $agene_num = pop @atmp;
    my $agenome_id = join( ".", @atmp );
    my (@btmp) = split /\./, $bid;
    my $bgene_num  = pop @btmp;
    my $bgenome_id = join( ".", @btmp );

    my $asuffix = " ";
    my $bsuffix = " ";
    if ( $agene_num =~ /[a-z]/i ) {
        $asuffix = $agene_num;
        $asuffix   =~ s/[^a-z]//gi;
        $agene_num =~ s/[a-z]//gi;
    }
    if ( $bgene_num =~ /[a-z]/i ) {
        $bsuffix = $bgene_num;
        $bsuffix   =~ s/[^a-z]//gi;
        $bgene_num =~ s/[a-z]//gi;
    }

    if ( $agenome_id ne $bgenome_id ) { return $agenome_id cmp $bgenome_id }
    elsif ( $agene_num != $bgene_num ) { return $agene_num <=> $bgene_num }
    else { return $asuffix cmp $bsuffix }
}

sub compare_gene_positions {
    my ( $a, $b, $test_pseudo ) = @_;
    if ( ! defined $test_pseudo ) { $test_pseudo = 0 }

    if ( $test_pseudo ) {
        if ( $$a{is_pseudogene} && ! $$b{is_pseudogene} ) {
            return 1;
        }
        elsif ( $$b{is_pseudogene} && ! $$a{is_pseudogene} ) {
            return -1;
        }
    }

    if ( $$a{start_site} < $$b{start_site} ) {
        return -1;
    }
    elsif ( $$a{start_site} > $$b{start_site} ) {
        return 1;
    }

    if ( $$a{stop_site} < $$b{stop_site} ) {
        return -1;
    }
    elsif ( $$a{stop_site} > $$b{stop_site} ) {
        return 1;
    }

    if ( $$a{orientation} < $$b{orientation} ) {
        return -1;
    }
    elsif ( $$a{orientation} > $$b{orientation} ) {
        return 1;
    }
    else {
        return 0;
    }
}


########################################
# sequences and fasta files

# extract subsequence
# reverse complements if begin>end
sub subsequence {
    my ( $sequence, $begin, $end ) = @_;

    # establish boundaries
    my $max = length($sequence);

    my ( $b, $e ) = ( $begin, $end );
    if ( $b > $e ) {
        ( $b, $e ) = ( $end, $begin );
    }

    # if subsequence requested is outside the known sequ8ence, return Ns
    if ( $e < 1 || $b > $max ) {
        return lpad( "N", $e - $b + 1, "N" );
    }

    # pad with Ns if subsequence starts before known sequence
    my $lpad = "";
    while ( $b < 1 ) { $lpad .= "N"; $b++ }

    # pad with Ns if subsequence ends after known sequence
    my $rpad = "";
    while ( $e > length($sequence) ) { $rpad .= "N"; $e-- }

    # extract subsequence from known portion of sequence
    my $subseq = $lpad . substr( $sequence, $b - 1, $e - $b + 1 ) . $rpad;

    # return reverse strand if end > begin
    if ( $begin > $end ) {
        return reverse_complement($subseq);
    }

    # otherwise return forward strand
    else {
        return $subseq;
    }
}

sub reverse_complement {
    my ($seq) = @_;

    if ( !defined $seq ) { return "" }

    my $revcomp = reverse $seq;
    $revcomp =~
      tr/AGCTMRWSYKVHDBNagctmrwsykvhdbn/TCGAKYWSRMBDHVNtcgakywsrmbdhvn/;

    return $revcomp;
}

# translate DNA sequence to amino acid sequence
sub DNA2AA {
    my ( $seq, $codon_start, $alternate_startcodon, $stopcodon_readthru,
        $alternate_stopcodon )
      = @_;
    if ( !defined $seq )         { return "" }
    if ( !defined $codon_start ) { $codon_start = 1 }

    my $len = length($seq);
    if ( !$len ) { return "" }

    my $start_site = $codon_start - 1;
    my $protein;
    while ( $start_site <= $len - 3 ) {
        my $codon_dna = substr( $seq, $start_site, 3 );
        my $AA = &codon2aa($codon_dna);
        $protein .= $AA;
        $start_site = $start_site + 3;
    }
    if ( defined $alternate_startcodon ) {
        $protein = "M" . substr( $protein, 1 );
    }
    if ( defined $stopcodon_readthru ) {
        my $pos = $$stopcodon_readthru{aa_position};
        my $aa  = $$stopcodon_readthru{aa};

        #print "protein $protein\n"
        #. "length " . length( $protein ) . "\n"
        #. "exception $pos: $aa\n";
        $protein =
          substr( $protein, 0, $pos - 1 ) . $aa . substr( $protein, $pos );
    }
    if ( defined $alternate_stopcodon ) {
        $protein = substr( $protein, 0, length($protein) - 1 ) . "*";
    }

    return $protein;
}

sub codon2aa {
    my ($codon) = @_;

    for my $key ( keys %codons ) {
        if ( $codon =~ /$key/i ) {
            return $codons{$key};
        }
    }

    # probably contain an IUPAC ambiguity code;
    return "X";
}

# read fasta file into hash;
sub loadFasta {
    my ($fasta) = @_;
    
    if ( ! -e $fasta ) {
        if ( -e "$fasta.pal" || -e "$fasta.pin" || -e "$fasta.nin" ) {
            my $vigorspace = get_parameter( "vigorspace" );
            my $tmpfasta = "$vigorspace/loadFasta.tmp";
            my $cmd = "$myBin/fastacmd -d $fasta -D 1 2> /dev/null | sed 's/\t/ /g' | sed 's/^>gnl|[^ ]* */>/' > $tmpfasta";
            &runCmd($vigorspace, $cmd);
            $fasta = $tmpfasta;            
        }
    }

    my %contents;
    open( FASTA, "<$fasta" ) || die "\nCould not read fasta file \"$fasta\"\n";

    my $sequence = "";
    my $defline;
    while ( my $line = <FASTA> ) {
        chomp $line;
        $line =~ s/\r//g;
        $line =~ s/\t/ /g;
        if ( $line =~ /^ *$/ ) { next }
        if ( $line =~ /^>/ )   {
            if ( length($sequence) ) {
                my %entry;
                $entry{defline}  = $defline;
                $entry{sequence} = $sequence;
                $entry{seqlen}   = length($sequence);
                if ( $sequence =~ /\*$/ ) { $entry{seqlen}-- }
                my ($id) = get_defline_id($defline);
                $entry{id} = $id;
                $contents{$id} = \%entry;
            }
            $defline = $line;
            $defline =~ s/^ *> *//;
            $sequence = "";
        }
        else {
            $sequence .= $line;
        }
    }
    close FASTA;
    if ( length($sequence) ) {
        my %entry;
        $entry{defline}  = $defline;
        $entry{sequence} = $sequence;
        $entry{seqlen}   = length($sequence);
        if ( $sequence =~ /\*$/ ) { $entry{seqlen}-- }
        my ($id) = get_defline_id($defline);
        $entry{id} = $id;
        $contents{$id} = \%entry;
    }

    return %contents;
}

# fetch next sequence from fasta file
sub next_sequence {
    if ($next_sequence_eof) { return undef }
    my ($FASTA) = @_;

    my %entry;
    $entry{id}          = "";
    $entry{defline}     = "";
    $entry{sequence}    = "";
    $entry{frag_offset} = 0;

    if ( defined $next_sequence_buffer && length($next_sequence_buffer) ) {
        $next_sequence_buffer =~ s/^ *> *//;
        $entry{defline} = $next_sequence_buffer;
        ( $entry{id} ) = split / /, $entry{defline};
    }

    while ( $next_sequence_buffer = <$FASTA> ) {
        chomp $next_sequence_buffer;
        $next_sequence_buffer =~ s/\r//g;
        $next_sequence_buffer =~ s/\t/ /g;
        if ( $next_sequence_buffer =~ /^ *$/ ) { next }
        if ( $next_sequence_buffer =~ /^>/ )   {
            my $seqlen = length( $entry{sequence} );
            if ($seqlen) {
                $entry{seqlen} = $seqlen;
                set_genome_gaps( \%entry );
                return \%entry;
            }
            $next_sequence_buffer =~ s/^ *> *//;
            $entry{defline} = $next_sequence_buffer;
            ( $entry{id} ) = split / /, $entry{defline};
            $entry{sequence} = "";
        }
        else {
            $entry{sequence} .= $next_sequence_buffer;
        }
    }

    $next_sequence_eof = 1;
    my $seqlen = length( $entry{sequence} );
    if ($seqlen) {
        $entry{seqlen} = $seqlen;
        set_genome_gaps( \%entry );
        return \%entry;
    }
    return undef;
}

sub set_genome_gaps {
    my ($genome) = @_;

    if ( $$genome{seqlen} ) {
        if ( !$$genome{is_complete} ) {
            for my $i ( -2 .. 0 ) {
                $$genome{in_gap}{$i} = -1;
            }
            for my $i ( $$genome{seqlen} + 1 .. $$genome{seqlen} + 3 ) {
                $$genome{in_gap}{$i} = -1;
            }
        }

        $$genome{gaps} = find_ngaps( $$genome{sequence} );

        if ( defined $$genome{gaps} ) {
            my $gapno = 0;
            for my $gap ( @{ $$genome{gaps} } ) {
                for my $i ( $$gap{begin} .. $$gap{end} ) {
                    $$genome{in_gap}{$i} = $gapno;
                }
                $gapno++;
            }
        }
    }
}

sub find_ngaps {
    my ($sequence) = @_;

    my @gaps;
    my $lastgap;
    for my $gap ( find_regexp( '[Nn]{'.MIN_GAP_LN.',}', $sequence ) ) {
        if ( defined $lastgap && $$gap{begin} < $$lastgap{end} + 6 ) {
            $$lastgap{end} = $$gap{end};
        }
        else {
            push @gaps, $gap;
            $lastgap = $gap;
        }
    }

    return \@gaps;
}

sub in_gap {
    my ( $genome, $begin, $end ) = @_;

    if ( $begin > $end ) {
        for my $pos ( $end .. $begin ) {
            if ( defined $$genome{in_gap}{$pos} ) { return 1 }
        }
    }
    else {
        for my $pos ( $begin .. $end ) {
            if ( defined $$genome{in_gap}{$pos} ) { return 1 }
        }
    }
    return 0;
}

sub find_gaps {
    my ( $genome, $beg, $end, $trim ) = @_;
    if ( ! defined $trim ) { $trim = 0 }

    my @gaps;
    if ( !defined $$genome{gaps} ) { return @gaps }

    my ( $left, $right ) = ( $beg, $end );
    if ( $left > $right ) { ( $left, $right ) = ( $right, $left ) }

    for my $gap ( @{ $$genome{gaps} } ) {
        if ( $trim ) {
            if ( $left <= $$gap{end} && $right >= $$gap{begin} ) {
                my $trimmedgap = {
                    begin => $left > $$gap{begin} ? $left : $$gap{begin},
                    end => $right < $$gap{end} ? $right : $$gap{end}
                };
                push @gaps, $trimmedgap;
            }
        }
        else {
            if ( $left <= $$gap{begin} && $right >= $$gap{end} ) {
                push @gaps, $gap;
            }
        }
    }
    return @gaps;
}

# parse id from defline
sub get_defline_id {
    my ($defline) = @_;

    my ($rawid) = split / /, $defline;
    $rawid =~ s/[,\| ]*$//;

    my @id = split( /\|/, $rawid );
    if ( scalar @id < 2 ) { return $rawid }

    if ( $id[0] eq "gnl" && @id > 2 ) {
        shift @id;
        shift @id;
        $rawid = join( "|", @id );
        if ( @id <= 2 ) { return $rawid }
    }

    my $i = 0;
    while ( $i < @id - 1 ) {
        if (
            index( ".gi.gb.rf.emb.dbj.pir.prf.sp.ref.",
                "." . lc( $id[$i] ) . "." ) >= 0
          )
        {
            return lc( $id[$i] ) . "|" . $id[ $i + 1 ];
        }
        $i++;
    }
    return $rawid;
}

# return blast db size
sub get_db_size {
    my ($dbpath) = @_;

    my $vigorspace = get_parameter("vigorspace");
    my $cmd = "$myBin/fastacmd -d $dbpath -I 2> /dev/null";
    my $results = &runCmdAndGetResults($vigorspace, $cmd);
    my @tmp = split(/\n/, $results);
     
    for my $dbsize (@tmp) {
        if ( $dbsize =~ /total letters/ ) {
            $dbsize =~ s/total letters.*//;
            $dbsize =~ s/^.*; //;
            $dbsize =~ s/[^0-9]//g;
            return $dbsize;
        }
    }

    return undef;
}

########################################
# blast reference da    tabase against genome
#
sub blast_sequence {
    my ( $genome, $database, $options, $dbg ) = @_;
    if ( !defined $dbg ) { $dbg = DEBUG }

    my $vigorspace = get_parameter("vigorspace");

    my @hsps;
    my $tmpfile = "$vigorspace/seq.fasta";
    my $xmlfile = "$vigorspace/seq.xml";
    unlink $tmpfile;
    unlink $xmlfile;

    my %uniqhits;
    open( TMP, ">$tmpfile" );
    print TMP ">$$genome{defline}\n$$genome{sequence}";
    close TMP;
    my $blastcmd = "$myBin/blastall -i $tmpfile -d $database $options -m 7 1> $xmlfile 2> /dev/null";
    &runCmd($vigorspace, $blastcmd);

#print "$blastcmd\n";
#print "\n" . `cat $tmpfile` . "\n";
    my @seqhsps = sort { order_left_to_right( $a, $b ) } parse_blastxml( $xmlfile, 0 );
    unlink $tmpfile;
    unlink $xmlfile;

    return @seqhsps;
}

# parse blast result xml
sub parse_blastxml {
    my ( $xmlfile, $get_alignment_strings ) = @_;
    my $dbg = DEBUG;

    #    my $euler = 2.718281828;
    my %xmlattr;
    $xmlattr{"Iteration_query-def"} = "query_definition";
    $xmlattr{"Iteration_query-len"} = "query_length";
    $xmlattr{"Hit_id"}              = "subject_acc";
    $xmlattr{"Hit_def"}             = "subject_definition";
    $xmlattr{"Hit_len"}             = "subject_length";
    $xmlattr{"Hsp_bit-score"}       = "bit_score";

    #    $xmlattr{"Hsp_score"} = "hsp_score";
    $xmlattr{"Hsp_evalue"}     = "evalue";
    $xmlattr{"Hsp_query-from"} = "query_left";
    $xmlattr{"Hsp_query-to"}   = "query_right";
    $xmlattr{"Hsp_hit-from"}   = "subject_begin";
    $xmlattr{"Hsp_hit-to"}     = "subject_end";

    $xmlattr{"Hsp_query-frame"} = "query_frame";
    $xmlattr{"Hsp_hit-frame"}   = "subject_frame";
    $xmlattr{"Hsp_identity"}    = "num_identical";
    $xmlattr{"Hsp_positive"}    = "num_similar";
    $xmlattr{"Hsp_gaps"}        = "num_gaps";
    $xmlattr{"Hsp_align-len"}   = "alignment_length";
    $xmlattr{"Hsp_qseq"}        = "query_alignstr";
    my $last_attribute = "query_alignstr";

    # bit flags signal which alignment strings to retrieve:
    # 0=none, 1 = query, 2 = subject, 4 = midline
    # (i.e. 3 = query & subject alignment strings)
    #
    # note: we always EXTRACT the query string so we can
    # count stop codons, but we only keep it if flagged
    if ($get_alignment_strings) {
        if ( $get_alignment_strings / 2 % 2 ) {
            $last_attribute = "subject_alignstr";
            $xmlattr{"Hsp_hseq"} = $last_attribute;
        }
        if ( $get_alignment_strings / 4 % 2 ) {
            $last_attribute = "midline_alignstr";
            $xmlattr{"Hsp_midline"} = $last_attribute;
        }
    }

    my @hits;
    my $hit;
    my $query_id;
    my $query_definition;
    my $query_length;
    my $subject_id;
    my $subject_name;
    my $subject_definition;
    my $subject_length;

    my $linenum = 0;
    open( XML, "<$xmlfile" );

    #print "$xmlfile\n";
    while ( my $line = <XML> ) {

        if ($dbg) { print $line }

        #print $line;
        $linenum++;
        chomp $line;
        $line =~ s/[\n\r]/ /g;
        $line =~ s/^\s+//;
        $line =~ s/\s+$//;
        $line =~ s/&lt;/</g;
        $line =~ s/&gt;/>/g;
        $line =~ s/&amp;/&/g;
        $line =~ s/&quot;/"/g;
        $line =~ s/&apos;/'/g;

        # check XML tag
        if ( substr( $line, 0, 1 ) eq "<" ) {
            $line = substr( $line, 1 );
            my ($xml_tag) = split( ">", $line );
            my $attrname = $xmlattr{$xml_tag};
            if ( !defined $attrname ) { next }

            # get tag value
            my $attrval;
            $line = substr( $line, length($xml_tag) + 1 );
            my $close_tag = "</" . $xml_tag . ">";
            my $eod = index( $line, $close_tag );
            if ( $eod == -1 ) {
                $attrval = $line;
            }
            else {
                $attrval = substr( $line, 0, $eod );
            }

            # save tag value
            if ( $attrname eq "subject_definition" ) {
                $subject_definition = $attrval;
                $subject_definition =~ s/^[ >]*//;
                $subject_definition =~ s/ *$//;
                if ( defined $$hit{subject_acc} ) {

                    #print "ACCESSION=$$hit{subject_acc}\n";
                    if ( substr( $$hit{subject_acc}, 0, 4 ) ne "gnl|" ) {
                        $subject_definition =
                          "$$hit{subject_acc} $subject_definition";
                    }

                    #print "DEFLINE=$subject_definition\n";
                }
                $$hit{subject_definition} = $subject_definition;
                $subject_id               = get_defline_id($subject_definition);
                $subject_name             = get_reference_name($subject_id);
                $$hit{subject_id}         = $subject_id;
                $$hit{subject_name}       = $subject_name;

            }
            elsif ( $attrname eq "subject_length" ) {
                $subject_length = $attrval;
                $$hit{subject_length} = $subject_length;
            }
            elsif ( $attrname eq "query_definition" ) {
                $query_definition = $attrval;
                $query_definition =~ s/^[ >]*//;
                $query_definition =~ s/ *$//;
                $$hit{query_definition} = $query_definition;
                $query_id               = get_defline_id($query_definition);
                $$hit{query_id}         = $query_id;
            }
            elsif ( $attrname eq "query_length" ) {
                $query_length = $attrval;
                $$hit{query_length} = $query_length;
            }
            else {
                $$hit{$attrname} = $attrval;
            }
            if ( !defined $attrval ) {
                my $tmpcmd = "tail +n " . ( $linenum - 5 ) . " | head -n 20 ";

                #print `$tmpcmd` . "\n";
            }

            # save hit
            if ( $attrname eq $last_attribute ) {
                my $orientation = 1;
                if ( defined $$hit{query_frame} ) {
                    if ( $$hit{query_frame} < 0 ) { $orientation = -1 }
                }
                elsif ( defined $$hit{subject_frame} ) {
                    if ( $$hit{subject_frame} < 0 ) { $orientation = -1 }
                }
                if ( $$hit{query_left} > $$hit{query_right} ) {
                    $orientation = -$orientation;
                    my $tmp = $$hit{query_left};
                    $$hit{query_left}  = $$hit{query_right};
                    $$hit{query_right} = $tmp;
                }
                if ( $$hit{subject_begin} > $$hit{subject_end} ) {
                    $orientation = -$orientation;
                    my $tmp = $$hit{subject_begin};
                    $$hit{subject_begin} = $$hit{subject_end};
                    $$hit{subject_end}   = $tmp;
                }
                $$hit{orientation} = $orientation;
                calculate_blast_percentages($hit);

                $$hit{has_frameshift} = 0;

                #print "BLAST HSP\n";
                #print_blasthits( 0,  $hit );
                # save this hit and start a new one
                push @hits, $hit;

                my %newhit;
                $newhit{query_id}           = $query_id;
                $newhit{query_definition}   = $query_definition;
                $newhit{query_length}       = $query_length;
                $newhit{subject_id}         = $subject_id;
                $newhit{subject_name}       = $subject_name;
                $newhit{subject_definition} = $subject_definition;
                $newhit{subject_length}     = $subject_length;
                $hit                        = \%newhit;
            }
        }
    }

    close XML;
    return @hits;
}


sub make_hsp {
    my ( $genome, $template, $sbjleft, $sbjright, $qryleft, $qryright, $numid, $numsim,  $frame ) = @_;

    my %new = %$template;
    delete $new{hsps};
    delete $new{query_alignstr};
    delete $new{subject_alignstr};
    delete $new{midline_alignstr};

    $new{subject_begin}    = $sbjleft;
    $new{subject_end}      = $sbjright;
    $new{subject_coverage} = $sbjright - $sbjleft + 1;
    $new{longest_orf}      = $new{subject_coverage};
    $new{subject_alignstr} = lpad( "X", $new{subject_coverage}, "X" );

    $new{query_left}     = $qryleft;
    $new{query_right}    = $qryright;
    $new{query_coverage} = $qryright - $qryleft + 1;
    $new{query_alignstr} =
      substr( $$genome{sequence}, $qryleft - 1, $new{query_coverage} );

    if ( defined $frame ) { $new{query_frame} = $frame }

    $new{alignment_length} =
      maxval( $new{subject_coverage}, $new{query_coverage} / 3 );
    $new{num_identical} = $numid;
    $new{num_similar}   = $numsim;

    $new{bit_score} =
      minval( $new{subject_coverage}, $new{query_coverage} / 3.0 ) / 4.0 +
      $numsim / 2.0 + $numid;
    if ( $new{bit_score} < 2 ) { $new{bit_score} = 2 }

    $new{evalue} = 0.5;

    $new{has_frameshift} = 0;

    calculate_blast_percentages( \%new );
    
    return \%new;
}

# calculate %identity, %coverage, etc.
sub calculate_blast_percentages {
    my ($hsp) = @_;

    my $frameshift = $$hsp{has_frameshift};
    if ($frameshift) {
        for my $attr (
            "subject_coverage", "query_coverage",
            "num_identical",    "num_similar",
            "longest_orf",      "bit_score"
          )
        {
            $$hsp{$attr} = int( 0.90 * $$hsp{$attr} );
        }
        $$hsp{evalue} = ( $$hsp{evalue} + 1e-70 ) / 0.90;
    }
    if ( !defined $$hsp{subject_coverage} ) {
        $$hsp{subject_coverage} = $$hsp{subject_end} - $$hsp{subject_begin} + 1;
    }
    if ( !defined $$hsp{query_coverage} ) {
        $$hsp{query_coverage} = $$hsp{query_right} - $$hsp{query_left} + 1;
    }

    $$hsp{pct_identity} =
      int( 1000. * $$hsp{num_identical} / ( $$hsp{alignment_length} ) ) / 10.;
    if ( $$hsp{pct_identity} > 100 ) { $$hsp{pct_identity} = 100 }

    $$hsp{pct_similarity} =
      int( 1000. * $$hsp{num_similar} / ( $$hsp{alignment_length} ) ) / 10.;
    if ( $$hsp{pct_similarity} > 100 ) { $$hsp{pct_similarity} = 100 }

    $$hsp{pct_qcoverage} =
      int( 1000. * $$hsp{query_coverage} / $$hsp{query_length} ) / 10.;
    if ( $$hsp{pct_qcoverage} > 100 ) { $$hsp{pct_qcoverage} = 100 }
    $$hsp{pct_scoverage} =
      int( 1000. * $$hsp{subject_coverage} / $$hsp{subject_length} ) / 10.;
    if ( $$hsp{pct_scoverage} > 100 ) { $$hsp{pct_scoverage} = 100 }
    $$hsp{pct_coverage} = maxval( $$hsp{pct_qcoverage}, $$hsp{pct_scoverage} );

    $$hsp{query_gaps} =
      int( $$hsp{alignment_length} - $$hsp{query_coverage} / 3.0 );
    $$hsp{subject_gaps} = $$hsp{alignment_length} - $$hsp{subject_coverage};

    if ( defined $$hsp{query_alignstr} ) {
        my @orfs = sort { length($b) <=> length($b) } split /\*+/,
          $$hsp{query_alignstr};
        $$hsp{query_stops} = @orfs - 1;
        $$hsp{longest_orf} = length( $orfs[0] );
    }
    else {
        $$hsp{query_stops} = 0;
        $$hsp{longest_orf} = $$hsp{query_coverage} / 3;
    }

    $$hsp{vigor_pctsimilarity} = calculate_vigor_pctsimilarity($hsp);
    $$hsp{vigor_matchwgt}      = $$hsp{vigor_pctsimilarity} / 100.0 *
      maxval( $$hsp{alignment_length}, $$hsp{subject_length} );
}


sub calculate_vigor_pctsimilarity {
    my ($hit) = @_;

    #    my $pct = ( $$hit{pct_identity} + $$hit{pct_similarity} ) / 2.0;

    # match weight
    my $mismatch = $$hit{subject_coverage} - $$hit{num_similar};
    my $similar  = $$hit{num_similar} - $$hit{num_identical};
    my $weight   = 0.10 * $mismatch + 0.60 * $similar +
      $$hit{num_identical};    # - $nomatch * 0.10;

    # as percentage of subject length
    my $pct = 100.0 * $weight /
      maxval( $$hit{alignment_length}, $$hit{subject_length} );

    if ( $pct > 100.0 ) { $pct = 100; }

    return $pct;
}

########################################
# find genomic regions likely to contain genes
sub find_candidates {
    my ($genome, $ignore_list ) = @_;
    my $dbg = DEBUG;
    my $verbose = get_parameter("verbose");

    #my $frag_blastopts = get_parameter("candidate_blastopts");
    my $frag_blastopts = candidate_blastopts( $$genome{seqlen} );
    my $candidate_db = get_parameter("reference_db");
    my $reference_dbsize = get_parameter("reference_dbsize");
    my $min_candidate_sbjcoverage = get_parameter("min_candidate_sbjcoverage");
    my $min_candidate_pctsimilarity = get_parameter("min_candidate_pctsimilarity");
    my $candidate_evalue_cutoff = get_parameter("candidate_evalue");
    if ( $candidate_evalue_cutoff =~ /\w/ ) {
        $frag_blastopts .= " -e $candidate_evalue_cutoff";
    }

    # blast genomic sequence against reference proteins
    my @hits;
    my $top;
    {
        my @rawhsps = blast_sequence( $genome, $candidate_db, $frag_blastopts );
        if ($verbose) {
            print "\nRAW CANDIDATE HSPS\n";
            print_blasthits( 0, @rawhsps );
        }        

        my @hsps;
        for my $hsp ( @rawhsps ) {
            if ( defined $ignore_list && index( ",$ignore_list,", ",$$hsp{subject_id}," ) >= 0 ) { next }
            if ( $$hsp{pct_similarity} >= $min_candidate_pctsimilarity ) { push @hsps, $hsp }            
        }
        
        # join hsps into hits
        if ( @hsps ) {
            my @rawhits = sort { order_left_to_right( $a, $b ) } best_gene_blasthits( $genome, \@hsps );
            if ($dbg) {
                print "\nRAW CANDIDATE HITS\n";
                print_blasthits( 0, @rawhits );
            }
            for my $hit (@rawhits) {

                # screen hits coverage
                if ( check_coverage( $min_candidate_sbjcoverage, $hit, $genome ) ) {
                    push @hits, $hit;
                }
            }
        }
    }
    if ($verbose) {
        print "\nCOVERED CANDIDATE HITS\n";
        print_blasthits( 0, @hits );
    }

    # resolve overlapping hits
    my $i = 0;
    while ( $i < @hits ) {
        my $ihit = $hits[$i];
        my $j    = $i + 1;
        while ( $j < @hits ) {
            my $jhit = $hits[$j];
            if ( $$ihit{bin} eq $$jhit{bin} || shared_cds( $ihit, $jhit ) )
            {
                $j++;
            }
            elsif ( candidates_overlap( $ihit, $jhit ) ) {
                my $comparison = default_candidate_comparison( $ihit, $jhit );
                if ( $comparison > 0 ) {
                    print "\n    Resolved Overlapping genes\n";
                    print_blasthits( 4, $ihit, $jhit );
                    splice @hits, $j, 1;
                    next;
                }
                elsif ( $comparison < 0 ) {
                    print "\n    Resolved Overlapping genes\n";
                    print_blasthits( 4, $jhit, $ihit );
                    splice @hits, $i, 1;
                    $i--;
                    last;
                }
                else {
                    print "\n    Resolved Overlapping genes\n";
                    print_blasthits( 4, $jhit, $ihit );
                    splice @hits, $i, 1;
                    $i--;
                    last;
                }
            }
            else {
                $j++;
            }
        }
        $i++;
    }

    @hits = sort { order_left_to_right( $a, $b ) } @hits;

    if ($verbose) {
        print "\nFINAL CANDIDATES\n";
        print_blasthits( 0, @hits );
    }

    my @rawgenes;
    {
        my $bins;
        for my $hit ( @hits ) {
            if ( ! exists $bins->{ $hit->{bin} } ) {
                my $rawgene;
                my $gene_name = get_reference_name( $hit->{subject_id} );
                $$rawgene{gene_name} = $gene_name;
                my @permutations = ( $hit );
                $$rawgene{query_id}     = $hit->{query_id};
                $$rawgene{orientation}  = $hit->{orientation};
                $$rawgene{query_left}   = $hit->{query_left};
                $$rawgene{query_right}  = $hit->{query_right};
                $bins->{ $hit->{bin} } = $rawgene;
            }
            push @{ $bins->{ $hit->{bin} }->{permutations} }, $hit
        }
        @rawgenes = values %$bins;
    }
    return sort { order_left_to_right( $a, $b ) } @rawgenes;
}

sub best_gene_blasthits {
    my ( $genome, $hsps ) = @_;
    my $dbg = DEBUG;
    #if ( get_parameter( "reference_db" ) =~ /yfv/i ) { $dbg = 1 }

    my @genehits = best_subject_hits( $genome, $hsps );
    my %bins;
    for my $genehit (@genehits) {
        my $geneid  = get_reference_name( $$genehit{subject_id} );
        my $binid = "$geneid|$$genehit{bin}";

        my $score = score_permutation($genehit);
        $$genehit{permutation_score} = $score;
        if ( !exists $bins{$binid} ) {
            my @permutations = ($genehit);
            $bins{$binid}{candidates} = \@permutations;
        }
        else {
            push @{ $bins{$binid}{candidates} }, $genehit;
        }
    }

    for my $binid ( keys %bins ) {
        my @candidates = @{ $bins{$binid}{candidates} };
        if ($dbg) {
            print "\nbin $binid candidates\n";
            print_blasthits( 0, @candidates );
        }

        # select best candidates
        my %selections;

        # select by similarity
        {
            if ($dbg) { print "\n  selecting by similarity\n" }
            my $reqdvalue;
            for my $candidate (
                sort { $$b{permutation_score} <=> $$a{permutation_score} }
                @candidates )
            {
                if ( !defined $reqdvalue ) {
                    $reqdvalue = 0.98 * $$candidate{permutation_score};
                    if ($dbg) { print "  required value: $reqdvalue\n" }
                }
                elsif ( $$candidate{permutation_score} < $reqdvalue ) {
                    last;
                }
                my $descriptor = get_candidate_descriptor( $genome, $candidate );
                if ( exists $selections{$descriptor} ) { next }
                $selections{$descriptor} = $candidate;
                if ($dbg) {
                    print
"  descriptor: $descriptor  candidate: $$candidate{subject_id}  value: $$candidate{permutation_score}  reqd: $reqdvalue\n";
                }
            }
        }

        # select by weight
        {
            if ($dbg) { print "\n  selecting by weight\n" }
            my $reqdvalue;
            for my $candidate (
                sort { $$b{vigor_matchwgt} <=> $$a{vigor_matchwgt} }
                @candidates )
            {
                if ( !defined $reqdvalue ) {
                    $reqdvalue = 0.98 * $$candidate{vigor_matchwgt};
                    if ($dbg) { print "required value: $reqdvalue\n" }
                }
                elsif ( $$candidate{vigor_matchwgt} < $reqdvalue ) {
                    last;
                }
                my $descriptor = get_candidate_descriptor( $genome, $candidate );
                if ( exists $selections{$descriptor} ) { next }
                $selections{$descriptor} = $candidate;
                if ($dbg) {
                    print
"  descriptor: $descriptor  candidate: $$candidate{subject_id}  value: $$candidate{permutation_score}  reqd: $reqdvalue\n";
                }
            }
        }

        # select by coverage
        {
            if ($dbg) { print "\n  selecting by coverage\n" }
            my $reqdvalue;
            for my $selected ( sort{ $$b{pct_scoverage} <=> $$a{pct_scoverage} } values %selections ) {
                $reqdvalue = $$selected{pct_scoverage};
                last;
            }
            if ($dbg) { print "  required value: $reqdvalue\n" }
            for my $candidate ( sort { select_by_coverage_sort( $a, $b ) } @candidates ) {
                if ( $$candidate{pct_scoverage} > $reqdvalue ) {
                    my $descriptor = get_candidate_descriptor( $genome, $candidate);
                    if ( ! exists $selections{$descriptor} ) {
                        $selections{$descriptor} = $candidate;
                        if ($dbg) {
                            print "  descriptor: $descriptor  candidate: $$candidate{subject_id}  value: $$candidate{permutation_score}  reqd: $reqdvalue\n";
                        }
                        last;
                    }
                }
                else {
                    last;
                }
            }
        }

        my @selected = values %selections;
        if ($dbg) {
            print "\n  final candidates\n";
            print_blasthits( 2, @selected );
        }
        $bins{$binid}{selected} = \@selected;
    }

    my @besthits;
    for my $binid ( sort { $a cmp $b } keys %bins ) {
        push @besthits, @{ $bins{$binid}{selected} };
    }
    return @besthits;
}

sub select_by_coverage_sort {
    my ( $a, $b ) = @_;
    
    if ( $$a{pct_scoverage} > $$b{pct_scoverage} ) {
        return -1;
    }    
    if ( $$a{pct_scoverage} < $$b{pct_scoverage} ) {
        return 1;
    }
    my $aval = $$a{permutation_score} * $$a{vigor_matchwgt};
    my $bval = $$b{permutation_score} * $$b{vigor_matchwgt};
    
    return $bval <=> $aval;
}

sub score_permutation {
    my ($permutation) = @_;

    my $score =
      $$permutation{pct_similarity} / 100.0 * $$permutation{pct_identity} /
      100.0 * sqrt( $$permutation{pct_scoverage} / 100.0 );

    return $score;
}

sub best_subject_hits {
    my ( $genome, $hsps, $query_ratio, $matpep ) = @_;
    if ( !defined $query_ratio ) { $query_ratio = 3 }
    if ( !defined $matpep )      { $matpep      = 0 }
    my $dbg = DEBUG;

#if ( `whoami` =~ /jhoover/ ) { $dbg = 1 }
#if ( get_parameter( "reference_db" ) =~ /yfv/i ) { $dbg = 1 }
#if ( get_reference_name( $$hsps[0]{subject_id} ) eq "SP" ) { $dbg = 1 }
#if ( get_reference_name( $$hsps[0]{subject_id} ) =~ /^(pTP|L1-3|E2B)$/ ) { $dbg = 1 }

    my @besthits = ();
    my @ngaps;
    if ( defined $$genome{ngaps} ) { push @ngaps, @{ $$genome{ngaps} } }

    my %subjects;
    my %bins;
    if ( !$matpep ) {
        my %genes;
        my $subject_id = "";
        my $structure;
        my $genename;
        my $product;
        #my $trims = 0;
        for my $hsp ( sort { $a->{subject_id} cmp $b->{subject_id} } @$hsps ) {

            # trim HSPs to exon structure (sometimes blast overextends the HSP)
            if ( $hsp->{subject_id} ne $subject_id ) {
                $subject_id = $hsp->{subject_id};
                $structure = get_gene_structure( $subject_id );
                $genename = get_reference_name( $subject_id );
                $product = get_reference_product( $subject_id );
            }
            my $expected = expected_exon( init_exon_from_hsp($hsp), $structure );
            if ( !defined $expected ) { next }
            if ( $hsp->{subject_begin} < $expected->{subject_begin} ) {
                if ( $hsp->{orientation} == 1 ) {
                    $hsp->{query_left} += 3 *
                      ( $expected->{subject_begin} - $hsp->{subject_begin} );
                }
                else {
                    $hsp->{query_right} -= 3 *
                      ( $expected->{subject_begin} - $hsp->{subject_begin} );
                }
                #if ( !$trims ) { print "\n" }
                #print "    trimmed 5' end of overextended HSP: $subject_id | $genename | $product from $hsp->{subject_begin} to $expected->{subject_begin}\n";
                $hsp->{subject_begin} = $expected->{subject_begin};
                #$trims++;
            }
            if ( $hsp->{subject_end} > $expected->{subject_end} ) {
                if ( $hsp->{orientation} == 1 ) {
                    $hsp->{query_right} -=
                      3 * ( $hsp->{subject_end} - $expected->{subject_end} );
                }
                else {
                    $hsp->{query_left} +=
                      3 * ( $hsp->{subject_end} - $expected->{subject_end} );
                }
                #if ( !$trims ) { print "\n" }
                #print "    trimmed 3' end of overextended HSP: $subject_id | $genename | $product from $hsp->{subject_end} to $expected->{subject_end}\n";
                $hsp->{subject_end} = $expected->{subject_end};
                #$trims++;
            }
            
            # PROJECT POSITION OF GENE BASED ON hsp
            ( $$hsp{proj_left}, $$hsp{proj_right} ) = project_gene($hsp);

            # group by subject/orientation
            my $genekey = "$genename\t$$hsp{orientation}";
            if ( exists $genes{$genekey} ) {
                push @{ $genes{$genekey} }, $hsp;
            }
            else {
                my @tmp = ($hsp);
                $genes{$genekey} = \@tmp;
            }
        }
        #if ($trims) { print "\n" }
        for my $genekey ( keys %genes ) {
            my $binno = 0;
            my %gbins;
#my ( $gsym ) = split /\t/, $genekey;
#if ( $gsym eq "EBNA-LP" ) { $dbg = 1 }
            for my $hsp ( sort { $$a{proj_left} <=> $$b{proj_left} }
                @{ $genes{$genekey} } )
            {
                my $assigned = 0;
                my $hspwidth = $$hsp{proj_right} - $$hsp{proj_left} + 1;
                for my $binid ( keys %gbins ) {
                    my $binwidth = ( $$hsp{proj_right} - $$hsp{proj_left} + 1 );
                    my $width = $binwidth > $hspwidth ? $binwidth : $hspwidth;
                    my $bin = $bins{$binid};
                    my $overlap = minval( $gbins{$binid}{right}, $$hsp{proj_right} )
                                    - maxval( $gbins{$binid}{left}, $$hsp{proj_left} ) + 1;
                    if ( $overlap >= $width / 3.0 ) {
                        if ( $$hsp{proj_right} > $gbins{$binid}{right} ) {
                            $gbins{$binid}{right} = $$hsp{proj_right};
                        }
                        push @{ $bins{$binid}{ $$hsp{subject_id} } }, $hsp;
                        $$hsp{bin} = $binid;
                        $assigned = 1;
                        if ($dbg) {
                            print "bin $binid  $gbins{$binid}{left}-$gbins{$binid}{right}  overlap=$overlap  ref $$hsp{subject_id}  sbj $$hsp{subject_begin}-$$hsp{subject_end}  qry $$hsp{query_left}-$$hsp{query_right}  proj $$hsp{proj_left}-$$hsp{proj_right}\n";
                        }
                        last;
                    }
                }
                if ( ! $assigned ) {
                    $binno++;
                    my $binid = "$genekey\t$binno";
                    $gbins{$binid}{left}  = $$hsp{proj_left};
                    $gbins{$binid}{right} = $$hsp{proj_right};
                    push @{ $bins{$binid}{ $$hsp{subject_id} } }, $hsp;
                    $$hsp{bin} = $binid;
                    if ($dbg) {
                        print "bin $binid  $gbins{$binid}{left}-$gbins{$binid}{right}  overlap=new  ref $$hsp{subject_id}  sbj $$hsp{subject_begin}-$$hsp{subject_end}  qry $$hsp{query_left}-$$hsp{query_right}  proj $$hsp{proj_left}-$$hsp{proj_right}\n";
                    }
                }
            }

            # adjustment for circular genomes
            if ( $genome->{is_circular} ) {
                for my $b1 ( keys %gbins ) {
                    if ( ! defined $gbins{$b1} ) { next }
                    if ( $gbins{$b1}{right} > $genome->{original_seqlen} ) {
#print "B1: $b1 $gbins{$b1}{left} - $gbins{$b1}{right}";
                        if ( $gbins{$b1}{left} > $genome->{original_seqlen} ) {
                            delete $gbins{$b1};
                            delete $bins{$b1};
#print "  dropped\n";
                            next;
                        }
                        my ( $l, $r ) = ( $gbins{$b1}{left}, $gbins{$b1}{right} );
                        $l -= $genome->{original_seqlen};
                        $r -= $genome->{original_seqlen};
#print "  adjusted: $l to $r\n";
                        my $w1 = $r - $l + 1;
                        for my $b2 ( keys %gbins ) {
                            if ( ! defined $gbins{$b2} ) { next }
                            if ( $b2 eq $b1 ) { next }
                            
                            my $w2 = $gbins{$b2}{right} - $gbins{$b2}{left} + 1;
#print "  B1: $b2 $gbins{$b2}{left} - $gbins{$b2}{right}";
                            my $w = $w1 < $w2 ? $w1 : $w2;
                            my $overlap = minval( $gbins{$b2}{right}, $r )
                                    - maxval( $gbins{$b2}{left}, $l ) + 1;
                            if ( $overlap >= $w / 3.0 ) {
                                delete $gbins{$b2};
                                delete $bins{$b2};
#print "  dropped\n";
                            }
                            else {
#print "  kept\n";
                            }
                        }
                    }
                }
            }
#$dbg = 0;
        }
        
        #
        for my $binid ( keys %bins ) {
            my ( $genename, $orientation, $binno ) = split /\t/, $binid;
            for my $subject_id ( keys %{ $bins{$binid} } ) {
                my $subjkey = "$subject_id\t$orientation\t$binno";
                $subjects{$subjkey} = $bins{$binid}{$subject_id};
            }
        }
    }
    else {
        for my $hsp (@$hsps) {
            my $subjkey = "$$hsp{subject_id}\t$$hsp{orientation}";
            if ( exists $subjects{$subjkey} ) {
                push @{ $subjects{$subjkey} }, $hsp;
            }
            else {
                my @tmp = ($hsp);
                $subjects{$subjkey} = \@tmp;
            }
        }
    }

    for my $subjkey ( keys %subjects ) {
        my ( $subject_id, $orientation ) = split /\t/, $subjkey;

        #if ( $subject_id eq "AP_000539.1" ) { $dbg = 1 } else { $dbg = 0 }
        my @tmp =
          sort { $$a{subject_begin} <=> $$b{subject_begin} }
          remove_conflicting_hsps2( @{ $subjects{$subjkey} } );
        my $hit = generate_hit( \@tmp, $query_ratio );
        push @besthits, $hit;
    }

    if ($dbg) {
        print "BEST HITS\n";
        print_blasthits( 0, @besthits );
    }

    return @besthits;
}

sub get_candidate_descriptor {
    my ( $genome, $hit ) = @_;

    my $ref       = get_reference_seq( $hit->{subject_id} );
    my $structure = get_gene_structure( $hit->{subject_id} );
    my $ori = $$hit{orientation};
    my $fs;
    my $qf;
    my @edges;
    my @tmp  = @{ $$hit{hsps} };
    my $nhsp = @tmp;
    my $margins;
    my $lf;
    my $genome_begin = $ori == 1 ? 1 : $$genome{seqlen};
    my $genome_end = $ori == 1 ? $$genome{seqlen} : 1;

    for my $i ( 0 .. @tmp - 1 ) {
        if ( !defined $qf ) {
            $qf = $tmp[$i]{query_frame};
            $lf = $qf;
            $fs = 0;
        }
        elsif ( $tmp[$i]{query_frame} ne $lf ) {
            $fs++;
            $lf = $tmp[$i]{query_frame};
            $qf = 0;
        }

        my $expected =
          expected_exon( init_exon_from_hsp( $tmp[$i] ), $structure );
        my $x         = $expected->{order};

        my $dbeg = $ori == 1 ? $tmp[$i]{query_left} : $tmp[$i]{query_right};
        my $dend = $ori == 1 ? $tmp[$i]{query_right} : $tmp[$i]{query_left};
        my $begmargin = $tmp[$i]->{subject_begin} - $expected->{subject_begin};
        if ( $begmargin > 1 && $ori * ( $dbeg - $genome_begin ) < 3 ) {
            $begmargin = 1;
        }
        my $endmargin = $expected->{subject_end} - $tmp[$i]->{subject_end};
        if ( $endmargin > 1 && $ori * ( $genome_end - $dend ) < 3 ) {
            $endmargin = 1;
        }
        if ( !defined $margins->{$x}->{beg}
            || abs($begmargin) < $margins->{$x}->{beg} )
        {
            $margins->{$x}->{beg} = $begmargin;
        }
        if ( !defined $margins->{$x}->{end}
            || abs($endmargin) < $margins->{$x}->{end} )
        {
            $margins->{$x}->{end} = $endmargin;
        }
    }
    for my $x ( sort { $a <=> $b } keys %$margins ) {
        my $edge = "$x:$margins->{$x}->{beg},$margins->{$x}->{end}";
        push @edges, $edge;
    }

    my $descriptor =
        $hit->{orientation} . "|"
      . $hit->{query_left} . "-"
      . $hit->{query_right} . "|"
      . $nhsp . "/"
      . $qf . "/"
      . $fs . "|"
      . join( "|", @edges );

    #print "descriptor $descriptor\n";
    #print_blasthits( 0, @{ $hit->{hsps} } );

    return $descriptor;
}

#sub get_blast_aligncoords {
#    my ($blasthit) = @_;
#    my @hsps = @{ $$blasthit{hsps} };
#
#    #my $alignment =  $$blasthit{subject_length};
#    my $alignment =
#      "$$blasthit{orientation}+$$blasthit{subject_begin}-"
#      . ( $$blasthit{subject_length} - $$blasthit{subject_end} );
#    for my $hsp (@hsps) {
#
##$alignment .=  ",$$hsp{query_left}-$$hsp{query_right}|$$hsp{orientation}|$$hsp{subject_begin}-$$hsp{subject_end}";
#        $alignment .= ",$$hsp{query_left}-$$hsp{query_right}";
#    }
#
#    return $alignment;
#}

sub generate_hit {
    my ( $hsps, $query_ratio ) = @_;
    if ( !defined $query_ratio ) {
        $query_ratio = 3;
    }

    my $dbg = DEBUG;

    #if ( $$hsps[0]{subject_id} eq "AP_000117.1" ) { $dbg = 1 }
    if ($dbg) {
        print "HSPS\n";
        print_blasthits( 0, @$hsps );
    }

    #    my $structure = "";
    #    if ( ! $matpep ) {
    #        get_gene_structure( $$hsps[0]{subject_id} );
    #    }
    my @hithsps = sort { $$a{subject_end} <=> $$b{subject_end} } @$hsps;

    my $score = 0;
    my $ori;
    my ( $send, $qend );

#    # remove conflicting HSPs
#    remove_conflicting_hsps2( @hithsps );
##    remove_conflicting_hsps( \@hithsps, "query_left", "query_right", 0.5,
##        "bit_score" );
#    if ($dbg) {
#        print "\nDe-conflicted\n";
#        print_blasthits( 0, @hithsps );
#    }

    #    my ( $smin, $smax );
    #    my ( $qmin, $qmax );
    #    my %exons;
    for my $i ( 0 .. @hithsps - 1 ) {
        if ( $query_ratio != 3 ) {
            $hithsps[$i]{query_frame} = 0;
            $hithsps[$i]{query_stops} = 0;
        }

        # add hsp to score
        my $weight =
          ( $hithsps[$i]{num_identical} + $hithsps[$i]{num_similar} ) / 2.0;
        $hithsps[$i]{weight} = $weight;

        my $similarity =
          ( $hithsps[$i]{pct_identity} + $hithsps[$i]{pct_similarity} ) / 200.0;
        my $scr = $similarity * $weight;

        my $sovr = 0;
        my $qovr = 0;
        if ( !defined $ori ) {
            $ori = $hithsps[$i]{orientation};
            if ( $ori == 1 ) {
                $qend = $hithsps[$i]{query_right};
            }
            else {
                $qend = $hithsps[$i]{query_left};
            }
            $send = $hithsps[$i]{subject_end};

            #            $qmin = $hithsps[$i]{query_left};
            #            $qmax = $hithsps[$i]{query_right};
            #            $smin = $hithsps[$i]{subject_begin};
            #            $smax = $hithsps[$i]{subject_end};
        }
        else {
            if ( $hithsps[$i]{subject_begin} < $send ) {
                $sovr = $send - $hithsps[$i]{subject_begin} + 1;
            }
            $send = $hithsps[$i]{subject_end};

            if ( $ori == 1 ) {
                if ( $hithsps[$i]{query_left} <= $qend ) {
                    $qovr = ( $qend - $hithsps[$i]{query_left} + 1 );
                }
                $qend = $hithsps[$i]{query_right};
            }
            else {
                if ( $hithsps[$i]{query_right} >= $qend ) {
                    $qovr = ( $hithsps[$i]{query_right} - $qend + 1 );
                }
                $qend = $hithsps[$i]{query_left};
            }

#            if ( $qmin > $hithsps[$i]{query_left} ) { $qmin = $hithsps[$i]{query_left} }
#            if ( $qmax < $hithsps[$i]{query_right} ) { $qmax = $hithsps[$i]{query_right} }
#            if ( $smin > $hithsps[$i]{subject_begin} ) { $smin = $hithsps[$i]{subject_begin} }
#            if ( $smax < $hithsps[$i]{subject_end} ) { $smax = $hithsps[$i]{subject_end} }
        }

        $score += ( $scr - $qovr / $query_ratio - $sovr );
        if ($dbg) {
            print "Q: $hithsps[$i]{query_left}-$hithsps[$i]{query_right}"
              . " S: $hithsps[$i]{subject_begin}-$hithsps[$i]{subject_end}"
              . " Scores: w: $weight s: $scr qovr: $qovr  sovr: $sovr  Total: $score\n";
        }

      #        # construct exons for checking gene structure
      #        my $exon = init_exon_from_hsp( $hithsps[$i] );
      #        my $expected = expected_exon( $exon, $structure );
      #        if ( ! defined $expected ) {
      #            $score -= $hithsps[$i]{alignment_length};
      #        }
      #        elsif ( exists $exons{ $$expected{order} } ) {
      #             $exons{ $$expected{order} }{subject_end} = $hithsps[$i]{subject_end};
      #             if ( $ori == 1 ) {
      #                 $exons{ $$expected{order} }{dna_end} = $hithsps[$i]{query_right};
      #             }
      #             else {
      #                 $exons{ $$expected{order} }{dna_end} = $hithsps[$i]{query_left};
      #             }
      #        }
      #        else {
      #            $exons{ $$expected{order} } = $exon;
      #            $exons{ $$expected{order} }{expected} =  $expected;
      #        }
    }

#    # score compliance with expected structure
#    my @exons = sort { $$a{expected}{order} <=> $$b{expected}{order} } values %exons;
#    if ( @exons ) {
#        set_cdna_coordinates( @exons );
#    }
#    if ( $ori == 1 ) {
#        my $offset = $qmin;
#
#    }
#    else {
#        my $offset = $qmax;
#    }

    # build hit
    my $hit = sum_hsps(@hithsps);
    $$hit{vigor_matchwgt} = $score;

    # return hit
    if ($dbg) {
        print "HIT $$hit{vigor_matchwgt}\n";
        print_blasthits( 0, $hit );
        print_blasthits( 0, @hithsps );
    }

    return $hit;
}

# calculate hit statistics from constituent hsps
sub sum_hsps {
    my (@arghsps) = @_;

    my $hit;
    my @hsps;
    push @hsps, sort { $$a{subject_begin} <=> $$b{subject_begin} } @arghsps;

    # initialize hit
    %$hit = %{ $hsps[0] };
    $$hit{hsps} = \@hsps;

    # stretch hit to span all hsps
    for my $i ( 1 .. @hsps - 1 ) {
        if ( ${ $hsps[$i] }{subject_begin} < $$hit{subject_begin} ) {
            $$hit{subject_begin} = $hsps[$i]{subject_begin};
        }
        if ( ${ $hsps[$i] }{subject_end} > $$hit{subject_end} ) {
            $$hit{subject_end} = $hsps[$i]{subject_end};
        }
        if ( ${ $hsps[$i] }{query_left} < $$hit{query_left} ) {
            $$hit{query_left} = $hsps[$i]{query_left};
        }
        if ( ${ $hsps[$i] }{query_right} > $$hit{query_right} ) {
            $$hit{query_right} = $hsps[$i]{query_right};
        }
        if ( ${ $hsps[$i] }{query_frame} != $$hit{query_frame} ) {
            $$hit{query_frame} = 0;
        }
        if ( ${ $hsps[$i] }{evalue} < $$hit{evalue} ) {
            $$hit{evalue} = $hsps[$i]{evalue};
        }
        if ( ${ $hsps[$i] }{bit_score} > $$hit{bit_score} ) {
            $$hit{bit_score} = $hsps[$i]{bit_score};
        }
    }

    # sum hsp alignment counts
    my @attributes = (
        "alignment_length", "num_identical",
        "num_similar",      "subject_gaps",
        "query_gaps",       "query_stops"
    );
    for my $attr (@attributes) {
        $$hit{$attr} = 0;
        for my $hsp (@hsps) {
            $$hit{$attr} += $$hsp{$attr};
        }
    }

    # join hsp alignment strings
    for
      my $alignment ( "query_alignstr", "midline_alignstr", "subject_alignstr" )
    {
        my $alignstr = "";
        for my $hsp ( sort { compare_qry_positions( $a, $b ) } @hsps ) {
            if ( !defined $$hsp{$alignment} ) {
                $alignstr = undef;
                last;
            }
            $alignstr .= $$hsp{$alignment};
        }
        if ( defined $alignstr ) {
            $$hit{$alignment} = $alignstr;
        }
    }

    # calculate coverage for each hsp
    for my $hsp (@hsps) {
        if ( !defined $$hsp{subject_coverage} ) {
            $$hsp{subject_coverage} =
              $$hsp{subject_end} - $$hsp{subject_begin} + 1;
        }
        if ( !defined $$hsp{query_coverage} ) {
            $$hsp{query_coverage} = $$hsp{query_right} - $$hsp{query_left} + 1;
        }
    }

    # calculate hit coverage without double counting overlapping hsps
    my $subject_coverage = 0;
    my $subject_end      = 0;
    for my $hsp ( sort { $$a{subject_begin} <=> $$b{subject_begin} } @hsps ) {
        if ( $$hsp{subject_begin} > $subject_end ) {
            $subject_end = $$hsp{subject_end};
            $subject_coverage += $subject_end - $$hsp{subject_begin} + 1;
        }
        elsif ( $$hsp{subject_end} > $subject_end ) {
            $subject_coverage += $$hsp{subject_end} - $subject_end;
            $subject_end = $$hsp{subject_end};
        }
    }
    $$hit{subject_coverage} = $subject_coverage;

    my $query_coverage = 0;
    my $query_right    = 0;
    for my $hsp ( sort { $$a{query_left} <=> $$b{query_left} } @hsps ) {
        if ( $$hsp{query_left} > $query_right ) {
            $query_right = $$hsp{query_right};
            $query_coverage += $query_right - $$hsp{query_left} + 1;
        }
        elsif ( $$hsp{query_right} > $query_right ) {
            $query_coverage += $$hsp{query_right} - $query_right;
            $query_right = $$hsp{query_right};
        }
    }
    $$hit{query_coverage} = $query_coverage;

    calculate_blast_percentages($hit);

    return $hit;
}

sub default_candidate_comparison {
    my ( $ihit, $jhit ) = @_;

    if ( $$ihit{vigor_matchwgt} > $$jhit{vigor_matchwgt} ) {
        return 1;
    }
    elsif ( $$jhit{vigor_matchwgt} > $$ihit{vigor_matchwgt} ) {
        return -1;
    }
    elsif ( $$ihit{bit_score} > $$jhit{bit_score} ) {
        return 1;
    }
    elsif ( $$jhit{bit_score} > $$ihit{bit_score} ) {
        return -1;
    }
    elsif ( $$ihit{vigor_pctsimilarity} > $$jhit{vigor_pctsimilarity} ) {
        return 1;
    }
    elsif ( $$jhit{vigor_pctsimilarity} > $$ihit{vigor_pctsimilarity} ) {
        return -1;
    }
    elsif ( $$ihit{subject_id} lt $$jhit{subject_id} ) {
        return 1;
    }
    elsif ( $$jhit{subject_id} lt $$ihit{subject_id} ) {
        return -1;
    }
    else {
        return 0;
    }
}

# determine if candidates overlap sufficiently
# to be mututally exclusive
sub candidates_overlap {
    my ( $hit1, $hit2, $overlap_threshold ) = @_;

    if ( !defined $overlap_threshold ) {
        $overlap_threshold = get_parameter("candidate_overlap_threshold");
    }

    # calculate span overlap
    my $dbg = DEBUG;

#if ( $$hit1{subject_id} eq "AAQ10563.1" || $$hit2{subject_id} eq "AAQ10563.1" ) { $dbg = 1 }
#if ( `whoami` =~ /jhoover/ && get_reference_name( $$hit1{subject_id} ) eq "7a" ) { $dbg = 1 }

    my $overlap;

    # calculate in-frame overlap between hsps
    if ($dbg) {
        print "\nOVERLAP1\n";
        print_blasthits( 0, @{ $$hit1{hsps} } );
        print "OVERLAP2\n";
        print_blasthits( 0, @{ $$hit2{hsps} } );
    }
    $overlap = blasthit_hsp_overlap( $hit1, $hit2 );

    # compare overlap to threshold
    my $name1       = get_reference_name( $$hit1{subject_id} );
    my $name2       = get_reference_name( $$hit2{subject_id} );
    my $overlapped1 = $overlap / $$hit1{query_coverage};
    my $overlapped2 = $overlap / $$hit2{query_coverage};
    if ($dbg) {
        print "\nOVERLAP ( threshold $overlap_threshold)\n"
          . "  1. $name1 $overlap / $$hit1{query_coverage} = $overlapped1\n"
          . "  2. $name2 $overlap / $$hit2{query_coverage} = $overlapped2\n";
    }
    if ( $overlap <= 0 ) {
        if ($dbg) { print "HITS DO NOT OVERLAP\n" }
        return 0;
    }

    if ( $overlap / $$hit1{query_coverage} >= $overlap_threshold ) {
        if ($dbg) { print "HITS OVERLAP\n" }
        return 1;
    }
    if ( $overlap / $$hit2{query_coverage} >= $overlap_threshold ) {
        if ($dbg) { print "HITS OVERLAP\n" }
        return 1;
    }
    if ($dbg) { print "HITS DO NOT OVERLAP\n" }
    return 0;
}

sub blasthit_hsp_overlap {
    my ( $hit1, $hit2 ) = @_;

    my $overlap = 0;
    for my $hsp1 ( sort { compare_qry_positions( $a, $b ) } @{ $$hit1{hsps} } )
    {
        for my $hsp2 ( sort { compare_qry_positions( $a, $b ) }
            @{ $$hit2{hsps} } )
        {
            if ( $$hsp2{query_frame} != $$hsp1{query_frame} ) { next }
            my $hsp_overlap =
              minval( $$hsp2{query_right}, $$hsp1{query_right} ) -
              maxval( $$hsp2{query_left}, $$hsp1{query_left} ) + 1;
            if ( $hsp_overlap > 0 ) { $overlap += $hsp_overlap }
        }
    }

    return $overlap;
}

sub remove_undefs {
    my (@in) = @_;
    if ( !@in ) { return () }

    my @out;
    for my $i ( 0 .. @in - 1 ) {
        if ( defined $in[$i] ) { push @out, $in[$i] }
    }

    #print "remove undefs in " . scalar @in . " out " . scalar @out . "\n";
    return @out;
}
########################################
# fill in pieces of the potential gene
sub fill_in_candidate_gaps {
    my ( $genome, $potentialgene ) = @_;

    # loop through possible permutations for this candidate
    # generate all possible small exons for this permutation
    my @permutations;
    for my $candidate ( @{ $$potentialgene{permutations} } ) {
        my @candidates = find_missing_pieces( $genome, $candidate );
        if (@candidates) {
            push @permutations, @candidates;
        }
    }

    # update candidate's list of permutations
    $$potentialgene{permutations} = \@permutations;
}

sub find_missing_pieces {
    my ( $genome, $rawgene ) = @_;
    my $dbg = DEBUG;

#if ( $$genome{id} =~ /^255$/ ) { $dbg = 1 } 
#if ( $$rawgene{subject_id} eq "AP_000539.1" ) { $dbg = 1 }
#if ( get_reference_name( $$rawgene{subject_id} ) =~ /^EBNALP$/ ) { $dbg = 1 }

    # attempt to locate any missing exons
    my %newexons;
    my $num_missing = 0;
    my $complete_genome = $$genome{is_complete};
    my $ori = $$rawgene{orientation};
    my $ref = get_reference_seq( $$rawgene{subject_id} );
    my $tiny5 = get_tiny_exon5( $$rawgene{subject_id} );
    my $tiny3 = get_tiny_exon3( $$rawgene{subject_id} );
    my $structure = get_gene_structure( $$rawgene{subject_id} );
    my $spliced = allow_splicing( $$rawgene{subject_id} );
    my @expected_exons = @{ $$structure{exons} };
    my $vigorspace = get_parameter("vigorspace");
    my %found;
    my @exons;

    for my $hsp ( @{ $$rawgene{hsps} } ) {
        my $rawexon = init_exon_from_hsp($hsp);
        push @exons, $rawexon;

        my $expectedexon = expected_exon( $rawexon, $structure );
        if ( exists $found{ $$expectedexon{order} } ) {
            push @{ $found{ $$expectedexon{order} } }, $hsp;
        }
        else {
            my @tmp = ($hsp);
            $found{ $$expectedexon{order} } = \@tmp;
        }
    }

    if ($dbg) {
        print "check for missing exons\n";
        print_blasthits( 0, @{ $$rawgene{hsps} } );

        for my $exonid ( sort { $a <=> $b } keys %found ) {
            print "exon #$exonid\n";
            print_blasthits( 0, @{ $found{$exonid} } );
        }
    }

    # check for each expected exon
    if ( $spliced ) {
        my $continue = 1;
        while ( $continue ) {
            $continue = 0;
            for my $expected ( 1 .. $$structure{num_exons} ) {
        
                # have we found it?
                if ( exists $found{$expected} ) { next }
        
                # exon is missing, where should we look for it?
                $num_missing++;
        
                if ($dbg) {
                    print_hash( "missing $expected",
                        $expected_exons[ $expected - 1 ] );
                }
        
                # estimate position based on preceding exon
                my @ranges;
                if ( exists $found{ $expected - 1 } ) {
                    my @befores = sort { $$b{subject_end} <=> $$a{subject_end} } @{ $found{ $expected - 1 } };
                    for my $hsp ( @befores ) {
                        my $before = init_exon_from_hsp( $hsp );
                        if ( $dbg ) {
                            print_hash( "preceding exon", $before );
                        }
                        my $expected_sz =  int( 1.1 * $expected_exons[$expected-1]{exon_size} );
                        my $before_proj_end = $$before{dna_end}
                            + $ori * 3 * ( $expected_exons[ $expected - 2 ]{subject_end} - $$before{subject_end} );
        
                        my ( $min_intron, $max_intron ) = get_intron_size( $$rawgene{subject_id}, $expected - 1 );
                            
                        if ($dbg) {
                            print "intron $min_intron-$max_intron  prior exon ends at s: $$before{subject_end}  d: $$before{dna_end}\n";
                        }
            
                        my $range;
                        $$range{begin} = $before_proj_end + $ori * $min_intron - $ori * 2;
                        $$range{end}   = $before_proj_end + $ori * $max_intron + $ori * $expected_sz + $ori * 2;
                        push @ranges, $range;
                        if ($dbg) {
                            print "uprange  $$range{begin} to $$range{end}\n";
                        }
                    }
                }
        
                # estimate position based on following exon
                if ( exists $found{ $expected + 1 } ) {
                    my @afters = sort { $$a{subject_begin} <=> $$b{subject_begin} } @{ $found{ $expected + 1 } };
                    for my $hsp( @afters ) {
                        my $after = init_exon_from_hsp( $hsp );
                        if ( $dbg ) {
                            print_hash( "following exon", $after );
                        }

                        my $expected_sz =  int( 1.1 * $expected_exons[$expected-1]{exon_size} );
                        my $after_proj_begin = $$after{dna_begin}
                            - $ori * 3 * ( $$after{subject_begin} - $expected_exons[$expected]{subject_begin} );
                        my ( $min_intron, $max_intron ) = get_intron_size( $$rawgene{subject_id}, $expected );
                        if ($dbg) {
                            print "intron $min_intron-$max_intron  following exon starts at s: $$after{subject_begin}  d: $$after{dna_begin}  proj: $after_proj_begin\n";
                        }
            
                        my $range;
                        $$range{begin} = $after_proj_begin - $ori * $max_intron - $ori * $expected_sz - $ori * 2;
                        $$range{end} = $after_proj_begin - $ori * $min_intron + $ori * 2;

                        push @ranges, $range;
                        if ($dbg) {
                            print "downrange  $$range{begin} to $$range{end}\n";
                        }
                    }
                }
        
                # combined range estimates
                my $dna_begin;
                my $dna_end;
                for my $range (@ranges) {
                    if ( !defined $dna_begin || $ori * $$range{begin} < $ori * $dna_begin ) {
                        $dna_begin = $$range{begin};
                    }
                    if ( !defined $dna_end || $ori * $$range{end} > $ori * $dna_end ) {
                        $dna_end = $$range{end};
                    }
                    if ($dbg) {
                        print "merging range $$range{begin} to $$range{end}\n";
                    }
                }
                if ( !defined $dna_begin ) { next }

                # search for HSP
                my $expectedexon = $expected_exons[ $expected - 1 ];
                my $expectedsz = 3 * ( $$expectedexon{subject_end} - $$expectedexon{subject_begin} + 1 );
                my $ref_begin = $$expectedexon{subject_begin};
                my $ref_end   = $$expectedexon{subject_end};
                if ( $ref_end > $$rawgene{subject_length} ) {
                    $ref_end = $$rawgene{subject_length};
                }
        
                # orient range in forward direction
                if ( $ori == -1 ) {
                    ( $dna_begin, $dna_end ) = ( $dna_end, $dna_begin );
                }
                if ($dbg) {
                    print "searching for subject $ref_begin-$ref_end in range $dna_begin to $dna_end\n";
                }
        
                # find N-gaps in range
                my $midrange = ( $dna_begin + $dna_end ) / 2.0;
                my $gapstart = $dna_begin < 1 ? 1 : $dna_begin;
                my $gapend = $dna_end > $$genome{seqlen} ? $$genome{seqlen} : $dna_end;
                my $gapmid = ( $dna_begin + $dna_end ) / 2.0;
                my @gaps = find_gaps( $genome, $gapstart, $gapend );
                if ( $dna_begin < 1 ) {
                    if ( ! $complete_genome ) {
                        push @gaps, { begin => $dna_begin, end => 0 };
                    }
                    $dna_begin = 1;
                }
                if ( $dna_end > $$genome{seqlen} ) {
                    if ( ! $complete_genome ) {
                        push @gaps, { begin => $$genome{seqlen} + 1, end => $dna_end };
                    }
                    $dna_end = $$genome{seqlen};
                }
                if ($dbg) {
                    print "#gaps " . @gaps . "  dna range $dna_begin to $dna_end\n";
                }
        
                # could exon be hidden in gap?
                my $ngaps = 0;
                my $gapscore;
                my $bestgap;
                for my $gap (@gaps) {
                    if ( $dbg ) { print "gap: $$gap{begin} to $$gap{end}\n" }
                    my $gapsize = $$gap{end} - $$gap{begin} + 1;
                    $ngaps += $gapsize;
                    my $score = abs( ( $$gap{end} + $$gap{begin} ) / 2.0 - $midrange );
                    if ( $gapsize < $expectedsz ) {
                        $score += ( $expectedsz - $gapsize );
                    }
                    if ( !defined $gapscore || $score < $gapscore ) {
                        $gapscore = $score;
                        $bestgap  = $gap;
                    }
                }
                if ( defined $bestgap && $expectedsz - $ngaps < 10 ) {
                    my $gaphsp = make_gap_hsp(
                        $genome, $midrange, $bestgap, $rawgene,
                        $$expectedexon{subject_begin},
                        $$expectedexon{subject_end}
                    );
                    push @{ $newexons{$expected} }, $gaphsp;
                    push @{ $found{$expected} }, $gaphsp;
                    $continue = 1;
                if ($dbg) {
                        print "potential exon in gap\n";
                        print_blasthits( 0, $gaphsp );
                    }
                }
        
                if ($dbg) {
                    print "exon# $expected  range: $dna_begin-$dna_end  gaps: $ngaps\n";
                }
    
                # search for HSP
                my @search;
                if ( $expected == 1 && $tiny5 =~ /\w/ ) {
                    @search = find_tiny_exon5( $genome, $dna_begin, $dna_end, $rawgene, $tiny5 );
                    if ( $dbg ) {
                        print "exon $expected tiny5 candidates\n";
                        print_blasthits( 0, @search );
                    }
                    if ( @search ) {
                        push @{ $newexons{$expected} }, @search;
                        push @{ $found{$expected} }, @search;
                        $continue = 1;
                        last;
                    }
                }
                elsif ( $expected == $structure->{num_exons} && $tiny3 =~ /\w/ ) {
                    @search = find_tiny_exon3( $genome, $dna_begin, $dna_end, $rawgene, $tiny3 );
                    if ( $dbg ) {
                        print "exon $expected tiny3 candidates\n";
                        print_blasthits( 0, @search );
                    }
                    if ( @search ) {
                        push @{ $newexons{$expected} }, @search;
                        push @{ $found{$expected} }, @search;
                        $continue = 1;
                        last;
                    }
                }
                else {
                    my $rl = $ref_end - $ref_begin + 1;
                    my $minsim = 75;
                    if ( $rl < 4 ) {
                        $minsim = 100;
                    }
                    @search = find_small_hsps(
                        $genome, $dna_begin, $dna_end, $rawgene,
                        0.0005,    5,          10,        $minsim,
                        $ref,    $ref_begin, $ref_end );
                    if ( $dbg ) {
                        print "exon $expected small candidates\n";
                        print_blasthits( 0, @search );
                    }
                    if ( @search ) {
                        push @{ $newexons{$expected} }, @search;
                        push @{ $found{$expected} }, @search;
                        $continue = 1;
                        last;
                    }
                }

            }
        }
    }
    if ($dbg) {
        for my $exonid ( sort { $a <=> $b } keys %found ) {
            print "exon #$exonid\n";
            print_blasthits( 0, @{ $found{$exonid} } );
        }
    }

    # look for missing (frameshifted?) pieces of reference
    # (skip if frameshift sensitivity has been set low)
    my @permutations = ( $rawgene );
    my $fs_sensitivity = get_parameter("frameshift_sensitivity");
    if ( $fs_sensitivity > 0 && ! $spliced ) {
        my $smallest = 12;
        if ( $fs_sensitivity > 1 ) { $smallest = 5 }
        my @fillhsps;

        # set up pieces missing from 5'/3' edges
        my $missing = $$rawgene{subject_begin} - 1;
        if ($dbg) {
            print "missing from start $missing\n";
        }
        if ( exists $found{1} && $missing >= 3 ) {
            my $exon;
            $$exon{is_edge}     = 1;
            $$exon{subject_end} = 0;
            $$exon{dna_end}     =
              $exons[0]{dna_begin} - $ori * 3 * $missing - $ori * 25;
            if ( $$exon{dna_end} < 0 ) {
                $$exon{dna_end} = 0;
            }
            elsif ( $$exon{dna_end} > $$genome{seqlen} + 1 ) {
                $$exon{dna_end} = $$genome{seqlen} + 1;
            }
            unshift @exons, $exon;
            if ($dbg) {
                print_hash( "\nprefix exon", $exon );
            }
        }

        $missing = $$rawgene{subject_length} - $$rawgene{subject_end};
        if ($dbg) {
            print "missing from ($$structure{num_exons}) end $missing\n";
        }
        if ( exists $found{ $$structure{num_exons} } && $missing >= $smallest )
        {
            my $exon;
            $$exon{is_edge}       = 1;
            $$exon{subject_begin} = $$rawgene{subject_length} + 1;
            $$exon{dna_begin}     =
              $exons[ @exons - 1 ]{dna_end} + $ori * 3 * $missing + $ori * 25;
            if ( $$exon{dna_begin} < 0 ) {
                $$exon{dna_begin} = 0;
            }
            elsif ( $$exon{dna_begin} > $$genome{seqlen} + 1 ) {
                $$exon{dna_begin} = $$genome{seqlen} + 1;
            }
            push @exons, $exon;
            if ($dbg) {
                print_hash( "\nsuffix exon", $exon );
            }
        }
        if ($dbg) {
            for my $exon (@exons) {
                print_hash( "\ncheck exons2", $exon );
            }
        }

        # check for missing portion of subject in the dna between the hsps
        if ( @exons > 1 ) {
            while ( @exons > 1 ) {
                my $lastexon = shift @exons;
                
                # ignore frameshifted regions below a certain length
                # (relax standard for missing start codon)
                my $limit = $smallest;
                if ( $$lastexon{subject_end} == 0 ) {
                    $limit = 3;
                }
                
                # don't bother to find pieces missing in-between
                # if frameshift sensitivity has been set low
                elsif ( $fs_sensitivity < 1 ) {
                    next;
                }
                if ($dbg) {
                    print_hash( "\nlastexon", $lastexon );
                    print_hash( "\nexons[0]", $exons[0] );
                }
                
                # establish the dna & ref ranges
                my $dna_begin = $$lastexon{dna_end} + $ori;
                my $dna_end   = $exons[0]{dna_begin} - $ori;
                if ( $ori == -1 ) {
                    ( $dna_begin, $dna_end ) = ( $dna_end, $dna_begin );
                }
                
                my $ref_begin = $$lastexon{subject_end} + 1;
                my $ref_end   = $exons[0]{subject_begin} - 1;
                my $missing   = $ref_end - $ref_begin + 1;
                if ($dbg) {
                    print "missing $missing at ref $ref_begin-$ref_end from dna $dna_begin-$dna_end\n";
                }
                
                # find the "fill" HSPs
                if ( $missing >= $limit )      {
                    push @fillhsps,
                      find_small_hsps(
                        $genome, $dna_begin, $dna_end, $rawgene,
                        1E-5,    1E-3,       10,       60,
                        $ref,    $ref_begin, $ref_end
                      );
                }
            }
        }

        # update the candidate with the "fill" HSPs
        if (@fillhsps) {
            if ($dbg) {
                print "\nfill hsps\n";
                print_blasthits( 0, @fillhsps );
            }
            push @fillhsps, @{ $$rawgene{hsps} };
            @fillhsps = remove_conflicting_hsps2( @fillhsps );
            my $hit = generate_hit( \@fillhsps );
            @permutations = ( $hit );
            if ($dbg) {
                print "\npermutations\n";
                print_blasthits( 0, @permutations );
            }
        }
    }

    # generate permutation for all combinations of expected exons
    for my $exonid ( sort { $a <=> $b } keys %newexons ) {
        my @candidates = @{ $newexons{$exonid} };
        if ( ! @candidates ) { next }
        if ($dbg) {
            print "candidates for exon $exonid: " . @candidates . "\n";
        }
        for my $i ( 0 .. @permutations - 1 ) {
            my $permutation = $permutations[$i];
            if ( !defined $permutation ) { next }
            $permutations[$i] = undef;
            for my $exonhsp (@candidates) {
                my @hsps = @{ $$permutation{hsps} };
                push @hsps, $exonhsp;
                @hsps = remove_conflicting_hsps2( @hsps );
                my $hit = generate_hit( \@hsps );
                push @permutations, $hit;
            }
        }

        @permutations = remove_undefs(@permutations);
    }

    # return results
    if ($dbg) {
        print "final permutations\n";
        print_blasthits( 0, @permutations );
    }
    return @permutations;
}

sub make_gap_hsp {
    my ( $genome, $midrange, $gap, $hit, $sbjleft, $sbjright ) = @_;

    my %new = %$hit;
    delete $new{hsps};
    delete $new{query_alignstr};
    delete $new{subject_alignstr};
    delete $new{midline_alignstr};
    $new{in_gap} = 1;

    $new{subject_begin}    = $sbjleft;
    $new{subject_end}      = $sbjright;
    $new{subject_coverage} = $sbjright - $sbjleft + 1;

    $new{query_left}     = $$gap{begin} + 1;
    $new{query_right}    = $$gap{end} - 1;
    $new{query_coverage} = $new{query_right} - $new{query_left} + 1;

    # adjust edges to preserve codon boundaries
    my $edge = 1;
#print "$new{query_left} to $new{query_right} ($new{query_coverage})\n";
    while ( $new{query_coverage} % 3 > 0 ) {
        if ($edge) {
            $new{query_left}++;
        }
        else {
            $new{query_right}--;
        }
        $new{query_coverage}--;
        $edge = 1 - $edge;
#print "$new{query_left} to $new{query_right} ($new{query_coverage})\n";
    }

    # adjust edges to exon size
#print "$new{query_left} to $new{query_right} $new{query_coverage} vs $new{subject_coverage}\n";
    while ( $new{query_coverage} > 3 * $new{subject_coverage} ) {
        if ($edge) {
            $new{query_left} += 2;
            $new{query_right}--;
        }
        else {
            $new{query_right} -= 2;
            $new{query_left}++;
        }
        $new{query_coverage} -= 3;
        $edge = 1 - $edge;

#print "$new{query_left} to $new{query_right} $new{query_coverage} vs $new{subject_coverage}\n";
    }
    
    # shift as close as possible to expected location
    while ( ( $new{query_left} + $new{query_right} ) / 2.0 < $midrange ) {
        if (  $new{query_right} >= $$gap{end} - 1 ) { last }
        $new{query_left}++;
        $new{query_right}++;
#print "$new{query_left} to $new{query_right} in gap $$gap{begin} to $$gap{end} ($midrange)\n";
    }
    while ( ( $new{query_left} + $new{query_right} ) / 2.0 > $midrange ) {
        if ( $new{query_left} <= $$gap{begin} + 1 ) { last }
        $new{query_left}--;
        $new{query_right}--;
#print "$new{query_left} to $new{query_right} in gap $$gap{begin} to $$gap{end} ($midrange)\n";
    }

    # set metrics
    $new{alignment_length} = $new{query_coverage};
    $new{num_identical}    =
      int( $$hit{pct_identity} / 100.0 * $new{query_coverage} );
    $new{num_similar} =
      int( $$hit{pct_similarity} / 100.0 * $new{query_coverage} );

    my $frame =
      $$hit{orientation} == 1
      ? ( $new{query_left} - 1 ) % 3 + 1
      : -( ( $$genome{seqlen} - $new{query_right} ) % 3 + 1 );
    $new{bit_score}      = 0;
    $new{evalue}         = 0.5;
    $new{has_frameshift} = 0;

    $new{vigor_pctsimilarity} = calculate_vigor_pctsimilarity( \%new ) * 0.9;
    $new{vigor_matchwgt}      = $new{vigor_pctsimilarity} / 100.0 *
      maxval( $new{alignment_length}, $new{subject_length} );

    return \%new;
}

sub find_small_hsps {
    my (
        $genome, $dna_begin, $dna_end, $rawgene,   $mineval, $maxeval,
        $xeval,  $minsim,    $ref,     $ref_begin, $ref_end
      )
      = @_;
    my $dbg = DEBUG;

    my $vigorspace = get_parameter("vigorspace");

    my $missing_aa = subsequence( $$ref{sequence}, $ref_begin, $ref_end );
    my $numaa = $ref_end - $ref_begin + 1;
    my $min_numsim = $numaa < 3 ? $numaa : 3;
    my $min_numid = $numaa < 3 ? $numaa : 2;
    open( AA, ">$vigorspace/aa" );
    print AA ">AA\n$missing_aa\n";
    close AA;
    if ($dbg) { print "missing exon: $missing_aa\n" }

    my $dna_seq = subsequence( $$genome{sequence}, $dna_begin, $dna_end );
    open( DNA, ">$vigorspace/dna" );
    print DNA ">DNA\n$dna_seq\n";
    close DNA;
    if ($dbg) { print "candidate dna: $dna_seq\n" }
    my $cmd = "$myBin/formatdb -i $vigorspace/dna -p F";
    &runCmd($vigorspace, $cmd);

    #print "cmd: $cmd\n" . `$cmd` . "\n";

    my $evalue = $mineval;
    my @hsps = ();
    while ( $evalue <= $maxeval && @hsps < 1 ) {
        $cmd = "$myBin/blastall -p tblastn -i $vigorspace/aa -d $vigorspace/dna -g F -W 1 -C 0 -F F -e $evalue -m 7 -o $vigorspace/xml 2> /dev/null";
        
        if (&runCmd($vigorspace, $cmd, 1)) {
            @hsps = parse_blastxml("$vigorspace/xml");
            $evalue = $evalue * $xeval;
        }
        else {
            warn "Blast failed while searching aminoacid sequence \"$missing_aa\" against nucleotide sequence \"$dna_seq\"\n";
            last;
        }
    }

    my @newhsps;
    if (@hsps) {
        my $dna_offset     = $dna_begin - 1;
        my $subject_offset = $ref_begin - 1;
        for my $hsp (@hsps) {
            if ( $$hsp{orientation} != $$rawgene{orientation} ) {
                next;
            }
            if ( $$hsp{pct_similarity} < $minsim ) { next }
            if ( $$hsp{num_similar} < $min_numsim ) { next }
            if ( $$hsp{num_identical} < $min_numid ) { next }
            $$hsp{query_id}         = $$rawgene{query_id};
            $$hsp{query_definition} = $$rawgene{query_definition};
            $$hsp{query_length}     = $$rawgene{query_length};

            $$hsp{subject_id}         = $$rawgene{subject_id};
            $$hsp{subject_definition} = $$rawgene{query_definition};
            $$hsp{subject_length}     = $$rawgene{subject_length};

            ( $$hsp{query_left}, $$hsp{subject_begin} ) =
              ( $$hsp{subject_begin}, $$hsp{query_left} );
            ( $$hsp{query_right}, $$hsp{subject_end} ) =
              ( $$hsp{subject_end}, $$hsp{query_right} );
            ( $$hsp{subject_coverage}, $$hsp{query_coverage} ) =
              ( $$hsp{query_coverage}, $$hsp{subject_coverage} );

            $$hsp{query_left}    += $dna_offset;
            $$hsp{query_right}   += $dna_offset;
            $$hsp{subject_begin} += $subject_offset;
            $$hsp{subject_end}   += $subject_offset;

            if ( $$hsp{orientation} == 1 ) {
                $$hsp{query_frame} = "" . ( ( $$hsp{query_left} - 1 ) % 3 + 1 );
            }
            else {
                $$hsp{query_frame} =
                  "-" . ( ( $$genome{seqlen} - $$hsp{query_right} ) % 3 + 1 );
            }

            calculate_blast_percentages($hsp);
            push @newhsps, $hsp;
        }
    }

    return @newhsps;
}

# find possible locations for tiny 3' exon
sub find_tiny_exon3 {
    my ( $genome, $left, $right, $hit, $tiny_exon3 ) = @_;

    my $dbg = DEBUG;
    #if ( $$hit{subject_id} =~ /YP_089657/ ) { $dbg = 1 }

    if ($dbg) {
        print "find small exon at 3' end ($tiny_exon3)\n";
        print_blasthits( 0, $hit );
    }

    # include possible splice acceptors to help limit the possibilities
    my $offset = 0;
    my $cdna = $tiny_exon3;
    if ( $cdna =~ /^([A-Z]+):([0-9]+)$/ ) {
        $offset = $2;
        $cdna = $1;
    }
    my @acceptors   = keys %{ get_splice_acceptors() };
    my $splice_exon = "(" . join( "|", @acceptors ) . ")" . $cdna;
    if ( $$hit{orientation} == -1 ) {
        $splice_exon = reverse_complement($splice_exon);
        $splice_exon =~ tr/[]()/][)(/;
    }

    my @hsps;
    for my $exon ( find_regexp( $splice_exon, $$genome{sequence}, $left, $right, 1 ) ) {
        if ($dbg) {
            print "stop $$hit{begin}-$$hit{end}: $$hit{string}\n";
        }
        
        # remove splice acceptor from location
        if (  $$hit{orientation} == 1 ) {
            $$exon{begin} += 2;
            $$exon{end} -= $offset;
        }
        else {
            $$exon{end} -= 2;
            $$exon{begin} += $offset;
        }

        # create and save as HSP
        push @hsps, make_hsp(
            $genome, $hit,
            $$hit{subject_length} + 1, $$hit{subject_length} + 1,
            $$exon{begin}, $$exon{end}, 1, 1
        );
    }

    return @hsps;
}

# find possible locations for tiny 5' exon
sub find_tiny_exon5 {
    my ( $genome, $left, $right, $hit, $tiny_exon5 ) = @_;

    my $dbg = DEBUG;
    #if ( $$hit{subject_id} =~ /YP_089657/ ) { $dbg = 1 }

    if ($dbg) {
        print "find small exon at 5' end ($tiny_exon5)\n";
        print_blasthits( 0, $hit );
    }

    # include possible splice donors to help limit the possibilities
    my $offset = 0;;
    my $cdna = $tiny_exon5;
    if ( $cdna =~ /^([A-Z]+):([0-9]+)$/ ) {
        $offset = $2;
        $cdna = $1;
    }
    my @donors = keys %{ get_splice_donors() };
    my $splice_exon = $cdna . "(" . join( "|", @donors ) . ")";
    if ( $$hit{orientation} == -1 ) {
        $splice_exon = reverse_complement($splice_exon);
        $splice_exon =~ tr/[]()/][)(/;
    }

    my @hsps;
    for my $exon ( find_regexp( $splice_exon, $$genome{sequence}, $left, $right, 1 ) ) {
        if ($dbg) {
            print "stop $$hit{begin}-$$hit{end}: $$hit{string}\n";
        }

        # remove splice donor from location
        if (  $$hit{orientation} == 1 ) {
            $$exon{end} -= 2;
            $$exon{begin} += $offset;
        }
        else {
            $$exon{begin} += 2;
            $$exon{end} -= $offset;
        }

        # create and save as HSP
        push @hsps, make_hsp(
            $genome, $hit,
            1, 1,
            $$exon{begin}, $$exon{end}, 1, 1
        );
    }

    return @hsps;
}


sub find_regexp {
    my ( $regexp, $string, $begin, $end, $allpossibilities ) = @_;
    if ( !defined $allpossibilities ) { $allpossibilities = 0 }

    my $dbg = DEBUG;
    if ( !defined $begin ) { $begin = 1 }
    if ( !defined $end )   { $end   = length($string) }

    if ($dbg) { print "search for $regexp from $begin-$end\n" }

    pos($string) = maxval( $begin - 1, 0 );
    my @matches;
    while ( $string =~ m/$regexp/g ) {
        my %match;
        $match{string} = $&;
        $match{end}    = pos($string);
        if ( $match{end} > $end ) { last }
        $match{begin} = $match{end} - length( $match{string} ) + 1;
        if ($allpossibilities) { pos($string) = $match{begin} }
        if ($dbg) {
            print
"$regexp matches \"$match{string}\" at $match{begin}-$match{end}\n";
        }
        push @matches, \%match;

        #pos($string) = $match{begin};
    }

    return @matches;
}



sub get_refstat_attributes {
    return (
        "pct_refidentity", "pct_refsimilarity", "pct_refsimilarity25",
        "pct_refcoverage", "num_refidentical",  "num_refsimilar",
        "num_refsim25",    "num_refcovered",    "pct_reftrunc5",
        "pct_refgap",      "pct_reftrunc3",     "num_reftrunc5",
        "num_refgap",      "num_reftrunc3"
    );
}

########################################
sub find_gene {
    my ( $genome, $rawgene ) = @_;

    my $dbg = DEBUG;
    #if ( get_reference_name( $$rawgene{subject_id} ) eq "Pp1ab" ) { $dbg = 1 } 

    my $min_gene_exon_size = get_parameter("min_gene_exon_size");

    my $sequence = $$genome{sequence};
    my $seq_len  = length($sequence);
    my $warning;

    # find cds in raw gene
    if ($dbg) { print "find_gene:find_cds\n" }
    my $cds = find_cds( $genome, $rawgene );
    if ($dbg) {
        print_hash( "find_gene CDS", $cds );
    }

    # get protein sequence
    my $protein = $$cds{protein};

    #print "FINAL PROT =$protein\n";
    my $protein_length = length($protein);
    if ( $protein =~ /\*$/ ) {
        $protein_length--;
    }

    # return gene
    my %gene;
    $gene{genome} = $genome;
    $gene{is_pseudogene} = 0;
    if ( is_pseudogene($cds) ) { $gene{is_pseudogene} = 1 }
    $gene{invalidsplicing} = $$cds{invalidsplicing};
    $gene{embeddedstops} = $$cds{embeddedstops};
    $gene{frameshifts} = $$cds{frameshifts};
    $gene{gaperrors} = $$cds{gaperrors};
    $gene{cdsnotes} = $$cds{cdsnotes};
    $gene{cdserrors} = $$cds{cdserrors};

    my $genomeID = $$rawgene{query_id};
    $genomeID =~ s/\|/_/g;
    $gene{gene_id} = "$genomeID.$$rawgene{gene_num}";

    $gene{genome} = $genome;
    $gene{ref_id} = $$rawgene{subject_id};
    $gene{ref_length} = $$rawgene{subject_length};
    my @exons = @{ $$cds{exons} };
    my $first = $exons[0];
    my $last  = $exons[ @exons - 1 ];
    $gene{ref_begin} = $$first{subject_begin};
    $gene{ref_end}   = $$last{subject_end};

    $gene{orientation}           = $$rawgene{orientation};
    $gene{exons}                 = $$cds{exons};
    $gene{noncanonical_splicing} = $$cds{noncanonical_splicing};

#if ( defined $gene{noncanonical_splicing} ) {
#    print "gene: noncanon: " . join( " | ", @{ $gene{noncanonical_splicing} } ) . "\n";
#}

    $gene{cdna}             = $$cds{cdna};
    $gene{start_site}       = $$cds{start_site};
    $gene{codon_start}      = $$cds{codon_start};
    $gene{start_truncation} = $$cds{start_truncation};
    $gene{stop_site}        = $$cds{stop_site};
    $gene{stop_truncation}  = $$cds{stop_truncation};
    if ( exists $$cds{alternate_startcodon} ) {
        $gene{alternate_startcodon} = $$cds{alternate_startcodon};
    }
    if ( exists $$cds{stopcodon_readthru} ) {
        $gene{stopcodon_readthru} = $$cds{stopcodon_readthru};
    }
    if ( exists $$cds{alternate_stopcodon} ) {
        $gene{alternate_stopcodon} = $$cds{alternate_stopcodon};
    }
    $gene{protein}        = $protein;
    $gene{protein_length} = $protein_length;

    $gene{gene_name} = get_reference_name( $gene{ref_id} );
    if ( !defined $gene{gene_name} ) {
        $gene{gene_name} = $gene{gene_id};
    }
    $gene{product_name} = get_reference_product( $gene{ref_id} );
    $gene{note}         = get_reference_note( $gene{ref_id} );

    #print "find: ref=$gene{ref_id}  note=$gene{note}\n";
    $gene{gene_synonym} = get_reference_gene_synonym( $gene{ref_id} );

    #print "find_gene out\n";
    # score gene versus reference
    score_gene( $genome, \%gene );

    # double-check 5' truncation
    if (   $gene{start_truncation}
        && $gene{protein} =~ /^M/i
        && $gene{num_reftrunc5} == 0 )
    {
        $gene{start_truncation} = 0;
        my @exons = @{ $gene{exons} };
        $exons[0]{fuzzy_begin} = 0;
    }

    # return gene
    return \%gene;
}

sub score_gene {
    my ( $genome, $gene ) = @_;

    my $dbg = DEBUG;

    #if ( $$gene{gene_name} eq "L1-3" ) { $dbg = 1 }
    #if ( $$gene{ref_id} eq "YP_002213832.1" ) { $dbg = 1 }
    #if ( allow_splicing( $$gene{ref_id} ) ) { $dbg = 1 }
    #if ( $$gene{gene_name}  eq "E2B" ) { $dbg = 1 }

    # validate gene vs reference
    my $ref    = get_reference_seq( $$gene{ref_id} );
    my $refpep = $$ref{sequence};

    score_pep_vs_reference( $$gene{gene_id}, $$gene{protein}, $$gene{ref_id},
        $refpep, $gene );
    round_refstats($gene);

    $$gene{ref_alignment} = bl2seq_alignment( $genome, $gene );

    $$gene{splicing_quality} = sum_splicing_quality( @{ $$gene{exons} } );

    my $fs_penalty = calculate_frameshift_penalty(
        $$gene{ref_id},    $$gene{start_site},
        $$gene{stop_site}, $$gene{frameshifts}
    );
    my $splicing_penalty =
      calculate_splicing_penalty( $$gene{ref_id}, $$gene{invalidsplicing} );
    my $stop_penalty =
      calculate_stops_penalty( $$gene{ref_id}, $$gene{protein} );
    my $cds_penalty =
      calculate_cds_penalty( $$gene{ref_id}, $$gene{cdserrors} );

    my $partial_penalty = 0;
    my $ref_begin       = 1;
    if ( $$gene{start_truncation} ) {
        $ref_begin       = $$gene{ref_begin};
        $partial_penalty =
          maxval( 5, ( $$gene{ref_end} - $$gene{ref_begin} + 1 ) / 50.0 );
    }
    my $ref_end = $$gene{ref_length};
    if ( $$gene{stop_truncation} ) {
        $ref_end         = $$gene{ref_end};
        $partial_penalty =
          maxval( 5, ( $$gene{ref_end} - $$gene{ref_begin} + 1 ) / 50.0 );
    }
    my $ref_size = $ref_end - $ref_begin + 1 - $$gene{num_refgap};

    my $match         = $$gene{match_quality};
    my $form_mismatch = gene_form_mismatch($gene);
    $match -= $form_mismatch / 2.0;
    my $mismatch = maxval( 0, $$gene{num_refcovered} - $match );

    $$gene{gene_quality} = $match - 1.5 * $mismatch + $$gene{splicing_quality} +
      $$gene{pct_refcoverage} / 100.0 - $partial_penalty;
    $$gene{gene_quality} =
      $$gene{gene_quality} * $fs_penalty * $splicing_penalty * $stop_penalty *
      $cds_penalty;
    if ($dbg) {
        print
"\nPENALTIES quality $$gene{gene_quality}  match/mismatch $match/$mismatch cds $cds_penalty  stop $stop_penalty  fs $fs_penalty  spl $splicing_penalty  partial $partial_penalty REF $$gene{ref_length}\n";
        for my $stat ( get_refstat_attributes() ) {
            print "$stat: $$gene{$stat}\n";
        }
        print "\n";
        print_genehits($gene);
    }
    return;
}

# find longest cds within raw gene
sub find_cds {
    my ( $genome, $rawgene ) = @_;
    my $dbg = DEBUG;
    #if ( get_reference_name( $$rawgene{subject_id} ) eq "Pp1ab" ) { $dbg = 1 }
 
    if ($dbg) {
        print "$$rawgene{subject_id}\n";
        print_blasthits( 0, @{ $$rawgene{hsps} } );
    }
    if ( $$genome{is_circular} && $$rawgene{query_left} > $$genome{original_seqlen} ) {
        $$rawgene{query_left} -= $$genome{original_seqlen};
        $$rawgene{query_right} -= $$genome{original_seqlen};
        for my $hsp (  @{ $$rawgene{hsps} } ) {
            $$hsp{query_left} -= $$genome{original_seqlen};
            $$hsp{query_right} -= $$genome{original_seqlen};
        }
    }

    # get reference information
    my $ref = get_reference_seq( $$rawgene{subject_id} );
    my $stop_readthru  = get_readthru_exception( $$rawgene{subject_id} );
    my $allow_splicing = allow_splicing( $$rawgene{subject_id} );
    my $allow_ribosomal_slippage = allow_ribosomal_slippage( $$rawgene{subject_id} );
    my $allow_rna_editing = allow_rna_editing( $$rawgene{subject_id} );
    my $is_optional    = is_optional( $$rawgene{subject_id} );
    my $is_required    = is_required( $$rawgene{subject_id} );

    # convert hsps into exons
    my @exons;
    my $ori = $$rawgene{orientation};
    my @missingexons;

    # adjust hsps according to splicing rules
    if ($allow_splicing) {
        if ( $dbg ) { print "extract_spliced_exons\n" }
        @exons = extract_spliced_exons( $genome, $rawgene );

        # diagnose missing exons
        my $structure = get_gene_structure( $$rawgene{subject_id} );
        my %found;
        my @missing;
        for my $exon (@exons) {
            my $expected = expected_exon( $exon, $structure );
            $found{ $$expected{order} } = $exon;
        }
        for my $i ( 1 .. $$structure{num_exons} ) {
            if ( !exists $found{$i} ) { push @missingexons, $i }
        }
    }

    # adjust hsps according to slippage rules
    elsif ( $allow_ribosomal_slippage ) {
        if ( $dbg ) { print "extract_ribosomalslippage_exons in\n" }
        @exons = extract_ribosomalslippage_exons( $genome, $rawgene );
        if ( $dbg ) { print "extract_ribosomalslippage_exons out\n" }
    }

    # adjust hsps according to rna editing rules
    elsif ( $allow_rna_editing ) {
        if ( $dbg ) { print "extract_rnaediting_exons\n" }
        @exons = extract_rnaediting_exons( $genome, $rawgene );
    }

    # adjust exons according to default rules
    else {
        if ( $dbg ) { print "extract_unspliced_exons\n" }
        @exons = extract_unspliced_exons( $genome, $rawgene );
    }

    # find possible starts
    if ( $dbg ) { print "find_possible_starts\n" }
    my @starts = find_possible_starts( $genome, $rawgene, @exons );

    # find possible stops
    if ( $dbg ) { print "find_possible_stops\n" }
    my @stops = find_possible_stops( $genome, $rawgene, @exons );

    # score all possible ORFs
    if ( $dbg ) { print "set_cdna_coordinates\n" }
    set_cdna_coordinates(@exons);
    if ( $dbg ) { print "cdna_from_exons\n" }
    my $cdna      = cdna_from_exons( $genome, @exons );
    my $cdnalen   = length($cdna);
    my $ref_begin = maxval( 1, $exons[0]{subject_begin} );
    my $ref_end   = minval( $$ref{seqlen}, $exons[ @exons - 1 ]{subject_end} );
    my $ref_structure = get_gene_structure( $$ref{id} );
    my $refseq = subsequence( $$ref{sequence}, $ref_begin, $ref_end );
    if ($dbg) {
        print "REF=$refseq\n";
    }
    my %ref_profile = profile_peptide($refseq);

    my $bestorf;
    for my $start (@starts) {
        for my $stop (@stops) {
            if ( $ori * $$stop{dnapos} < $ori * $$start{dnapos} ) { next }
            my ( $cdna_begin, $cdna_end ) =
              dna_coords_to_cdna( $$start{dnapos}, $$stop{dnapos}, \@exons );
            my $orf = subsequence( $cdna, $cdna_begin, $cdna_end );
            my $aa = DNA2AA( $orf, 1, $$start{translation_exception},
                undef, $$stop{translation_exception} );
            my $aascore = score_profile( $aa, %ref_profile );

            my $pseudo_penalties = 0;
            my $lastorder        = 0;
            my $lastframe;
            for my $exon (@exons) {
                if ( $ori * $$exon{dna_begin} > $ori * $$start{dnapos} ) {
                    next;
                }
                if ( $ori * $$exon{dna_end} < $ori * $$stop{dnapos} ) { next }
                my $expected = expected_exon( $exon, $ref_structure );
                if ( !defined $expected ) {
                    $pseudo_penalties += 3;
                    if ($dbg) {
                        print "pseudo penalty: unexpected exon\n";
                    }
                }
                elsif ($$expected{order} eq $lastorder
                    && $$exon{frame} != $lastframe )
                {
                    $pseudo_penalties++;
                    if ($dbg) {
                        print "pseudo penalty: frameshift\n";
                    }
                }
                $lastorder = $$expected{order};
                $lastframe = $$exon{frame};
            }
            my @internalorfs = split /\*/, substr( $aa, 0, length($aa) - 1 );
            if ( @internalorfs > 1
                && allow_stopcodon_readthru( $$rawgene{subject_id} ) )
            {
                pop @internalorfs;
            }
            if ( @internalorfs > 1 ) {
                $pseudo_penalties += ( @internalorfs - 1 );
                if ($dbg) {
                    print "pseudo penalty: internalorfs ("
                      . ( @internalorfs - 1 ) . ")\n";
                }
            }
            if ( $$start{error} ) {
                $pseudo_penalties++;
                if ($dbg) {
                    print "pseudo penalty: bad start\n";
                }
            }
            if ( $$stop{error} ) {
                $pseudo_penalties++;
                if ($dbg) {
                    print "pseudo penalty: bad stop\n";
                }
            }

            $pseudo_penalties = sqrt($pseudo_penalties) * $$ref{seqlen} / 20.0;

            my $score =
              $aascore - $$start{penalty} - $$stop{penalty} - $pseudo_penalties;

            if ($dbg) {
                print
"orf dna $$start{dnapos}..$$stop{dnapos}  cdna $cdna_begin-$cdna_end  score $score  aascore $aascore  penalties start $$start{penalty}  stop $$stop{penalty}  pseudo $pseudo_penalties\n";
                print "AA=$aa\n";
                print "orf=$orf\n";
            }

            if ( !defined $$bestorf{score} || $score > $$bestorf{score} ) {
                $$bestorf{score}   = $score;
                $$bestorf{aascore} = $aascore;
                $$bestorf{start}   = $start;
                $$bestorf{stop}    = $stop;
                $$bestorf{cdna}    = $orf;
                $$bestorf{aa}      = $aa;
            }
        }
    }

    if ($dbg) {
        print
"bestorf $$bestorf{start}{dnapos}..$$bestorf{stop}{dnapos}  score $$bestorf{score}  aascore $$bestorf{aascore}\n";
        print_hash( "\nbestorf", $bestorf );
    }

    $$bestorf{start}{codon_start} = 1;
    if (   $$bestorf{start}{start_truncation}
        && $$bestorf{start}{gap_truncation} )
    {
        while ( !defined $$genome{in_gap}{ $$bestorf{start}{dnapos} - $ori } ) {
            $$bestorf{start}{dnapos} -= $ori;
            $$bestorf{start}{codon_start} =
              $$bestorf{start}{codon_start} % 3 + 1;
        }
        $exons[0]{dna_begin} = $$bestorf{start}{dnapos};
    }

    # format and return best CDS
    my $cds;

    # trim exons to best orf and calculate cdna
    while (@exons
        && $ori * $exons[0]{dna_end} < $ori * $$bestorf{start}{dnapos} )
    {
        shift @exons;
    }
    $exons[0]{dna_begin} = $$bestorf{start}{dnapos};
    $exons[0]{subject_begin} = maxval( 1, $$bestorf{start}{refpos} );
    if ( $$bestorf{start}{start_truncation} ) {
        $exons[0]{fuzzy_begin} = 1;
    }
    else {
        $exons[0]{fuzzy_begin} = 0;
    }

    while (@exons
        && $ori * $exons[ @exons - 1 ]{dna_begin} >
        $ori * $$bestorf{stop}{dnapos} )
    {
        pop @exons;
    }
    $exons[ @exons - 1 ]{dna_end}     = $$bestorf{stop}{dnapos};
    $exons[ @exons - 1 ]{subject_end} =
      minval( $$rawgene{subject_length} + 1, $$bestorf{stop}{refpos} );
    if ( $$bestorf{stop}{stop_truncation} ) {
        $exons[ @exons - 1 ]{fuzzy_end} = 1;
    }
    else {
        $exons[ @exons - 1 ]{fuzzy_end} = 0;
    }

    set_cdna_coordinates(@exons);
    $$cds{exons} = \@exons;
    $$cds{cdna}  = cdna_from_exons( $genome, @{ $$cds{exons} } );

    # start codon
    $$cds{start_site}       = $$bestorf{start}{dnapos};
    $$cds{codon_start}      = $$bestorf{start}{codon_start};
    $$cds{start_truncation} = $$bestorf{start}{start_truncation};
    if ( defined $$bestorf{start}{translation_exception} ) {
        $$cds{alternate_startcodon} = $$bestorf{start}{translation_exception};
        $$cds{alternate_startcodon}{aa_position} = int( length($cdna) / 3 );
    }

    # stop codon and embedded stops
    $$cds{stop_site}       = $$bestorf{stop}{dnapos};
    $$cds{stop_truncation} = $$bestorf{stop}{stop_truncation};
    if ( defined $$bestorf{stop}{translation_exception} ) {
        $$cds{alternate_stopcodon} = $$bestorf{stop}{translation_exception};
    }
    {
        my @internalstops;
        my @orfs = split /(\*)/, $$bestorf{aa};
        if ( $$bestorf{aa} =~ /\*$/ ) { pop @orfs }
        if ( @orfs > 1 ) {
            my $stop = 0;
            for my $orf (@orfs) {
                if ( $orf =~ /\*/ ) {
                    push @internalstops, $stop + 1;
                }
                $stop += length($orf);
            }

            # stop codon read thru?
            if (@internalstops) {
                if ( allow_stopcodon_readthru( $$rawgene{subject_id} ) ) {
                    my $readthru = pop @internalstops;
                    my $translation_exception;
                    $$translation_exception{aa} =
                      get_readthru_exception( $$rawgene{subject_id} );
                    $$translation_exception{aa_position} = $readthru;
                    (
                        $$translation_exception{dna_begin},
                        $$translation_exception{dna_end}
                      )
                      = pep_coords_to_dna( $readthru, $readthru,
                        $$cds{codon_start}, $$cds{exons} );
                    $$cds{stopcodon_readthru} = $translation_exception;
                    $$cds{protein}            = DNA2AA(
                        $$cds{cdna},
                        $$cds{codon_start},
                        $$cds{alternate_startcodon},
                        $$cds{stopcodon_readthru},
                        $$cds{alternate_stopcodon}
                    );
                    $$cds{protein_length} = length( $$cds{protein} );
                    if ( $$cds{protein} =~ /\*$/ ) { $$cds{protein_length}-- }
                }
            }

            # cds interrupted by stops
            if (@internalstops) {
                my @embeddedstops;
                for my $stop (@internalstops) {
                    my ( $begin, $end ) =
                      pep_coords_to_dna( $stop, $stop, $$cds{codon_start},
                        $$cds{exons} );
                    push @embeddedstops, "$begin..$end";
                }
                $$cds{embeddedstops} = \@embeddedstops;
            }
        }
    }

    # protein
    $$cds{protein} = DNA2AA(
        $$cds{cdna}, $$cds{codon_start},
        $$cds{alternate_startcodon},
        $$cds{stopcodon_readthru},
        $$cds{alternate_stopcodon}
    );
    $$cds{protein_length} = length( $$cds{protein} );
    if ( $$cds{protein} =~ /\*$/ ) { $$cds{protein_length}-- }

    # cds notes
    {
        my @notes;
        if ( defined $$bestorf{start}{note} ) {
            push @notes, $$bestorf{start}{note};
        }
        if ( defined $$bestorf{stop}{note} ) {
            push @notes, $$bestorf{stop}{note};
        }
        if (@notes) {
            $$cds{cdsnotes} = \@notes;
        }
    }

    # cds errors
    {
        my @cdserrors;
        if ( $$bestorf{start}{error} ) {
            push @cdserrors, "could not find valid start codon";
        }
        if ( $$bestorf{stop}{error} ) {
            push @cdserrors, "could not find valid stop codon";
        }
        if ( @missingexons ) {
            my $message = "missing exons " . join( ", ", @missingexons );
            if ( @missingexons == 1 ) { $message =~ s/exons/exon/ }
            else { $message =~ s/,([^,]+$)/ and$1/ }
            push @cdserrors, $message;
        }
        if (@cdserrors) {
            $$cds{cdserrors} = \@cdserrors;
        }
    }

    # splicing exceptions
    if ($allow_splicing) {
        my @noncanon;
        my @badsplice;
        my $num_inferred = 0;
        for my $i ( 1 .. @exons - 1 ) {
            if ( $exons[$i]{spliced} > 2 ) { $num_inferred++ }
            if ( $exons[$i]{spliced} == 2 || $exons[$i]{spliced} == 4 ) {
                my @intron = ( $exons[ $i - 1 ]{dna_end} + $ori, $exons[$i]{dna_begin} - $ori );
                if ( $$genome{is_circular} ) {
                    if ( $intron[0] > $$genome{original_seqlen} ) {
                        $intron[0] -= $$genome{original_seqlen};
                    }
                    if ( $intron[1] > $$genome{original_seqlen} ) {
                        $intron[1] -= $$genome{original_seqlen};
                    }
                }
                push @noncanon, join( "..", @intron );
            }
            elsif ( $exons[$i]{spliced} == -1 ) {
                my @intron = ( $exons[ $i - 1 ]{dna_end} + $ori, $exons[$i]{dna_begin} - $ori );
                if ( $$genome{is_circular} ) {
                    if ( $intron[0] > $$genome{original_seqlen} ) {
                        $intron[0] -= $$genome{original_seqlen};
                    }
                    if ( $intron[1] > $$genome{original_seqlen} ) {
                        $intron[1] -= $$genome{original_seqlen};
                    }
                }
                push @badsplice, join( "..", @intron );
            }
        }
        if (@noncanon)  { $$cds{noncanonical_splicing} = \@noncanon }
        if (@badsplice) { $$cds{invalidsplicing}       = \@badsplice }
        if ( $num_inferred > 0 ) {
            my $message =
              "splice sites within gaps inferred from reference sequence";
            if ( $num_inferred == 1 ) {
                $message =~ s/sites/site/;
                $message =~ s/gaps/gap/;
            }
            if ( defined $$cds{cdsnotes} && @{ $$cds{cdsnotes} } ) {
                push @{ $$cds{cdsnotes} }, $message;
            }
            else {
                my @tmp = ($message);
                $$cds{cdsnotes} = \@tmp;
            }
        }
    }

    # gap and frameshift errors
    my @gaperrors;
    my @frameshifts;
    my $has_slippage = 0;
    for my $i ( 1 .. @exons - 1 ) {
        if ( defined $exons[$i]{spliced} && $exons[$i]{spliced} != 0 ) { next }
        if ( defined $exons[$i]{slippage} ) {
            $has_slippage = 1;
            if ( $exons[$i]{slippage} == 2 ) {
                my $message =
                  "position of slippage inferred from reference sequence";
                if ( defined $$cds{cdsnotes} && @{ $$cds{cdsnotes} } ) {
                    push @{ $$cds{cdsnotes} }, $message;
                }
                else {
                    my @tmp = ($message);
                    $$cds{cdsnotes} = \@tmp;
                }
            }
            next;
        }
        if ( $exons[ $i - 1 ]{frame} eq $exons[$i]{frame} ) { next }
        if ( $exons[ $i - 1 ]{in_gap} ) {
            my $gaperror;
            $$gaperror{gap_dna_begin} = $exons[ $i - 1 ]{gap_dna_begin};
            $$gaperror{gap_dna_end}   = $exons[ $i - 1 ]{gap_dna_end};
            $$gaperror{gapsize}       = $exons[ $i - 1 ]{gapsize};
            $$gaperror{expectedsize}  = $exons[ $i - 1 ]{expectedsize};
            push @gaperrors, $gaperror;
        }
        else {
            my $frameshift;
            $$frameshift{dna_begin}     = $exons[ $i - 1 ]{dna_end};
            $$frameshift{dna_end}       = $exons[$i]{dna_begin};
            $$frameshift{orientation}   = $exons[$i]{orientation};
            $$frameshift{subject_begin} = $exons[ $i - 1 ]{subject_end};
            $$frameshift{subject_end}   = $exons[$i]{subject_begin};
            $$frameshift{frame_begin}   = $exons[ $i - 1 ]{frame};
            $$frameshift{frame_end}     = $exons[$i]{frame};
            push @frameshifts, $frameshift;
        }
    }
    if ( $allow_ribosomal_slippage && ! $has_slippage ) {
        my $structure = get_gene_structure( $$rawgene{subject_id} );
        my $expected = ${ $$structure{exons} }[0];
        my $slippage_pos = $$expected{subject_end} + 0.5;
        my ( undef, $slippage_shift ) = describe_ribosomal_slippage( $$rawgene{subject_id} );
        my $bestf;
        my $bestdistance = 30;
        for my $f ( 0..@frameshifts-1 ) {
            my $fs = $frameshifts[$f];
            my $shift = $ori * ( $$fs{frame_end} - $$fs{frame_begin} );
            if ( $shift == $slippage_shift ) {
                my $distance = abs( ( $$fs{subject_begin} + $$fs{subject_end} ) / 2.0 - $slippage_pos );
                if ( $distance < $bestdistance ) {
                    $bestdistance = $distance;
                    $bestf = $f;
                }
            }
        }
        if ( defined $bestf ) {
            $frameshifts[$bestf] = undef;
            @frameshifts = remove_undefs( @frameshifts );
            push @{ $$cds{cdserrors} }, "slippage motif not present";
        }
    }
    if ( @gaperrors && !get_parameter("jcvi_rules") ) {
        $$cds{gaperrors} = \@gaperrors;
    }
    if (@frameshifts) {
        $$cds{frameshifts} = \@frameshifts;
    }

    # return CDS
    if ($dbg) {
        print_hash( "CDS", $cds );
    }
    return $cds;
}

# find possible starts
sub find_possible_starts {
    my ( $genome, $rawgene, @exons ) = @_;
    my $dbg = DEBUG;

    #if ( get_reference_name( $$rawgene{subject_id} ) =~ /^(IVa2|pTP|E2B)$/ ) { $dbg = 1 }

    my $jcvi_rules = get_parameter("jcvi_rules");
    my @starts;

    if ($dbg) {
        print
"--- find possible starts ----------------------------------------------------------\n";
        print_blasthits( 0, @{ $$rawgene{hsps} } );
    }

    # genome boundaries
    my $genome_begin = 1;
    my $genome_end   = $$genome{seqlen};
    my $ori          = $$rawgene{orientation};
    if ( $ori == -1 ) {
        ( $genome_begin, $genome_end ) = ( $genome_end, $genome_begin );
    }

    if ($dbg) {
        print "\ngenome $genome_begin..$genome_end\n";
        print "  gaps:";
        if ( defined $$genome{gaps} ) {
            for my $gap ( @{ $$genome{gaps} } ) {
                print "  $$gap{begin}..$$gap{end}";
            }
        }
        print "\n";
    }

    # reference details
    my $start_codons   = get_startcodons( $$rawgene{subject_id} );
    my $structure      = get_gene_structure( $$rawgene{subject_id} );
    my $gene_variation = get_gene_variation( $$rawgene{subject_id} );
    my $ref            = get_reference_seq( $$rawgene{subject_id} );

    if ($dbg) {
        print_hash( "structure", $structure );
        print "\nexpected\n";
        for my $exon ( @{ $$structure{exons} } ) {
            print_hash( "expected$$exon{order}", $exon );
        }
        print "\nfound\n";
        my $i = 0;
        for my $exon (@exons) {
            $i++;
            print_hash( "exon$i", $exon );
        }
    }

    # initial cdna
    my $cdna          = cdna_from_exons( $genome, @exons );
    my $cdnalen       = length($cdna);
    my $dna_begin     = $exons[0]{dna_begin};
    my $subject_begin = $exons[0]{subject_begin};

# if the initial exon is missing or in a gap, the gene is truncated on the 5' end
    my $first = expected_exon( $exons[0], $structure );
    if ( $$first{order} > 1 || defined $exons[0]{in_gap} ) {
        while ( defined $exons[0]{in_gap} ) { shift @exons }
        my $start;
        $$start{start_truncation} = 1;
        $$start{gap_truncation}   = 0;
        $$start{refpos}           = $exons[0]{subject_begin};
        $$start{dnapos}           = $exons[0]{dna_begin};
        $$start{penalty}          = 1;
        $$start{error}            = 0;
        @starts                   = ($start);

        if ($dbg) {
            print "\n";
            print_hash( "missingexon start", $start );
        }

        return @starts;
    }

    # do we have perfect match to start of reference?
    elsif ( $subject_begin == 1
        && is_start_codon( substr( $cdna, 0, 3 ), $start_codons ) == 1 )
    {
        my $start;
        $$start{start_truncation} = 0;
        $$start{refpos}           = $subject_begin;
        $$start{dnapos}           = $dna_begin;
        $$start{penalty}          = 0;
        $$start{error}            = 0;
        @starts                   = ($start);

        # alternate start codon?
        if ( codon2aa( substr( $cdna, 0, 3 ) ) ne "M" ) {
            my $translation_exception;
            $$translation_exception{aa}          = "M";
            $$translation_exception{aa_position} = 1;
            (
                $$translation_exception{dna_begin},
                $$translation_exception{dna_end}
              )
              = cdna_coords_to_dna( 1, 3, \@exons );
            $$start{translation_exception} = $translation_exception;
        }

        if ($dbg) {
            print "\n";
            print_hash( "perfect start", $start );
        }
    }

    # search for potential starts
    else {

        # calculate expected and allowed range for searching
        my $stretch = $cdnalen /
          ( $exons[ @exons - 1 ]{subject_end} - $exons[0]{subject_begin} + 1.0 )
          / 3.0;
        my $expansion = maxval( $stretch, 1.05 );
        my $shrinkage = minval( $stretch, 0.95 );

        my $lo = int( ( $subject_begin - 1 ) * $expansion + 0.5 );
        my $hi = int( ( $subject_begin - 1 ) * $shrinkage - 0.5 );
        my $expected_subject_lo = $subject_begin - $lo;
        my $min_subject         = $expected_subject_lo - 20;
        my $expected_subject_hi = $subject_begin - $hi;
        my $max_subject         = $expected_subject_hi + 6;
        if ( $gene_variation == 0 ) {
            my $min_subject = maxval( -1, $min_subject );
            my $max_subject = minval( 2,  $max_subject );
        }

        if ($dbg) {
            print
"ori $ori  genome begin $genome_begin  dna begin $dna_begin  subject begin $subject_begin\n";
            print
"stretch: $stretch  range: $min_subject-$expected_subject_lo-$expected_subject_hi-$max_subject\n";
        }

        # default start (from raw alignment)
        {
            my $start;
            $$start{refpos}  = $subject_begin;
            $$start{dnapos}  = $dna_begin;
            $$start{penalty} = 5;
            $$start{error}   = 1;
            @starts          = ($start);

            if ($dbg) {
                print "\n";
                print_hash( "default start", $start );
            }
        }

        # extend cdna
        my $gaplimit = 3;
        if ( !$jcvi_rules ) { $gaplimit = $expected_subject_lo }
        while ($ori * ( $dna_begin - 3 * $ori ) >= $ori * $genome_begin
            && $subject_begin - 1 >= $min_subject )
        {

            if ( $subject_begin <= 1 ) {
                if (
                    is_start_codon(
                        subsequence(
                            $$genome{sequence}, $dna_begin,
                            $dna_begin + 2 * $ori
                        ),
                        $start_codons
                    ) == 1
                  )
                {
                    $min_subject = $subject_begin;
                    last;
                }
            }

            # possible cdna truncation
            if ( defined $$genome{in_gap}{ $dna_begin - 3 * $ori } ) {
                if ( $subject_begin >= $expected_subject_lo ) {
                    my $start;
                    $$start{start_truncation} = 1;
                    $$start{gap_truncation}   = 1;
                    $$start{refpos}           = $subject_begin;
                    $$start{dnapos}           = $dna_begin;
                    if ( $subject_begin < 1 ) {
                        $$start{penalty} = 1.0 + 1.5 * ( 1 - $subject_begin );
                    }
                    else {
                        $$start{penalty} = 1;
                    }
                    $$start{error} = 0;
                    push @starts, $start;

                    if ($dbg) {
                        print "\n";
                        print_hash( "gap start", $start );
                    }
                }

                my $pos = $dna_begin - 3 * $ori;
                if ($dbg) { print "skipping gap: $pos ($subject_begin)" }
                while ( defined $$genome{in_gap}{ $dna_begin - 3 * $ori }
                    && $subject_begin >= $min_subject )
                {
                    $dna_begin -= 3 * $ori;
                    $subject_begin--;
                }
                if ($dbg) { print " - $dna_begin ($subject_begin)\n" }
                if ( $subject_begin < $gaplimit ) {
                    $subject_begin = $exons[0]{subject_begin};
                    $min_subject   = $subject_begin;
                    $dna_begin     = $exons[0]{dna_begin};
                    last;
                }
            }

            # add codon to cdna
            else {
                $dna_begin -= 3 * $ori;
                $subject_begin--;
                $exons[0]{dna_begin}     = $dna_begin;
                $exons[0]{subject_begin} = $subject_begin;

                if ($dbg) {
                    print "extend dna: $dna_begin  subj: $subject_begin\n";
                }
            }
        }
        if (   $ori * ( $dna_begin - 3 * $ori ) < $ori * $genome_begin
            && $subject_begin - 1 >= $min_subject && ! $$genome{is_complete} )
        {
            my $start;
            $$start{start_truncation} = 1;
            $$start{gap_truncation}   = 1;
            $$start{refpos}           = $exons[0]{subject_begin};
            $$start{dnapos}           = $exons[0]{dna_begin};
            if ( $exons[0]{subject_begin} < 1 ) {
                $$start{penalty} = 1 + 1.5 * ( 1 - $exons[0]{subject_begin} );
            }
            else {
                $$start{penalty} = 1;
            }
            $$start{error} = 0;
            push @starts, $start;

            if ($dbg) {
                print "\n";
                print_hash( "genome start", $start );
            }
        }

        # recompute cdna
        set_cdna_coordinates(@exons);
        $cdna          = cdna_from_exons( $genome, @exons );
        $cdnalen       = length($cdna);
        $subject_begin = $exons[0]{subject_begin};
        my $cdnapos = 1;

        if ($dbg) {
            print
"cdnalen: $cdnalen  subject begin: $subject_begin  min: $min_subject  max: $max_subject\n";
        }

        # find potential start codons
        while ( $subject_begin <= $max_subject ) {
            my $codon = subsequence( $cdna, $cdnapos, $cdnapos + 2 );
            my $startcheck = is_start_codon( $codon, $start_codons );

            if ($dbg) {
                print
"cdnapos $cdnapos  subjpos $subject_begin  codon: $codon  check $startcheck\n";
            }

            # known start codon
            if ( $startcheck == 1 ) {
                my $start;
                $$start{start_truncation} = 0;
                $$start{refpos}           = $subject_begin;
                my ( $beg, $end ) =
                  cdna_coords_to_dna( $cdnapos, $cdnapos + 2, \@exons );
                $$start{dnapos} = $beg;
                if ( $subject_begin > $expected_subject_hi ) {
                    $$start{penalty} =
                      2.0 * ( $subject_begin - $expected_subject_hi );
                }
                elsif ( $subject_begin < $expected_subject_lo ) {
                    $$start{penalty} =
                      1.5 * ( $expected_subject_lo - $subject_begin );
                }
                else {
                    $$start{penalty} = 0;
                }
                $$start{error} = 0;
                push @starts, $start;

                # alternate start codon
                if ( codon2aa($codon) ne "M" ) {
                    my $translation_exception;
                    $$translation_exception{aa}          = "M";
                    $$translation_exception{aa_position} = 1;
                    $$translation_exception{dna_begin}   = $beg;
                    $$translation_exception{dna_end}     = $end;
                    $$start{translation_exception} = $translation_exception;
                }

                if ($dbg) {
                    print "\n";
                    print_hash( "codon start", $start );
                }
            }

            # possible start codon (ambiguity codes)
            elsif ($startcheck == 2
                && $codon !~ /N.*N/i
                && !$jcvi_rules
                && $subject_begin >= 0
                && $subject_begin <= 2 )
            {

                my ( $beg, $end ) =
                  cdna_coords_to_dna( $cdnapos, $cdnapos + 2, \@exons );
                my $start;
                $$start{start_truncation} = 0;
                $$start{refpos}           = $subject_begin;
                $$start{dnapos}           = $beg;
                $$start{codon_start}      = 1;
                if ( $subject_begin > $expected_subject_hi ) {
                    $$start{penalty} =
                      1 + 2.0 * ( $subject_begin - $expected_subject_hi );
                }
                elsif ( $subject_begin < $expected_subject_lo ) {
                    $$start{penalty} =
                      1 + 1.5 * ( $expected_subject_lo - $subject_begin );
                }
                else {
                    $$start{penalty} = 1 + ( $subject_begin - 1 ) / 10.0;
                }
                $$start{penalty} += 1;
                $$start{error} = 0;
                $$start{note}  =
                  "probable start codon inferred from reference sequence";
                my $translation_exception;
                $$translation_exception{aa}          = "M";
                $$translation_exception{aa_position} = 1;
                $$translation_exception{dna_begin}   = $beg;
                $$translation_exception{dna_end}     = $end;
                $$start{translation_exception}       = $translation_exception;
                push @starts, $start;

                if ($dbg) {
                    print "\n";
                    print_hash( "inferred start", $start );
                }
            }
            $cdnapos += 3;
            $subject_begin++;
        }
    }
    if ($dbg) {
        for my $start (@starts) {
            print_hash( "\nstart", $start );
        }
        print
"--- end find_possible_starts ----------------------------------------------\n";
    }

    my $refaa = uc( substr( $$ref{sequence}, 0, 30 ) );
    my $aa = DNA2AA($cdna);
    for my $start (@starts) {
        if (   $$start{refpos} > 2
            && $$start{error} == 0
            && $$start{start_truncation} == 0 )
        {
            my ($cdnapos) =
              dna_coords_to_cdna( $$start{dnapos}, $$start{dnapos}, \@exons );
            my $aapos    = ( $cdnapos - 1 ) / 3;
            my $qryaa    = substr( $aa, $aapos, 12 );
            my $indexpos = index( $refaa, substr( $qryaa, 0, 3 ) );
            if ($dbg) {
                print
"start ref: $$start{refpos}  dna: $$start{dnapos}  cdna: $cdnapos  aa: $aapos  index: $indexpos\n";
                print "qryaa: $qryaa  refaa: $refaa\n";
            }
            if ( $indexpos > 0 && abs( $indexpos + 1 - $$start{refpos} ) < 2 ) {
                $$start{penalty} = 5;
                $$start{error}   = 1;
                if ($dbg) {
                    print "\n";
                    print_hash( "internal start", $start );
                }
            }
        }
    }
    return @starts;
}

# find possible stops
sub find_possible_stops {
    my ( $genome, $rawgene, @exons ) = @_;
    my $dbg = DEBUG;

    #if ( get_reference_name( $$rawgene{subject_id} ) =~ /^(IVa2|pTP|E2B)$/ ) { $dbg = 1 }
    if ($dbg) {
        print
"--- find possible stops ----------------------------------------------------------\n";
        for my $exon (@exons) {
            print_hash( "\nexon", $exon );
        }
        print "\n";
    }

    my $jcvi_rules = get_parameter("jcvi_rules");
    my @stops;

    # genome boundaries
    my $genome_begin = 1;
    my $genome_end   = $$genome{seqlen};
    my $ori          = $$rawgene{orientation};
    if ( $ori == -1 ) {
        ( $genome_begin, $genome_end ) = ( $genome_end, $genome_begin );
    }

    # reference details
    my $structure = get_gene_structure( $$rawgene{subject_id} );

    # initial cdna
    my $cdna           = cdna_from_exons( $genome, @exons );
    my $cdnalen        = length($cdna);
    my $dna_end        = $exons[ @exons - 1 ]{dna_end};
    my $subject_end    = $exons[ @exons - 1 ]{subject_end};
    my $subject_length = $$rawgene{subject_length};

    # calculate expected and allowed range for searching
    my $stretch = $cdnalen /
      ( $exons[ @exons - 1 ]{subject_end} - $exons[0]{subject_begin} + 1.0 ) /
      3.0;
    my $expansion = maxval( $stretch, 1.05 );
    my $shrinkage = minval( $stretch, 0.95 );

    my $expected_subject_lo =
      minval( $subject_length - 5, int( $subject_length * 0.95 ) );
    my $expected_subject_hi =
      maxval( $subject_length + 10, int( $subject_length * 1.05 + 0.5 ) );
    my $min_subject = minval( $subject_length - 5,
        int( $shrinkage * 0.95 * $subject_length - 0.5 ) );
    my $max_subject = maxval( $subject_length + 10,
        int( $expansion * 1.05 * $subject_length + 0.5 ) );

    if ($dbg) {
        print
"stretch: $stretch  sublen $subject_length  range: $min_subject-$expected_subject_lo-$expected_subject_hi-$max_subject\n";
    }

    # if the last exon is missing or in a gap, the gene is truncated on the 3' end
    my $partial3 = is_reference_partial3( $$rawgene{subject_id} );
    if ( !$partial3 ) {
        my $last = expected_exon( $exons[ @exons - 1 ], $structure );
        if ( defined $last && $$last{order} < $$structure{num_exons} ) {
            $partial3 = 1;
        }
        elsif ( $exons[ @exons - 1 ]{in_gap} ) {
            $partial3 = 1;
        }
    }
    if ($partial3) {
        while ( @exons > 1 && defined $exons[@exons-1]{in_gap} ) { pop @exons }
        my $stop;
        $$stop{stop_truncation} = 1;
        $$stop{refpos}          = $exons[@exons-1]{subject_end};
        $$stop{dnapos}          = $exons[@exons-1]{dna_end};
        $$stop{penalty}         = 1;
        $$stop{error}           = 0;
        push @stops, $stop;

        if ($dbg) {
            print "\n";
            print_hash( "missingexon stop", $stop );
        }
        return @stops;
    }

    # default stop (from raw alignment)
    {
        my $stop;
        $$stop{refpos}  = $subject_end;
        $$stop{dnapos}  = $dna_end;
        $$stop{penalty} = 5 + 1.5 * ( $subject_length - $subject_end );
        $$stop{error}   = 1;
        @stops          = ($stop);

        if ($dbg) {
            print "\n";
            print_hash( "default stop", $stop );
        }
    }

    # extend cdna
    my $gaplimit = $subject_length - 2;
    if ( !$jcvi_rules ) { $gaplimit = $expected_subject_hi }
    while ($ori * ( $dna_end + 3 * $ori ) <= $ori * $genome_end
        && $subject_end + 1 <= $max_subject )
    {

        if ($dbg) {
            print
"subject_length $subject_length  max_subject $max_subject  subject_end $subject_end  genome_end $genome_end  dna_end $dna_end\n";
        }
        if ( $subject_end > $subject_length - 5 ) {
            if (
                codon2aa(
                    subsequence(
                        $$genome{sequence}, $dna_end - 2 * $ori, $dna_end
                    )
                ) eq "*"
              )
            {
                $max_subject = $subject_end;
                last;
            }
        }

        # possible cdna truncation
        if ( defined $$genome{in_gap}{ $dna_end + 3 * $ori } ) {
            if ( $subject_end <= $expected_subject_hi ) {
                my $stop;
                $$stop{stop_truncation} = 1;
                $$stop{refpos}          = $subject_end;
                $$stop{dnapos}          = $dna_end;
                if ( $subject_end > $expected_subject_hi ) {
                    $$stop{penalty} = ( $subject_end - $expected_subject_hi );
                }
                else {
                    $$stop{penalty} = 0.5;
                }
                $$stop{error} = 0;
                push @stops, $stop;

                if ($dbg) {
                    print "\n";
                    print_hash( "gap stop", $stop );
                }
            }

            my $pos = $dna_end + 3 * $ori;
            if ($dbg) { print "skipping gap: $pos ($subject_end)" }
            while ( defined $$genome{in_gap}{ $dna_end + 3 * $ori }
                && $subject_end <= $max_subject )
            {
                $dna_end += 3 * $ori;
                $subject_end++;
            }
            if ($dbg) { print " - $dna_end ($subject_end)\n" }
            if ( $subject_end > $gaplimit ) {
                $subject_end = $exons[ @exons - 1 ]{subject_end};
                $max_subject = $subject_end;
                $dna_end     = $exons[ @exons - 1 ]{dna_end};
                last;
            }
        }

        # add codon to cdna
        else {
            $dna_end += 3 * $ori;
            $subject_end++;
            $exons[ @exons - 1 ]{dna_end}     = $dna_end;
            $exons[ @exons - 1 ]{subject_end} = $subject_end;
            if ($dbg) {
                print "extend to DNA $dna_end  Subject $subject_end\n";
            }
        }
    }
    if (   $ori * ( $dna_end + 3 * $ori ) > $ori * $genome_end
        && $subject_end + 1 <= $max_subject && ! $$genome{is_complete} )
    {
        my $stop;
        $$stop{stop_truncation} = 1;
        $$stop{refpos}          = $exons[ @exons - 1 ]{subject_end};
        $$stop{dnapos}          = $exons[ @exons - 1 ]{dna_end};

        if ( $exons[ @exons - 1 ]{subject_end} > $expected_subject_hi ) {
            $$stop{penalty} = ( $subject_end - $expected_subject_hi );
        }
        else {
            $$stop{penalty} = 0.5;
        }

        $$stop{error} = 0;
        push @stops, $stop;

        if ($dbg) {
            print "\n";
            print_hash( "genome stop", $stop );
        }
    }

    # recompute cdna
    set_cdna_coordinates(@exons);
    $cdna        = cdna_from_exons( $genome, @exons );
    $cdnalen     = length($cdna);
    $subject_end = $exons[ @exons - 1 ]{subject_end};
    if ($dbg) {
        print "cdnalen = $cdnalen  subject_end=$subject_end\n";
        for my $exon (@exons) {
            print_hash( "\nextended exon", $exon );
        }
        print "\n";
    }

    my $cdnapos = $cdnalen;
    if ($dbg) {
        print "stop cdna=$cdna\n";
    }

    # find potential stop codons
    while ( $subject_end >= $min_subject + 1 ) {
        my $codon = subsequence( $cdna, $cdnapos - 2, $cdnapos );
        my $stopcheck = is_stop_codon($codon);
        if ($dbg) {
            print "cdnapos $cdnapos  subject $subject_end  codon $codon\n";
        }

        # known stop codon
        if ( $stopcheck == 1 ) {
            my $stop;
            $$stop{stop_truncation} = 0;
            $$stop{refpos}          = $subject_end;
            my ( $beg, $end ) =
              cdna_coords_to_dna( $cdnapos - 2, $cdnapos, \@exons );
            $$stop{dnapos} = $end;
            if ( $subject_end > $expected_subject_hi ) {
                $$stop{penalty} = ( $subject_end - $expected_subject_hi );
            }
            elsif ( $subject_end < $expected_subject_lo ) {
                $$stop{penalty} = 1.5 * ( $expected_subject_lo - $subject_end );
            }
            else {
                $$stop{penalty} = 0;
            }
            $$stop{error} = 0;
            push @stops, $stop;

            if ($dbg) {
                print "\n";
                print_hash( "codon stop cdnapos=$cdnapos", $stop );
            }
        }

        # possible stop codon (ambiguity codes)
        elsif ($stopcheck == 2
            && $codon !~ /N.*N/
            && !$jcvi_rules
            && $subject_end >= $subject_length - 1
            && $subject_end <= $subject_length + 1 )
        {

            my ( $beg, $end ) =
              cdna_coords_to_dna( $cdnapos - 2, $cdnapos, \@exons );
            my $stop;
            $$stop{stop_truncation} = 0;
            $$stop{refpos}          = $subject_end;
            $$stop{dnapos}          = $end;
            if ( $subject_end > $expected_subject_hi ) {
                $$stop{penalty} = ( $subject_end - $expected_subject_hi );
            }
            elsif ( $subject_end < $expected_subject_lo ) {
                $$stop{penalty} = 1.5 * ( $expected_subject_lo - $subject_end );
            }
            else {
                $$stop{penalty} = 0;
            }
            $$stop{penalty} += 1;
            $$stop{error} = 0;
            $$stop{note}  =
              "probable stop codon inferred from reference sequence";

            my $translation_exception;
            $$translation_exception{aa}          = "*";
            $$translation_exception{aa_position} = undef;
            $$translation_exception{dna_begin}   = $beg;
            $$translation_exception{dna_end}     = $end;
            $$stop{translation_exception}        = $translation_exception;

            push @stops, $stop;

            if ($dbg) {
                print "\n";
                print_hash( "inferred stop", $stop );
            }
        }
        $cdnapos -= 3;
        $subject_end--;
    }
    if ($dbg) {
        for my $stop (@stops) {
            print_hash( "\nstop", $stop );
        }
        print
"--- end find_possible_stops ----------------------------------------------\n";
    }

    return @stops;
}

sub resolve_ambiguous_splicing {
    my ($selected) = @_;
    my $dbg = DEBUG;

    if ( !defined $selected ) { return }
    if ( !@$selected )        { return }

    my @scores;
    for my $selected (@$selected) {
        my $location   = get_gene_location($selected);
        my $spliceform = get_splice_form( $$selected{ref_id} );
        if ($dbg) {
            print "location: $location\n";
            print "spliceform: $spliceform\n";
        }
        if ( $spliceform gt "" && $location !~ /[><]/ ) {
            my @reftmp = split /[EeIi]/, $spliceform;
            shift @reftmp;
            my @seltmp;
            my $prev;
            for my $exon ( @{ $$selected{exons} } ) {
                if ( defined $prev ) {
                    push @seltmp,
                      abs( $$exon{dna_begin} - $$prev{dna_end} ) - 1;
                }
                push @seltmp, $$exon{cdna_end} - $$exon{cdna_begin} + 1;
                $prev = $exon;
            }
            my $ns    = @seltmp;
            my $nr    = @reftmp;
            my $n     = minval( $ns, $nr );
            my $error = 0;
            for my $i ( 0 .. $n - 1 ) {
                $error += abs( $seltmp[$i] - $reftmp[$i] );
            }
            for my $i ( $n .. $ns - 1 ) {
                $error += $seltmp[$i];
            }
            for my $i ( $n .. $nr - 1 ) {
                $error += $reftmp[$i];
            }
            if ($dbg) {
                print "ref: " . join( "-", @reftmp ) . "\n";
                print "vig: " . join( "-", @seltmp ) . "\n";
                print "err: $error\n";
            }
            my $score;
            $$score{error} = $error;
            $$score{gene}  = $selected;
            push @scores, $score;
        }
    }
    if (@scores) {
        @scores = sort { $$a{error} <=> $$b{error} } @scores;
        my $gene  = $scores[0]{gene};
        my $least = $scores[0]{error};
        shift @scores;
        while ( @scores && $scores[0]{error} == $least ) {
            my $score    = shift @scores;
            my $location = get_gene_location( $$score{gene} );
            $$gene{equivalent_splicing}{$location} = 1;
        }
        @$selected = ($gene);
    }
    else {
        my $gene = shift @$selected;
        while (@$selected) {
            my $tmpgene  = shift @$selected;
            my $location = get_gene_location($tmpgene);
            $$gene{equivalent_splicing}{$location} = 1;
        }
        @$selected = ($gene);
    }

    if ($dbg) {
        print "\nSELECTED SPLCING\n";
        print_genehits(@$selected);
    }
}

sub calculate_frameshift_penalty {
    my ( $ref_id, $dna_begin, $dna_end, $frameshifts ) = @_;

    my $left  = minval( $dna_begin, $dna_end );
    my $right = maxval( $dna_begin, $dna_end );

    my $penalty = 1;
    if ( defined $frameshifts ) {
        my $variation_penalty = get_variation_penalty($ref_id);
        for my $fs (@$frameshifts) {
            my ( $fsleft, $fsright ) = ( $$fs{dna_begin}, $$fs{dna_end} );
            if ( $fsleft > $fsright ) {
                ( $fsleft, $fsright ) = ( $$fs{dna_end}, $$fs{dna_begin} );
            }
            if ( $fsleft >= $left && $fsright <= $right ) {
                $penalty = $variation_penalty * $penalty;
            }
            $variation_penalty = sqrt($variation_penalty);
            $variation_penalty = sqrt($variation_penalty);
        }
    }

    return $penalty;
}

sub calculate_splicing_penalty {
    my ( $ref_id, $badsplices ) = @_;

    my $penalty = 1;
    if ( defined $badsplices ) {
        my $variation_penalty = get_variation_penalty($ref_id);
        for my $spl (@$badsplices) {
            $penalty           = $penalty * $variation_penalty;
            $variation_penalty = sqrt( sqrt($variation_penalty) );
        }
    }

    return $penalty;
}

sub calculate_stops_penalty {
    my ( $ref_id, $aa ) = @_;

    my $penalty = 1;

    my $tmpaa = $aa;
    $tmpaa =~ s/\*$//;
    my $full_length = length($tmpaa);
    $tmpaa =~ s/\*//g;
    my $stopcount = $full_length - length($tmpaa);

    if ( $stopcount == 0 ) { return $penalty }

    my $variation_penalty = get_variation_penalty($ref_id);
    for my $i ( 1 .. $stopcount ) {
        $penalty           = $variation_penalty * $penalty;
        $variation_penalty = sqrt($variation_penalty);
    }

    return $penalty;
}

sub calculate_cds_penalty {
    my ( $ref_id, $cdserrors ) = @_;

    my $penalty = 1;

    if ( !defined $cdserrors ) { return $penalty }
    if ( !@$cdserrors )        { return $penalty }

    my $variation_penalty = get_variation_penalty($ref_id);
    for my $error (@$cdserrors) {
        $penalty           = $variation_penalty * $penalty;
        $variation_penalty = sqrt($variation_penalty);
    }
    return $penalty;
}

sub sum_splicing_quality {
    my (@exons) = @_;

    my $quality = 0;
    for my $exon (@exons) {
        if ( exists $$exon{splice_quality} ) {
            $quality += $$exon{splice_quality};
        }
    }

    return $quality;
}

# derive exons from hsps by starting with the hsp boundaries and
# "polishing" the boundaries to splice sites while maintaining the proper framing
sub extract_spliced_exons {
    my ( $genome, $rawgene ) = @_;
    my $dbg = DEBUG;

#if ( $$rawgene{subject_id} eq "gi_391886505-NS2" ) { $dbg = 1 }
#if ( `whoami` =~ /jhoover/ && get_reference_name( $$rawgene{subject_id} ) =~ /EBNALP/ ) { $dbg = 1 }

    if ($dbg) {
        print "\nRAWGENE HSPS\n";
        print_blasthits( 0, @{ $$rawgene{hsps} } );
    }

    my $structure = get_gene_structure( $$rawgene{subject_id} );
    my $splicesite_extend_range = get_parameter("splicesite_extend_range");

    my $ori      = $$rawgene{orientation};
    my $sequence = $$genome{sequence};
    my $seqlen   = length($sequence);
    my ( $genome_begin, $genome_end ) = ( 1, $seqlen );
    if ( $ori == -1 ) {
        ( $genome_begin, $genome_end ) = ( $seqlen, 1 );
    }

    my $ref              = get_reference_seq( $$rawgene{subject_id} );
    my $gene             = get_reference_gene( $$rawgene{subject_id} );
    my $gene_splicepairs = get_gene_splicepairs($gene);
    my $tiny3 = get_tiny_exon3( $$rawgene{subject_id} );
    my $tiny5 = get_tiny_exon5( $$rawgene{subject_id} );
    
    # reformat hsps as exons
    my @hsps = sort { compare_qry_positions( $a, $b ) } @{ $$rawgene{hsps} };
    my @exons;
    if ($dbg) { print "\nRAWEXONS\n" }
    for my $hsp (@hsps) {
        my $exon = init_exon_from_hsp($hsp);
        if (@exons) {
            my $expected     = expected_exon( $exon, $structure );
            my $prev         = pop @exons;
            my $prevexpected = expected_exon( $prev, $structure );
            if ( $$prevexpected{order} == $$expected{order} ) {
                my @fused =
                  fuse_exons( $genome, $ref, $prev, $exon );
                push @exons, @fused;
            }
            else {
                push @exons, $prev, $exon;
            }
        }
        else {
            push @exons, $exon;
        }
        if ($dbg) {
            print
"  F: $$exon{frame}  Q: $$exon{dna_begin}-$$exon{dna_end}  S: $$exon{subject_begin}-$$exon{subject_end}\n";
        }
    }
    

    # try to adjust exon boundaries to splice sites
    $exons[0]{spliced}          = 0;
    $exons[0]{splicing_quality} = 0;
    my @spliced = ( shift @exons );
    while (@exons) {
        my $exon5 = pop @spliced;
        my $exon3 = shift @exons;
        $$exon3{spliced} = -1;
        $$exon3{splicing_quality} = -5;
        

        my $expected5 = expected_exon( $exon5, $structure );
        my $expected3 = expected_exon( $exon3, $structure );
        my $expected_intron_size = $$expected3{dna_begin} - $$expected5{dna_end} + 1;
        my $intron_size = $ori * ( $$exon3{dna_begin} - $$exon5{dna_end} ) + 1;

        if ($dbg) {
            print "\nSPLICING\n";
            print_hash( "exon5",  $exon5 );
            print_hash( "exon3",  $exon3 );
        }

        # fragments of same exon?
        if (   $$expected5{order} == $$expected3{order} )
        {
            $$exon3{spliced}          = 0;
            $$exon3{splicing_quality} = 0;
            push @spliced, $exon5, $exon3;
            if ( $dbg ) {
                print "continuation of same exon\n";
            }
        }

        # missing exon in between?
        elsif ( $$expected5{order} + 1 != $$expected3{order} ) {
            if ($dbg) {
                print "missing exon\n";
                print_hash( "MISSING X5", $expected5 );
                print_hash( "MISSING X3", $expected3 );
            }

            $$exon3{spliced}          = 0;
            $$exon3{splicing_quality} = -5;
            push @spliced, $exon5, $exon3;
        }

        # exon in gap?
        elsif ( defined $exon5->{in_gap} ) {
            if ($dbg) {
                print "gap exon\n";
                print_hash( "gap X5", $expected5 );
                print_hash( "gap X3", $expected3 );
            }
            $$exon3{spliced}          = 1;
            $$exon3{splicing_quality} = 0;
            push @spliced, $exon5, $exon3;
        }

        # find splice sites

        # establish expected protein by which we will evaluate potential splicing
        else {
            my $refseq = subsequence( $$ref{sequence}, $$rawgene{subject_begin},
                $$exon3{subject_end} );
            my $reflen     = length($refseq);
            my %refprofile = profile_peptide($refseq);
            if ($dbg) {
                print "reference AA: $refseq\n";
            }

            # establish potential range for intron
            my ( $min_intron_size, $max_intron_size ) = get_intron_size( $$ref{id}, $$expected5{order} );
            if ($dbg) {
                print "expected $expected_intron_size  limits $min_intron_size-$max_intron_size\n";
            }

            my $begin5  = $$exon5{dna_end} - 15 * $ori;
            my $sbegin5 = $$exon5{subject_end} - 5;
            while ( $sbegin5 > maxval( 1, $$exon3{subject_begin} - 5 ) ) {
                $begin5 -= 3 * $ori;
                $sbegin5--;
            }
            if ( $begin5 < 1 ) { $begin5 = 1 }
            elsif ( $begin5 > $$genome{seqlen} ) { $begin5 = $$genome{seqlen} }
            if ( $ori * $begin5 < $ori * $$exon5{dna_begin} ) {
                $begin5 = $$exon5{dna_begin};
            }

            my $end3  = $$exon3{dna_begin} + 15 * $ori;
            my $send3 = $$exon3{subject_begin} + 5;
            while (
                $send3 < minval( $$ref{seqlen} + 1, $$exon5{subject_end} + 5 ) )
            {
                $end3 += 3 * $ori;
                $send3++;
            }
            if ( $end3 < 1 ) { $end3 = 1 }
            elsif ( $end3 > $$genome{seqlen} ) { $end3 = $$genome{seqlen} }
            if ( $ori * $end3 > $ori * $$exon3{dna_end} ) {
                $end3 = $$exon3{dna_end};
            }

            # try different splicings
            my $bestscore;
            my @bestexons;
            my $proj5 = $$exon5{dna_end} + $ori * 3 *
              ( $$expected5{subject_end} - $$exon5{subject_end} );
            my $proj3 = $$exon3{dna_begin} - $ori * 3 *
              ( $$exon3{subject_begin} - $$expected3{subject_begin} );
            my $tolerance = 25;
            if ($dbg) {
                print_hash( "exon5",     $exon5 );
                print_hash( "expected5", $expected5 );
                print_hash( "exon3",     $exon3 );
                print_hash( "expected3", $expected3 );
                print "expected intron $expected_intron_size\n";
                print
"gap projection: $$exon5{dna_begin}..$proj5 , $proj3..$$exon3{dna_end}\n";
            }

            my $d5_start = $begin5;
            my $d5_end   = $end3 - $min_intron_size;
            if ( $ori == -1 ) {
                $d5_start = $end3 + $min_intron_size;
                $d5_end   = $begin5;
            }
            if ($dbg) {
                print
"X5 $$exon5{dna_end} | $$exon5{subject_end}  X3 $$exon3{dna_begin} | $$exon3{subject_begin}  "
                  . "Structure $$structure{spliceform}  Ori $ori  Intron $min_intron_size-$max_intron_size  "
                  . "range $begin5..$end3  d5 $d5_start..$d5_end\n";
            }
            if ( $d5_start < $proj5 - $tolerance ) {
                $d5_start = $proj5 - $tolerance;
                if ($dbg) {
                    print "adjusting d5 range to $d5_start-$d5_end\n";
                }
            }
            if ( $d5_end > $proj5 + $tolerance ) {
                $d5_end = $proj5 + $tolerance;
                if ($dbg) {
                    print "adjusting d5 range to $d5_start-$d5_end\n";
                }
            }
            if ( $$expected5{order} ==1 && $tiny5 =~ /\w/ ) {
                $d5_start = $$exon5{dna_end};
                $d5_end = $d5_start;
                if ($dbg) {
                    print "tiny_exon5 adjusting d5 range to $d5_start-$d5_end\n";
                }
            }

            my $d5 = $d5_start - 1;
            while ( $d5 < $d5_end ) {
                $d5++;

                if ( in_gap( $genome, $d5 + $ori, $d5 + 2 * $ori )
                    && abs( $d5 + $ori - $proj5 ) > 2 )
                {
                    if ($dbg) {
                        print "gap: d5=$d5 proj5=$proj5\n";
                    }
                    next;
                }

                my $donor = subsequence( $$genome{sequence}, $d5 + $ori, $d5 + 2 * $ori );
                if ( $dbg && $donor =~ /N/i ) {
                    print
" *** DONOR $donor is_gap? d5=$d5 $$genome{in_gap}{$d5-1}|$$genome{in_gap}{$d5}|$$genome{in_gap}{$d5+1}\n";
                    print "gaps @ "
                      . join( ",",
                        sort { $a <=> $b } keys %{ $$genome{in_gap} } )
                      . "\n";
                }

                if ( validate_splicing( $gene_splicepairs, $donor ) == 0 ) {
                    if ($dbg) {
                        print "$d5: $donor not a valid donor\n";
                    }
                    next;
                }

                my %tmp5 = %$exon5;
                adjust_exon_end( \%tmp5, $d5, $ref );

                my $d3_start = $d5 + $ori * $min_intron_size + $ori;
                my $d3_end   = $d5 + $ori * $max_intron_size + $ori;
                if ( $ori * $d3_end > $ori * $end3 ) { $d3_end = $end3 }
                if ( $ori == -1 ) {
                    ( $d3_start, $d3_end ) = ( $d3_end, $d3_start );
                }
                if ($dbg) {
                    print
                      "D5 $d5  Range $begin5..$end3  D3 $d3_start..$d3_end\n";
                }
                if ( $d3_start < $proj3 - $tolerance ) {
                    $d3_start = $proj3 - $tolerance;
                    if ($dbg) {
                        print "adjusting d3 range to $d3_start-$d3_end\n";
                    }
                }
                if ( $d3_end > $proj3 + $tolerance ) {
                    $d3_end = $proj3 + $tolerance;
                    if ($dbg) {
                        print "adjusting d3 range to $d3_start-$d3_end\n";
                    }
                }
                if ( $$expected3{order} == $$structure{num_exons} && $tiny3 =~ /\w/ ) {
                    $d3_start = $$exon3{dna_begin};
                    $d3_end = $d3_start;
                    if ($dbg) {
                        print "tiny_exon3 adjusting d3 range to $d3_start-$d3_end\n";
                    }
                }
                my $d3 = $d3_start - 1;
                while ( $d3 < $d3_end ) {
                    $d3++;
                    if ($dbg) {
                        print "\nintron between $d5-$d3\n";
                    }

                    # are we in a gap
                    if ( in_gap( $genome, $d3 - 2 * $ori, $d3 - $ori )
                        && abs( $d3 - $ori - $proj3 ) > 2 )
                    {
                        next;
                    }

                    # do we have a valid donor/acceptor pair?
                    my $acceptor = subsequence( $$genome{sequence}, $d3 - 2 * $ori, $d3 - $ori );
                    my ( $splicing, $splice_quality ) = validate_splicing( $gene_splicepairs, $donor, $acceptor );
                    if ($dbg) {
                        print "$d5-$d3  $donor+$acceptor = $splicing | $splice_quality";
                        if ( defined $bestscore ) {
                            print "  $d5-$d3 bestscore $bestscore\n";
                        }
                        else {
                            print "\n";
                        }
                    }
                    if ( $splicing == 0 ) { next }

                    # generate protein resulting from this splicing
                    my %tmp3 = %$exon3;
                    adjust_exon_begin( \%tmp3, $d3 );

                    $tmp3{spliced}        = $splicing;
                    $tmp3{splice_quality} = $splice_quality;

                    my @workingexons = @spliced;
                    push @workingexons, \%tmp5;
                    push @workingexons, \%tmp3;
                    my $splicing_quality = sum_splicing_quality(@workingexons);

                    set_cdna_coordinates(@workingexons);
                    my $cdna = cdna_from_exons( $genome, @workingexons );

                    # adjust expected edges for deletions/insertions
                    my $adjust5 =
                      int( $tmp5{cdna_end} / 3.0 ) - $tmp5{subject_end};
                    my $expedge5 = $$expected5{cdna_end} + $adjust5;
                    my $adjust3  =
                      int( ( $tmp3{cdna_begin} - 1 ) / 3.0 ) + 1 -
                      $tmp3{subject_begin};
                    my $expedge3 = $$expected3{cdna_begin} + $adjust3;

                    # calculate edge mismatches (exon sizing errors)
                    my $error5 = abs( $tmp5{cdna_end} - $expedge5 );
                    $error5 += 10.0 * $error5 / $$expected5{exon_size};
                    my $error3 = abs( $tmp3{cdna_begin} - $expedge3 );
                    $error3 += 10.0 * $error3 / $$expected3{exon_size};
                    my $errori =
                      abs( $intron_size - $expected_intron_size ) / 75.0;
                    my $sizing_error =
                      ( maxval( $error5, $error3 ) + $errori ) / 3.0;

           #if ( $sizing_error > $reflen / 10.0 && $sizing_error > 50 ) { next }

                    # score resulting protein
                    my $aa = DNA2AA($cdna);
                    $aa =~ s/\*$//;
                    my $aascore = score_profile( $aa, %refprofile );

                    my $stops_penalty =
                      calculate_stops_penalty( $$ref{id}, $aa );

                    my $score = $aascore * $stops_penalty + $splicing_quality -
                      $sizing_error;

                    if ($dbg) {

                        $aa =~ s/(.{1,60})/$1\n/g;
                        print
"-- $tmp5{dna_end}-$tmp3{dna_begin} ---------------------\n";
                        print "AA=$aa\n";

                        print
"exons=$tmp5{dna_begin}..$tmp5{dna_end}-$tmp3{dna_begin}..$tmp3{dna_end}\n"
                          . "  cdna=$tmp5{cdna_end}-$tmp3{cdna_begin} | $tmp5{subject_end}-$tmp3{subject_begin}\n"
                          . "  expd=$$expected5{cdna_end}-$$expected3{cdna_begin} | $$expected5{subject_end}-$$expected3{subject_begin}\n"
                          . "  adj=$adjust5 -> $expedge5..$expedge3 <- $adjust3\n"
                          . "  splicing_score= quality: $splicing_quality  sizing: $sizing_error 5=$error5  3=$error3  i=$errori\n"
                          . "  score=$score  stp penalty=$stops_penalty  aascore=$aascore\n"
                          . "  spl quality=$splicing_quality  donor_acceptor=$donor+$acceptor\n";
                        print
"--------------------------------------------------------\n";
                    }

                    # keep best scoring splice
                    if ( !defined $bestscore || $score > $bestscore ) {
                        $bestscore = $score;
                        @bestexons = ( \%tmp5, \%tmp3 );
                    }
                }
            }

            # use exons from best scoring splice
            if ($dbg) {
                print "NEW EXONS:\n";
                my $exonid = @spliced;
                print "\n--------------------------------------------\n";
                for my $exon (@bestexons) {
                    $exonid++;
                    print_hash( "exon$exonid", $exon );
                }
                print "--------------------------------------------\n";
            }

            if (@bestexons) {
                push @spliced, @bestexons;
            }
            else {
                $$exon3{spliced}        = -1;
                $$exon3{splice_quality} = -5;
                push @spliced, $exon5, $exon3;
            }

            if ($dbg) {
                print "\n--------------------------------------------\n";
                print "\nSPLICED:\n";
                my $exonid = 1;
                for my $exon (@spliced) {
                    print_hash( "exon$exonid", $exon );
                    $exonid++;
                }
                print "--------------------------------------------\n";
                print "\nREMAINING:\n";
                for my $exon (@exons) {
                    print_hash( "exon$exonid", $exon );
                    $exonid++;
                }
                print "--------------------------------------------\n";
            }
        }
    }

    # finalize and return exons
    set_cdna_coordinates(@spliced);

    if ($dbg) {
        print "\n--- final exons -----------------------\n";
        for my $exon (@spliced) {
            print_hash( "\nexon", $exon );
        }
        print "\n---------------------------------------\n";
    }

    return @spliced;
}

# returns gene exon/intron structure
# note: dna positions are offsets from wherever thstoprt of gene is found
sub get_gene_structure {
    my ($ref_id) = @_;
    my $dbg = DEBUG;

    #if ( get_reference_gene( $ref_id ) eq "CM1" ) { $dbg = 1 }

    my $ref = get_reference_seq($ref_id);
    my %structure;

    if ( !allow_splicing($ref_id) && !allow_ribosomal_slippage($ref_id) ) {
        $structure{spliceform} = "e" . 3 * ( $$ref{seqlen} + 1 );
        if ($dbg) { print "unspliced length $$ref{seqlen}\n" }
        my $exon;
        $$exon{fake_exon}     = 1;
        $$exon{subject_begin} = 1;
        $$exon{subject_end}   = $$ref{seqlen} + 1;
        $$exon{dna_begin}     = 0;
        $$exon{dna_end}       = 3 * $$ref{seqlen} + 2;
        $$exon{cdna_begin}    = 1;
        $$exon{cdna_end}      = 3 * $$ref{seqlen} + 3;
        $$exon{exon_size}     = $$exon{dna_end} - $$exon{dna_begin} + 1;

        $$exon{order} = 1;

        if ($dbg) {
            print_hash( "EXP X1", $exon );
        }
        my @exons = ($exon);
        $structure{exons}       = \@exons;
        $structure{num_exons}   = 1;
        $structure{cdna_length} = $$exon{cdna_end};
        $structure{dna_length}  = $$exon{dna_end};
    }

    else {
        $structure{spliceform} = get_splice_form($ref_id);
        if ( !defined $structure{spliceform} || $structure{spliceform} !~ /i/ )
        {
            die
"\nSpliced reference missing/malformed splice_form tag: $ref_id\n";
        }
        if ($dbg) { print "spliceform $structure{spliceform}\n" }
        my @tmp              = split /(i-{0,1}[0-9]+)/, $structure{spliceform};
        my $subject_position = 0;
        my $cdna_position    = 0;
        my $dna_offset       = 0;
        my $order            = 0;
        for my $f ( 0 .. @tmp - 1 ) {
            my $nucleotides = substr( $tmp[$f], 1 );
            if ( $tmp[$f] =~ /i/ ) {
                my $intron;
                $$intron{intron_size}   = $nucleotides;
                $$intron{subject_begin} = int( $subject_position + 0.7 );
                $$intron{subject_end}   = int($subject_position) + 1;
                $$intron{cdna_begin}    = $cdna_position + 1;
                $$intron{cdna_end}      = $cdna_position;
                $$intron{dna_begin}     = $dna_offset;
                $dna_offset += $nucleotides;
                $$intron{dna_end}     = $dna_offset - 1;
                $$intron{intron_size} =
                  $$intron{dna_end} - $$intron{dna_begin} + 1;

                $$intron{order} = $order;

                if ( defined $structure{introns} ) {
                    push @{ $structure{introns} }, $intron;
                }
                else {
                    my @introns = ($intron);
                    $structure{introns} = \@introns;
                }

                if ($dbg) {
                    print_hash( "EXP I$order", $intron );
                }
            }
            else {
                my $exon;
                $$exon{fake_exon}     = 1;
                $$exon{exon_size}     = $nucleotides;
                $$exon{subject_begin} = int($subject_position) + 1;
                $subject_position += $nucleotides / 3.0;
                if ( $subject_position - int($subject_position) > 0.7 ) {
                    $subject_position = int($subject_position) + 1;
                }
                if ( $subject_position < $$exon{subject_begin} ) { $subject_position = $$exon{subject_begin} }
                $$exon{subject_end} = int($subject_position);
                $$exon{cdna_begin}  = $cdna_position + 1;
                $cdna_position += $nucleotides;
                $$exon{cdna_end}  = $cdna_position;
                $$exon{dna_begin} = $dna_offset;
                $dna_offset += $nucleotides;
                $$exon{dna_end}   = $dna_offset - 1;
                $$exon{exon_size} = $$exon{dna_end} - $$exon{dna_begin} + 1;

                $order++;
                $$exon{order} = $order;

                if ( defined $structure{exons} ) {
                    push @{ $structure{exons} }, $exon;
                }
                else {
                    my @exons = ($exon);
                    $structure{exons} = \@exons;
                }

                if ($dbg) {
                    print_hash( "EXP X$order", $exon );
                }
            }
        }
        if ( is_reference_partial3($ref_id) ) {
            pop @{ $structure{exons} };
            pop @{ $structure{introns} };
        }
        $structure{num_exons}   = scalar @{ $structure{exons} };
        $structure{cdna_length} = $cdna_position;
        $structure{dna_length}  = $dna_offset;
    }

    return \%structure;
}

sub project_gene {
    my ($hsp) = @_;
    my $dbg = DEBUG;
    #if ( $$hsp{subject_id} =~ /^AFY97832.1$/ ) { $dbg = 1 }

    my $structure = get_gene_structure( $$hsp{subject_id} );
    my $expected = expected_exon( init_exon_from_hsp($hsp), $structure );
    if ( !defined $expected ) {
        return ( $$hsp{query_left}, $$hsp{query_right} );
    }

    if ( $dbg ) {
        #print_hash( "structure", $structure );
        print_hash( "expected", $expected );
        print_blasthits( 0, $hsp );
    }
    my $offset5 = $$expected{dna_begin} + 3 *
      ( $$hsp{subject_begin} - $$expected{subject_begin} );
    my $offset3 = $$structure{dna_length} - ( $$expected{dna_end} + 3 *
          ( $$hsp{subject_end} - $$expected{subject_end} ) );
    if ( $dbg ) {
        print "offsets 5': $offset5  3': $offset3\n";
    }
    if ( $$hsp{orientation} == 1 ) {
        if ( $dbg ) {
            print "return ( " . join( ", ", ( $$hsp{query_left} - $offset5, $$hsp{query_right} + $offset3 ) ) . " )\n";
        }
        return ( $$hsp{query_left} - $offset5, $$hsp{query_right} + $offset3 );
    }
    else {
        if ( $dbg ) {
            print "return ( " . join( ", ", ( $$hsp{query_left} - $offset3, $$hsp{query_right} + $offset5 ) ) . " )\n";
        }
        return ( $$hsp{query_left} - $offset3, $$hsp{query_right} + $offset5 );
    }
}

sub expected_exon {
    my ( $exon, $structure ) = @_;

    my $exonsize = $$exon{subject_end} - $$exon{subject_begin} + 1;
    my @expected = @{ $$structure{exons} };

   #print "find exon $exonsize | $$exon{subject_end} - $$exon{subject_begin}\n";

    my $bestoverlap = 0;
    my $bestexon;
    for my $i ( 1 .. @expected ) {
        my $exp     = $expected[ $i - 1 ];
        my $left    = maxval( $$exp{subject_begin}, $$exon{subject_begin} );
        my $right   = minval( $$exp{subject_end}, $$exon{subject_end} );
        my $overlap = $right - $left + 1;

        #print_hash( "$i: $overlap", $exp );

        if ( $overlap > $bestoverlap ) {
            $bestexon    = $exp;
            $bestoverlap = $overlap;
        }
    }

    #if ( ! defined $bestexon ) {
    #    print_hash( "exon", $exon );
    #    print_hash( "structure", $structure );
    #}

    return $bestexon;
}

sub gene_form_mismatch {
    my ($gene) = @_;

    my $structure = get_gene_structure( $$gene{ref_id} );
    my $subjmin   = $$gene{ref_begin};
    my $subjmax   = $$gene{ref_end};

    my $mismatch  = 0;
    my $lastorder = 0;
    for my $exon ( @{ $$gene{exons} } ) {
        my $expected = expected_exon( $exon, $structure );
        if ( !defined $expected ) {
            $mismatch +=
              ( $$exon{subject_end} - $$exon{subject_begin} + 1 ) / 2.0;

#print "exon ?: $$exon{subject_begin}-$$exon{subject_end}  expected: ?-?  min-max=$subjmin-$subjmax\n";
            next;
        }

        if ( $$expected{order} == $lastorder ) { next }
        $lastorder = $$expected{order};

        my $subject_end = $$exon{subject_end};
        if ( defined $$exon{fused_subject_end} ) {
            $subject_end = $$exon{fused_subject_end};
        }
        my $subject_begin = $$exon{subject_begin};
        if ( defined $$exon{fused_subject_begin} ) {
            $subject_begin = $$exon{fused_subject_begin};
        }

        else {
            $mismatch +=
              abs( $subject_end - minval( $subjmax, $$expected{subject_end} ) )
              / 2.0;
            $mismatch +=
              abs(
                $subject_begin - maxval( $subjmin, $$expected{subject_begin} ) )
              / 2.0;
        }

#print "exon $lastorder: $subject_begin-$subject_end  expected: $$expected{subject_begin}-$$expected{subject_end}  min-max=$subjmin-$subjmax\n";
    }

    return $mismatch;
}

#
sub extract_ribosomalslippage_exons {
    my ( $genome, $rawgene ) = @_;
    my $dbg = DEBUG;
    #if ( $$rawgene{subject_id} =~ /NP_598309.1/ ) { $dbg = 1 }
    #if ( `whoami` =~ /jhoover/ )  { $dbg = 1 }

    # initialization
    my $structure = get_gene_structure( $$rawgene{subject_id} );
    if ($dbg) {
        print_hash( "structure", $structure );
    }
    my $ref   = get_reference_seq( $$rawgene{subject_id} );
    my $ori   = $$rawgene{orientation};
    my @exons = extract_unspliced_exons( $genome, $rawgene );
    if ($dbg) {
        my $x = 1;
        for my $exon (@exons) {
            print_hash( "in$x", $exon );
            $x++;
        }
    }

    my ( $exon5, $exon3 );
    my $x = 0;
    for my $exon ( sort { $$a{subject_end} <=> $$b{subject_end} } @exons ) {
        $x++;
        my $expected = expected_exon( $exon, $structure );
        if ($dbg) {
            print_hash( "exon$x",     $exon );
            print_hash( "structure",  $structure );
            print_hash( "expected$x", $expected );
        }
        if ( defined $expected ) {
            if ( $$expected{order} == 1 ) {
                $exon5 = $exon;
                $$exon5{proj_dna} = $$exon5{dna_end} + 3 * $ori *
                  ( $$expected{subject_end} - $$exon{subject_end} );
                $$exon5{proj_sbj} = $$expected{subject_end};
                if ($dbg) {
                    print_hash( "exon5", $exon5 );
                }
            }
            elsif ( $$expected{order} == 2 ) {
                $exon3 = $exon;
                $$exon3{proj_dna} = $$exon3{dna_begin} - 3 * $ori *
                  ( $$exon{subject_begin} - $$expected{subject_begin} );
                $$exon3{proj_sbj} = $$expected{subject_begin};
                if ($dbg) {
                    print_hash( "exon3", $exon3 );
                }
                last;
            }
        }
    }

    # find position of ribosomal slippage
    if ( defined $exon5 && defined $exon3 ) {
        my ( $slippage_motif, $slippage_shift, $slippage_offset ) =
          describe_ribosomal_slippage( $$rawgene{subject_id} );
        my $slippage_motiflen = length($slippage_motif);
        if ($dbg) {
            print "slippage motif: $slippage_motif  shift: $slippage_shift  offset: $slippage_offset\n";
            print_hash( "exon5", $exon5 );
            print_hash( "exon3", $exon3 );
        }

        my $slip_point;
        my $slipframe = $$exon5{frame} + $ori * $slippage_shift;
        if ( $ori > 0 ) {
            if ( $slipframe < 1 ) { $slipframe += 3 }
            elsif ( $slipframe > 3 ) { $slipframe -= 3 }
        }
        else {
            if ( $slipframe > -1 ) { $slipframe -= 3 }
            elsif ( $slipframe < -3 ) { $slipframe += 3 }
        }

        if ($dbg) {
            print
"startframe: $$exon5{frame}  nextexpected: $slipframe  nextfound: $$exon3{frame}\n";
        }

        # determine range to search for slippage
        if ( $$exon3{frame} = $slipframe ) {
            my $bufleft = REF_BUFFER;
            if ( $$exon5{proj_sbj} - $bufleft < 1 ) {
                $bufleft = $$exon5{proj_sbj} - 1;
            }
            my $refleft = $$exon5{proj_sbj} - $bufleft;

            my $bufright = REF_BUFFER;
            if ( $$exon3{proj_sbj} - 1 + $bufright > $$rawgene{subject_length} )
            {
                $bufright = $$rawgene{subject_length} - $$exon3{proj_sbj} + 1;
            }
            my $refright = $$exon3{proj_sbj} - 1 + $bufright;

            my $bufbegin  = $$exon5{proj_dna} + $ori - $ori * 3 * $bufleft;
            my $slipbegin = $$exon5{proj_dna} + $ori - $ori * 3 * SLIP_RANGE;
            if ( $ori == 1 ) {
                while ( $bufbegin < 1 ) {
                    $bufbegin += 3;
                    $refleft++;
                    if ( $bufbegin > $slipbegin ) { $slipbegin = $bufbegin }
                }
            }
            else {
                while ( $bufbegin > $$genome{seqlen} ) {
                    $bufbegin -= 3;
                    $refleft++;
                    if ( $bufbegin < $slipbegin ) { $slipbegin = $bufbegin }
                }
            }

            my $bufend  = $$exon3{proj_dna} - $ori + $ori * 3 * $bufright;
            my $slipend = $$exon3{proj_dna} - $ori + $ori * 3 * SLIP_RANGE;
            if ( $ori == 1 ) {
                while ( $bufend > $$genome{seqlen} ) {
                    $bufend -= 3;
                    $refright--;
                    if ( $bufend < $slipend ) { $slipend = $bufend }
                }
            }
            else {
                while ( $bufend < 1 ) {
                    $bufend += 3;
                    $refright--;
                    if ( $bufend > $slipend ) { $slipend = $bufend }
                }
            }
            if ($dbg) {
                print "slippage $slipbegin-$slipend ($ori)\n";
                print "  buffer $bufbegin-$bufend ($ori)\n";
                print " subject $refleft-$refright\n";
            }

            # try the possible points of slippage,
            # keep the one that yields the closest match to the reference
            my $subref =
              substr( $$ref{sequence}, $refleft, $refright - $refleft + 1 );
            my %refprofile = profile_peptide($subref);

            $slipbegin -= $ori * $slippage_motiflen;
            $slipend += $ori * $slippage_motiflen;
            my $slippos   = $slipbegin;
            my $subgenome =
              subsequence( $$genome{sequence}, $slipbegin, $slipend );
            my $sublen = length($subgenome);
            for my $slippage ( find_regexp( $slippage_motif, $subgenome, 1, $sublen, 1 ) ) {
                my $slippos = $$slippage{end} + $ori * $slippage_offset;
                if ( $dbg ) { print "range: $slipbegin-$slipend  motif: $$slippage{begin}-$$slippage{end} $$slippage{string}  slippos $slippos = $$slippage{end} - $ori * $slippage_offset\n" }
                if ( $ori == 1 ) {
                    $slippos = $slippos + $slipbegin - 1;
                }
                else {
                    $slippos = $slipend + 1 - $slippos;
                }
                my $slip2 = $slippos - $ori * ( $slippage_shift + 1 );
                my $framecheck = ( abs( $slip2 - $bufbegin ) + 1 ) % 3;
                if ( $dbg ) { print "bufbegin $bufbegin  slippos $slippos  slip2 $slip2  framecheck $framecheck\n" }
                if ( $framecheck == 0 ) {
                    if ( $slip2 <= $$genome{seqlen} && $slip2 >= 1 ) {
                        if ($dbg) {
                            print
    "($slipbegin..$slipend) $$slippage{begin}-$$slippage{end} = $$slippage{string} =~ /$slippage_motif$/ => $slippos,$slip2\n";
                        }
                        my $penalty = abs( $$exon5{proj_dna} - $slip2 ) / 30.0;
                        my $reward  = $$slippage{string};
                        $reward =~ s/N//gi;
                        $reward = length($reward);
                        my $pep =
                          DNA2AA(
                                subsequence( $$genome{sequence}, $bufbegin, $slip2 )
                              . subsequence( $$genome{sequence}, $slippos, $bufend )
                          );
                        my $aascore = score_profile( $pep, %refprofile );
                        my $score   = $aascore - $penalty + $reward;
    
                        if ( !defined $slip_point || $score > $$slip_point{score} )
                        {
                            $$slip_point{score}    = $score;
                            $$slip_point{position} = $slippos;
                        }
                        if ($dbg) {
                            print
    "slippos=$bufbegin..$slippos,$slip2..$bufend  score=$score aascore=$aascore penalty=$penalty reward=$reward\npep=$pep\nref=$subref\n";
                        }
                    }
                }
            }
        }
        delete $$exon5{proj_dna};
        delete $$exon5{proj_sbj};
        delete $$exon3{proj_dna};
        delete $$exon3{proj_sbj};
        if ( defined $slip_point ) {
            my $in_gap = in_gap( $genome, $$slip_point{position} - 4 * $ori,
                $$slip_point{position} );
            adjust_exon_end( $exon5,
                $$slip_point{position} - $ori * ( $slippage_shift + 1 ), $ref );
            $$exon3{slippage} = 1 + $in_gap;
            adjust_exon_begin( $exon3, $$slip_point{position} );
        }
    }

    # format and return exons
    set_cdna_coordinates(@exons);
    if ($dbg) {
        my $x = 1;
        for my $exon (@exons) {
            print_hash( "out$x", $exon );
            $x++;
        }
    }
    return @exons;
}

# extract exons for gene with RNA editing
sub extract_rnaediting_exons {
    my ( $genome, $rawgene ) = @_;
    my $dbg = DEBUG;
    #if ( `whoami` =~ /jhoover/ && $$rawgene{subject_id} eq "298388487-I" ) { $dbg = 1 }

    # initialization
    my $rna_edit = get_rna_edit( $$rawgene{subject_id} );
    my ( $sz, $old, $new, $note ) = split /\//, $rna_edit;
    my $fs = rna_editing_frameshift( $$rawgene{subject_id} );

    my $ref        = get_reference_seq( $$rawgene{subject_id} );
    my %refprofile = profile_peptide(
        subsequence(
            $$ref{sequence}, $$rawgene{subject_begin},
            $$rawgene{subject_end}
        )
    );

    # convert HSPs to exons
    my $ori          = $$rawgene{orientation};
    my $genome_begin = 1;
    my $genome_end   = $$genome{seqlen};
    if ( $ori == -1 ) {
        my ( $genome_begin, $genome_end ) = ( $genome_end, $genome_begin );
    }
    my @exons = extract_unspliced_exons( $genome, $rawgene );
    {
        my @tmp;
        for my $exon (@exons) {
            if ( $$exon{in_gap} ) {
                my $prev = pop @tmp;
                $$exon{dna_begin}          = $$prev{dna_begin};
                $$exon{subject_begin}      = $$prev{subject_begin};
                $$exon{orig_dna_begin}     = $$prev{orig_dna_begin};
                $$exon{orig_subject_begin} = $$prev{prog_subject_begin};
                delete $$exon{in_gap};
                while ( ( $ori * ( $$exon{dna_end} - $$exon{dna_begin} ) + 1 ) %
                    3 != 0 )
                {
                    $$exon{dna_end} += $ori;
                }
                push @tmp, $exon;
            }
            else {
                push @tmp, $exon;
            }
        }
        set_cdna_coordinates(@tmp);
        @exons = @tmp;
    }

    if ($dbg) {
        print
          "\nRNA EDITING: \"$rna_edit\" fs=$fs old=$old new=$new note=$note\n";
        print_blasthits( 0, @{ $$rawgene{hsps} } );
        for my $exon (@exons) {
            print_hash( "\nin exon", $exon );
        }
    }

    # extend exons to projected start and stop codons
    my @extexons;
    for my $exon (@exons) {
        my %tmp = %$exon;
        push @extexons, \%tmp;
    }
    if ( $extexons[0]{subject_begin} > 1 ) {
        my $dnabegin = $extexons[0]{dna_begin} - 3 * $ori *
          ( $extexons[0]{subject_begin} - 1 );
        my $sbjbegin = 1;
        while ( $ori * $dnabegin < $genome_begin ) {
            $dnabegin += 3 * $ori;
            $sbjbegin++;
        }
        $extexons[0]{dna_begin}     = $dnabegin;
        $extexons[0]{subject_begin} = $sbjbegin;
    }
    if ( $extexons[ @extexons - 1 ]{subject_end} < $$ref{seqlen} ) {
        my $dnaend = $extexons[0]{dna_end} + 3 * $ori *
          ( $$ref{seqlen} - $extexons[ @extexons - 1 ]{subject_end} );
        my $sbjend = $$ref{seqlen};
        while ( $ori * $dnaend > $genome_end ) {
            $dnaend -= 3 * $ori;
            $sbjend--;
        }
        $extexons[0]{dna_end}     = $dnaend;
        $extexons[0]{subject_end} = $sbjend;
    }

    # find position where editing occurs
    my $edit_point;
    my $bestscore =
      0.9 * score_profile( DNA2AA( cdna_from_exons( $genome, @exons ) ),
        %refprofile );
    for my $i ( 0 .. @extexons - 1 ) {

        # check framing
        if ( $fs != 0 && $i < @extexons - 1 ) {
            my $beginframe = $extexons[$i]{frame};
            my $editframe  = $beginframe - $ori * $fs;
            if ( $beginframe > 0 ) {
                if ( $editframe < 1 ) { $editframe += 3 }
                elsif ( $editframe > 3 ) { $editframe -= 3 }
            }
            else {
                if ( $editframe > -1 ) { $editframe -= 3 }
                elsif ( $editframe < -3 ) { $editframe += 3 }
            }
            if ($dbg) {
                print
"beginframe: $beginframe  nextexpected: $editframe  nextfound: $extexons[$i+1]{frame}\n";
            }
            if ( $extexons[ $i + 1 ]{frame} ne $editframe ) { next }
        }

      # fuse with frameshifted exon (RNA editing will compensate for frameshift)
        my @tmp = @extexons;
        if ( $i < @extexons - 1 ) {
            my %tmpexon = %{ $extexons[$i] };
            $tmpexon{dna_end}     = $extexons[ $i + 1 ]{dna_end};
            $tmpexon{subject_end} = $extexons[ $i + 1 ]{subject_end};
            $tmp[$i]              = \%tmpexon;
            $tmp[ $i + 1 ]        = undef;
            @tmp                  = remove_undefs(@tmp);
        }
        set_cdna_coordinates(@tmp);
        my $cdna = cdna_from_exons( $genome, @tmp );

        # try the possible points of edit,
        # keep the one that yields the closest match to the reference
        my $begin = $tmp[$i]{cdna_begin};

        #        if ( defined $tmp[$i]{in_gap} ) {
        #            $begin = $tmp[ $i - 1 ]{cdna_begin};
        #        }
        my $end = $tmp[$i]{cdna_end};
        if ($dbg) {
            print "edit point range: $begin-$end\n";
            print "cdna=$cdna\n";
        }
        my @matches = find_regexp( $old, $cdna, $begin, $end, 1 );
        my @gapmatches = find_regexp( "(N{20,})", $cdna, $begin, $end );
        for my $match (@gapmatches) {
            $$match{in_gap} = 1;
            if ( $sz < 0 ) {
                my $s = $sz;
                while ( $s < 0 ) {
                    $$match{begin}++;
                    $s++;
                    if ( $s < 0 ) {
                        $$match{end}--;
                        $s++;
                    }
                }
            }
            $$match{string} =
              lpad( "N", $$match{end} - $$match{begin} + 1, "N" );
            push @matches, $match;
        }

        for my $match (@matches) {
            if ($dbg) {
                print_hash( "\nedit site?", $match );
            }
            my $newseq;
            my $align;
            my $dna_xref;
            if ( defined $$match{in_gap} ) {
                $newseq =
                  lpad( $$match{string}, length( $$match{string} ) + $sz, "N" );
                if ( $sz > 0 ) {
                    my $half = int( length( $$match{string} ) / 2.0 + 0.5 );
                    my $old  =
                      lpad( lpad( "N", $half, "N" ) . lpad( "-", $sz, "-" ),
                        length($newseq), "N" );
                    $align = "$old\n$newseq";
                }
                elsif ( $sz < 0 ) {
                    my $half = int( length($newseq) / 2.0 + 0.5 );
                    my $new  =
                      lpad( lpad( "N", $half, "N" ) . lpad( "-", $sz, "-" ),
                        length( $$match{string} ), "N" );
                    $align = "$$match{string}\n$new";
                }
                else {
                    $align = "$$match{string}\n$newseq";
                }
                if ($dbg) {
                    print "gap match\n$align\n";
                }
            }
            else {
                if ( $$match{string} !~ /$old/ ) { next }
                my @submatches = ( $1, $2, $3, $4, $5 );
                $newseq = $new;
                $newseq =~ s/\$//g;
                $dna_xref = $newseq;
                $dna_xref =~ s/[a-z]/0/gi;
                if ( $newseq =~ /oldseq/i ) {
                    $newseq =~ s/oldseq/$$match{string}/gi;
                    if ($dbg) {
                        print
"oldseq|\"$$match{string}\": edit from $$match{string} to $newseq\n";
                    }
                }
                for my $i ( 1 .. 5 ) {
                    if ( $newseq =~ /$i/ ) {
                        $newseq =~ s/$i/$submatches[$i-1]/g;
                        my $dnamatch =
                          lpad( "1", length( $submatches[ $i - 1 ] ), "1" );
                        $dna_xref =~ s/$i/$dnamatch/g;
                        if ($dbg) {
                            print
"$i|\"$submatches[$i-1]\": edit from $$match{string} to $newseq\n";
                        }
                    }
                }
                $align = align_edit_strings( $$match{string}, $newseq );
                if ($dbg) {
                    print "edit match\n$align\n";
                }
            }

            my $newcdna =
                substr( $cdna, 0, $$match{begin} - 1 ) . $newseq
              . substr( $cdna, $$match{end} );

            my $pep = DNA2AA($newcdna);
            my $score = score_profile( $pep, %refprofile );
            if ( $$match{string} =~ /N/i ) {
                my $matlen  = length( $$match{string} );
                my $nstring = $$match{string};
                $nstring =~ s/[^Nn]//g;
                my $nlen  = length($nstring);
                my $scale = sqrt( 1.0 - $nlen / $matlen / 15.0 );

#print "match $$match{string}  N discount $scale  orig score $score  discounted " . ( $scale * $score ) . "\n";
                $score = 0.95 * $score;
            }

            if ($dbg) {
                print
"edit $rna_edit match $$match{begin}-$$match{end} change $$match{string} to $newseq\n";
                print "ref=$$ref{sequence}\n";
                print "pep=$pep\n";
                print "i=$i  score=$score  best=$bestscore\n";
            }

            if ( $score > $bestscore ) {
                $bestscore               = $score;
                $$edit_point{score}      = $score;
                $$edit_point{i}          = $i;
                $$edit_point{cdna_begin} = $$match{begin};
                $$edit_point{cdna_end}   =
                  $$match{end} + length($newseq) - length( $$match{string} );
                ( $$edit_point{dna_begin}, $$edit_point{dna_end} ) =
                  cdna_coords_to_dna( $$match{begin}, $$match{end}, \@tmp );
                $$edit_point{newseq} = $newseq;
                $$edit_point{in_gap} = $$match{in_gap};

                if ($dbg) {
                    print_hash( "\nedit_point", $edit_point );
                    print "\n";
                }
                $$edit_point{align} = $align;
            }
        }
    }

    if ( defined $edit_point ) {
        my $i = $$edit_point{i};
        if ($dbg) {
            print
"\n*** i=$i  editpos D=$$edit_point{dna_begin}-$$edit_point{dna_end}  C=$$edit_point{cdna_begin}-$$edit_point{cdna_end}  score=$$edit_point{score}\n";
        }
        $exons[$i]{rna_edit}{cdna_begin} = $$edit_point{cdna_begin};
        $exons[$i]{rna_edit}{cdna_end}   = $$edit_point{cdna_end};
        $exons[$i]{rna_edit}{dna_begin}  = $$edit_point{dna_begin};
        $exons[$i]{rna_edit}{dna_end}    = $$edit_point{dna_end};
        $exons[$i]{rna_edit}{newseq}     = $$edit_point{newseq};
        $exons[$i]{rna_edit}{align}      = $$edit_point{align};
        $exons[$i]{rna_edit}{note}       = $note;
        if ( defined $$edit_point{in_gap} ) {

            if ( defined $exons[$i]{rna_edit}{note} ) {
                $exons[$i]{rna_edit}{note} .=
                  "; likely site of rna editing falls within gap";
            }
            else {
                $exons[$i]{rna_edit}{note} =
                  "likely position of rna editing falls within gap";
            }
            $exons[$i]{rna_edit}{in_gap} = 1;
        }
        if ( $i < @exons - 1 ) {
            $exons[$i]{dna_end}     = $exons[ $i + 1 ]{dna_end};
            $exons[$i]{subject_end} = $exons[ $i + 1 ]{subject_end};
            $exons[ $i + 1 ]        = undef;
            @exons                  = remove_undefs(@exons);
            if ($dbg) {
                for my $exon (@exons) {
                    print_hash( "\nmid exon", $exon );
                }
            }
        }

        if ($dbg) {
            print "$i. $exons[$i]{dna_begin}-"
              . ( $$edit_point{dna_begin} - $ori ) . "="
              . subsequence(
                $$genome{sequence},
                $exons[$i]{dna_begin},
                $$edit_point{dna_begin} - $ori
              )
              . "\n"
              . "   $$edit_point{dna_begin}-$$edit_point{dna_end}=$$edit_point{newseq}\n"
              . ( $$edit_point{dna_end} + $ori )
              . "-$exons[$i]{dna_end}="
              . subsequence(
                $$genome{sequence},
                $$edit_point{dna_end} + $ori,
                $exons[$i]{dna_end}
              )
              . "\n";
            print "align=\n$$edit_point{align}\n";
        }
    }

    # format and return exons
    if ($dbg) {
        for my $exon (@exons) {
            print_hash( "\nout exon", $exon );
        }
    }
    set_cdna_coordinates(@exons);
    if ( $exons[ @exons - 1 ]{cdna_end} % 3 != 0 ) {
        while ( $exons[ @exons - 1 ]{cdna_end} % 3 != 0 ) {
            $exons[ @exons - 1 ]{cdna_end}++;
            $exons[ @exons - 1 ]{dna_end} += $ori;
        }
        $exons[ @exons - 1 ]{subject_end}++

    }
    return @exons;
}

sub align_edit_strings {
    my ( $old, $new ) = @_;

    my $lo = length($old);
    my $ln = length($new);

    # substitution
    if ( $lo == $ln ) {
    }

    # insertion
    elsif ( $lo < $ln ) {
        my $po = 0;
        while ( $po < $lo && substr( $old, $po, 1 ) eq substr( $new, $po, 1 ) )
        {
            $po++;
        }
        if ( $po > 0 ) {
            my $o = substr( $old, 0, $po ) . lpad( "-", $ln - $lo, "-" );
            if ( $po < $lo ) {
                $old = $o . substr( $old, $po + 1 );
            }
            else {
                $old = $o;
            }
        }
        else {
            $old = lpad( "-", $ln - $lo, "-" ) . $old;
        }
    }

    # deletion
    else {
        my $pn = 0;
        while ( $pn < $ln && substr( $old, $pn, 1 ) eq substr( $new, $pn, 1 ) )
        {
            $pn++;
        }
        if ( $pn > 0 ) {
            my $n = substr( $new, 0, $pn ) . lpad( "-", $lo - $ln, "-" );
            if ( $pn < $ln ) {
                $new = $n . substr( $new, $pn + 1 );
            }
            else {
                $new = $n;
            }
        }
        else {
            $new = lpad( "-", $lo - $ln, "-" ) . $new;
        }
    }

    # return result
    return "$old\n$new";
}

# return unspliced exons (frameshift may cause return of multiple exons (separated by the shifts)
sub extract_unspliced_exons {
    my ( $genome, $rawgene ) = @_;
    my $dbg = DEBUG;
    #if ( $$rawgene{subject_id} eq "NP_062883.2" ) { $dbg = 1 }

    # get reference sequence for comparison
    my $ref       = get_reference_seq( $$rawgene{subject_id} );
    my $structure = get_gene_structure( $$rawgene{subject_id} );

    # covert hsps to exons, fuse exons
    my $ori = $$rawgene{orientation};
    my @exons;
    my $lastexon;
    for my $hsp ( sort { compare_qry_positions( $a, $b ) } @{ $$rawgene{hsps} } ) {
        if ( $dbg ){
            print "next HSP\n";
            print_blasthits( 0, $hsp );
        }
        my $exon = init_exon_from_hsp($hsp);
        my $expected = expected_exon( $exon, $structure );
        if ($dbg) {
            print "exon F: $$exon{frame}  D: $$exon{dna_begin}-$$exon{dna_end}  S: $$exon{subject_begin}-$$exon{subject_end}\n";
            print_hash( "expected", $expected );
        }
        if ( !@exons ) {
            push @exons, $exon;
            $lastexon = $exon;
            if ($dbg) {
                print "  add exon #"
                  . ( @exons + 1 )
                  . " FR1: $$exon{frame} DNA: $$exon{dna_begin}-$$exon{dna_end} ($ori) REF: $$exon{subject_begin}-$$exon{subject_end}\n";
            }
        }
        else {
            my $lastexon = pop @exons;
            my $lastexpected = expected_exon( $lastexon, $structure );
            if ($dbg) {
                print_hash( "lastexpected", $lastexpected );
            }
            if (   defined $lastexpected
                && defined $expected
                && $$lastexpected{order} == $$expected{order} )
            {
                my @fused = fuse_exons( $genome, $ref, $lastexon, $exon );
                for my $exon (@fused) {
                    push @exons, $exon;
                    if ($dbg) {
                        print "  add exon #"
                          . ( @exons + 1 )
                          . " FR2: $$exon{frame} DNA: $$exon{dna_begin}-$$exon{dna_end} ($ori) REF: $$exon{subject_begin}-$$exon{subject_end}\n";
                    }
                }
            }
            else {
                push @exons, $lastexon, $exon;
                if ($dbg) {
                    print "  add exon #"
                      . ( @exons + 1 )
                      . " FR3: $$exon{frame} DNA: $$exon{dna_begin}-$$exon{dna_end} ($ori) REF: $$exon{subject_begin}-$$exon{subject_end}\n";
                }
            }
        }
    }
    if ( $dbg ) { print "set_cdna_coordinates and return\n" }
    # format and return exons
    set_cdna_coordinates(@exons);
    return @exons;
}

sub init_exon_from_hsp {
    my ($hsp) = @_;

    my %exon;
    $exon{subject_begin} = $$hsp{subject_begin};
    $exon{subject_end}   = $$hsp{subject_end};
    if ( $$hsp{orientation} == 1 ) {
        $exon{dna_begin} = $$hsp{query_left};
        $exon{dna_end}   = $$hsp{query_right};
    }
    else {
        $exon{dna_begin} = $$hsp{query_right};
        $exon{dna_end}   = $$hsp{query_left};
    }
    $exon{orig_dna_begin}     = $exon{dna_begin};
    $exon{orig_dna_end}       = $exon{dna_end};
    $exon{orig_subject_begin} = $exon{subject_begin};
    $exon{orig_subject_end}   = $exon{subject_end};
    $exon{frame}              = $$hsp{query_frame};
    $exon{frame} =~ s/\+//;
    $exon{orientation} = $$hsp{orientation};

    if ( defined $$hsp{in_gap} ) { $exon{in_gap} = $$hsp{in_gap} }

    return \%exon;
}

sub set_cdna_coordinates {
    my (@exons) = @_;
    my $cdnalen = 0;
    for my $exon (@exons) {
        $$exon{cdna_begin} = $cdnalen + 1;
        $cdnalen += abs( $$exon{dna_end} - $$exon{dna_begin} ) + 1;
        $$exon{cdna_end} = $cdnalen;
        if ( defined $$exon{rna_edit} ) {
            my $editspan =
              $$exon{rna_edit}{cdna_end} - $$exon{rna_edit}{cdna_begin};
            $$exon{rna_edit}{cdna_begin} =
              abs( $$exon{rna_edit}{dna_begin} - $$exon{dna_begin} ) +
              $$exon{cdna_begin};
            $$exon{rna_edit}{cdna_end} =
              $$exon{rna_edit}{cdna_begin} + $editspan;
            $cdnalen += $editspan -
              abs( $$exon{rna_edit}{dna_end} - $$exon{rna_edit}{dna_begin} );
            $$exon{cdna_end} = $cdnalen;

        }

#        if ( $verbose ) {
#            print "EXON D: $$exon{dna_begin}-$$exon{dna_end}  C: $$exon{cdna_begin}-$$exon{cdna_end}\n";
#        }
    }
}

sub cdna_from_exons {
    my ( $genome, @exons ) = @_;

    my $cdna = "";
    my $rna_edit;
    for my $exon (@exons) {
        if ( defined $$exon{rna_edit} ) {
            $rna_edit = $$exon{rna_edit};
            my $seqbefore =
              subsequence( $$genome{sequence}, $$exon{dna_begin},
                $$exon{rna_edit}{dna_begin} - $$exon{orientation} );
            my $seqafter =
              subsequence( $$genome{sequence},
                $$exon{rna_edit}{dna_end} + $$exon{orientation},
                $$exon{dna_end} );
            $cdna .= $seqbefore . $$exon{rna_edit}{newseq} . $seqafter;

#print "$$exon{dna_begin}-" . ($$exon{rna_edit}{dna_begin}-$$exon{orientation} ) . "=$seqbefore\n"
#    . "$$exon{rna_edit}{dna_begin}-$$exon{rna_edit}{dna_end}=$$exon{rna_edit}{newseq}\n"
#    . ($$exon{rna_edit}{dna_end}+$$exon{orientation}) . "-$$exon{dna_end}=$seqafter\n";
#print "\nAA=" . DNA2AA($cdna) . "\n";
        }
        else {
            my $exon_seq =
              subsequence( $$genome{sequence}, $$exon{dna_begin},
                $$exon{dna_end} );

            #print "$$exon{dna_begin}-$$exon{dna_end}=$exon_seq\n";
            $cdna .= $exon_seq;
        }
    }
    return $cdna;
}

sub adjust_exon_begin {
    my ( $exon, $newbegin ) = @_;
    if ( $newbegin == $$exon{dna_begin} ) { return }

    my $ori = $$exon{orientation};
    my $sbjdelta = int( $ori * ( $newbegin - $$exon{dna_begin} ) / 3.0 + 0.5 );
    $$exon{dna_begin} = $newbegin;
    $$exon{subject_begin} += $sbjdelta;
    if ( $$exon{subject_begin} < 1 ) { $$exon{subject_begin} = 1 }
}

sub adjust_exon_cdnabegin {
    my ( $exon, $newbegin ) = @_;
    if ( $newbegin == $$exon{cdna_begin} ) { return }

    my $ori = $$exon{orientation};
    my $sbjdelta = int( ( $newbegin - $$exon{cdna_begin} ) / 3.0 + 0.5 );
    $$exon{dna_begin} += $ori * ( $newbegin - $$exon{cdna_begin} );
    $$exon{cdna_begin} = $newbegin;
    $$exon{subject_begin} += $sbjdelta;
    if ( $$exon{subject_begin} < 1 ) { $$exon{subject_begin} = 1 }
    $$exon{fuzzy_begin} = 0;
}

sub adjust_exon_end {
    my ( $exon, $newend, $ref ) = @_;
    if ( $newend == $$exon{dna_end} ) { return }

    my $ori = $$exon{orientation};
    my $sbjdelta = int( $ori * ( $newend - $$exon{dna_end} ) / 3.0 + 0.5 );
    $$exon{dna_end} = $newend;
    $$exon{subject_end} += $sbjdelta;
    if ( $$exon{subject_end} > $$ref{seqlen} ) {
        $$exon{subject_end} = $$ref{seqlen};
    }
}

sub adjust_exon_cdnaend {
    my ( $exon, $newend, $reflen ) = @_;
    if ( $newend == $$exon{cdna_end} ) { return }

    my $ori = $$exon{orientation};
    my $sbjdelta = int( ( $newend - $$exon{cdna_end} ) / 3.0 + 0.5 );
    $$exon{dna_end} += $ori * ( $newend - $$exon{cdna_end} );
    $$exon{cdna_end} = $newend;
    $$exon{subject_end} += $sbjdelta;
    if ( $$exon{subject_end} > $reflen ) {
        $$exon{subject_end} = $reflen;
    }
    $$exon{fuzzy_end} = 0;
}

# fuse exons into one (unless there is a frameshift, then grow th edges till they meet)
# currently handles only TWO exons
sub fuse_exons {
    my ( $genome, $ref, @exons ) = @_;
    my $dbg = DEBUG;
    #if ( $$ref{id} eq "NP_062883.2" ) { $dbg = 1 }

    if ($dbg) {
        print
"\n\n-----------------------------------------------------\nFUSE_EXONS\n";
    }

    for my $exon (@exons) {
        if ($dbg) {
            print_hash( "in", $exon );
        }
        if ( !defined $$exon{orig_dna_begin} ) {
            $$exon{orig_dna_begin} = $$exon{dna_begin};
        }
        if ( !defined $$exon{orig_dna_end} ) {
            $$exon{orig_dna_end} = $$exon{dna_end};
        }
        if ( !defined $$exon{orig_subject_begin} ) {
            $$exon{orig_subject_begin} = $$exon{subject_begin};
        }
        if ( !defined $$exon{orig_subject_end} ) {
            $$exon{orig_subject_end} = $$exon{subject_end};
        }
    }
    my %tmp      = %{ $exons[0] };
    my $lastexon = \%tmp;
    my %tmp2     = %{ $exons[1] };
    my $exon     = \%tmp2;
    my $ori      = $exons[0]{orientation};

    # exons in same frame - merge
    if ( $$exon{frame} eq $$lastexon{frame} ) {
        $$lastexon{orig_dna_end} = $$exon{orig_dna_end};
        $$lastexon{dna_end}      = $$exon{dna_end};
        $$lastexon{subject_end}  = $$exon{subject_end};
        return ($lastexon);
    }

    # exons out of frame
    else {

        # re-blast the two exons
        # (sometimes blast over-extends HSPs)
        {
            my $dbeg = $$lastexon{dna_end};
            my $sbeg = $$lastexon{subject_end};
            my $move = 40;
            while ( $move && $sbeg > $$lastexon{subject_begin} ) {
                $sbeg--;
                $dbeg -= 3 * $ori;
                $move--;
                if ( $dbg ) { print "re-blast $$lastexon{subject_begin}-$$exon{subject_end}: sbeg $sbeg  dbeg $dbeg\n" }
            }

            my $dend = $$exon{dna_begin};
            my $send = $$exon{subject_begin};            
            $move = 40;
            while ( $move && $send < $$exon{subject_end} ) {
                $send++;
                $dend += 3 * $ori;
                $move--;
                if ( $dbg ) { print "re-blast $$lastexon{subject_begin}-$$exon{subject_end}: send $send  dend $dend\n" }
            }

            my $rsz = ( $send - $sbeg + 1 );
            if ( $rsz > ( abs( $dend - $dbeg ) + 1 ) / 3.0 ) {
                $rsz = ( abs( $dend - $dbeg ) + 1 ) / 3.0;
            }
            my $rawgene = {
                orientation => $ori,
                query_id => $$genome{id},
                query_definition => $$genome{id},
                query_length => $$genome{seqlen},
                subject_id => $$ref{id},
                subject_definition => $$ref{defline},
                subject_length => $$ref{seqlen}
            };
            my @tst = sort { $$a{subject_begin} <=> $$b{subject_begin} } 
                remove_conflicting_hsps2 ( 
                    find_small_hsps(
                        $genome, $dbeg, $dend, $rawgene, 1E-5, 1E-5,
                        10,  60, $ref, $sbeg, $send ) );
            if ( $dbg ) {
                print "possible re-blast HSPs D: $dbeg-$dend ($ori)  S: $sbeg-$send\n";
                print_blasthits( 0, @tst );
            }
            while ( @tst && $tst[0]{query_frame} == $$lastexon{frame} ) {
                $$lastexon{subject_end} = $tst[0]{subject_end};
                $$lastexon{dna_end} = $ori == 1 ? $tst[0]{query_right} : $tst[0]{query_left};
                shift @tst;
            }
            while ( @tst && $tst[@tst-1]{query_frame} == $$exon{frame} ) {
                $$exon{subject_begin} = $tst[@tst-1]{subject_begin};
                $$exon{dna_begin} = $ori == 1 ? $tst[@tst-1]{query_left} : $tst[0]{query_right};
                pop @tst
            }
#            if ( $dbg ) {
#                print_hash( "adjusted lastexon", $lastexon );
#                print_hash( "adjusted exon",     $exon );
#            }
        }

        # gaps between exons?
        my $gbegin = $$lastexon{dna_end} - 10 * $ori;
        while ( defined $$genome{in_gap}{$gbegin} ) {
            $gbegin -= $ori;
        }
        my $gend = $$exon{dna_begin} + 10 * $ori;
        while ( defined $$genome{in_gap}{$gend} ) {
            $gend += $ori;
        }
        my @gaps = find_gaps( $genome, $gbegin, $gend );
        if ($dbg) {
            print "find gaps in $gbegin-$gend = " . @gaps . "\n";
            for my $gap (@gaps) {
                print_hash( "gap", $gap );
            }
        }
        if ( $dbg ) {
            print_hash( "out-of-frame lastexon", $lastexon );
            print_hash( "out-of-frame exon",     $exon );
            print "\n" . @gaps . " gaps between exons\n";
        }

        # no gaps
        # try to figure where the frameshift occurs
        if ( ! @gaps ) {
            
            
            # calculate the frameshift offset between exons
             my $adj = ( abs( $$exon{frame} ) - abs( $$lastexon{frame} ) + 4 ) % 3;
             $adj = $ori * $adj;
            if ( $adj == 2 ) {
                $adj = -1;
            }
            elsif ( $adj == -2 ) {
                $adj = 1;
            }
            
            # get the reference substring for scoring the potential frameshift positions
            my $subref = subsequence( $$ref{sequence}, $$lastexon{subject_begin}, $$exon{subject_end} );
            my %subpro = profile_peptide( $subref );
            my $subscore;
            my $subpos;
            
            # determine the range of dna positions to test as point of frameshift
            my $maxlast = $$lastexon{dna_end};
            my $minnext = $$exon{dna_begin};
            while ( $ori * $maxlast > $ori * $$exon{dna_begin} ) {
                $maxlast -= $ori * 3;
            }
            while ( $ori * $minnext < $ori * $$lastexon{dna_end} ) {
                $minnext += $ori * 3;
            }
            if ( $dbg ) { print "test range $maxlast to $minnext  adjustment: $adj\n" }
            my $pos = $maxlast;
            while ( $ori * $pos <= $ori * $minnext ) {
                my $lastcds = subsequence( $$genome{sequence}, $$lastexon{dna_begin}, $pos );
                my $nextcds = subsequence( $$genome{sequence}, $pos + $adj, $$exon{dna_end} );
                my $subaa = DNA2AA( $lastcds . $nextcds );
                my $score = score_profile( $subaa, %subpro );
                my $stop = $subaa;
                while ( $stop =~ /\*/ ) {
                    $stop =~ s/\*//;
                    $score--;
                }
                if ( $dbg ) {
                    print "$$lastexon{dna_begin}-$pos+" . ( $pos + $adj ) . "-$$exon{dna_end}\n"
                        . "subaa: $subaa\n"
                        . "refaa: $subref\n"
                        . "score: $score";
                }
                if ( ! defined $subscore || $score > $subscore ) {
                    if ( $dbg ) { print "*" }
                    $subpos = $pos;
                    $subscore = $score;
                }
                if ( $dbg ) { print "\n" }
                $pos += 3 * $ori
            }
            
            # adjust the exon edges to the highest scoring position
            $$lastexon{subject_end} += int( $ori * ( $subpos - $$lastexon{dna_end} ) / 3.0 + $ori * 0.5 );
            $$lastexon{dna_end} = $subpos;
    
            $subpos += $adj;
            $$exon{subject_begin} -= int( $ori * ( $$exon{dna_begin} - $subpos ) / 3.0 + $ori * 0.5 );
            $$exon{dna_begin} = $subpos;
    
            if ($dbg) {
                print_hash( "unoverlapped lastexon", $lastexon );
                print_hash( "unoverlapped exon",     $exon );
            }
            return ( $lastexon, $exon );
        }

        # gaps exist - assume frameshift caused by mis-sized gap
        # add an exon to represent the gap
        else {
            my $gapexon;
            $$gapexon{in_gap}      = 1;
            $$gapexon{frame}       = $$lastexon{frame};
            $$gapexon{orientation} = $ori;
            if ( $ori == 1 ) {
                $$gapexon{gap_dna_begin} = $gaps[0]{begin};
                $$gapexon{gap_dna_end}   = $gaps[ @gaps - 1 ]{end};
            }
            else {
                $$gapexon{gap_dna_begin} = $gaps[ @gaps - 1 ]{end};
                $$gapexon{gap_dna_end}   = $gaps[0]{begin};
            }

            # grow the original exons to meet the gap
            while (
                $ori * $$lastexon{dna_end} > $ori * $$gapexon{gap_dna_begin} )
            {
                $$lastexon{dna_end} -= 3 * $ori;
                $$lastexon{subject_end}--;
            }
            while (
                $ori * $$lastexon{dna_end} < $ori * $$gapexon{gap_dna_begin} )
            {
                $$lastexon{dna_end} += 3 * $ori;
                $$lastexon{subject_end}++;
            }
            while ( $ori * $$exon{dna_begin} < $ori * $$gapexon{gap_dna_end} ) {
                $$exon{dna_begin} += 3 * $ori;
                $$exon{subject_begin}++;
            }
            while ( $ori * $$exon{dna_begin} > $ori * $$gapexon{gap_dna_end} ) {
                $$exon{dna_begin} -= 3 * $ori;
                $$exon{subject_begin}--;
            }
            while ( $$lastexon{subject_end} >= $$exon{subject_begin} ) {
                $$lastexon{subject_end}--;
                $$exon{subject_begin}++;
            }
            if ($dbg) {
                print_hash( "togap lastexon", $lastexon );
                print_hash( "togap exon",     $exon );
            }

            # define the "gap exon"
            $$gapexon{dna_begin}     = $$lastexon{dna_end} + $ori;
            $$gapexon{dna_end}       = $$exon{dna_begin} - $ori;
            $$gapexon{subject_begin} = $$lastexon{subject_end} + 1;
            $$gapexon{subject_end}   = $$exon{subject_begin} - 1;

            $$gapexon{gapsize} =
              $ori * ( $$gapexon{gap_dna_end} - $$gapexon{gap_dna_begin} ) + 1;
            if ($dbg) {
                print_hash( "init gapexon", $gapexon );
            }

            # adjust the "gap exon" for the frameshift
            my $dnasize =
              $ori * ( $$gapexon{dna_end} - $$gapexon{dna_begin} ) + 1;
            my $expected_dnasize =
              3 * ( $$gapexon{subject_end} - $$gapexon{subject_begin} + 1 );
            my $expected_gapsize =
              $expected_dnasize - $dnasize + $$gapexon{gapsize};
            while ( $dnasize % 3 != 0 ) {
                $dnasize--;
                $$gapexon{dna_end} -= $ori;
                if ($dbg) {
                    print
                      "gap framing: size $dnasize  end $$gapexon{dna_end} \n";
                }
            }
            $$gapexon{expectedsize} = $expected_gapsize;
            my $fused_dna_begin     = $$lastexon{dna_begin};
            my $fused_dna_end       = $$exon{dna_end};
            my $fused_subject_begin = $$lastexon{subject_begin};
            my $fused_subject_end   = $$exon{subject_end};
            for my $x ( $lastexon, $gapexon, $exon ) {
                $$x{fused_dna_begin}     = $fused_dna_begin;
                $$x{fused_dna_end}       = $fused_dna_end;
                $$x{fused_subject_begin} = $fused_subject_begin;
                $$x{fused_subject_end}   = $fused_subject_end;
            }
            if ($dbg) {
                print_hash( "unoverlapped lastexon", $lastexon );
                print_hash( "gap exon", $gapexon );
                print_hash( "unoverlapped exon",     $exon );
            }
            return ( $lastexon, $gapexon, $exon );
        }

    }
}

# quick and dirty routines for comparing two proteins based on kmer counts
# generate a protein kmer profile
sub profile_peptide {
    my ($peptide) = @_;
    my $binsz = 3;

    my $pep = uc($peptide);
    $pep =~ s/\*$//;
    my $peplen = length($pep);

    my %profile;
    $profile{aa}     = $pep;
    $profile{aa_len} = $peplen;

    for my $sz ( 1 .. $binsz ) {
        for my $i ( 0 .. length($pep) - $sz ) {
            my $bin = substr( $pep, $i, $sz );
            $profile{bins}{$bin}++;
        }
    }

    return %profile;
}

# compare a protein to a previous computed kmer profile
sub score_profile {
    my ( $peptide, %refprofile ) = @_;

    my $dbg = $refprofile{dbg};
    if ( !defined $dbg ) { $dbg = 0 }

    my %refbins;
    for my $bin ( keys %{ $refprofile{bins} } ) {
        $refbins{$bin} = $refprofile{bins}->{$bin};
    }

    my %pepprofile = profile_peptide($peptide);
    my %pepbins    = %{ $pepprofile{bins} };

    if ($dbg) {
        print "pep=$pepprofile{aa}\nref=$refprofile{aa}\n";
        print "peplen=$pepprofile{aa_len}  ref=$refprofile{aa_len}\n";
    }

    my %matches;
    for my $bin ( keys %pepbins ) {
        if ( exists $refbins{$bin} ) {
            my $match =
              $pepbins{$bin} < $refbins{$bin} ? $pepbins{$bin} : $refbins{$bin};
            $matches{$bin} = $match;
            $pepbins{$bin} -= $match;
            $refbins{$bin} -= $match;
        }
    }

    # try to match-up kmers from gaps (X) with reference kmers
    my $continue = 1;
    while ($continue) {
        $continue = 0;
        for my $pep ( sort { $pepbins{$b} <=> $pepbins{$a} } keys %pepbins ) {
            if ( $pepbins{$pep} < 1 ) { last }
            if ( $pep !~ /X/i )       { next }
            my $Xs   = $pepbins{$pep};
            my $test = $pep;
            $test =~ s/X/./gi;
            $test =~ s/\*/\\*/g;
            for my $ref ( sort { $refbins{$b} <=> $refbins{$a} } keys %refbins )
            {
                if ( $refbins{$ref} < 1 ) { last }
                if ( $ref !~ /$test/ )    { next }
                my $match =
                    $pepbins{$pep} < $refbins{$ref}
                  ? $pepbins{$pep}
                  : $refbins{$ref};
                $pepbins{$pep} -= $match;
                $refbins{$ref} -= $match;
                if ( $pepbins{$pep} > 0 ) { $continue = 1 }
                last;
            }
        }
    }

    my $matchcount = 0;
    for my $value ( values %matches ) {
        $matchcount += 2 * $value;
    }
    my $mismatchcount = 0;
    for my $value ( values %pepbins ) {
        $mismatchcount += $value;
    }
    for my $value ( values %refbins ) {
        $mismatchcount += $value;
    }

    my $len_mismatch = abs( $pepprofile{aa_len} - $refprofile{aa_len} ) / 2.0;

    my $score = $refprofile{aa_len} * $matchcount /
      ( $matchcount + $mismatchcount + $len_mismatch + 0.000001 );
    if ($dbg) {
        print
"1. match=$matchcount  mismatch=$mismatchcount  len_mismatch=$len_mismatch  score=$score\n";
    }

#    if ( uc substr( $peptide, 0, 1 ) eq "M" ) {
#        $len_mismatch -= 5;
#        if ( $len_mismatch < 0 ) { $len_mismatch = 0 }
#        $match += 5;
#        $score = 100.0 * $match / ( $match + $mismatch + $len_mismatch );
#        if ( $dbg ) { print "2. match=$match  mismatch=$mismatch  len_mismatch=$len_mismatch  score=$score\n" }
#    }

    my $rescore = int( 2.0 * $score + 0.5 ) / 2.0;
    $rescore += ( $score - $rescore ) / 100.0;
    return $rescore;
}

# calculate quality of alignment between reference and gene
sub score_pep_vs_reference {
    my ( $pep_id, $peptide, $ref_id, $ref_pep, $stats ) = @_;
    my $dbg = DEBUG;

    if ($dbg) {
        print
          "\n----------------------------------\nscore_pep_vs_reference in\n";
    }
    my $ref_length = length($ref_pep);
    if ( substr( $ref_pep, $ref_length - 1, 1 ) eq "*" ) {
        $ref_length--;
    }
    my $pep = uc $peptide;

    my $alignment =
      align_pep_to_reference( "predicted", $pep, "reference", $ref_pep );
    if ( substr( $pep, length($pep) - 1, 1 ) eq "*" ) {
        $pep = substr( $pep, 0, length($pep) - 1 );
    }
    $pep =~ s/X{2,}//g;
    my $pep_length = length($pep);

    if ($dbg) { print "Q $pep_id vs S $ref_id alignment=\n$alignment\n" }
    my (
        $alignlen,    $num_covered,   $num_identical,
        $num_similar, $num_sim25,     $num_reftrunc5,
        $num_refgap,  $num_reftrunc3, $num_overextended
      )
      = get_alignment_counts( "predicted", "reference", $alignment );
    if ( $num_reftrunc5 == 0 && $num_overextended > 0 ) {
        $num_reftrunc5 = -$num_overextended;
    }
    if ($dbg) {
        print
"align len $alignlen  identical $num_identical  similar $num_similar  ref len $ref_length  covered $num_covered\n";
    }
    (
        $$stats{alignlen},      $$stats{num_covered},
        $$stats{num_identical}, $$stats{num_similar},
        $$stats{num_sim25},     $$stats{num_reftrunc5},
        $$stats{num_refgap},    $$stats{num_reftrunc3},
        $$stats{num_overextended}
      )
      = (
        $alignlen,    $num_covered,   $num_identical,
        $num_similar, $num_sim25,     $num_reftrunc5,
        $num_refgap,  $num_reftrunc3, $num_overextended
      );

    if ( $alignlen == 0 ) {
        $$stats{num_refcovered}      = 0;
        $$stats{num_refsimilar}      = 0;
        $$stats{num_refsim25}        = 0;
        $$stats{num_refidentical}    = 0;
        $$stats{num_pepcoverage}     = 0;
        $$stats{pct_refcoverage}     = 0;
        $$stats{pct_refsimilarity}   = 0;
        $$stats{pct_refsimilarity25} = 0;
        $$stats{pct_refidentity}     = 0;
        $$stats{pct_reftrunc5}       = 0;
        $$stats{pct_refgap}          = 0;
        $$stats{pct_reftrunc3}       = 0;
        $$stats{match_quality}       = 0;
    }
    else {
        $$stats{num_refcovered}      = $num_covered;
        $$stats{pct_refcoverage}     = 100.0 * $num_covered / $ref_length;
        $$stats{num_refsimilar}      = $num_similar;
        $$stats{pct_refsimilarity}   = 100.0 * $num_similar / $alignlen;
        $$stats{num_refsim25}        = $num_sim25;
        $$stats{pct_refsimilarity25} =
          100.0 * $num_sim25 / minval( 25.0, $alignlen );
        $$stats{num_refidentical} = $num_identical;
        $$stats{pct_refidentity}  = 100.0 * $num_identical / $alignlen;
        $$stats{num_reftrunc5}    = $num_reftrunc5;
        $$stats{pct_reftrunc5}    = 100.0 * $num_reftrunc5 / $ref_length;
        $$stats{num_refgap}       = $num_refgap;
        $$stats{pct_refgap}       = 100.0 * $num_refgap / $ref_length;
        $$stats{num_reftrunc3}    = $num_reftrunc3;
        $$stats{pct_reftrunc3}    = 100.0 * $num_reftrunc3 / $ref_length;
        $$stats{match_quality}    =
          1.00 * $num_identical + 0.50 * ( $num_similar - $num_identical ) +
          0.15 * ( $num_covered - $num_similar ) - 0.15 *
          ( $pep_length - $num_covered );
    }

    $$stats{pep_alignment} = $alignment;

    #print "score_pep_vs_reference out\n";
}

sub round_refstats {
    my ($stats) = @_;
    for my $statname ( get_refstat_attributes() ) {
        if ( $statname eq "alignment" ) { next }
        if ( defined $$stats{$statname} ) {
            $$stats{$statname} = int( 10.0 * $$stats{$statname} + 0.5 ) / 10.0;
            if ( index( $$stats{$statname}, "." ) < 0 ) {
                $$stats{$statname} .= ".0";
            }
        }
    }
}

sub get_alignment_counts {
    my ( $gene_id, $ref_id, $alignment ) = @_;
    my $dbg = DEBUG;

    #if ( getlogin() eq "jhoover" ) { $dbg = 1 }

    # parse the genome:reference alignment
    if ( !defined $alignment ) { $alignment = " \n \n " }
    my @tmp = split /\n/, $alignment;
    shift @tmp;
    shift @tmp;
    shift @tmp;

    my $i       = 0;
    my $genestr = "";
    my $refstr  = "";
    my $compstr = "";
    my $prevstr = "";
    my $margin;
    my $field = 0;

    for my $line (@tmp) {
        if ($dbg) {
            print "gene $gene_id ref $ref_id line $line\n";
        }
        my $seqid = $line;
        $seqid =~ s/ .*$//;
        if ( $seqid eq $gene_id ) {
            $field = 1;
            my $alignstr = $line;
            $alignstr =~ s/^[^ ]*  *//;
            $genestr .= $alignstr;
            $margin = index( $line, $alignstr );
            $prevstr = $alignstr;
        }
        elsif ( $seqid eq $ref_id ) {
            $field = 2;
            my $alignstr = $line;
            $alignstr =~ s/^[^ ]*  *//;
            $refstr .= $alignstr;
            $margin = index( $line, $alignstr );
            $prevstr = $alignstr;
        }
        else {
            $field++;
            if ( $field == 3 ) {
                if ( $line !~ /^ *$/ ) {
                    my $alignstr = substr( $line, $margin );
                    $compstr .= $alignstr;
                }
                else {
                    my $alignstr = $prevstr;
                    $alignstr =~ s/[^ ]/ /g;
                    $compstr .= $alignstr;
                }
            }
            else {
                $field = 4;
                next;
            }
        }
    }
    if ($dbg) {
        print "\nalignment\ngene=$genestr\ncomp=$compstr\n ref=$refstr\n";
    }
    my $tmp = $refstr;
    $tmp =~ s/-//g;
    my $reflen = length($tmp);

    while ( length($compstr) < length($refstr) ) {
        $compstr .= "          ";
    }
    if ($dbg) { print "reflen=$reflen\n" }

    # check for 5' and 3' truncation
    my $num_reftrunc5 = 0;
    if ( $genestr =~ /^(\-+)/ ) {
        $num_reftrunc5 = length($1);
        $genestr       = substr( $genestr, $num_reftrunc5 );
        $compstr       = substr( $compstr, $num_reftrunc5 );
        $refstr        = substr( $refstr, $num_reftrunc5 );
    }
    my $num_reftrunc3 = 0;
    if ( $genestr =~ /(\-+)$/ ) {
        $num_reftrunc3 = length($1);
        $genestr = substr( $genestr, 0, length($genestr) - $num_reftrunc3 );
        $compstr = substr( $compstr, 0, length($compstr) - $num_reftrunc3 );
        $refstr  = substr( $refstr, 0, length($refstr) - $num_reftrunc3 );
    }
    if ($dbg) {
        print
          "\nafter truncation\ngene=$genestr\ncomp=$compstr\n ref=$refstr\n";
    }

    # count and remove internal gaps on gene
    my $pos     = 0;
    my $tmpgene = $genestr;
    $genestr = "";
    my $tmpcomp = $compstr;
    $compstr = "";
    my $tmpref = $refstr;
    $refstr = "";
    my $num_refgap = 0;

    while ( $tmpgene =~ /(XX-+XX)/i ) {
        my $x1 = $1;
        my $x2 = $x1;
        $x2      =~ s/-/X/g;
        $tmpgene =~ s/$x1/$x2/;
    }
    for my $frag ( split /([xX]{2,})/, $tmpgene ) {

        if ( $frag !~ /[xX]{2,}/ ) {
            $genestr .= substr( $tmpgene, $pos, length($frag) );
            $compstr .= substr( $tmpcomp, $pos, length($frag) );
            $refstr  .= substr( $tmpref,  $pos, length($frag) );
        }
        else {
            $num_refgap += length($frag);
        }
        $pos += length($frag);
    }
    if ($dbg) {
        print "\nafter gene gaps\ngene=$genestr\ncomp=$compstr\n ref=$refstr\n";
    }
    my $alignlen = length($refstr);

    # calculate coverage and matches
    my $num_overextended = 0;
    if ( $refstr =~ /^(-+)M/i ) {
        $num_overextended = length($1);
        if ( substr( $refstr, $num_overextended, 1 ) !~ /M/i ) {
            $num_overextended = 0;
        }
        if ($dbg) {
            print "gene: "
              . substr( $genestr, 0, $num_overextended + 5 ) . "\n";
            print " ref: " . substr( $refstr, 0, $num_overextended + 5 ) . "\n";
        }
    }

    $refstr =~ s/-//g;
    my $num_covered = length($refstr);

    my $comp25 = substr( $compstr, 0, 25 );
    $comp25 =~ s/ //g;
    my $num_similar25 = length($comp25);
    $comp25 =~ s/\.//g;
    my $lowsim = $num_similar25 - length($comp25);
    $num_similar25 -= $lowsim / 4.0;

    $compstr =~ s/ //g;
    my $num_similar = length($compstr);
    $compstr =~ s/\.//g;
    $lowsim = $num_similar - length($compstr);
    $num_similar -= $lowsim / 4.0;

    $compstr =~ s/[^*]//g;
    my $num_identical = length($compstr);

    if ($dbg) {
        print
"alignlen=$alignlen  num_reftrunc5=$num_reftrunc5  num_refgap=$num_refgap  num_reftrunc3=$num_reftrunc3\n";
        print
"num_covered=$num_covered  num_similar=$num_similar  num_identical=$num_identical\n";
    }

    #print "get_alignment_counts out\n";
    return (
        $alignlen,    $num_covered,   $num_identical,
        $num_similar, $num_similar25, $num_reftrunc5,
        $num_refgap,  $num_reftrunc3, $num_overextended
    );
}


sub mapto_polyprotein {
    my ( $gene, $polyprotein, $matpepdb, $min_gene_size, $min_gene_coverage ) = @_;
    my $dbg = DEBUG;
    my $verbose = get_parameter("verbose");

    if ( $verbose || $dbg ) {
        print "polyprotein gene\n";
        print_genehits( $gene );
    }

    # initialization
    if ( !defined $polyprotein ) { $polyprotein = $$gene{protein} }
    if ( ! defined $polyprotein ) { return undef }

    my $polylen = length($polyprotein);
    if ( substr( $polyprotein, $polylen - 1, 1 ) eq "*" ) { $polylen-- }

    if ( !defined $matpepdb ) {
        $matpepdb = get_matpep_db( $$gene{ref_id} );
        if ( !defined $matpepdb ) {
            my @empty = ();
            return \@empty;
        }
    }
    if ( $matpepdb ne $refmatdb ) {
        %refmat   = loadFasta($matpepdb);
        $refmatdb = $matpepdb;
    }
    my $refdb = basename( $matpepdb );

    my $matpep_mincoverage   = get_parameter("mature_pep_mincoverage");
    my $matpep_minsimilarity = get_parameter("mature_pep_minsimilarity");
    my $matpep_minidentity   = get_parameter("mature_pep_minidentity");
    my $vigorspace           = get_parameter("vigorspace");
    my $jcvi_rules           = get_parameter("jcvi_rules");


    # use reference sequence to fill in gaps
    my $ref = get_reference_seq( $$gene{ref_id} );
    if ($dbg) {
        my $poly = $polyprotein;
        $poly =~ s/(.{1,60})/$1\n/g;
        print "\n>POLY length=$polylen\n$poly\n";
        my $refpoly = $$ref{sequence};
        $refpoly =~ s/(.{1,60})/$1\n/g;
        print "\n>$$ref{id} length=$$ref{seqlen}\n$refpoly\n";
    }
    my @gaps;
    for my $gap ( find_regexp( "[xX]{2,}", $polyprotein ) ) {
        my ( $ref_begin, $ref_end ) =
          pep_coords_to_ref( $$gap{begin}, $$gap{end}, $$gene{codon_start},
            $$gene{exons}, $ref );
        if ($dbg) {
            print "GAP $$gap{begin}-$$gap{end} = REF $ref_begin, $ref_end\n";
        }
        if ( !defined $ref_begin ) { next }

        my $old    = subsequence( $polyprotein, $$gap{begin}, $$gap{end} );
        my $oldlen = length($old);
        my $new    = lc( subsequence( $$ref{sequence}, $ref_begin, $ref_end ) );
        my $newlen = length($new);
        while ( $newlen < $oldlen ) {
            my $half = int( $newlen / 2 ) + 1;
            $new = substr( $new, 0, $half ) . "X" . substr( $new, $half );
            $newlen++;
        }
        while ( $newlen > $oldlen ) {
            my $half = int( $newlen / 2 ) + 1;
            $new = substr( $new, 0, $half - 1 ) . substr( $new, $half );
            $newlen--;
        }

        $polyprotein =
            subsequence( $polyprotein, 1, $$gap{begin} - 1 ) . $new
          . subsequence( $polyprotein, $$gap{end} + 1, $polylen );
        push @gaps, "$$gap{begin}|$$gap{end}|$old";
        if ($dbg) {
            print
              "GAP FILL: $$gap{begin}-$$gap{end}\n  >OLD\n$old\n  >NEW\n$new\n";
        }
    }

    # use reference sequence to extend polyprotein fragment
    my $fragoffset = 0;
    my $fraglen    = $polylen;
    if ( $$gene{start_truncation} ) {
        $polyprotein =
          lc( substr( $$ref{sequence}, 0, $$gene{ref_begin} - 1 ) )
          . $polyprotein;
        $polylen    = length($polyprotein);
        $fragoffset = $$gene{ref_begin} - 1;
    }
    if ( $$gene{stop_truncation} ) {
        $polyprotein .= lc( substr( $$ref{sequence}, $$gene{ref_end} ) );
        $polylen = length($polyprotein);
    }
    if ( $dbg && $polyprotein ne uc($polyprotein) ) {
        my $tmp = $polyprotein;
        $tmp =~ s/(.{1,60})/$1\n/g;
        print "updated polyprotein\n$tmp\n";
    }

    my $mpsbj = "$vigorspace/sbj";
    unlink $mpsbj;
    open( MPSBJ, ">$mpsbj" );
    print MPSBJ ">pol\n$polyprotein\n";
    close MPSBJ;

 #    system
 #"cd $vigorspace; $myBin/formatdb -i $mpsbj -l $vigorspace/formatdb_mpsbj.log";

    my $mpxml = "$vigorspace/xml";
    unlink $mpxml;

    &runCmd($vigorspace, "$myBin/blastall -p blastp -e 10 -W 1 -C 0 -F \"\" -v 0 -b 2000 -m 7 -d $matpepdb -i $mpsbj > $mpxml 2> /dev/null");

#"cd $vigorspace; $myBin/blastall -p blastp -e 10 -M BLOSUM45 -g F -F \"\" -v 0 -b 1 -m 7 -i $mpqry -d $mpsbj > $mpxml";
    my @tmp = parse_blastxml( $mpxml, 0 );

    #    if ($verbose) {
    #        print "SUBJECT HSPS\n";
    #        print_blasthits( 0, @tmp );
    #    }

    #    {
    #        my %best;
    #        for my $hit ( sort { $$a{evalue} <=> $$b{evalue} } @tmp ) {
    #            if ( !exists $best{ $$hit{subject_id} } ) {
    #                $best{ $$hit{subject_id} } = $hit;
    #            }
    #        }
    #        @tmp = sort { $$a{subject_begin} <=> $$b{subject_begin} } values %best;
    #    }
    #    push @tmp, @tmpend;
    my @subjs =
      sort { $$a{query_left} <=> $$b{query_left} }
      best_subject_hits( undef, \@tmp, 1, 1 );

    #    if ($verbose) {
    #        print "SUBJECT HITS\n";
    #        print_blasthits( 0, @subjs );
    #    }

    # invert subject and query
    my @rawmatpeps;
    for my $hit (@subjs) {

        #        ( $$hit{query_id}, $$hit{subject_id} ) =
        #          ( $$hit{subject_id}, $$hit{query_id} );
        #        ( $$hit{query_definition}, $$hit{subject_definition} ) =
        #          ( $$hit{subject_definition}, $$hit{query_definition} );
        #        ( $$hit{query_left}, $$hit{subject_begin} ) =
        #          ( $$hit{subject_begin}, $$hit{query_left} );
        #        ( $$hit{query_right}, $$hit{subject_end} ) =
        #          ( $$hit{subject_end}, $$hit{query_right} );
        #        ( $$hit{query_length}, $$hit{subject_length} ) =
        #          ( $$hit{subject_length}, $$hit{query_length} );
        #        ( $$hit{pct_scoverage}, $$hit{pct_qcoverage} ) =
        #          ( $$hit{pct_qcoverage}, $$hit{pct_scoverage} );

        # project coverage to end of reference
        $$hit{orig_query_left}    = $$hit{query_left};
        $$hit{orig_query_right}   = $$hit{query_right};
        $$hit{orig_subject_begin} = $$hit{subject_begin};
        $$hit{orig_subject_end}   = $$hit{subject_end};
        if ( $$hit{subject_begin} > 1 ) {
            my $orig_left = $$hit{query_left};
            my $adjust    = $$hit{subject_begin} - 1;
            if ( $adjust < 20 ) {
                $$hit{query_left}    -= $adjust;
                $$hit{subject_begin} -= $adjust;
            }
            if ( $$hit{query_left} < 1 ) {
                $$hit{subject_begin} += ( 1 - $$hit{query_left} );
                $$hit{query_left} = 1;
            }
            $$hit{alignment_length} += ( $orig_left - $$hit{query_left} );
        }
        if ( $$hit{subject_end} < $$hit{subject_length} ) {
            my $orig_right = $$hit{query_right};
            my $adjust     =
              minval( 10, $$hit{subject_length} - $$hit{subject_end} );
            if ( $adjust < 20 ) {
                $$hit{query_right} += $adjust;
                $$hit{subject_end} += $adjust;
            }
            if ( $$hit{query_right} > $polylen ) {
                $$hit{subject_end} -= ( $$hit{query_right} - $polylen );
                $$hit{query_right} = $polylen;
            }
            $$hit{alignment_length} += ( $$hit{query_right} - $orig_right );
        }
        $$hit{query_coverage}   = $$hit{query_right} - $$hit{query_left} + 1;
        $$hit{subject_coverage} = $$hit{subject_end} - $$hit{subject_begin} + 1;
        calculate_blast_percentages($hit);

        if (   $$hit{pct_identity} < $matpep_minidentity
            || $$hit{pct_similarity} < $matpep_minsimilarity )
        {
            next;
        }
        if ( !check_coverage( $matpep_mincoverage, $hit, $gene ) ) { next }

        # init "fuzzy" flags
        $$hit{fuzzy_begin} = 0;
        $$hit{fuzzy_end}   = 0;
        push @rawmatpeps, $hit;
    }
    if ($verbose) {
        print "init candidates\n";
        print_blasthits( 0, @rawmatpeps );
    }
    if ( !@rawmatpeps ) { return undef }

    # score edges
    my %products;
    for my $hit (@rawmatpeps) {
        my $product = get_defline_tag( $$hit{subject_definition}, "product" );
        my $edgeL   = $$hit{query_left};
        my $edgeR   = $$hit{query_right};
        my $weightL = 1.0;
        if ( $$hit{orig_query_left} != $$hit{query_left} ) {
            $weightL =
              ( $$hit{subject_length} - $$hit{orig_query_left} +
                  $$hit{query_left} ) / $$hit{subject_length} *
              $$hit{pct_identity} / 125.0;
        }
        my $weightR = 1.0;
        if ( $$hit{orig_query_left} != $$hit{query_left} ) {
            $weightR =
              ( $$hit{subject_length} - $$hit{query_right} +
                  $$hit{orig_query_right} ) / $$hit{subject_length} *
              $$hit{pct_identity} / 125.0;
        }
        $products{$product}{left}{$edgeL}  += $weightL;
        $products{$product}{right}{$edgeR} += $weightR;
        $products{$product}{total}++;
    }
    if ($dbg) {
        print "edges\n";
        for my $product ( sort { $a cmp $b } keys %products ) {
            print "  $product ($products{$product}{total})\n";
            print "    L:";
            for my $left ( sort { $a <=> $b }
                keys %{ $products{$product}{left} } )
            {
                print "  $left ($products{$product}{left}{$left})";
            }
            print "\n";
            print "    R:";
            for my $right ( sort { $a <=> $b }
                keys %{ $products{$product}{right} } )
            {
                print "  $right ($products{$product}{right}{$right})";
            }
            print "\n";
        }
    }

    # remove redundant hits
    @rawmatpeps = sort { $$a{query_left} <=> $$b{query_left} } @rawmatpeps;
    for my $i ( 0 .. @rawmatpeps - 2 ) {
        if ( !defined $rawmatpeps[$i] ) { next }
        for my $j ( $i + 1 .. @rawmatpeps - 1 ) {
            if ( !defined $rawmatpeps[$j] ) { next }
            if (
                abs(
                    $rawmatpeps[$i]{query_left} - $rawmatpeps[$j]{query_left}
                ) < 30
                && abs(
                    $rawmatpeps[$i]{query_right} - $rawmatpeps[$j]{query_right}
                ) < 30
              )
            {
                if (
                    (
                        $rawmatpeps[$i]{num_similar} +
                        $rawmatpeps[$i]{num_identical} / 10.0
                    ) / 1.1 > (
                        $rawmatpeps[$j]{num_similar} +
                          $rawmatpeps[$j]{num_identical} / 10.0
                    ) / 1.1
                  )
                {
                    $rawmatpeps[$j] = undef;
                }
                else {
                    $rawmatpeps[$i] = undef;
                    last;
                }
            }
            elsif (
                $rawmatpeps[$j]{query_left} - $rawmatpeps[$i]{query_left} > 30 )
            {
                last;
            }
        }
    }
    @rawmatpeps = remove_undefs(@rawmatpeps);

    #print "nr candidates\n";
    #print_blasthits( 0,  @rawmatpeps );

    # find linking edges
    my $maxgap = 30;
    for my $i ( 0 .. @rawmatpeps - 2 ) {
        for my $j ( $i + 1 .. @rawmatpeps - 1 ) {
            if (
                abs(
                    $rawmatpeps[$i]{query_right} + 1 -
                      $rawmatpeps[$j]{query_left}
                ) < $maxgap
              )
            {
                push @{ $rawmatpeps[$i]{right_links} }, $j;
                push @{ $rawmatpeps[$j]{left_links} },  $i;
            }
        }
    }

    # find best path through tree for each potential starting point
    # (mat pep without a neighbor to the left)
    #print "FINDING PATHS\n";
    my @bestpaths;
    for my $start (@rawmatpeps) {
        if ( exists $$start{left_links} ) { next }

        my @path;
        my @stack;
        my $bestpath;
        my $node = $start;
        while ( defined $node ) {
            push @path, $node;
            if ( defined $$node{right_links} && @{ $$node{right_links} } ) {
                my @links = @{ $$node{right_links} };
                $node = $rawmatpeps[ pop @links ];
                push @stack, \@links;
            }
            else {
                my $score = 0;
                my $right = $path[0]{query_left} - 1;
                for my $n (@path) {
                    $score +=
                      ( $$n{num_similar} + $$n{num_identical} / 10.0 ) / 1.1;
                    $score -= abs( $right + 1 - $$n{query_left} );
                    $right = $$n{query_right};
                }

#print "Path  $path[0]{query_left}-$path[@path-1]{query_right}  Size " . @path . "  Score $score\n";
#print_blasthits( 0,  @path );
                if ( !defined $bestpath || $score > $$bestpath{score} ) {
                    $$bestpath{left}  = $path[0]{query_left};
                    $$bestpath{right} = $path[ @path - 1 ]{query_right};
                    $$bestpath{width} =
                      $$bestpath{right} - $$bestpath{left} + 1;
                    $$bestpath{size}  = scalar @path;
                    $$bestpath{score} = $score;
                    my @savepath = @path;
                    $$bestpath{path} = \@savepath;
                }
                $node = undef;
                while ( @stack && !defined $node ) {
                    pop @path;
                    my $links = pop @stack;
                    if ( defined $links && @$links ) {
                        $node = $rawmatpeps[ pop @$links ];
                        push @stack, $links;
                    }
                }
            }
        }
        push @bestpaths, $bestpath;
    }
    @bestpaths = sort { $$a{left} <=> $$b{left} } @bestpaths;
    if ($verbose) {
        print "\nBEST PATHS\n";
        for my $bestpath (@bestpaths) {
            print
"Bestpath $$bestpath{left}-$$bestpath{right}  Size $$bestpath{size}  Score $$bestpath{score}\n";
            print_blasthits( 0, @{ $$bestpath{path} } );
        }
    }

    # remove redundant paths
    my $right = 0;
    for my $i ( 0 .. @bestpaths - 1 ) {
        if ( !defined $bestpaths[$i] ) { next }
        for my $j ( $i + 1 .. @bestpaths - 1 ) {
            if ( !defined $bestpaths[$j] ) { next }
            if ( $bestpaths[$j]{left} >= $bestpaths[$i]{right} ) { last }
            my $overlap =
              minval( $bestpaths[$i]{right}, $bestpaths[$j]{right} ) -
              maxval( $bestpaths[$i]{left}, $bestpaths[$j]{left} ) + 1;
            if ( $overlap >
                minval( $bestpaths[$i]{width}, $bestpaths[$j]{width} ) / 2.0 )
            {
                if ( $bestpaths[$i]{score} -
                    abs( $right + 1 - $bestpaths[$i]{left} ) >
                    $bestpaths[$j]{score} -
                    abs( $right + 1 - $bestpaths[$j]{left} ) )
                {
                    $bestpaths[$j] = undef;
                    next;
                }
                else {
                    $bestpaths[$i] = undef;
                    last;
                }
            }
        }
        if ( defined $bestpaths[$i] ) {
            $bestpaths[$i]{score} -= abs( $right + 1 - $bestpaths[$i]{left} );
            $right = $bestpaths[$i]{right};
        }
    }
    @bestpaths = remove_undefs(@bestpaths);
    if ($verbose) {
        print "\nFINAL PATHS\n";
        for my $bestpath (@bestpaths) {
            print
"Bestpath $$bestpath{left}-$$bestpath{right}  Size $$bestpath{size}  Score $$bestpath{score}\n";
            print_blasthits( 0, @{ $$bestpath{path} } );
        }
    }

    # adjust edges to eliminate gaps
    {
        my $maxgap = 15;
        @rawmatpeps = ();
        my @tmp;
        for my $bestpath (@bestpaths) {
            push @tmp, @{ $$bestpath{path} };
        }
        if ($verbose) {
            print "\nSELECTED MATURE PEPS\n";
            print_blasthits( 0, @tmp );
        }

        # make sure mature peptides start at first aa
        {
            my $raw = $tmp[0];

            # large gap implies a missing mature peptide
            my $maxgap = 1.05 * ( $$raw{subject_begin} - 1 ) + $maxgap;
            if ( $$raw{query_left} > $maxgap ) {
                $$raw{fuzzy_begin} = 1;
            }
            else {
                $$raw{query_left} = 1;
            }
        }

        # mature peptides should be adjacent
        # adjust alignment boundaries appropriately
        while ( @tmp > 1 ) {
            my $raw = shift @tmp;
            push @rawmatpeps, $raw;

            # large gap implies a missing mature peptide
            my $rawnext = $tmp[0];
            my $maxgap  =
              1.05 * ( $$raw{subject_length} - $$raw{subject_end} ) + 1.05 *
              ( $$rawnext{subject_begin} - 1 ) + $maxgap;
            my $gap = $$rawnext{query_left} - $$raw{query_right};
            if ( $gap > $maxgap ) {
                $$raw{fuzzy_end}       = 1;
                $$rawnext{fuzzy_begin} = 1;
            }
            else {
                adjust_matpep_edges( $polyprotein, \%refmat, $raw, $rawnext,
                    %products );
            }
        }

        # make sure mature peptides end at last aa
        {
            my $raw = $tmp[0];
            push @rawmatpeps, $raw;

            # large gap implies a missing mature peptide
            my $maxgap =
              1.05 * ( $$raw{subject_length} - $$raw{subject_end} ) + $maxgap;
            if ( $polylen + 1 - $$raw{query_right} > $maxgap ) {
                $$raw{fuzzy_end} = 1;
            }
            else {
                $$raw{query_right} = $polylen;
            }
        }
    }

    if ($verbose) {
        print "\nPOLISHED MATURE PEPS\n";
        print_blasthits( 0, @rawmatpeps );
    }

    # if partial polyprotein, adjust back to fragment
    if ( $fragoffset > 0 || $fraglen < $polylen ) {
        $polyprotein = substr( $polyprotein, $fragoffset, $fraglen );
        $polylen     = $fraglen;

        #print "frag offset-$fragoffset  length=$fraglen\n$polyprotein\n";
        my @tmp;
        for my $pep (@rawmatpeps) {

#print " in: $$pep{query_left}-$$pep{query_right} ($$pep{fuzzy_begin},$$pep{fuzzy_end})\n";
            $$pep{query_left}  -= $fragoffset;
            $$pep{query_right} -= $fragoffset;
            if ( $$pep{query_right} < 1 )       { next }
            if ( $$pep{query_left} > $fraglen ) { next }
            if ( $$pep{query_left} < 1 )        {
                $$pep{query_left}  = 1;
                $$pep{fuzzy_begin} = $$gene{start_truncation};

#print "md1: $$pep{query_left}-$$pep{query_right} ($$pep{fuzzy_begin},$$pep{fuzzy_end})\n";
            }
            if ( $$pep{query_right} > $fraglen ) {
                $$pep{query_right} = $fraglen;
                $$pep{fuzzy_end}   = $$gene{stop_truncation};

#print "md2: $$pep{query_left}-$$pep{query_right} ($$pep{fuzzy_begin},$$pep{fuzzy_end})\n";
            }
            push @tmp, $pep;

#print "out: $$pep{query_left}-$$pep{query_right} ($$pep{fuzzy_begin},$$pep{fuzzy_end})\n";
        }
        @rawmatpeps = @tmp;

        if ($verbose) {
            print "\nFRAGMENT MATURE PEPS\n";
            print_blasthits( 0, @rawmatpeps );
        }
    }

    # if gaps filled from reference, restore the gaps
    if (@gaps) {
        for my $gap (@gaps) {
            my ( $begin, $end, $old ) = split /\|/, $gap;
            my $cur = subsequence( $polyprotein, $begin, $end );

            $polyprotein =
                subsequence( $polyprotein, 1, $begin - 1 ) . $old
              . subsequence( $polyprotein, $end + 1, $polylen );

    #print "restore gap $gap $begin-$end\nfrom $cur\n  to $old\n$polyprotein\n";
        }
    }

    # return the final list of mature peptides
    my $pepno = 0;
    my @matpeps;
    for my $raw (@rawmatpeps) {
        $pepno++;

        my $pep;
        $$pep{pep_id}       = "$$gene{gene_id}.$pepno";
        $$pep{ref_id}       = $$raw{subject_id};
        $$pep{ref_db}       = $refdb;
        $$pep{product_name} =
          get_reference_product( $$raw{subject_definition} );
        $$pep{pep_note} = get_reference_note( $$raw{subject_definition} );

        #        $$pep{product_name} = $$raw{subject_definition};
        #        $$pep{product_name} =~ s/^[^ ]* *//;
        $$pep{pep_begin}    = $$raw{query_left};
        $$pep{pep_end}      = $$raw{query_right};
        $$pep{pep_length}   = $$pep{pep_end} - $$pep{pep_begin} + 1;
        $$pep{pep_sequence} =
          substr( $polyprotein, $$pep{pep_begin} - 1, $$pep{pep_length} );
        $$pep{fuzzy_begin} = $$raw{fuzzy_begin};
        $$pep{fuzzy_end}   = $$raw{fuzzy_end};

        # eliminate small fragments of mature peptides
        if ( $$pep{pep_sequence} =~ /^X{3,}/i ) {
            $$pep{pep_sequence} =~ s/^X+//gi;
            my $newlen = length( $$pep{pep_sequence} );
            if ( $newlen < 10 ) { next }
            my $delta = $$pep{pep_length} - $newlen;
            $$pep{fuzzy_begin} = 1;
            $$pep{pep_begin} += $delta;
            $$pep{pep_length} -= $delta;
        }
        if ( $$pep{pep_sequence} =~ /X{3,}$/i ) {
            $$pep{pep_sequence} =~ s/X+$//gi;
            my $newlen = length( $$pep{pep_sequence} );
            if ( $newlen < 10 ) { next }
            my $delta = $$pep{pep_length} - $newlen;
            $$pep{fuzzy_end} = 1;
            $$pep{pep_end}    -= $delta;
            $$pep{pep_length} -= $delta;
        }
        if ( $$pep{fuzzy_begin} || $$pep{fuzzy_end} ) {
            if ( $$pep{pep_length} < 10 ) { next }
        }

        # compare to reference
        $$pep{ref_id} = $$raw{subject_id};
        if ( $$pep{ref_id} ne "?" ) {
            my $refstart = $$raw{subject_begin};
            $$pep{ref_begin} = maxval( $refstart, 1 );
            my $refend = $$raw{subject_end};
            $$pep{ref_end} = minval( $refend, $$raw{subject_length} );
            $$pep{ref_length} = $$raw{subject_length};
            my $refpep = $refmat{ $$pep{ref_id} }{sequence};
            score_pep_vs_reference( $$pep{pep_id}, $$pep{pep_sequence},
                $$pep{ref_id}, $refpep, $pep );
            if ( $$pep{num_reftrunc5} == 0 ) { $$pep{fuzzy_begin} = 0 }
            if ( $$pep{num_reftrunc3} == 0 ) { $$pep{fuzzy_end}   = 0 }
            #delete( $$pep{pep_alignment} );
            round_refstats($pep);
        }

        if ( $pepno == 1 && $$gene{start_truncation} && $$pep{ref_begin} > 1 ) {
            $$pep{fuzzy_begin} = 1;
        }
        elsif ($pepno == @rawmatpeps
            && $$gene{stop_truncation}
            && $$pep{ref_end} < $$pep{ref_length} )
        {
            $$pep{fuzzy_end} = 1;
        }

        if ( defined $min_gene_size && $$pep{pep_length} >= $min_gene_size ) {
        }
        elsif ( defined $min_gene_coverage
            && $$pep{pct_refcoverage} >= $min_gene_coverage )
        {
        }
        elsif ( defined $min_gene_size || defined $min_gene_coverage ) {
            next;
        }

        push @matpeps, $pep;
    }

    #    if ( ! $matpeps[0]{fuzzy_begin} ) {
    #        $matpeps[0]{fuzzy_begin} = $$gene{start_truncation};
    #    }
    #    if ( ! $matpeps[@matpeps-1]{fuzzy_end} ) {
    #        $matpeps[@matpeps-1]{fuzzy_end} = $$gene{stop_truncation};
    #    }
    return \@matpeps;
}


=comment Not used anywhere in Vigor package
sub matpep_grouping {
    my ( $a, $b ) = @_;

    return $$a{query_left} + $$a{query_right} <=> $$b{query_left} +
      $$b{query_right};
}
=cut

sub adjust_matpep_edges {
    my ( $polyprotein, $refmat, $mpleft, $mpright, %products ) = @_;
    my $dbg = DEBUG;

    if ( $$mpleft{query_right} + 1 == $$mpright{query_left} ) { return }

    if ($dbg) {
        print
"adjust edge between $$mpleft{query_left}-$$mpleft{query_right} and $$mpright{query_left}-$$mpright{query_right}\n";
    }

    my $polylen = length($polyprotein);
    if ( $polyprotein =~ /\*$/ ) { $polylen-- }

    my $lenl = 0;
    my $refl;
    my %profl;
    if ( $$mpleft{subject_id} ne "?" ) {
        $refl  = $$refmat{ $$mpleft{subject_id} }{sequence};
        %profl = profile_peptide($refl);
        $lenl  = length($refl);
    }
    my $lenr = 0;
    my $refr;
    my %profr;
    if ( $$mpright{subject_id} ne "?" ) {
        $refr  = $$refmat{ $$mpright{subject_id} }{sequence};
        %profr = profile_peptide($refr);
        $lenr  = length($refr);
    }
    if ( !$lenl && !$lenr ) {
        $$mpleft{query_right} =
          int( ( $$mpright{query_left} - 1 + $$mpleft{query_right} ) / 2.0 );
        $$mpright{query_left}  = $$mpleft{query_right} + 1;
        $$mpleft{fuzzy_end}    = 0;
        $$mpright{fuzzy_begin} = 0;
        if ($dbg) {
            print
"new edge $$mpleft{query_left}-$$mpleft{query_right} | $$mpright{query_left}-$$mpright{query_right}\n";
        }
        return;
    }

    my $bestedge;
    my $bestscore = 0;
    my $lowscore  = 999999999;
    my $begin = minval( $$mpleft{query_right}, $$mpright{query_left} - 1 ) - 1;
    my $end = maxval( $$mpleft{query_right}, $$mpright{query_left} - 1 ) + 1;
    if ( $$mpleft{subject_id} eq "?" ) { $begin -= 3 }
    if ( $$mpright{subject_id} eq "?" ) { $end += 3 }
    if ( $begin < 1 )      { $begin = 1 }
    if ( $end > $polylen ) { $end   = $polylen }

    if ($dbg) {
        print
"L: $$mpleft{subject_id}  R: $$mpright{subject_id}  segment: $$mpleft{query_left}-$$mpright{query_right}  edge=$$mpleft{query_right}|$$mpright{query_left} $begin-$end\n";
    }

    my $prodleft = get_defline_tag( $$mpleft{subject_definition}, "product" );
    for my $edge ( $begin .. $end ) {

        my $pepscoreL  = 0;
        my $edgescoreL = -1;
        if ($lenl) {
            my $mp = substr(
                $polyprotein,
                $$mpleft{query_left} - 1,
                $edge - $$mpleft{query_left} + 1
            );
            $pepscoreL = score_profile( $mp, %profl );
            my $product =
              get_defline_tag( $$mpleft{subject_definition}, "product" );
            if ( exists $products{$product}{right}{$edge} ) {
                $edgescoreL = 6.0 / $products{$product}{total} *
                  $products{$product}{right}{$edge};
            }
        }

        my $pepscoreR  = 0;
        my $edgescoreR = -1;
        if ($lenr) {
            my $mp =
              substr( $polyprotein, $edge, $$mpright{query_right} - $edge );
            $pepscoreR = score_profile( $mp, %profr );
            my $product =
              get_defline_tag( $$mpright{subject_definition}, "product" );
            if ( exists $products{$product}{left}{ $edge + 1 } ) {
                $edgescoreR = 6.0 / $products{$product}{total} *
                  $products{$product}{left}{ $edge + 1 };
            }
        }

        my $score = $pepscoreL + $edgescoreL + $pepscoreR + $edgescoreR;
        if ( $score > $bestscore ) {
            $bestscore = $score;
            $bestedge  = $edge;
        }
        elsif ( $score < $lowscore ) {
            $lowscore = $score;
        }
        if ($dbg) {
            print
"edge=$edge  L($lenl)=$pepscoreL+$edgescoreL  R($lenr)=$pepscoreR+$edgescoreR  score=$score  range=$lowscore-$bestscore\n";
        }
    }
    if ($dbg) {
        print "  bestedge=$bestedge  scorerange=$lowscore-$bestscore\n";
    }

    if ($dbg) { print "$$mpleft{subject_id} $$mpleft{query_right}\n" }
    if ($dbg) { print "$$mpright{subject_id} $$mpright{query_left}\n" }
    if ( $$mpleft{subject_id} ne "?" && $$mpright{subject_id} ne "?" ) {
        if ( $$mpleft{query_right} == $$mpright{query_left} - 1 ) {
            $$mpleft{fuzzy_end}    = 0;
            $$mpright{fuzzy_begin} = 0;
            if ($dbg) {
                print "return 2. $$mpleft{query_right}|$$mpright{query_left}\n";
            }
            return;
        }
    }
    if ( $$mpleft{subject_id} eq "?" ) {
        $$mpleft{query_right} = $$mpright{query_left} - 1;
        if ($dbg) {
            print "return 6. $$mpleft{query_right}|$$mpright{query_left}\n";
        }
    }
    elsif ( $$mpright{subject_id} eq "?" ) {
        $$mpright{query_left} = $$mpleft{query_right} + 1;
        if ($dbg) {
            print "return 7. $$mpleft{query_right}|$$mpright{query_left}\n";
        }
    }
    else {
        $$mpleft{query_right} = $bestedge;
        $$mpright{query_left} = $bestedge + 1;
        if ($dbg) {
            print "return 8. $$mpleft{query_right}|$$mpright{query_left}\n";
        }
    }
    $$mpleft{fuzzy_end}    = 0;
    $$mpright{fuzzy_begin} = 0;

    return;
}

sub remove_redundant_matpeps {
    my ($genes) = @_;

    for my $i ( 0 .. @$genes - 1 ) {
        my $big = $$genes[$i];
        if ( $$big{is_pseudogene} )        { next }
        if ( !defined $$big{mature_peps} ) { next }

        for my $j ( 0 .. @$genes - 1 ) {
            if ( $j == $i ) { next }
            my $small = $$genes[$j];
            if ( $$small{is_pseudogene} )                          { next }
            if ( !defined $$small{mature_peps} )                   { next }
            if ( $$small{orientation} != $$big{orientation} )      { next }
            if ( $$small{start_site} != $$big{start_site} )        { next }
            if ( $$small{protein_length} > $$big{protein_length} ) { next }
            if (   $$small{protein_length} == $$big{protein_length}
                && $$small{gene_name} gt $$big{gene_name} )
            {
                next;
            }

            my @pepSmall = @{ $$small{mature_peps} };
            my @pepBig   = @{ $$big{mature_peps} };
            for my $a ( 0 .. @pepSmall - 1 ) {
                for my $b ( 0 .. @pepBig - 1 ) {
                    if (   $pepSmall[$a]{pep_begin} == $pepBig[$b]{pep_begin}
                        && $pepSmall[$a]{pep_end} == $pepBig[$b]{pep_end} )
                    {
                        $pepSmall[$a] = undef;
                        last;
                    }
                }
            }

            my @tmp = remove_undefs(@pepSmall);
            $$small{mature_peps} = \@tmp;
        }
    }
}

sub cleanup_genes {
    my ( $RPT, $genome, $genes ) = @_;

    my %names;
    my %groups;
    my @dups;
    my $min_pseudogene_identity   = get_parameter("min_pseudogene_identity");
    my $min_pseudogene_similarity = get_parameter("min_pseudogene_similarity");
    my $min_pseudogene_coverage   = get_parameter("min_pseudogene_coverage");
    my @keep;

    # remove mutually exclusive genes
    @$genes = sort { compare_gene_positions( $a, $b ) } remove_undefs(@$genes);
    my %groupid;
    for my $g ( 0 .. @$genes - 1 ) {
        my $gene = $$genes[$g];
        if ( !defined $gene ) {
            next;
        }
        my %exclusions = get_reference_exclusions( $$gene{ref_id} );
        
        if (%exclusions) {
            for my $x ( 0 .. @$genes - 1 ) {
                my $excluded = $$genes[$x];
                
                if ( !defined $excluded ) { 
                    next;
                }
                if ( exists $exclusions{ $$excluded{gene_name} } ) {
                    if ( $$gene{is_pseudogene} <= $$excluded{is_pseudogene} ) {
                        $$genes[$x] = undef;
                    }
                    else {
                        $$genes[$g] = undef;
                        last;
                    }
                }
            }
        }
    }
    @$genes = remove_undefs(@$genes);

    # remove failed genes that do not qualify as pseudogenes
    my $geneno = 0;
    my %groupnum;
    my %genecount;
    for my $gene (@$genes) {
        $genecount{ $$gene{gene_id} }++;
    }
    for my $gene (@$genes) {
        if (   defined($$gene{frameshifts}) && @{ $$gene{frameshifts} } 
            || defined ($$gene{embeddedstops}) && @{ $$gene{embeddedstops} } 
            || defined ($$gene{cdserrors})     && @{ $$gene{cdserrors} }  )
        {
            if ( is_optional( $$gene{ref_id} ) ) {
                print "\ndropping pseudogene $$gene{gene_id} $$gene{gene_name} $$gene{product_name}: gene is optional\n";
                print_genehits($gene);
                print $RPT "\ndropping pseudogene $$gene{gene_id} $$gene{gene_name} $$gene{product_name}: gene is optional\n";
                my $STDOUT = *STDOUT;
                *STDOUT = $RPT;
                print_genehits($gene);
                *STDOUT = $STDOUT;

                $genecount{ $$gene{gene_id} }--;
                next;
            }
            else {
                my @fails;

#print "min  id $min_pseudogene_identity  sim $min_pseudogene_similarity  cov $min_pseudogene_coverage\n";
#print "gene id $$gene{pct_refidentity}  sim $$gene{pct_refsimilarity}  cov $$gene{pct_refcoverage}\n";
                if ( $$gene{pct_refidentity} < $min_pseudogene_identity ) {
                    push(@fails, "%id $$gene{pct_refidentity} < $min_pseudogene_identity");
                }
                if ( $$gene{pct_refsimilarity} < $min_pseudogene_similarity ) {
                    push(@fails,"%sim $$gene{pct_refsimilarity} < $min_pseudogene_similarity");
                }
                if ( !check_coverage( $min_pseudogene_coverage, $gene, $genome )) {
                    push(@fails, "%cov $$gene{pct_refcoverage} < $min_pseudogene_coverage");
                }

                #                if (   @fails == 0
                #                    && defined $$gene{embeddedstops}
                #                    && @{ $$gene{embeddedstops} } > 3 )
                #                {
                #                    push @fails, "more than three stop codons";
                #                }
                if (@fails) {
                    if ( is_required( $$gene{ref_id} ) && $genecount{ $$gene{gene_id} } == 1 ) {
                        next;
                    }
                    my $message =
                        "\ndropping pseudogene $$gene{gene_id}, "
                      . "fails match requirements: "
                      . join( ", ", @fails ) . "\n";
                    print $message;
                    print_genehits($gene);
                    print $RPT $message;
                    my $STDOUT = *STDOUT;
                    *STDOUT = $RPT;
                    print_genehits($gene);
                    *STDOUT = $STDOUT;

                    $genecount{ $$gene{gene_id} }--;
                    next;
                }
            }
        }

        # adjust identifiers to remove gaps caused by dropped genes
        $geneno++;
        my $gene_id = $$gene{gene_id};
        my $postfix;
        if ( $gene_id =~ /([a-z])$/ ) {
            $postfix = $1;
        }
        $gene_id =~ s/\.[^.]+$//;
        if ( defined $postfix ) {
            my $origid = $$gene{gene_id};
            $origid =~ s/$postfix$//;
            if ( exists $groupnum{$origid} ) {
                $gene_id .= ".$groupnum{$origid}$postfix";
                $geneno--;
            }
            else {
                $groupnum{$origid} = $geneno;
                $gene_id .= ".$geneno$postfix";
            }
        }
        else {
            $gene_id .= ".$geneno";
        }

        if ( $gene_id ne $$gene{gene_id} ) {
            if ( defined $$gene{mature_peps} ) {
                my $pepno = 0;
                for my $matpep ( @{ $$gene{mature_peps} } ) {
                    $pepno++;
                    $$matpep{pep_id} = "$gene_id.$pepno";
                }
            }
            $$gene{gene_id} = $gene_id;
        }
        push @keep, $gene;
    }

    # find genes with duplicated names
    for my $gene ( sort { compare_gene_positions( $a, $b, 1 ) } @keep ) {
        if ( $$gene{gene_id} !~ /[b-z]$/ ) {
            if ( exists $names{ $$gene{gene_name} } ) {
                push @dups, $gene;
            }
            else {
                $names{ $$gene{gene_name} } = $gene;
            }
        }
        else {
            my $group = $$gene{gene_id};
            $group =~ s/[b-z]$//;
            if ( exists $groups{$group} ) {
                push @{ $groups{$group} }, $gene;
            }
            else {
                my @tmp = ($gene);
                $groups{$group} = \@tmp;
            }
        }
    }

    # rename genes with duplicated names
    for my $dup (@dups) {
        my $origname = $$dup{gene_name};
        my $origprod = get_reference_product( $$dup{ref_id} );
        my $modifier = 2;
        my $newname  = modify_gene_name( $$dup{ref_id}, $origname, $modifier );
        while ( exists $names{$newname} ) {
            $modifier++;
            $newname = modify_gene_name( $$dup{ref_id}, $origname, $modifier );
        }
        my $newproduct = get_reference_product( $$dup{ref_id}, undef, $modifier );
        if ( !defined $newproduct ) { $newproduct = $origprod }

        my $newnote;
        my $symbol = gene_name_to_symbol($origname);
        if ( !defined $symbol || index( $origprod, $symbol ) >= 0 ) {
            if ( $origprod ne "hypothetical protein" ) {
                $newnote = "duplication of $origprod";
            }
        }
        else {
            $newnote = "duplication of $symbol $origprod";
        }

        $$dup{gene_name}    = $newname;
        $$dup{note}         = $newnote;
        $names{$newname}    = $dup;

        my $group = $$dup{gene_id};
        $group =~ s/a$//;
        if ( exists $groups{$group} ) {
            for my $grp ( @{ $groups{$group} } ) {
                $$grp{gene_name}    = $newname;
                $$grp{note}         = $newnote;
            }
        }
    }
    
    # remove redundant mature peptides
    remove_redundant_matpeps( $genes );    

    # update gene list
    @$genes = sort { compare_gene_ids( $a, $b ) } @keep;
    
    ## Scanning all the genes stretching the 3' end of partial genes till the end of the molecule
    my $seg_ln = $genome->{seqlen};
    
    foreach my $gene (@{$genes}) {
        if (!$gene->{stop_truncation} || $seg_ln - $gene->{stop_site} >= CODON_LN) {
            next;
        }
        my $current_end = $gene->{stop_site};
        $gene->{stop_site} = $seg_ln;
        
        foreach my $exon (@{$gene->{exons}}) { 
            if ($exon->{dna_end} == $current_end) {
                $exon->{dna_end} = $seg_ln;
            }
=comment For now we ignore these two coordinates. One need to streamline the gene data structure removing elements that are no longer used 
            if ($exon->{cdna_end} == $current_end) {
                $exon->{cdna_end} = $seg_ln;
            }
            if ($exon->{orig_dna_end} == $current_end) {
                $exon->{orig_dna_end} = $seg_ln;
            }
=cut            
        }
    }
}
sub modify_gene_name {
    my ( $ref_id, $name, $modifier ) = @_;

    my $newname = get_reference_name( $ref_id, undef, $modifier );

    #print "modified $name with $modifier producing $newname\n";
    if ( defined $newname ) { return $newname }

    $newname = $name;
    my $roman = $roman_numerals[ $modifier - 1 ];
    if ( $newname =~ /-/ ) {
        $newname =~ s/-/-$roman-/;
    }
    else {
        $newname .= "-$roman";
    }

    #print "modified $name with $modifier producing $newname\n";
    return $newname;
}

####################################
# pairwise alignments
# align peptide to reference
sub align_pep_to_reference {
    my ( $pepid, $pepseq, $refid, $refseq ) = @_;

    my $vigorspace = get_parameter("vigorspace");

    my $tmppepseq = $pepseq;
    $tmppepseq =~ s/\*$//;
    if ( $tmppepseq =~ /\*/ ) {

        #print "seq $tmppepseq\n";
        $tmppepseq =~ s/\*/x/g;

        #print " => $tmppepseq\n";
    }

    open( TMP, ">$vigorspace/tmp_align_fasta" )
      || die "Could not write to the tmp_align_fasta file\.\n";
    print TMP ">$pepid\n$tmppepseq\n>$refid\n$refseq\n";

#if ( getlogin() eq "jhoover" ) { print ">$tmppepid\n$tmppepseq\n>$tmprefid\n$refseq\n" }
    close TMP;

    my $cmd = "$myBin/clustalw -align -GAPOPEN=6.0 -GAPEXT=0.3 -infile=$vigorspace/tmp_align_fasta &> $vigorspace/tmp_align_fasta.log";

#"muscle -in $vigorspace/tmp_align_fasta -out $vigorspace/tmp_align_fasta.aln -clwstrict &> $vigorspace/tmp_align_fasta.log";
    &runCmd($vigorspace, $cmd);

    my $alignment = &runCmdAndGetResults($vigorspace, "cat tmp_align_fasta.aln");
    $alignment =~ s/[\n\r\s]+$//;

    #print $alignment . "\n" . `cat $vigorspace/tmp_align_fasta.log` . "\n";

    #unlink "$vigorspace/tmp_align_fasta";
    #unlink "$vigorspace/tmp_align_fasta.aln";
    #unlink "$vigorspace/tmp_align_fasta.log";

    return $alignment;
}

# align gene span to reference
sub bl2seq_alignment {
    my ( $frag, $gene, $bl2seq_blastopts ) = @_;
    my $dbg = DEBUG;

    #if ( `whoami` =~ /jhoover/ ) { $dbg = 1 }

    my $genome = $frag;
    if ( defined $$frag{full_genome} ) { $genome = $$frag{full_genome} }
    my $strand = $$gene{orientation};

    my $genespan = bl2seq_makespan( $frag, $genome, $gene );
    my $refaa = get_reference_seq( $$gene{ref_id} );
    if ( !defined $bl2seq_blastopts ) {
        $bl2seq_blastopts = bl2seq_blastopts( $$refaa{seqlen} );
    }
    my $output = bl2seq_runner( $genespan, $refaa, $bl2seq_blastopts );

    # calculate offset adjustment for the coordinates
    my $offset = $$genespan{left} - 1;
    if ( $strand == -1 ) { $offset = $$genespan{right} + 1 }

# calculate the largest possible coordinate so we can size the margin appropriately
    my $coordsz = $$gene{start_site} + $$frag{frag_offset};
    if ( $$gene{stop_site} + $$frag{frag_offset} > $coordsz ) {
        $coordsz = $$gene{stop_site} + $$frag{frag_offset};
    }
    if ( length( $reference_seqs{ $$gene{ref_id} }{sequence} ) > $coordsz ) {
        $coordsz = length( $reference_seqs{ $$gene{ref_id} }{sequence} );
    }
    $coordsz = length($coordsz);

    # parse bl2seq output
    my @hsps =
      bl2seq_parser( $output, $coordsz, $offset, $$genome{seqlen}, $strand );

    # sort HSPs by position instead of score
    @hsps = sort { $$a{subject_begin} <=> $$b{subject_begin} } @hsps;

    # build report
    my $alignment = "$$genespan{fasta}\n";
    $alignment .= "Length = $$refaa{seqlen}\n";

    for my $i ( 0 .. @hsps - 1 ) {
        if ( !defined $hsps[$i]{body} ) { next }
        if (   $i > 0
            && $hsps[$i]{query_frame} ne $hsps[ $i - 1 ]{query_frame}
            && sign( $hsps[$i]{query_frame} ) ==
            sign( $hsps[ $i - 1 ]{query_frame} ) )
        {
            $alignment .=
              bl2seq_frameshift( $genome, $strand, $refaa, $hsps[ $i - 1 ],
                $hsps[$i] );
        }
        $alignment .= "\n$hsps[$i]{body}";
    }

    #print "bl2seq_alignment out\n";
    return $alignment;
}

sub bl2seq_makespan {
    my ( $frag, $genome, $gene ) = @_;

    #print "bl2seq_makespan in\n";

#print "gene " . ( $$gene{start_site} + $$frag{frag_offset} ) . "-" . ( $$gene{stop_site} + $$frag{frag_offset} ) . "\n";
    my $extend_span    = 150;
    my $allow_splicing = allow_splicing( $$gene{ref_id} );
    my $gene_span      = "";

    my $first =
      minval( $$gene{start_site}, $$gene{stop_site} ) + $$frag{frag_offset};
    my $left = $first - $extend_span;
    if ( $left < 1 ) {
        $left = 1;
    }

    #print "left=$left\n";

    my $right =
      $left + $extend_span + abs( $$gene{stop_site} - $$gene{start_site} ) +
      $extend_span;

    #print "right=$right   max=$$genome{seqlen}\n";
    if ( $right > $$genome{seqlen} ) {
        $left -= $right - $$genome{seqlen};
        if ( $left < 1 ) { $left = 1 }
        $right = $$genome{seqlen};
    }

    #print "left-right $left-$right\n";

    my $end = $left - 1;
    for
      my $exon ( sort { $$a{dna_begin} <=> $$b{dna_begin} } @{ $$gene{exons} } )
    {

#print "exon " . ( $$exon{dna_begin} + $$frag{frag_offset} ) . "-" . ( $$exon{dna_end} + $$frag{frag_offset} ) . "\n";
        my $begin =
          minval( $$exon{dna_begin}, $$exon{dna_end} ) + $$frag{frag_offset};

        #print "end=$end  begin=$begin\n";
        if ( $begin > $end + 1 ) {

            #print "intron " . ( $end+1 ) . "-" . ( $begin-1) . "\n";
            my $intron =
              lc subsequence( $$genome{sequence}, $end + 1, $begin - 1 );
            $gene_span .= $intron;

            #print "span=$gene_span\n";
            $end = $begin - 1;
        }
        elsif ( $begin <= $end ) {
            $begin = $end + 1;
        }
        $end =
          maxval( $$exon{dna_begin}, $$exon{dna_end} ) + $$frag{frag_offset};
        if ( $end >= $begin ) {
            $gene_span .= uc subsequence( $$genome{sequence}, $begin, $end );
        }

        #print "exon $begin-$end\n";
        #print "span=$gene_span\n";
    }

    if ( $right > $end ) {
        $gene_span .= lc subsequence( $$genome{sequence}, $end + 1, $right );

        #print "intron " . ( $end+1 ) . "-" . ( $right ) . "\n";
        #print "span=$gene_span\n";
    }

    if ( $$gene{orientation} == -1 ) {
        $gene_span = reverse_complement($gene_span);
    }

    $gene_span =~ s/(.{1,60})/$1\n/g;
    
    my $startcodon = $$gene{start_site} + $$frag{frag_offset};
    my $stopsite   = $$gene{stop_site} + $$frag{frag_offset};

    my $span;
    $$span{left}  = $left;
    $$span{right} = $right;
    $$span{fasta} = ">$$gene{gene_id}\t$startcodon\t$stopsite\n$gene_span";

    #print "bl2seq_makespan out\n";
    return $span;
}

sub bl2seq_runner {

    #print "bl2seq_runner in\n";
    my ( $span, $ref, $bl2seq_blastopts ) = @_;
    my $vigorspace = get_parameter("vigorspace");

    # use bl2seq to align extended gene span to reference sequence
    unlink "$vigorspace/mutantspan.fasta";
    open( TMP, ">$vigorspace/mutantspan.fasta" )
      || die "could not write mutantspan.fasta to disk";
    print TMP $$span{fasta};
    close TMP;

    unlink "$vigorspace/mutantref.fasta";
    my $ref_sequence = $$ref{sequence};
    open( TMP, ">$vigorspace/mutantref.fasta" )
      || die "could not write mutantref.fasta to disk";
    print TMP ">$$ref{defline}\n$ref_sequence";
    close TMP;

    my $cmd = "$myBin/bl2seq $bl2seq_blastopts -i $vigorspace/mutantspan.fasta -j $vigorspace/mutantref.fasta";
    my $output = &runCmdAndGetResults($vigorspace, $cmd);
   return $output;
}

sub bl2seq_blastopts {
    my ($protsz) = @_;

    my $rawopts = get_parameter("candidate_blastopts");
    my @opts    = split /  */, $rawopts;

    my $bl2seqopts = "";

    #    my $evalue = 0.001;
    #    while (@opts) {
    #        my $opt = shift @opts;
    #        if ( $opt eq "-e" ) {
    #            $bl2seqopts .= "-e $evalue";
    #        }
    #        elsif ($opt eq "-p" || $opt eq "-M" )
    #        {
    #            $bl2seqopts .= $opt . " " . ( shift @opts ) . " ";
    #        }
    #    }
    $bl2seqopts .= " -p blastx -g F -e 0.001 -F \"\"";

    #print "   raw $rawopts\nbl2seq $bl2seqopts\n";
    return $bl2seqopts;
}

sub candidate_blastopts {
    my ($fragsz) = @_;

    my $expected_size    = 5000;
    my $rawopts          = get_parameter("candidate_blastopts");
    my $candidate_evalue = get_parameter("candidate_evalue");

    my @opts = split /  */, $rawopts;

    my $evalue;
    my $dbsize;
    for my $i ( 0 .. @opts - 1 ) {
        if ( $opts[$i] eq "-e" ) {
            $evalue = $opts[ $i + 1 ];
        }
    }
    if ( !defined $evalue ) {
        $evalue = 10;
        push @opts, "-e", $evalue;
    }
    if ( $candidate_evalue ne "" ) {
        $evalue = $candidate_evalue;
    }
    else {
        my $scale = $expected_size / $fragsz;
        if ( $scale < 1.0 ) { $scale = 1.0 }
        $evalue = $evalue * $scale * $scale;
    }

    my (@tmp) = split /[eE]/, $evalue;

    #print "  fragsz=$fragsz new value=$evalue";
    if ( @tmp == 1 ) {
        $evalue = int( 100000.0 * $evalue + 0.5 ) / 100000.0;
    }
    elsif ( @tmp < 3 ) {
        $tmp[0] = int( 10.0 * $tmp[0] + 0.5 ) / 10.0;
        $evalue = join( "E", @tmp );
    }

    #print "  reformatted=$evalue\n";

    #print "$evalue\n";

    my $blastopts = "";
    while (@opts) {
        my $opt = shift @opts;
        if ( $opt eq "-e" ) {
            $blastopts .= "-e $evalue ";
            shift @opts;
        }
        else {
            $blastopts .= $opt . " " . ( shift @opts ) . " ";
        }
    }

    #print "   raw $rawopts\nfragsz $fragsz\nfrag $blastopts\n";
    return $blastopts;
}

sub bl2seq_parser {
    my ( $input, $coordsz, $offset, $genomelen, $strand ) = @_;
    #print "bl2seq_parser $subject_id in\n";

    #print "INPUT=\n$input\n";
    # break input into HSP blocks
    my @blks = split /Score =/, $input;
    if ( !@blks ) { return undef }

    # trim header and tail
    my $hdr  = shift @blks;
    my $tail = pop @blks;
    if ( !defined $tail ) { return undef }
    my @tmp = split /Lambda     K      H/, $tail;
    push @blks, $tmp[0];

    # parse HSPs and
    # adjust coordinates from sub-span aligned to full genome
    my $oldmargin;
    my $newmargin;
    my @hsps;
    for my $blk (@blks) {
        $blk = "Score = $blk";

        #print "BLK=\n$blk\n";
        my $hsp;
        my @body;
        my $blank = 1;
        for my $line ( split /\n/, $blk ) {

            # skip consecutive blank lines
            if ( $line =~ /^[ \t]*$/ ) {
                if ($blank) { next }
                $blank = 1;
                push @body, $line;
                next;
            }
            $blank = 0;

            # query (genome) string
            if ( $line =~ /(Query: +([0-9]+) +)([^ ]+) +([0-9]+) *$/ ) {
                my ( $prefix, $begin, $seq, $end ) = ( $1, $2, $3, $4 );
                $oldmargin = length($prefix);

                # adjust coordinates
                if ( $strand == 1 ) {
                    $begin = $offset + $begin;
                    $end   = $offset + $end;
                }
                else {
                    $begin = $offset - $begin;
                    $end   = $offset - $end;
                }

                # reformat line
                $prefix =
                  "Query: " . substr( "$begin            ", 0, $coordsz ) . " ";
                if ( !defined $newmargin ) { $newmargin = length($prefix) }
                $line = $prefix . $seq . " " . $end;

                # update HSP ranges
                if ( !defined $$hsp{query_begin} ) {
                    $$hsp{query_begin} = $begin;
                }
                $$hsp{query_end}  = $end;
                $$hsp{query_left} =
                  minval( $$hsp{query_begin}, $$hsp{query_end} );
                $$hsp{query_right} =
                  maxval( $$hsp{query_begin}, $$hsp{query_end} );
                $$hsp{query_coverage} =
                  $$hsp{query_right} - $$hsp{query_left} + 1;
            }

            # subject (reference protein) string
            elsif ( $line =~ /(Sbjct: +([0-9]+) +)([^ ]+) +([0-9]+) *$/ ) {
                my ( $prefix, $begin, $seq, $end ) = ( $1, $2, $3, $4 );
                $oldmargin = length($prefix);

                # reformat line
                $prefix =
                  "Sbjct: " . substr( "$begin            ", 0, $coordsz ) . " ";
                if ( !defined $newmargin ) { $newmargin = length($prefix) }
                $line = $prefix . $seq . " " . $end;

                # update HSP ranges
                if ( !defined $$hsp{subject_begin} ) {
                    $$hsp{subject_begin} = $begin;
                }
                $$hsp{subject_end}   = $end;
                $$hsp{subject_begin} =
                  minval( $$hsp{subject_begin}, $$hsp{subject_end} );
                $$hsp{subject_end} =
                  maxval( $$hsp{subject_begin}, $$hsp{subject_end} );
                $$hsp{subject_coverage} =
                  $$hsp{subject_end} - $$hsp{subject_begin} + 1;
            }

            # score/identities/frame
            elsif ( $line =~ /^ *Score = *([^ ]*)/ ) {

                #print "line=$line score=$1\n";
                $$hsp{score} = $1;
                $line =~ s/^  *//;
            }
            elsif ( $line =~ /^ *Identities = / ) {
                $line =~ s/^  *//;
            }
            elsif ( $line =~ /^ * Frame = ([^ ]+) *$/ ) {
                $$hsp{query_frame} = $1;
                $line =~ s/^  *//;
            }

            # comparison string;
            else {
                if ( $oldmargin > $newmargin ) {
                    $line = substr( $line, $oldmargin - $newmargin );
                }
                elsif ( $oldmargin < $newmargin ) {
                    $line =
                      substr( "               ", 0, $newmargin - $oldmargin )
                      . $line;
                }
            }
            push @body, $line;
        }

        if ($blank) { pop @body }
        $$hsp{body} = join( "\n", @body ) . "\n";

        my $adjusted_frame;
        if ( $strand == 1 ) {
            $adjusted_frame = "+" . ( ( $$hsp{query_begin} - 1 ) % 3 + 1 );
        }
        else {
            $adjusted_frame =
              "-" . ( ( $genomelen - $$hsp{query_begin} ) % 3 + 1 );
        }
        if ( $$hsp{query_frame} ne $adjusted_frame ) {
            my $oldexp = "Frame = \\" . $$hsp{query_frame};
            my $newexp = "Frame = " . $adjusted_frame;
            $$hsp{body} =~ s/$oldexp/$newexp/;
            $$hsp{query_frame} = $adjusted_frame;
        }

        if ( $strand == 1 ) {
            if ( $$hsp{query_frame} < 0 ) { next }
        }
        else {
            if ( $$hsp{query_frame} =~ /\+/ ) { next }
        }
        push @hsps, $hsp;

        #print_hash( "HSP", $hsp );
    }

    # remove conflicting HSPs
    remove_conflicting_hsps( \@hsps, "query_left", "query_right", 0.5,
        "score" );
    remove_conflicting_hsps( \@hsps, "subject_begin", "subject_end", 0.5,
        "score" );

    #print "bl2seq_parser out\n";
    return @hsps;
}

sub remove_conflicting_hsps2 {
    my ( @rawhsps) = @_;
    my $dbg = DEBUG;
    #if ( $rawhsps[0]{subject_id} eq "NP_062883.2" ) { $dbg = 1 }

    if ( !@rawhsps ) { return @rawhsps }

    my @tmp = sort { $a->{subject_begin} <=> $b->{subject_begin} } @rawhsps;
    if ($dbg) {
        print "\nremove conflicting HSPs from\n";
        print_blasthits( 0, @tmp );
    }
    my $deletes = 0;

    # for each HSP
    for my $i ( 0 .. @tmp - 1 ) {
        my $ihsp = $tmp[$i];
        if ( !defined $ihsp ) { next }
        my $isize = $ihsp->{subject_begin} - $ihsp->{subject_end} + 1;

        # compare to other HSPs
        for my $j ( $i + 1 .. @tmp - 1 ) {
            my $jhsp = $tmp[$j];
            if ( !defined $jhsp ) { next }
            if ( $ihsp->{subject_id} ne $jhsp->{subject_id} ) { next }

            # is there a conflict between these HSPs?
            my $conflict = hsps_conflict( $ihsp, $jhsp );
            if ( $conflict == 0 ) { next }

            # HSPs conflict, discard lower scoring HSP
            if ( $ihsp->{evalue} < $jhsp->{evalue} ) {
                if ($dbg) {
                    print "    S1 keep/remove ($conflict)\n";
                    print_blasthits( 4, $ihsp, $jhsp );
                }
                $tmp[$j] = undef;
                $deletes++;
                next;
            }
            if ( $ihsp->{evalue} > $jhsp->{evalue} ) {
                if ($dbg) {
                    print "    S2 keep/remove ($conflict)\n";
                    print_blasthits( 4, $jhsp, $ihsp );
                }
                $tmp[$i] = undef;
                $deletes++;
                last;
            }
            if ( $ihsp->{vigor_matchwgt} >= $jhsp->{vigor_matchwgt} ) {
                if ($dbg) {
                    print "    S3 keep/remove ($conflict)\n";
                    print_blasthits( 4, $ihsp, $jhsp );
                }
                $tmp[$j] = undef;
                $deletes++;
                next;
            }
            elsif ( $ihsp->{vigor_matchwgt} < $jhsp->{vigor_matchwgt} ) {
                if ($dbg) {
                    print "    S4 keep/remove ($conflict)\n";
                    print_blasthits( 4, $jhsp, $ihsp );
                }
                $tmp[$i] = undef;
                $deletes++;
                last;
            }
            elsif ( $ihsp->{num_identical} >= $jhsp->{num_identical} ) {
                if ($dbg) {
                    print "    S5 keep/remove ($conflict)\n";
                    print_blasthits( 4, $ihsp, $jhsp );
                }
                $tmp[$j] = undef;
                $deletes++;
                next;
            }
            elsif ( $ihsp->{num_identical} < $jhsp->{num_identical} ) {
                if ($dbg) {
                    print "    S6 keep/remove ($conflict)\n";
                    print_blasthits( 4, $jhsp, $ihsp );
                }
                $tmp[$i] = undef;
                $deletes++;
                last;
            }
            else {
                if ($dbg) {
                    print "    S7 keep/remove ($conflict)\n";
                    print_blasthits( 4, $jhsp, $ihsp );
                }
                $tmp[$j] = undef;
                $deletes++;
                next;
            }
        }
    }

    # remove the discarded HSPs
    if ($deletes) {
        my @keep;
        for my $hsp (@tmp) {
            if ( defined $hsp ) { push @keep, $hsp }
        }
        if ($dbg) {
            print "  final HSPs (1)\n";
            print_blasthits( 2, @keep );
        }
        return @keep;
    }

    # return the surviving HSPs
    if ($dbg) {
        print "  final HSPs (2)\n";
        print_blasthits( 2, @rawhsps );
    }

    return @rawhsps;
}

sub hsps_conflict {
    my ( $ihsp, $jhsp ) = @_;
    my $dbg = DEBUG;
    #if ( $$ihsp{subject_id} eq "NP_062883.2" ) { $dbg = 1 }
    
    if ( $dbg ) {
        print "HSPs conflict?\n";
        print_blasthits( 0, $ihsp, $jhsp );
    }

    # inconsistent strand
    if ( $ihsp->{orientation} != $jhsp->{orientation} ) {
        return 1;
    }
    
    # inconsistent relative position
    my $ori = $ihsp->{orientation};
    if ( $ori == 1 ) {
        if ( $$ihsp{subject_end} > $$jhsp{subject_end} ) {
            if ( $$ihsp{query_right} <= $$jhsp{query_right} ) { return 2 }
        }
        elsif ( $$ihsp{subject_end} < $$jhsp{subject_end} ) {
            if ( $$ihsp{query_right} >= $$jhsp{query_right} ) { return 3 }
        }
        if ( $$ihsp{subject_begin} > $$jhsp{subject_begin} ) {
            if ( $$ihsp{query_left} <= $$jhsp{query_left} ) { return 4 }
        }
        elsif ( $$ihsp{subject_begin} < $$jhsp{subject_begin} ) {
            if ( $$ihsp{query_left} >= $$jhsp{query_left} ) { return 5 }
        }

        if ( $$ihsp{query_right} > $$jhsp{query_right} ) {
            if ( $$ihsp{subject_end} <= $$jhsp{subject_end} ) { return 6 }
        }
        elsif ( $$ihsp{query_right} < $$jhsp{query_right} ) {
            if ( $$ihsp{subject_end} >= $$jhsp{subject_end} ) { return 7 }
        }
        if ( $$ihsp{query_left} > $$jhsp{query_left} ) {
            if ( $$ihsp{subject_begin} <= $$jhsp{subject_begin} ) { return 8 }
        }
        if ( $$ihsp{query_left} < $$jhsp{query_left} ) {
            if ( $$ihsp{subject_begin} >= $$jhsp{subject_begin} ) { return 9 }
        }
    }
    else {
        if ( $$ihsp{subject_end} > $$jhsp{subject_end} ) {
            if ( $$ihsp{query_left} >= $$jhsp{query_left} ) { return 2 }
        }
        elsif ( $$ihsp{subject_end} < $$jhsp{subject_end} ) {
            if ( $$ihsp{query_left} <= $$jhsp{query_left} ) { return 3 }
        }
        if ( $$ihsp{subject_begin} > $$jhsp{subject_begin} ) {
            if ( $$ihsp{query_right} >= $$jhsp{query_right} ) { return 4 }
        }
        elsif ( $$ihsp{subject_begin} < $$jhsp{subject_begin} ) {
            if ( $$ihsp{query_right} <= $$jhsp{query_right} ) { return 5 }
        }

        if ( $$ihsp{query_left} > $$jhsp{query_left} ) {
            if ( $$ihsp{subject_end} >= $$jhsp{subject_end} ) { return 6 }
        }
        elsif ( $$ihsp{query_left} < $$jhsp{query_left} ) {
            if ( $$ihsp{subject_end} <= $$jhsp{subject_end} ) { return 7 }
        }
        if ( $$ihsp{query_right} > $$jhsp{query_right} ) {
            if ( $$ihsp{subject_begin} >= $$jhsp{subject_begin} ) { return 8 }
        }
        if ( $$ihsp{query_right} < $$jhsp{query_right} ) {
            if ( $$ihsp{subject_begin} <= $$jhsp{subject_begin} ) { return 9 }
        }
    }
#    my @iedges;
#    push @iedges, { dna => $ori == 1 ? $$ihsp{query_left} : $$ihsp{query_right}, sbj => $$ihsp{subject_begin} };
#    push @iedges, { dna => $ori == 1 ? $$ihsp{query_right} : $$ihsp{query_left}, sbj => $$ihsp{subject_end} };
#    my @jedges;
#    push @jedges, { dna => $ori == 1 ? $$jhsp{query_left} : $$jhsp{query_right}, sbj => $$jhsp{subject_begin} };
#    push @jedges, { dna => $ori == 1 ? $$jhsp{query_right} : $$jhsp{query_left}, sbj => $$jhsp{subject_end} };
#    my $case = 1;
#    for my $iedge ( @iedges ) {
#        for my $jedge ( @jedges ) {
#            if ( $$iedge{sbj} < $$jedge{sbj} ) {
#                if ( $ori * $$iedge{dna} >= $ori * $$jedge{dna} ) { return $case+1 }
#            }
#            elsif ( $$iedge{sbj} > $$jedge{sbj} ) {
#                if ( $ori * $$iedge{dna} <= $ori * $$jedge{dna} ) { return $case+2 }
#            }
#            if ( $ori * $$iedge{dna} < $ori * $$jedge{dna} ) {
#                if ( $$iedge{sbj} >= $$jedge{sbj} ) { return $case+3 }
#            }
#            elsif ( $ori * $$iedge{dna} > $ori * $$jedge{dna} ) {
#                if ( $$iedge{sbj} <= $$jedge{sbj} ) { return $case+4 }
#            }
#            $case += 4;
#        }
#    }
    
    # significant overlap on subject sequence
    {
        my $isize = $ihsp->{subject_end} - $ihsp->{subject_begin} + 1;
        my $jsize = $jhsp->{subject_end} - $jhsp->{subject_begin} + 1;
        my $oleft =
            $ihsp->{subject_begin} > $jhsp->{subject_begin}
          ? $ihsp->{subject_begin}
          : $jhsp->{subject_begin};
        my $oright =
            $ihsp->{subject_end} < $jhsp->{subject_end}
          ? $ihsp->{subject_end}
          : $jhsp->{subject_end};
        my $osize = $oright - $oleft + 1;

        #        if ( $osize / $ihsp->{subject_unit} >= 30.0 ) {
        #            return 8;
        #        }
        #        els
        if ( $osize > 0.6 * $isize || $osize > 0.6 * $jsize ) {
            return 10;
        }
        elsif ( $osize > 5.0
            && ( $osize >= 0.4 * $isize || $osize >= 0.4 * $jsize ) )
        {
            return 11;
        }
    }

    # significant overlap on query sequence
    {
        my $isize = $ihsp->{query_right} - $ihsp->{query_left} + 1;
        my $jsize = $jhsp->{query_right} - $jhsp->{query_left} + 1;
        my $oleft =
            $ihsp->{query_left} > $jhsp->{query_left}
          ? $ihsp->{query_left}
          : $jhsp->{query_left};
        my $oright =
            $ihsp->{query_right} < $jhsp->{query_right}
          ? $ihsp->{query_right}
          : $jhsp->{query_right};
        my $osize = $oright - $oleft + 1;
        if ( $osize > 0.6 * $isize || $osize > 0.6 * $jsize ) {
            return 12;
        }
        elsif ( $osize / 3.0 > 5.0
            && ( $osize >= 0.4 * $isize || $osize >= 0.4 * $jsize ) )
        {
            return 13;
        }
    }

    # HSPs are consistent with one-another
    return 0;
}

sub remove_conflicting_hsps {
    my ( $hsps, $left, $right, $fraction, $score ) = @_;
    my $dbg = DEBUG;
    #if ( get_reference_name( $$hsps[0]{subject_id} ) eq "NS3c" ) { $dbg = 1 }

    if ($dbg) {
        print
          "----------------------------------------------------------------\n";
        print "DECONFLICT $left, $right, $fraction, $score\n";
        print_blasthits( 0, @$hsps );
    }

    for my $i ( 0 .. @$hsps - 2 ) {
        if ( !defined $$hsps[$i] ) { next }

        #print_blasthits( 4, $$hsps[$i] );
        for my $j ( $i + 1 .. @$hsps - 1 ) {
            if ( !defined $$hsps[$j] ) { next }

            #print_blasthits( 4, $$hsps[$j] );
            #            if ( $$hsps[$j]{subject_id} ne $$hsps[$i]{subject_id} ) { next }
            my $overlap =
              minval( $$hsps[$i]{$right}, $$hsps[$j]{$right} ) -
              maxval( $$hsps[$i]{$left}, $$hsps[$j]{$left} ) + 1;
            my $allowed = $fraction * minval(
                $$hsps[$i]{$right} - $$hsps[$i]{$left} + 1,
                $$hsps[$j]{$right} - $$hsps[$j]{$left} + 1
            );
            if ($dbg) {
                print "  overlap $overlap  allowed $allowed\n";
                print_blasthits( 2, $$hsps[$i], $$hsps[$j] );
            }
            if ( $overlap > $allowed ) {
                if ( $$hsps[$i]{$score} >= $$hsps[$j]{$score} ) {
                    $$hsps[$j] = undef;
                }
                else {
                    $$hsps[$i] = undef;
                    last;
                }
            }
        }
    }
    if ($dbg) {

    }
    @$hsps = remove_undefs(@$hsps);
    if ($dbg) {
        print "RESULT\n";
        print_blasthits( 0, @$hsps );
        print
          "----------------------------------------------------------------\n";
    }
}

sub bl2seq_frameshift {

    #print "bl2seq_frameshift in\n";
    my ( $genome, $strand, $ref, $hsp1, $hsp2 ) = @_;

    #print "hsp1\n$$hsp1{body}\nhsp2\n$$hsp2{body}\n";

    if ( allow_splicing( $$ref{id} ) ) {
        my $qleft  = minval( $$hsp2{query_begin},   $$hsp1{query_end} );
        my $qright = maxval( $$hsp2{query_begin},   $$hsp1{query_end} );
        my $sleft  = minval( $$hsp2{subject_begin}, $$hsp1{subject_end} );
        my $sright = maxval( $$hsp2{subject_begin}, $$hsp1{subject_end} );
        if ( $qright - $qleft + 1 >= 3 * ( $sright - $sleft + 1 ) + 50 ) {
            return "";
        }
    }

    my $bufsz  = 80;
    my $qbegin = $$hsp1{query_end} + $strand - 3 * $bufsz * $strand;
    my $sbegin = $$hsp1{subject_end} + 1 - $bufsz;
    my $qend   = $$hsp2{query_begin} - $strand + 3 * $bufsz * $strand;
    my $send   = $$hsp2{subject_begin} - 1 + $bufsz;

    #print "bl2seq_frameshift 1\n";
    while ( $qbegin < 1 ) {
        $qbegin += 3;
        $sbegin += $strand;
    }

    #print "bl2seq_frameshift 2\n";
    while ( $qbegin > $$genome{seqlen} ) {
        $qbegin -= 3;
        $sbegin -= $strand;
    }

    #print "bl2seq_frameshift 3\n";
    while ( $qend < 1 ) {
        $qend += 3;
        $send += $strand;
    }

    #print "bl2seq_frameshift 4\n";
    while ( $qend > $$genome{seqlen} ) {
        $qend -= 3;
        $send -= $strand;
    }

    #print "bl2seq_frameshift 5\n";
    while ( $strand * ( $qend - $qbegin ) + 1 < 18 ) {
        my $nextbegin    = $qbegin - 3 * $strand;
        my $infiniteloop = 1;
        if ( $nextbegin >= 1 && $nextbegin <= $$genome{seqlen} ) {
            $qbegin = $nextbegin;
            $sbegin--;
            $infiniteloop = 0;
        }
        my $nextend = $qend - 3 * $strand;
        if ( $nextend >= 1 && $nextend <= $$genome{seqlen} ) {
            $qend         = $nextend;
            $infiniteloop = 0;
            $send++;
        }
        if ($infiniteloop) { last }
    }
    if ( $sbegin < 1 )           { $sbegin = 1 }
    if ( $send > $$ref{seqlen} ) { $send   = $$ref{seqlen} }

    my $qryseq = subsequence( $$genome{sequence}, $qbegin, $qend );
    my $sbjseq = substr( $$ref{sequence}, $sbegin, $send - $sbegin + 1 );

    my @frameshifts;
    my $frame = eval $$hsp1{query_frame};

    #print "bl2seq_frameshift 6\n";
    for my $f ( 0 .. 2 ) {
        if ( $f >= length($qryseq) ) { last }
        my $aa = DNA2AA( substr( $qryseq, $f ) );
        $aa = format_frameaa( $sbjseq, $aa );

        my $label = " frame +$frame";
        if ( $strand == -1 ) {
            $label = " frame $frame";
        }

        if ( $f > 0 ) {

            #            $label = substr( "  ", 0, $f ) . $label;
            $aa = substr( "  ", 0, $f ) . $aa;
        }

        $aa =~ s/ +$//;
        $aa =~ s/(.{1,60})/$1\n/g;
        
        my $frameshift;
        $$frameshift{id}    = $f;
        $$frameshift{frame} = $label;
        $$frameshift{aa}    = $aa;

        push @frameshifts, $frameshift;

        $frame += $strand;
        if ( $frame == 4 ) {
            $frame = 1;
        }
        elsif ( $frame == -4 ) {
            $frame = -1;
        }
    }

    $qryseq =~ s/(.{1,60})/$1\n/g;
    
    my $ori    = $qend < $qbegin ? -1 : 1;
    my $qb     = $qbegin;
    my @qrytmp = split /\n/, $qryseq;
    for my $i ( 0 .. @qrytmp - 1 ) {
        my $qe = $qb + $ori * 59;
        if ( $ori * $qe > $ori * $qend ) { $qe = $qend }
        $qrytmp[$i] =
          lpad( $qb, 6 ) . " " . rpad( $qrytmp[$i], 60 ) . " " . rpad( $qe, 6 );
        $qb += $ori * 60;
    }
    $qryseq = join( "\n", @qrytmp );

    $sbjseq =~ s/(.{1,60})/$1\n/g;
    
    my $report;
    if ( allow_splicing( $$ref{id} ) ) {
        $report =
            "\n------------------------------------------------------------\n"
          . "frameshifted or splicing from $$hsp1{query_frame} to $$hsp2{query_frame} between $qbegin-$qend\n"
          . "$qryseq\n";
    }
    else {
        $report =
            "\n------------------------------------------------------------\n"
          . "frameshifted from $$hsp1{query_frame} to $$hsp2{query_frame} between $qbegin-$qend\n"
          . "$qryseq\n";
    }
    for my $frameshift ( sort { $$a{id} <=> $$b{id} } @frameshifts ) {
        my $qb = $qbegin;
        my @frmtmp = split /\n/, $$frameshift{aa};
        for my $i ( 0 .. @frmtmp - 1 ) {
            my $qe = $qb + $ori * 59;
            if ( $ori * $qe > $ori * $qend ) { $qe = $qend }
            my $hilite = "";
            if ( $frmtmp[$i] =~ /-[A-Z]--[A-Z]--[A-Z]-/ ) { $hilite = " |" }
            $frmtmp[$i] =
                lpad( $qb, 6 ) . " "
              . rpad( $frmtmp[$i], 60 ) . " "
              . rpad( $qe, 6 )
              . $hilite;
            $qb += $ori * 60;
        }
        my $frmseq = join( "\n", @frmtmp );

        $report .= "$$frameshift{frame}\n$frmseq\n";
    }
    $report .= "\nreference sequence $sbegin-$send\n$sbjseq";
    $report .= "------------------------------------------------------------\n";

    #print "bl2seq_frameshift out\n";
    return $report;

}

sub format_frameaa {
    my ( $sbjaa, $qryaa ) = @_;

    my @sbj = split /|/, lc $sbjaa;
    my $sbjlen = @sbj;

    my @qry;
    if ( defined $qryaa ) {
        @qry = split /|/, lc $qryaa;
    }
    my $qrylen = @qry;

    my %adjustments;
    for my $adj ( -$sbjlen .. +$sbjlen ) {
        my $matches   = 0;
        my $formatted = "";
        for my $q ( 0 .. $qrylen - 1 ) {
            my $s  = $q + $adj;
            my $ch = " " . $qry[$q] . " ";
            if ( $s >= 0 && $s < $sbjlen ) {
                if ( $qry[$q] eq $sbj[$s] ) {
                    $matches++;
                    $ch = "-" . uc $sbj[$s] . "-";
                }
            }
            $formatted .= $ch;
        }
        $adjustments{$adj}{formatted} = $formatted;
        $adjustments{$adj}{matches}   = $matches;
    }
    for my $adj (
        sort { $adjustments{$b}{matches} <=> $adjustments{$a}{matches} }
        keys %adjustments
      )
    {
        return $adjustments{$adj}{formatted};
        last;
    }
    return " " . join( "  ", @qry ) . " ";
}


########################################
# reports

# write standard TBL file
sub std_tbl_report {
    my ( $TBL, $genome, $inCDSs ) = @_;
    my $dbg = DEBUG;

    if ( !defined $TBL )    { return }
    if ( !defined $inCDSs ) { return }

    # adjust partials to edge of gap
    my $CDSs;
    for my $cds (@$inCDSs) {
        my $genesym = gene_name_to_symbol( $$cds{gene_name} );

 #print "name: $$cds{gene_name}  pseudo: $$cds{is_pseudogene}  sym: $genesym\n";
        if ( $$cds{is_pseudogene} && !defined $genesym ) {
            next;
        }

        # copy the gene
        # initialize codon start
        my %new;
        %new = %$cds;
        my @exons;
        for my $exon ( @{ $new{exons} } ) {
            my %tmp = %$exon;
            push @exons, \%tmp;
        }

        push @$CDSs, \%new;
    }

    # establish gene boundaries from CDSs
    my %genes;
    my $order = 0;
    for my $cds (@$CDSs) {
        my $strand = 1;
        if ( $$cds{start_site} > $$cds{stop_site} ) { $strand = -1 }
        my $gene_symbol = gene_name_to_symbol( $$cds{gene_name} );
        my $locus_tag   = gene_name_to_locus( $$cds{gene_name} );
#print "transcript: $$cds{gene_name}  symbol: $gene_symbol";
#if ( defined $locus_tag ) { print "  locus: $locus_tag" }
#print "\n";
        my $gene_tag = $gene_symbol;
        if ( defined $locus_tag ) {
            $gene_tag = $locus_tag;
        }
        $$cds{gene_tag} = $gene_tag;

        if ( exists $genes{$gene_tag} && $strand == $genes{$gene_tag}{strand} )
        {
            push @{ $genes{$gene_tag}{transcripts} }, $cds;

            #            my @exons = @{ $$cds{exons} };
            #            my $start = $exons[0]{dna_begin};
            #            my $stop  = $exons[ @exons - 1 ]{dna_end};
            my $start = $$cds{start_site};
            my $stop  = $$cds{stop_site};
            if ( $$cds{orientation} * $start <
                $$cds{orientation} * $genes{$gene_tag}{dna_begin} )
            {
                $genes{$gene_tag}{dna_begin}        = $start;
                $genes{$gene_tag}{start_truncation} = $$cds{start_truncation};
            }
            elsif ($genes{$gene_tag}{dna_begin} == $start
                && $$cds{start_truncation} )
            {
                $genes{$gene_tag}{start_truncation} = $$cds{start_truncation};
            }
            if ( $$cds{orientation} * $stop >
                $$cds{orientation} * $genes{$gene_tag}{dna_end} )
            {
                $genes{$gene_tag}{dna_end}         = $stop;
                $genes{$gene_tag}{stop_truncation} = $$cds{stop_truncation};
            }
            elsif ($genes{$gene_tag}{dna_end} == $stop
                && $$cds{stop_truncation} )
            {
                $genes{$gene_tag}{stop_truncation} = $$cds{stop_truncation};
            }
        }
        else {

#            if ( exists $genes{ $gene_tag } ) {
#                for my $i ( "II", "III", "IV", "V", "VI", "VII", "VII", "VIII", "IX", "X" ) {
#                    my $new_tag;
#                    if ( defined $gene_symbol ) {
#                        $gene_symbol .= "-$i";
#                        $new_tag = $gene_symbol;
#                    }
#                    if ( defined $locus_tag ) {
#                        $locus_tag .= "v$i";
#                        $new_tag = $locus_tag;
#                    }
#                    if ( ! exists $genes{ $new_tag } ) { last }
#                }
#            }
            $genes{$gene_tag}{order}  = $order++;
            $genes{$gene_tag}{strand} = $strand;
            $genes{$gene_tag}{genome} = $$cds{genome};
            my @tmp = ($cds);
            $genes{$gene_tag}{transcripts}      = \@tmp;
            $genes{$gene_tag}{orientation}      = $$cds{orientation};
            $genes{$gene_tag}{gene_synonym}     = $$cds{gene_synonym};
            $genes{$gene_tag}{gene_symbol}      = $gene_symbol;
            $genes{$gene_tag}{locus_tag}        = $locus_tag;
            $genes{$gene_tag}{dna_begin}        = $$cds{start_site};
            $genes{$gene_tag}{start_truncation} = $$cds{start_truncation};
            $genes{$gene_tag}{dna_end}          = $$cds{stop_site};
            $genes{$gene_tag}{stop_truncation}  = $$cds{stop_truncation};
            $genes{$gene_tag}{pseudogenes}      = 0;
        }

        if ( $$cds{is_pseudogene} ) {
            $genes{$gene_tag}{pseudogenes}++;
        }
    }

    # genomic sequence
    print $TBL ">Features $$genome{id}\n";

    # features
    for my $gene_tag ( sort { $genes{$a}{order} <=> $genes{$b}{order} }
        keys %genes )
    {

        # gene
        my $start = $genes{$gene_tag}{dna_begin};
        my $end = $genes{$gene_tag}{dna_end};
        my $start2;
        my $end2;
        if ( $genes{$gene_tag}{genome}{is_circular} ) {
            if ( $genes{$gene_tag}{strand} == 1 ) {
                if ( $end > $genes{$gene_tag}{genome}{original_seqlen} ) {
                    $end2 = $end - $genes{$gene_tag}{genome}{original_seqlen};
                    $end = $genes{$gene_tag}{genome}{original_seqlen};                    
                    $start2 = 1;
                    if ( $genes{$gene_tag}{start_truncation} ) { $start = "<$start" }
                    if ( $genes{$gene_tag}{stop_truncation} ) { $end2 = ">$end2" }
                }
            }
            else {
                if ( $start > $genes{$gene_tag}{genome}{original_seqlen} ) {
                    $start2 = $start - $genes{$gene_tag}{genome}{original_seqlen};
                    $start = $genes{$gene_tag}{genome}{original_seqlen};
                    $end2 = 1;
                    if ( $genes{$gene_tag}{start_truncation} ) { $start2 = "<$start2" }
                    if ( $genes{$gene_tag}{stop_truncation} ) { $end = ">$end" }
                }
            }
        }
        else {
            if ( $genes{$gene_tag}{start_truncation} ) { $start = "<$start" }
            if ( $genes{$gene_tag}{stop_truncation} ) { $end = ">$end" }
        }
        print $TBL "$start\t$end\tgene\n";
        if ( defined $start2 ) {
        print $TBL "$start2\t$end2\tgene\n";
        }
        if ( defined $genes{$gene_tag}{locus_tag} ) {
            print $TBL "\t\t\tlocus_tag\t$genes{$gene_tag}{locus_tag}\n";
        }
        if ( defined $genes{$gene_tag}{gene_symbol} ) {
            print $TBL "\t\t\tgene\t$genes{$gene_tag}{gene_symbol}\n";
        }
        if ( defined $genes{$gene_tag}{gene_synonym} ) {
            print $TBL "\t\t\tgene_syn\t$genes{$gene_tag}{gene_synonym}\n";
        }
        if (   $genes{$gene_tag}{pseudogenes} == 1
            && $genes{$gene_tag}{transcripts} == 1 )
        {

        }

        my @transcripts = @{ $genes{$gene_tag}{transcripts} };

        # if all CDSs are pseudogenes flag gene as pseudogene
        if ( $genes{$gene_tag}{pseudogenes} == @transcripts ) {
            print $TBL "\t\t\tpseudogene\tunitary\n";
        }

    #        # pseudogene within single-transcript gene, add product and notes to gene
    #        if ( $genes{ $gene_tag }{pseudogenes} == 1 && @transcripts == 1 ) {
    #            print $TBL "\t\t\tproduct\t$transcripts[0]{product_name}\n";
    #            my $pseudogene_note = is_pseudogene( $transcripts[0], 1 );
##print "psnote: $pseudogene_note\n";
       #            if ( defined $transcripts[0]{note} && $transcripts[0]{note} gt "" ) {
##print "note: $transcripts[0]{note}\n";
        #                if ( defined $pseudogene_note && $pseudogene_note gt "" ) {
        #                    $pseudogene_note = "$pseudogene_note; $transcripts[0]{note}";
        #                }
        #                else {
        #                    $pseudogene_note = $transcripts[0]{note};
        #                }
        #            }
##print "finalnote: $pseudogene_note\n";
        #            if ( defined $pseudogene_note && $pseudogene_note gt "" ) {
        #                print $TBL "\t\t\tnote\t$pseudogene_note\n";
        #                tbl_error_features( $TBL, $transcripts[0] );
        #            }
        #            next;
        #        }

        # transcripts
        for my $cds (@transcripts) {

            # pseudogene within multi-transcript gene
            if ( $$cds{is_pseudogene} ) {
                my $start = $$cds{start_site};
                if ( $$cds{start_truncation} ) { $start = "<$start" }
                my $end = $$cds{stop_site};
                if ( $$cds{stop_truncation} ) { $end = ">$end" }
                print $TBL "$start\t$end\tmisc_feature\n";
                if ( defined $genes{$gene_tag}{locus_tag} ) {
                    print $TBL
                      "\t\t\tlocus_tag\t$genes{$gene_tag}{locus_tag}\n";
                }
                if ( defined $genes{$gene_tag}{gene_symbol} ) {
                    print $TBL "\t\t\tgene\t$genes{$gene_tag}{gene_symbol}\n";
                }
                print $TBL "\t\t\tproduct\t$$cds{product_name}\n";

                my $pseudogene_note = is_pseudogene( $cds, 1 );
                if ( defined $$cds{note} && $$cds{note} gt "" ) {
                    if ( defined $pseudogene_note && $pseudogene_note gt "" ) {
                        $pseudogene_note = "$pseudogene_note; $$cds{note}";
                    }
                    else {
                        $pseudogene_note = $$cds{note};
                    }
                }
                if ( defined $pseudogene_note && $pseudogene_note gt "" ) {
                    print $TBL "\t\t\tnote\t$pseudogene_note\n";
                    tbl_error_features( $TBL, $cds,
                        $genes{$gene_tag}{locus_tag} );
                }
                next;
            }

            # CDS
            my $location = get_gene_location($cds);
            my @exons    = split /,/, $location;

            for my $i ( 0 .. @exons - 1 ) {
                my $exon = $exons[$i];
                my ( $exstart, $exend ) = split /\.\./, $exon;
                if ( $i == 0 ) {
                    print $TBL "$exstart\t$exend\tCDS\n";
                }
                else {
                    print $TBL "$exstart\t$exend\n";
                }
            }

            print $TBL "\t\t\tcodon_start\t$$cds{codon_start}\n";

            if ( defined $$cds{alternate_startcodon} ) {
                print $TBL "\t\t\ttransl_except\t"
                  . "(pos:$$cds{alternate_startcodon}{dna_begin}..$$cds{alternate_startcodon}{dna_end},"
                  . "aa:$$cds{alternate_startcodon}{aa})\n";
            }
            if ( defined $$cds{stopcodon_readthru} ) {
                print $TBL "\t\t\ttransl_except\t"
                  . "(pos:$$cds{stopcodon_readthru}{dna_begin}..$$cds{stopcodon_readthru}{dna_end},"
                  . "aa:$$cds{stopcodon_readthru}{aa})\n";
            }
            if ( defined $$cds{alternate_stopcodon} ) {
                print $TBL "\t\t\ttransl_except\t"
                  . "(pos:$$cds{alternate_stopcodon}{dna_begin}..$$cds{alternate_stopcodon}{dna_end},"
                  . "aa:$$cds{alternate_stopcodon}{aa})\n";
            }

            # identification
            print $TBL "\t\t\tprotein_id\t$$cds{gene_id}\n";
            if ( defined $genes{$gene_tag}{locus_tag} ) {
                print $TBL "\t\t\tlocus_tag\t$genes{$gene_tag}{locus_tag}\n";
            }
            if ( defined $genes{$gene_tag}{gene_symbol} ) {
                print $TBL "\t\t\tgene\t$genes{$gene_tag}{gene_symbol}\n";
            }
            print $TBL "\t\t\tproduct\t$$cds{product_name}\n";

            # ribosomal slippage
            if ( allow_ribosomal_slippage( $$cds{ref_id} ) ) {
                for my $exon ( @{ $$cds{exons} } ) {
                    if ( defined $$exon{slippage} ) {
                        print $TBL "\t\t\tribosomal_slippage\n";
                        last;
                    }
                }
            }

            # RNA editing
            my $rnaedit_note;
            my $rnaedit;
            my $rna_editing = 0;
            for my $exon ( @{ $$cds{exons} } ) {
                if ( defined $$exon{rna_edit} ) {
                    if ( !$rna_editing ) {
                        print $TBL "\t\t\texception\tRNA editing\n";
                    }
                    $rna_editing++;
                    $rnaedit      = $$exon{rna_edit};
                    $rnaedit_note = $$exon{rna_edit}{note};
                }
            }

            # notes
            my $note;
            if ( defined $$cds{note} ) {
                $note = $$cds{note};
            }

            #            if ( $location =~ /[><]/ ) {
            #                if ( defined $note ) {
            #                    if ( $note !~ /; *$/ ) {
            #                        $note =~ s/  *$//;
            #                        $note .= "; ";
            #                    }
            #                    $note .= "incomplete protein sequence";
            #                }
            #                else {
            #                    $note = "incomplete protein sequence";
            #                }
            #            }
            if ( defined $$cds{cdsnotes} && @{ $$cds{cdsnotes} } ) {
                my $message = join( "; ", @{ $$cds{cdsnotes} } );
                if ( defined $note ) {
                    if ( $note !~ /; *$/ ) {
                        $note =~ s/  *$//;
                        $note .= "; ";
                    }
                    $note .= $message;
                }
                else {
                    $note = $message;
                }
            }
            if ( defined $rnaedit_note ) {
                if ( defined $note ) {
                    if ( $note !~ /; *$/ ) {
                        $note =~ s/  *$//;
                        $note .= "; ";
                    }
                    $note .= $rnaedit_note;
                }
                else {
                    $note = $rnaedit_note;
                }
            }
            if ( defined $$cds{noncanonical_splicing}
                && @{ $$cds{noncanonical_splicing} } )
            {
                my $message = "non-canonical splicing";

                if ( defined $note ) {
                    if ( $note !~ /, *$/ ) {
                        $note =~ s/  *$//;
                        $note .= "; ";
                    }
                    $note .= $message;
                }
                else {
                    $note = $message;
                }
            }
            if ( defined $$cds{equivalent_splicing} ) {
                my @tmp     = keys %{ $$cds{equivalent_splicing} };
                my $message =
"ambiguous results from splicing algorithm, other possibilities exist";

                if ( defined $note ) {
                    if ( $note !~ /; *$/ ) {
                        $note =~ s/  *$//;
                        $note .= "; ";
                    }
                    $note .= $message;
                }
                else {
                    $note = $message;
                }
            }
            if ( defined $note ) {
                print $TBL "\t\t\tnote\t$note\n";
            }

            if ( defined $rnaedit ) {
                print $TBL
                  "$$rnaedit{dna_begin}\t$$rnaedit{dna_end}\tmisc_feature\n";
                if ( defined $genes{$gene_tag}{locus_tag} ) {
                    print $TBL
                      "\t\t\tlocus_tag\t$genes{$gene_tag}{locus_tag}\n";
                }
                my $oldseq =
                  subsequence( $$genome{sequence}, $$rnaedit{dna_begin},
                    $$rnaedit{dna_end} );
                print $TBL
"\t\t\tnote\tlocation of RNA editing ($oldseq => $$rnaedit{newseq}) in $$cds{product_name}\n";
            }
        }
    }

    # mature peptides
    for my $cds (@$CDSs) {
        if ( defined $$cds{mature_peps} && @{ $$cds{mature_peps} } ) {
            if ( is_pseudogene($cds) ) { next }

            my $gene_tag = $$cds{gene_tag};
            print $TBL "\n>Features $$cds{gene_id}\n";

            for my $pep ( @{ $$cds{mature_peps} } ) {
                my $start = $$pep{pep_begin};
                if ( $$pep{fuzzy_begin} ) {
                    $start = "<$start";
                }
                my $end = $$pep{pep_end};
                if ( $$pep{fuzzy_end} ) {
                    $end = ">$end";
                }
                if ( $$pep{product_name} =~ /\bsignal +peptide\b/i ) {
                    print $TBL "$start\t$end\tsig_peptide\n";
                }
                else {
                    print $TBL "$start\t$end\tmat_peptide\n";
                    print $TBL "\t\t\tproduct\t$$pep{product_name}\n";
                }
                if ( defined $genes{$gene_tag}{locus_tag} ) {
                    print $TBL
                      "\t\t\tlocus_tag\t$genes{$gene_tag}{locus_tag}\n";
                }
                if ( defined $genes{$gene_tag}{gene_symbol} ) {
                    print $TBL "\t\t\tgene\t$genes{$gene_tag}{gene_symbol}\n";
                }
                my $note;
                if ( defined $$pep{pep_note} ) { $note = $$pep{pep_note} }
                if ( defined $note ) {
                    print $TBL "\t\t\tnote\t$note\n";
                }
            }
        }
    }
}

sub tbl_error_features {
    my ( $TBL, $cds, $genetag ) = @_;

    if ( defined $$cds{frameshifts} && @{ $$cds{frameshifts} } ) {
        for my $fs ( @{ $$cds{frameshifts} } ) {
            my $fb = $$fs{frame_begin};
            if ( $fb > 0 ) { $fb = "+$fb" }
            my $fe = $$fs{frame_end};
            if ( $fe > 0 ) { $fe = "+$fe" }
            print $TBL "$$fs{dna_begin}\t$$fs{dna_end}\tmisc_feature\n";
            if ( defined $genetag ) {
                print $TBL "\t\t\tlocus_tag\t$genetag\n";
            }
            print $TBL
"\t\t\tnote\tapproximate location of $fb => $fe frameshift in $$cds{product_name}\n";
        }
    }
}

sub is_pseudogene {
    my ( $gene, $tbl_note ) = @_;
    if ( !defined $tbl_note ) { $tbl_note = 0 }

    my $pseudonote;
    if ( defined $$gene{cdserrors} && @{ $$gene{cdserrors} } ) {
        $pseudonote = join( "; ", @{ $$gene{cdserrors} } );
    }
    if ( defined $$gene{embeddedstops} && @{ $$gene{embeddedstops} } ) {
        my $message;
        if ( @{ $$gene{embeddedstops} } == 1 ) {
            $message = "CDS interrupted by stop codon";
        }
        elsif ( @{ $$gene{embeddedstops} } > 3 ) {
            $message = "CDS interrupted by many stop codons";
        }
        else {
            $message = "CDS interrupted by stop codons";
        }
        if ( !$tbl_note && @{ $$gene{embeddedstops} } <= 3 ) {
            $message .= " at " . join( ", ", @{ $$gene{embeddedstops} } );
            $message =~ s/, ([^,]+)$/ and $1/;
        }

        if ( defined $pseudonote ) {
            $pseudonote .= "; $message";
        }
        else {
            $pseudonote = $message;
        }
    }

    if ( defined $$gene{invalidsplicing} && @{ $$gene{invalidsplicing} } ) {
        my $message;
        if ( @{ $$gene{invalidsplicing} } == 1 ) {
            $message = "invalid splicing for intron";
        }
        else {
            $message = "invalid splicing for introns";
        }
        if ( !$tbl_note ) {
            $message .= " at " . join( ", ", @{ $$gene{invalidsplicing} } );
            $message =~ s/, ([^,]+)$/ and $1/;
        }
        if ( defined $pseudonote ) {
            $pseudonote .= "; $message";
        }
        else {
            $pseudonote = $message;
        }
    }

    if ( defined $$gene{frameshifts} && @{ $$gene{frameshifts} } ) {
        my $message;
        if ($tbl_note) {
            if ( @{ $$gene{frameshifts} } == 1 ) {
                $message = "frameshift in genome";
            }
            else {
                $message = "frameshifts in genome";
            }
        }
        else {
            $message = frameshift_message( @{ $$gene{frameshifts} } );
        }
        if ( defined $message ) {
            if ( defined $pseudonote ) {
                $pseudonote .= "; $message";
            }
            else {
                $pseudonote = $message;
            }
        }
    }

    if ( defined $$gene{gaperrors} && @{ $$gene{gaperrors} } ) {
        my $message;
        if ($tbl_note) {
            if ( @{ $$gene{gaperrors} } == 1 ) {
                $message = "gap sizing error in genome";
            }
            else {
                $message = "gap sizing errors in genome";
            }
        }
        else {
            $message = gaperror_message( @{ $$gene{gaperrors} } );
        }
        if ( defined $message ) {
            if ( defined $pseudonote ) {
                $pseudonote .= "; $message";
            }
            else {
                $pseudonote = $message;
            }
        }
    }

    return $pseudonote;
}

# write cds fasta
sub std_cds_report {
    my ( $CDS, $genome, $genes ) = @_;

    if ( !defined $CDS ) { return }

    my $ref_db = basename( get_parameter("reference_db") );

    for my $gene (@$genes) {

        my $location = get_gene_location($gene);
        my $cdna     = $$gene{cdna};

        my $defline = ">$$gene{gene_id}";
        if ( $$gene{is_pseudogene} ) { $defline .= " pseudogene" }
        $defline .= " location=$location";
        if ( defined $$gene{codon_start} ) {
            $defline .= " codon_start=$$gene{codon_start}";
        }
        else {
            $defline .= " codon_start=1";
        }
        if ( exists $$gene{translation_exception} ) {
            $defline .=
                " translation_exception="
              . "\"$$gene{translation_exception}{dna_begin}..$$gene{translation_exception}{dna_end}:"
              . "$$gene{translation_exception}{aa}\"";
        }
        $defline .=
" gene=\"$$gene{gene_name}\" product=\"$$gene{product_name}\" ref_db=\"$ref_db\" ref_id=\"$$gene{ref_id}\"\n";

        print $CDS $defline;

        $cdna =~ s/(.{1,60})/$1\n/g;
       print $CDS $cdna;

        if ( defined $$gene{mature_peps} ) {
            for my $pep ( @{ $$gene{mature_peps} } ) {

                my $location = get_matpep_location( $gene, $pep );

                my @tmp_description;
                if ( defined $$pep{pep_name} ) {
                    push @tmp_description, $$pep{pep_name};
                }
                if ( defined $$pep{product_name} ) {
                    push @tmp_description, $$pep{product_name};
                }
                my $pep_description = join( " ", @tmp_description );
                if ( $pep_description !~ /mat[ure]* *pep[tide]*/i ) {
                    $pep_description = "mature peptide, $pep_description";
                }

                my $defline = ">$$pep{pep_id}";
                if ( $$gene{is_pseudogene} ) { $defline .= " pseudogene" }
                $defline .=
" mat_peptide location=$location gene=\"$$gene{gene_name}\" product=\"$pep_description\"  ref_db=\"$$pep{ref_db}\" ref_id=\"$$pep{ref_id}\"\n";

                print $CDS $defline;

                my $cdna = "";
                for my $exon ( split /,/, $location ) {
                    $exon =~ s/[><]//g;
                    my ( $start, $end ) = split /\.\./, $exon;
                    $cdna .= subsequence( $$genome{sequence}, $start, $end );
                }
                $cdna =~ s/(.{1,60})/$1\n/g;
                print $CDS $cdna;
            }
        }
    }
}

# write peptide fasta
sub std_pep_report {
    my ( $PEP, $genome, $genes ) = @_;

    if ( !defined $PEP ) { return }

    my $ref_db = basename( get_parameter("reference_db") );

    for my $gene (@$genes) {

        my $location = get_gene_location($gene);
        my $protein  = $$gene{protein};
        my $numaa    = length($protein);
        if ( substr( $protein, $numaa - 1, 1 ) eq "*" ) {
            $numaa--;
        }

        my $defline = ">$$gene{gene_id}";
        if ( $$gene{is_pseudogene} ) { $defline .= " pseudogene" }
        $defline .= " location=$location";
        if ( defined $$gene{codon_start} ) {
            $defline .= " codon_start=$$gene{codon_start}";
        }
        else {
            $defline .= " codon_start=1";
        }
        if ( exists $$gene{translation_exception} ) {
            $defline .=
                " translation_exception="
              . "\"$$gene{translation_exception}{dna_begin}..$$gene{translation_exception}{dna_end}:"
              . "$$gene{translation_exception}{aa}\"";
        }
        $defline .=
" gene=\"$$gene{gene_name}\" product=\"$$gene{product_name}\" ref_db=\"$ref_db\" ref_id=\"$$gene{ref_id}\"\n";

        print $PEP $defline;

        $protein =~ s/(.{1,60})/$1\n/g;
        print $PEP $protein;

        if ( defined $$gene{mature_peps} ) {
            for my $pep ( @{ $$gene{mature_peps} } ) {

                my $protein = $$pep{pep_sequence};
                my $numaa   = length($protein);
                if ( substr( $protein, $numaa - 1, 1 ) eq "*" ) {
                    $numaa--;
                }

                my $location = get_matpep_location( $gene, $pep );

                my @tmp_description;
                if ( defined $$pep{pep_name} ) {
                    push @tmp_description, $$pep{pep_name};
                }
                if ( defined $$pep{product_name} ) {
                    push @tmp_description, $$pep{product_name};
                }
                my $pep_description = join( " ", @tmp_description );

                my $defline = ">$$pep{pep_id}";
                if ( $$gene{is_pseudogene} ) { $defline .= " pseudogene" }
                $defline .=
" mat_peptide location=$location gene=\"$$gene{gene_name}\" product=\"$pep_description\"  ref_db=\"$$pep{ref_db}\" ref_id=\"$$pep{ref_id}\"\n";

                print $PEP $defline;

                $protein =~ s/(.{1,60})/$1\n/g;
                print $PEP $protein;
            }
        }
    }
}

# write standard ALIGN report
sub std_align_report {
    my ( $ALIGN, $FS, $genome, $genes ) = @_;

    if ( !defined $ALIGN && !defined $FS ) { return }

    for my $gene (@$genes) {
        if ( defined $ALIGN ) { gene_align_report( $ALIGN, $gene ) }
        if ( defined $FS && is_pseudogene($gene) ) {
            gene_align_report( $FS, $gene );
        }
        if ( defined $$gene{mature_peps} ) {
            for my $pep ( @{ $$gene{mature_peps} } ) {
                if ( defined $ALIGN ) { gene_align_report( $ALIGN, $pep ) }
            }
        }
    }
    return;
}

sub std_autotasker_report {
    my ( $AT, $genome, $genes ) = @_;

    for my $gene (@$genes) {
        my $location = get_gene_location($gene);
        my $baseline =
            "$$genome{id}\t"
          . basename( get_parameter("reference_db") ) . "\t"
          . "$location\t"
          . "$$gene{ref_id}\t"
          . "$$gene{gene_name}\t"
          . "$$gene{product_name}";
        my $reported = 0;
        if ( $location =~ /^</ ) {
            my $truncation = maxval( 3, 3 * $$gene{num_reftrunc5} );
            print $AT "$baseline\tT5\t$truncation\n";
            $reported = 1;
        }
        if ( $location =~ />[0-9]+$/ ) {
            my $truncation = 3 * $$gene{num_reftrunc3} + 3;
            print $AT "$baseline\tT3\t$truncation\n";
            $reported = 1;
        }

        if ( defined $$gene{frameshifts} ) {
            for my $fs ( @{ $$gene{frameshifts} } ) {
                my $fb = $$fs{frame_begin};
                if ( $fb > 0 ) { $fb = "+$fb" }
                my $fe = $$fs{frame_end};
                if ( $fe > 0 ) { $fe = "+$fe" }
                print $AT
                  "$baseline\tFS\t$$fs{dna_begin}..$$fs{dna_end},$fb=>$fe\n";
                $reported = 1;
            }
        }

        if ( defined $$gene{gaperrors} ) {
            for my $ge ( @{ $$gene{gaperrors} } ) {
                print $AT
"$baseline\tGS\t$$ge{gap_dna_begin}..$$ge{gap_dna_end},$$ge{gapsize}=>$$ge{expectedsize}\n";
                $reported = 1;
            }
        }

        if ( defined $$gene{embeddedstops} ) {
            for my $stop ( @{ $$gene{embeddedstops} } ) {
                print $AT "$baseline\tES\t$stop\n";
                $reported = 1;
            }
        }

        if ( defined $$gene{cdserrors} ) {
            for my $error ( @{ $$gene{cdserrors} } ) {
                if ( $error =~ /start codon/i ) {
                    print $AT "$baseline\tC5\t\n";
                    $reported = 1;
                }
                elsif ( $error =~ /stop codon/i ) {
                    print $AT "$baseline\tC3\t\n";
                    $reported = 1;
                }
            }
        }

        if ( defined $$gene{invalidsplicing} ) {
            for my $intron ( @{ $$gene{invalidsplicing} } ) {
                print $AT "$baseline\tIS\t$intron\n";
                $reported = 1;
            }
        }

        if ( !$reported ) {
            print $AT "$baseline\tOK\t\n";
        }
    }
}

sub print_genehits {
    my (@genes) = @_;
    print format_genehits(@genes);
}

sub format_genehits {
    my (@genes) = @_;

    my $blanks50 = "                                                  ";
    my $header   =
        substr( "gene_id$blanks50", 0, 18 ) . " "
      . substr( "%id$blanks50",         0, 6 ) . " "
      . substr( "%sim$blanks50",        0, 6 ) . " "
      . substr( "%cov$blanks50",        0, 6 ) . " "
      . substr( "%t5$blanks50",         0, 5 ) . " "
      . substr( "%gap$blanks50",        0, 5 ) . " "
      . substr( "%t3$blanks50",         0, 5 ) . " "
      . substr( "start..stop$blanks50", 0, 16 ) . " "
      . substr( "pep sz$blanks50",      0, 6 ) . " "
      . substr( "ref sz$blanks50",      0, 6 ) . " "
      . substr( "ref id$blanks50",      0, 17 ) . " "
      . "definition";

    my $hits = "$header\n";

    for my $gene (@genes) {
        my $gene_id = substr( "$$gene{gene_id}$blanks50", 0, 18 );

        my $ident = "      ";
        if ( defined $$gene{pct_refidentity} ) {
            $ident = substr( "$$gene{pct_refidentity}      ", 0, 6 );
        }

        my $sim = "      ";
        if ( defined $$gene{pct_refsimilarity} ) {
            $sim = substr( "$$gene{pct_refsimilarity}      ", 0, 6 );
        }

        my $cov = "      ";
        if ( defined $$gene{pct_refcoverage} ) {
            $cov = substr( "$$gene{pct_refcoverage}      ", 0, 6 );
        }

        my $t5 = "     ";
        if ( defined $$gene{pct_reftrunc5} ) {
            $t5 = substr( "$$gene{pct_reftrunc5}      ", 0, 5 );
        }

        my $gap = "     ";
        if ( defined $$gene{pct_refgap} ) {
            $gap = substr( "$$gene{pct_refgap}      ", 0, 5 );
        }

        my $t3 = "     ";
        if ( defined $$gene{pct_reftrunc3} ) {
            $t3 = substr( "$$gene{pct_reftrunc3}      ", 0, 5 );
        }

        my $query_range = get_gene_location($gene);
        my @qranges = split /,/, $query_range;
        $query_range = shift @qranges;
        $query_range = substr( "$query_range$blanks50", 0, 16 );

        my $pepsz = "";
        if ( defined $$gene{protein_length} ) {
            $pepsz = $$gene{protein_length};
        }
        $pepsz = substr( "$pepsz$blanks50", 0, 6 );

        my $refsz = "";
        if ( defined $$gene{ref_length} ) { $refsz = $$gene{ref_length} }
        $refsz = substr( "$refsz$blanks50", 0, 6 );

        my $subject_id = $$gene{ref_id};
        $subject_id = substr( "$subject_id$blanks50", 0, 17 );

        my $product_name = $$gene{gene_name} . " | " . $$gene{product_name};
        my @prodwrapped = wraptext( $product_name, 50 );
        $product_name = shift @prodwrapped;
        my $line = join(
            " ",
            (
                $gene_id, $ident, $sim,        $cov,
                $t5,      $gap,   $t3,         $query_range,
                $pepsz,   $refsz, $subject_id, $product_name
            )
        );
        $hits .= "$line\n";

        for my $qrange (@qranges) {
            $product_name = "";
            if (@prodwrapped) { $product_name = shift @prodwrapped }
            $hits .= substr( $blanks50, 0, 18 ) . " "    # gene_id
              . substr( $blanks50,          0, 6 ) . " "     # %id
              . substr( $blanks50,          0, 6 ) . " "     # %sim
              . substr( $blanks50,          0, 6 ) . " "     # %cov
              . substr( $blanks50,          0, 5 ) . " "     # %t5
              . substr( $blanks50,          0, 5 ) . " "     # %gap
              . substr( $blanks50,          0, 5 ) . " "     # %t3
              . substr( "$qrange$blanks50", 0, 16 ) . " "    # dna start..stop
              . substr( "$blanks50",        0, 6 ) . " "     # pep sz
              . substr( "$blanks50",        0, 6 ) . " "     # ref sz
              . substr( $blanks50,          0, 17 ) . " "    # ref id
              . $product_name . "\n";
        }
        for $product_name (@prodwrapped) {
            $hits .= substr( $blanks50, 0, 18 ) . " "        # gene_id
              . substr( $blanks50,   0, 6 ) . " "            # %id
              . substr( $blanks50,   0, 6 ) . " "            # %sim
              . substr( $blanks50,   0, 6 ) . " "            # %cov
              . substr( $blanks50,   0, 5 ) . " "            # %t5
              . substr( $blanks50,   0, 5 ) . " "            # %gap
              . substr( $blanks50,   0, 5 ) . " "            # %t3
              . substr( $blanks50,   0, 16 ) . " "           # dna start..stop
              . substr( "$blanks50", 0, 6 ) . " "            # pep sz
              . substr( "$blanks50", 0, 6 ) . " "            # ref sz
              . substr( $blanks50,   0, 17 ) . " "           # ref id
              . $product_name . "\n";
        }

        # warn of frameshifts and gap errors in genomic sequence
        if ( defined $$gene{frameshifts} && @{ $$gene{frameshifts} } ) {
            my $message = frameshift_message( @{ $$gene{frameshifts} } );
            $hits .= substr( $blanks50, 0, 20 ) . $message . "\n";
        }
        if ( defined $$gene{gaperrors} && @{ $$gene{gaperrors} } ) {
            my $message = gaperror_message( @{ $$gene{gaperrors} } );
            $hits .= substr( $blanks50, 0, 20 ) . $message . "\n";
        }

        # warn of other embedded stop codons
        if ( defined $$gene{embeddedstops} && @{ $$gene{embeddedstops} } ) {
            my $message = join( ", ", @{ $$gene{embeddedstops} } );
            if ( @{ $$gene{embeddedstops} } > 5 ) {
                $message = "many stop codons embedded in CDS";
            }
            elsif ( @{ $$gene{embeddedstops} } > 1 ) {
                $message =~ s/, ([^,]*)$/ and $1/;
                $message = "stop codons embedded in CDS at $message";
            }
            else {
                $message = "stop codon embedded in CDS at $message";
            }
            $hits .= substr( $blanks50, 0, 20 ) . $message . "\n";
        }

        # warn about invalid or non-canonical splicing
        if ( defined $$gene{invalidsplicing}
            && @{ $$gene{invalidsplicing} } )
        {
            my $message =
              "invalid splicing for introns at "
              . join( ", ", @{ $$gene{invalidsplicing} } );
            if ( @{ $$gene{invalidsplicing} } > 1 ) {
                $message =~ s/, ([^,]*)$/ and $1/;
            }
            else {
                $message =~ s/ introns/ intron /;
            }

            $hits .= substr( $blanks50, 0, 20 ) . $message . "\n";
        }
        if ( defined $$gene{noncanonical_splicing}
            && @{ $$gene{noncanonical_splicing} } )
        {
            my $message =
              "non-canonical splicing for introns at "
              . join( ", ", @{ $$gene{noncanonical_splicing} } );
            if ( @{ $$gene{noncanonical_splicing} } > 1 ) {
                $message =~ s/, ([^,]*)$/ and $1/;
            }
            else {
                $message =~ s/ introns / intron /;
            }
            $hits .= substr( $blanks50, 0, 20 ) . $message . "\n";
        }

        # warn about translation exception
        if ( defined $$gene{alternate_startcodon} ) {
            my $message =
                "translation exception "
              . "$$gene{alternate_startcodon}{dna_begin}..$$gene{alternate_startcodon}{dna_end}:$$gene{alternate_startcodon}{aa}";
            $hits .= substr( $blanks50, 0, 20 ) . $message . "\n";
        }
        if ( defined $$gene{stopcodon_readthru} ) {
            my $message =
                "translation exception "
              . "$$gene{stopcodon_readthru}{dna_begin}..$$gene{stopcodon_readthru}{dna_end}:$$gene{stopcodon_readthru}{aa}";
            $hits .= substr( $blanks50, 0, 20 ) . $message . "\n";
        }
        if ( defined $$gene{alternate_stopcodon} ) {
            my $message =
                "translation exception "
              . "$$gene{alternate_stopcodon}{dna_begin}..$$gene{alternate_stopcodon}{dna_end}:$$gene{alternate_stopcodon}{aa}";
            $hits .= substr( $blanks50, 0, 20 ) . $message . "\n";
        }
        for my $exon ( @{ $$gene{exons} } ) {
            if ( defined $$exon{rna_edit} ) {
                my $message;
                if ( !defined $$exon{rna_edit}{in_gap} ) {
                    $message =
                        "transcription exception "
                      . "$$exon{rna_edit}{dna_begin}..$$exon{rna_edit}{dna_end}:$$exon{rna_edit}{newseq}";
                }
                else {
                    $message =
                      "transcription exception likely falls within gap";
                }
                $hits .= substr( $blanks50, 0, 20 ) . $message . "\n";
            }
        }

        # other warnings
        if ( defined $$gene{cdsnotes} && @{ $$gene{cdsnotes} } ) {
            for my $message ( @{ $$gene{cdsnotes} } ) {
                $hits .= substr( $blanks50, 0, 20 ) . $message . "\n";
            }
        }
        if ( defined $$gene{cdserrors} && @{ $$gene{cdserrors} } ) {
            my $message = join( ", ", @{ $$gene{cdserrors} } );
            if ( @{ $$gene{cdserrors} } > 1 ) {
                $message =~ s/, ([^,]*)$/ and $1/;
            }
            $hits .= substr( $blanks50, 0, 20 ) . $message . "\n";
        }

        if ( defined $$gene{equivalent_splicing} ) {
            my @tmp = keys %{ $$gene{equivalent_splicing} };
            my $message = join( "; ", @tmp );
            if ( @tmp > 1 ) {
                $message =~ s/; ([^;]*)$/ and $1/;
                $message =
"ambiguous results from splicing algorithm, other possibilities: $message";
            }
            else {
                $message =
"ambiguous results from splicing algorithm, other possibility: $message";
            }
            $hits .= substr( $blanks50, 0, 20 ) . $message . "\n";
        }

        # miscellaneous warnings
        if ( $$gene{warning} ) {
            $hits .= substr( $blanks50, 0, 20 ) . $$gene{warning} . "\n";
        }

        if ( defined $$gene{mature_peps} ) {
            for my $pep ( @{ $$gene{mature_peps} } ) {
                my $location = get_matpep_location( $gene, $pep );
                my @qranges = split /,/, $location;

                my @tmp_description;
                if ( defined $$pep{pep_name} ) {
                    push @tmp_description, $$pep{pep_name};
                }
                if ( defined $$pep{product_name} ) {
                    push @tmp_description, $$pep{product_name};
                }
                my $product_name = join( " ", @tmp_description );
                my @prodwrapped = wraptext( $product_name, 50 );
                $product_name = shift @prodwrapped;

                my $qrange = shift @qranges;

                if ( defined $$pep{ref_id} ) {
                    my $pep_id = substr( "$$pep{pep_id}$blanks50", 0, 18 );

                    $ident = " ";
                    if ( defined $$pep{pct_refidentity} ) {
                        $ident = $$pep{pct_refidentity};
                    }
                    $ident = substr( "$ident$blanks50", 0, 6 );

                    $sim = " ";
                    if ( defined $$pep{pct_refsimilarity} ) {
                        $sim = $$pep{pct_refsimilarity};
                    }
                    $sim = substr( "$sim$blanks50", 0, 6 );

                    $cov = " ";
                    if ( defined $$pep{pct_refcoverage} ) {
                        $cov = $$pep{pct_refcoverage};
                    }
                    $cov = substr( "$cov$blanks50", 0, 6 );

                    $t5 = " ";
                    if ( defined $$pep{pct_reftrunc5} ) {
                        $t5 = $$pep{pct_reftrunc5};
                    }
                    $t5 = substr( "$t5$blanks50", 0, 5 );

                    $gap = " ";
                    if ( defined $$pep{pct_refgap} ) {
                        $gap = $$pep{pct_refgap};
                    }
                    $gap = substr( "$gap$blanks50", 0, 5 );

                    $t3 = " ";
                    if ( defined $$pep{pct_reftrunc3} ) {
                        $t3 = $$pep{pct_reftrunc3};
                    }
                    $t3 = substr( "$t3$blanks50", 0, 5 );

                    $query_range = substr( "$qrange$blanks50", 0, 16 );

                    my $pepsz = "";
                    if ( defined $$pep{pep_length} ) {
                        $pepsz = $$pep{pep_length};
                    }
                    $pepsz = substr( "$pepsz$blanks50", 0, 6 );

                    my $refsz = "";
                    if ( defined $$pep{ref_length} ) {
                        $refsz = $$pep{ref_length};
                    }

                    $refsz = substr( "$refsz$blanks50", 0, 6 );

                    $subject_id = substr( "$$pep{ref_id}$blanks50", 0, 17 );

                    my $line = join(
                        " ",
                        (
                            $pep_id, $ident,       $sim,
                            $cov,    $t5,          $gap,
                            $t3,     $query_range, $pepsz,
                            $refsz,  $subject_id,  $product_name
                        )
                    );

                    $hits .= "$line\n";
                }
                else {
                    $hits .= substr( $blanks50, 0, 18 ) . " "    # gene_id
                      . substr( $blanks50,          0, 6 ) . " "   # %id
                      . substr( $blanks50,          0, 6 ) . " "   # %sim
                      . substr( $blanks50,          0, 6 ) . " "   # %cov
                      . substr( $blanks50,          0, 5 ) . " "   # %t5
                      . substr( $blanks50,          0, 5 ) . " "   # %gap
                      . substr( $blanks50,          0, 5 ) . " "   # %t3
                      . substr( "$qrange$blanks50", 0, 16 ) . " "  # start..stop
                      . substr( "$blanks50",        0, 6 ) . " "   # pep sz
                      . substr( "$blanks50",        0, 6 ) . " "   # ref sz
                      . substr( $blanks50,          0, 17 ) . " "  # ref id
                      . $product_name    # product_name
                      . "\n";
                }

                for $qrange (@qranges) {
                    $product_name = "";
                    if (@prodwrapped) { $product_name = shift @prodwrapped }
                    $hits .= substr( $blanks50, 0, 18 ) . " "    # gene_id
                      . substr( $blanks50,          0, 6 ) . " "   # %cov
                      . substr( $blanks50,          0, 6 ) . " "   # %id
                      . substr( $blanks50,          0, 6 ) . " "   # %sim
                      . substr( $blanks50,          0, 5 ) . " "   # %t5
                      . substr( $blanks50,          0, 5 ) . " "   # %gap
                      . substr( $blanks50,          0, 5 ) . " "   # %t3
                      . substr( "$qrange$blanks50", 0, 16 ) . " "  # start..stop
                      . substr( "$blanks50",        0, 6 ) . " "   # pep sz
                      . substr( "$blanks50",        0, 6 ) . " "   # ref sz
                      . substr( $blanks50,          0, 17 ) . " "  # ref id
                      . $product_name . "\n";
                }
                for $product_name (@prodwrapped) {
                    $hits .= substr( $blanks50, 0, 18 ) . " "      # gene_id
                      . substr( $blanks50,   0, 6 ) . " "     # %cov
                      . substr( $blanks50,   0, 6 ) . " "     # %id
                      . substr( $blanks50,   0, 6 ) . " "     # %sim
                      . substr( $blanks50,   0, 5 ) . " "     # %t5
                      . substr( $blanks50,   0, 5 ) . " "     # %gap
                      . substr( $blanks50,   0, 5 ) . " "     # %t3
                      . substr( $blanks50,   0, 16 ) . " "    # dna start..stop
                      . substr( "$blanks50", 0, 6 ) . " "     # pep sz
                      . substr( "$blanks50", 0, 6 ) . " "     # ref sz
                      . substr( $blanks50,   0, 17 ) . " "    # ref id
                      . $product_name . "\n";
                }
            }
        }
    }

    return $hits;
}

sub gene_align_report {
    my ( $OUT, $gene ) = @_;
    if ( !defined $OUT ) { return }


    if ( defined $$gene{is_pseudogene} ) {
        for my $i ( 1 .. 24 ) { print $OUT "=====" }
        print $OUT "\n\n";
        print $OUT format_genehits($gene) . "\n";
    }
    else {
        for my $i ( 1 .. 24 ) { print $OUT "-----" }
        print $OUT "\n\n";
        print $OUT "$$gene{ref_id}  $$gene{product_name}\n\n";
    }
    if ( defined $$gene{pep_alignment} ) {
        print $OUT "$$gene{pep_alignment}\n";
        print $OUT "\n";
    }

    if ( defined $$gene{ref_alignment} ) {
        print $OUT "$$gene{ref_alignment}\n";
        print $OUT "\n";
    }

    if ( defined $$gene{mp_alignment} ) {
        print $OUT "$$gene{mp_alignment}\n";
        print $OUT "\n";
    }
}

sub print_blasthits {
    my ( $margin, @hits ) = @_;

    my $show_refdb = 0;
    my @blasthits  = @hits;
    if ( @blasthits > 0 && $blasthits[ @blasthits - 1 ] eq "refdb" ) {
        $show_refdb = 1;
        pop @blasthits;
    }

    my $padmargin = lpad( "", $margin );

    my $blanks40 = "                                        ";
    my $header   =
        substr( "subject_id$blanks40", 0, 17 ) . " "
      . substr( "sbjlen,$blanks40",       0, 6 ) . " "
      . substr( "subject_range$blanks40", 0, 14 ) . " "
      . substr( "ori$blanks40",           0, 3 ) . " "
      . substr( "query_range$blanks40",   0, 14 ) . " "
      . substr( "hsp",                    0, 3 ) . " "
      . substr( "qf$blanks40",            0, 2 ) . " "
      . substr( "fs$blanks40",            0, 2 ) . " "
      . substr( "gap$blanks40",           0, 3 ) . " "
      . substr( "%id$blanks40",           0, 6 ) . " "
      . substr( "%sim$blanks40",          0, 6 ) . " "
      . substr( "%cov$blanks40",          0, 6 ) . " "
      . substr( "evalue$blanks40",        0, 8 ) . " ";
    my $desc_size = 56;
    if ($show_refdb) {
        $header .= substr( "ref DB$blanks40", 0, 12 ) . " ";
        $desc_size -= 13;
    }
    $header .= "subject definition";

    print "$padmargin$header\n";

    my $oldsid = "";
    for my $blasthit (@blasthits) {

        my @hits = ($blasthit);
        if ( exists $$blasthit{permutations} ) {
            push @hits, @{ $$blasthit{permutations} };
        }

        for my $hit (@hits) {
            my $subject_id = substr( "$$hit{subject_id}$blanks40", 0, 17 );

            my $sbjlen = substr( "$$hit{subject_length}$blanks40", 0, 6 );

            my $subject_range =
              substr( "$$hit{subject_begin}-$$hit{subject_end}$blanks40",
                0, 14 );

            my $ori = substr( "$$hit{orientation}$blanks40", 0, 3 );

            my $query_range =
              substr( "$$hit{query_left}-$$hit{query_right}$blanks40", 0, 14 );

            my $hsps;
            my $fs;
            my $qf;
            if ( defined $$hit{hsps} ) {
                my @tmp = @{ $$hit{hsps} };
                $hsps = @tmp;
                $qf   = $tmp[0]{query_frame};
                $fs   = 0;
                my $lf = $qf;
                for my $i ( 1 .. @tmp - 1 ) {
                    if ( $tmp[$i]{query_frame} ne $lf ) {
                        $fs++;
                        $lf = $tmp[$i]{query_frame};
                        $qf = 0;
                    }
                }
            }
            else {
                $hsps = 1;
                $qf   = $$hit{query_frame};
                $fs   = 0;
            }
            $hsps = substr( "$hsps$blanks40", 0, 3 );

            my $query_frame = substr( "$qf$blanks40", 0, 2 );

            my $frameshifts = substr( "$fs$blanks40", 0, 2 );
            
            my $gap = "N  ";
            if ( ! defined $$hit{hsps} ) {
                if ( defined $$hit{in_gap} ) { $gap = "Y  " }
            }
            else {
                for my $hsp ( @{ $$hit{hsps} } ) {
                    if ( defined $$hsp{in_gap} ) {
                        $gap = "Y  ";
                        last;
                    }
                }
            }
            
            my $pct_identity = $$hit{pct_identity};
            if ( index( $pct_identity, "." ) < 0 ) { $pct_identity .= ".0" }
            $pct_identity = substr( "$pct_identity$blanks40", 0, 6 );

            my $pct_similarity = int( 10.0 * $$hit{pct_similarity} ) / 10.0;
            if ( index( $pct_similarity, "." ) < 0 ) { $pct_similarity .= ".0" }
            $pct_similarity = substr( "$pct_similarity$blanks40", 0, 6 );

            my $pct_scoverage = $$hit{pct_scoverage};
            if ( index( $pct_scoverage, "." ) < 0 ) { $pct_scoverage .= ".0" }
            $pct_scoverage = substr( "$pct_scoverage$blanks40", 0, 6 );

            my $evalue =
              substr( format_evalue( $$hit{evalue} ) . $blanks40, 0, 8 );

            my $ref = get_reference_seq( $$hit{subject_id} );

            my $definition;
            if ( defined $ref ) {
                $definition =
                    get_reference_name( $$hit{subject_id} ) . " | "
                  . get_reference_product( $$hit{subject_id} );
            }
            elsif ( !defined $$hit{subject_id} ) {
                $definition = "predicted protein";
            }
            elsif ( $$hit{subject_definition} =~ /product="([^"]+)"/ ) {
                $definition = $1;
            }
            elsif ( $$hit{subject_definition} =~ /product=([^ ])+ / ) {
                $definition = $1;
            }
            else {
                $definition = $$hit{subject_definition};
                $definition =~ s/^[^ ]* *//;
            }

            if ($show_refdb) {
                my $refdb = get_defline_tag( $$hit{subject_definition}, "db" );
                $refdb = substr( "$refdb$blanks40", 0, 12 );
                my $line = join(
                    " ",
                    (
                        $subject_id,     $sbjlen,         $subject_range,
                        $ori,            $query_range,    $hsps,
                        $query_frame,    $frameshifts,    $gap,
                        $pct_identity,   $pct_similarity, $pct_scoverage,
                        $evalue,         $refdb,          $definition
                    )
                );
                print "$padmargin$line\n";
            }
            else {
                my $line = join(
                    " ",
                    (
                        $subject_id,     $sbjlen,         $subject_range,
                        $ori,            $query_range,    $hsps,
                        $query_frame,    $frameshifts,    $gap,
                        $pct_identity,   $pct_similarity, $pct_scoverage,
                        $evalue,         $definition
                    )
                );
                print "$padmargin$line\n";
            }
        }
    }
    return;
}

sub stats_report {
    my ( $OUT, $tab_delimited, $genome, $genelist ) = @_;

    my $mutations   = 0;
    my $genes       = 0;
    my $cdsbases    = 0;
    my $refbases    = 0;
    my $pepbases    = 0;
    my $coveredaa   = 0;
    my $identicalaa = 0;
    my $similaraa   = 0;

    for my $gene (@$genelist) {
        if ( is_pseudogene($gene) ) {
            $mutations++;
        }
        else {
            $genes++;
            $cdsbases  += length( $$gene{cdna} );
            $refbases  += $$gene{ref_length};
            $pepbases  += $$gene{protein_length};
            $coveredaa +=
              int( $$gene{pct_refcoverage} / 100.0 * $$gene{ref_length} + 0.5 );
            $identicalaa +=
              int( $$gene{pct_refidentity} / 100.0 * $$gene{ref_length} + 0.5 );
            $similaraa +=
              int(
                $$gene{pct_refsimilarity} / 100.0 * $$gene{ref_length} + 0.5 );
        }
    }
    my $coverage   = "-";
    my $similarity = "-";
    my $identity   = "-";
    if ($refbases) {
        $coverage = int( 1000.0 * $coveredaa / $refbases + 0.5 ) / 10.0;
        if ( index( $coverage, "." ) < 0 ) { $coverage .= ".0" }
        $identity = int( 1000.0 * $identicalaa / $refbases + 0.5 ) / 10.0;
        if ( index( $identity, "." ) < 0 ) { $identity .= ".0" }
        $similarity = int( 1000.0 * $similaraa / $refbases + 0.5 ) / 10.0;
        if ( index( $similarity, "." ) < 0 ) { $similarity .= ".0" }
    }
    my $refdb = basename( get_parameter("reference_db") );

    if ($tab_delimited) {
        print $OUT join(
            "\t",
            (
                $$genome{id}, length( $$genome{sequence} ),
                $genes,       $mutations,
                $cdsbases,    $pepbases,
                $identity,    $similarity,
                $coverage,    $refdb
            )
          )
          . "\n";
    }
    else {
        print $OUT join(
            "  ",
            (
                rpad( $$genome{id}, 20 ),
                lpad( length( $$genome{sequence} ), 7 ),
                lpad( $genes,                       length("Genes") ),
                lpad( $mutations,                   length("Pseudogenes") ),
                lpad( $cdsbases,                    length("CDS Bases") ),
                lpad( $pepbases,                    length("Peptide Bases") ),
                lpad( $identity,                    length("%Ref Identity") ),
                lpad( $similarity,                  length("%Ref Similarity") ),
                lpad( $coverage,                    length("%Ref Coverage") ),
                " $refdb"
            )
          )
          . "\n";
          print $OUT "gaps: ";
          my @gaps = @{ $$genome{gaps} };
          for my $g ( 0..@gaps-1 ) {
              if ( $g > 0 ) { print $OUT ", " }
              print $OUT "$gaps[$g]{begin}-$gaps[$g]{end}"
          }
          print $OUT "\n";
    }
}

sub stats_header {
    my ( $OUT, $tab_delimited ) = @_;

    if ($tab_delimited) {
        print $OUT join(
            "\t",
            (
                "Sequence",
                "Length",
                "Genes",
                "Pseudogenes",
                "CDS Bases",
                "Peptide Bases",
                "%Ref Identity",
                "%Ref Similarity",
                "%Ref Coverage",
                "Ref DB"
            )
          )
          . "\n";
    }
    else {
        print $OUT join(
            "  ",
            (
                rpad( "Sequence", 20 ),
                lpad( "Length",          7 ),
                lpad( "Genes",           length("Genes") ),
                lpad( "Pseudogenes",     length("Pseudogenes") ),
                lpad( "CDS Bases",       length("CDS Bases") ),
                lpad( "Peptide Bases",   length("Peptide Bases") ),
                lpad( "%Ref Identity",   length("%Ref Identity") ),
                lpad( "%Ref Similarity", length("%Ref Similarity") ),
                lpad( "%Ref Coverage",   length("%Ref Coverage") ),
                " Ref DB"
            )
          )
          . "\n";
    }
}

sub frameshift_message {
    my (@frameshifts) = @_;

    my $fs_message = "";
    for my $frameshift (@frameshifts) {
        if ( length($fs_message) ) { $fs_message .= ", " }
        my $fb = $$frameshift{frame_begin};
        if ( $fb > 0 ) { $fb = "+$fb" }
        my $fe = $$frameshift{frame_end};
        if ( $fe > 0 ) { $fe = "+$fe" }
        $fs_message .=
          "$$frameshift{dna_begin}..$$frameshift{dna_end} ($fb=>$fe)";
    }
    if ( @frameshifts > 1 ) {
        $fs_message = "frameshifts in genome at or near positions $fs_message";
        $fs_message =~ s/, ([^,]*)$/ and $1/;
    }
    elsif ( @frameshifts == 1 ) {
        $fs_message = "frameshift in genome at or near position $fs_message";
    }

    return $fs_message;
}

sub gaperror_message {
    my (@gaperrors) = @_;

    my $ge_message = "";
    for my $gaperror (@gaperrors) {
        if ( length($ge_message) ) { $ge_message .= ", " }
        $ge_message .=
"$$gaperror{gap_dna_begin}..$$gaperror{gap_dna_end} ($$gaperror{gapsize}=>$$gaperror{expectedsize})";
    }
    if ( @gaperrors > 1 ) {
        $ge_message = "gap sizing errors in genome at $ge_message";
        $ge_message =~ s/, ([^,]*)$/ and $1/;
    }
    elsif ( @gaperrors == 1 ) {
        $ge_message = "gap sizing error at $ge_message";
    }

    return $ge_message;
}

sub wraptext {
    my ( $text, $width ) = @_;
    if ( !defined $text ) { $text = "" }
    $text =~ s/\t/ /g;
    $text =~ s/  */ /g;
    $text =~ s/^  *//;
    $text =~ s/  *$//;

    my @words = split /( )/, $text;
    my $phrase = ( shift @words );
    my @wrapped;
    for my $word (@words) {
        if ( $word eq " " ) {
            if ( length($phrase) < $width - 1 ) {
                if ( $phrase ne "" ) { $phrase .= " " }
            }
            else {
                $phrase =~ s/ $//;
                push @wrapped, $phrase;
                $phrase = "";
            }
        }
        elsif ( $phrase eq "" ) {
            $phrase = $word;
        }
        elsif ( length( $phrase . $word ) <= $width ) {
            $phrase .= $word;
        }
        else {
            $phrase =~ s/ $//;
            push @wrapped, $phrase;
            $phrase = $word;
        }
    }

    if ( $phrase ne "" ) { push @wrapped, $phrase }
    return @wrapped;
}

sub get_gene_location {
    my ($gene) = @_;
    my $dbg = DEBUG;

    my @exons = @{ $$gene{exons} };
    my $genome_length = $$gene{genome}->{seqlen};
    my $location = "";
    for my $exon (@exons) {
        if ( $location ne "" ) { $location .= "," }
        if ( $$exon{fuzzy_begin} ) { $location .= "<" }
        my $xbegin = $$exon{dna_begin};
        my $xend = $$exon{dna_end};
        if ( $gene->{genome}->{is_circular} ) {
            if ( $xbegin > $gene->{genome}->{original_seqlen} ) {
                $xbegin -= $gene->{genome}->{original_seqlen};
                $xend -= $gene->{genome}->{original_seqlen};
            }
        }
        $location .= $xbegin;
        $location .= "..";
        if ( $$exon{fuzzy_end} ) { $location .= ">" }
        $location .= $xend;
    }

    return $location;
}

sub get_matpep_location {
    my ( $gene, $pep ) = @_;

    my ( $poly_begin, $poly_end ) = ( $$pep{pep_begin}, $$pep{pep_end} );
    my $location =
      pep_coords_to_dna( $poly_begin, $poly_end, $$gene{codon_start},
        $$gene{exons}, 1 );

    if ( $$pep{fuzzy_begin} ) {
        if ( $location !~ /^</ ) {
            $location = "<$location";
        }
    }
    else {
        $location =~ s/<//;
    }

    if ( $$pep{fuzzy_end} ) {
        my @tmp = split /\.\./, $location;
        my $dnaend = pop @tmp;
        if ( $dnaend !~ /^>/ ) {
            $dnaend = ">$dnaend";
        }
        push @tmp, $dnaend;
        $location = join( "..", @tmp );
    }
    else {
        $location =~ s/>//;
    }

    return $location;
}

sub pep_coords_to_dna {
    my ( $pep_begin, $pep_end, $codon_start, $exons, $formatted ) = @_;

    my $cdna_begin = 3 * ( $pep_begin - 1 ) + $codon_start;
    my $cdna_end   = 3 * $pep_end + $codon_start - 1;

    return cdna_coords_to_dna( $cdna_begin, $cdna_end, $exons, $formatted );
}

sub pep_coords_to_ref {
    my ( $pep_begin, $pep_end, $codon_start, $exons, $ref ) = @_;
    my $dbg = DEBUG;

    my $cdna_begin = 3 * ( $pep_begin - 1 ) + $codon_start;
    my $cdna_end   = 3 * $pep_end + $codon_start - 1;

    my $begin_exon;
    for my $exon (@$exons) {

        if (   $$exon{cdna_begin} <= $cdna_begin
            && $$exon{cdna_end} >= $cdna_begin )
        {
            $begin_exon = $exon;
            last;
        }
    }

    my $end_exon;
    for my $exon (@$exons) {

        if ( $$exon{cdna_begin} <= $cdna_end && $$exon{cdna_end} >= $cdna_end )
        {
            $end_exon = $exon;
            last;
        }
    }

    if ( !defined $begin_exon ) {
        if ( !defined $end_exon ) {
            if ($dbg) {
                print "endE: $cdna_begin-$cdna_end not defined\n";
            }
            return undef;
        }
        $begin_exon = $end_exon;
    }
    elsif ( !defined $end_exon ) {
        $end_exon = $begin_exon;
    }
    if ($dbg) {
        print_hash( "begin=$cdna_begin", $begin_exon );
        print_hash( "end=$cdna_end",     $end_exon );
    }

    my $subject_begin;
    if (
        abs( $cdna_begin - $$begin_exon{cdna_begin} ) <=
        abs( $$begin_exon{cdna_end} - $cdna_begin ) )
    {
        $subject_begin = $$begin_exon{subject_begin} +
          int( ( $cdna_begin - $$begin_exon{cdna_begin} + 2 ) / 3 );
        if ($dbg) {
            print "beginFromB: $cdna_begin -> $subject_begin\n";
        }
    }
    else {
        $subject_begin = $$begin_exon{subject_end} - 1 -
          int( ( $$begin_exon{cdna_end} - $cdna_begin ) / 3 );
        if ($dbg) {
            print "beginFromE: $cdna_begin -> $subject_begin\n";
        }
    }

    my $subject_end;
    if (
        abs( $$end_exon{cdna_end} - $cdna_end ) <=
        abs( $cdna_end - $$end_exon{cdna_begin} ) )
    {
        $subject_end = $$end_exon{subject_end} - 1 -
          int( ( $$end_exon{cdna_end} - $cdna_end - 2 ) / 3 );
        if ($dbg) {
            print "endFromE: $cdna_end -> $subject_end\n";
        }
    }
    else {
        $subject_end = $$end_exon{subject_begin} +
          int( ( $cdna_end - $$end_exon{cdna_begin} ) / 3 );
        if ($dbg) {
            print "endFromB: $cdna_end -> $subject_end\n";
        }
    }

    return ( $subject_begin, $subject_end );
}

sub cdna_coords_to_dna {
    my ( $cdna_begin, $cdna_end, $exons, $formatted ) = @_;
    my $dbg = DEBUG;

    #if ( `whoami` =~ /jhoover/ ) { $dbg = 1 }
    if ($dbg) {
        print "return dna coordinates for cdna $cdna_begin-$cdna_end\n";
    }
    my @tmpexons;
    for my $exon ( @{$exons} ) {
        push @tmpexons, $exon;
        if ($dbg) {
            print
"exon $$exon{dna_begin}-$$exon{dna_end}|$$exon{cdna_begin}-$$exon{cdna_end}\n";
        }
    }

    my $exon = shift @tmpexons;
    my $ori  = $$exon{orientation};
    while ( $$exon{cdna_end} < $cdna_begin && @tmpexons ) {
        if ($dbg) {
            print
"skip exon $$exon{dna_begin}-$$exon{dna_end}|$$exon{cdna_begin}-$$exon{cdna_end}\n";
        }
        $exon = shift @tmpexons;
    }

    my $cdnapos = $cdna_begin;
    my $fuzzy   = 0;
    if (   defined $$exon{fuzzy_begin}
        && $$exon{fuzzy_begin}
        && $cdnapos == $$exon{cdna_begin} )
    {
        $fuzzy = 1;
    }
    my $adjust = 0;
    if ( defined $$exon{rna_edit} ) {
        if (   $cdnapos >= $$exon{rna_edit}{cdna_begin}
            && $cdnapos <= $$exon{rna_edit}{cdna_end} )
        {
            my ( $old, $new ) = split /\n/, $$exon{rna_edit}{align};
            my $editpos  = $cdnapos - $$exon{rna_edit}{cdna_begin};
            my $newpos   = $$exon{rna_edit}{cdna_begin} - 1;
            my $dnapos   = $$exon{rna_edit}{dna_begin} - $$exon{orientation};
            my $alignpos = 0;
            my $alignlen = length($new);
            $fuzzy = 0;
            if ($dbg) {
                print_hash( "rnaedit", $$exon{rna_edit} );
                print_hash( "exon",    $exon );
                print "cdnapos=$cdnapos\n";
            }
            while ( $alignpos < $alignlen && $newpos < $cdnapos ) {
                if (   substr( $old, $alignpos, 1 ) ne "-"
                    && substr( $new, $alignpos, 1 ) ne "-" )
                {
                    $fuzzy = 0;
                }
                else {
                    $fuzzy = 1;
                }
                if ( substr( $old, $alignpos, 1 ) eq "-" ) {
                    $adjust--;
                }
                else {
                    $dnapos += $$exon{orientation};
                }
                if ( substr( $new, $alignpos, 1 ) ne "-" ) {
                    $newpos++;
                }
                $alignpos++;
                if ($dbg) {
                    print
"$alignpos: cdnapos: $newpos  dnapos: $dnapos  fuzzy $fuzzy  adjust: $adjust\n";
                }
            }
            if ($fuzzy) {
                $cdnapos++;
            }
        }
        elsif ( $cdnapos > $$exon{rna_edit}{cdna_end} ) {
            $adjust =
              abs( $$exon{rna_edit}{dna_end} - $$exon{rna_edit}{dna_begin} ) -
              ( $$exon{rna_edit}{cdna_end} - $$exon{rna_edit}{cdna_begin} );
            if ($dbg) {
                print
"adjust begin $adjust = abs( $$exon{rna_edit}{dna_end} - $$exon{rna_edit}{dna_begin} ) - ( $$exon{rna_edit}{cdna_end} - $$exon{rna_edit}{cdna_begin} )\n";
            }
        }
    }
    if ($dbg) {
        print
"dna_begin = $$exon{dna_begin} + $ori * ( $cdnapos + $adjust - $$exon{cdna_begin} )\n";
    }
    my $dna_begin =
      $$exon{dna_begin} + $ori * ( $cdnapos + $adjust - $$exon{cdna_begin} );
    my $formatted_range = "";
    if ($fuzzy) {
        $formatted_range = "<";
    }
    $formatted_range .= "$dna_begin..";

    while ( $$exon{cdna_end} < $cdna_end && @tmpexons ) {
        if ( defined $$exon{fuzzy_end} && $$exon{fuzzy_end} ) {
            $formatted_range .= ">";
        }
        $formatted_range .= "$$exon{dna_end},";
        $exon = shift @tmpexons;
        if ( defined $$exon{fuzzy_begin} && $$exon{fuzzy_begin} ) {
            $formatted_range .= "<";
        }
        $formatted_range .= "$$exon{dna_begin}..";
    }
    $cdnapos = $cdna_end;
    $fuzzy   = 0;
    if (   defined $$exon{fuzzy_end}
        && $$exon{fuzzy_end}
        && $cdnapos == $$exon{cdna_end} )
    {
        $fuzzy = 1;
    }
    $adjust = 0;
    if ( defined $$exon{rna_edit} ) {
        if (   $cdnapos >= $$exon{rna_edit}{cdna_begin}
            && $cdnapos <= $$exon{rna_edit}{cdna_end} )
        {
            my ( $old, $new ) = split /\n/, $$exon{rna_edit}{align};
            my $editpos  = $cdnapos - $$exon{rna_edit}{cdna_begin};
            my $newpos   = $$exon{rna_edit}{cdna_begin} - 1;
            my $dnapos   = $$exon{rna_edit}{dna_begin} - $$exon{orientation};
            my $alignpos = 0;
            my $alignlen = length($new);
            $fuzzy = 0;

            if ($dbg) {
                print_hash( "rnaedit", $$exon{rna_edit} );
                print_hash( "exon",    $exon );
                print "cdnapos=$cdnapos\n";
            }
            while ( $alignpos < $alignlen && $newpos < $cdnapos ) {
                if (   substr( $old, $alignpos, 1 ) ne "-"
                    && substr( $new, $alignpos, 1 ) ne "-" )
                {
                    $fuzzy = 0;
                }
                else {
                    $fuzzy = 1;
                }
                if ( substr( $old, $alignpos, 1 ) eq "-" ) {
                    $adjust--;
                }
                else {
                    $dnapos += $$exon{orientation};
                }
                if ( substr( $new, $alignpos, 1 ) ne "-" ) {
                    $newpos++;
                }
                $alignpos++;
                if ($dbg) {
                    print
"$alignpos: cdnapos: $newpos  dnapos: $dnapos  fuzzy $fuzzy  adjust: $adjust\n";
                }
            }
        }
        elsif ( $cdnapos > $$exon{rna_edit}{cdna_end} ) {
            $adjust =
              abs( $$exon{rna_edit}{dna_end} - $$exon{rna_edit}{dna_begin} ) -
              ( $$exon{rna_edit}{cdna_end} - $$exon{rna_edit}{cdna_begin} );

#print "adjust end $adjust = abs( $$exon{rna_edit}{dna_end} - $$exon{rna_edit}{dna_begin} ) - ( $$exon{rna_edit}{cdna_end} - $$exon{rna_edit}{cdna_begin} )\n";
        }
    }
    my $dna_end =
      $$exon{dna_begin} + $ori * ( $cdnapos + $adjust - $$exon{cdna_begin} );
    if ($fuzzy) {
        $formatted_range .= ">";
    }
    $formatted_range .= $dna_end;

    # return results
    if ($dbg) {
        print
"$dna_end = $$exon{dna_begin} + $ori * ( $cdna_end - $$exon{cdna_begin} )\n";
    }

    if ($formatted) {
        return $formatted_range;
    }
    else {
        return ( $dna_begin, $dna_end );
    }
}

sub dna_coords_to_cdna {
    my ( $dnabegin, $dnaend, $exons ) = @_;
    my $dbg = DEBUG;

    #if ( `whoami` =~ /jhoover/ ) { $dbg = 1 }

    my $ori = $$exons[0]{orientation};

    my $cdnabegin;
    for my $exon (@$exons) {
        if (   $ori * $dnabegin <= $ori * $$exon{dna_end}
            && $ori * $dnabegin >= $ori * $$exon{dna_begin} )
        {
            if ($dbg) {
                print_hash( "dnabegin=$dnabegin", $exon );
            }
            my $adjust = 0;
            if ( defined $$exon{rna_edit} ) {
                if ( $ori * $dnabegin > $ori * $$exon{rna_edit}{dna_end} ) {
                    $adjust =
                      abs( $$exon{rna_edit}{cdna_end} -
                          $$exon{rna_edit}{cdna_begin} ) -
                      abs( $$exon{rna_edit}{dna_end} -
                          $$exon{rna_edit}{dna_begin} );
                    if ($dbg) {
                        print "post edit begin adjustment: $adjust\n";
                    }
                }
                elsif ( $ori * $dnabegin >= $ori * $$exon{rna_edit}{dna_begin} )
                {
                    my ( $old, $new ) = split /\n/, $$exon{rna_edit}{align};
                    my $cdnapos  = $$exon{rna_edit}{cdna_begin};
                    my $dnapos   = $$exon{rna_edit}{dna_begin};
                    my $newpos   = $cdnapos - 1;
                    my $alignpos = 0;
                    while ( $ori * $dnapos <= $dnabegin ) {
                        if ( substr( $new, $alignpos, 1 ) ne "-" ) {
                            if ( substr( $old, $alignpos, 1 ) ne "-" ) {
                                $newpos = $cdnapos;
                            }
                            $cdnapos++;
                        }
                        if ( substr( $old, $alignpos, 1 ) ne "-" ) {
                            $dnapos += $ori;
                        }
                        $alignpos++;
                    }
                    $cdnabegin = $newpos;
                    if ($dbg) {
                        print "mid edit begin estimate: $newpos\n";
                    }
                    last;
                }
            }
            my $cdnaoffset = $ori * ( $dnabegin - $$exon{dna_begin} ) + $adjust;
            $cdnabegin = $$exon{cdna_begin} + $cdnaoffset;
            if ($dbg) {
                print "cdnabegin: $cdnabegin\n";
            }
            last;
        }
    }

    my $cdnaend;
    for my $exon (@$exons) {
        if (   $ori * $dnaend <= $ori * $$exon{dna_end}
            && $ori * $dnaend >= $ori * $$exon{dna_begin} )
        {
            if ($dbg) {
                print_hash( "dnaend=$dnaend", $exon );
            }
            my $adjust = 0;
            if ( defined $$exon{rna_edit} ) {
                if ( $ori * $dnaend > $ori * $$exon{rna_edit}{dna_end} ) {
                    $adjust =
                      abs( $$exon{rna_edit}{cdna_end} -
                          $$exon{rna_edit}{cdna_begin} ) -
                      abs( $$exon{rna_edit}{dna_end} -
                          $$exon{rna_edit}{dna_begin} );
                    if ($dbg) {
                        print "post edit end adjustment: $adjust\n";
                    }
                }
                elsif ( $ori * $dnaend >= $ori * $$exon{rna_edit}{dna_begin} ) {
                    my ( $old, $new ) = split /\n/, $$exon{rna_edit}{align};
                    my $cdnapos  = $$exon{rna_edit}{cdna_begin};
                    my $dnapos   = $$exon{rna_edit}{dna_begin};
                    my $newpos   = $cdnapos - 1;
                    my $alignpos = 0;
                    while ( $ori * $dnapos <= $dnaend ) {
                        if ( substr( $new, $alignpos, 1 ) ne "-" ) {
                            if ( substr( $old, $alignpos, 1 ) ne "-" ) {
                                $newpos = $cdnapos;
                            }
                            $cdnapos++;
                        }
                        if ( substr( $old, $alignpos, 1 ) ne "-" ) {
                            $dnapos += $ori;
                        }
                        $alignpos++;
                    }
                    $cdnaend = $newpos;
                    if ($dbg) {
                        print "mid edit end estimate: $newpos\n";
                    }
                    last;
                }
            }
            my $cdnaoffset = $ori * ( $dnaend - $$exon{dna_begin} ) + $adjust;
            $cdnaend = $$exon{cdna_begin} + $cdnaoffset;
            if ($dbg) {
                print "cdnaend: $cdnaend\n";
            }
            last;
        }
    }

    if ( !defined $cdnabegin || !defined $cdnaend ) {
        print "could not convert location $dnabegin-$dnaend\n";
        for my $exon (@$exons) {
            print_hash( "exon", $exon );
        }
    }

    return ( $cdnabegin, $cdnaend );
}

sub format_evalue {
    my ($evalue) = @_;

    if ( defined $evalue ) {
        my ( $a, $b ) = split /[eE]\-/, $evalue;
        if ( !defined $b || !length($b) ) { $b = 0 }
        if ( $a < 0.0 ) { $a = 0.0 }
        if ( $a > 0 ) {
            while ( $a < 1.0 ) {
                $a = 10. * $a;
                $b++;
            }
            while ( $a >= 10.0 ) {
                $a = $a / 10.0;
                $b--;
            }
            $a = int( 10.0 * $a + 0.5 ) / 10.;
            if ( $a == 10.0 ) {
                $a = 1.0;
                $b++;
            }
            if ( $b >= 1 ) {
                $b =~ s/^0*//;
                $b = lpad( $b, 3, "0" );
            }
            else {
                $b = "";
            }
        }
        if ($b) {
            if ( index( $a, "." ) < 0 ) { $a .= ".0" }
            $evalue = "$a-e$b";
        }
        else {
            $evalue = $a;
        }
    }

    return $evalue;
}

# current date-time
sub now {
    my ( $sec, $min, $hour, $mday, $mon, $year, $wday, $yday, $isdst ) =
      localtime(time);

    $sec  = &lpad( $sec,         2, "0" );
    $min  = &lpad( $min,         2, "0" );
    $hour = &lpad( $hour,        2, "0" );
    $mday = &lpad( $mday,        2, "0" );
    $mon  = &lpad( $mon + 1,     2, "0" );
    $year = &lpad( $year + 1900, 4, "0" );
    my $now = "$year-$mon-$mday $hour:$min:$sec";
    return $now;
}

#
# right pad string to fixed length
sub rpad {
    my ( $text, $pad_len, $pad_char ) = @_;

    if ( $pad_len <= 0 ) {
        return "";
    }
    elsif ( $pad_len <= length($text) ) {
        return substr( $text, 0, $pad_len );
    }

    if ( !defined $pad_char ) {
        $pad_char = " ";
    }
    elsif ( length($pad_char) > 1 ) {
        $pad_char = substr( $pad_char, 0, 1 );
    }

    if ( $pad_len > length($text) ) {
        $text .= $pad_char x ( $pad_len - length($text) );
    }

    return "$text";
}

#
# left pad string to fixed length
sub lpad {
    my ( $text, $pad_len, $pad_char ) = @_;

    if ( $pad_len <= 0 ) {
        return "";
    }
    elsif ( $pad_len < length($text) ) {
        return substr( $text, 0, $pad_len );
    }

    if ( !defined $pad_char ) {
        $pad_char = " ";
    }
    elsif ( length($pad_char) > 1 ) {
        $pad_char = substr( $pad_char, 0, 1 );
    }

    if ( $pad_len > length($text) ) {
        $text = $pad_char x ( $pad_len - length($text) ) . $text;
    }

    return "$text";
}

sub print_hash {
    my ( $hashname, $hash ) = @_;

    our %hashes;

    #if ( exists $hashes{$hashname} ) { return }
    $hashes{$hashname} = $hash;

    print "\n$hashname\n";
    for my $hashkey ( sort keys %$hash ) {
        if ( !defined $$hash{$hashkey} ) {
            print "$hashkey=\n";
            next;
        }
        my $hashval = "$$hash{$hashkey}";
        if ( $hashval =~ /HASH/ ) {
            print_hash( "$hashname.$hashkey", $$hash{$hashkey} );
        }
        elsif ( $hashval =~ /ARRAY/ ) {
            for my $element ( @{ $$hash{$hashkey} } ) {
                if ( "$element" =~ /HASH/ ) {
                    print_hash( "$hashname.$hashkey.element", $element );
                }
                else {
                    print "$hashname=>ARRAY: $$hash{$hashkey}\n";
                }
                last;
            }
        }
        else {
            print "$hashkey=$hashval\n";
        }
    }
    return;
}

# ---- standardize reference data ----
# return standardized reference name (name is NOT required)
sub get_reference_name {
    my ( $reference, $refs, $N ) = @_;
    return get_reference_gene( $reference, $refs, $N );
}

sub gene_name_to_symbol {
    my ($gene_name) = @_;

    my $symbol = $gene_name;
    if ( get_parameter("use_locus_tags") ) {
        if ( $symbol =~ /^HYP\b/ ) {
            return undef;
        }
        else {
            $symbol =~ s/-[0-9].*//;
            return $symbol;
        }
    }
    else {
        return $symbol;
    }
}

sub gene_name_to_locus {
    my ( $gene_name, $is_pseudogene ) = @_;
    if ( !defined $is_pseudogene ) { $is_pseudogene = 0 }

    if ( get_parameter("use_locus_tags") ) {
        my $symbol = $gene_name;
        $symbol =~ s/-($roman_regexp)/d$1/;
        if ( $is_pseudogene || $symbol =~ /^HYP\b/ ) {
            $symbol =~ s/-/p/g;
        }
        else {
            $symbol =~ s/-[0-9].*$//;
        }
        $symbol =~ s/[^A-Za-z0-9-]//g;
        return "vigor_$symbol";
    }
    else {
        return undef;
    }
}

sub get_reference_gene {
    my ( $reference, $refs, $N ) = @_;
    if ( !defined $N ) { $N = "" }

    my $definition;
    my $seq = get_reference_seq( $reference, $refs );
    if ( defined $seq ) {
        $definition = $$seq{defline};
    }
    else {
        $definition = $reference;
    }

    #if ( $N =~ /\w/ ) { print "modiifier=$N  definition=$definition\n" }
    $definition =~ s/\s\s+/ /g;
    $definition =~ s/^[^ ]* */ /;
    $definition =~ s/^\s+//;
    $definition =~ s/\s+$//;
    my ($refid) = split / /, $definition;
    $definition = " $definition ";
    if ( $definition =~ / gene$N="([^"]+)"/i ) {
        my $name = $1;
        $name =~ s/^ +//;
        $name =~ s/ +$//;

#if ( $N =~ /\w/ ) { print "modifier=$N  definition=$definition  name=$name\n" }
        return $name;
    }
    elsif ( $definition =~ / gene$N=([^ ]+) /i ) {
        my $name = $1;
        $name =~ s/^ +//;
        $name =~ s/ +$//;

#if ( $N =~ /\w/ ) { print "modiifier=$N  definition=$definition  name=$name\n" }
        return $name;
    }
    else {
        return undef;
    }
}

sub get_reference_exclusions {
    my ($reference) = @_;

    my %exclusions;
    my $seq = get_reference_seq($reference);
    if ( !defined $seq ) { return %exclusions }

    my @list;
    if ( " $$seq{defline} " =~ / excludes_gene=N / ) {
        return %exclusions;
    }
    elsif ( " $$seq{defline} " =~ / excludes_gene="([^"]+)"/ ) {
        @list = split /,/, $1;
    }
    elsif ( " $$seq{defline} " =~ / excludes_gene=([^ ]+) / ) {
        @list = split /,/, $1;
    }

    for my $excluded (@list) {
        $exclusions{$excluded} = 1;
    }

    return %exclusions;
}

# ---- standardize reference data ----
# return genbank note
sub get_reference_note {
    my ( $reference, $refs, $N ) = @_;
    if ( !defined $N ) { $N = "" }

    my $definition;
    my $seq = get_reference_seq( $reference, $refs );
    if ( defined $seq ) {
        $definition = $$seq{defline};
    }
    else {
        $definition = $reference;
    }

    #print "definition=$definition\t";
    $definition =~ s/\s\s+/ /g;
    $definition =~ s/^[^ ]* */ /;
    $definition =~ s/^\s+//;
    $definition =~ s/\s+$//;
    my ($refid) = split / /, $definition;
    $definition = " $definition ";

    if ( $definition =~ / note$N="([^"]+)"/i ) {
        my $note = $1;
        $note =~ s/^ +//;
        $note =~ s/ +$//;
        return $note;
    }
    elsif ( $definition =~ / note$N=([^ ]+) /i ) {
        my $note = $1;
        $note =~ s/^ +//;
        $note =~ s/ +$//;
        return $note;
    }
    else {
        return undef;
    }
}

# ---- standardize reference data ----
# return genbank "note""s for

sub get_reference_gene_synonym {
    my ( $reference, $refs, $N ) = @_;
    if ( !defined $N ) { $N = "" }

    my $definition;
    my $seq = get_reference_seq( $reference, $refs );
    if ( defined $seq ) {
        $definition = $$seq{defline};
    }
    else {
        $definition = $reference;
    }

    $definition =~ s/\s\s+/ /g;
    $definition =~ s/^[^ ]* */ /;
    $definition =~ s/^\s+//;
    $definition =~ s/\s+$//;
    my ($refid) = split / /, $definition;
    $definition = " $definition ";

    if ( $definition =~ / gene_synonym$N="([^"]+)"/i ) {
        my $gene_synonym = $1;
        $gene_synonym =~ s/^ +//;
        $gene_synonym =~ s/ +$//;
        return $gene_synonym;
    }
    elsif ( $definition =~ / gene_synonym$N=([^ ]+) /i ) {
        my $gene_synonym = $1;
        $gene_synonym =~ s/^ +//;
        $gene_synonym =~ s/ +$//;
        return $gene_synonym;
    }
    else {
        return undef;
    }
}

sub is_reference_partial3 {
    my ($reference) = @_;

    my $definition;
    my $seq = get_reference_seq($reference);
    if ( defined $seq ) {
        $definition = $$seq{defline};
    }
    else {
        $definition = $reference;
    }

    $definition =~ s/\s\s+/ /g;
    $definition =~ s/^[^ ]* */ /;
    $definition =~ s/^\s+//;
    $definition =~ s/\s+$//;
    my ($refid) = split / /, $definition;
    $definition = " $definition ";
    if ( $definition =~ / tiny_exon3="i[]100000:0-9]+XXX"/ ) {
        return 1;
    }
    else {
        return 0;
    }
}

sub get_reference_product {
    my ( $reference, $refs, $N ) = @_;
    if ( !defined $N ) { $N = "" }

    my $definition;
    my $seq = get_reference_seq( $reference, $refs );
    if ( defined $seq ) {
        $definition = $$seq{defline};
    }
    else {
        $definition = $reference;
    }

    $definition =~ s/\s\s+/ /g;
    $definition =~ s/^[^ ]* */ /;
    $definition =~ s/^\s+//;
    $definition =~ s/\s+$//;
    my ($refid) = split / /, $definition;
    $definition = " $definition ";

    if ( " $definition " =~ / product$N="([^"]+)"/i ) {
        $definition = $1;
        $definition =~ s/^ +//;
        $definition =~ s/ +$//;
    }
    elsif ( " $definition " =~ / product$N=([^ ]+) /i ) {
        $definition = $1;
        $definition =~ s/^ +//;
        $definition =~ s/ +$//;
    }
    else {
        $definition = undef;
    }

    return $definition;
}

sub check_coverage {
    my ( $min_coverage, $feature, $context ) = @_;
    my $dbg = DEBUG;
    #if ( `whoami` =~ /jhoover/ && defined $$feature{subject_id} && $$feature{subject_id} eq "NP_067150.1" ) { $dbg = 1 }

    if ( exists $$feature{subject_id} ) {
        # check coverage of mature peptide
        if ( exists $$context{protein} ) {
            if ( !exists $$context{exons} ) { return 1 }
            if ($dbg) {
                print "CHECK COVERAGE OF MATPEP\n";
                print_blasthits( 2, $feature );
            }
            my $cdna_begin = 3 * $$feature{query_left} - 2;
            my $cdna_end   = 3 * $$feature{query_right};

            # extract exon alignment covering mature peptide
            my @tmpexons;
            my $prev;
            my $next;
            for my $exon ( @{ $$context{exons} } ) {
                if ( $$exon{cdna_end} < $cdna_begin ) {
                    $prev = $exon;
                    next;
                }
                elsif ( $$exon{cdna_begin} > $cdna_end ) {
                    my $next = $exon;
                    last;
                }
                my %tmp = %$exon;
                adjust_exon_cdnabegin( \%tmp,
                    maxval( $$exon{cdna_begin}, $cdna_begin ) );
                adjust_exon_cdnaend(
                    \%tmp,
                    minval( $$exon{cdna_end}, $cdna_end ),
                    $$feature{subject_length}
                );
                push @tmpexons, \%tmp;
            }

            # find how much missing coverage is due to genome gaps
            my $hidden5         = 0;
            my $hidden_internal = 0;
            my $hidden3         = 0;
            if ( $tmpexons[0]{fuzzy_begin} ) {
                $hidden5 = $$feature{subject_begin} - 1;
                if ($dbg)            { print "  hidden5=$hidden5\n" }
                if ( defined $prev ) {
                    $hidden5 =
                      minval( $hidden5,
                        $tmpexons[0]{subject_begin} - $$prev{subject_end} - 1 );
                    if ( $hidden5 < 0 ) { $hidden5 = 0 }
                    if ($dbg) { print "  hidden5=$hidden5\n" }
                }
            }
            for my $i ( 1 .. @tmpexons - 1 ) {
                my $hidden = maxval( 0,
                    $tmpexons[$i]{subject_begin} -
                      $tmpexons[ $i - 1 ]{subject_end} - 1 );
                $hidden_internal += $hidden;
                if ($dbg) { print "  hidden_internal=$hidden_internal\n" }
            }
            if ( $tmpexons[ @tmpexons - 1 ]{fuzzy_end} ) {
                $hidden3 = $$feature{subject_length} - $$feature{subject_end};
                if ($dbg)            { print "  hidden3=$hidden3\n" }
                if ( defined $next ) {
                    $hidden3 = minval( $hidden3,
                        $$next{subject_begin} -
                          $tmpexons[ @tmpexons - 1 ]{subject_end} - 1 );
                    if ( $hidden3 < 0 ) { $hidden3 = 0 }
                    if ($dbg) { print "  hidden3=$hidden3\n" }
                }
            }

            # adjust coverage for genome gaps;
            my $pct_hidden =
              100.0 * ( $hidden5 + $hidden_internal + $hidden3 ) /
              $$feature{subject_length};
            if ($dbg) { print "  pct_hidden=$pct_hidden\n" }

            # check coverage
            if ($dbg) {
                print
"if ( $$feature{pct_scoverage} + $pct_hidden >= $min_coverage ) { return 1 }\n";
            }
            if ( $$feature{pct_scoverage} + $pct_hidden >= $min_coverage ) {
                return 1;
            }
            return 0;
        }

        # check coverage of blast hit
        else {
            if ($dbg) {
                print "CHECK COVERAGE OF BLAST HIT\n";
                print_blasthits( 2, $feature );
            }
            my $genome_begin  = 1;
            my $genome_end    = $$context{seqlen};
            my $is_complete   = $$context{is_complete};
            my $feature_begin = $$feature{query_left};
            my $feature_end   = $$feature{query_right};
            my $proj_begin = $feature_begin;
            my $proj_end = $feature_end;
            for my $hsp ( @{ $$feature{hsps} } ) {
                my ( $left, $right ) = project_gene( $hsp );
                if ( $left < $proj_begin ) { $proj_begin = $left }
                if ( $right > $proj_end ) { $proj_end = $right }
            }

            # determine how much is missing
            my $missing5 = $$feature{subject_begin} - 1;
            my $missing3 = $$feature{subject_length} - $$feature{subject_end};
            my $missing_internal = 0;
            push my @hsps, @{ $$feature{hsps} };
            my $last = shift @hsps;
            for my $hsp (@hsps) {
                if ( $$hsp{subject_begin} > $$last{subject_end} ) {
                    $missing_internal += $$hsp{subject_begin} - $$last{subject_end};
                }
                $last = $hsp;
            }
            if ($dbg) {
                print "  missing5=$missing5  missing3=$missing3  missing_internal=$missing_internal\n";
            }

            # try to account for missing portion by gaps / genome edges
            if ( $$feature{orientation} == -1 ) {
                ( $genome_begin, $genome_end ) = ( $genome_end, $genome_begin );
                ( $feature_begin, $feature_end ) = ( $feature_end, $feature_begin );
                ( $proj_begin, $proj_end ) = ( $proj_end, $proj_begin );
            }
            
            # how much is due to internal gaps
            my $hidden_internal = 0;
            if ( $missing_internal > 0 ) {
                if ( $dbg ) { print "gaps between $feature_begin and $feature_end\n" }
                my @gaps = find_gaps( $context, $feature_begin, $feature_end, 1 );
                my $gaplength = 0;
                for my $gap (@gaps) {
                    $gaplength += abs( $gap->{end} - $gap->{begin} ) + 1;
                    if ( $dbg ) { print "  gap $$gap{begin} - $$gap{end} => $gaplength\n" }
                }
                $gaplength = int( $gaplength / 3.0 );
                $hidden_internal = $missing_internal <= $gaplength ? $missing_internal : $gaplength;
            }
            
            # how much is due to 5' gaps/edge
            my $hidden5 = 0;
            if ( $missing5 > 0 ) {
                my @gaps = find_gaps( $context, $proj_begin, $feature_begin, 1 );
                if ( ! $is_complete ) {
                    if (  $$feature{orientation} * $proj_begin < $$feature{orientation} * $genome_begin ) {
                        push @gaps, { begin => $proj_begin, end => $genome_begin - $$feature{orientation} };
                    }
                }
                if ( $dbg ) { print "gaps between $proj_begin and $feature_begin\n" }
                my $gaplength = 0;
                for my $gap (@gaps) {
                    $gaplength += abs( $gap->{end} - $gap->{begin} ) + 1;
                    if ( $dbg ) { print "  gap $$gap{begin} - $$gap{end} => $gaplength\n" }
                }
                $gaplength = int( $gaplength / 3.0 );
                $hidden5 = $missing5 <= $gaplength ? $missing5 : $gaplength;
            }
            
            # how much is due to 3' gaps/edge
            my $hidden3 = 0;
            if ( $missing3 > 0 ) {
                my @gaps = find_gaps( $context, $feature_end, $proj_end, 1 );
                if ( ! $is_complete ) {
                    if (  $$feature{orientation} * $proj_end > $$feature{orientation} * $genome_end ) {
                        push @gaps, { begin => $genome_end + $$feature{orientation}, end => $proj_end };
                    }
                }
                if ( $dbg ) { print "gaps between $feature_end and $proj_end\n" }
                my $gaplength = 0;
                for my $gap (@gaps) {
                    $gaplength += abs( $gap->{end} - $gap->{begin} ) + 1;
                    if ( $dbg ) { print "  gap $$gap{begin} - $$gap{end} => $gaplength\n" }
                }
                $gaplength = int( $gaplength / 3.0 );
                $hidden3 = $missing3 <= $gaplength ? $missing3 : $gaplength;
            }

            if ($dbg) {
                print
"  hidden5=$hidden5  hidden3=$hidden3  hidden_internal=$hidden_internal\n";
            }

            # adjust coverage for pieces hidden by gaps/edge
            my $pct_hidden =
              100.0 * ( $hidden5 + $hidden3 + $hidden_internal ) /
              $$feature{subject_length};
            if ($dbg) { print "  pct_hidden=$pct_hidden\n" }

            # check coverage
            if ($dbg) {
                print
"if ( $$feature{pct_scoverage} + $pct_hidden >= $min_coverage ) { return 1 }\n";
            }
            if ( $$feature{pct_scoverage} + $pct_hidden >= $min_coverage ) {
                return 1;
            }
            return 0;
        }

    }

    # check coverage of gene
    elsif ( exists $$feature{ref_id} ) {
        if ($dbg) {
            print "CHECK COVERAGE OF GENE\n";
            print_genehits($feature);
        }

        # find how much missing coverage is due to genome gaps
        my $hidden5 = 0;
        if ( $$feature{start_truncation} ) {
            $hidden5 = $$feature{num_reftrunc5};
        }
        my $hidden_internal = $$feature{num_refgap};
        my $hidden3         = 0;
        if ( $$feature{stop_truncation} ) {
            $hidden3 = $$feature{num_reftrunc3};
        }

        #        my @tmpexons         = @{ $$feature{exons} };
        #        if ( $$feature{start_truncation} ) {
        #            $hidden5 = $$feature{ref_begin} - 1;
        #        }
        #
        #        if ($dbg) {
        #            print "  hidden5=$hidden5\n";
        #            print "  hidden_internal=$hidden_internal\n";
        #        }
        #        for my $i ( 1 .. @tmpexons - 1 ) {
        #            my $hidden = maxval( 0,
        #                $tmpexons[$i]{subject_begin} -
        #                  $tmpexons[ $i - 1 ]{subject_end} -
        #                  1 );
        #            $hidden_internal += $hidden;
        #            if ($dbg) { print "  hidden_internal=$hidden_internal\n" }
        #        }
        #        if ( $$feature{stop_truncation} ) {
        #            $hidden3 = $$feature{ref_length} - $$feature{ref_end};
        #        }
        #        if ($dbg) { print "  hidden3=$hidden3\n" }

        # adjust coverage for genome gaps;
        my $pct_hidden =
          100.0 * ( $hidden5 + $hidden_internal + $hidden3 ) /
          $$feature{ref_length};
        if ($dbg) { print "  pct_hidden=$pct_hidden\n" }

        # check coverage
        if ($dbg) {
            print
"if ( $$feature{pct_refcoverage} + $pct_hidden >= $min_coverage ) { return 1 }\n";
        }
        if ( $$feature{pct_refcoverage} + $pct_hidden >= $min_coverage ) {
            return 1;
        }
        return 0;
    }

    # unexpected case
    else {
        return 0;
    }
}

sub allow_stopcodon_readthru {
    my ($ref_id) = @_;

    my $ref     = get_reference_seq($ref_id);
    my $defline = $$ref{defline};
    if ( " $defline " =~ / stopcodon_readthru=N /i ) {
        return 0;
    }
    elsif ( " $defline " =~ / stopcodon_readthru="Y/i ) {
        return 1;
    }
    elsif ( " $defline " =~ / stopcodon_readthru=Y/i ) {
        return 1;
    }
    return 0;
}

sub get_readthru_exception {
    my ($ref_id) = @_;

    my $ref     = get_reference_seq($ref_id);
    my $defline = $$ref{defline};
    if ( " $defline " =~ / stopcodon_readthru=N /i ) {
        return "*";
    }
    elsif ( " $defline " =~ / stopcodon_readthru="Y:([A-Z])/i ) {
        return $1;
    }
    elsif ( " $defline " =~ / stopcodon_readthru=Y:([A-Z])/ ) {
        return $1;
    }
    elsif ( " $defline " =~ / stopcodon_readthru="Y/i ) {
        return "X";
    }
    elsif ( " $defline " =~ / stopcodon_readthru=Y/i ) {
        return "X";
    }
    else {
        return "*";
    }
}

sub allow_ribosomal_slippage {
    my ($ref_id) = @_;

    my $ref     = get_reference_seq($ref_id);
    my $defline = $$ref{defline};
    if ( " $defline " =~ / ribosomal_slippage=N /i ) {
        return 0;
    }
    elsif ( " $defline " =~ / ribosomal_slippage="{0,1}Y/ ) {
        return 1;
    }
    return 0;
}

sub allow_rna_editing {
    my ($ref_id) = @_;

    my $ref     = get_reference_seq($ref_id);
    my $defline = $$ref{defline};

    if ( " $defline " =~ / rna_editing=N /i ) {
        return 0;
    }
    elsif ( " $defline " =~ / rna_editing=/i ) {
        return 1;
    }
    return 0;
}

sub rna_editing_frameshift {
    my ($ref_id) = @_;

    my $edit = get_rna_edit($ref_id);
    if ( !defined $edit ) { return 0 }

    my ($sz) = split /\//, $edit;
    while ( $sz >= 3 ) { $sz -= 3 }
    while ( $sz <= -3 ) { $sz += 3 }

    return $sz;
}

sub get_rna_edit {
    my ($ref_id) = @_;

    my $ref     = get_reference_seq($ref_id);
    my $defline = $$ref{defline};
    if ( " $defline " =~ / rna_editing=N /i ) {
        return undef;
    }
    elsif ( " $defline " =~ / rna_editing="(.*\/.*\/.*\/.*\/)" / ) {
        return $1;
    }
    elsif ( " $defline " =~ / rna_editing=(.*\/.*\/.*\/.*\/) / ) {
        return $1;
    }

    return undef;
}

sub describe_ribosomal_slippage {
    my ($ref_id) = @_;

    my $slippage_motif = get_parameter("slippage_motif");
    my $slippage_shift = get_parameter("slippage_frameshift");
    my $slippage_offset = get_parameter("slippage_offset");

    my $ref = get_reference_seq($ref_id);
    if ( defined $$ref{defline} ) {
        if ( $$ref{defline} =~ / slippage_motif="([^"]+)"/i ) {
            $slippage_motif = $1;
        }
        if ( "$$ref{defline} " =~ /slippage_frameshift=([^ ]+) /i ) {
            $slippage_shift = $1;
        }
        if ( "$$ref{defline} " =~ /slippage_offset=([^ ]+) /i ) {
            $slippage_offset = $1;
        }
    }
    return ( $slippage_motif, $slippage_shift, $slippage_offset );
}

sub get_startcodons {
    my ($ref_id) = @_;

    my $ref     = get_reference_seq($ref_id);
    my $defline = $$ref{defline};

    my %starts;
    $starts{ATG} = 1;

    if ( " $defline " =~ / alternate_startcodon=N /i ) {
    }
    if ( " $defline " =~ / alternate_startcodon="([^"]+)"/i ) {
        for my $codon ( split /,/, uc $1 ) {
            $starts{$codon} = 1;
        }
    }
    elsif ( " $defline " =~ / alternate_startcodon=([^ ]+) /i ) {
        for my $codon ( split /,/, uc $1 ) {
            $starts{$codon} = 1;
        }
    }

    my @list = keys %starts;
    return \@list;
}

sub is_start_codon {
    my ( $codon, $codonlist ) = @_;

    my $test = $codon;
    $test =~ s/N/./gi;
    for my $start (@$codonlist) {
        if ( $start =~ /$test/i ) {
            if ( $codon =~ /N/ ) {
                return 2;
            }
            else {
                return 1;
            }
        }
    }
    return 0;
}

sub is_stop_codon {
    my ($codon) = @_;

    my @stoplist = ( "TAA", "TAG", "TAR", "TGA" );
    my $test = $codon;
    $test =~ s/N/./gi;
    for my $stop (@stoplist) {
        if ( $stop =~ /$test/i ) {
            if ( $codon =~ /N/ ) {
                return 2;
            }
            else {
                return 1;
            }
        }
    }
    return 0;
}

sub allow_splicing {
    my ($ref_id) = @_;

    my $ref     = get_reference_seq($ref_id);
    my $defline = $$ref{defline};
    if ( " $defline " =~ / spliced=N /i ) {
        return 0;
    }
    elsif ( " $defline " =~ / spliced="{0,1}Y/i ) {
        return 1;
    }
    return 0;
}

sub get_splice_form {
    my ($ref_id) = @_;

    my $ref     = get_reference_seq($ref_id);
    my $defline = $$ref{defline};
    if ( " $defline " =~ / splice_form="([^"]+)"/ ) {
        return $1;
    }
    elsif ( " $defline " =~ / splice_form=([^ ]+) / ) {
        return $1;
    }
    return "";
}

sub get_intron_size {
    my ( $ref_id, $intron_id ) = @_;
    if ( !defined $intron_id ) { $intron_id = 1 }

    my $ref     = get_reference_seq($ref_id);
    my $defline = $$ref{defline};

#    my $intron_range;
#    if ( " $defline " =~ / intron_size="([^"]+)" /i ) {
#
#        #print "get_intron_size1( $ref_id, $1, $2 )\n";
#        $intron_range = $1;
#    }
#    elsif ( " $defline " =~ / intron_size=([^ ]+) /i ) {
#        $intron_range = $1;
#    }
#
#    if ( defined $intron_range ) {
#        my @ranges = split /,/, $intron_range;
#        if ( $intron_id < 1 || $intron_id > @ranges ) {
#            return (
#                get_parameter("min_intron_size"),
#                get_parameter("max_intron_size")
#            );
#        }
#
#        my $range = $ranges[ $intron_id - 1 ];
#        my ( $lo, $hi ) = split /-/, $range;
#        if ( !defined $hi ) {
#            $hi = $lo + 100;
#            $lo = $lo - 100;
#        }
#
#        return ( $lo, $hi );
#    }
#    else {
        my $structure = get_gene_structure($ref_id);
        if (   $intron_id < 1 || !defined $$structure{introns} || $intron_id > @{ $$structure{introns} } ) {
            return ( get_parameter("min_intron_size"), get_parameter("max_intron_size") );
        }

        my @introns = @{ $$structure{introns} };
#print "introns:";
#for my $i ( 0..@introns-1 ) { print "  " . ($i+1) . ": $introns[$i]{intron_size}" }
#print "\n";
        my $size    = $introns[ $intron_id - 1 ]{intron_size};
        my $lo      = int( 0.9 * $size );
        my $hi      = int( 1.1 * $size );
        return ( $lo, $hi );
#    }
}

sub is_required {
    my ($ref_id) = @_;

    my $gene_name = get_reference_name($ref_id);
    if ( exists $required_genes{$gene_name} ) { return 1 }
    if ( exists $optional_genes{$gene_name} ) { return 0 }

    my $ref     = get_reference_seq($ref_id);
    my $defline = $$ref{defline};
    if ( " $defline " =~ / is_required /i ) {
        $required_genes{$gene_name} = 1;
        return 1;
    }
    if ( " $defline " =~ / is_optional /i ) {
        $optional_genes{$gene_name} = 1;
        return 0;
    }

    return get_parameter("default_gene_required");
}

sub is_optional {
    my ($ref_id) = @_;
    my $gene_name = get_reference_name($ref_id);
    if ( exists $optional_genes{$gene_name} ) { return 1 }
    if ( exists $required_genes{$gene_name} ) { return 0 }

    my $ref     = get_reference_seq($ref_id);
    my $defline = $$ref{defline};
    if ( " $defline " =~ / is_optional /i ) {
        $optional_genes{$gene_name} = 1;
        return 1;
    }
    if ( " $defline " =~ / is_required /i ) {
        $required_genes{$gene_name} = 1;
        return 0;
    }

    return get_parameter("default_gene_optional");
}

sub get_tiny_exon3 {
    my ($ref_id) = @_;

    my $ref = get_reference_seq($ref_id);
    if ( !defined $ref ) { return "" }

    my $defline = $$ref{defline};
    if ( " $defline " =~ / tiny_exon3="([^"]+)"/i ) {
        return $1;
    }
    elsif ( " $defline " =~ / tiny_exon3=([^ ]+) /i ) {
        return $1;
    }

    return "";
}

sub get_tiny_exon5 {
    my ($ref_id) = @_;

    my $ref = get_reference_seq($ref_id);
    if ( !defined $ref ) { return "" }

    my $defline = $$ref{defline};
    if ( " $defline " =~ / tiny_exon5="([^"]+)"/i ) {
        return $1;
    }
    elsif ( " $defline " =~ / tiny_exon5=([^ ]+) /i ) {
        return $1;
    }

    return "";
}

sub get_gene_variation {
    my ($ref_id) = @_;

    my $gene_name = get_reference_name($ref_id);
    my $gene_variation;

    if ( exists $gene_variations{$gene_name} ) {
        $gene_variation = $gene_variations{$gene_name};
    }
    else {
        my $ref     = get_reference_seq($ref_id);
        my $defline = $$ref{defline};
        if ( !defined $defline ) {
            print_hash( "ref $ref_id", $ref );
        }
        if ( " $defline " =~ / gene_variation=([0-9]) /i ) {
            $gene_variation = $1;
            $gene_variations{$gene_name} = $gene_variation;
        }
        else {
            $gene_variation = get_parameter("default_gene_variation");
        }
    }
    if    ( $gene_variation < 0 ) { $gene_variation = 0 }
    elsif ( $gene_variation > 5 ) { $gene_variation = 5 }

    return $gene_variation;
}

sub get_variation_penalty {
    my ($ref_id) = @_;

    my $gene_variation = get_gene_variation($ref_id);

    my $variation = get_parameter("variation");

    my $penalty = 1;
    for my $i ( 1 .. $gene_variation ) {
        $penalty = $variation * $penalty;
    }

    return $penalty;
}

sub validate_splicing {
    my ( $splicepairs, $donor, $acceptor ) = @_;

    if ( defined $donor ) {
        my $dtmp = $donor;
        $dtmp =~ s/N/./gi;
        my $splicetype;
        for my $dkey ( keys %$splicepairs ) {
            if ( $dkey =~ /$dtmp/i ) {
                if ( defined $acceptor ) {
                    my $atmp = $acceptor;
                    $atmp =~ s/N/./gi;
                    for my $akey ( keys %{ $$splicepairs{$dkey} } ) {
                        if ( $akey =~ /$atmp/i ) {
                            if ( $acceptor =~ /N/i || $donor =~ /N/i ) {
                                if ( "AG" =~ /$atmp/i && "GT" =~ /$dtmp/i ) {
                                    return ( 3, 0 );
                                }
                                else {
                                    return ( 4, 0 );
                                }
                            }
                            elsif ( $akey =~ /AG/i && $dkey =~ /GT/i ) {
                                return ( 1, $$splicepairs{$dkey}{$akey} );
                            }
                            else {
                                return ( 2, $$splicepairs{$dkey}{$akey} );
                            }
                        }
                    }
                }
                else {
                    if ( $donor =~ /N/i ) {
                        if ( "GT" =~ /$dtmp/i ) {
                            return 3;
                        }
                        else {
                            return 4;
                        }
                    }
                    elsif ( $donor =~ /GT/i ) {
                        return 1;
                    }
                    else {
                        return 2;
                    }
                }
            }
        }
        if ( defined $acceptor ) {
            return ( 0, -5 );
        }
        else {
            return 0;
        }
    }
    if ( defined $acceptor ) {
        my $atmp = $acceptor;
        $atmp =~ s/N/./gi;
        my $best = 0;
        for my $dkey ( keys %$splicepairs ) {
            for my $akey ( keys %{ $$splicepairs{$dkey} } ) {
                if ( $akey =~ /$atmp/i ) {
                    if ( $acceptor =~ /N/i ) {
                        if ( "AG" =~ /$atmp/i ) {
                            return 3;
                        }
                        else {
                            return 4;
                        }
                    }
                    elsif ( $akey =~ /AG/i ) {
                        return 1;
                    }
                    else {
                        $best = 2;
                    }
                }
            }
        }
        return $best;
    }

    return undef;
}

sub get_gene_splicepairs {
    my ($gene) = @_;

    my %gene_splicepairs;
    for my $donor ( keys %splice_pairs ) {
        for my $acceptor ( keys %{ $splice_pairs{$donor} } ) {
            $gene_splicepairs{ uc $donor }{ uc $acceptor } =
              $splice_pairs{$donor}{$acceptor};
        }
    }
    if ( exists $noncanonical_splicing_genes{$gene} ) {
        for my $donor ( keys %{ $noncanonical_splicing_genes{$gene} } ) {
            for my $acceptor (
                keys %{ $noncanonical_splicing_genes{$gene}{$donor} } )
            {
                $gene_splicepairs{ uc $donor }{ uc $acceptor } = 0;
            }
        }
    }

    return \%gene_splicepairs;
}


sub get_matpep_db {
    my ($ref_id) = @_;

    my $ref     = get_reference_seq($ref_id);
    my $defline = $$ref{defline};
    if ( " $defline " =~ / matpepdb=N /i ) {
        return undef;
    }
    elsif ( " $defline " =~ / matpepdb=default /i ) {
        return get_parameter("mature_pep_refdb");
    }
    elsif ( " $defline " =~ / matpepdb="default"/i ) {
        return get_parameter("mature_pep_refdb");
    }
    elsif ( " $defline " =~ / matpepdb="([^"]+)"/i ) {
        my $db = $1;
        if ( $db =~ /^<vigordata>/ ) {
            $db =~ s/^<vigordata>//;
            $db = "$myData/$db";
        }
        return $db;
    }
    elsif ( " $defline " =~ / matpepdb=([^ ]+) /i ) {
        my $db = $1;
        if ( $db =~ /^<vigordata>/ ) {
            $db =~ s/^<vigordata>//;
            $db = "$myData/$db";
        }
        return $db;
    }
    return undef;
}

sub shared_cds {
    my ( $hit1, $hit2 ) = @_;

    my $name1 = get_reference_name( $$hit1{subject_id} );
    my $name2 = get_reference_name( $$hit2{subject_id} );

    if ( $name1 ne $name2 ) {
        my $sym1 = gene_name_to_symbol( $name1 );
        my $sym2 = gene_name_to_symbol( $name2 );
        if ( defined $sym1 && defined $sym2 && $sym1 =~ /\w/ && $sym1 eq $sym2  ) {
            return 1;
        }
    }

    my $ref1 = get_reference_seq( $$hit1{subject_id} );
    my $ref2 = get_reference_seq( $$hit2{subject_id} );
    if ( " $$ref1{defline} " =~ / shared_cds=N /i ) {
        return 0;
    }
    elsif ( " $$ref2{defline} " =~ / shared_cds=N /i ) {
        return 0;
    }
    if ( $$ref1{defline} =~ / shared_cds="([^"]*)"/i ) {
        my $sharedcds = $1;
        my @sharedList = split /,/, $sharedcds;
        for my $shared (@sharedList) {
            if ( $shared eq $name2 ) { return 1 }
        }
    }
    if ( $$ref2{defline} =~ / shared_cds="([^"]*)"/i ) {
        my $sharedcds = $1;
        my @sharedList = split /,/, $sharedcds;
        for my $shared (@sharedList) {
            if ( $shared eq $name1 ) { return 1 }
        }
    }
    return 0;
}

sub make_flavor_reference {
    my ($flavor) = @_;

    my $refname = $flavor;
    my $flavordb;
    my $concatenated = 0;

    my $f = lc($flavor);
    if ( exists $$flavorDB{$f} ) {
        $flavordb = $$flavorDB{$f};
        if ( -e "$flavordb.pal" ) { $concatenated = 1 }
    }
    elsif ( -e "$flavor.pin" ) {
        $flavordb = realpath( $flavor );
    }
    elsif ( -e "$myData/$flavor.pin" ) {
        $flavordb = "$myData/$flavor";
    }
    elsif ( -e "$myData/$flavor.pal" ) {
        $flavordb     = "$myData/$flavor";
        $concatenated = 1;
    }
    else {
        die "\nCannot find reference database $flavor\n";
    }

    if ( !$concatenated ) { setup_reference_db($flavordb) }

    return ( $flavordb, $concatenated );
}

sub choose_reference {
    my ( $genome, $virusdb ) = @_;
    if ( !defined $virusdb ) { $virusdb = "$myData/virus_db" }

    my $dbg = get_parameter("verbose");

    #if ( `whoami` =~ /jhoover/ ) { $dbg = 1 }

    my $vigorspace = get_parameter("vigorspace");
    my $tmpfile    = "$vigorspace/chooseref.fasta";
    my $xmlfile    = "$vigorspace/chooseref.xml";
    my $logfile    = "$vigorspace/chooseref.log";
    open( TMP, ">$tmpfile" );
    my $seq = $$genome{sequence};
    $seq =~ s/(.{1,60})/$1\n/g;
    print TMP ">$$genome{id}\n$seq";
    close TMP;

#    my $blastcmd =
#"$myBin/blastall -p blastx -i $tmpfile -d $virusdb -e 1e-5 -M BLOSUM45 -X 5 -Z 12 -f 14 -F \"\" -v 0 -b 100 -m 7 1> $xmlfile 2> $logfile";
    my $blastcmd = "$myBin/blastall -p blastx -i $tmpfile -d $virusdb -e 1e-5 -M BLOSUM45 -g F -F \"\" -z 3000000 -v 0 -b 100 -m 7 1> $xmlfile 2> $logfile";
    &runCmd($vigorspace, $blastcmd);
    #print "cd $vigorspace; $blastcmd\n";
    #    my @hsps = parse_blastxml( $xmlfile, 0 );
    my @hsps = parse_blastxml( $xmlfile, 0 );

    #print "\nDATABASE HAPS\n";
    #print_blasthits( 0, @hsps );
    unlink $tmpfile;
    unlink $xmlfile;
    unlink $logfile;
    if ( !@hsps ) {
        my $db = get_parameter("reference_db");
        return basename($db);
    }

    my @hits = sort {
        get_reference_name( $$a{subject_id} )
          cmp get_reference_name( $$b{subject_id} )
    } best_subject_hits( $genome, \@hsps );
    if ($dbg) {
        print "\nDATABASE HITS\n";
        print_blasthits( 0, @hits, "refdb" );
    }

    my %ref_db;
    my $bestpct;
    for my $hit ( sort { score_hit($b) <=> score_hit($b) } @hits ) {
        my $gene = get_defline_tag( $$hit{subject_definition}, "gene" );
        my $db   = get_defline_tag( $$hit{subject_definition}, "db" );

        my $wgt = score_hit($hit);
        if ( $wgt < 0 ) { next }
        $wgt = int( 10.0 * $wgt + 0.5 ) / 10.0;
        if ( $wgt == int($wgt) ) { $wgt .= ".0" }

        my $pct =
          ( $$hit{pct_identity} / 100.0 + $$hit{pct_similarity} / 100.0 ) / 2.0;
        $pct = int( 1000.0 * $pct + 0.5 ) / 10.0;
        if ( $pct == int($pct) ) { $pct .= ".0" }

        if ( !exists $ref_db{$db}{$gene}
            || $pct > $ref_db{$db}{$gene}{pct} )
        {

            $ref_db{$db}{$gene}{wgt} = $wgt;
            $ref_db{$db}{$gene}{pct} = $pct;
            if ( !defined $bestpct || $pct > $bestpct ) { $bestpct = $pct }
        }
    }
    my @bestgenes;
    if ($dbg) {
        print rpad( "db", 8 ) . " "
          . rpad( "gene", 8 ) . " "
          . lpad( "pct",    5 ) . " "
          . lpad( "weight", 8 ) . "\n";
    }
    for my $db ( sort { $a cmp $b } keys %ref_db ) {
        for my $gene (
            sort { $ref_db{$db}{$b}{wgt} <=> $ref_db{$db}{$a}{wgt} }
            keys %{ $ref_db{$db} }
          )
        {
            if ($dbg) {
                print rpad( $db, 8 ) . " "
                  . rpad( $gene, 8 ) . " "
                  . lpad( $ref_db{$db}{$gene}{pct}, 5 ) . " "
                  . lpad( $ref_db{$db}{$gene}{wgt}, 8 ) . "\n";
            }
            if ( $ref_db{$db}{$gene}{pct} > 0.80 * $bestpct ) {
                my $best;
                $$best{db}   = $db;
                $$best{gene} = $gene;
                $$best{wgt}  = $ref_db{$db}{$gene}{wgt};
                $$best{pct}  = $ref_db{$db}{$gene}{pct};
                push @bestgenes, $best;
            }
        }
    }

    my $newpath;
    if ( !@bestgenes ) {
        $newpath = $virusdb;
    }
    else {
        my %dbstats;
        for my $best ( sort { $$b{wgt} <=> $$a{wgt} } @bestgenes ) {
            if ( !exists $dbstats{ $$best{db} } ) {
                $dbstats{ $$best{db} }{db}    = $$best{db};
                $dbstats{ $$best{db} }{scale} = 1.0;
                $dbstats{ $$best{db} }{wgt}   = $$best{wgt};
            }
            else {
                $dbstats{ $$best{db} }{scale} =
                  $dbstats{ $$best{db} }{scale} / 2.0;
                $dbstats{ $$best{db} }{wgt} +=
                  $dbstats{ $$best{db} }{scale} * $$best{wgt};
            }
        }

        my @tmp = sort { $$b{wgt} <=> $$a{wgt} } values %dbstats;
        if ($dbg) {
            print rpad( "db", 12 ) . " " . lpad( "wgt", 7 ) . "\n";
            for my $db (@tmp) {
                my $wgt = int( 100.0 * $$db{wgt} + 0.5 ) / 100.0;
                if ( $wgt == int($wgt) ) { $wgt .= ".0" }
                print rpad( $$db{db}, 12 ) . " " . lpad( $wgt, 7 ) . "\n";
            }
        }
        $newpath = "$myData/$tmp[0]{db}";
    }

    my $curpath = get_parameter("reference_db");
    my $newdb   = basename($newpath);
    if ( $newpath ne $curpath ) { setup_reference_db($newpath) }

    return $newdb;
}


sub setup_reference_db {
    my ($newpath) = @_;
    my $newdb = basename($newpath);

    my %allrefs = %{ get_reference_seqs() };
    initialize_defaults();
    set_parameter( "reference_db", $newpath );
    load_reference_db();
    set_reference_seqs(%allrefs);

    for my $config (@$configData) {
        if ( $$config{db} eq $newdb && $$config{type} eq "parameter" ) {
            set_parameter( $$config{col1}, $$config{col2} );
        }
    }
}

sub get_defline_tag {
    my ( $defline, $tag ) = @_;

    if ( " $defline " =~ / $tag="([^"]*)"/ ) {
        return $1;
    }
    elsif ( " $defline " =~ / $tag=([^ ]*) / ) {
        return $1;
    }
    else { return "" }
}

sub score_hit {
    my ($hit) = @_;

    my $match    = 0.5 * $$hit{num_identical} + 0.5 * $$hit{num_similar};
    my $mismatch =
      maxval( 0, $$hit{subject_end} - $$hit{subject_begin} + 1 - $match );
    my $score = $match - 1.5 * $mismatch;

    return $score;
}

sub split_gene_at_gaps {
    my ( $genome, $gene ) = @_;
    my $dbg = DEBUG;
    #if ( $$gene{gene_name} =~ /LMP2A/ ) { $dbg = 1 }

    if ( $$gene{is_pseudogene} ) {
        $$gene{product_name} = modify_product_name( $$gene{product_name}, $$gene{start_truncation}, $$gene{stop_truncation} );    
        return ($gene);
    }

    my @chunks;
    my $skipgap     = 0;
    my $codon_start = $$gene{codon_start};

    my $chunk;
    $$chunk{dna_begin}        = $$gene{start_site};
    $$chunk{codon_start}      = $codon_start;
    $$chunk{start_truncation} = $$gene{start_truncation};
    $$chunk{stop_truncation}  = 1;

    for my $exon ( @{ $$gene{exons} } ) {
        my $i = $$exon{dna_begin};
        while ($$gene{orientation} * $i <= $$gene{orientation} * $$exon{dna_end} ) {
            if ( $dbg ) { 
                print "i=$i  skipgap=$skipgap  codon_start=$codon_start is_gap=" 
            }
            if ( defined $$genome{in_gap}{$i} ) {

                if ( $dbg ) {
                     print "Y\n" 
                }
                if ( !$skipgap ) {
                    my %tmp = %$chunk;
                    push @chunks, \%tmp;

                    if ( $dbg ) { 
                        print_hash( "chunk", \%tmp ) 
                    }
                    $skipgap = 1;
                }
            }
            else {
                if ( $dbg ) { 
                    print "N\n" 
                }
                if ($skipgap) {
                    $$chunk{dna_begin}        = $i;
                    $$chunk{codon_start}      = $codon_start;
                    $$chunk{start_truncation} = 1;
                    $$chunk{stop_truncation}  = 1;
                    $skipgap                  = 0;
                }
                $$chunk{dna_end} = $i;
            }
            $i += $$gene{orientation};
            $codon_start--; # ?!
            
            if ( defined($$exon{rna_edit}) && $i == $$exon{rna_edit}{dna_begin} ) {
                $i           = $$exon{rna_edit}{dna_end};
                $codon_start = 3;
            }
            if ( $codon_start <= 0 ) { 
                $codon_start += 3 
            }
        }
    }
    if ( !@chunks ) {
        $$gene{product_name} = modify_product_name( $$gene{product_name}, $$gene{start_truncation}, $$gene{stop_truncation} );    
        return ($gene);
    }

    if ( !$skipgap ) {
        my %tmp = %$chunk;
        $tmp{stop_truncation} = $$gene{stop_truncation};
        push @chunks, \%tmp;

        if ( $dbg ) { print_hash( "chunk", \%tmp ) }
    }

    my $fragid = 0;
    my @frags;
    for my $chunk (@chunks) {
        $fragid++;
        $$chunk{fragid} = $fragid;
        my $frag = subgene( $genome, $gene, $chunk );
        push @frags, $frag;

        if ( $dbg ) {
            print_hash( "chunk", $chunk );
            print_genehits( $frag );
        }
    }

    if ( get_parameter("verbose") ) {
        if ( $dbg ) {
            print "GENE->FRAGS\n";
            print_genehits( $gene, @frags );
        }
    }
    return @frags;
}

sub modify_product_name {
    my ( $product_name, $trunc5, $trunc3 ) = @_;

    my $putative_name = $product_name;
    if ( $putative_name !~ /putative/i ) {
        $putative_name = "putative $putative_name";
    }

    if ($trunc5) {
        if ($trunc3) {
            return "$putative_name, fragment";
        }
        else {
            return "$putative_name, C-terminal";
        }
    }
    elsif ($trunc3) {
        return "$putative_name, N-terminal";
    }
    else {
        return $product_name;
    }
}

sub subgene {
    my ( $genome, $gene, $chunk ) = @_;
    my $dbg = DEBUG;

    my $sub;
    my $suffix = substr( "abcdefghijklmnopqrstuvwxyz", $$chunk{fragid} - 1, 1 );
    $$sub{gene_id} = "$$gene{gene_id}$suffix";
    $$sub{is_pseudogene} = $$gene{is_pseudogene};
    $$sub{genome} = $$gene{genome};
    $$sub{orientation} = $$gene{orientation};
    $$sub{start_site} = $$chunk{dna_begin};
    $$sub{codon_start} = $$chunk{codon_start};
    if ( $$sub{start_site} == $$gene{start_site} ) {
        $$sub{start_truncation} = $$gene{start_truncation};
        if ( defined $$gene{alternate_startcodon} ) {
            $$sub{alternate_startcodon} = $$gene{alternate_startcodon};
        }
    }
    else {
        $$sub{start_truncation} = 1;
    }

    $$sub{stop_site} = $$chunk{dna_end};
    if ( $$sub{stop_site} == $$gene{stop_site} ) {
        $$sub{stop_truncation} = $$gene{stop_truncation};
        if ( defined $$gene{alternate_stopcodon} ) {
            $$sub{alternate_stopcodon} = $$gene{alternate_stopcodon};
        }
    }
    else {
        $$sub{stop_truncation} = 1;
    }

    my @exons;

    if ( $dbg ) { print_hash("chunk", $chunk ) }
    for my $exon ( @{ $$gene{exons} } ) {

        if ( $dbg ) { print_hash( "  exon in", $exon ) }
        if ( $$gene{orientation} * $$exon{dna_end} <
            $$gene{orientation} * $$chunk{dna_begin} )
        {
            next;
        }
        if ( $$gene{orientation} * $$exon{dna_begin} >
            $$gene{orientation} * $$chunk{dna_end} )
        {
            last;
        }

        my %tmp = %$exon;
        if ( $$gene{orientation} * $$exon{dna_begin} <
            $$gene{orientation} * $$chunk{dna_begin} )
        {
            $tmp{subject_begin} +=
              int( abs( $$chunk{dna_begin} - $$exon{dna_begin} ) / 3 );
            $tmp{dna_begin} = $$chunk{dna_begin};
        }
        if ( $$gene{orientation} * $$exon{dna_end} >
            $$gene{orientation} * $$chunk{dna_end} )
        {
            $tmp{subject_end} -=
              int( abs( $$exon{dna_end} - $$chunk{dna_end} ) / 3 );
            $tmp{dna_end} = $$chunk{dna_end};
        }
        if ( defined $$exon{rna_edit} ) {
            if ( $$gene{orientation} * $$exon{rna_edit}{dna_begin} >
                   $$gene{orientation} * $$chunk{dna_end}
                || $$gene{orientation} * $$exon{rna_edit}{dna_end} <
                $$gene{orientation} * $$chunk{dna_begin} )
            {
                delete $tmp{rna_edit};
            }
        }

        if ( $dbg ) { print_hash("    exon out", \%tmp) }
        push @exons, \%tmp;
    }
    $exons[0]{fuzzy_begin} = $$sub{start_truncation};
    $exons[ @exons - 1 ]{fuzzy_end} = $$sub{stop_truncation};

    set_cdna_coordinates(@exons);
    $$sub{exons} = \@exons;
    $$sub{cdna}  = cdna_from_exons( $genome, @exons );

    if ( $dbg ) { print "cdna start $$sub{codon_start}  seq: $$sub{cdna}\n" }

    $$sub{protein} =
      DNA2AA( $$sub{cdna}, $$sub{codon_start}, $$sub{alternate_startcodon},
        undef, $$sub{alternate_stopcodon} );

    if ( $dbg ) { print "protein $$sub{protein}\n" }
    if ( defined $$gene{stopcodon_readthru}
        && $$gene{orientation} * $$gene{stopcodon_readthru}{dna_begin} >=
        $$gene{orientation} * $$chunk{dna_begin}
        && $$gene{orientation} * $$gene{stopcodon_readthru}{dna_end} <=
        $$gene{orientation} * $$chunk{dna_end} )
    {
        my $aa_pos = index( $$sub{protein}, "*" ) + 1;
        if ( $aa_pos > 0 ) {
            $$sub{stopcodon_readthru}         = $$gene{stopcodon_readthru};
            $$sub{stopcodon_readthru}{aa_pos} = $aa_pos;
            $$sub{protein}                    =
                substr( $$sub{protein}, 0, $aa_pos - 1 )
              . $$sub{stopcodon_readthru}{aa}
              . substr( $$sub{protein}, $aa_pos );

            if ( $dbg ) { print "pos $aa_pos  aa $$sub{stopcodon_readthru}{aa}  protein $$sub{protein}\n" }
        }
    }
    $$sub{protein_length} = length( $$sub{protein} );
    if ( ! $$sub{stop_truncation} ) { $$sub{protein_length}-- }

    $$sub{ref_id}     = $$gene{ref_id};
    $$sub{ref_length} = $$gene{ref_length};
    $$sub{ref_begin}  = $exons[0]{subject_begin};
    $$sub{ref_end}    = $exons[ @exons - 1 ]{subject_end};

    $$sub{gene_name}    = $$gene{gene_name};
    $$sub{gene_synonym} = $$gene{gene_synonym};
    $$sub{product_name} = modify_product_name( $$gene{product_name}, $$sub{start_truncation}, $$sub{stop_truncation} );
    $$sub{note}         = $$gene{note};
    if ( $$sub{start_site} != $$gene{start_site} || $$sub{stop_site} != $$gene{stop_site} ) {
        if ( defined $$sub{note} ) {
            $$sub{note} .= "; coding region disrupted by sequencing gap";
        }
        else {
            $$sub{note} = "coding region disrupted by sequencing gap";
        }
    }

    if ( defined $$gene{cdsnotes} && @{ $$gene{cdsnotes} } ) {
        my @cdsnotes;
        for my $note ( $$gene{cdsnotes} ) {
            if ( $note =~ /start/i && !$$sub{start_truncation} ) {
                push @cdsnotes, $note;
            }
            elsif ( $note =~ /stop/i && !$$sub{stop_truncation} ) {
                push @cdsnotes, $note;
            }
        }
        if (@cdsnotes) {
            $$sub{cdsnotes} = \@cdsnotes;
        }
    }

    if ( defined $$gene{noncanonical_splicing}
        && @{ $$gene{noncanonical_splicing} } )
    {
        my @noncanon;
        for my $note ( @{ $$gene{noncanonical_splicing} } ) {
            if ( $note =~ /([0-9]+)\.\.([0-9]+)/ ) {
                my ( $begin, $end ) = ( $1, $2 );
                if ( $$gene{orientation} * $begin >=
                       $$gene{orientation} * $$chunk{dna_begin}
                    && $$gene{orientation} * $end <=
                    $$gene{orientation} * $$chunk{dna_end} )
                {
                    push @noncanon, $note;
                }
            }
        }
        if (@noncanon) {
            $$sub{noncanonical_splicing} = \@noncanon;
        }
    }

    score_gene( $genome, $sub );

    if ( $dbg ) {
        print "\n\n";
        print_genehits( $gene, $sub );
        print_hash( "sub", $sub );
    }
    return $sub;
}

sub runCmd {
    my ($work_dir, $cmd, $ignore_crashes) = @_;
    my $current_dir = getcwd();
    my $success = 1;
    
    unless ($work_dir eq '.' || $work_dir eq $current_dir) {
        if (-d $work_dir) {
            chdir($work_dir) || die "\n\nImpossible to change directory to \"$work_dir\".\n\n";
        }
        else {
            die "\n\nImpossible to find working directory \"$work_dir\".\n\n";
        }
    }
    if (system($cmd)) {
        if ($ignore_crashes) {
            warn "\n\nProblems executing the following command: \"$cmd\" (Error Number: \"$!\" - Child error: \"$?\")\n\n";
            $success = 0;
        }
        else {
            die "\n\nProblems executing the following command: \"$cmd\" (Error Number: \"$!\" - Child error: \"$?\")\n\n";
        }
    }
    unless ($work_dir eq '.' || $work_dir eq $current_dir) {
        chdir($current_dir) || die "\n\nImpossible to change directory back to \"$current_dir\".\n\n";
    }
    return($success);
}

sub runCmdAndGetResults {
    my ($work_dir, $cmd, $ignore_crashes) = @_;
    my $current_dir = getcwd();
    my $results = undef;
    
    unless ($work_dir eq '.' || $work_dir eq $current_dir) {
        if (-d $work_dir) {
            chdir($work_dir) || die "\n\nImpossible to change directory to \"$work_dir\".\n\n";
        }
        else {
            die "\n\nImpossible to find working directory \"$work_dir\".\n\n";
        }
    }
    $results = `$cmd`;
    
    if($! || $?) {
        if ($ignore_crashes) {
            warn "\n\nProblems executing the following command: \"$cmd\" (Error Number: \"$!\" - Child error: \"$?\")\n\n";
        }
        else {
            die "\n\nProblems executing the following command: \"$cmd\" (Error Number: \"$!\" - Child error: \"$?\")\n\n";
        }
    }
    unless ($work_dir eq '.' || $work_dir eq $current_dir) {
        chdir($current_dir) || die "\n\nImpossible to change directory back to \"$current_dir\".\n\n";
    }
    return($results);
}
1;
