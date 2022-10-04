#!/usr/bin/env perl 

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
use Getopt::Std;
use Cwd ('realpath', 'abs_path');
use File::Basename;
use File::Path;

my $program = realpath($0);
our $myBin = dirname($program);
our $myData = dirname ( $myBin ) . "/data3";
our $myConf = "$myBin/conf";
our $flavorData;

require "$myBin/VIGOR3.pm";
require "$myBin/db.pm";
*STDERR = *STDOUT;
$|++;

# initialize VIGOR
initialize_defaults();
set_parameter( "command_line", $0 . " " . join( " ", @ARGV ) );
my $vigorspace = create_workspace();

load_config_db();

# get command line parameters
my %args;
&getopts( 'aAc:Cd:D:e:f:G:hi:I:jlLmo:O:P:s:vx:0', \%args );
my ( $opt_a,   $opt_A,   $opt_c,   $opt_C,   $opt_d,   $opt_D,   $opt_e,   $opt_f,   $opt_G,   $opt_h,   $opt_i,   $opt_I,   $opt_j,   $opt_l,   $opt_L,   $opt_m,   $opt_o,   $opt_O,   $opt_P,   $opt_s,   $opt_v,   $opt_x,   $opt_0 ) =
   ( $args{a}, $args{A}, $args{c}, $args{C}, $args{d}, $args{D}, $args{e}, $args{f}, $args{G}, $args{h}, $args{i}, $args{I}, $args{j}, $args{l}, $args{L}, $args{m}, $args{o}, $args{O}, $args{P}, $args{s}, $args{v}, $args{x}, $args{0} );

if ( $opt_h ) {
    print_help();
    exit(0);
}

if ( ! defined $opt_i ) {
    if ( defined $opt_I ) {
        $opt_i = $opt_I;
    }
    else {
        die "\nNo input file specified (-i)\n";
    }
}
elsif ( defined $opt_I && $opt_i ne $opt_I ) {
    die "\nAmbigous input fasta specified -i \"$opt_i\" and -I \"$opt_I\"\n";
}

if ( ! -r $opt_i ) {
    die "\nCannot read input file specified (-i $opt_i)\n";
}

if ( ! defined $opt_O ) {
    if ( defined $opt_o ) {
        $opt_O = $opt_o;
    }
    else {
        die "\nOutput prefix is required (-O)\n";
    }
}
elsif ( defined $opt_o && $opt_o ne $opt_O ) {
    die "\nAmbigous output prefix specified -O \"$opt_O\" and -o \"$opt_o\"\n";
}

if ( ! defined $opt_v ) { $opt_v = 0 }
my $verbose = $opt_v;
set_parameter( "verbose", $verbose );

# default to -A (autoselect genome) if no genbank file (-G) or refernce database (-D) specified
if ( ! defined $opt_D ) {
    if ( defined $opt_d ) {
        $opt_D = $opt_d;
    }
}
elsif ( defined $opt_d && $opt_d ne $opt_D ) {
    die "\nAmbigous reference database specified -d \"$opt_d\" and -D \"$opt_D\"\n";
}

if ( ! defined $opt_A ) { $opt_A = $opt_a }

if ( ! defined $opt_G && ! defined $opt_D && ! defined $opt_A ) {
    $opt_A = 1;
}

my $autoselect = 0;
if ( $opt_A ) {
    set_parameter( "reference_db", "$myData/virus_db" );
    $autoselect = 1;
}

# use user specified reference database
if ( defined $opt_D ) {
    ( my $flavordb, $autoselect ) = make_flavor_reference( $opt_D );
    set_parameter( "reference_db", $flavordb );
}

# use genbank file as reference database
elsif ( defined $opt_G ) {
    if ( $opt_G !~ /\.gbk/i ) { die "-G $opt_G must specify a genbank file.\n" }
    set_parameter( "reference_db", $opt_G );
    set_parameter( "mature_pep_refdb", $opt_G );
    set_parameter( "default_gene_required", 1 );
    set_parameter( "candidate_blastopts", "-p blastx -G 8 -E 2 -M BLOSUM80 -e 1e-5 -F \"\" -g F -Z 3000000" );
    set_parameter( "candidate_overlap_threshold", 1.0 );
}

# is genome complete sequence?
my $genome_is_complete = 0;
my $genome_is_circular;
if ( defined $opt_C ) {
    $genome_is_complete = 1;
}
if ( defined $opt_0 ) {
    $genome_is_complete = 1;
    $genome_is_circular = 1;
}

# gene coverage requirements (in TBL report)
my $min_gene_size;
my $min_gene_coverage;
if ( defined $opt_s ) {
    $min_gene_size = $opt_s;
}
if ( defined $opt_c ) {
    $min_gene_coverage = $opt_c;
}

# override frameshift sesnsitivity
if ( defined $opt_f ) {
    if ( $opt_f == 0 || $opt_f == 1 || $opt_f == 2 ) { set_parameter( "frameshift_sensitivity", $opt_f ) }
    else { die "\n invalid value for framehshift sensitivity (-f $opt_f), expected 0, 1, or 2\n" }
}

# turn of locus tags
if ( defined $opt_l ) {
    if ( defined $opt_L ) { die "\n incompatible switches, use -l or -L but not both\n" }
    set_parameter( "use_locus_tags", 0 );
}
elsif ( defined $opt_L ) {
    set_parameter( "use_locus_tags", 1 );
}

# command line parameter overrides
if ( defined $opt_m ) {
    turnoff_refmatch_requirements();
}
if ( defined $opt_P ) {
    parse_parameters( $opt_P );
}

if ( defined $opt_e ) {
    if ( $opt_e > 0 ) {
        set_parameter( "candidate_evalue", $opt_e );
    }
    else {
        die "\nE-value must be > 0 (-e $opt_e)\n";
    }
}
if ( defined $opt_j ) {
    set_parameter( "jcvi_rules", 0 );
}

# open report file (log file for run)
open( my $RPT, ">$opt_O.rpt" ) || die "cannot write to the report file \"$opt_O.rpt\"\.$!\n";

# display parameters
print now() . "\n";
print $RPT now() . "\n\n";
show_parameters( $RPT ); 

# load reference sequences
load_reference_db();

# open input and output files
open(my $INPUT, "<$opt_i" ) || die "cannot open the file -i \"$opt_i\"\.$!\n";

open( my $ALIGN, ">$opt_O.aln" ) || die "cannot write to the alignment report \"$opt_O.aln\"\.$!\n";
open( my $FS, ">$opt_O.fs" ) || die "cannot write to the assembly report \"$opt_O.fs\"\.$!\n";
open( my $AT, ">$opt_O.at" ) || die "cannot write to the auto-tasker file \"$opt_O.asm\"\.$!\n";
open( my $TBL, ">$opt_O.tbl" ) || die "cannot write to the TBL file \"$opt_O.tbl\"\.$!\n";
open( my $CDS, ">$opt_O.cds" ) || die "cannot write to the CDS fasta \"$opt_O.cds\"\.$!\n";
open( my $PEP, ">$opt_O.pep" ) || die "cannot write to the protein fasta \"$opt_O.pep\"\.$!\n";
open( my $STATS, ">$opt_O.stats" ) || die "cannot write to the stats file \"$opt_O.stats\"\.$!\n";
stats_header( $STATS, 1 );

# for each genomic sequence
my $reference_db = get_parameter( "reference_db" );
while ( my $genome = next_sequence( $INPUT ) ) {
    
    my @genes;
    my $gene_num = 0;
    my $fragcount = 0;
    print "\ngenomic sequence $$genome{id}, length=". length( $$genome{sequence} ) ."\n";

    $genome->{is_complete} = $genome_is_complete;
    $genome->{is_circular} = $genome_is_circular;
    if ( $genome->{is_circular} ) {
        $genome->{original_seqlen} = $genome->{seqlen};
        $genome->{sequence} .= $genome->{sequence};
        $genome->{seqlen} += $genome->{seqlen};
    }

    if ( $autoselect ) {
        my $newdb = choose_reference( $genome, $reference_db );
        print "  selecting database " . basename( $newdb ) . "\n";
    
        if ( defined $opt_m ) {
            turnoff_refmatch_requirements();
        }
        if ( defined $opt_P ) {
            parse_parameters( $opt_P );
        }
        if ( defined $opt_e ) {
            set_parameter( "candidate_evalue", $opt_e );
        }
        if ( defined $opt_j ) {
            set_parameter( "jcvi_rules", 0 );
        }
        if ( defined $opt_l ) {
            set_parameter( "use_locus_tags", 0 );
        }
        elsif ( defined $opt_L ) {
            set_parameter( "use_locus_tags", 1 );
        }
        if ( defined $opt_f ) {
            set_parameter( "frameshift_sensitivity", $opt_f );
        }
    }

    # find candidate genes
    my @candidates = find_candidates( $genome, $opt_x );
    
    # process each candidate gene
    for my $candidate ( @candidates ) {

        if ( $verbose ) {
            print "===========================================================\n";
        }

        # fill in missing pieces of gene
        # (missing exons, gaps in the alignment)
        fill_in_candidate_gaps( $genome, $candidate );
        my @permutations = @{ $$candidate{permutations} }; 

        if ( @{$$candidate{permutations}} > 1 ) {
            print "\n    found " . @permutations . " permutations for gene $$candidate{gene_name}\n";
        }
        else {
            print "\n    found " . @permutations . " permutation for gene $$candidate{gene_name}\n";
        }
        print_blasthits( 4, @permutations );

        # find best permutation of gene
        $gene_num++;
        my @selected;
        my $topquality;
        my $permuteid = 0;
        for my $permutation ( @permutations ) {
            $permuteid++;
            if ( $verbose ) {
                print "\n-------------------------------------------------------\n"
                    . "PERMUTATION\n";
                print_blasthits( 0,  $permutation, @{$$permutation{hsps}} );
            }
            $$permutation{gene_num} = $gene_num;
            my $tmpgene = find_gene( $genome, $permutation );
            if ( defined $min_gene_size && $$tmpgene{protein_length} >= $min_gene_size ) {
#print "size $$tmpgene{protein_length} >= $min_gene_size\n";
            }
            elsif ( defined $min_gene_coverage && $$tmpgene{pct_refcoverage} >= $min_gene_coverage ) {
#print "coverage $$tmpgene{pct_refcoverage} >= $min_gene_coverage\n";
            }
            elsif ( defined $min_gene_size || defined $min_gene_coverage ) {
                next;
            }

            if ( $verbose && @permutations > 1 ) {
                print "\nMATCH QUALITY = $$tmpgene{match_quality}  PSEUDOGENE $$tmpgene{is_pseudogene}  GENE QUALITY = $$tmpgene{gene_quality}\n";
                print_genehits( $tmpgene );
            }

            if ( ! defined $topquality ) {
                @selected = ( $tmpgene );
                $topquality = $$tmpgene{gene_quality};
            }
            elsif ( $selected[0]{is_pseudogene} && ! $$tmpgene{is_pseudogene} ) {
                @selected = ( $tmpgene );
                $topquality = $$tmpgene{gene_quality};
            }
            elsif ( ! $selected[0]{is_pseudogene} && $$tmpgene{is_pseudogene} ) {
            }
            elsif ( $$tmpgene{gene_quality} > $topquality ) {
                @selected = ( $tmpgene );
                $topquality = $$tmpgene{gene_quality};
            }
            elsif ( $$tmpgene{gene_quality} == $topquality ) {
                if ( get_gene_location( $tmpgene ) ne get_gene_location( $selected[0] ) ) {
                    push @selected, $tmpgene;
                }
                else {
                }
            }
        }
            
        # final cut
        if ( @selected > 1 ) {
            resolve_ambiguous_splicing( \@selected );                
        }
        if ( @selected == 1 ) {
            print "        found " . @selected . " gene\n";
        }
        elsif ( @selected > 1 ) {
            print "        found " . @selected . " genes\n";
        }
        else {
            $gene_num--;
            print "        found no genes\n";
            print "===========================================================\n";
            next;
        }

        # split genes at gaps (JCVI rules)
        my $jcvi_rules = get_parameter( "jcvi_rules" );
        if ( $jcvi_rules ) {
            my @tmp;
            for my $selected ( @selected ) {
                my @frags = split_gene_at_gaps( $genome, $selected );
                for my $frag ( @frags ) {
                    if ( defined $min_gene_size && $$frag{protein_length} >= $min_gene_size ) {
                    }
                    elsif ( defined $min_gene_coverage && $$frag{pct_refcoverage} >= $min_gene_coverage ) {
                    }
                    elsif ( defined $min_gene_size || defined $min_gene_coverage ) {
                        next;
                    }
                    push @tmp, $frag;                    
                }
            }
            @selected = @tmp;
        }
    
        # annotate mature peptides            
        my $printmsg = 1;
        for my $selected ( @selected ) {
            my $mpdb = get_matpep_db( $$selected{ref_id} );
            
            if ( defined $mpdb ) {
                $$selected{mature_peps} = mapto_polyprotein($selected, $$selected{protein}, $mpdb, $min_gene_size, $min_gene_coverage );
                
                if ( $jcvi_rules && defined $$selected{mature_peps} && @{ $$selected{mature_peps} } ) {
                    for my $pep ( @{ $$selected{mature_peps} } ) {
                        $$pep{product_name} = modify_product_name( $$pep{product_name}, $$pep{fuzzy_begin}, $$pep{fuzzy_end} );
                    }
                }
            }
        }

        if ( $verbose ) {
            print "\nSELECTED PERMUTATIONS\n";
            print_genehits( @selected );
            print "===========================================================\n";
        }
        
        push @genes, @selected;
    }

    # clean up final gene set
    cleanup_genes( $RPT, $genome, \@genes );

    # print summary of results
    print "\n";
    print $RPT "\n";
    stats_header( *STDOUT, 0 );
    stats_report( *STDOUT, 0, $genome, \@genes );
    stats_header( $RPT, 0 );
    stats_report( $RPT, 0, $genome, \@genes );
    stats_report( $STATS, 1, $genome, \@genes );
    print "\n";
    print $RPT "\n";
    print_genehits( @genes );
    {
        my $STDOUT = *STDOUT;
        *STDOUT = $RPT;
        print_genehits( @genes );
        *STDOUT = $STDOUT;    
    }
    print "\n";
    print $RPT "\n";

    # produce standard reports
    std_tbl_report( $TBL, $genome, \@genes );
    std_cds_report( $CDS, $genome, \@genes );
    std_pep_report( $PEP, $genome, \@genes );
    std_align_report( $ALIGN, $FS, $genome, \@genes );
    std_autotasker_report( $AT, $genome, \@genes );
}

# close files and clean up workspace
if ( defined $INPUT ) { close $INPUT }
if ( defined $STATS ) { close $STATS }
if ( defined $TBL ) { close $TBL }
if ( defined $CDS ) { close $CDS }
if ( defined $PEP ) { close $PEP }
if ( defined $ALIGN ) { close $ALIGN }
if ( defined $FS ) { close $FS }
if ( defined $AT ) { close $AT }
close $RPT;

if ( $vigorspace =~ /\/[0-9]+\.[0-9]+$/ ) {
    system "rm -rf $vigorspace";
}
exit(0);

sub print_help {
    print "Usage:\n"
        . "  -- allow VIGOR to choose the reference database\n"
        . "  \$./VIGOR3.pl -i inputfasta -o outputprefix\n"
        . "\n"
        . "  -- tell VIGOR which reference database to use\n"
        . "  \$./VIGOR3.pl -d refdb -i inputfasta -o outputprefix\n"
        . "\n";
        
    print "Command Line Options\n"
        . "  -a auto-select the reference database, equivalent to \"-d any\", default behavior unless\n"
        . "      overridden by -d or -G, (-A is a synonym for this option)\n"
        . "  -d <ref db>, specify the reference database to be used, (-D is a synonym for this option)\n"
        . "  -e <evalue>, override the default evalue used to identify potential genes, the default\n"
        . "     is usually 1E-5, but varies by reference database\n"
        . "  -c <pct ref> minimum coverage of reference product (0-100) required to report a gene, by\n"
        . "     default coverage is ignored\n"
        . "  -C complete (linear) genome (do not treat edges as gaps)\n"
        . "  -0 complete circular genome (allows gene to span origin)\n"
        . "  -f <0, 1, or 2>, frameshift sensitivity, 0=ignore frameshifts, 1=normal (default), 2=sensitive\n"
        . "  -G <genbank file>, use a genbank file as the reference database, caution: VIGOR genbank\n"
        . "     parsing is fairly rudimentary and many genbank files are unparseable.  Partial genes will\n"
        . "     be ignored. Note: genbank files do not record enough information to handle RNA editing\n"
        . "  -i <input fasta>, path to fasta file of genomic sequences to be annotated, (-I is a synonym\n"
        . "      for this option)\n"
        . "  -l do NOT use locus_tags in TBL file output (incompatible with -L)\n"
        . "  -L USE locus_tags in TBL file output (incompatible with -l)\n"
        . "  -o <output prefix>, prefix for outputfile files, e.g. if the ouput prefix is /mydir/anno\n"
        . "     VIGOR will create output files /mydir/anno.tbl, /mydir/anno.stats, etc., (-O is a synonym\n"
        . "     for this option)\n"
        . "  -P <parameter=value~~...~~paramaeter=value>, override default values of VIGOR parameters\n"
        . "  -j turn off JCVI rules, JCVI rules treat gaps and ambiguity codes conservatively, use\n"
        . "     this option to relax these constraints and produce a more speculative annotation\n"
        . "  -m ignore reference match requirements (coverage/identity/similarity), sometimes useful\n"
        . "     when running VIGOR to evaluate raw contigs and rough draft sequences\n"
        . "  -s <gene size> minimum size (aa) of product required to report a gene, by default size is\n"
        . "     ignored\n"
        . "  -x <ref_id,...,ref_id> comma separated list of reference sequence IDs to ignore\n"
        . "     (useful when debugging a reference database)\n"
        . "\n";

    print "Outputs:\n"
        . "  outputprefix.rpt - summary of program results\n"
        . "  outputprefix.stats - run statistics (per genome sequence) in tab-delimited format\n"
        . "  outputprefix.cds - fasta file of predicted CDSs\n"
        . "  outputprefix.pep - fasta file of predicted proteins\n"
        . "  outputprefix.tbl - predicted features in GenBank tbl format\n"
        . "  outputprefix.aln - alignment of predicted protein to reference, and reference protein to genome\n"
        . "  outputprefix.fs - subset of aln report for those genes with potential sequencing issues\n"
        . "  outputprefix.at - potential sequencing issues in tab-delimited format\n"
        . "\n";

    print "Reference Databases:\n";
    {
        my %dbs;
        for my $flavor ( sort { $$a{flavor} cmp $$b{flavor} } @$flavorData ) {
            if ( defined $$flavor{description} && $$flavor{description} !~ /^ *$/ ) {
                $$flavor{show} = 1;
                if ( ! exists $dbs{$$flavor{db}} ) { $dbs{$$flavor{db}} = $flavor }
            }
            else {
                $$flavor{show} = 0;
                $$flavor{description} = "(no description provided)";
            }
        }

        for my $flavor ( sort { $$a{flavor} cmp $$b{flavor} } @$flavorData ) {
            if ( ! $$flavor{show} ) {
                if ( exists $dbs{$$flavor{db}} ) {
                    my $parent = $dbs{$$flavor{db}};
                    if ( exists $$parent{synonyms} ) {
                        $$parent{synonyms} .= ", $$flavor{flavor}";
                    }
                    else {
                        $$parent{synonyms} = $$flavor{flavor};
                    }
                }
                else {
                    $$flavor{show} = 1;
                }
            }
        }
    }
    my $note = 0;
    print "  " . rpad( "Name", 12 ) . rpad( "Description", 45 ) . " (Synonyms)\n";
    for my $flavor ( sort { $$a{description} cmp $$b{description} } @$flavorData ) {
        if ( $$flavor{show} ) {
            my $synonyms = $$flavor{synonyms};
            if ( ! defined $synonyms ) { $synonyms = "" }
            my $description = $$flavor{description};
            $description =~ s/^ +//;
            if ( $description =~ /\*$/ ) { $note = 1 }
            if ( $synonyms =~ /\w/ ) { $synonyms = "($synonyms)" }
            print "  " . rpad( $$flavor{flavor}, 12 ) . rpad( $description, 45 ) . " $synonyms\n";
        }
    }
    if ( $note ) {
        print "\n* non-standard grouping, must be invoked directly, not included in \"any virus\" via -A or as a\n"
            . "  subset of other -D specifications\n";
    }
}
