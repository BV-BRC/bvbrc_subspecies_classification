#!/usr/local/bin/perl -w

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
use Cwd 'realpath';
use File::Basename;
use File::Copy;

my $program = realpath($0);
our $myBin = dirname($program);
our $myData = "$myBin/data";
our $myConf = "$myBin/conf";
require "$myBin/VIGOR3.pm";
$|++;

my ( $genomeFasta, $tblFile, $utrDB ) = @ARGV;
if ( ! -e $genomeFasta ) { die "cannot find fasta $genomeFasta\n" };
if ( ! -r $genomeFasta ) { die "cannot read fasta $genomeFasta\n" };
if ( ! -e $tblFile ) { die "cannot find TBL file $tblFile\n" };
if ( ! -r $tblFile ) { die "cannot read TBL file $tblFile\n" };
if ( ! -w $tblFile ) { die "cannot update TBL file $tblFile\n" };
if ( ! -e $utrDB ) { die "cannot find UTR database $utrDB\n" };
if ( ! -r $utrDB ) { die "cannot read UTR database $utrDB\n" };

initialize_defaults();
my $vigorspace = create_workspace();
set_parameter( "selectivity", 1.0 );
set_reference_seqs( loadFasta( $utrDB) );

my $saveId = 0;
my $saveFile = "$tblFile.save$saveId";
while ( -e $saveFile ) {
    $saveId++;
    $saveFile = "$tblFile.save$saveId";
}

copy($tblFile, $saveFile);
print "original TBL file saved to $saveFile\n";

my @tblInfo = parseTBL( $tblFile );

open( my $FASTA, "<$genomeFasta" ) || die "could not open fasta file $genomeFasta\n"; 
while ( my $genome = next_sequence( $FASTA ) ) {
    print "\n>$$genome{id}\n";
    open( GENOME, ">$vigorspace/genome") || die "could not open $vigorspace/genome for output\n";
    print GENOME ">$$genome{id}\n$$genome{sequence}\n";
    close GENOME;
    
    my @features = findFeatures( \@tblInfo, $$genome{id} );

    my $cmd = "$myBin/blastall -p blastn -d $utrDB -i genome -e 1e-30 -g F -F \"\" -b 1 -m 7 -o blast.xml > blast.log";
    runCmd($vigorspace, $cmd);

    my @hsps = sort { blastOrder( $a, $b ) } parse_blastxml( "$vigorspace/blast.xml" );
#    print "\n$$genome{id} raw data\n";
#    print_blasthits( 0, @hsps );    
#    print "\n";

    my @hits = sort{ blastOrder( $a, $b ) } best_gene_blasthits( $genome, \@hsps, 1.0 );
    for my $i ( 0..@hits-1 ) {
        my $hit = $hits[$i];
        if ( $$hit{subject_end} < $$hit{subject_length} && $$hit{subject_end} > $$hit{subject_length} - 6 ) {
            my $end =  $$hit{query_right};
            if ( $$hit{orientation} == -1 ) {
                $end =  $$hit{query_left};
            }
            my $nextend = $end + $$hit{orientation};
            while ( $$hit{subject_end} < $$hit{subject_length} && ! defined $$genome{is_gap}{$nextend} ) {
                $end = $nextend;
                $nextend += $$hit{orientation};
                $$hit{subject_end}++;
            }
            if ( $$hit{orientation} == -1 ) {
                $$hit{query_left} = $end;
            }
            else {
                $$hit{query_right} = $end;
            }
        }
        if ( $$hit{subject_begin} > 1 && $$hit{subject_begin} < 5 ) {
            my $begin =  $$hit{query_left};
            if ( $$hit{orientation} == -1 ) {
                $begin =  $$hit{query_right};
            }
            my $nextbegin = $begin - $$hit{orientation};
            while ( $$hit{subject_begin} > 1 && ! defined $$genome{is_gap}{$nextbegin} ) {
                $begin = $nextbegin;
                $nextbegin -= $$hit{orientation};
                $$hit{subject_begin}--;
            }
            if ( $$hit{orientation} == -1 ) {
                $$hit{query_right} = $begin;
            }
            else {
                $$hit{query_end} = $begin;
            }
        }
        my $qsize = $$hit{query_right} - $$hit{query_left} + 1;
        my $ssize = $$hit{subject_end} - $$hit{subject_begin} + 1;
        if ( $qsize > 1.1 * $ssize ) { $hits[$i] = undef }
    }
    @hits = sort { blastOrder( $a, $b ) } remove_undefs( @hits );
    
    print "\nbest hits\n";
    print_blasthits( 0, @hits );    
    print "\n";
#    print "\n$$genome{id} Blast Results\n";
#    print_blasthits( 0, @hits );
    
    my $gene;
    my $genesize;
    for my $hit ( @hits ) {
        my $hitgene = get_reference_name( $$hit{subject_id} );
        if ( ! defined $gene || $$gene{name} ne $hitgene ) {
            $gene = findGene( \@features, $hitgene );
            if ( ! defined $gene ) { next }
            $genesize = abs( $$gene{stop} - $$gene{start} ) + 1;
            if ( ! defined $$gene{utr_start} ) {
                $$gene{utr_start} = $$gene{start};
                $$gene{utr_stop} = $$gene{stop};
                $$gene{orig_trunc5} = $$gene{trunc5};
                $$gene{orig_trunc3} = $$gene{trunc3};
                $$gene{trunc5} = "<";
                $$gene{trunc3} = ">";
            }
        }
        
        if ( $$gene{start} > $$gene{stop} ) {
            if ( $$hit{orientation} != -1 ) { next }
            my $overlap = minval( $$hit{query_right}, $$gene{start} )
                - maxval( $$hit{query_left}, $$gene{stop} ) + 1;
            my $updated = 0;
#print "gene $$gene{name}  location $$gene{start}..$$gene{stop}  size $genesize  overlap $overlap\n";
            if ( $overlap >= 0.75 * $genesize ) {
                my $oldtrunc5 = $$gene{trunc5};
                my $oldtrunc3 = $$gene{trunc3};
                if ( $$gene{orig_trunc3} ne ">" && $$hit{query_left} < $$gene{utr_stop} ) {
                    $$gene{utr_stop} = $$hit{query_left};
                    if ( $$hit{subject_end} == $$hit{subject_length} ) { $$gene{trunc3} = "" }
                    $updated = 1;
                }
                if ( $$gene{orig_trunc5} ne "<" && $$hit{query_right} > $$gene{utr_start} ) {
                    $$gene{utr_start} = $$hit{query_right};
                    if ( $$hit{subject_begin} == 1 ) { $$gene{trunc5} = "" }
                    $updated = 1;
                }
                if ( $updated ) {
                    print "  updated gene $$gene{name} from $oldtrunc5$$gene{start}..$oldtrunc3$$gene{stop} to $$gene{trunc5}$$gene{utr_start}..$$gene{trunc3}$$gene{utr_stop}\n";
                }
            }
        }
        else {
            if ( $$hit{orientation} != 1 ) { next }
            my $overlap = minval( $$hit{query_right}, $$gene{stop} )
                - maxval( $$hit{query_left}, $$gene{start} ) + 1;
#print "gene $$gene{name}  location $$gene{start}..$$gene{stop}  size $genesize  overlap $overlap\n";
            if ( $overlap >= 0.75 * $genesize ) {
                my $updated = 0;
                my $oldtrunc5 = $$gene{trunc5};
                my $oldtrunc3 = $$gene{trunc3};
                if ( $$gene{orig_trunc5} ne "<" && $$hit{query_left} < $$gene{utr_start} ) {
                    $$gene{utr_start} = $$hit{query_left};
                    if ( $$hit{subject_begin} == 1 ) { $$gene{trunc5} = "" }
                    $updated = 1;
                }
                if ( $$gene{orig_trunc3} ne ">" && $$hit{query_right} > $$gene{utr_stop} ) {
                    $$gene{utr_stop} = $$hit{query_right};
                    if ( $$hit{subject_end} == $$hit{subject_length} ) { $$gene{trunc3} = "" }
                    $updated = 1;
                }                
                if ( $updated ) {
                    print "  updated gene $$gene{name} from $oldtrunc5$$gene{start}..$oldtrunc3$$gene{stop} to $$gene{trunc5}$$gene{utr_start}..$$gene{trunc3}$$gene{utr_stop}\n";
                }
            }
        }
    }    
}
close $FASTA;
if ( -e $vigorspace && -d $vigorspace ) { system "rm -rf $vigorspace" } 

rewriteTBL( \@tblInfo, $tblFile );
#print `cat $tblFile` . "\n";
exit(0);    

sub blastOrder {
    my ( $a, $b ) = @_;
    
    if ( get_reference_name( $$a{subject_id} ) lt get_reference_name( $$b{subject_id} ) ) { return -1 }
    if ( get_reference_name( $$a{subject_id} ) gt get_reference_name( $$b{subject_id} ) ) { return 1 }
    return $$b{pct_scoverage} <=> $$a{pct_scoverage};
}

sub parseTBL {
    my ( $tblFile ) = @_;
    
    my @tblinfo;
    my $cmd = "cat $tblFile | sed 's/\\r/\\n/g'";
    my $tblData = runCmdAndGetResults('.', $cmd);
    my @seqs = split />Features /, $tblData;
    shift @seqs;
    while ( @seqs ) {
        my $seqbody = shift @seqs;
        chomp $seqbody;
        my @seqdata = split /(([><]*[0-9]+)\t([><]*[0-9]+)\t(gene|mat_peptide))/, $seqbody;
        my $seqid = shift @seqdata;
        chomp $seqid;
        my $sequence;
        $$sequence{id} = $seqid;
        push @tblinfo, $sequence;
        
        while ( @seqdata ) {
            shift @seqdata;
            my $start = shift @seqdata;
            my $stop = shift @seqdata;
            my $type = shift @seqdata;
            my $body = substr( shift @seqdata, 1 );
            if ( substr( $body, length( $body ) - 1) ne "\n" ) { $body .= "\n" }

            my $feature;
            $$feature{type} = $type;
            $$feature{trunc5} = "";
            if ( $start =~ /^[><]/ ) {
                $$feature{trunc5} = substr( $start, 0, 1 );
                $start = substr( $start, 1 );
            }
            $$feature{start} = $start;
#print "START=\"$$feature{trunc5}\" + \"$$feature{start}\"\n";
            
            $$feature{trunc3} = "";            
            if ( $stop =~ /^[><]/ ) {
                $$feature{trunc3} = substr( $stop, 0, 1 );
                $stop = substr( $stop, 1 );
            }
            $$feature{stop} = $stop;

            $$feature{body} = $body;
            if ( $$feature{type} eq "gene" ) {
                my @genelines = split /\n/, $$feature{body};
                for my $geneline ( @genelines ) {
                    chomp $geneline;
                    if ( $geneline =~ /^\t\t\tgene\t([^ ]+)$/ ) {
                        $$feature{gene} = $1;
                    }
                    if ( $geneline =~ /^\t\t\tlocus_tag\t([^ ]+)$/ ) {
                        $$feature{locus_tag} = $1;
                    }
                }
                if ( defined $$feature{gene} ) {
                    $$feature{name} = $$feature{gene};
                }
                else {
                    $$feature{name} = $$feature{locus_tag};
                }
            }

            if ( defined $$sequence{features} ) {
                push @{$$sequence{features}}, $feature;
            }
            else {
                my @features = ( $feature );
                $$sequence{features} = \@features;
            } 
        }
    }
    
    return @tblinfo;
}

sub findFeatures {
    my ( $tblinfo, $genomeid ) = @_;
    for my $genome ( @$tblinfo ) {
        if ( $$genome{id} eq $genomeid ) {
            if ( defined $$genome{features} ) { return @{ $$genome{features} } }
        }
    }
    return undef;
}

sub findGene {
    my ( $features, $geneid ) = @_;
    
    for my $feature ( @$features ) {
        if ( index ( $$feature{body}, "\tpseudogene\t" ) < 0  ) {
            if ( defined $$feature{name} && $$feature{name} eq $geneid ) { return $feature }
        }
    }
    return undef;
}

sub rewriteTBL {
    my ( $tblinfo, $tblfile ) = @_;
    
    open( TBL, ">$tblfile" ) || die "cannot write to $tblfile\n";
    
    for my $ts ( @$tblinfo ) {
        print TBL ">Features $$ts{id}\n";
        for my $feature ( @{$$ts{features}} ) {
            if ( defined $$feature{utr_start} ) {
#print "trunc5=\"$$feature{trunc5}\"\n";
#print "utr_start=\"$$feature{utr_start}\"\n";
                print TBL "$$feature{trunc5}$$feature{utr_start}\t$$feature{trunc3}$$feature{utr_stop}\t$$feature{type}\n$$feature{body}";
            }
            else {
                print TBL "$$feature{trunc5}$$feature{start}\t$$feature{trunc3}$$feature{stop}\t$$feature{type}\n$$feature{body}";
            }
            if ( defined $$feature{utr_start} && $$feature{utr_start} != $$feature{start} ) {
                my $utr5_end = $$feature{start};
                if ( $$feature{start} < $$feature{stop} ) {
                    $utr5_end--;
                }
                else {
                    $utr5_end++;
                }
                print TBL "$$feature{trunc5}$$feature{utr_start}\t$utr5_end\t5'UTR\n";
                if ( defined $$feature{locus_tag} ) {
                    print TBL "\t\t\tlocus_tag\t$$feature{locus_tag}\n";
                }
                if ( defined $$feature{gene} ) {
                    print TBL "\t\t\tgene\t$$feature{gene}\n";
                }
            }
            if ( defined $$feature{utr_stop} && $$feature{utr_stop} != $$feature{stop} ) {
                my $utr3_start = $$feature{stop};
                if ( $$feature{start} < $$feature{stop} ) {
                    $utr3_start++;
                }
                else {
                    $utr3_start--;
                }
                print TBL "$utr3_start\t$$feature{trunc3}$$feature{utr_stop}\t3'UTR\n";
                if ( defined $$feature{locus_tag} ) {
                    print TBL "\t\t\tlocus_tag\t$$feature{locus_tag}\n";
                }
                if ( defined $$feature{gene} ) {
                    print TBL "\t\t\tgene\t$$feature{gene}\n";
                }
            }
        }
    }
    close TBL;
}

