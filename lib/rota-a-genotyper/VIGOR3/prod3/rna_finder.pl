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

my $program = realpath($0);
our $myBin = dirname($program);
our $myData = "$myBin/data";
our $myConf = "$myBin/conf";

require "$myBin/VIGOR3.pm";
*STDERR = *STDOUT;
$|++;

initialize_defaults();
my $vigorspace = create_workspace();

my $use_locus_tags = get_parameter( "use_locus_tags" ); 
if ( $ARGV[0] eq "-l" ) {
    $use_locus_tags = 0;
    shift @ARGV;
}
elsif ( $ARGV[0] eq "-L" ) {
    $use_locus_tags = 1;
    shift @ARGV;
}
my $fastaFile = $ARGV[0];
my $tblFile = $ARGV[1];
my $rnadb = $ARGV[2];

open ( my $FASTA, "<$fastaFile" );
open ( my $TBL, ">$tblFile" );
my %refs = loadFasta( $rnadb );
set_reference_seqs( %refs );
while ( my $genome = next_sequence( $FASTA ) ) {
    annotate_rna( $genome, $TBL );
}
close $TBL;
exit(0);

sub    annotate_rna {
    my ( $genome, $TBL ) = @_;
    
    my $seq = $$genome{sequence};
    $seq =~ s/(.{60})_/$1\n/g;
    open( FASTA, ">$vigorspace/genome.fasta" );
    print FASTA ">$$genome{id}\n$seq\n";
    close FASTA;
    
    my $blastxml = "$vigorspace/blastxml";
    my $cmd = "$myBin/blastall -p blastn -i $vigorspace/genome.fasta -F \"\" -b 250 -e 1e-70 -m 7 -d $rnadb > $blastxml";
    runCmd($vigorspace, $cmd);
    
    my @hits = sort { $$b{vigor_matchwgt} <=> $$a{vigor_matchwgt}} parse_blastxml( $blastxml, 0 );
    print_blasthits( 3, @hits );

    my %rna;
    for my $hit ( @hits ) {
        my $subjid = $$hit{subject_id};
        my $qstart = $$hit{query_left};
        my $qend = $$hit{query_right};
        my $ori = $$hit{orientation};
        if ( $ori == -1 ) {
            ( $qstart, $qend ) = ( $qend, $qstart );
        }

        my $sstart = $$hit{subject_begin};
        my $send = $$hit{subject_end};
        
        my $t5 = $sstart - 1;
        my $t3 = $$hit{subject_length} - $send;
        
        my $rna_id;
        if ( $$hit{subject_definition} =~ / rna_id="([^"]+)"/ ) {
            $rna_id = $1;            
        }
        else {
            die "\nreference RNA sequence \"$subjid\" does not define an rna_id\n";
        }
        my $gene;
        if ( $$hit{subject_definition} =~ / gene="([^"]+)"/ ) {
            $gene = $1;            
        }
        my $product;
        if ( $$hit{subject_definition} =~ / product="([^"]+)"/ ) {
            $product = $1;            
        }
        my $note;
        if ( $$hit{subject_definition} =~ / note="([^"]+)"/ ) {
            $note = $1;            
        }
        my $rna_type = "misc_RNA";
        if ( $$hit{subject_definition} =~ / rna_type="([^"]+)"/ ) {
            $rna_type = $1;            
        }
        elsif ( "$$hit{subject_definition} " =~ / rna_type=([^ ]+) / ) {
            $rna_type = $1;            
        }
        if ( !defined $rna{$rna_id} ) {
            my %miscrna;
            $miscrna{qstart} = $qstart;
            $miscrna{qend} = $qend;
            $miscrna{sstart} = $sstart;
            $miscrna{send} = $send;
            $miscrna{t5} = $t5;
            $miscrna{t3} = $t3;
            $miscrna{gene} = $gene;
            $miscrna{product} = $product;
            $miscrna{rna_type} = $rna_type;
            $miscrna{note} = $note;
            $miscrna{orientation} = $ori;
            $rna{$rna_id} = \%miscrna;
        }
        else {
            if ( $qstart > $rna{$rna_id}{qstart} - 50
                    && $qstart <  $rna{$rna_id}{qstart}
                    && $t5 <= $rna{$rna_id}{t5} ) { 
                $rna{$rna_id}{qstart} = $qstart;
                $rna{$rna_id}{t5} = $t5;
            }
            if ( $qend < $rna{$rna_id}{qend} + 50
                    && $qend > $rna{$rna_id}{qend}
                    && $t3 <= $rna{$rna_id}{t3} ) {
                $rna{$rna_id}{qend} = $qend;
                $rna{$rna_id}{t3} = $t3;
            }
        }
    }
    unlink "$vigorspace/genome.fasta";
    if ( keys %rna ) {
        print $TBL ">Features $$genome{id}\n";
    }
    for my $rna_id ( sort { $rna{$a}{qstart} <=> $rna{$b}{qstart} } keys %rna ) {
        if ( $rna{$rna_id}{t5} > 0 ) {
            $rna{$rna_id}{qstart} = "<$rna{$rna_id}{qstart}";
        }
        if ( $rna{$rna_id}{t3} > 0 ) {
            $rna{$rna_id}{qend} = ">$rna{$rna_id}{qend}";
        }
        print $TBL "$rna{$rna_id}{qstart}\t$rna{$rna_id}{qend}\t$rna{$rna_id}{rna_type}\n";
        if ( $use_locus_tags ) {
            my $tag = $rna{$rna_id}{gene};
            if ( ! defined $tag ) { $tag = $rna_id }
            $tag =~ s/-/p/gi;
            $tag =~ s/_/x/g;
            print $TBL "\t\t\tlocus_tag\tvigor_$tag\n";
        }
        for my $attr ( "gene", "product", "note" ) {
            if ( defined $rna{$rna_id}{$attr} ) {
                print $TBL "\t\t\t$attr\t$rna{$rna_id}{$attr}\n";
            }
        }
    }
}
