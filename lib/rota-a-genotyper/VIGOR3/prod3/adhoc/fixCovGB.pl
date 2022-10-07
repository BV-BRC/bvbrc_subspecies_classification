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
use Bio::SeqIO;
use Getopt::Std;
use Cwd 'realpath';
use File::Basename;

my $program = realpath($0);
our $myBin = dirname( dirname($program) );
our $myData = "$myBin/data";
our $myConf = "$myBin/conf";
require "$myBin/VIGOR3.pm";

set_default_splice_pairs();
my %splicePairs = %{get_splice_pairs()};

my $mature_peptides = 0;
if ( $ARGV[0] eq "-m" ) {
    $mature_peptides = 1;
    shift @ARGV;
}

# for each genbank file
while ( @ARGV ) {

# get features from genbank file (gene, protein, mat_peptide)
    my $gbFile = shift @ARGV;
    print STDERR "$gbFile\n";
    
    my $genbank = Bio::SeqIO->new( -file => "<$gbFile", -format => "Genbank" );
    my $virus = $genbank->next_seq;
    my $genome;
    $$genome{id} = $virus->display_id;
    $$genome{sequence} = $virus->seq();
    $$genome{seqlen} = length( $$genome{sequence} );

    my $gppfx = $$genome{id};
    $gppfx =~ s/\..*$//;
    my $gpno = 0;
         
    for my $feature ( $virus->get_SeqFeatures ) {
        if ( $feature->primary_tag ne "CDS" ) { next }
        my $exons = feature_exons( $genome, $feature, 2 );
        if ( @$exons != 2 ) { next }
        my $orfa = shift @$exons;
        my $slippage = subsequence( $$genome{sequence}, $$orfa{end}-6, $$orfa{end} );
        print "$$genome{id}: $$orfa{end}: $slippage\n";
        if ( $slippage =~ /tttaaac/i ) { next }
        
        my $gseq = $$genome{sequence};
        $gseq =~ s/(.{60})/$1\n/g;
        $gseq .= "\n";
        open FASTA, ">$gppfx.fasta";
        print FASTA ">$gppfx\n$gseq\n";
        close FASTA;
        
        my $outPrefix = $gbFile;
        $outPrefix = dirname( $gbFile ) . "/$gppfx";
#        system "/local/devel/ANNOTATION/VIRAL/VIGOR/new/VIGOR3.pl -G $gbFile -i $gppfx.fasta -O $outPrefix";
        my $cmd = "/local/devel/VIRIFX/software/VIGOR3/prod3/VIGOR3.pl -D covall_db -i $gppfx.fasta -O $outPrefix";
        runCmd('.', $cmd);
        last;
    }
}
exit(0);

sub feature_exons {
    my ( $genome, $feature, $flag ) = @_;
    if ( ! defined $flag ) { $flag = -1 }

#   for my $tag ($feature->get_all_tags) {             
#      print "  tag: ", $tag, "\n";             
#      for my $value ($feature->get_tag_values($tag)) {                
#         print "    value: ", $value, "\n";             
#      }          
#   }       

#    my @tmp = $feature->get_tag_values( "protein_id" );
#    my $id = shift @tmp;
#print "FEATURE $id\n";
    
    my $frameadj = 0;
    if ( $feature->has_tag( "codon_start" ) ) {
        my @tmp = $feature->get_tag_values( "codon_start" );
        $frameadj = shift @tmp;
        $frameadj--;
    }

    my @locations;
    if ( $feature->location->isa('Bio::Location::SplitLocationI') ) {
        for my $location ( $feature->location->sub_Location ) {
            push @locations, $location;
        }
    }
    else {
        push @locations, $feature->location;
    }
    @locations = sort { $a->strand * $a->start <=> $b->strand * $b->start } @locations;
    
    my @exons;
    for my $location ( @locations ) {

        my $start = $location->start;
        my $end = $location->end;
        my $strand = $location->strand;
        if ( $strand == -1 ) {
            ( $start, $end ) = ( $end, $start );
        }
        
#print "frameadj=$frameadj  location=$start-$end ($strand)\n" ;
        my $exon;
        $$exon{start} = $start;
        $$exon{codon_start} = $start + $frameadj;
        $$exon{end} = $end;
        $$exon{codon_end} = $end;
        $$exon{strand} = $strand;
        if ( $strand == -1 ) {
            ( $$exon{codon_left}, $$exon{codon_right} ) = ( $$exon{codon_end}, $$exon{codon_start} );
            $$exon{frame} = ( $$genome{seqlen} + 1 - $$exon{codon_start} ) % 3;
            $$exon{frame} = -$$exon{frame};
            if ( $$exon{frame} == 0 ) { $$exon{frame} = -3 }
        }
        else {
            ( $$exon{codon_left}, $$exon{codon_right} ) = ( $$exon{codon_start}, $$exon{codon_end} );
            $$exon{frame} = $$exon{codon_start} % 3;
            if ( $$exon{frame} == 0 ) { $$exon{frame} = 3 }
        }

        $frameadj = -( $$exon{codon_right} - $$exon{codon_left} + 1 ) % 3;
        if ( $frameadj < 0 ) { $frameadj += 3 }

#print "exon=$$exon{codon_left}-$$exon{codon_right} ($$exon{frame})\n";
        while( ( $$exon{codon_right} - $$exon{codon_left} + 1 ) % 3 > 0 ) {
            $$exon{codon_end} -= $strand;
            if ( $strand == 1 ) {
                $$exon{codon_right} -= 1;
            }
            else {
                $$exon{codon_left} += 1;
            }
#print "exon=$$exon{codon_left}-$$exon{codon_right} ($$exon{frame})\n";
        }
        push @exons, $exon;
    }

    return \@exons;    
}
