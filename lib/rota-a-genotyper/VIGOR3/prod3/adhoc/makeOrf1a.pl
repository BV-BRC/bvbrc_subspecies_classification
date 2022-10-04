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
require "$myBin/VIGOR.pm";

initialize_defaults();
my $vigorspace = create_workspace();

# for each genbank file
while ( @ARGV ) {

# get features from genbank file (gene, protein, mat_peptide)
    my $gbFile = shift @ARGV;
    
    my $genbank = Bio::SeqIO->new( -file => "<$gbFile", -format => "Genbank" );
    my $virus = $genbank->next_seq;
    my $genome;
    $$genome{id} = $virus->display_id;
    $$genome{sequence} = $virus->seq();
    $$genome{seqlen} = length( $$genome{sequence} );
    
    my @classification = $virus->species->classification;
    $$genome{organism} = shift @classification;
         
    my @genes;
    my @proteins;
    my @matpeps;
    my %protgene;
    for my $feature ( $virus->get_SeqFeatures ) {
        if ( $feature->location->start_pos_type ne 'EXACT'
                || $feature->location->end_pos_type ne 'EXACT' 
                || $feature->location->location_type ne 'EXACT' ) {
            next;
        }
        if ( $feature->primary_tag eq "CDS" ) {
            push @proteins, $feature;
        }
    }

    
    for my $protein ( @proteins ) {
        
        if ( $protein->location =~ /[><]/ ) { next }
        if ( ! $protein->has_tag( "protein_id" ) ) { next }

        
        if ( ! $protein->has_tag( "gene" ) ) { next }
        my @tmp = $protein->get_tag_values( "gene" );
        my $gene = shift @tmp;
        if ( $gene !~ /1ab/i ) { next }

        @tmp = $protein->get_tag_values( "protein_id" );
        my $id = shift @tmp;
        $id =~ s/\..*$//;
        
        my $start = $protein->location->start;
        my $end = $protein->location->end;
        
        for my $slide ( find_regexp( "AAATTTC", $$genome{sequence}, $start+9000 ) ) {
            #print "$$slide{begin}-$$slide{end} = $$slide{string}\n";
            my $gp1a = subsequence( $$genome{sequence}, $start, $$slide{end} + 300 );
            my $gp1a_aa = DNA2AA( $gp1a, 1 );
            my $stop = index( $gp1a_aa, "*" );
            $gp1a_aa = substr( $gp1a_aa, 0, $stop );
            my $gp1alen = length( $gp1a_aa );
            my $gp1aend = $start + 3 * $gp1alen + 2;
            $gp1a_aa =~ s/(.{60})/$1\n/g;
                    
            print ">$id" . "a position=$start..$end length=$gp1alen gene=\"orf1a\" product=\"polyprotein 1a\" shared_cds=\"orf1ab\" matpepdb=\"default\" [$$genome{organism}]\n$gp1a_aa\n";
        }        
    }
}
exit(0);
