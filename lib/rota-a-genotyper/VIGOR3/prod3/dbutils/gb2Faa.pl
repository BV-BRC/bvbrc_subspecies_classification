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

# for each genbank file
while ( @ARGV ) {

    my $gbFile = shift @ARGV;
    print STDERR "$gbFile\n";

    # get genome information
    my $genbank = Bio::SeqIO->new( -file => "<$gbFile", -format => "Genbank" );
    my $genome = $genbank->next_seq;
    my $genome_id = basename( $gbFile );
    $genome_id =~ s/\..*$//;

    # get CDS features from genbank file
    my @proteins;
    for my $feature ( $genome->get_SeqFeatures ) {
        if ( $feature->primary_tag eq "CDS" ) {
            push @proteins, $feature;
        }
    }

    # output CDS
    @proteins = sort { $a->location->start <=> $b->location->start } @proteins;
    my $protein_number = 0;
    for my $protein ( @proteins ) {
        my $id;
        if ( $protein->has_tag( "protein_id" ) ) {
            for my $value ( $protein->get_tag_values( "protein_id" ) ) {
                $id = $value;
                last;
            }
        }
        else {
            $protein_number++;
            $id = "CDS#$protein_number";
        }
        
        my $gene;
        if ( $protein->has_tag( "gene" ) ) {
            for my $value ( $protein->get_tag_values( "gene" ) ) {
                $gene = $value;
                last;
            }
        }
        else {
            $gene="unknown";
        }

        my $product = "uncharacterized protein";
        if ( $protein->has_tag( "product" ) ) {
            for my $value ( $protein->get_tag_values( "product" ) ) {
                $product = $value;
                last;
            }
        }
        
        my $location = $protein->location->to_FTstring();        

        my $aa = "";
        for my $value ( $protein->get_tag_values( "translation" ) ) {
            $aa .= $value;
        }
        $aa =~ s/(.{60})/$1\n/g;
            
        my $geneid = $gene;
        $geneid =~ s/\W+/_/g;
        my $prodid = $product;
        $prodid =~ s/\W+/_/g;
        
        open( OUT, ">$genome_id.$geneid.$prodid.faa" );
        print OUT ">$id  gene=\"$gene\"  product=\"$product\"  location=\"$location\"\n$aa\n";
        close OUT;
    }
}

exit(0);
