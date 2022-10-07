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
our $myBin = dirname( dirname($program) );
our $myData = "$myBin/data3";
our $myConf = "$myBin/conf";

require "$myBin/VIGOR3.pm";
require "$myBin/db.pm";
use Bio::SeqIO;
$|++;

set_default_codons();
my $vigorspace = create_workspace();
my $outfile = shift @ARGV;
if ( -e $outfile ) { die "\nRequested output file \"$outfile\" already exists.\n" }
my $cmd = "touch $outfile";
runCmd('.', $outfile);

# for each genbank file
my $tmpfile = "$vigorspace/matpep.tmp";
open( OUT, ">$outfile" ) || die "\n could not open output file \"$outfile\"\n";
while ( @ARGV ) {

# get features from genbank file (gene, protein, mat_peptide)
    my $gbFile = shift @ARGV;
    my $genbank = Bio::SeqIO->new( -file => "<$gbFile", -format => "Genbank" );
    my $protein = $genbank->next_seq;
    my $accession = $protein->display_id;
    my $taxon = $protein->species->species;
    my @matpeps;
    my %protgene;
    for my $feature ( $protein->get_SeqFeatures ) {
        my $tag = $feature->primary_tag;
        if ( $feature->primary_tag eq "mat_peptide" ) {
            push @matpeps, $feature;
        }
    }
    my $pepno = 0;
    for my $mp ( @matpeps ) {
        $pepno++;
        if ( $mp->location->start_pos_type ne 'EXACT'
                || $mp->location->end_pos_type ne 'EXACT' 
                || $mp->location->location_type ne 'EXACT' ) {
            print "fuzzy location: accession $accession matpep #$pepno\n";
            next;
        }
        my $protein_id = $accession . "_mp$pepno";
        if ( $mp->has_tag( "protein_id" ) ) {
            for my $value ( $mp->get_tag_values( "protein_id" ) ) {
                $protein_id = $value;
                last;
            }
        }
        my $product_id = "mature peptide";
        if ( $mp->has_tag( "product" ) ) {
            for my $value ( $mp->get_tag_values( "product" ) ) {
                $product_id = $value;
                last;
            }
        }
        my $sequence = DNA2AA( $mp->spliced_seq->seq, 1 );

        print OUT ">$protein_id organism=\"$taxon\" product=\"$product_id\"\n$sequence\n";
        print ">$protein_id organism=\"$taxon\" product=\"$product_id\"\n";
    }
}
close OUT;

exit(0);

