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

# for each genbank file

# get features from genbank file (gene, protein, mat_peptide)
while ( @ARGV ) {
    my $gbFile = shift @ARGV;
        
    my $genbank = Bio::SeqIO->new( -file => "<$gbFile", -format => "Genbank" );
    my $genome = $genbank->next_seq;
    my $accession = $genome->display_id;
    my $no = 0;
    my @repeats;
    for my $feature ( $genome->get_SeqFeatures ) {
        #print "type = " . $feature->primary_tag . "\n";
        if ( $feature->primary_tag !~ /misc_RNA/i ) { next }
        
        $no++;
        my $start = $feature->location->start;
        my $end =  $feature->location->end;
        my $defline = ">$accession" . "_rna$no location=$start..$end";
        for my $tag ($feature->get_all_tags ) {
            $defline .= " $tag=\"";
            for my $value ($feature->get_tag_values($tag)) {                
                $defline .= $value . " ";
            }
            $defline =~ s/ +$//;
            $defline .= "\"";
        }
        my $sequence = $feature->spliced_seq->seq;
        print "$defline\n$sequence\n";
    }
}
exit(0);
