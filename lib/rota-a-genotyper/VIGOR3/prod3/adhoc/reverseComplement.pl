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
our $myData = "$myBin/data";
our $myConf = "$myBin/conf";
require "$myBin/VIGOR.pm";
$|++;

set_default_codons();

my %fasta = loadFasta( $ARGV[0] );
for my $id ( sort keys %fasta ) {
    my $seq = $fasta{$id};
    $$seq{id} = "RV-$$seq{id}";
    $$seq{defline} = "RV-$$seq{defline}";
    $$seq{sequence} = reverse_complement( $$seq{sequence} );
    $$seq{sequence} =~ s/(.{60})/$1\n/g;
    if ( substr( $$seq{sequence}, length( $$seq{sequence} ) - 1, 1 ) ne "\n" ) {
        $$seq{sequence} .= "\n";
    }
    print ">$$seq{defline}\n$$seq{sequence}\n";
}
exit(0);