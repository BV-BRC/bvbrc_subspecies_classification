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
require "$myBin/VIGOR3.pm";
$|++;
initialize_defaults();
my $vigorspace = create_workspace();
set_parameter( "verbose", 1 );

my %mpref = loadFasta( $ARGV[1] );
my %refseqs = loadFasta( $ARGV[0] );
set_reference_seqs( %mpref );

# write results, sorted by gene
for my $seqid ( sort keys %refseqs ) {
    my $gene;
    $$gene{start_truncation} = 0;
    $$gene{stop_truncation} = 0;
    $$gene{gene_id} = $seqid;
    $$gene{ref_id} = $seqid;
    $$gene{protein} = $refseqs{$seqid}{sequence};
    $$gene{protein_length} = $refseqs{$seqid}{seqlen};
    my $matpeps = mapto_polyprotein( $gene, undef , $ARGV[1] );
    
    for my $mp ( @$matpeps ) {
        my $pepseq = subsequence( $$gene{protein}, $$mp{pep_start}, $$mp{pep_end} );
        print ">$$mp{pep_id} $$mp{pep_start}-$$mp{pep_end}\n$pepseq\n";
    }
}
exit(0);