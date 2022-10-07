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

my %refseqs = loadFasta( $ARGV[0] );
my $tag = $ARGV[1];
my $curated;
if ( defined $ARGV[2] ) { $curated = $ARGV[2] }
set_reference_seqs( %refseqs );

# write results, sorted by gene
for my $seqid ( sort { get_defline_tag( $refseqs{$a}{defline}, $tag ) cmp get_defline_tag( $refseqs{$b}{defline}, $tag ) } keys %refseqs ) {
    if (  defined $curated && $refseqs{$seqid}{defline} !~ / curated=$curated/i ) { next }
#print STDERR "$tag=" . get_defline_tag( $refseqs{$seqid}{defline}, $tag ) . "\n";
    my $ref = get_reference_seq( $seqid );
    my $seq = $$ref{sequence};
    $seq =~ s/\*$//;
    my $len = length( $seq );
    $seq =~ s/(.{60})/$1\n/g;
    $seq =~ s/\n$//;
    if ( $$ref{defline} =~ / length=/ ) {
        print "\n>$$ref{defline}\n$seq\n";
    }
    else {
        print "\n>$$ref{defline} length=$len\n$seq\n";
    }
}
exit(0);