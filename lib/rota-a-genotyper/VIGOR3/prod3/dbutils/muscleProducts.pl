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

if ( ! -e $ARGV[0] || ! -r $ARGV[0] ) { die "\ncannot read input from \"$ARGV[0]\"\n" }
my %refseqs = loadFasta( $ARGV[0] );
my $outprefix = $ARGV[1];
if ( ! defined $outprefix ) {
    $outprefix = "./" . basename( $ARGV[0] );
}
else {
    $outprefix =~ s/\\/\//g;
}
set_reference_seqs( %refseqs );

# write results, sorted by gene
my %products;
for my $seqid ( keys %refseqs ) {
    my $ref = get_reference_seq( $seqid );
    my $product = get_reference_product( $seqid );
    if ( ! defined $products{$product}{members} ) {
        my @tmp = ( $ref );
        my $product = get_reference_product( $seqid );
        $products{$product}{members} = \@tmp;
    }
    else {
        push @{$products{$product}{members}}, $ref;
    }
}

my $productid = 0;
for my $product ( sort keys %products ) {
    my @members = @{$products{$product}{members}};
    my $prod = $product;
    $prod =~ s/\W/_/g;
    $productid++;
    my $fastafile = $outprefix . "." . $prod . ".fasta";
    my $musclefile = $outprefix . "." . $prod . ".aln";
    open( OUT, ">$fastafile" ) || die "\ncannot write to \"$fastafile\"\n" ;
    for my $member ( @members ) {
        print OUT ">$$member{defline}\n$$member{sequence}\n";
    }
    close OUT;
    print "\n============================================================================================\n";
    if ( @members > 1 ) {
        my $cmd = "$myBin/muscle -clw -in $fastafile -out $musclefile";
        print "product # $productid  product $product\n$cmd\n";
        $cmd .= "&> /dev/null"; 
        runCmd('.', $cmd);
        $cmd = "cat $musclefile | sed 's/>/\\n>/'";
        my $results = runCmdAndGetResults('.', $cmd);
        print "$results\n";
    }
}
exit(0);