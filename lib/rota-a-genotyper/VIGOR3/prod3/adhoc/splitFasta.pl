#!/usr/local/bin/perl

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
use Cwd 'realpath';
use File::Basename;
our $errorMessage;

require "getopts.pl";

#-----------------------------------------------------------------------
# get options from command line
my ( $fasta, $maxseqcnt, $maxseqlen, $outdir ) = &initialize;
print "SPLITFASTA: $fasta, $maxseqcnt, $maxseqlen, $outdir \n";
my $gzip = 0;
if ( substr( $fasta, length( $fasta ) - 3 ) eq ".gz" ) { $gzip = 1 }
my $partition_root = basename( $fasta );
if ( $gzip ) { $partition_root =~ s/\.gz$// }
$partition_root .= ".p";
$partition_root = $outdir . "/" . $partition_root;

if ( $gzip ) {
    open( IN, "gunzip -c $fasta |" );
} else {
    open( IN, "<$fasta" );
}
my $partition_maxcnt = $maxseqcnt;
my $partition_maxlen = $maxseqlen;
my $partition_id = 1;
my $partition_cnt = 0;
my $partition_len = 0;
my $partition_file = $partition_root . $partition_id . ".fasta";
open( OUT, ">$partition_file" );

my $prevlen = 0;
while ( my $line = <IN> ) {
    if ( substr( $line, 0, 1 ) eq ">" ) {
        if ( ( defined $partition_maxcnt && $partition_cnt >= $partition_maxcnt )
                || ( defined $partition_maxlen && $partition_len+$prevlen > $partition_maxlen ) ) {
            close( OUT );
            print "$partition_file: $partition_cnt $partition_len\n";
            $partition_id++;
            $partition_file = $partition_root . $partition_id . ".fasta";
            open( OUT, ">$partition_file" );
            $partition_cnt = 1;
            $partition_len = 0;
        } else {
            $partition_cnt++;
        }
        $prevlen = 0;
    } else {
        $partition_len += length( $line );
        $prevlen += length( $line );        
    }
    print OUT $line;             
}
close( OUT );
print "$partition_file: $partition_cnt $partition_len\n";
exit(0);

sub initialize {
    use vars qw( $opt_F $opt_N $opt_L $opt_h $opt_d );
    &Getopts('F:N:L:d:h');

    if ( $opt_h ) {
        print
"
usage: ./split.pl -F fasta -N number of seqs

-F <fasta> path to fasta file to be split
-N <number of seqs> maximum number of seqs per split
-L <total length of seq>   maximum total length of seqs per split
-d <output directory> optional, defaults to current directory
";
        exit(0);
    }

    if ( !$opt_F ) {
        die "\nYou must specify the protein fasta file (-F)."
    } elsif ( ! -e $opt_F) {
        die "\nThe specified fasta file does not exist: \"-F " . $opt_F . "\".";
    }
    
    if ( $opt_N ) {
        if ( $opt_N < 1 ) {
            die "\n-N must greater than 0\n";
        }
        
    }
    
    if ( $opt_L ) {
        if ( $opt_L < 100000 ) {
            die "\n-L must be >= 100000\n";
        }
        
    }
    
    if ( !$opt_N && ! $opt_L ) {
        die "\nEither -N or -L or both must be specified\n" 
    }
    
    my $outdir = ".";
    if ( defined $opt_d ) { $outdir = $opt_d }
    $outdir = realpath( $outdir );
    if ( ! -d $outdir ) {
        die "\n\"$outdir\" is not a directory.\n";
    } elsif ( ! -e $outdir ) {
        die "\nCannot find output directory \"$outdir\".\n";
    } elsif ( ! -w $outdir ) {
        die "\nOutput directory \"$outdir\" is not write-enabled.\n";
    }
    if ( substr( $outdir, length($outdir)-1, 1) eq "/" ) {
        $outdir = substr( $outdir, 0, length($outdir)-1 );
    }
    
    return ( $opt_F, $opt_N, $opt_L, $outdir );
}
