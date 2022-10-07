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
require "$myBin/VIGOR3.pm";

my $namefile = shift @ARGV;
my $dbfile = shift @ARGV;

if ( ! -r $dbfile ) { die "\ncannot read $dbfile\n" }
my $dbname = basename( $dbfile );

my %data;
if ( $namefile ne "" ) {
    if ( ! -r $namefile ) { die "\ncannot read $namefile\n" }
    my $cmd = "cat $namefile | sed 's/\r//g'";
    my $namedata = runCmdAndGetResults('.', $cmd);
    my @namedata = split /\n/, $namedata;
    for my $datum ( @namedata  ) {
        my @values = split /\t/, $datum;
        if ( @values ) {
            for my $i ( 0..@values-1 ) {
                if ( defined $values[$i] && $values[$i] =~ /^"(.*)"$/ ) {
                    $values[$i] = $1;
                }
                if ( defined $values[$i] ) {
                    $values[$i] =~ s/^ +//;
                    $values[$i] =~ s/ +$//;
                }
            }
            if ( $values[0] gt "" ) {
                if ( ! defined $values[1] || $values[1] le "" ) {
                    die "\n no product name define on line $datum\n";
                }
                my $gene = $values[0];
                my %info;
                $info{product} = $values[1];
                if ( defined $values[2] && $values[2] gt "" ) {
                    $values[2] =~ s/;/,/g;
                    if ( $values[2] =~ /,/ ) {
                        $values[2] =~ s/, */, /g;
                        $info{note} = "synonyms: $values[2]";
                    }
                    else {
                        $info{note} = "synonym: $values[2]";
                    }
                }
                if ( defined $values[3] && $values[3] gt "" ) {
                    if ( exists $info{note} ) {
                        $info{note} .= "; $values[3]";
                    }
                    else {
                        $info{note} = $values[3];
                    }                
                }
                $data{$gene} = \%info;
    #print_hash( $gene, \%info );
            }
        }
    }
}

my %fasta = loadFasta( $dbfile );
set_reference_seqs( %fasta );
for my $refid ( sort { get_reference_gene( $a ) cmp get_reference_gene( $b ) } keys %fasta ) {
    my $gene = get_reference_gene( $refid );
    my %ref = %{ $fasta{$refid} };
    if ( $ref{sequence} =~ /\*$/ ) { $ref{sequence} =~ s/\*+$// }
    $ref{defline} .= " ";
    $ref{defline} =~ s/ curated=.//g;
    $ref{defline} =~ s/ original[^=]*="[^"]*"//g;
    $ref{defline} =~ s/ original[^=]*=[^ ]* / /g;
    $ref{defline} =~ s/ besthit[^=]*="[^"]*"//g;
    $ref{defline} =~ s/ besthit[^=]*=[^ ]* / /g;
    $ref{defline} =~ s/ length=[^ ]* / /g;
    $ref{defline} =~ s/ db="[^"]*"//g;
    $ref{defline} =~ s/ db=[^ ] / /g;
    $ref{defline} =~ s/ +$//;
    if ( exists $data{$gene} ) {
        $ref{defline} =~ s/ product=\"[^"]*\"//;
        $ref{defline} .= " product=\"$data{$gene}{product}\"";
        $ref{defline} =~ s/ note=\"[^"]*\"//;
        if ( exists $data{$gene}{note} ) {
            $ref{defline} .= " note=\"$data{$gene}{note}\"";
        }
    }
    $ref{defline} .= " length=" . length( $ref{sequence} ) . " db=\"$dbname\"";
    $ref{defline} = clean_defline( $ref{defline} );
    $ref{sequence} =~ s/(.{60})/$1\n/g;
    if ( substr( $ref{sequence}, length( $ref{sequence} ) -1 ) ne "\n" ) { $ref{sequence} .= "\n" }
    print ">$ref{defline}\n$ref{sequence}\n";
}

exit( 0 );

sub clean_defline {
    my ( $defline ) = @_;
    
    my %tags;
    my @misc;
    
    my $tmp = $defline;
    
    my $id = $tmp;
    $id =~ s/ .*//;
    $tmp =~ s/^[^ ]* *//;
    
    while ( $tmp ) {
        if ( $tmp =~ /^([^= ]+)=("[^"]*")/ ) {
            my ( $tag, $value ) = ( $1, $2 );
            if ( ! exists $tags{$tag} ) { $tags{$tag} = $value }
            $tmp =~ s/^([^= ]+)=("[^"]*") *//;
        }
        elsif ( "$tmp " =~ /^([^= ]+)=([^ ]*) / ) {
            my ( $tag, $value ) = ( $1, $2 );
            if ( ! exists $tags{$tag} ) { $tags{$tag} = $value }
            $tmp =~ s/^([^= ]+)=([^ ]*) *//;
        }
        elsif ( "$tmp " =~ /^(is_[^ ]*) / ) {
            my $tag = $1;
            $tags{$tag} = undef;
            $tmp =~ s/^(is_[^ ]*) *//;
        }
        elsif ( "$tmp " =~ /^([^ ]*) / ) {
            my $value = $1;
            push @misc, $value;
            $tmp =~ s/^([^ ]*) *//;
        }
    }
    
    $tmp = $id;
    for my $tag ( sort { $a cmp $b } keys %tags ) {
        if ( defined $tags{$tag} ) {
            $tmp .= " $tag=$tags{$tag}";
        }
        else {
            $tmp .= " $tag";
        }
    }
    if ( @misc ) {
        $tmp .= " " . join( " ", @misc );
    }
    
    return $tmp;
}