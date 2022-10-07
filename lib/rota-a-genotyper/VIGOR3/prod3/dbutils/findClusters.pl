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
$|++;

use Cwd 'realpath';
use File::Basename;
use Getopt::Std;
my $program = realpath($0);
our $myBin = dirname( dirname($program) );
require "$myBin/VIGOR3.pm";
$|++;

# get user options (input and output file paths, coverage and similarity)
my %args;
&getopts( 'c:i:o:s:u:', \%args );

my $infasta = $args{i};
if ( ! defined $infasta ) {
    die "\nNo input file specified (-i)\n";
}
elsif ( ! -e $infasta || ! -r $infasta ) {
    die "\nCannot read sequence from $infasta\n";
}
my %inseqs = loadFasta( $infasta );
    
my $outfasta = $args{o};
if ( ! defined $outfasta ) {
    die "\nNo output file specified (-o)\n";
}
elsif ( ! open OUT, ">$outfasta" ) {
    die "\nCannot write to \"$outfasta\"\n";
}

my $uniqfile = $args{u};
open( UNIQ, ">$uniqfile" ) || die "couild not open $uniqfile for output\n";

my $coverage = $args{c};
if ( ! defined $coverage ) { $coverage = 0.80 }
elsif ( $coverage < 0.3 && $coverage > 1.0 ) { die "\nThe value of -c (coverage) must be between 0.3 and 1.0\n" }
my $similarity = $args{s};
if ( ! defined $similarity ) { $similarity = 0.80 }
elsif ( $similarity < 0.55 && $similarity > 1.0 ) { die "\nThe value of -s (similarity) must be between 0.55 and 1.0\n" }

# create a working directory
my $workdir = dirname( $outfasta );
{
    my $tmpdir = "$workdir/FC_" . rand(1000000);
    while ( -e $tmpdir ) {
        my $tmpdir = "$workdir/FC_" . rand(1000000);
    }
    $workdir = $tmpdir;
}
mkdir $workdir || die "\nCannot create working directory \"$workdir\"\n";

print "working directory: $workdir\n";

# use cd-hit to reduce redundancy ( 95% identity 95 %coverage)
my $cmd = "$myBin/cd-hit -i $infasta -o $workdir/nr -aL $coverage -aS $coverage -G 0 -c $similarity -g 1 -B 1 -t 5 2>&1 | tee $workdir/cluster.log";
print "clustering: $cmd\n";
runCmd('.', $cmd);
print "flushing clustalw buffers...\n";
sleep 60;    # give cd-hit time to flush buffers to disk

# find clusters
print "loading cluster definitions\n";
my %multiclusters;
open( CL, "<$workdir/nr.clstr" ) || die "\nCannot read cluster deflines\n";
my $cluster_id = 0;
my %clusters;
while ( my $member = <CL> ) {
    chomp $member;
    if ( $member =~ /^>Cluster (.*)$/ ) {        
        $cluster_id = $1;
    }
    else {
        $member =~ s/^.*>//;
        $member =~ s/\.\.\..*$//;
        $member =~ s/[\t ].*//;
        push @{$clusters{$cluster_id}{members}}, $member;
    }
}
close CL;

# output clusters
for my $cluster_id ( sort { $a <=> $b } keys %clusters ) {
    $clusters{$cluster_id}{id} = $cluster_id;
    for my $member_id ( @{ $clusters{$cluster_id}{members} } ) {
        if ( ! exists $inseqs{$member_id} ) {
            print "skipping: $member_id\n";
        }
        else {
            my $member = $inseqs{$member_id};
            my $def = clean_defline( $$member{defline} . " length=$$member{seqlen}" );
            $def .= " cluster_id=$cluster_id";
            print OUT ">$def\n";
            print ">$def\n";
            my $seq = $$member{sequence};
            $seq =~ s/(.{60})/$1\n/g;
            if ( $seq !~ /\n$/ ) { $seq .= "\n" }
            print OUT "$seq\n";
            $clusters{$cluster_id}{representative} = $member;
            my $product = get_defline_tag( $$member{defline}, "product" );
            $clusters{$cluster_id}{products}{$product}++;
            my $gene = get_defline_tag( $$member{defline}, "gene" );
            $clusters{$cluster_id}{genes}{$gene}++;
            $clusters{$cluster_id}{size}++;
        }
    }
    print OUT "\n";    
}
close OUT;

# output file with one representative per cluster
for my $c ( sort { $a <=> $b } keys %clusters ) {
    my $cluster = $clusters{$c};
    my @products = sort { $cluster->{products}->{$b} <=> $cluster->{products}->{$a} } keys %{ $cluster->{products} };
    my @genes = sort { $cluster->{genes}->{$b} <=> $cluster->{genes}->{$a} } keys %{ $cluster->{genes} };
    my $seq = $cluster->{representative};
    $seq->{defline} =~ s/ gene="[^"]*"//;
    $seq->{defline} .= " gene=\"$genes[0]\"";
    $seq->{defline} =~ s/ product="[^"]*"//;
    $seq->{defline} .= " product=\"$products[0]\"";
    $seq->{defline} .= " length=$seq->{seqlen} cluster_id=$cluster->{id} cluster_size=$cluster->{size}";
    my $sequence = $seq->{sequence};
    $sequence =~ s/(.{60})/$1\n/g;
    print UNIQ ">$seq->{defline}\n$sequence\n";
}
close UNIQ;

# delete workspce and exit
print "cleaning up\n";
system "rm -rf $workdir";
print "done\n";

exit(0);


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
        if ( $tag =~ /^besthit/i ) { next }
        elsif ( $tag =~ /^original/i ) { next }
        elsif ( $tag =~ /(curated|cluster_size|related_proteins)/i ) { next }
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
