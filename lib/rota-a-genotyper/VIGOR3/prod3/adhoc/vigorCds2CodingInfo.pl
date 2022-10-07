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

my @part_list = ( "I", "II", "III", "IV", "V", "VI", "VII", "VIII", "IX", "X" );
my %transcript_labels;
my %transcripts;
my @matpeps;

# load CDS file (fasta format)
my $fastafile = shift @ARGV;
if ( ! -r $fastafile ) { die "\ncannot read fasta $fastafile\n" }
my %fasta = loadFasta( $fastafile );

# parse defline
for my $cds ( values %fasta ) {
    my ( $cdsid, $deftags ) = parse_defline( $$cds{defline} );
#print "$$cds{defline}\n";

    # skip pseudogenes
    if ( exists $$deftags{pseudogene} ) { next }
    $$cds{tags} = $deftags;
        
    # skip partial genes
    my $location = $$cds{tags}{location};
    $location =~ s/[><]//g;

    # convert location into exons;
#print "  location=$location\n";
    my @exons;
    my ( $left, $right, $ori );
    my $exon_num = 0;
    for my $exon ( split /,/, $location ) {
        my ( $start, $stop ) = split /\.\./, $exon;
        my $x;
        $exon_num++;
        $$x{exon_num} = $exon_num;
        $$x{left} = $start;
        $$x{right} = $stop;
        $$x{ori} = 1;
        push @exons, $x;
        
        # update span
        if ( ! defined $left || minval( $start, $stop ) < $left ) {
            $left = minval( $start, $stop );
        }
        if ( ! defined $right || maxval( $start, $stop ) > $right ) {
            $right = maxval( $start, $stop );
        }
    }
        
    # adjust exons reverse orientated transcript
    if ( $left != $exons[0]{left} ) {
        for my $x ( @exons ) {
            $$x{ori} = -1;
            ( $$x{left}, $$x{right} ) = ( $$x{right}, $$x{left} );
        }
        $ori = -1;
    }
    else {
        $ori = 1;
    }

    # convert coordinates to 0-space (from 1-base)
    for my $x ( @exons ) {
        $$x{left}--;
#print "    exon $$x{exon_num}: $$x{left}-$$x{right} ($$x{ori})\n";
    }
    $left--;

    # adjust for partial codons
    my $cdslen = $$cds{seqlen};
    if ( defined $$cds{tags}{codon_start} && $$cds{tags}{codon_start} > 1 ) {
#print "codon_start=$$cds{tags}{codon_start}  left-right=$left-$right  cdslen=$cdslen\n";
        my $adjust5 = $$cds{tags}{codon_start} - 1;
        $cdslen -= $adjust5;
        if ( $ori == 1 ) {
            $left += $adjust5;
            $exons[0]{left} = $left; 
        }
        else {
            $right -= $adjust5;
            $exons[0]{right} = $right;
        }
#print "   adjust5=$adjust5 => $left-$right  cdslen=$cdslen\n";
    }
    if ( $cdslen % 3 > 0 ) {
        my $adjust3 = $cdslen % 3;
        if ( $ori == 1 ) {
            $right -= $adjust3;
            $exons[@exons-1]{right} = $right;
        }
        else {
            $left += $adjust3;
            $exons[@exons-1]{left} = $left;
        }
        $cdslen -= $adjust3;
#print "   adjust3=$adjust3 => $left-$right\n";
    }

    # save span and exons
    $$cds{left} = $left;
    $$cds{right} = $right;
    $$cds{ori} = $ori;
    $$cds{cdslen} = $cdslen;
#print "    span $left-$right ($ori)\n";
    $$cds{exons} = \@exons;

    # save entry
    if ( exists $$deftags{mat_peptide} ) {
#print "  save as mat_peptide\n";
        my @tmp = split /\./, $cdsid;
        pop @tmp;
        my $transcript_id = join( ".", @tmp );
#print "cdsid=$cdsid  transcript_id=$transcript_id\n";
        $$cds{transcript_id} = $transcript_id;
        $$cds{tags}{product} =~ s/^mature peptide, //;
        push @matpeps, $cds;
    }
    else {
#print "  save as transcript $$cds{tags}{gene}\n";
        my $transcript_label = $$cds{tags}{gene};
        $$cds{transcript_label} = $transcript_label;
        $transcript_labels{$transcript_label}{$cdsid} = $cds;
        $$cds{transcript_id} = $cdsid;
        $transcripts{$cdsid} = $cds;
    }
}

# adjust names of split transcripts to gurantee uniquenes
for my $name ( keys %transcript_labels ) {
    my @parts = values %{ $transcript_labels{$name} };
    if ( @parts > 1 ) {
        my $part_num = 0;
        for my $transcript ( sort {  $$a{ori} * $$a{left} <=> $$b{ori} * $$b{left} } @parts ) {
            $$transcript{transcript_label} .= "-$part_list[$part_num]";
            $part_num++;
        }
    }
}

# convert matpep coordinates to peptide, assign a label, and associate with parent transcript
for my $mp ( sort { $$a{ori} * $$a{left} <=> $$b{ori} * $$b{left} } @matpeps ) {
    my $tid = $$mp{transcript_id};
    my $transcript = $transcripts{$tid};
    
#print "  mp $$mp{tags}{gene} $$mp{left}-$$mp{right} ($$mp{ori})\n";
#print "    transcript $$transcript{tags}{gene} $$transcript{left}-$$transcript{right} ($$transcript{ori})\n";
    my @mpexons = @{ $$mp{exons} };
    my @texons = @{ $$transcript{exons} };
    
    my $offset = 0;
    if ( $texons[0]{ori} == 1 ) {
        while ( $texons[0]{right} < $mpexons[0]{left} ) {
            $offset += ( $texons[0]{right} - $texons[0]{left} );
            shift @texons;
        }
        $offset += ( $mpexons[0]{left} - $texons[0]{left} );
    }
    else {
        while ( $texons[@texons-1]{left} > $mpexons[@mpexons-1]{right} ) {
            $offset += ( $texons[@texons-1]{right} - $texons[@texons-1]{left} );
            pop @texons;
        }
        $offset += ( $texons[@texons-1]{right} - $mpexons[@mpexons-1]{right} );
    }
    $$mp{left} = $offset / 3;
    $$mp{right} = $$mp{left} + $$mp{seqlen} / 3;
#print "    adjusted $$mp{left}-$$mp{right}\n";
    push @{ $$transcript{matpeps} }, $mp;
    $$mp{matpep_label} = "$$transcript{transcript_label}.p". @{ $$transcript{matpeps} };
}

# output information;
for my $transcript ( sort { $$a{left} <=> $$b{left} } values %transcripts ) {
    my ( $gene_id ) = split /-/, $$transcript{tags}{gene};
    for my $exon ( sort{ $$a{ori} * $$a{left} <=> $$b{ori} * $$b{left} } @{ $$transcript{exons} } ) {
        print "EXON\t$gene_id\t$$transcript{transcript_label}\t$$transcript{transcript_label}.x$$exon{exon_num}\t$$exon{left}\t$$exon{right}\t$$exon{ori}\n";
    }
    my $aalen = $$transcript{cdslen} / 3;
    if ( $$transcript{tags}{location} !~ />[0-9]+$/ ) { $aalen-- }
    print "CDS\t$gene_id\t$$transcript{transcript_label}\t$$transcript{left}\t$$transcript{right}\n";
    print "PROTEIN\t$$transcript{transcript_label}\t$$transcript{transcript_label}\t$$transcript{tags}{product}\t0\t$aalen\n";
    if ( defined $$transcript{matpeps} ) {
        for my $mp ( sort { $$a{left} <=> $$b{left} } @{ $$transcript{matpeps} } ) {
            print "PROTEIN\t$$transcript{transcript_label}\t$$mp{matpep_label}\t$$mp{tags}{product}\t$$mp{left}\t$$mp{right}\n";
        }
    }
}
exit( 0 );

sub parse_defline {
    my ( $defline ) = @_;
    
    my %tags;
    my @misc;
    
    my $tmp = $defline;
    
    my $id = $tmp;
    $id =~ s/ .*//;
    $tmp =~ s/^[^ ]* *//;
    
    while ( $tmp ) {
        if ( $tmp =~ /^([^= ]+)="([^"]*)"/ ) {
            my ( $tag, $value ) = ( $1, $2 );
            if ( ! exists $tags{$tag} ) { $tags{$tag} = $value }
            $tmp =~ s/^([^= ]+)="([^"]*)" *//;
        }
        elsif ( "$tmp " =~ /^([^= ]+)=([^ ]*) / ) {
            my ( $tag, $value ) = ( $1, $2 );
            if ( ! exists $tags{$tag} ) { $tags{$tag} = $value }
            $tmp =~ s/^([^= ]+)=([^ ]*) *//;
        }
        elsif ( "$tmp " =~ /^([^ ]*) / ) {
            my $tag = $1;
            $tags{$tag} = undef;
            $tmp =~ s/^([^ ]*) *//;
        }
    }
    
    for my $tag ( sort { $a cmp $b } keys %tags ) {
        if ( $tag =~ /^(curated|besthit|original|cluster|related)/i ) { delete $tags{$tag} }
        elsif ( defined $tags{$tag} ) {
            $tags{$tag} =~ s/^ +//;
            $tags{$tag} =~ s/ +$//;
            $tags{$tag} =~ s/  +/ /g;
        }
    }
    if ( @misc ) {
        $tmp .= " " . join( " ", @misc );
    }
    
    return ( $id, \%tags );
}