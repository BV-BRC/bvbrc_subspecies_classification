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
set_parameter( "default_gene_variation", 0 );
set_parameter( "selectivity", 1 );
my @commonTags = sort { $a cmp $b } (
    "gene", "product", "note",
    "shared_cds",
    "spliced", 
    "matpepdb",
    "stop_codon_readthru", "alternate_startcodon",
    "ribosomal_slippage",
    "rna_editing"
    );
my @allTags = sort { $a cmp $b } (
    "gene", "product", "note",
    "shared_cds", "excludes_gene",
    "is_required", "is_optional",
    "spliced", "intron_size", "splice_form", "noncanonical_splicing", "tiny_exon3",
    "matpepdb",
    "stop_codon_readthru", "alternate_startcodon",
    "ribosomal_slippage", "slippage_motif", "slippage_frameshift",
    "rna_editing"
    );
    
my $rawFasta = shift @ARGV;
if ( ! -r $rawFasta || ! -e $rawFasta ) {
    die "\nCannot read input fasta \"$rawFasta\"\n";
}
my $annotationRef = shift @ARGV;
if ( ! -r $annotationRef || ! -e $annotationRef ) {
    die "\nCannot read annotation reference \"$annotationRef\"\n";
}
my $coverageCutoff = shift @ARGV;
if ( $coverageCutoff < 30 || $coverageCutoff > 100 ) {
    die "\nInvalid %coverage ($coverageCutoff), expected 30-100\n"
}
my $identityCutoff = shift @ARGV;
if ( $identityCutoff < 30 || $identityCutoff > 100 ) {
    die "\nInvalid %identity ($identityCutoff), expected 30-100\n"
}

my @transferTags;
{
    for my $tag ( @ARGV ) {
        if ( lc $tag eq "common") {
            push @transferTags, @commonTags;    
        }
        elsif ( lc $tag eq "all" ) {
            push @transferTags, @allTags;    
        }
        else {
            push @transferTags, lc $tag ;
        }
    }
    my %uniq;
    for my $tag ( @transferTags ) {
        $uniq{lc $tag} = 1;
    }
    @transferTags = sort { $a cmp $b } keys %uniq;
}
if ( ! @transferTags ) {
    die "\nNo tags specified for transfer\n";
}

my %rawSeqs = loadFasta( $rawFasta );

my $tmpRef = $vigorspace . "/" . basename( $annotationRef );
my $cmd = "cp $annotationRef $tmpRef";
runCmd($vigorspace, $cmd);
$cmd = "formatdb -i $tmpRef -p T -l $vigorspace/formatdb.log";
runCmd($vigorspace, $cmd);
my %templateSeqs = loadFasta( $tmpRef );
set_reference_seqs( %templateSeqs );

# write results, sorted by gene
my $queryFasta = "$vigorspace/query.fasta";
my $blastXML = "$vigorspace/blast.xml";
for my $rawSeq ( sort { $$a{id} cmp $$b{id} } values %rawSeqs ) {
    $$rawSeq{defline} = cleanDefline( $$rawSeq{defline} );
    unlink $queryFasta;
    unlink $blastXML;
    writeFasta( $queryFasta, $rawSeq );
    $cmd = "blastall -p blastp -b 5 -e 0.01 -F \"\" -i $queryFasta -d $tmpRef -m 7 -o $blastXML";
    runCmd($vigorspace, $cmd);
    my @hsps = parse_blastxml( $blastXML, 0 );
    my @hits = sort { $$b{vigor_matchwgt} <=> $$a{vigor_matchwgt} } best_subject_hits( undef, \@hsps );

#print "\n> $$rawSeq{defline}\n";
#print "\n     HSPS\n";
#print_blasthits( 5, @hsps ); 
#print "\n     HITS\n";
#print_blasthits( 5, @hits );

    my %original;
    for my $tag ( @transferTags ) {
        if ( $$rawSeq{defline} =~ /^(.*) $tag=("[^"]*")(.*)$/i ) {
            $original{$tag} = $2;
            $$rawSeq{defline} = "$1$3";
        }
        elsif ( lc( $tag ) eq "rna_editing" && "$$rawSeq{defline} " =~ /^(.*) rna_editing=([0-9]+\/[^\/]+\/[^\/]+\/[^\/]+\/) (.*)/i ) {
            $original{$tag} = "\"$2\"";
            $$rawSeq{defline} = "$1$3";
        }
        elsif (  "$$rawSeq{defline} " =~ /^(.*) $tag=([^ ]*)( .*)/i ) {
            $original{$tag} = $2;
            $$rawSeq{defline} = "$1$3";
        }
        elsif ( $tag =~ /^is_/ && "$$rawSeq{defline} " =~ /^(.*) ($tag)( .*)/i ) {
            $original{$tag} = $2;
            $$rawSeq{defline} = "$1$3";
        }
    }
 
    my %best;
    my $curated = 0;
    my $besthit;
    if ( @hits ) {
        my $hit = shift @hits;
        $besthit = "$$hit{subject_id}/$$hit{pct_scoverage}/$$hit{pct_identity}";
        if ( $$hit{pct_identity} >= $identityCutoff
                && $$hit{pct_scoverage} >= $coverageCutoff
                && $$hit{pct_qcoverage} >= $coverageCutoff ) {
            $curated = 1;
        }
        my $refSeq = get_reference_seq( $$hit{subject_id} );
        for my $tag ( @transferTags ) {
#            if ( $tag eq "rna_editing" && $$refSeq{defline} =~ /rna_editing/ ) {
#                my $value = "";
#                if ( $$refSeq{defline} =~ / rna_editing=(.*\/.*\/.*\/.*\/)/i ) {
#                    $value = $1;
#                }
#                print "defline=$$refSeq{defline}\n$tag=$value\n";
#            }
            if ( $$refSeq{defline} =~ /^(.*) $tag=("[^"]*")(.*)$/i ) {
                $best{$tag} = $2;
            }
            elsif ( lc( $tag ) eq "rna_editing" && "$$refSeq{defline} " =~ /^(.*) rna_editing=(.*\/.*\/.*\/.*\/) (.*)/i ) {
                $best{$tag} = "\"$2\"";
            }
            elsif (  "$$refSeq{defline} " =~ /^(.*) $tag=([^ ]*)( .*)/i ) {
                $best{$tag} = $2;
            }
            elsif ( $tag =~ /^is_/ && "$$refSeq{defline} " =~ /^(.*) ($tag)( .*)/i ) {
                $best{$tag} = $2;
            }
        }
    }
    else {
        $besthit = "no hit";
        %best = %original;
    }
    
    my %annotation;
    my %rejected;
    my $rej_prefix; 
    if ( $curated ) {
        %annotation = %best;
        %rejected = %original;
        $rej_prefix = "original_";
    }
    else {
        %annotation = %original;
        %rejected = %best;
        $rej_prefix = "besthit_";
    }
    
    for my $tag ( @transferTags ) {
        if ( defined $annotation{$tag} ) {
            if ( $tag =~ /^is_/ ) {
                $$rawSeq{defline} .= " $tag";
            }
            else {
                $$rawSeq{defline} .= " $tag=$annotation{$tag}";
            }
        }
    }
    if ( $curated ) {
        $$rawSeq{defline} .= " curated=Y besthit=\"$besthit\"";
    }
    else {
        $$rawSeq{defline} .= " curated=N besthit=\"$besthit\"";
    }
    for my $tag ( @transferTags ) {
        if ( defined $rejected{$tag} ) {
            if ( ! defined $annotation{$tag} || $annotation{$tag} ne $rejected{$tag} ) {
                if ( $tag =~ /^is_/ ) {
                    $$rawSeq{defline} .= " $rej_prefix$tag";
                }
                else {
                    $$rawSeq{defline} .= " $rej_prefix$tag=$rejected{$tag}";
                }
            }
        }
        elsif ( defined $annotation{$tag} ) {
            $$rawSeq{defline} .= " $rej_prefix$tag=\"\"";
        }
    }
}

set_reference_seqs( %rawSeqs );
my @seqs = sort { get_reference_gene( $$a{id} ) cmp get_reference_gene( $$b{id} ) } values %rawSeqs; 
for my $seq  ( @seqs ) {
    my $sequence = $$seq{sequence};
    $sequence =~ s/(.{60})/$1\n/g;    
    if ( substr( $sequence, length( $sequence ) - 1  ) ne "\n" ) { $sequence .= "\n" }
    
    print ">$$seq{defline}\n$sequence\n";
}
exit(0);

sub writeFasta {
    my ( $fasta, @seqs ) = @_;
    unlink $fasta;
    open( FASTA, ">$fasta" ) || die "\nCould not write to \"$fasta\"\n";
    for my $seq ( @seqs ) {
        my $sequence = $$seq{sequence};
        $sequence =~ s/(.{60})/$1\n/g;
        if ( substr( $sequence, length( $sequence ) - 1  ) ne "\n" ) { $sequence .= "\n" }
        print FASTA ">$$seq{defline}\n$sequence\n";
    }
    close FASTA;
}

sub cleanDefline {
    my ( $defline ) = @_;

    $defline = " $defline ";
    while ( $defline =~ / curated=[YN] / ) {
        $defline =~ s/curated=[YN] / /;
    }
    while ( $defline =~ / besthit[^=]*="[^"]+" / ) {
         $defline =~ s/ besthit[^=]*="[^"]+" / /;
    }
    while ( $defline =~ / besthit[^=]*=[^ ]+ / ) {
         $defline =~ s/ besthit[^=]*=[^ ]+ / /;
    }
    while ( $defline =~ / original[^=]+="[^"]*" / ) {
        $defline =~ s/ original[^=]+="[^"]*" / /;
    }
    while ( $defline =~ / original[^=]+=[^ ]* / ) {
        $defline =~ s/ original[^=]+=[^ ]* / /;
    }
    $defline =~ s/^ +//;
    $defline =~ s/ +$//;
    $defline =~ s/  +/ /g;
    
    return $defline;
}