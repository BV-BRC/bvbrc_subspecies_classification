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
use Getopt::Std;
use Cwd 'realpath';
use File::Basename;

my $program = realpath($0);
our $myBin = dirname( dirname($program) );
our $myData = "$myBin/data";
our $myConf = "$myBin/conf";
require "$myBin/VIGOR3.pm";

set_default_splice_pairs();
my %splicePairs = %{get_splice_pairs()};

my $mature_peptides = 0;
if ( $ARGV[0] eq "-m" ) {
    $mature_peptides = 1;
    shift @ARGV;
}

# for each genbank file
while ( @ARGV ) {

# get features from genbank file (gene, protein, mat_peptide)
    my $gbFile = shift @ARGV;
    print STDERR "$gbFile\n";
    
    my $genbank = Bio::SeqIO->new( -file => "<$gbFile", -format => "Genbank" );
    my $virus = $genbank->next_seq;
    my $genome;
    $$genome{id} = $virus->display_id;
    $$genome{sequence} = $virus->seq();
    $$genome{seqlen} = length( $$genome{sequence} );

    my $gppfx = $$genome{id};
    $gppfx =~ s/\..*$//;
    my $gpno = 0;
    
    my @classification = $virus->species->classification;
    $$genome{organism} = shift @classification;
         
    my @genes;
    my @proteins;
    my @matpeps;
    my %protgene;
    for my $feature ( $virus->get_SeqFeatures ) {
        my $locstr = $feature->location->to_FTstring();
        if ( $locstr =~ /[><]/ ) { next }

        if ( $feature->primary_tag eq "gene" ) {
            push @genes, $feature;
        }
        else {
            if ( $feature->has_tag( "translation" ) ) {
                if ( $feature->primary_tag eq "CDS" ) {
                    push @proteins, $feature;
                }
                elsif ( $feature->primary_tag eq "mat_peptide" ) {
                    my %mp;
                    $mp{matpep} = $feature;
                    $mp{polypep} = $proteins[@proteins-1];
                    push @matpeps, \%mp;
                }
            }
        }
    }

    # order genes and make sure each has a unique name
    @genes = sort { $a->location->strand * $a->location->start <=> $b->location->strand * $b->location->start } @genes;
    my %geneNames;
    for my $gene ( @genes ) {
        if ( $gene->has_tag( "gene" ) ) {
            for my $value ( $gene->get_tag_values( "gene" ) ) {
                $geneNames{$value} = $gene;
            }
        }
        else {
            my $locus = "";
            if ( $gene->has_tag( "locus_tag" ) ) {
                for my $value ( $gene->get_tag_values( "locus_tag" ) ) {
                    $locus = $value;
                    #print STDERR "locus=$locus\n";
                    my @tmp = split /[-_]/, $locus;
                    $locus = pop @tmp;
                    #print STDERR "  locus=$locus\n";
                    if ( length( $locus ) > 0 && length( $locus ) <= 8 ) {
                        $gene->add_tag_value( "gene", $locus );
                        $geneNames{$locus} = $gene;
                        #print STDERR "  gene=$locus\n";
                    }
                    last;
                }
            }
        }
    }

    my @tmp;
    for my $gene ( @genes ) {
        if ( $gene->has_tag( "gene" ) ) {  push @tmp, $gene }
    }
    @genes = @tmp;
    
    # assign unique gene id to each splice variant
    @proteins = sort { $a->location->start <=> $b->location->start } @proteins;
#    my $geneNo = 0;
    for my $protein ( @proteins ) {
        $gpno++;
        
        if ( ! $protein->has_tag( "protein_id" ) ) {
            $protein->add_tag_value( "protein_id", $gppfx . "-" . $gpno );
        }
        # associate protein with gene
        if ( ! $protein->has_tag( "gene" ) ) {
            for my $gene ( @genes ) {
                if ( $gene->location->strand == $protein->location->strand
                        && $gene->location->contains($protein->location ) ) {
                    my @tmp = $gene->get_tag_values( "gene" );
                    $protein->add_tag_value( "gene", $tmp[0] );
                    last;
                }
            }
        }
        if ( ! $protein->has_tag( "gene" ) ) {
            my $product = "";
            if ( $protein->has_tag( "product" ) ) {
                for my $value ( $protein->get_tag_values( "product" ) ) {
                    $product = $value;
                    last;
                }
            }
            #print STDERR "product=$product\n";
            $product =~ s/ *[^ ]*protein[^ ]* *//gi;
            $product =~ s/hypothetical//gi;
            $product =~ s/non[- ]*structural/NSP/gi;
            $product =~ s/ //g;
            #print STDERR "  product=$product\n";
            if ( length( $product ) > 0 && length( $product ) <= 6 ) {
                $protein->add_tag_value( "gene", $product );
                #print STDERR "  gene=$product\n";
                $geneNames{$product} = $protein;
            }
        }
        if ( ! $protein->has_tag( "gene" ) ) {
            my $geneId = $$genome{id} . "p$gpno";
            $protein->add_tag_value( "gene", $geneId );
            #print STDERR "  gene=$geneId\n";
            $geneNames{$geneId} = $protein;
        }            

        my $geneid = "";
        for my $value ( $protein->get_tag_values( "gene" ) ) {
            $geneid = $value;
            last;
        }
        if ( exists $protgene{$geneid} ) {
            my $variant = 2;
            while ( exists $protgene{"$geneid-$variant"} ) {
                $variant++;
            }
            $geneid = "$geneid-$variant";
            $protein->remove_tag( "gene" );
            $protein->add_tag_value( "gene", $geneid );            
        }
        $protgene{$geneid} = $protein;
    }
    
    # look for shared CDS, noncanonical splicing, and intron sizing
    my %sharedCDS;
    for my $i ( 0..@proteins-1 ) {
        my $prot1 = $proteins[$i];
        if ( ! $prot1->has_tag( "gene" ) ) { next }
        my @tmp1 = $prot1->get_tag_values( "gene" );
        my $id1 = shift @tmp1;
#print "I=$id1\n";
        my $exons1 = feature_exons( $genome, $prot1, 1 );
        for my $j ( 0..@proteins-1 ) {
            if ( $j == $i ) { next }
            my $prot2 = $proteins[$j];
            if ( ! $prot2->has_tag( "gene" ) ) { next }
            my @tmp2 = $prot2->get_tag_values( "gene" );
            my $id2 = shift @tmp2;
#print "J=$id2\n";
            if ( $prot2->location->strand ne $prot1->location->strand ) { next }
            my $exons2 = feature_exons( $genome, $prot2, 2 );
            
            if ( exons_overlap( $exons1, $exons2 ) >= 60 ) {
                $sharedCDS{$id1}{$id2} = 1;
                $sharedCDS{$id2}{$id1} = 1;
            }
        }
    }
    
    # output protein fasta
    if ( ! $mature_peptides ) {
        for my $protein ( sort { order_proteins( $a, $b ) } @proteins ) {
        
            if ( $protein->location =~ /[><]/ ) { next }
            if ( ! $protein->has_tag( "protein_id" ) ) { next }

            my @tmp = $protein->get_tag_values( "protein_id" );
            my $id = shift @tmp;
            
            # check for special cases
            my ( $min_intron, $max_intron );
            my $riboSlippage = "N";
            if ( $protein->has_tag( "ribosomal_slippage" ) ) {
                $riboSlippage = "Y";
            }
            elsif (  $protein->has_tag( "note" ) ) {
                my $note = " ";
                for my $value ( $protein->get_tag_values( "note" ) ) {
                    $note .= $value . " ";
                }
                if ( $note =~ /ribosomal/i && ( $note =~ /frame[ -]*shift/i || $note =~ / slipp/ ) ) {
                    $riboSlippage = "Y";
                }
            }
            my $stopReadThru = "N";
            if ( $protein->has_tag( "transl_except" ) && $protein->has_tag( "note" ) ) {
                my $note = " ";
                for my $value ( $protein->get_tag_values( "note" ) ) {
                    $note .= $value . " ";
                }
                if ( $note =~ /stop/i && $note=~/reads* *thr/ ) {
                    $stopReadThru = "Y";
                }
            }
            #my $rnaEditing = "N";    #not parseable from genbank files
            my $spliceform;
            my $noncanon;
            my $splicing = "N";
            if ( $protein->location->isa('Bio::Location::SplitLocationI') ) {
                if ( $riboSlippage eq "N" ) { $splicing = "Y" }
                my $tmp = feature_exons( $genome, $protein, 0 );
                my @exons = @$tmp;
                $spliceform = "e"
                    . ( $exons[0]{strand} * ( $exons[0]{end} - $exons[0]{start} ) + 1 );
                for my $intron ( 1..@exons-1 ) {
                    my $intron_size =  $exons[$intron]{strand}
                        * ( $exons[$intron]{start} - $exons[$intron-1]{end} ) - 1;
                    $spliceform .= "i$intron_size"
                        . "e"
                        . ( $exons[$intron]{strand} * ( $exons[$intron]{end} - $exons[$intron]{start} ) + 1 );
                        
                    if ( ! defined $min_intron || $intron_size < $min_intron ) { $min_intron = $intron_size }
                    if ( ! defined $max_intron || $intron_size > $max_intron ) { $max_intron = $intron_size }            
                }
                if ( $splicing eq "Y" && @exons > 1 ) {
                    my %nc_splices;
#print_hash( "genome", $genome );
                    for my $i ( 1..@exons-1 ) {
                        my $s = $exons[$i-1]{end} + $exons[$i-1]{strand};
                        my $e = $exons[$i-1]{end} + 2 * $exons[$i-1]{strand};
                        my $donor = subsequence( $$genome{sequence}, $s, $e );
#print "intron start[$i]=$exons[$i-1]{end}  s-e=$s-$e  donor=$donor\n";

                        $s = $exons[$i]{start} - 2 * $exons[$i]{strand};
                        $e = $exons[$i]{start} - $exons[$i]{strand};
                        my $acceptor = subsequence( $$genome{sequence}, $s, $e );
#print "intron end[$i]=$exons[$i]{end}  s-e=$s-$e  acceptor=$acceptor\n";

                        if ( ! exists $splicePairs{$donor}{$acceptor} ) {
                            $nc_splices{$donor}{$acceptor} = 1
                        }
                    }
                    for my $donor ( keys %nc_splices ) {
                        for my $acceptor ( keys %{ $nc_splices{$donor} } ) {
                            if ( defined $noncanon ) {
                                $noncanon .= ";$donor+$acceptor";
                            }
                            else {
                                $noncanon = "$donor+$acceptor";
                            }
                        }
                    }
                }
            }
            
            my $geneid = "";
            for my $value ( $protein->get_tag_values( "gene" ) ) {
                if ( length( $geneid ) ) { $geneid .= " " }
                $geneid .= $value;
            }
            
            my $product = "";
            if ( $protein->has_tag( "product" ) ) {
                for my $value ( $protein->get_tag_values( "product" ) ) {
                    if ( length( $product ) ) { $product .= " " }
                    $product .= $value;
                }
            }
            else {
                $product = "hypothetical protein";
            }
            
            my $shared_cds;
            if ( exists $sharedCDS{$geneid} ) {
                $shared_cds = join( ",", keys %{$sharedCDS{$geneid}} );
            }
            
            my $matpepdb;
            if ( $product =~ /polyprotein/i ) {
                $matpepdb = "default";
            }
        
            # format defline            
            my $defline = ">$id gene=\"$geneid\" product=\"$product\"";
            if ( $splicing eq "Y" ) {
                $defline .= " spliced=Y splice_form=\"$spliceform\"";
                $min_intron -= 50;
                if ( $min_intron < 40 ) { $min_intron = 40 }
                $max_intron += 50;
                $defline .= " intron_size=$min_intron-$max_intron";
                if ( defined $noncanon ) {
                    $defline .= " noncanonical_splicing=\"$noncanon\"";
                }
            }
            if ( $riboSlippage eq "Y" ) { $defline .= " ribosomal_slippage=Y splice_form=\"$spliceform\"" }
            if ( $stopReadThru eq "Y" ) { $defline .= " stopcodon_readthru=Y" }
            if ( defined $shared_cds ) { $defline .= " shared_cds=\"$shared_cds\"" }
            if ( defined $matpepdb ) { $defline .= " matpepdb=\"$matpepdb\"" }
            
            $defline .= " organism=\"$$genome{organism}\"";
            
#print "$geneid  $id  $product\n";        
            # format sequence
            my $aa = "";
            for my $value ( $protein->get_tag_values( "translation" ) ) {
                $aa .= $value;
            }
            if ( $aa =~ /xxxxx/i ) { next }
            $aa =~ s/(.{60})/$1\n/g;
            if ( substr( $aa, length( $aa ) - 1 ) ne "\n" ) { $aa .= "\n" }
            
            # write entry    
            print "$defline\n$aa\n";
        }
    }
    
    # output mat_peptide fasta
    else {
        my %mpref;
        for my $mp ( @matpeps ) {
            my $matpep = $$mp{matpep};
            my $polypep = $$mp{polypep};
            $mpref{$polypep}++;
        
            # format defline    
            my $id = "";
            if ( $matpep->has_tag( "protein_id" ) ) {
                for my $value ( $matpep->get_tag_values( "protein_id" ) ) {
                    if ( length( $id ) ) { $id .= " " }
                    $id .= $value;
                }
            }
            else {
                for my $value ( $polypep->get_tag_values( "protein_id" ) ) {
                    if ( length( $id ) ) { $id .= " " }
                    $id .= $value;
                }
                $id .= "_MP$mpref{$polypep}";
            }
            
            my $geneid = "";
            for my $value ( $polypep->get_tag_values( "gene" ) ) {
                if ( length( $geneid ) ) { $geneid .= " " }
                $geneid .= $value;
            }
            
            my $product = "";
            if ( $matpep->has_tag( "product" ) ) {
                for my $value ( $matpep->get_tag_values( "product" ) ) {
                    if ( length( $product ) ) { $product .= " " }
                    $product .= $value;
                }
            }
            else {
                $product = $id;
            }
            
            my $defline = ">$id gene=\"$geneid\" product=\"$product\"";
        
            # format sequence
            my $polyaa = "";
            for my $value ( $polypep->get_tag_values( "translation" ) ) {
                $polyaa .= $value;
            }
            
            my $polystart = $polypep->location->start;
            my $polyend = $polypep->location->end;
            my $mpstart = $matpep->location->start;
            my $mpend = $matpep->location->end;
            my $start_offset = abs( $mpstart - $polystart ) / 3;
            my $end_offset = abs( $polyend - $mpend ) / 3 - 1;
    #print "$defline\n";
    #print "poly $polystart-$polyend\n";
    #print "poly $mpstart-$mpend\n";
    #print "off +$start_offset -$end_offset\n";
            
            my $mpaa = substr( $polyaa, 0, length( $polyaa ) - $end_offset );
            $mpaa = substr( $mpaa, $start_offset );
            if ( $mpaa =~ /xxxxx/i ) { next }
            
    #print "poly aa=$polyaa\n";
    #print "  mp aa=$mpaa\n";
            $mpaa =~ s/(.{60})/$1\n/g;
            if ( substr( $mpaa, length( $mpaa ) - 1 ) ne "\n" ) { $mpaa .= "\n" }
                
            # write entry    
            print "$defline\n$mpaa\n";
        }
    }
}
exit(0);

sub feature_exons {
    my ( $genome, $feature, $flag ) = @_;
    if ( ! defined $flag ) { $flag = -1 }

#   for my $tag ($feature->get_all_tags) {             
#      print "  tag: ", $tag, "\n";             
#      for my $value ($feature->get_tag_values($tag)) {                
#         print "    value: ", $value, "\n";             
#      }          
#   }       

    my @tmp = $feature->get_tag_values( "protein_id" );
    my $id = shift @tmp;
#print "FEATURE $id\n";
    
    my $frameadj = 0;
    if ( $feature->has_tag( "codon_start" ) ) {
        my @tmp = $feature->get_tag_values( "codon_start" );
        $frameadj = shift @tmp;
        $frameadj--;
    }

    my @locations;
    if ( $feature->location->isa('Bio::Location::SplitLocationI') ) {
        for my $location ( $feature->location->sub_Location ) {
            push @locations, $location;
        }
    }
    else {
        push @locations, $feature->location;
    }
    @locations = sort { $a->strand * $a->start <=> $b->strand * $b->start } @locations;
    
    my @exons;
    for my $location ( @locations ) {

        my $start = $location->start;
        my $end = $location->end;
        my $strand = $location->strand;
        if ( $strand == -1 ) {
            ( $start, $end ) = ( $end, $start );
        }
        
#print "frameadj=$frameadj  location=$start-$end ($strand)\n" ;
        my $exon;
        $$exon{start} = $start;
        $$exon{codon_start} = $start + $frameadj;
        $$exon{end} = $end;
        $$exon{codon_end} = $end;
        $$exon{strand} = $strand;
        if ( $strand == -1 ) {
            ( $$exon{codon_left}, $$exon{codon_right} ) = ( $$exon{codon_end}, $$exon{codon_start} );
            $$exon{frame} = ( $$genome{seqlen} + 1 - $$exon{codon_start} ) % 3;
            $$exon{frame} = -$$exon{frame};
            if ( $$exon{frame} == 0 ) { $$exon{frame} = -3 }
        }
        else {
            ( $$exon{codon_left}, $$exon{codon_right} ) = ( $$exon{codon_start}, $$exon{codon_end} );
            $$exon{frame} = $$exon{codon_start} % 3;
            if ( $$exon{frame} == 0 ) { $$exon{frame} = 3 }
        }

        $frameadj = -( $$exon{codon_right} - $$exon{codon_left} + 1 ) % 3;
        if ( $frameadj < 0 ) { $frameadj += 3 }

#print "exon=$$exon{codon_left}-$$exon{codon_right} ($$exon{frame})\n";
        while( ( $$exon{codon_right} - $$exon{codon_left} + 1 ) % 3 > 0 ) {
            $$exon{codon_end} -= $strand;
            if ( $strand == 1 ) {
                $$exon{codon_right} -= 1;
            }
            else {
                $$exon{codon_left} += 1;
            }
#print "exon=$$exon{codon_left}-$$exon{codon_right} ($$exon{frame})\n";
        }
        push @exons, $exon;
    }

    return \@exons;    
}

sub order_proteins {
    my ( $a, $b ) = @_;

    my @tmp = $a->get_tag_values( "gene" );
    my $genea = shift @tmp;

    @tmp = $b->get_tag_values( "gene" );
    my $geneb = shift @tmp;

    if ( $genea ne $geneb ) { return $genea cmp $geneb }
    
    return $a->location->strand * $b->location->start <=> $a->location->strand * $b->location->start; 
}

sub exons_overlap {
    my ( $exons1, $exons2 ) = @_;
    
    my $overlap = 0;
#print "start overlap\n";
    for my $exon1 ( @$exons1 ) {
        for my $exon2 ( @$exons2 ) {
            if ( $$exon1{frame} == $$exon2{frame} ) {
                my $overlapR = minval( $$exon1{codon_right}, $$exon2{codon_right} );
                my $overlapL = maxval( $$exon1{codon_left}, $$exon2{codon_left} );
                if ( $overlapR >= $overlapL ) {
                    $overlap += ( $overlapR - $overlapL + 1 );
#print "$$exon1{codon_left}-$$exon1{codon_right} ($$exon1{frame}) | $$exon2{codon_left}-$$exon2{codon_right} ($$exon2{frame}) | $overlapL-$overlapR | $overlap\n";
                }
            }
        }
    }
    return $overlap;
}