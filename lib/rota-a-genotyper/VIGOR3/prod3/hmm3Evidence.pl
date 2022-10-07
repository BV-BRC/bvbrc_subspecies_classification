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
our $myBin = dirname($program);
require "$myBin/VIGOR3.pm";
our $myData = dirname( $myBin ) . "/data3";

*STDERR = *STDOUT;
$|++;

# get command line parameters
my %args;
&getopts( 'p:h:v:', \%args );
my ( $opt_p, $opt_h, $opt_v) =    ( $args{p}, $args{h}, $args{v} );

if ( ! -e $opt_p || ! -r $opt_p ) {
    die "\nCannot read peptide fasta -p \"$opt_p\"\n";
}

my $cmd = "touch $opt_h";
runCmd('.', $cmd);

if ( -d $opt_h || ! -w $opt_h) {
    die "\nCannot read/write htab output -h \"$opt_h\"\n";
}
else {
    unlink $opt_h;
}


if ( defined $opt_v ) {
    if ( ! -e $opt_v || ! -w $opt_v ) {
        die "\nCannot read/write VIGOR tbl -v \"$opt_v\"\n";
    }

    # save original tbl file
    my $saveId = 0;
    my $saveFile = "$opt_v.save$saveId";
    while ( -e $saveFile ) {
        $saveId++;
        $saveFile = "$opt_v.save$saveId";
    }
    my $cmd = "cp $opt_v $saveFile\n";
    runCmd('.', $cmd);
    print "original TBL file saved to $saveFile\n";
}


# create workspace on scratch
my $workspace = create_workspace();

# load query sequences
my %proteins = loadFasta( $opt_p );

# run HMM3 search for each protein
my $evidenceHits;
for my $seq ( sort { compare_ids( $$a{id}, $$b{id} ) } values %proteins ) {
    if ( $$seq{defline} =~ /(mat_peptide|pseudogene)/i ) { next }
    my $qryfasta = "$workspace/qry.fasta";
    unlink $qryfasta;
    open QRY, ">$qryfasta" || die "\nCould not write query fasta \"$qryfasta\"\n";
    print QRY ">$$seq{defline}\n$$seq{sequence}\n";
    close QRY;
    
    my $hmm3out = "$workspace/hmm3.out";
    unlink $hmm3out;
    
    my $hmm3log = "$workspace/hmm3.log";
    unlink $hmm3log;

    my $cmd = "/usr/local/packages/hmmer-3.0/bin/hmmscan --cut_nc -o $hmm3out /usr/local/db/common/ergatis/PFAM_TIGR/current/ALL_LIB.HMM $qryfasta &> $hmm3log";
    runCmd($workspace, $cmd);

# format htab results
    my $htablog = "$workspace/htab.log";
    unlink $htablog;

    my $htabout = "$workspace/htab.out";
    unlink $htabout;

    $cmd = "/usr/local/devel/ANNOTATION/ard/ergatis-v2r13b2/bin/hmmpfam2htab --input_file $hmm3out --mldbm_file /usr/local/db/common/ergatis/PFAM_TIGR/current/ALL_LIB.HMM.db --output_htab $workspace/htab.tmp # &> $htablog";
    runCmd($workspace, $cmd);
    open( TMP, "<$workspace/htab.tmp" ) || die "\nCould not read $workspace/htab.tmp\n";
    open( HTAB, ">$htabout" ) || die "\nCould not write $htabout\n";
    for my $line ( <TMP> ) {
        my @data = split /\t/, $line;
        if ( $data[12] >= $data[17] ) { print HTAB $line }
    }
    close TMP;
    close HTAB;
    $cmd = "cat $htabout >> $opt_h";
    runCmd($workspace, $cmd);
    
# extact results for this protein and add to final results
    $cmd = "cut -f1,16 $htabout | sort | uniq";
    my $results = runCmdAndGetResults($workspace, $cmd);
    my @hits = split(/\n/, $results);
    
    for my $hit ( @hits ) {
        my @tmp = split /\t/, $hit;
        my $result;
        print "$$seq{id}\t$hit\n";
        my $hmm = $tmp[0];
        my $fam = $tmp[1];
        $hmm =~ s/ //g;
        $fam =~ s/^ +//;
        if ( $fam =~ /: / ) {
            $fam =~ s/: .*$//;
        }
        else {
            $fam =~ s/ .*//;
        }
        if ( $fam =~ /^ *$/ ) { $fam = $tmp[1] }
        my $note = "$fam ($hmm)";
        if ( exists $$evidenceHits{$$seq{id}} ) {
            push @{ $$evidenceHits{$$seq{id}} }, $note;
        }
        else {
            my @tmphits = ( $note );
            $$evidenceHits{$$seq{id}} = \@tmphits;
        }
    }
}

# update VIGOR table file
if ( defined $opt_v ) {    
    my @tblinfo = parseTBL( $opt_v );
    for my $protein_id ( keys %$evidenceHits ) {
        my $note = join( ", ", sort { $a cmp $b } @{ $$evidenceHits{$protein_id} } );
        if ( index( $note, "," ) < 0 ) {
            $note = "identified by match to protein family $note";
        }
        else {
            $note = "identified by matches to protein families $note";
            $note =~ s/, ([^,]+)$/ and $1/;
        }
        my $feature = findProtein( \@tblinfo, $protein_id );
        if ( defined $feature ) {
            if ( defined $$feature{note} ) {
                $$feature{note} .= "; $note";
            }
            else {
                $$feature{note} = $note;
            }
        }
    } 
    
    rewriteTBL( \@tblinfo, $opt_v );
}

# delete workspace
#system "rm -rf $workspace\n";
exit(0);

sub parseTBL {
    my ( $tblFile ) = @_;
    
    my @tblinfo;
    my $cmd = "cat $tblFile | sed 's/\\r/\\n/g'";
    my $tblData = runCmdAndGetResults('.', $cmd);
    my @seqs = split />Features /, $tblData;
    shift @seqs;
    while ( @seqs ) {
        my $seqbody = shift @seqs;
        chomp $seqbody;
        my @seqdata = split /(([><]*[0-9]+)\t([><]*[0-9]+)\t([^\n]+))\n/, $seqbody;
        my $seqid = shift @seqdata;
        chomp $seqid;
        my $sequence;
        $$sequence{id} = $seqid;
        push @tblinfo, $sequence;
        
        while ( @seqdata ) {
            shift @seqdata;
            my $start = shift @seqdata;
            my $stop = shift @seqdata;
            my $type = shift @seqdata;
            my $body = shift @seqdata;
            if ( substr( $body, length( $body ) - 1) ne "\n" ) { $body .= "\n" }

            my $feature;
            $$feature{type} = $type;
            $$feature{trunc5} = "";
            if ( $start =~ /^[><]/ ) {
                $$feature{trunc5} = substr( $start, 0, 1 );
                $start = substr( $start, 1 );
            }
            $$feature{start} = $start;
#print "START=\"$$feature{trunc5}\" + \"$$feature{start}\"\n";
            
            $$feature{trunc3} = "";            
            if ( $stop =~ /^[><]/ ) {
                $$feature{trunc3} = substr( $stop, 0, 1 );
                $stop = substr( $stop, 1 );
            }
            $$feature{stop} = $stop;

            $$feature{body} = $body;
            if ( $$feature{type} eq "CDS" ) {
                my $newbody = "";
                my @cdslines = split /\n/, $$feature{body};
                for my $cdsline ( @cdslines ) {
                    if ( $cdsline =~ /^\t\t\tprotein_id\t([^ ]+)$/ ) {
                        $$feature{protein_id} = $1;
                    }
                    elsif ( $cdsline =~ /^\t\t\tnote\t(.+)$/ ) {
                        $$feature{note} = $1;
                    }
                    else {
                        $newbody .= "$cdsline\n";
                    }
                }
                $$feature{body} = $newbody;
            }

            if ( defined $$sequence{features} ) {
                push @{$$sequence{features}}, $feature;
            }
            else {
                my @features = ( $feature );
                $$sequence{features} = \@features;
            } 
        }
    }
    
    return @tblinfo;
}

sub findProtein {
    my ( $tblinfo, $protein_id ) = @_;

#print "find |$protein_id|\n";    
    for my $sequence ( @{ $tblinfo} ) {
#print "$$sequence{id}\n";
        for my $feature ( @{ $$sequence{features} } ) {
#print_hash( "  feature", $feature );
            if ( defined $$feature{protein_id} && $$feature{protein_id} eq $protein_id ) { return $feature }
        }
    }
#print "\n  NOT FOUND\n";
    return undef;
}

sub rewriteTBL {
    my ( $tblinfo, $tblfile ) = @_;
    
    open( TBL, ">$tblfile" ) || die "cannot write to $tblfile\n";
    
    for my $ts ( @$tblinfo ) {
        print TBL ">Features $$ts{id}\n";
        for my $feature ( @{$$ts{features}} ) {
            print TBL "$$feature{trunc5}$$feature{start}\t$$feature{trunc3}$$feature{stop}\t$$feature{type}\n$$feature{body}";
            if ( defined $$feature{protein_id} ) { print TBL "\t\t\tprotein_id\t$$feature{protein_id}\n" }
            if ( defined $$feature{note} ) { print TBL "\t\t\tnote\t$$feature{note}\n" }
        }
    }
    close TBL;
}