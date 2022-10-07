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
use warnings;
use DBI;
use DBD::SQLite;
use Cwd 'realpath';
our $errorMessage;

my %dbSessions;
my %sequences;
my $sybaseServer = "SYBTIGR";
    
sub connectSQLite {
    my ( $dbfile, $autocommit ) = @_;
    if ( !defined $autocommit ) { $autocommit = 0 }
    $errorMessage = undef;
    
    my $dsn = "dbi:SQLite:$dbfile";
       my %session;
       $session{dsn} = $dsn;
       $session{dbfile} = realpath($dbfile);
       $session{username} = "";
       $session{password} = "";
       $session{schema} = "";
    $session{rdbms} = "sqlite";
    $session{autocommit} =    $autocommit;
    $dbSessions{attachments} = undef;
    
    #my $dbh = DBI->connect( $dsn, "", "", {PrintError=>0,RaiseError=>0,AutoCommit=>$autocommit } );
    my $dbh = DBI->connect( $dsn, "", "", {PrintError=>1,RaiseError=>1, sqlite_open_flags => DBD::SQLite::OPEN_READONLY} );
    if ( !$dbh ) {
        $errorMessage = "connectSQLite connect to $dbfile failed: $DBI::errstr";
        die "$errorMessage\n";
    }

    $dbSessions{$dbh} = \%session;
    
    return $dbh;
}

sub querySQLArrayHash {
    my $result = &querySQLinternal(1,@_);
    if ( !defined $result ) {
        $errorMessage = "querySQLArrayHash: " . $errorMessage;
        return undef;
    }
    return $result;
}

sub querySQLinternal {
    my @args = @_;
    my $nargs = @args;

    my $fetchMode = $args[0];
    my $dbh = $args[1];
    my $sql = $args[2];

    my $ibind = 3;
    my $keycol;
    if ( $fetchMode == 2 || $fetchMode==3 ) {
        $keycol = $args[3];
        $ibind++;
    }
    
    my @binds = ();
    while ( $ibind<$nargs ) {
        if ( length($args[$ibind]) ) {
            push(@binds,$args[$ibind]);
        } else {
            push(@binds,undef);
        }
        $ibind++;
    }
    $errorMessage = undef;
    
 # prepare SQL
     my $stmt = $dbh->prepare("$sql");
    if ( !defined $stmt) {
        $errorMessage = "querySQLinternal could not prepare SQL $sql: $DBI::errstr";
        return undef;
    }
    

# execute statement
    if ( !$stmt->execute(@binds) ) {
        $errorMessage = "querySQLinternal could not execute SQL $sql: $DBI::errstr";
        return undef;
    }
    
    
# fetch the rows back from the SELECT statement;

# fetch array of arrays
    if ( $fetchMode==0 ) {
        my @results = ();
        while ( my @row = $stmt->fetchrow ) {
            push(@results,\@row);
        }
        $stmt->finish;
        return \@results;

    } elsif ( $fetchMode==1 ) {
        my @results = ();
        while ( my $row = $stmt->fetchrow_hashref() ) {
            push(@results,$row);
        }
        $stmt->finish;
        return \@results;
        
    } elsif ( $fetchMode==2 ) {
        my %results;
        while ( my $row = $stmt->fetchrow_hashref() ) {
            my $hashkey = $$row{$keycol};
            $results{$hashkey} = $row;
        }
        $stmt->finish;
        if ( !defined \%results ) { return 0 }
        return \%results;
        
    } elsif ( $fetchMode==3 ) {
        my %results;
        while ( my @row = $stmt->fetchrow ) {
            my $hashkey = $row[$keycol];
            $results{$hashkey} = \@row;
        }
        $stmt->finish;
        return \%results;

    } else {
        $errorMessage = "querySQLinternal: unkonw fetchMode: $fetchMode.";
        return undef;
    }
}
1;

