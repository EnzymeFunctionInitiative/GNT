#!/usr/bin/env perl

BEGIN {
    die "Please load efishared before runing this script" if not $ENV{EFI_SHARED};
    use lib $ENV{EFI_SHARED};
}


use strict;
use warnings;


use Getopt::Long;
use Data::Dumper;


use EFI::Database;


#my ($uniprotFile, $uniref50File, $uniref90File, $configFile);
my ($uniprotFile, $outFile, $configFile, $uniRefVer, $useV2);
my $result = GetOptions(
    "uniprot-ids=s"             => \$uniprotFile, # input is list of UniProt IDs
    "uniref-mapping=s"          => \$outFile,
    #"uniref50=s"                => \$uniref50File,
    #"uniref90=s"                => \$uniref90File,
    "config=s"                  => \$configFile,
    "uniref-version=i"          => \$uniRefVer, # assume input is list of UniRef IDs, and we get all of the uniprot IDs in the input cluster IDs
    "v2"                        => \$useV2,
);

my $defaultNbSize = 10;
$uniRefVer = 0 if not $uniRefVer or ($uniRefVer != 50 and $uniRefVer != 90);


my $usage = <<USAGE;
usage: $0 --uniprot <input_uniprot_id_file> --uniref50 <output_file> --uniref90 <output_file>
    --config            configuration file to use; if not present looks for EFI_CONFIG env. var
USAGE


die "No --uniprot-ids file provided\n$usage" if not $uniprotFile or not -f $uniprotFile;
die "No --uniref-mapping file provided\n$usage" if not $outFile;
#die "No --uniref50 file provided\n$usage" if not $uniref50File;
#die "No --uniref90 file provided\n$usage" if not $uniref90File;
die "No configuration file found in environment or as argument: \n$usage" if (not $configFile and not $ENV{EFI_CONFIG}) or ($configFile and not -f $configFile) or ($ENV{EFI_CONFIG} and not -f $ENV{EFI_CONFIG});

$configFile = $ENV{EFI_CONFIG} if not $configFile or not -f $configFile;

my %dbArgs;
$dbArgs{config_file_path} = $configFile;

my $mysqlDb = new EFI::Database(%dbArgs);
my $mysqlDbh = $mysqlDb->getHandle();

my @uniprotIds = getIds($uniprotFile);
#my ($ur50Ids, $ur90Ids, $ur50Size, $ur90Size) = getUniRefSeedIds($mysqlDbh, @uniprotIds);
my $map;
if ($useV2) {
    ($map) = getUniRefSeedIds2($mysqlDbh, @uniprotIds);
} else {
    ($map) = getUniRefSeedIds($mysqlDbh, @uniprotIds);
}

writeMapFile($outFile, $map);

#writeFile($uniref50File, $ur50Ids, $ur90Size, $ur50Size);
#writeFile($uniref90File, $ur90Ids, $ur90Size);



sub writeMapFile {
    my $file = shift;
    my $ids = shift;

    open my $fh, ">", $file or die "Unable to open --uniref file $file: $!";

    foreach my $row (@$ids) {
        my @parts = @$row;
        $fh->print(join("\t", @parts) . "\n");
    }

    close $fh;
}


#sub writeFile {
#    my $file = shift;
#    my $ids = shift;
#    my $ur90Size = shift || {};
#    my $ur50Size = shift || {};
#
#    open my $fh, ">", $file or die "Unable to open --uniref file $file: $!";
#
#    foreach my $id (@$ids) {
#        my @parts = ($id);
#        if ($uniref90File or $uniref50File) {
#            push @parts, ($ur90Size->{$id} ? $ur90Size->{$id} : "");
#            push @parts, ($ur50Size->{$id} ? $ur50Size->{$id} : "");
#        }
#        $fh->print(join("\t", @parts) . "\n");
#    }
#
#    close $fh;
#}


sub getIds {
    my $file = shift;

    open my $fh, "<", $file or die "Unable to open --uniprot file $file: $!";

    my @ids;
    while (<$fh>) {
        chomp;
        s/^\s*(.*?)\s*$/$1/;
        next if not $_;
        my @parts = split(m/,/);
        push @ids, @parts;
    }

    close $fh;

    return @ids;
}


sub getUniRefSeedIds {
    my $dbh = shift;
    my @ids = @_;

    my @map;
    my (%ur50, %ur90);

    my $col = $uniRefVer ? "uniref${uniRefVer}_seed" : "accession";
    foreach my $id (@ids) {
        my @p = split(m/\|/, $id);
        my $qid = $id;
        $qid = $p[0] if scalar @p > 1;
        my $sql = "SELECT * FROM uniref WHERE $col = '$qid'";
        my $sth = $dbh->prepare($sql);
        $sth->execute;
        while (my $row = $sth->fetchrow_hashref) {
            $ur50{$row->{uniref50_seed}}++;
            $ur90{$row->{uniref90_seed}}++;
            push @map, [$row->{accession}, $row->{uniref90_seed}, $row->{uniref50_seed}];
        }
    }

    return \@map, \%ur50, \%ur90;
}


sub getUniRefSeedIds2 {
    my $dbh = shift;
    my @ids = @_;
    my %ids;
    map { my @a = split(m/\t/, $_); $ids{$a[$#a]} = 1; } @ids;

    my @map;
    my (%ur50, %ur90);

    my $col = $uniRefVer ? "uniref${uniRefVer}_seed" : "accession";

    my $sql = "SELECT * FROM uniref";
    my $sth = $dbh->prepare($sql);
    $sth->execute;

    my $c = 0;
    while (my $row = $sth->fetchrow_hashref) {
        next if not $ids{$row->{$col}};
        $ur50{$row->{uniref50_seed}}++;
        $ur90{$row->{uniref90_seed}}++;
        push @map, [$row->{accession}, $row->{uniref90_seed}, $row->{uniref50_seed}];
        if ($c++ > 1000) {
            print "$c\r";
            $c = 0;
        }
    }
    print "\n";

    return \@map, \%ur50, \%ur90;
}


#sub getUniRefSeedIds {
#    my $dbh = shift;
#    my @ids = @_;
#
#    my (%ur50, %ur90);
#    my (@ur50, @ur90);
#
#    foreach my $id (@ids) {
#        my $sql = "SELECT * FROM uniref WHERE accession = '$id'";
#        my $sth = $dbh->prepare($sql);
#        $sth->execute;
#        my $row = $sth->fetchrow_hashref;
#        if ($row) {
#            if (not $ur50{$row->{uniref50_seed}}) {
#                push @ur50, $row->{uniref50_seed};
#            }
#            $ur50{$row->{uniref50_seed}}++;
#            if (not $ur90{$row->{uniref90_seed}}) {
#                push @ur90, $row->{uniref90_seed};
#            }
#            $ur90{$row->{uniref90_seed}}++;
#        } else {
#            push @ur50, $id;
#            push @ur90, $id;
#        }
#    }
#
#    return \@ur50, \@ur90, \%ur50, \%ur90;
#}

