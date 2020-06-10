#!/usr/bin/env perl

BEGIN {
    die "Please load efishared before runing this script" if not $ENV{EFISHARED};
    use lib $ENV{EFISHARED};
}


use strict;
use warnings;


use Getopt::Long;


use EFI::Database;


my ($uniprotFile, $uniref50File, $uniref90File, $configFile);
my $result = GetOptions(
    "uniprot=s"             => \$uniprotFile,
    "uniref50=s"            => \$uniref50File,
    "uniref90=s"            => \$uniref90File,
    "config=s"              => \$configFile,
);

my $defaultNbSize = 10;

my $usage = <<USAGE;
usage: $0 --uniprot <input_uniprot_id_file> --uniref50 <output_file> --uniref90 <output_file>
    --config            configuration file to use; if not present looks for EFI_CONFIG env. var
USAGE


die "No --uniprot file provided\n$usage" if not $uniprotFile or not -f $uniprotFile;
die "No --uniref50 file provided\n$usage" if not $uniref50File;
die "No --uniref90 file provided\n$usage" if not $uniref90File;
die "No configuration file found in environment or as argument: \n$usage" if (not $configFile or -f $configFile) and not exists $ENV{EFI_CONFIG} and not -f $ENV{EFI_CONFIG};

$configFile = $ENV{EFI_CONFIG} if not $configFile or not -f $configFile;

my %dbArgs;
$dbArgs{config_file_path} = $configFile;

my $mysqlDb = new EFI::Database(%dbArgs);
my $mysqlDbh = $mysqlDb->getHandle();

my @uniprotIds = getIds($uniprotFile);
my ($ur50Ids, $ur90Ids) = getUniRefSeedIds($mysqlDbh, @uniprotIds);

writeFile($uniref50File, $ur50Ids);
writeFile($uniref90File, $ur90Ids);




sub writeFile {
    my $file = shift;
    my $ids = shift;

    open my $fh, ">", $file or die "Unable to open --uniref file $file: $!";

    foreach my $id (@$ids) {
        $fh->print("$id\n");
    }

    close $fh;
}


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

    my (%ur50, %ur90);
    my (@ur50, @ur90);

    foreach my $id (@ids) {
        my $sql = "SELECT * FROM uniref WHERE accession = '$id'";
        my $sth = $dbh->prepare($sql);
        $sth->execute;
        my $row = $sth->fetchrow_hashref;
        if ($row) {
            if (not $ur50{$row->{uniref50_seed}}) {
                $ur50{$row->{uniref50_seed}} = 1;
                push @ur50, $id;
            }
            if (not $ur90{$row->{uniref90_seed}}) {
                $ur90{$row->{uniref90_seed}} = 1;
                push @ur90, $id;
            }
        } else {
            push @ur50, $id;
            push @ur90, $id;
        }
    }

    return \@ur50, \@ur90;
}


