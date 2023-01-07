#!/usr/bin/env perl

BEGIN {
    die "Please load efishared before runing this script" if not $ENV{EFISHARED};
    use lib $ENV{EFISHARED};
}

use strict;

use FindBin;
use Getopt::Long qw(:config pass_through);
use lib $FindBin::Bin . "/lib";

use EFI::SchedulerApi;
use EFI::Util qw(usesSlurm);
use EFI::GNN::Arrows;


my ($diagramZipFile, $blastSeq, $evalue, $maxNumSeq, $outputFile, $scheduler, $queue, $dryRun,
    $legacy, $title, $nbSize, $idFile, $jobType, $fastaFile, $jobId, $seqDbType, $reverseUniRef);
my ($taxFile, $taxTreeId, $taxIdType, $jobDir, $resultsDirName);
my $result = GetOptions(
    "zip-file=s"            => \$diagramZipFile,

    "blast=s"               => \$blastSeq,
    "evalue=n"              => \$evalue,
    "max-seq=n"             => \$maxNumSeq,
    "nb-size=n"             => \$nbSize, # neighborhood size

    "id-file=s"             => \$idFile,
    "fasta-file=s"          => \$fastaFile,

    "tax-file=s"            => \$taxFile,
    "tax-tree-id=s"         => \$taxTreeId,
    "tax-id-type=s"         => \$taxIdType,

    "output=s"              => \$outputFile,
    "title=s"               => \$title,
    "job-type=s"            => \$jobType,

    "job-id=s"              => \$jobId,
    "scheduler=s"           => \$scheduler,
    "queue=s"               => \$queue,
    "dryrun"                => \$dryRun,
    "legacy"                => \$legacy,
    "seq-db-type=s"         => \$seqDbType,
    "reverse-uniref"        => \$reverseUniRef, # Assume that the input ID list is UniProt IDs, not UniRef
    "job-dir=s"             => \$jobDir,
    "results-dir-name=s"    => \$resultsDirName,
);

my $usage = <<USAGE
usage: $0 --diagram-file <filename> --zip-file <filename> [--zip-file <filename> OR
    --blast <seq_or_file> OR --id-file <filename> OR --fasta-file <filename>]

    # Mode 1
    --zip-file          the zip'ed file to unzip to a SQLite file (simple unzip job submit)

    # Mode 2
    --blast             the sequence for Option A, which uses BLAST to get similar sequences
    --evalue            the evalue to use for BLAST
    --max-seq           the maximum number of sequences to return from the BLAST
    --nb-size           the neighborhood window on either side of the query sequence

    # Mode 3
    --id-file           file containing a list of IDs to use to generate the diagrams

    # Mode 4 (works like Mode 3)
    --fasta-file        file containing FASTA sequences with headers; we extract the IDs from
                        the headers and use those IDs to generate the diagrams

    # Mode 5
    --tax-file          path to the taxonomy json file
    --tax-tree-id       node ID
    --tax-id-type       ID type to use (uniprot|uniref90|uniref50)

    --output            output sqlite file for Options A-D
    --title             the job title to save in the output file
    --job-type          the string to put in for the job type (used by the web app)

    --scheduler         scheduler type (default to torque, but also can be slurm)
    --queue             the cluster queue to use
    --dryrun            if this flag is present, the jobs aren't executed but the job scripts
                        are output to the terminal
    --legacy            if this flag is present, the legacy modules are used
    --seq-db-type       uniprot (default), uniprot-nf, uniref##/uniref##-nf
    --reverse-uniref    if --seq-db-type is uniref##[-nf], then assume input ID list is
                        UniProt; otherwise assume input ID list are UniRef cluster IDs
USAGE
;

my $diagramVersion = $EFI::GNN::Arrows::Version;


if (not -f $diagramZipFile and not $blastSeq and not -f $idFile and not -f $fastaFile and not -f $taxFile) {
    die "$usage";
}

die "The efignt module must be loaded." if not $ENV{EFI_GNN};
die "The efidb module must be loaded." if not $ENV{EFI_DB_MOD};

my $blastMod = $legacy ? "blast" : "BLAST";
if ($blastSeq and $outputFile) {
    if (not $ENV{"BLASTDB"}) {
        die "The $blastMod module must be loaded.";
    } elsif (not $ENV{"EFIDBPATH"}) {
        die "The efidb module must be loaded.";
    }
}


my $toolpath = $ENV{EFI_GNN};
my $efiGnnMod = $ENV{EFI_GNN_MOD};
my $dbMod = $ENV{EFI_DB_MOD};

$jobDir = $ENV{PWD} if not $jobDir;
$resultsDirName = "output" if not $resultsDirName;
my $outputDir = "$jobDir/$resultsDirName";
mkdir $outputDir;

$outputFile = "$outputDir/$outputFile" if $outputFile !~ m%^/%;

$diagramZipFile = "$outputDir/$diagramZipFile"  if $diagramZipFile and $diagramZipFile !~ /^\//;
$queue = "efi"                                  unless $queue =~ /\w/;
$evalue = 5                                     if not $evalue;
$maxNumSeq = 200                                if not $maxNumSeq;
$title = ""                                     if not $title;
$nbSize = 10                                    if not $nbSize;
$jobId = ""                                     if not defined $jobId;


my $stderrFile = "$outputDir/stderr.log";
my $jobCompletedFile = "$outputDir/job.completed";
my $jobErrorFile = "$outputDir/job.error";
my $jobNamePrefix = $jobId ? "${jobId}_" : "";

#my $ur50IdFile = "$outputDir/uniref50.ids";
#my $ur90IdFile = "$outputDir/uniref90.ids";
my $urIdMapFile = "$outputDir/uniref_ids.map";
my $outputUr50File = "$outputFile.uniref50";
my $outputUr90File = "$outputFile.uniref90";


my $schedType = "torque";
$schedType = "slurm" if (defined($scheduler) and $scheduler eq "slurm") or (not defined($scheduler) and usesSlurm());
my $SS = new EFI::SchedulerApi(type => $schedType, queue => $queue, resource => [1, 1, "50GB"], dryrun => $dryRun);


my $titleArg = $title ? "-title \"$title\"" : "";


my $B = $SS->getBuilder();
$B->addAction("rm -f $stderrFile");
$B->addAction("touch $stderrFile");
$B->addAction("module load $efiGnnMod");
$B->addAction("module load $dbMod");

my $jobId;


###################################################################################################
# This job runs a BLAST on the input sequence, then extracts the sequence IDs from the output BLAST
# and then finds all of the neighbors for those IDs and creates the sqlite database from that.
if ($blastSeq) {
    $jobType = "BLAST" if not $jobType;

    my $blastDb = $ENV{EFI_DB_DIR} . "/" . getBlastDbName($seqDbType);
    my $seqFile = "$outputDir/query.fa";
    my $blastOutFile = "$outputDir/blast.raw";
    my $blastIdListFile = "$outputDir/blast.ids";

    open QUERY, "> $seqFile" or die "Unable to open $outputDir/query.fa for writing: $!";
    print QUERY $blastSeq;
    close QUERY;

    $B->resource(1, 1, "70gb");
    $B->addAction("module load $blastMod");
    $B->addAction("blastall -p blastp -i $seqFile -d $blastDb -m 8 -e $evalue -b $maxNumSeq -o $blastOutFile");
    #$B->addAction("grep -v '#' $blastOutFile | cut -f 2,11,12 | sort -k3,3nr | cut -d'|' -f2 > $blastIdListFile");
    #$B->addAction("grep -v '#' $blastOutFile | cut -f 2,11,12 | sort -k3,3nr | sed 's/[\t ]\\{1,\\}/|/g' | cut -d'|' -f2,4 > $blastIdListFile");
    $B->addAction("grep -v '#' $blastOutFile | cut -f 2,11,12 | sort -k3,3nr | sed 's/^[^|]\\+|\\([^|]\\+\\)|[^\t ]\\+\\(.*\\)\$/\\1\\2/' | cut -f1,2 | sed 's/[	]/|/g' > $blastIdListFile");
#    $B->addAction("create_diagram_db.pl -id-file $blastIdListFile -db-file $outputFile -blast-seq-file $seqFile -job-type $jobType $titleArg -nb-size $nbSize");
    outputCreateScript($blastIdListFile, $jobType, "-blast-seq-file $seqFile");
    $B->addAction("echo $diagramVersion > $outputDir/diagram.version");

    addBashErrorCheck($B, 1, $outputFile);
}

elsif ($idFile or ($taxFile and defined $taxTreeId and $taxIdType)) {
    $jobType = "ID_LOOKUP" if not $jobType;

    my $inputFile = "$jobDir/$jobId.txt";
    $B->resource(1, 1, "15gb");
    if ($taxFile and defined $taxTreeId and $taxIdType) {
        $inputFile = "$outputDir/ids_from_tax_tree.txt";
        $B->addAction("extract_taxonomy_tree.pl --json-file $taxFile --output-file $inputFile --id-type $taxIdType --tree-id $taxTreeId");
    } else {
        $B->addAction("cp $idFile $inputFile");
    }
    outputCreateScript($inputFile, $jobType, "--do-id-mapping");
    $B->addAction("echo $diagramVersion > $outputDir/diagram.version");

    addBashErrorCheck($B, 0, $outputFile);
}

elsif ($fastaFile) {
    $jobType = "FASTA" if not $jobType;

    my $inputFile = "$jobDir/$jobId.fasta";
    my $tempIdFile = "$outputDir/temp-ids.txt";

    $B->resource(1, 1, "10gb");

    $B->addAction("cp $fastaFile $inputFile");
    $B->addAction("extract_ids_from_fasta.pl --fasta-file $inputFile --output-file $tempIdFile");
    outputCreateScript($tempIdFile, $jobType, "--do-id-mapping");
    $B->addAction("rm $tempIdFile");
    $B->addAction("echo $diagramVersion > $outputDir/diagram.version");
    
    addBashErrorCheck($B, 0, $outputFile);
}

else {
    $jobType = "unzip";
    $B->resource(1, 1, "5gb");
    ###################################################################################################
    # This job simply unzips the file.
    if ($diagramZipFile =~ m/\.zip$/i) {
        $B->addAction("$toolpath/unzip_file.pl -in $diagramZipFile -out $outputFile -out-ext sqlite 2> $stderrFile");
    } else {
        $B->addAction("cp $diagramZipFile $outputFile");
    }
    $B->addAction("$toolpath/check_diagram_version.pl -db-file $outputFile -version $diagramVersion -version-file $outputDir/diagram.version");
    addBashErrorCheck($B, 1, $outputFile);
}



$jobType = lc $jobType;

my $jobName = "${jobNamePrefix}diagram_$jobType";
my $jobScript = "$outputDir/$jobName.sh";

$B->jobName($jobName);
$B->renderToFile($jobScript);
$jobId = $SS->submit($jobScript);
chomp $jobId;

print "Diagram job ($jobType) is :\n $jobId";



sub outputCreateScript {
    my $idFile = shift;
    my $jobType = shift;
    my $extra = shift || "";

    if ($seqDbType =~ m/uniref([59]0)/) {
        my $uv = $1;
        my $revArg = $reverseUniRef ? "" : "--uniref-version $uv";
        $B->addAction("get_uniref_ids.pl --uniprot-ids $idFile --uniref-mapping $urIdMapFile $revArg");
        $extra .= " --uniref $uv";
        $idFile = $urIdMapFile;
    }
    $B->addAction("create_diagram_db.pl --id-file $idFile $extra --db-file $outputFile --job-type $jobType $titleArg --nb-size $nbSize");
}


sub addBashErrorCheck {
    my ($B, $markAbort, $outputFile) = @_;

    if ($markAbort) {
        $B->addAction("if [ \$? -ne 0 ]; then");
        $B->addAction("    touch $jobErrorFile");
        $B->addAction("fi");
    }
    $B->addAction("if [ ! -f \"$outputFile\" ]; then");
    $B->addAction("    touch $jobErrorFile");
    $B->addAction("fi");
    $B->addAction("touch $jobCompletedFile");

    $B->addAction("");
}


sub getBlastDbName {
    my $seqType = shift;
    my $suffix = $seqType =~ m/\-nf$/ ? "_nf" : "";
    my $name = $seqType =~ m/(uniref[59]0)/ ? "$1$suffix" : "combined$suffix";
    return "$name.fasta";
}


