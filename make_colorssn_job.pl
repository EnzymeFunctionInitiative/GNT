#!/usr/bin/env perl

BEGIN {
    die "The efishared and efignt environments must be loaded before running this script" if not exists $ENV{EFISHARED} or not exists $ENV{EFIGNN};
    use lib $ENV{EFISHARED};
}

use strict;
use warnings;

use Getopt::Long qw(:config pass_through);
use FindBin;
use File::Basename;
use lib $FindBin::Bin . "/lib";

use EFI::Database;
use EFI::SchedulerApi;
use EFI::Util qw(usesSlurm checkNetworkType);
use EFI::Config;
use EFI::GNN::Base;
use EFI::HMM::Job;


my ($ssnIn, $nbSize, $ssnOut, $cooc, $resultsDirName, $scheduler, $dryRun, $queue, $jobId, $jobDir);
my ($statsFile, $clusterSizeFile, $clusterNumMapFile, $swissprotClustersDescFile, $swissprotSinglesDescFile);
my ($jobConfigFile, $domainMapFileName, $mapFileName, $extraRam, $cleanup);
my ($optMsaOption, $optAaThreshold, $optAaList, $optMinSeqMsa, $optMaxSeqMsa);
my ($efiRefVer, $efiRefDb, $skipFasta, $convRatioFile);
my $result = GetOptions(
    "ssn-in=s"                  => \$ssnIn,
    "ssn-out=s"                 => \$ssnOut,
    "results-dir-name=s"        => \$resultsDirName, # name of results sub dir
    "job-dir=s"                 => \$jobDir,
    "scheduler=s"               => \$scheduler,
    "dry-run"                   => \$dryRun,
    "queue=s"                   => \$queue,
    "job-id=s"                  => \$jobId,
    "map-file-name=s"           => \$mapFileName,
    "domain-map-file-name=s"    => \$domainMapFileName,
    "stats=s"                   => \$statsFile,
    "conv-ratio=s"              => \$convRatioFile,
    "cluster-sizes=s"           => \$clusterSizeFile,
    "cluster-num-map=s"         => \$clusterNumMapFile,
    "sp-clusters-desc=s"        => \$swissprotClustersDescFile,
    "sp-singletons-desc=s"      => \$swissprotSinglesDescFile,
    "job-config=s"              => \$jobConfigFile,
    "extra-ram:i"               => \$extraRam,
    "cleanup"                   => \$cleanup,
    "opt-msa-option=s"          => \$optMsaOption,
    "opt-aa-threshold=s"        => \$optAaThreshold,
    "opt-aa-list=s"             => \$optAaList,
    "opt-min-seq-msa=s"         => \$optMinSeqMsa,
    "opt-max-seq-msa=s"         => \$optMaxSeqMsa,
    "efiref-ver=i"              => \$efiRefVer,
    "efiref-db=s"               => \$efiRefDb,
    "skip-fasta"                => \$skipFasta,
);

my $usage = <<USAGE
usage: $0 -ssnin <filename>

    -ssn-in             path to file of original ssn network to process
    -ssn-out            output filename (not path) for colorized sequence similarity network
    -out-dir            output directory
    -scheduler          scheduler type (default to torque, but also can be slurm)
    -dry-run            only generate the scripts, don't submit to queue
    -queue              the cluster queue to use
    -map-dir-name       the name of the sub-directory to use to output the list of uniprot IDs for each cluster (one file per cluster)
    -map-file-name      the name of the map file that lists uniprot ID, cluster #, and cluster color
    -stats              file to output tabular statistics to
    -clusters-sizes     file to output cluster sizes to
    -sp-clusters-desc   file to output SwissProt descriptions to for any clusters/nodes that have
                        SwissProt annotations (or any metanodes with children that have
                        SwissProt annotations)
    -sp-singletons-desc file to output SwissProt descriptions to for any singleton node that have
                        SwissProt annotations
    -job-config         file specifying the parameters and files to use as input output for this job
    -extra-ram          use increased amount of memory
    -generate-hmm       generate HMMs for each cluster using the FASTA files that are retreived
                        set to 'fast' to generate normal and fast MUSCLE runs, otherwise use 'normal'

The only required argument is -ssnin, all others have defaults.
USAGE
;

my $dbModule = $ENV{EFIDBMOD};
my $gntPath = $ENV{EFIGNN};
my $gntModule = $ENV{EFIGNNMOD};
my $extraPerl = "$ENV{EFI_PERL_ENV}";

#my $JC = LoadJobConfig($jobConfigFile, getDefaults());

if (not exists $ENV{EFICONFIG}) {
    die "Either the configuration file or the EFICONFIG environment variable must be set\n$usage";
}
my $configFile = $ENV{EFICONFIG};

$ssnIn = "$ENV{PWD}/$ssnIn" if $ssnIn and $ssnIn !~ m/^\//;
if (not $ssnIn or not -s $ssnIn) {
    $ssnIn = "" if not $ssnIn;
    die "-ssnin $ssnIn does not exist or has a zero size\n$usage";
}

die "$usage\nERROR: missing -queue parameter" if not $queue;



$jobDir = $ENV{PWD} if not $jobDir;
$resultsDirName = "output" if not $resultsDirName;
my $outputPath = "$jobDir/$resultsDirName";


my ($fn, $fp, $fx) = fileparse($ssnOut, ".zip", ".xgmml", ".xgmml.zip");
my $ssnOutZip = "$fn.zip";
$ssnOut =~ s/\.zip$/.xgmml/;

my $ssnInZip = $ssnIn;
if ($ssnInZip =~ /\.zip$/i) {
    my ($fn, $fp, $fx) = fileparse($ssnIn);
    my $fname = "$fn.xgmml";
    $ssnIn = "$jobDir/$fname";
}

my ($ssnType, $isDomain) = checkNetworkType($ssnInZip);

if (not $optMsaOption) {
    $optMsaOption = 0;
}
$optAaList = $optAaList ? $optAaList : "";
if ($optMsaOption =~ m/CR/ and $optAaList !~ m/^[A-Z,]+$/) {
    $optMsaOption = 0;
}


$skipFasta = 0 if $optMsaOption;


$mapFileName                = "mapping_table.txt"               if not $mapFileName;
$domainMapFileName          = "domain_mapping_table.txt"        if not $domainMapFileName;
$jobId                      = ""                                if not $jobId;
$statsFile                  = "stats.txt"                       if not $statsFile;
$convRatioFile              = "conv_ratio.txt"                  if not $convRatioFile;
$clusterSizeFile            = "cluster_sizes.txt"               if not $clusterSizeFile;
$clusterNumMapFile          = "cluster_num_map.txt"             if not $clusterNumMapFile;
$swissprotClustersDescFile  = "swissprot_clusters_desc.txt"     if not $swissprotClustersDescFile;
$swissprotSinglesDescFile   = "swissprot_singletons_desc.txt"   if not $swissprotSinglesDescFile;

my $jobNamePrefix = $jobId ? "${jobId}_" : "";

$cleanup = defined($cleanup);


my $clusterDataPath             = "cluster-data";
my $uniprotNodeDataDir          = "$clusterDataPath/uniprot-nodes";
my $uniprotDomainNodeDataDir    = "$clusterDataPath/uniprot-domain-nodes";
my $uniRef50NodeDataDir         = "$clusterDataPath/uniref50-nodes";
my $uniRef50DomainNodeDataDir   = "$clusterDataPath/uniref50-domain-nodes";
my $uniRef90NodeDataDir         = "$clusterDataPath/uniref90-nodes";
my $uniRef90DomainNodeDataDir   = "$clusterDataPath/uniref90-domain-nodes";
my $fastaUniProtDataDir         = "$clusterDataPath/fasta";
my $fastaUniProtDomainDataDir   = "$clusterDataPath/fasta-domain";
my $fastaUniRef90DataDir        = "$clusterDataPath/fasta-uniref90";
my $fastaUniRef90DomainDataDir  = "$clusterDataPath/fasta-uniref90-domain";
my $fastaUniRef50DataDir        = "$clusterDataPath/fasta-uniref50";
my $fastaUniRef50DomainDataDir  = "$clusterDataPath/fasta-uniref50-domain";
my $hmmDataDir                  = "$clusterDataPath/hmm";

my $uniprotIdZip = "$outputPath/UniProt_IDs.zip";
my $uniprotDomainIdZip = "$outputPath/UniProt_Domain_IDs.zip";
my $uniRef50IdZip = "$outputPath/UniRef50_IDs.zip";
my $uniRef50DomainIdZip = "$outputPath/UniRef50_Domain_IDs.zip";
my $uniRef90IdZip = "$outputPath/UniRef90_IDs.zip";
my $uniRef90DomainIdZip = "$outputPath/UniRef90_Domain_IDs.zip";
my $fastaZip = "$outputPath/FASTA.zip";
my $fastaDomainZip = "$outputPath/FASTA_Domain.zip";
my $fastaUniRef90Zip = "$outputPath/FASTA_UniRef90.zip";
my $fastaUniRef90DomainZip = "$outputPath/FASTA_UniRef90_Domain.zip";
my $fastaUniRef50Zip = "$outputPath/FASTA_UniRef50.zip";
my $fastaUniRef50DomainZip = "$outputPath/FASTA_UniRef50_Domain.zip";
my $hmmZip = "$outputPath/HMMs.zip";

# The if statements apply to the mkdir cmd, not the die().
my $mkPath = sub {
   my $dir = "$outputPath/$_[0]";
   if ($dryRun) {
       print "mkdir $dir\n";
   } else {
       mkdir $dir or die "Unable to create output dir $dir: $!" if not -d $dir;
   }
};

mkdir $outputPath or die "Unable to create output directory $outputPath: $!" if not $dryRun and not -d $outputPath;
&$mkPath($clusterDataPath);

&$mkPath($uniprotNodeDataDir);
&$mkPath($uniprotDomainNodeDataDir) if not $ssnType or $ssnType eq "UniProt" and $isDomain;
&$mkPath($fastaUniProtDataDir);
&$mkPath($fastaUniProtDomainDataDir) if not $ssnType or $ssnType eq "UniProt" and $isDomain;

&$mkPath($uniRef50NodeDataDir);
&$mkPath($uniRef50DomainNodeDataDir) if not $ssnType or $ssnType eq "UniRef50" and $isDomain;
&$mkPath($fastaUniRef50DataDir);
&$mkPath($fastaUniRef50DomainDataDir) if not $ssnType or $ssnType eq "UniRef50" and $isDomain;

&$mkPath($uniRef90NodeDataDir);
&$mkPath($uniRef90DomainNodeDataDir) if not $ssnType or $ssnType eq "UniRef90" and $isDomain;
&$mkPath($fastaUniRef90DataDir);
&$mkPath($fastaUniRef90DomainDataDir) if not $ssnType or $ssnType eq "UniRef90" and $isDomain;

&$mkPath($hmmDataDir) if $optMsaOption;


my $fileSize = 0;
if ($ssnInZip !~ m/\.zip/) { # If it's a .zip we can't predict apriori what the size will be.
    $fileSize = -s $ssnIn;
}

# Y = MX+B, M=emperically determined, B = safety factor; X = file size in MB; Y = RAM reservation in GB
my $ramReservation = defined $extraRam ? (length $extraRam ? $extraRam : 800) : 150;
if (!$extraRam && $fileSize) {
    my $ramPredictionM = 0.02;
    my $ramSafety = 10;
    $fileSize = $fileSize / 1024 / 1024; # MB
    $ramReservation = $ramPredictionM * $fileSize + $ramSafety;
    $ramReservation = int($ramReservation + 0.5);
    $ramReservation = 1000 if $ramReservation > 1000;
}

my $schedType = "torque";
$schedType = "slurm" if (defined($scheduler) and $scheduler eq "slurm") or (not defined($scheduler) and usesSlurm());
my $SS = new EFI::SchedulerApi(type => $schedType, queue => $queue, resource => [1, 1], dryrun => $dryRun);

my $absPath = sub {
    return $_[0] =~ m/^\// ? $_[0] : "$outputPath/$_[0]";
};


my $fileInfo = {
    color_only => 1,
    config_file => $configFile,
    tool_path => $gntPath,
    fasta_tool_path => "$gntPath/get_fasta.pl",
    cat_tool_path => "$gntPath/cat_files.pl",
    ssn_out => "$outputPath/$ssnOut",
    ssn_out_zip => "$outputPath/$ssnOutZip",

    domain_map_file => "$domainMapFileName",
    map_file => "$mapFileName",

    input_seqs_file => "ssn-sequences.fa",
};

$fileInfo->{uniprot} = {
    id_dir => &$absPath($uniprotNodeDataDir),
    fasta_dir => &$absPath($fastaUniProtDataDir),
    id_zip => $uniprotIdZip,
    fasta_zip => $fastaZip,
};


# The 'not $ssnType or' statement ensures that this happens.
if (not $ssnType or $ssnType eq "UniProt" and $isDomain) {
    $fileInfo->{uniprot_domain} = {
        id_dir => &$absPath($uniprotDomainNodeDataDir),
        fasta_dir => &$absPath($fastaUniProtDomainDataDir),
        id_zip => $uniprotDomainIdZip,
        fasta_zip => $fastaDomainZip,
    };
}

if (not $ssnType or $ssnType eq "UniRef90" or $ssnType eq "UniRef50") {
    $fileInfo->{uniref90} = {
        id_dir => &$absPath($uniRef90NodeDataDir),
        fasta_dir => &$absPath($fastaUniRef90DataDir),
        id_zip => $uniRef90IdZip,
        fasta_zip => $fastaUniRef90Zip,
    };
    if (not $ssnType or $isDomain and $ssnType eq "UniRef90") {
        $fileInfo->{uniref90_domain} = {
            id_dir => &$absPath($uniRef90DomainNodeDataDir),
            fasta_dir => &$absPath($fastaUniRef90DomainDataDir),
            id_zip => $uniRef90DomainIdZip,
            fasta_zip => $fastaUniRef90DomainZip,
        };
    }
}

if (not $ssnType or $ssnType eq "UniRef50") {
    $fileInfo->{uniref50} = {
        id_dir => &$absPath($uniRef50NodeDataDir),
        fasta_dir => &$absPath($fastaUniRef50DataDir),
        id_zip => $uniRef50IdZip,
        fasta_zip => $fastaUniRef50Zip,
    };
    if (not $ssnType or $isDomain) {
        $fileInfo->{uniref50_domain} = {
            id_dir => &$absPath($uniRef50DomainNodeDataDir),
            fasta_dir => &$absPath($fastaUniRef50DomainDataDir),
            id_zip => $uniRef50DomainIdZip,
            fasta_zip => $fastaUniRef50DomainZip,
        };
    }
}

if ($optMsaOption) {
    $fileInfo->{hmm_tool_path} = "$gntPath/build_hmm.pl";
    $fileInfo->{hmm_tool_dir} = "$gntPath/hmm";
    $fileInfo->{hmm_data_dir} = &$absPath($hmmDataDir);
    $fileInfo->{hmm_zip} = $hmmZip;
    $fileInfo->{hmm_logo_list} = "$outputPath/hmm_logos.txt";
    $fileInfo->{hmm_weblogo_list} = "$outputPath/weblogos.txt";
    $fileInfo->{hmm_histogram_list} = "$outputPath/histograms.txt";
    $fileInfo->{hmm_alignment_list} = "$outputPath/alignments.txt";
    $fileInfo->{hmm_consensus_residue_info_list} = "$outputPath/consensus_residue.txt";
    $fileInfo->{hmm_rel_path} = $hmmDataDir;
    $fileInfo->{hmm_count_aa_tool_path} = "$gntPath/count_aa.pl";
    $fileInfo->{hmm_collect_id_tool_path} = "$gntPath/collect_aa_hmm_ids.pl";

    $optAaThreshold = ($optAaThreshold and $optAaThreshold =~ m/^[0-9,\.]+$/) ? $optAaThreshold : 0;
    $fileInfo->{hmm_consensus_threshold} = $optAaThreshold;
    $fileInfo->{hmm_option} = $optMsaOption;
    $fileInfo->{hmm_amino_acids} = [split(m/,/, $optAaList)];
    my @colors = ("red", "blue", "orange", "DarkGreen", "Magenta", "Gray");
    $fileInfo->{hmm_weblogo_colors} = \@colors;
    $fileInfo->{hmm_max_seq_msa} = $optMaxSeqMsa // 700;

    $optMinSeqMsa = ($optMinSeqMsa and $optMinSeqMsa >= 1) ? $optMinSeqMsa : 5;
    $fileInfo->{hmm_min_seq_msa} = $optMinSeqMsa;

    $fileInfo->{output_path} = $outputPath;
    $fileInfo->{cluster_size_file} = $clusterSizeFile;
    $fileInfo->{ssn_type} = $ssnType;
    $fileInfo->{hmm_zip_prefix} = ""; #${ssnName}";

    $fileInfo->{compute_pim} = 1;

    $fileInfo->{gnt_module} = $gntModule;
    $fileInfo->{weblogo_path} = "$ENV{EFI_GROUP_HOME}/bin/weblogo"; 
    $fileInfo->{weblogo_lib} = "$ENV{EFI_GROUP_HOME}/lib/python3.7"; #TODO: get this from the environment
    $fileInfo->{extra_perl} = $extraPerl;
}


my $scriptArgs = 
    " -output-dir $outputPath" .
    " -ssnin $ssnIn" .
    " -ssnout $ssnOut" .
    " -uniprot-id-dir $uniprotNodeDataDir" .
    " -uniprot-domain-id-dir $uniprotDomainNodeDataDir" .
    " -uniref50-id-dir $uniRef50NodeDataDir" .
    " -uniref50-domain-id-dir $uniRef50DomainNodeDataDir" .
    " -uniref90-id-dir $uniRef90NodeDataDir" .
    " -uniref90-domain-id-dir $uniRef90DomainNodeDataDir" .
    " -id-out $mapFileName" .
    " -id-out-domain $domainMapFileName" .
    " -config $configFile" .
    " -stats \"$statsFile\"" .
    " -conv-ratio \"$convRatioFile\"" .
    " -cluster-sizes \"$clusterSizeFile\"" .
    " -cluster-num-map \"$clusterNumMapFile\"" .
    " -sp-clusters-desc \"$swissprotClustersDescFile\"" .
    " -sp-singletons-desc \"$swissprotSinglesDescFile\"" .
    ""
    ;
if ($efiRefVer) {
    $scriptArgs .= " --efiref-ver $efiRefVer --efiref-db $efiRefDb";
    $scriptArgs .= " --efiref-50-dir $fileInfo->{efiref50}->{id_dir}" if $efiRefVer == 50;
    $scriptArgs .= " --efiref-70-dir $fileInfo->{efiref70}->{id_dir}" if $efiRefVer <= 70;
}

my $B = $SS->getBuilder();

$B->resource(1, 1, "${ramReservation}gb");
$B->addAction("module load $dbModule");
$B->addAction("module load $gntModule");
$B->addAction("source $extraPerl");
$B->addAction("cd $outputPath");
$B->addAction("$gntPath/unzip_file.pl -in $ssnInZip -out $ssnIn") if $ssnInZip =~ /\.zip/i;
#TODO: remove this hack
#$B->addAction("/home/n-z/noberg/dev/EST/fix_xgmml.pl --input $ssnIn --rename");
$B->addAction("$gntPath/cluster_gnn.pl $scriptArgs");
EFI::GNN::Base::addFileActions($B, $fileInfo, $skipFasta);
if (not $optMsaOption) {
    $B->addAction("touch $outputPath/1.out.completed");
}

my $jobName = "${jobNamePrefix}color_ssn";
my $jobScript = "$outputPath/$jobName.sh";
$B->jobName($jobName);
$B->renderToFile($jobScript);
$jobId = $SS->submit($jobScript);
print "Color SSN job is:\n $jobId";

if ($optMsaOption) {
    $B = EFI::HMM::Job::makeJob($SS, $fileInfo, $jobId);
    $B->addAction("touch $outputPath/1.out.completed");
    $jobName = "${jobNamePrefix}hmm_and_stuff";
    $jobScript = "$outputPath/$jobName.sh";
    $B->jobName($jobName);
    $B->renderToFile($jobScript);
    $jobId = $SS->submit($jobScript);
    print "\nHMM and stuff job is:\n $jobId";
}

if ($cleanup) {
    $B = $SS->getBuilder(); 
    $B->resource(1, 1, "1GB");
    $B->dependency(0, $jobId);
    $B->addAction("cd $outputPath");
    my $paths = join(" ", map { "$outputPath/$_" } (
        $uniprotNodeDataDir, $uniprotDomainNodeDataDir, $fastaUniProtDataDir, $fastaUniProtDomainDataDir,
        $uniRef50NodeDataDir, $uniRef50DomainNodeDataDir, $fastaUniRef50DataDir, $fastaUniRef50DomainDataDir,
        $uniRef90NodeDataDir, $uniRef90DomainNodeDataDir, $fastaUniRef90DataDir, $fastaUniRef90DomainDataDir
    ));
    $B->addAction("rm -rf $paths");
    $jobName = "${jobNamePrefix}cleanup";
    $jobScript = "$outputPath/$jobName.sh";
    $B->jobName($jobName);
    $B->renderToFile($jobScript);
    $jobId = $SS->submit($jobScript);
}


sub getDefaults {
    my $defaults = {
        mapDirName                  => "cluster-data",
        mapFileName                 => "mapping_table.txt",
        jobId                       => "",
        statsFile                   => "stats.txt",
        clusterSizeFile             => "cluster_sizes.txt",
        swissprotClustersDescFile   => "swissprot_clusters_desc.txt",
        swissprotSinglesDescFile    => "swissprot_singletons_desc.txt",
    };
    return $defaults;
}


