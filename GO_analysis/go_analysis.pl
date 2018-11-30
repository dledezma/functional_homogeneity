#################################################################################
###### 23-11-18                                                            ######
###### Execute GO analysis                                                 ######
###### Executes GO analysis without levels 1,2. Creates random             ######
###### regulons keeping gene size and calculates their GO fraction.        ######
###### Argument usage:                                                     ######
###### [0] - path to go_analysis library. Must include:                    ######
######      a) regulons_in_gos_v3.pl                                       ######
######      b) random_regulons.pl                                          ######
######      c) golevels_of_random_sets_v2.pl                               ######
######      d) genes_gos_extended_no12.txt                                 ######
######      e) go_uniq_levels_extended_no12.txt                            ######
###### [1] - path to data sets dir (check regulon format)                  ######
###### [2] - output dir (subdirs will be created)                          ######
###### Output: three directories:                                          ######
######   - random_sets.A folder per data set including 100 randomizations. ######
######   - go_results. A file per data set.                                ######
######   - go_random_results. A file per data set.                         ######
###### History:                                                            ######
###### Notes:                                                              ######
#################################################################################

use strict;
use warnings;
my $time=localtime();

#Get args
my $library=$ARGV[0];
my $data_dir=$ARGV[1];
my $outdir=$ARGV[2];
chomp($outdir);

#Retrieve files from library
my $calculate_gos=$library . "regulons_in_gos_v3.pl";
my $create_random=$library . "random_regulons.pl";
my $cal_random_gos=$library . "golevels_of_random_sets_v3.pl";
my $f_genes_gos=$library . "genes_gos_extended_no12.txt";
my $f_go_levels=$library . "go_uniq_levels_extended_no12.txt";

#Create dirs for summaries
my $summ_dir=$outdir . "go_results/"; # Dir for real connectivity outputs
system("mkdir $summ_dir");
my $rsum_dir=$outdir . "go_random_results/"; # Dir for random connectivity outputs.
system("mkdir $rsum_dir");
my $rlevels_dir=$rsum_dir . "random_levels/"; # Dir for LEVELS random connectivity outputs.
system("mkdir $rlevels_dir");
my $rfraction_dir=$rsum_dir . "random_fractions/"; # Dir for FRACTION random connectivity outputs.
system("mkdir $rfraction_dir");
my $random_dir=$outdir . "random_sets/"; # Dir for randomized data sets
system("mkdir $random_dir");

# Extract datasets
opendir(DIR,$data_dir) || die "Cannot open $data_dir.\n";
my @datasets=readdir(DIR);
foreach my $file (sort @datasets){
  if($file=~/^(.+)\.txt/){
    my $dataset=$1;
    
    $time=localtime();
    print"\nDATASET: $dataset\nOn: $time\n";
    my $setpath=$data_dir . $dataset . ".txt"; # path to regulons file

    # Call calculation on real set
    my $outfile=$summ_dir . $dataset . "_golevels.txt";
    system("perl $calculate_gos $setpath $f_genes_gos $f_go_levels $outfile");

    # Create output dir for randomized datasets
    my $randir=$random_dir . $dataset . "/"; # path to output file
    system("mkdir $randir");

    # Create random data sets
    system("perl $create_random $setpath $randir 100");

    # Call calculations for random regulons
    my $outlevels=$rlevels_dir . $dataset . "_golevels.txt";
    my $outfraction=$rfraction_dir . $dataset . "_gofraction.txt";
    # Create random dir - will be deleted
    my $tempdir=$outdir . "temp/";
    system("mkdir $tempdir");
    system("perl $cal_random_gos $randir $tempdir $library $outlevels $outfraction");
    system("rm -r $tempdir");
  }
}
closedir(DIR);

