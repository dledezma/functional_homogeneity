#################################################################################
###### 22-11-18                                                            ######
###### Execute connectivity calculations (with enzymes)                    ######
###### Calls the creation of GUs, addition of SRs and connectivity         ######
###### calculation for all regulons and their random sets.                 ######
###### Argument usage:                                                     ######
###### [0] - path to connectivity library                                  ######
###### [1] - path to directory with regulon/go/pathways data sets          ######
###### [2] - output dir (subdirs will be created for each data set)        ######
###### Output: one directory per dataset with GU files and a directory     ######
###### with all summaries.                                                 ######
###### History:                                                            ######
###### Notes:                                                              ######
###### - Output directories will be named as data sets.                    ######
###### - add_super_reactions_v4.pl limit is 25 for all data sets           ######
###### - All directory paths in arguments must end in "/".                 ######
###### - Files inside connectivity library must have their original names. ######
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
my $gu_assembly=$library . "gu_assembly_groups_enz_fromfile.pl";
my $add_srs=$library . "add_super_reactions_v5.pl";
my $cal_conn=$library . "connectivity_from_gufiles_permissive_enx_v1.pl";
my $cal_rncon=$library . "conn_of_random_sets_v4.pl";
my $create_random=$library . "random_regulons_same_enzymes.pl";
my $get_enzymes=$library . "all_TRN_enzymes.pl";

#Create dirs for summaries
my $summ_dir=$outdir . "summaries/"; # Dir for real connectivity outputs
system("mkdir $summ_dir");
my $rsum_dir=$outdir . "summaries_random/"; # Dir for random connectivity outputs.
system("mkdir $rsum_dir");

# Extract datasets
opendir(DIR,$data_dir) || die "Cannot open $data_dir.\n";
my @datasets=readdir(DIR);
foreach my $file (sort @datasets){
  if($file=~/^(.+)\.txt/){
    my $dataset=$1;
    
    $time=localtime();
    print"\nDATASET: $dataset\nOn: $time\n";
    my $setpath=$data_dir . $dataset . ".txt"; # path to regulons file

    # Create output dir for dataset
    my $setdir=$outdir . $dataset . "/"; # path to output file
    system("mkdir $setdir");

    # Create inside directories
    my $rawdir=$setdir . "/" . "raw_gus/";
    system("mkdir $rawdir");
    my $comsetdir=$setdir . "/" . "complete_gus/";
    system("mkdir $comsetdir");  
    my $anasetdir=$setdir . "/" . "analysis/";
    system("mkdir $anasetdir");      

    # Call GU assembly
    system("perl $gu_assembly $setpath $rawdir $library");

    # Call addition of complementary pathway reactions
    my $sr_output=$anasetdir . "CPR_count.txt"; # Path to output file
    system("perl $add_srs $library $rawdir $comsetdir $sr_output 25");
    
    # Call Connectivity calculation
    my $summary=$summ_dir . $dataset . ".txt";
    system("perl $cal_conn $comsetdir $summary");

    ### Create random sets with same enzyme number ###
    # Obtain total enzymes per regulon from connectivity summary
    my $total_enz=$anasetdir . "total_enzymes.txt";
    system("grep -v '\#' $summary | cut -f1,2 > $total_enz");

    # Get all enzymes present in dataset as bag for randomization
    my $all_enz=$anasetdir . "all_enzymes_in_set.txt";
    system("perl $get_enzymes $setpath $all_enz");

    # Create directory to save random regulons
    my $random_dir=$setdir . "random_regulons/";
    system("mkdir $random_dir");

    # Create random regulons
    system("perl $create_random $total_enz $random_dir $all_enz 100");

    # Calculate random connectivity
    my $temp=$setdir . "temporal/";
    system("mkdir $temp");	# Creates temporal directory for random connectivity calculation.
    my $fransumm=$rsum_dir . $dataset . "_random.txt";
    system("perl $cal_rncon $random_dir $library $fransumm");
    system("rm -r $temp");	# Removes temporal directory.
  }
}
