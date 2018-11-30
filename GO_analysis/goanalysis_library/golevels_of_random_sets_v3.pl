#################################################################################
###### 18-05-18                                                            ######
###### Obtain GO levels of random gene groups                              ######
###### Retrieves random sets one at a time, creates a proper file for      ######
###### regulons_in_gos.pl and calculates GO levels and percentage of genes ######
###### in GOs for complete set. Saves GO level and perc values in a big table. ##
###### [0] - random regulons dir                                           ######
###### [1] - temporal directory                                            ######     
###### [2] - path to library including the following files:                ######    
######      a) regulons_in_gos_v3.pl                                       ######
######      b) genes_gos_extended_no12.txt                                 ######
######      c) go_uniq_levels_extended_no12.txt                            ######   
###### [3] - outfile for GO levels                                         ######
###### [4] - outfile for higher percentage of genes in GO                  ######
###### Output: matrix of dim num of GUs x num of simulations               ######
###### Notes:                                                              ######
###### - Check which scripts are being used. Re-check versions.            ######
###### History:                                                            ######
###### v2 > 040618 > Omits groups of 1 or less genes to avoid values of    ######
######      1 that are not informative. Returns NA to avoid missing values ######
######      and shifting of columns.                                       ######
###### v3 > 241118 > Retrieves scripts and datasets from library to avoid  ######
######               outdated versions.                                    ######
#################################################################################

use warnings;
use strict;
my $time=localtime();

#Get args
my $dirpath=$ARGV[0];
my $tempdir=$ARGV[1];
my $library=$ARGV[2];
my $outfile_lev=$ARGV[3];
my $outfile_perc=$ARGV[4];
chomp($outfile_perc);

#Retrieve files from library
my $calculate_gos=$library . "regulons_in_gos_v3.pl";
my $f_genes_gos=$library . "genes_gos_extended_no12.txt";
my $f_go_levels=$library . "go_uniq_levels_extended_no12.txt";

#Create path of temporal file for regulons_in_gos to receive
my $temp=$tempdir . "temp_result.txt";

# Create OUT file for GO levels
open(OLV,">$outfile_lev") || die "Cannot open OLV at $outfile_lev.\n";
print OLV"# From $0 on $time\n# Using: $dirpath\n# Percentages of genes of same set are in: $outfile_perc \n# Columns= Group name\n# Rows= Random set\n";

# Create OUT file for gene percentages
open(OPC,">$outfile_perc") || die "Cannot open OPC at $outfile_perc.\n";
print OPC"# From $0 on $time\n# Using: $dirpath\n# GO levels of same set are in: $outfile_lev \n# Columns= Group name\n# Rows= Random set\n";

# Set flag so that first random set prints column names in outfile
my $out_flag=0;
my %group_names=();

# Open directory
opendir(DIR,$dirpath) || die "Cannot open DIR at $dirpath.\n";
my @files=readdir(DIR);
foreach my $set (sort @files){
  if($set=~/^\./){
    next;
  }
  
  print"\n**$set\n";

  # Get random set file path
  my $set_path=$dirpath . $set;
 
  # Get all group names if first iteration
  if(!($out_flag)){
     
    # Open random set file
    open(RNS,$set_path) || die "Cannot open RNS at $set_path.\n";
    while(<RNS>){
      if($_=~/^\#/){
	next;
      }
      # Find an individual GU in set
      if($_=~/^([^\t]+)\t([^\t]+)$/){
	my $tf=$1;
	$group_names{$tf}++;
      }
    }
    close(RNS);

    # Print column names in both OUTFILEs if first random set
    foreach my $tf (sort keys %group_names){
      print OLV"\t$tf";
      print OPC"\t$tf";
    }
    $out_flag=1; 		# Value will change after first random set
    print OLV"\n";
    print OPC"\n";
  }

  # Call regulons_in_gos.pl for the current random set - calling version 2 and up implies that groups of 1 gene or smaller will be omitted. Version 3 implies that genes not present in the GO will be omitted from total fraction calculation.
  print"Calling GO level calculator...  \n$set_path\t$temp\n";
  system("perl $calculate_gos $set_path $f_genes_gos $f_go_levels $temp");

  # Print row name in outfiles
  print OLV"$set";
  print OPC"$set";

  # Parse results for random set. Begins with group names to avoid repeated lines in $temp outfile. Groups of one gene or less will not appear in file and NAs will be printed to avoid shifting of the columns.
  foreach my $gp (sort keys %group_names){
    my $qgp=quotemeta($gp);
    my $found=0;
    open(TMP,$temp) || die "Cannot open TMP at $temp.\n";
    while(<TMP>){
      if($_=~/^$qgp\t[^\t]+\t([^\t]+)\t[^\t]+\t[^\t]+\t([^\t]+)\t/i){
	my $level=$1;
	my $perc=$2;
	print OLV"\t$level";
	print OPC"\t$perc";
	print"$gp\t$perc\n";
	$found=1;
	last;
      }
    }
    close(TMP);

    if(!($found)){		# NAs are printed so that R can ignore those values and does not calculate them as 0. If nothing is printed there will be missing columns.
      print OLV"\tNA";
      print OPC"\tNA";
      print"$gp\tPRINTED NA\n";
    }
  }

  print OLV"\n";
  print OPC"\n";
}

closedir(DIR);

# Delete temporal file
system("rm $temp");

print"From: $time\n";
$time=localtime();
print"To: $time\n";
