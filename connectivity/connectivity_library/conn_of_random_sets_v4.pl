#################################################################################
###### 31-01-18                                                            ######
###### Obtain connectivity of random regulons                              ######
###### Retrieves random sets one at a time, creates a proper file for      ######
###### gu_assembly_for_groups.pl and calculates connectivity for complete  ######
###### set. Saves conn values in a big table.                              ######
###### [0] - random regulons file                                          ######
###### [1] - path to directory containing:                                 ######  
######       - gu_assembly_groups_enz_fromfile.pl                          ######
######       - add_super_reactions_v4.pl                                   ######
######       - connectivity_from_gufiles_permissive_enx_v1.pl              ######  
######       - gu_assembly library (check documentation)                   ######     
###### [2] - outfile                                                       ######
###### Output: matrix of dim num of GUs x num of simulations               ######
###### Notes:                                                              ######
###### - Check which scripts are being used for assembly and connectivity. ######
###### History:                                                            ######
###### 140518 > v2 > Adds super-reactions for a better negative control    ######
###### 250518 > v3 > Prints NA for sets with less than 2 enzymatic rxns    ######
######               for a better negative control.                        ######
###### 221118 > v4 > User provides scripts to be used. Formats were        ######
######               changed to use:                                       ######
######   - GU assembly: gu_assembly_groups_enz_fromfile.pl*                ######
######   - Add CPRs: add_super_reactions_v4.pl                             ######
######   - Calculate conn: connectivity_from_gufiles_permissive_enx_v1.pl* ######
######   * Changed from previous version.                                  ######
#################################################################################

use warnings;
use strict;
my $time=localtime();

#Get args
my $dirpath=$ARGV[0];
my $library=$ARGV[1];
my $outfile=$ARGV[2];
chomp($outfile);

#Retrieve files from library
my $gu_assembly=$library . "gu_assembly_groups_enz_fromfile.pl";
my $add_srs=$library . "add_super_reactions_v4.pl";
my $cal_conn=$library . "connectivity_from_gufiles_permissive_enx_v1.pl";

#Create temporal dir in same dir as outfile. Will be removed at the end.
my @outpath=split(/\//,$outfile);
my $tempdir="";
my $outsize=@outpath;
$outsize--;
for(my $i=0;$i<$outsize;$i++){
  $tempdir=$tempdir . $outpath[$i] . "/";
}
$tempdir=$tempdir . "temporal/";
system("mkdir $tempdir");

# Create OUT file
open(OUT,">$outfile") || die "Cannot open OUT at $outfile.\n";
print OUT"# From $0 on $time\n# Using: $dirpath\n# Columns= GU\n# Rows= Random set\n";

# Set flag so that first random set prints column names in outfile
my $out_flag=0;

# Open directory
opendir(DIR,$dirpath) || die "Cannot open DIR at $dirpath.\n";
my @files=readdir(DIR);
foreach my $set (sort @files){
  if($set=~/^\./){
    next;
  }
  
  # Create temporal directories for gu_assembly output.
  my $gudir=$tempdir . "tempi_raw/";
  my $srdir=$tempdir . "tempi_sr/";
  system("mkdir $gudir");
  system("mkdir $srdir");

  print"\n**$set\n";

  # Get random set file path
  my $set_path=$dirpath . $set;

  # If first random set, open it to print column names.
  if(!($out_flag)){
    # Open random set file
    open(RNS,$set_path) || die "Cannot open RNS at $set_path.\n";
    while(<RNS>){
      if($_=~/^\#/){
	next;
      }
      # Find an individual GU in set
      if($_=~/^([^\t]+)\t[^\t]+$/){
	my $tf=$1;
	print OUT"\t$tf";
      }
    }
    close(RNS);
    print OUT"\n";
    $out_flag=1;
  }
  
  # Call gu_assembly for the current set
  print"Creating raw GU...  \n";
  system("perl $gu_assembly $set_path $gudir $library");

  # Add super-reactions to all GUs in set
  # Create temp file for SR count - will be deleted
  my $srfile=$tempdir . "sr_temp.txt";
  print"Adding SRs... \n ";
  system("perl $add_srs $gudir $srdir $srfile 25");
  system("rm $srfile");

  print"Calculating connectivities for $set \n ";

  # Create second temp file for connectivity
  my $temp=$tempdir . "conn_temp.txt";

  # Calculate connectivity - results will be saved in $temp file
  system("perl $cal_conn $srdir $temp");
  
  # Parse connectivity and append to outfile.
  print OUT"$set"; 		# Print name of random set
  open(CNN,$temp) || die "Cannot open CNN at $temp.\n";
  while(<CNN>){
    if($_=~/^\#/){
      next;
    }
    if($_=~/^[^\t]+\t(\d+)\t[^\t]*\t[^\t]+\t[^\t]*\t[^\t]+\t[^\t]*\t[^\t]+\t[^\t]+\t([^\t]+)$/){
      my $tenz=$1;
      my $conn=$2;
      chomp($conn);
      if($tenz >= 2){	     # Only if random set has two or more enzymes to avoid artificial values of 0.
	print OUT"\t$conn"; 	# Should align with previous column names
      }else{
	print OUT"\tNA";
      }
    }
    print OUT"\n";
  }
  close(CNN);

  # Remove temporal files
  system("rm -r $gudir");
  system("rm -r $srdir");

}
closedir(DIR);

# Delete temporal dir
system("rm -r $tempdir");

print"From: $time\n";
$time=localtime();
print"To: $time\n";
