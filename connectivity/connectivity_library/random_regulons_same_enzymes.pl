#################################################################################
###### 15-02-18                                                            ######
###### Create random regulons with same number of enzymes                  ######
###### Creates all_regulon files for groups of genes mantaining the number ######
###### of enzymes in each group. Takes random enzymes from a list of       ######
###### enzymes (file provided by user). Resampling IS allowed.             ######
###### [0] - File with name of GU and number of enzymes                    ######
###### [1] - outdir                                                        ######
###### [2] - Dataset of all enzymes (all_TRN_enzymes.pl output)            ######
###### [3] - Number of repetitions (bootstraps)                            ######
###### Output: One file per bootstrap, same format as all_regulon files:   ######
###### [1] > GU                                                            ######
###### [2] > genes (gene1//gene2//geneN//)                                 ######
###### History:                                                            ######
###### Notes:                                                              ######
###### - Originally written for extended GUs, same code, different         ######
######   name and description for easier interpretation of pipelines.      ######
#################################################################################

use strict;
use warnings;
my $time=localtime();

#Get args
my $infile=$ARGV[0];
my $outdir=$ARGV[1];
my $enzymes=$ARGV[2];
my $bootstraps=$ARGV[3];
chomp($bootstraps);

my %tfs=();
my %enzs=();

# Create array of all participating genes (genes that are enzymes only).
open(ENZ,$enzymes) || die "Cannot open ENZ file at $enzymes.\n";
while(<ENZ>){
  if($_=~/^\#/){
    next;
  }
  if($_=~/^[^\t]+\t([^\t]+)\t/){
    my $gen=$1;
    chomp($gen);
    $enzs{$gen}++;
  }
}
close(ENZ);

# Turn %enzs into a vector for ran function.
my @atrn=keys(%enzs);

# Retrieve number of enzymes per GU
open(IN,$infile) || die "Cannot open IN at $infile.\n";
while(<IN>){
  if($_=~/^\#/){
    next;
  }
  if($_=~/^([^\t]+)\t([^\t]+)/){
    my $tf=$1;
    my $tenz=$2;
    chomp($tenz);

    $tfs{$tf}=$tenz;
    print"$tf\t$tenz\n";
  }
}
close(IN);

#Loop for number of boostraps
for(my $i=1;$i<=$bootstraps;$i++){
  
  # Create out file
  my $outfile=$outdir . "rn_cluster_" . $i . ".txt";  
  open(OUT,">$outfile") || die "Cannot open OUT file at $outfile.\n";
  print OUT "# From $0 on $time using: $infile\n# $enzymes\n";
  
  foreach my $tf (sort keys %tfs){
    print OUT"$tf\t";
    for(my $g=0;$g<$tfs{$tf};$g++){
      my $ran_gn=$atrn[rand @atrn];
      print OUT"$ran_gn//";
    }
    print OUT"\n";
  }
  close(OUT);
}
