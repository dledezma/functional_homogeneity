#################################################################################
###### 13-08-17                                                            ######
###### Create random regulons                                              ######
###### Creates all_regulon files with same size regulons using random genes #####
###### from original regulon set. Resample IS allowed                      ######
###### [0] - original regulon file                                         ######
###### [1] - outdir                                                        ######
###### [2] - Number of repetitions (bootstraps)                            ######
###### Output: One file per bootstrap, same format as all_regulon files:   ######
###### [1] > GU                                                            ######
###### [2] > genes (gene1//gene2//geneN//)                                 ######
###### History:                                                            ######
###### Notes:                                                              ######
###### - Randomization happens by putting all genes to randomize in an     ######
###### array and then taking random positions from it to recreate regulons ######
###### Taken positions are not deleted.                                    ######
#################################################################################

use strict;
use warnings;
my $time=localtime();

#Get args
my $infile=$ARGV[0];
my $outdir=$ARGV[1];
my $bootstraps=$ARGV[2];
chomp($bootstraps);

my %tfs=();
my %trn=();

# Create array of all participating genes (duplications are eliminated because hash is used to gather genes).
open(IN,$infile) || die "Cannot open IN file at $infile.\n";
while(<IN>){
  if($_=~/^\#/){
    next;
  }
  if($_=~/^([^\t]+)\t([^\t]+)/){
    my $tf=$1;
    my $genes=$2;
    chomp($genes);
    my $gn_count=0;
    my @tf_genes=split(/\/\//,$genes);
    foreach my $gn (@tf_genes){
      $trn{$gn}++;
      $gn_count++;
    }
    $tfs{$tf}=$gn_count;
    print"$tf\t$gn_count\n";
  }
}
close(IN);

# Turn %trn into a vector for ran function.
my @atrn=keys(%trn);


#Loop for number of boostraps
for(my $i=1;$i<=$bootstraps;$i++){
  
  # Create out file
  my $outfile=$outdir . "ran_" . $i . ".txt";  
  open(OUT,">$outfile") || die "Cannot open OUT file at $outfile.\n";
  print OUT "#From $0 on $time using $infile\n";
  
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

%trn=();
@atrn=();
%tfs=();
