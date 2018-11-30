#################################################################################
###### 25-05-17                                                            ######
###### Turn regulons into bnumbers                                         ######
###### Turns get_genes_from_regulons_vX.pl output into bnumbers in same format ##
###### Input:                                                              ######
###### [0] - get_genes_from_regulons_vX.pl output                          ######
###### [1] - gene links file                                               ######
###### [2] - Output file                                                   ######
###### Output: same format                                                 ######
###### History:                                                            ######
#################################################################################

use Sys::Hostname;
$host = hostname;
$time=localtime();
print"Begin: $time\n";

#Get args
$infile=$ARGV[0];
$gene_links=$ARGV[1];
$outfile=$ARGV[2];
chomp($outfile);

#Open OUT file & print header
open(OUT,">$outfile") || die "Cannot open OUT at $outfile.\n";
print OUT"# From $0 on $time @ $host.\n# Input = $infile\n\t$gene_links\n";

open(IN,$infile) || die "Cannot open IN file at $infile.\n";
while(<IN>){
  if($_=~/^([^\t]+)\t([^\t]+)$/){
    $TF=$1;
    $genes=$2;
    chomp($genes);

    print OUT"$TF\t";
    
    @gene=split(/\/\//,$genes);
    foreach $gn (@gene){
      
      $bnum=&IDS($gn);
      print OUT"$bnum//";
    }
    print OUT"\n";
  }
}
close(IN);

###########################
### IDS Function
###
### Turn gene names into bnumbers.
### -Input arguments:
###  [0] - common name
### -Returned values: bnumber.
###########################

sub IDS {
 
  local($name)=($_[0]);
  local($nm);

  open(IDS, $gene_links) || die "Cannot open IDS at $gene_links.\n";
  while(<IDS>){
    if($_=~/^[^\t]*\t[^\t]*\t([^\t]*)\t[^\t]*\t[^\t]*\t[^\t]*\t$name$/i){
      $nm=$1;
      chomp($nm);
      close(IDS);
      return($nm);
    }
  }
  close(IDS);
  print"Warning! $name was not found at $gene_links in &IDS.\n";
  return($name);
}
