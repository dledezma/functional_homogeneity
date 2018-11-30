#################################################################################
###### 07-05-18                                                            ######
###### Create data sets                                                    ######
###### Creates 5 data sets of gene groups for GU analysis:                 ######
###### a) General regulons (all TFs).                                      ######
###### b) Strict regulons (all TFs, two groups for dual TFs)               ######
###### c) Regulons by TF conformation (only Tfs with known conformation)   ######
###### d) Regulons by TF and effect (only Tfs with known conformation)     ######
###### e) Complex regulons (all TFs)                                       ######
###### f) Simple regulons (all TFs)                                        ######
###### Argument usage:                                                     ######
###### [0] - RegulonDB TF-gene data set                                    ######
###### [1] - TF_active                                                     ######
###### [2] - TF_inactive                                                   ######
###### [3] - RegulonDB TU-gene dataset                                     ######
###### [4] - Output dir                                                    ######
###### Output: 5 files will be created in output dir. Format is:           ######
###### [1] - GU name                                                       ######
###### [2] - Genes separated by //                                         ######
###### History:                                                            ######
#################################################################################

use strict;
use warnings;
my $time=localtime();

#Get args
my $regdb=$ARGV[0];
my $active=$ARGV[1];
my $inactive=$ARGV[2];
my $tulink=$ARGV[3];
my $outdir=$ARGV[4];
chomp($outdir);

# Create outfiles
my $out_a=$outdir . "general_regulons.txt";
my $out_b=$outdir . "strict_regulons.txt";
my $out_c=$outdir . "regulons_by_conformation.txt";
my $out_d=$outdir . "regulons_by_conformation_effect.txt";
my $out_e=$outdir . "complex_regulons.txt";
my $out_f=$outdir . "simple_regulons.txt";


###### (a) GENERAL REGULONS #####
# All genes that have a binding site for a TF, regardless of effect, conformation, num. of binding sites, etc.
#####

print"*Creating data set (a) - General regulons\n";

#Create and print header of OUT file
open(OUT,">$out_a") || die "Cannot open OUT file at $out_a\n";
print OUT"# From $0 on $time\n# Input = $regdb\n# $active\n# $inactive\n#GU name\tGenes\n";

# Get groups from data set
my %trn=();
open(RDB,$regdb) || die "Cannot open RDB at $regdb for dataset A.\n";
while(<RDB>){
  if($_=~/^([^\t]+)\t([^\t]+)\t([^\t]+)/){
    my $tf=$1;
    my $gene=$2;
    $trn{$tf}{$gene}++;
  }
}
close(RDB);

# Print output
foreach my $ky (sort keys %trn){
  print OUT"$ky\t";
  foreach my $gn (sort keys %{$trn{$ky}}){
    print OUT"$gn//";
  }
  print OUT"\n";
}

close(OUT);

###### (b) STRICT REGULONS #####
# All genes that have a binding site for a TF, with a defined effect. Duals are considered in both groups (TF+ & TF-). Omits interactions with unknown effect (?).
#####

print"*Creating data set (b) - Strict regulons\n";

#Create and print header of OUT file
open(OUT,">$out_b") || die "Cannot open OUT file at $out_b\n";
print OUT"# From $0 on $time\n# Input = $regdb\n# $active\n# $inactive\n#GU name\tGenes\n";

# Get groups from data set
%trn=();
open(RDB,$regdb) || die "Cannot open RDB at $regdb for dataset B.\n";
while(<RDB>){
  if($_=~/^([^\t]+)\t([^\t]+)\t([^\t]+)/){
    my $tf=$1;
    my $gene=$2;
    my $effect=$3;
    if($effect=~/^\+$/){
      $tf=$tf . "_pos";
      $trn{$tf}{$gene}++;
    }elsif($effect=~/^\-$/){
      $tf=$tf . "_neg";
      $trn{$tf}{$gene}++;
    }elsif($effect=~/^\+-$/){
      my $tf1=$tf . "_pos";
      $trn{$tf1}{$gene}++;
      $tf=$tf . "_neg";
      $trn{$tf}{$gene}++;
    }
  }
}
close(RDB);

# Print output
foreach my $ky (sort keys %trn){
  print OUT"$ky\t";
  foreach my $gn (sort keys %{$trn{$ky}}){
    print OUT"$gn//";
  }
  print OUT"\n";
}

close(OUT);


###### (c) REGULONS by TF Conformation #####
# All genes that have a binding site for a TF, and the effect (whichever) happens under a specific conformation of the TF. Inactive conformations are omitted. 
#####

print"*Creating data set (C) - REGULONS by TF Conformation\n";

#Create and print header of OUT file
open(OUT,">$out_c") || die "Cannot open OUT file at $out_c\n";
print OUT"# From $0 on $time\n# Input = $regdb\n# $active\n# $inactive\n# $tulink\n#GU name\tGenes\n";

%trn=();
open(ACT,$active) || die "Cannot open ACT at $active for dataset C.\n";
while(<ACT>){
  if($_=~/^[^\t]+\t([^\t]+)\t[^\t]+\t[^\t]+\t([^\t]+)/){
    my $conf=$1;
    $conf=&xml($conf);
    my $tu=$2;

    # Find genes in TU
    open(TUL,$tulink) || die "Cannot open TUL at $tulink for data set (c).\n";
    while(<TUL>){
      if($_=~/^[^\t]+\t$tu\t[^\t]+\t([^\t]+)/i){
	my $genes=$1;
	my @gns=split(/,/,$genes);
	foreach my $gn (@gns){
	  $trn{$conf}{$gn}++;
	}
      }
    }
    close(TUL);
  }
}
close(ACT);


      
# Print output
foreach my $ky (sort keys %trn){
  print OUT"$ky\t";
  foreach my $gn (sort keys %{$trn{$ky}}){
    print OUT"$gn//";
  }
  print OUT"\n";
}

close(OUT);

###### (d) REGULONS by TF CONFORMATION & EFFECT #####
# All genes that have a binding site for a TF, and the RI happens under a specific conformation of the TF with the same effect. Inactive conformations are omitted. 
#####

print"*Creating data set (d) - REGULONS by TF CONFORMATION & EFFECT\n";


#Create and print header of OUT file
open(OUT,">$out_d") || die "Cannot open OUT file at $out_d\n";
print OUT"# From $0 on $time\n# Input = $regdb\n# $active\n# $inactive\n# $tulink\n#GU name\tGenes\n";

%trn=();
open(ACT,$active) || die "Cannot open ACT at $active for dataset D.\n";
while(<ACT>){
  if($_=~/^[^\t]+\t([^\t]+)\t[^\t]+\t([^\t]+)\t([^\t]+)/){
    my $conf=$1;
    $conf=&xml($conf);
    my $effect=$2;
    my $tu=$3;

    if($effect=~/^\+$/){
      $conf=$conf . "_pos";
    }elsif($effect=~/^\-$/){
      $conf=$conf . "_neg";
    }elsif($effect=~/^\+\-$/){
      my $conf2=$conf . "_pos";

      # Add genes for positive effect - genes for negative effect will be done outside of if.
      # Find genes in TU
      open(TUL,$tulink) || die "Cannot open TUL at $tulink for data set D.\n";
      while(<TUL>){
	if($_=~/^[^\t]+\t$tu\t[^\t]+\t([^\t]+)/i){
	  my $genes=$1;
	  my @gns=split(/,/,$genes);
	  foreach my $gn (@gns){
	    $trn{$conf2}{$gn}++;
	  }
	}
      }
      close(TUL);
      $conf=$conf . "_neg"; 	# To find genes with negative effect outside if loop.
    }else{
      next;
    }

    # Find genes in TU
    open(TUL,$tulink) || die "Cannot open TUL at $tulink for data set D.\n";
    while(<TUL>){
      if($_=~/^[^\t]+\t$tu\t[^\t]+\t([^\t]+)/i){
	my $genes=$1;
	my @gns=split(/,/,$genes);
	foreach my $gn (@gns){
	  $trn{$conf}{$gn}++;
	}
      }
    }
    close(TUL);
  }
}
close(ACT);


      
# Print output
foreach my $ky (sort keys %trn){
  print OUT"$ky\t";
  foreach my $gn (sort keys %{$trn{$ky}}){
    print OUT"$gn//";
  }
  print OUT"\n";
}

close(OUT);


###### (e) COMPLEX REGULONS #####
# All genes that have a binding site for a TF, and the RI happens under a specific conformation of the TF with the same effect. Inactive conformations are omitted. 
#####

print"*Creating data set  (e) COMPLEX REGULONS\n";

#Create and print header of OUT file
open(OUT,">$out_e") || die "Cannot open OUT file at $out_e\n";
print OUT"# From $0 on $time\n# Input = $regdb\n# $active\n# $inactive\n#GU name\tGenes\n";

# # Get all combinations of TFs
my %tus=();
my %complex=();
open(RDB,$regdb) || die "Cannot open RDB at $regdb for dataset E.\n";
while(<RDB>){
  if($_=~/^([^\t]+)\t([^\t]+)\t[^\t]+/){
    my $tf=$1;
    my $gene=$2;

    my $flag=0; 		# Signals if TF is already on array - for cases when an RI is duplicated.

    foreach my $tfa (@{$tus{$gene}}){
      if($tfa eq $tf){
	$flag=1;
	last;
      }
    }

    if(!($flag)){		# If TF is not already in array
      push(@{$tus{$gene}},$tf);
    }
  }
}
close(RDB);

# Create new hash with all combinations of TFs that actually have a RI and the genes the combination regulates.
foreach my $gn (keys %tus){
  @{$tus{$gn}}=sort @{$tus{$gn}};
  # Create name of regulon
  my $name=$tus{$gn}[0];
  my $size=@{$tus{$gn}};
  for(my $i=1;$i<$size;$i++){
    $name=$name . "_" . $tus{$gn}[$i];
  }
  # Fill array of genes regulated by combination of TFs
  my $flag=0; 			# Signals if gene is already in array.
  foreach my $rgns (@{$complex{$name}}){
    if($gn eq $rgns){
      $flag=1;
      last;
    }
  }
  
  if(!($flag)){ 		# If gene is not yet in array of complex regulon.
    push(@{$complex{$name}},$gn);
  }
}


# Print output
foreach my $ky (sort keys %complex){
  print OUT"$ky\t";
  foreach my $gen (sort @{$complex{$ky}}){
    print OUT"$gen//";
  }
  print OUT"\n";
}

close(OUT);

###### (f) SIMPLE REGULONS #####
# All genes that have a binding site for ONLY one TF and not other. 
#####

print"*Creating data set (f) SIMPLE REGULONS\n";

#Create and print header of OUT file
open(OUT,">$out_f") || die "Cannot open OUT file at $out_f\n";
print OUT"# From $0 on $time\n# Input = $regdb\n# $active\n# $inactive\n#GU name\tGenes\n";

# # Get all combinations of TFs
%tus=();
my %simple=();
open(RDB,$regdb) || die "Cannot open RDB at $regdb for dataset F.\n";
while(<RDB>){
  if($_=~/^([^\t]+)\t([^\t]+)\t[^\t]+/){
    my $tf=$1;
    my $gene=$2;

    my $flag=0; 		# Signals if TF is already on array - for cases when an RI is duplicated.

    foreach my $tfa (@{$tus{$gene}}){
      if($tfa eq $tf){
	$flag=1;
	last;
      }
    }

    if(!($flag)){		# If TF is not already in array
      push(@{$tus{$gene}},$tf);
    }
  }
}
close(RDB);

# Create new hash with all combinations of TFs that actually have a RI and the genes the combination regulates.
foreach my $gn (keys %tus){
  # Filter all genes that are regulated by more than one TF.
  my $size=@{$tus{$gn}};
  if($size > 1){
    next;
  }
  push(@{$simple{$tus{$gn}[0]}},$gn);
}


# Print output
foreach my $ky (sort keys %simple){
  print OUT"$ky\t";
  foreach my $gn (sort @{$simple{$ky}}){
    print OUT"$gn//";
  }
  print OUT"\n";
}

close(OUT);


###########################
### XML Function
###
### Remove XML code from molecule names
### -Input arguments:
###  [0] - line
### -Returned values: line in a single char value
###########################

sub xml{

  my($line)=$_[0];
  
  chomp($line);

  while($line=~/(.*)\<.+\>(.*)/){
    $line=$1 . $2;
  }
  while($line=~/(.*)&alpha;(.*)/i){
    $line=$1 . "alpha" . $2;
  }
  while($line=~/(.*)&beta;(.*)/i){
    $line=$1 . "beta" . $2;
  }
  while($line=~/(.*)&gamma;(.*)/i){
    $line=$1 . "gamma" . $2;
  }
  while($line=~/(.*)&delta;(.*)/i){
    $line=$1 . "delta" . $2;
  }
  while($line=~/(.*)&sigma;(.*)/i){
    $line=$1 . "sigma" . $2;
  }
  while($line=~/(.*)&omega;(.*)/i){
    $line=$1 . "omega" . $2;
  }
  return($line);
  
}
