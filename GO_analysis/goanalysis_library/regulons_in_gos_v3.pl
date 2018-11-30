#################################################################################
###### 15-05-18                                                            ######
###### Lowest GO level for regulon format                                  ######
###### Finds GO that includes the most genes from groups in a file with    ######
###### format:                                                             ######
###### [1] - group name                                                    ######
###### [2] - Genes separated by //                                         ######
###### Argument usage:                                                     ######
###### [0] - File with group of genes                                      ######
###### [1] - gene_gos.txt                                                  ######
###### [2] - GO_levels file (get_uniq_lowest_levels.pl)                    ######
###### [3] - Output file                                                   ######
###### Output: one GO per line                                             ######
###### [1] - group name                                                    ######
###### [2] - GO ID                                                         ######
###### [3] - GO level                                                      ######
###### [4] - Genes in group                                                ######
###### [5] - Genes in GO                                                   ######
###### [6] - Percentage of genes in group in GO                            ######
###### [7] - GO name                                                       ######
###### History:                                                            ######
###### v2 > 040618 > Omits groups with 1 or less genes to omit artificial  ######
######      values of 1.                                                   ######
###### v3 > 060618 > Omits genes that are not present in GO to avoid       ######
######      values of 0 that are not due to a low coverage, but due to     ######
######      incomplete annotation of the genome.                           ######
######      - Adds a column of group size, since it might not be the same  ######
######      than the original, because of genes no in the GO.              ######
######      = Modified regexp to read go_uniq_levels.txt to cpmply with    ######
######      get_uniq_lowest_levels_v2.pl output (total genes column added) ######
###### Notes:                                                              ######
###### - Be aware of types of GOs considered (BP, CC, MF)                  ######
#################################################################################

use strict;
use warnings;
my $time=localtime();

#Get args
my $infile=$ARGV[0];
my $fgogns=$ARGV[1];
my $gotable=$ARGV[2];
my $outfile=$ARGV[3];
chomp($outfile);

#Create and print header of OUT file
open(OUT,">$outfile") || die "Cannot open OUT file at $outfile\n";
print OUT"# From $0 on $time\n# Input = $infile\n#$fgogns\n# $gotable\nGroup name\tGO ID\tGO level\tGenes in group\tGenes in GO\tPercentage\tGO name\n";

# Open gene groups file
open(IN, $infile) || die "Cannot open IN at $infile.\n";
while(<IN>){
  if($_=~/^\#/){
    next;
  }
  if($_=~/^([^\t]+)\t([^\t]+)$/){
    my $name=$1;
    my $genes=$2;
    chomp($genes);

    print"*$name\n";

    chop($genes); 		# Eliminate diagonals at EoL.
    chop($genes);

    my @allgns=split(/\/\//,$genes);
    my @gns=();			# Saves the genes that are in GO.

    # Find genes that are in GO BP and remove those that are not.
    foreach my $gn (@allgns){
      my $gflag=0;
      open(GGO,$fgogns) || die "Cannot open GGO at $fgogns.\n";
      while(<GGO>){
	if($_=~/^$gn\tGO\:\d+\tP/){
	  $gflag=1;
	  last;
	}
      }
      close(GGO);
      if($gflag){
	push(@gns,$gn);
      }
    }

    my $grsize=@gns;

    # Specific to v2 and up. Skips groups of 1 gene to omit uninformative values of 1.
    if($grsize <= 1){
      next;
    }

    my %high=();			# Saves the gos with highest overlap. A hash because ties are possible.
    my $highest=0;		# To avoid going through the complete hash %lowest for comparisons of overlap.
    
    open(GOS,$gotable) || die "Cannot open GOs at $gotable.\n";
    while(<GOS>){
      if($_=~/^([^\t]+)\t(\d+)\t\d+\t([^\t]+)\t([^\t]+)/){
	my $goid=$1;
	my $level=$2;
	my $gogns=$3;

	chop($gogns); 		# Eliminate last comma in line.
	
	my @gons=split(/,/,$gogns);
	
	my $cont=0; 		# Counts the genes in group that are in GO
	foreach my $gn (@gns){
	  foreach my $ggn (@gons){

	    if($gn eq $ggn){
	      #print"$goid - $gn eq $ggn\n";
	      $cont++;
	    }
	  }
	}
	
	if(($cont == 0) && ($highest == 0)){ # If no overlap was found, omit following steps to avoid $cont==$highest
	  next;
	}

	if($cont > $highest){	# If overlap is higher than highest, reset hash and save new highest.
	  %high=();
	  $high{$goid}=$level;
	  $highest=$cont;
	  #print"$goid\t$level\t$cont\n";
	}elsif($cont == $highest){ # If overlap is the same as the highest, add go to %highest.
	  $high{$goid}=$level;
	  #print"$goid\t$level\t$cont\n";
	}else{
	  next;
	}
      }
    }
    close(GOS);

    my $branch=0;		# Saves the lowest branch in the tree (highest GO level) present in hash.
    # Only find highest GO level in hash. Print happens next to permit ties.
    foreach my $ky (sort keys %high){
      if($high{$ky} > $branch){
	$branch=$high{$ky};
      }
    }
  
    if($branch == 0){
      print OUT"$name\tNO GOs\t0.0\t$grsize\tNA\t0.0\tNA\n";
      next;
    }

    ### Print all GOs in has with highest level.
    my $perc=sprintf "%.2f", ($highest/$grsize);
    foreach my $ky (sort keys %high){
      if($high{$ky} == $branch){
	print OUT"$name\t$ky\t$high{$ky}\t$grsize\t$highest\t$perc\t";
	# Find GO name
	open(GOS,$gotable) || die "Cannot open GOs at $gotable.\n";
	while(<GOS>){
	  if($_=~/^$ky\t\d+\t\d+\t[^\t]+\t([^\t]+)/){
	    my $godef=$1;
	    chomp($godef);
	    print OUT"$godef\n";
	    last;
	  }
	}
	close(GOS);
      }
    }
  }
}
close(IN);
close(OUT);





      
