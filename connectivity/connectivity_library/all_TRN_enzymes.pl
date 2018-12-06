#################################################################################
###### 19-02-18                                                            ######
###### Get all enzymes from Ecocyc that are in TRN                         ######
###### Retrieves all enzymes from Ecocyc and only prints those present     ######    
###### in the TRN                                                          ######     
###### [0] - all regulons file                                             ######
###### [1] - outfile                                                       ######
###### Output:                                                             ######
###### [f1] > Enzyme gene Ecocyc ID                                        ######
###### [f2] > Enzyme gene name                                             ######
###### [f3] > Reactions Ecocyc ID                                          ######
###### Notes:                                                              ######
###### - Calls all reactions, then enzymes and then genes.                 ######
###### History:                                                            ######
#################################################################################

use warnings;
use perlcyc;
use strict;			
my $time=localtime();
my $cyc = perlcyc -> new("ECOLI");
my $temp="";
my %genenz=();

###Get args
my $regulons=$ARGV[0];
my $outfile=$ARGV[1];
chomp($outfile);


#Print header on outfile
open(OUT,">$outfile") || die "Cannot open OUT at $outfile.\n";
print OUT"# From $0 on $time\n# f1 > Enzyme Ecocyc ID\n# f2 > Enzyme gene name\n# f3 > Reaction Ecocyc ID\n";


my @rxns=$cyc->all_rxns(); 	# retrieve all reactions
foreach my $rx (@rxns){
  my @enzymes=$cyc->enzymes_of_reaction($rx);
  foreach my $enz (@enzymes){
    my @genes=$cyc->genes_of_protein($enz); # Note that enzyme might be a complex and several genes will be returned.
      foreach my $gen (@genes){
	push(@{$genenz{$gen}},$rx); # Use of hash eliminates repeated genes.
      }
  }
}

foreach my $gn (sort keys %genenz){
  my @name=$cyc->get_slot_values($gn,"COMMON-NAME");

  # Check if gene name is present in the TRN
  my $qname=quotemeta($name[0]);
  open(TRN,$regulons) || die "Cannot open TRN at $regulons.\n";
  while(<TRN>){
    if($_=~/^[^\t]+\t.*$qname\/\//){
      print OUT"$gn\t$name[0]\t";
      foreach my $rxn (@{$genenz{$gn}}){
	print OUT"$rxn//";
      }
      print OUT"\n";
      last;
    }
  }
  close(TRN);
}
close(OUT);
