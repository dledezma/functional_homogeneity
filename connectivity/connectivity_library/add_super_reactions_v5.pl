#################################################################################
###### 19-10-15                                                            ######
###### Add super-reactions                                                 ######
###### Finds super-reactions at most N reactions away. Copies gu_assembly_vX.pl #
###### output and includes them.                                           ######
###### [0] - gu_assembly_vX.pl output directory                            ######
###### [1] - GU library including files:                                   ######
######     > rxn_rpd_links.txt                                             ######
######     > pwy_rxn_links.txt                                             ######
######     > gu_compound_dictionary.txt                                    ######
###### [2] - Output directory                                              ######
###### [3] - Counts output file path                                       ######
###### [4] - Number of allowed intermediate metabolites (default=9)        ######
###### History:                                                            ######
###### 21/10/15 > v1 > Print as reactants and products only molecules      ######
######              common to reaction in the GU and first not-GU reaction ######
###### 27/10/15 > v2 > Outputs file with super-reactions added per GU      ######
###### 12/11/15 > v2 > Input number of allowed intermediate reactions      ######
###### 04/12/15 > v3 > If no intersection between reactions in the border  ######
######                 of the path is found, print the metabolite already  ######
######                 present in RPD. If both, print first.               ######
###### 05/12/15 > v3 > Excluded "_Ext" suffix in &reactants and &products subs ##
######                 All reactions are now treated the same. Reversible  ######
######                 reactions return both sides.                        ######
###### 04/04/16 > v4 > Included search for reactions from effectors to all ######
######                 other reactions.                                    ######
###### Required files:                                                     ######
#################################################################################

use perlcyc;
$cyc = perlcyc -> new("ECOLI");

$time=localtime();

#Get args
$in_dir=$ARGV[0];
$library=$ARGV[1];
$out_dir=$ARGV[2];
$outfile=$ARGV[3];
$rxns_max=$ARGV[4];
chomp($rxns_max);

#Get name of necessary files.
$rxn_link=$library . "rxn_rpd_links.txt";
$pwy_link=$library . "pwy_rxn_rpd_links.txt";
$dictionary=$library . "gu_compound_dictionary.txt";

#Open outfile
open(OUT,">$outfile") || die "Cannot open OUT file at $outfile.\n";
print OUT"# From $0 on $time\n# Indir = $in_dir\n# Outdir = $out_dir\n# Maximum intermediate reactions = $rxns_max\n";

if(!($rxns_max)){
  $rxns_max=10;
}else{
  $rxns_max++;
}

#Open IN dir
opendir(DIR,$in_dir) || die "Cannot open DIR at $in_dir.\n";
@GU_files=readdir(DIR);
foreach $file (sort @GU_files){
  if(!($file=~/^\./)){

    print"\n**$file\n";
    $TF=$file;

    @all_rxns=();
    $c_rxn=0;
    undef %r_rxns;
    @rpt_rpds=();
    
    #Get file names from GU folder
    $o_objects=$in_dir . $file . "/objects.txt";
    $o_products=$in_dir . $file . "/reactants_products.txt";
    $o_reactions=$in_dir . $file . "/reactions.txt";
    $o_complexes=$in_dir . $file . "/complexes.txt";
    $o_modification=$in_dir . $file . "/modification.txt";

    #Create new dir and copy files
    $new_dir=$out_dir . $file;
    $objects=$new_dir . "/objects.txt";
    $products=$new_dir . "/reactants_products.txt";
    $reactions=$new_dir . "/reactions.txt";
    $complexes=$new_dir . "/complexes.txt";
    $modification=$new_dir . "/modification.txt";
    system("mkdir $new_dir");
    system("cp $o_objects $objects");
    system("cp $o_products $products");
    system("cp $o_reactions $reactions");
    system("cp $o_complexes $complexes");
    system("cp $o_modification $modification");

    #Save all reactions
    open(RXN,$o_reactions) || die "Cannot open RXN file at $reactions.\n";
    while(<RXN>){
      #Get last reaction for adding super_reactions
      if($_=~/^re(\d+)\t/){
	$c_rxn_i=$1;
      }
      if($_=~/^re\d+\t[^\t]+\t[^\t]+\t([^\t]+)$/){   #get ECOCYC reaction ID
	$rx=$1;
	chomp($rx);
	push(@all_rxns,$rx);
      }
    }
    close(RXN);

    #Set reaction counter to print next reaction
    $c_rxn_i++;
    $c_rxn=$c_rxn_i;

    foreach $rx (@all_rxns){      #foreach ECOCYC rxn ID taken from reactions.txt
      @pathways=();
      #Find reaction pathways
      open(LNK,$pwy_link) || die "Cannot open LNK at $pwy_link.\n";
      while(<LNK>){
	if($_=~/^([^\t]+)\t$rx$/){   
	  push(@pathways,$1);      #get pathways of reaction
	}
      }
      close(LNK);

      #Open outfiles
      open(RXN,">>$reactions") || die "Cannot open RXN at $reactions.\n";
      open(RPD,">>$products") || die "Cannot open RPD at $products.\n";

      foreach $pwy (@pathways){
	$threshold=0;        #threshold controls limit length of SR specified in args
	&before($rx,$rx);   #find reactions in pathway $pwy before $rx
	$threshold=0;
	&after($rx,$rx);     #find reactions in pathway $pwy after $rx

      }
	  
    }
   
    ###SRs that connect effector
    ###Find effector
    @effectors=();
    #Open complex file
    open(CPX,$o_complexes) || die "Cannot open CPX at $o_complexes.\n";
    while(<CPX>){
      if($_=~/^csa\d+\t$TF-[^\t]+\t([^\t]+)/i){
	$eff=$1;
	chomp($eff);
	if(!(($eff=~/^phosphate/) || ($eff=~/^ATP/))){
	  if($eff ne $TF){
	    push(@effectors,$eff);
	  }
	}
      }
    }
    close(CPX);

    #Find reactions where effector is present
    foreach $eff (@effectors){
      undef %eff_rxns;
      $q_eff=quotemeta($eff);
      
      open(RLN,$rxn_link) || die "Cannot open RLN at $rxn_link.\n";
      while(<RLN>){
	if($_=~/[^\t]+\t([^\t]+)\t([^\t]+)\t.*\/\/$q_eff\/\//){
	  $rxn=$1;
	  $dir=$2;
	  
	  #Omit reactions already present in GU
	  $flag_eff=0;
	  foreach $rx (@all_rxns){
	    if($rxn eq $rx){
	      $flag_eff=1;
	      last;
	    }
	  }
	  
	  if(!($flag_eff)){
	    $eff_rxns{$rxn}=$dir;
	  }
	}
      }
      close(RLN);
      #Find pathway of reactions found
      undef %eff_pwys;
      foreach $rx (keys %eff_rxns){
	open(PLN,$pwy_link) || die "Cannot open PLN on $pwy_link.\n";
	while(<PLN>){
	  if($_=~/^([^\t]+)\t$rx$/){
	    $pwy=$1;
	    $threshold=0;
	    &eff_before($rx,$rx);
	    $threshold=0;
	    &eff_after($rx,$rx);
	  }
	}
	close(PLN);
      }
    }
	    
    close(RXN);
    close(RPD);	  

    $total_srxns=$c_rxn-$c_rxn_i;
    print OUT"$file\t$total_srxns\n";
    print "$file\t$total_srxns\n";

  }
}
close(DIR);
close(OUT);

 undef %r_rxns;
undef %eff_pwys;
 undef %eff_rxns;
@effectors=();
@pathways=();


####################
#### Subroutine: &before
#### Input> rxn, reaction path
#### Returns> 0
####################

sub before{

  #First reaction
  if($threshold==0){
    $threshold++;                  #threshold increases in 1 because $rx is one reaction

    local($rxn,$path,@before,$bf,$flag_b,$ky)=($_[0],$_[1],(),"",0,"");
    
    @before=$cyc->get_predecessors($rxn,$pwy);    #get reactions before $rxn in same pathway
    foreach $bf (@before){              #foreach reaction before $rxn
    
      #Omit first reactions in GU
      $flag_b=0;                        
      foreach $ky (@all_rxns){        #foreach rxn of all ECOCYC rxn ids in GU
	if($ky eq $bf){               #If immediate rxn is part of the GU, omit - this step is omitted when looking from SRs from effectors
	  $flag_b=1;
	  last;
	}
      }
      if($flag_b){
	next;
      }else{
	&before($bf,$path);       #if immediate before reaction is not in the GU, keep looking in the immediate before reaction
      }
    }


  ##If first reaction was not part of the GU and is second iteration of the function
  #Note: threshold is number of intermediate reactions in path +1 (because it is increased after, making it the intermediate reactions + 2 rxns in the GU:first and last)
  }elsif(($threshold<$rxns_max) && ($threshold>0)){
    $threshold++;

    local($rxn,$path,@before,$bf,$line,$ky,$flag_b,$tpath)=($_[0],$_[1],(),"","","",0,"");
    $path=$rxn . "//" . $path;    #Add new rxn to path

    @before=$cyc->get_predecessors($rxn,$pwy);
    foreach $bf (@before){

      #Find reaction in GU
      $flag_b=0;
      foreach $ky (@all_rxns){                  ##Find if new immediate reaction is part of GU to print SR
	if(($ky eq $bf) && ($ky ne $rx)){       ##(each rxn in GU eq immediate predecessor) && (each rxn in GU ne present $rx)
	  $flag_b=1;
	  $tpath=$bf . "//" . $path;            ##Possible reactions path


	  #Omit repeated super-reactions
	  $flag2=0;                              
	  foreach $rp (keys %r_rxns){           ##Find if path has already been annotated
	    if($rp eq $tpath){
	      $flag2=1;
	    }
	  }

	  if(!($flag2)){	    
	    $r_rxns{$tpath}++;

	    #Get & print substrates and products
	    &Re_Pr($tpath);
	  }
	  last;
	  
	}
      }
      if($flag_b){
	next;
      }else{
	&before($bf,$path);
      }
    }
  }#elsif($threshold==4){
    #print"before: $path\n";
  #}
   
}

####################
#### Subroutine: after
#### Input> rxn, reaction path
#### Returns> 0
####################

sub after{

  #First reaction
  if($threshold==0){
    $threshold++;

    local($rxn,$path,@after,$af,$flag_a,$ky)=($_[0],$_[1],(),"",0,"");
    
    @after=$cyc->get_successors($rxn,$pwy);
    foreach $af (@after){
    
      #Omit first reactions in GU
      $flag_a=0;
      foreach $ky (@all_rxns){
	if($ky eq $af){
	  $flag_a=1;
	  last;
	}
      }
      if($flag_a){
	next;
      }else{
	&after($af,$path);
      }
    }


  #If first reaction is not part of the GU
  }elsif(($threshold<$rxns_max) && ($threshold>0)){
    $threshold++;

    local($rxn,$path,@after,$af,$ky,$flag_a,$tpath)=($_[0],$_[1],(),"","",0,"");
    $path=$path . "//" . $rxn;

    @after=$cyc->get_successors($rxn,$pwy);
    foreach $af (@after){

      #Find reaction in GU
      $flag_a=0;
      foreach $ky (@all_rxns){
	if(($ky eq $af) && ($ky ne $rx)){
	  $flag_a=1;
	  $tpath=$path . "//" . $af;


	  #Omit repeated super-reactions
	  $flag2=0;
	  foreach $rp (keys %r_rxns){
	    if($rp eq $tpath){
	      $flag2=1;
	    }
	  }

	  #If it is not repeated
	  if(!($flag2)){
	    $r_rxns{$tpath}++;	    

	    #Get & print RXN, substrates and products
	    &Re_Pr($tpath);


	  }
	  last;
	}
      }
      if($flag_a){
	next;
      }else{
	&after($af,$path);
      }
    }
  }#elsif($threshold==4){
    #print"after: $path\n";
  #}
    


}



####################
#### Subroutine: &eff_before
#### Same as &before but omits the filtering step in case that first immediate reaction is in GU.These cases are printed since they link effector to the rest of the GU
#### Input> rxn, reaction path
#### Returns> 0
####################

sub eff_before{

if($threshold<$rxns_max){
    $threshold++;

    local($rxn,$path,@before,$bf,$line,$ky,$flag_b,$tpath)=($_[0],$_[1],(),"","","",0,"");


    if($threshold>1){                ##Avoid repetition of first rxn in path
      $path=$rxn . "//" . $path;    #Add new rxn to path
    }
    @before=$cyc->get_predecessors($rxn,$pwy);
    foreach $bf (@before){

      #Find reaction in GU
      $flag_b=0;
      foreach $ky (@all_rxns){                  ##Find if new immediate reaction is part of GU to print SR
	if(($ky eq $bf) && ($ky ne $rx)){       ##(each rxn in GU eq immediate predecessor) && (each rxn in GU ne present $rx)
	  $flag_b=1;
	  $tpath=$bf . "//" . $path;            ##Possible reactions path


	  #Omit repeated super-reactions
	  $flag2=0;                              
	  foreach $rp (keys %r_rxns){           ##Find if path has already been annotated
	    if($rp eq $tpath){
	      $flag2=1;
	    }
	  }

	  if(!($flag2)){	    
	    $r_rxns{$tpath}++;

	    #Get & print substrates and products
	    &eff_Re_Pr("before",$tpath);
	    
	  }
	  last;
	  
	}
      }
      if($flag_b){
	next;
      }else{
	&before($bf,$path);
      }
    }
  }#elsif($threshold==4){
    #print"before: $path\n";
  #}
   
}


####################
#### Subroutine: eff_after
#### Same as &before but omits the filtering step in case that first immediate reaction is in GU.These cases are printed since they link effector to the rest of the GU
#### Input> rxn, reaction path
#### Returns> 0
####################

sub eff_after{

  if($threshold<$rxns_max){
    $threshold++;

    local($rxn,$path,@after,$af,$ky,$flag_a,$tpath)=($_[0],$_[1],(),"","",0,"");

    if($threshold>1){              ##Avoid repetition of first rxn in path
      $path=$path . "//" . $rxn;
    }

    @after=$cyc->get_successors($rxn,$pwy);
    foreach $af (@after){

      #Find reaction in GU
      $flag_a=0;
      foreach $ky (@all_rxns){
	if(($ky eq $af) && ($ky ne $rx)){
	  $flag_a=1;
	  $tpath=$path . "//" . $af;


	  #Omit repeated super-reactions
	  $flag2=0;
	  foreach $rp (keys %r_rxns){
	    if($rp eq $tpath){
	      $flag2=1;
	    }
	  }

	  #If it is not repeated
	  if(!($flag2)){
	    $r_rxns{$tpath}++;

	    #Get & print RXN, substrates and products
	    &eff_Re_Pr("after",$tpath);


	  }
	  last;
	}
      }
      if($flag_a){
	next;
      }else{
	&after($af,$path);
      }
    }
  }#elsif($threshold==4){
    #print"after: $path\n";
  #}
    


}



###########################
### Subroutine Re_Pr
###
### Finds metabolites common to 1st and 2nd reactions, prints them as substrates.
### Finds metabolites common to last and 2nd to last reactions,
### prints them as products.
### -Input arguments:
###  [0] - path
### -Returned values: none
###########################

sub Re_Pr{

  local($path)=($_[0]);

  #Initialize rpd variables (to avoid different paths, same RPD)
  $rxn_rpd="";
  undef %rxn_rpds;

  #Get individual reactions
  @s_rxns=();
  @s_rxns=split(/\/\//,$path);
  
  #Obtain reactions to compare.
  $first=$s_rxns[0];
  $second=$s_rxns[1];
  $last=pop(@s_rxns);
  $pre_last=pop(@s_rxns);

  ###Get first and second intersection
  #Get first reaction products
  @fst_prods=();
  $ref=&products($first);
  @fst_prods=@{$ref};
  
  #Get second reaction reactants
  @snd_rects=();
  $ref=&reactants($second);
  @snd_rects=@{$ref};

  #Find & print intersections
  $f_not=0;
  foreach $r (@snd_rects){
    foreach $p (@fst_prods){
      if($r eq $p){
	$f_not=1;
	#Save reaction met to compare and hash to later print
	$rxn_rpd=$rxn_rpd . "//" . $r;
	push(@{$rxn_rpds{reactants}},$r);
      }
    }
  }

  ###If no intersection was found, print the metabolite already in OBJ
  if(!($f_not)){
    print"++WARNING! Initial intersection not found in path:\n\t$path\n\tWill print already present metabolites\n";
    #Find all reactants and products in objects
    foreach $r (@snd_rects){
      $f_found=0;
      $q_r=quotemeta($r);
      #Open OBJ
      open(FBJ,$objects) || die "Cannot open FBJ at $objects.\n";
      while(<FBJ>){
	if($_=~/^$q_r\tSIMPLE/i){
	  #Save reaction met to compare and hash to later print
	  $rxn_rpd=$rxn_rpd . "//" . $r;
	  push(@{$rxn_rpds{reactants}},$r);
	  $f_found=1;
	  last;
	}
      }
      close(FBJ);
      if($f_found){
	last;
      }
    }

    if(!($f_found)){
      foreach $p (@fst_prods){
	$q_p=quotemeta($p);
	#Open OBJ
	open(FBJ,$objects) || die "Cannot open FBJ at $objects.\n";
	while(<FBJ>){
	  if($_=~/^$q_p\tSIMPLE/i){
	    #Save reaction met to compare and hash to later print
	    $rxn_rpd=$rxn_rpd . "//" . $p;
	    push(@{$rxn_rpds{reactants}},$p);
	    $f_found=1;
	    last;
	  }
	}
	close(FBJ);
	if($f_found){
	  last;
	}
      }
    }

    if(!($f_found)){
      ##If neither is already present
      print"WARNING!\n\tPath:\n\t$path\n\t does not have an overlap on its initial borders. Reactant or product will be missing in RPD.\n";
    }
  }

  ###Get last and second to last intersection
  #Get second to last reaction products
  @fst_prods=();
  $ref=&products($pre_last);
  @fst_prods=@{$ref};
  
  #Get last reaction reactants
  @snd_rects=();
  $ref=&reactants($last);
  @snd_rects=@{$ref};

  #Find & print intersections
  $f_not=0;
  foreach $r (@snd_rects){
    foreach $p (@fst_prods){
      if($r eq $p){
	$f_not=1;
	#Save reaction met to compare and hash to later print
	$rxn_rpd=$rxn_rpd . "//" . $r;
	push(@{$rxn_rpds{products}},$r);
      }
    }
  } 
  
  ###If no intersection was found, print the metabolite already in OBJ
  if(!($f_not)){
    print"++WARNING! Final intersection not found in path:\n\t$path\n\tWill print already present metabolites\n";
    #Find all reactants and products in objects
    foreach $r (@snd_rects){
      $f_found=0;
      $q_r=quotemeta($r);
      #Open OBJ
      open(FBJ,$objects) || die "Cannot open FBJ at $objects.\n";
      while(<FBJ>){
	if($_=~/^$q_r\tSIMPLE/i){
	  #Save reaction met to compare and hash to later print
	  $rxn_rpd=$rxn_rpd . "//" . $r;
	  push(@{$rxn_rpds{products}},$r);
	  $f_found=1;
	  last;
	}
      }
      close(FBJ);
      if($f_found){
	last;
      }
    }

    if(!($f_found)){
      foreach $p (@fst_prods){
	$q_p=quotemeta($p);
	#Open OBJ
	open(FBJ,$objects) || die "Cannot open FBJ at $objects.\n";
	while(<FBJ>){
	  if($_=~/^$q_p\tSIMPLE/i){
	    #Save reaction met to compare and hash to later print
	    $rxn_rpd=$rxn_rpd . "//" . $p;
	    push(@{$rxn_rpds{products}},$p);
	    $f_found=1;
	    last;
	  }
	}
	close(FBJ);
	if($f_found){
	  last;
	}
      }
    }
      

    if(!($f_found)){
      ##If neither is already present
      print"WARNING!\n\tPath:\n\t$path\n\tdoes not have an overlap on its final borders. Reactant or product will be missing in RPD.\n";
    }
  }


  #Ignore repeated RPD's
  $flag_rp=0;
  foreach $rp (@rpt_rpds){
    if($rxn_rpd eq $rp){
      $flag_rp=1;
    }
  }

  #Check if reactants=products
  $f_eq_rpd=0;
  $r_size=(@{$rxn_rpds{reactants}});
  $p_size=(@{$rxn_rpds{products}});
  $t_size=0;
  foreach $r (@{$rxn_rpds{reactants}}){
    foreach $p (@{$rxn_rpds{products}}){
      if($r eq $p){
	$t_size++;
      }
    }
  }
  #If all reactants are a subgroup of products or viceversa, omit (not informative).
  if(($t_size == $p_size) || ($t_size == $r_size)){
    $f_eq_rpd=1;
  }
  
  #If reaction has not been printed or reactants=products
  if((!($flag_rp)) && (!($f_eq_rpd))){
    #Print reaction
    print RXN"re$c_rxn\tSUPER_REACTION\t$path\n";

    #Print reactants & products
    foreach $r (@{$rxn_rpds{reactants}}){
      print RPD"re$c_rxn\t$r\treactant\n";
    }
    foreach $r (@{$rxn_rpds{products}}){
      print RPD"re$c_rxn\t$r\tproduct\n";
    }
    push(@rpt_rpds,$rxn_rpd);
    $c_rxn++;
  }



}



###########################
### Subroutine eff_Re_Pr
###
### Differs in the fact that one border is not present in the GU so intersection of reactants is not the reactant
### Finds metabolites common to 1st and 2nd reactions, prints them as substrates.
### Finds metabolites common to last and 2nd to last reactions,
### prints them as products.
### -Input arguments:
###  [0] - path
### -Returned values: none
###########################

sub eff_Re_Pr{

  local($type,$path,$error)=($_[0],$_[1],1);
  print"$path\n";

  #Initialize rpd variables (to avoid different paths, same RPD)
  $rxn_rpd="";
  undef %rxn_rpds;

  #Get individual reactions
  @s_rxns=();
  @s_rxns=split(/\/\//,$path);
  
  #Obtain reactions to compare.
  $first=$s_rxns[0];
  $second=$s_rxns[1];
  $last=pop(@s_rxns);
  $pre_last=pop(@s_rxns);

  if($type eq "before"){
    #Find effector in products of last reaction
    @fst_prods=();
    $ref=&products($last);
    @fst_prods=@{$ref};
    #Find effector of this iteration in products
    foreach $p (@fst_prods){
      if($p eq $eff){
	$error=0;
	#Save reaction met to compare and hash to later print
	$rxn_rpd=$rxn_rpd . "//" . $p;
	push(@{$rxn_rpds{products}},$p);
      }
    }
  
    #Warning for unidentified effectors
    if($error){
      print"WARNING! Effector not identified in effector SR. Will not be printed.\n";
    }
    
    ###Get first and second intersection
    #Get first reaction products
    @fst_prods=();
    $ref=&products($first);
    @fst_prods=@{$ref};
  
    #Get second reaction reactants
    @snd_rects=();
    $ref=&reactants($second);
    @snd_rects=@{$ref};

    #Find & print intersections
    $f_not=0;
    foreach $r (@snd_rects){
      foreach $p (@fst_prods){
	if($r eq $p){
	  $f_not=1;
	  #Save reaction met to compare and hash to later print
	  $rxn_rpd=$rxn_rpd . "//" . $r;
	  push(@{$rxn_rpds{reactants}},$r);
	}
      }
    }

    ###If no intersection was found, print the metabolite already in OBJ
    if(!($f_not)){
      print"++WARNING! Initial intersection not found in path:\n\t$path\n\tWill print already present metabolites\n";
      #Find all reactants and products in objects
      foreach $r (@snd_rects){
	$f_found=0;
	$q_r=quotemeta($r);
	#Open OBJ
	open(FBJ,$objects) || die "Cannot open FBJ at $objects.\n";
	while(<FBJ>){
	  if($_=~/^$q_r\tSIMPLE/i){
	    #Save reaction met to compare and hash to later print
	    $rxn_rpd=$rxn_rpd . "//" . $r;
	    push(@{$rxn_rpds{reactants}},$r);
	    $f_found=1;
	    last;
	  }
	}
	close(FBJ);
	if($f_found){
	  last;
	}
      }
      
      if(!($f_found)){
	foreach $p (@fst_prods){
	  $q_p=quotemeta($p);
	  #Open OBJ
	  open(FBJ,$objects) || die "Cannot open FBJ at $objects.\n";
	  while(<FBJ>){
	    if($_=~/^$q_p\tSIMPLE/i){
	      #Save reaction met to compare and hash to later print
	      $rxn_rpd=$rxn_rpd . "//" . $p;
	      push(@{$rxn_rpds{reactants}},$p);
	      $f_found=1;
	      last;
	    }
	  }
	  close(FBJ);
	  if($f_found){
	    last;
	  }
	}
      }
      
      if(!($f_found)){
	##If neither is already present
	print"WARNING!\n\tPath:\n\t$path\n\t does not have an overlap on its initial borders. Reactant or product will be missing in RPD.\n";
      }
    }
    
    #Ignore repeated RPD's
    $flag_rp=0;
    foreach $rp (@rpt_rpds){
      if($rxn_rpd eq $rp){
	$flag_rp=1;
      }
    }

    #Check if reactants=products
    $r_size=(@{$rxn_rpds{reactants}});
    $p_size=(@{$rxn_rpds{products}});
    $t_size=0;
    foreach $r (@{$rxn_rpds{reactants}}){
      foreach $p (@{$rxn_rpds{products}}){
	if($r eq $p){
	  $t_size++;
	}
      }
    }
    #If all reactants are a subgroup of products or viceversa, omit (not informative).
    if(($t_size == $p_size) || ($t_size == $r_size)){
      $flag_rp=1;
    }

    #If reaction has not been printed or reactants=products
    if((!($flag_rp)) && (!($error))){
      
      #Print reaction
      print RXN"re$c_rxn\tSUPER_REACTION\t$path\n";
      
      #Print reactants & products
      foreach $r (@{$rxn_rpds{reactants}}){
	print RPD"re$c_rxn\t$r\treactant\n";
      }
      foreach $r (@{$rxn_rpds{products}}){
	print RPD"re$c_rxn\t$r\tproduct\n";
      }
      push(@rpt_rpds,$rxn_rpd);
      $c_rxn++;
    }
  
  }


  if($type eq "after"){
    #Find effector in reactant of first reaction
    #Get last reaction reactants
    @fst_rects=();
    $ref=&reactants($first);
    @fst_rects=@{$ref};
    #Find effector of this iteration in reactants
    foreach $r (@fst_rects){
      if($r eq $eff){
	 print"$r - $eff\n";
	$error=0;
	#Save reaction met to compare and hash to later print
	$rxn_rpd=$rxn_rpd . "//" . $r;
	push(@{$rxn_rpds{reactants}},$r);
      }
    }
    
    #Warning for unidentified effectors
    if($error){
      print"WARNING! Effector $eff not identified in effector SR: $path.\n !!!Will not be printed.\n";
    }

    ###Get last and second to last intersection
    #Get second to last reaction products
    @fst_prods=();
    $ref=&products($pre_last);
    @fst_prods=@{$ref};
  
    #Get last reaction reactants
    @snd_rects=();
    $ref=&reactants($last);
    @snd_rects=@{$ref};

    #Find & print intersections
    $f_not=0;
    foreach $r (@snd_rects){
      foreach $p (@fst_prods){
	if($r eq $p){
	  $f_not=1;
	  #Save reaction met to compare and hash to later print
	  $rxn_rpd=$rxn_rpd . "//" . $r;
	  push(@{$rxn_rpds{products}},$r);
	}
      }
    } 
    
    ###If no intersection was found, print the metabolite already in OBJ
    if(!($f_not)){
      print"++WARNING! Final intersection not found in path:\n\t$path\n\tWill print already present metabolites\n";
      #Find all reactants and products in objects
      foreach $r (@snd_rects){
	$f_found=0;
	$q_r=quotemeta($r);
	#Open OBJ
	open(FBJ,$objects) || die "Cannot open FBJ at $objects.\n";
	while(<FBJ>){
	  if($_=~/^$q_r\tSIMPLE/i){
	    #Save reaction met to compare and hash to later print
	    $rxn_rpd=$rxn_rpd . "//" . $r;
	    push(@{$rxn_rpds{products}},$r);
	    $f_found=1;
	    last;
	  }
	}
	close(FBJ);
	if($f_found){
	  last;
	}
      }
      
      if(!($f_found)){
	foreach $p (@fst_prods){
	  $q_p=quotemeta($p);
	  #Open OBJ
	  open(FBJ,$objects) || die "Cannot open FBJ at $objects.\n";
	  while(<FBJ>){
	    if($_=~/^$q_p\tSIMPLE/i){
	      #Save reaction met to compare and hash to later print
	      $rxn_rpd=$rxn_rpd . "//" . $p;
	      push(@{$rxn_rpds{products}},$p);
	      $f_found=1;
	      last;
	    }
	  }
	  close(FBJ);
	  if($f_found){
	    last;
	  }
	}
      }
      
      
      if(!($f_found)){
	##If neither is already present
	print"WARNING!\n\tPath:\n\t$path\n\tdoes not have an overlap on its final borders. Reactant or product will be missing in RPD.\n";
      }
    }
    
    
    #Ignore repeated RPD's
    $flag_rp=0;
    foreach $rp (@rpt_rpds){
      if($rxn_rpd eq $rp){
	$flag_rp=1;
      }
    }
    
    #Check if reactants=products
    $r_size=(@{$rxn_rpds{reactants}});
    $p_size=(@{$rxn_rpds{products}});
    $t_size=0;
    foreach $r (@{$rxn_rpds{reactants}}){
      foreach $p (@{$rxn_rpds{products}}){
	if($r eq $p){
	  $t_size++;
	}
      }
    }
    #If all reactants are a subgroup of products or viceversa, omit (not informative).
    if(($t_size == $p_size) || ($t_size == $r_size)){
      $flag_rp=1;
    }

    #If reaction has not been printed or reactants=products
    if((!($flag_rp)) && (!($error))){

      #Print reaction
      print RXN"re$c_rxn\tSUPER_REACTION\t$path\n";
      
      #Print reactants & products
      foreach $r (@{$rxn_rpds{reactants}}){
	print RPD"re$c_rxn\t$r\treactant\n";
      }
      foreach $r (@{$rxn_rpds{products}}){
	print RPD"re$c_rxn\t$r\tproduct\n";
      }
      push(@rpt_rpds,$rxn_rpd);
      $c_rxn++;
    }

  }


}


###########################
### Subroutine reactants
###
### Finds reaction reactants
### -Input arguments:
###  [0] - reaction
### -Returned values: ref to array with reactants
### Notes:
###   - "_Ext" suffix is never added since only metabolites inside the cell can be compared.
###   - REVERSIBLE transport reactions are included since direction of the flux is unknown. 
###########################

sub reactants{
  
  local($rx)=($_[0]);
  @all_rects=();

  @f_phyr=$cyc->get_slot_values($rx,"REACTION-DIRECTION");
  if($f_phyr[0]=~/LEFT-TO-RIGHT|REVERSIBLE/){
    @subs=();
    @subs=$cyc->get_slot_values($rx,"LEFT");
    foreach $s (@subs){
      if(!($s=~/^PROTON$|^AMP$|^ATP$|^ADP$|^WATER$|^NAD.*$|^\|*Pi\|*$/i)){
	$s_nm=&IDS($s);
	push(@all_rects,$s_nm);
      }
    }
  }
  if($f_phyr[0]=~/RIGHT-TO-LEFT|REVERSIBLE/){
    @subs=();
    @subs=$cyc->get_slot_values($rx,"RIGHT");
    foreach $s (@subs){
      if(!($s=~/^PROTON$|^AMP$|^ATP$|^ADP$|^WATER$|^NAD.*$|^\|*Pi\|*$/i)){
	$s_nm=&IDS($s);
	push(@all_rects,$s_nm);
      }
    }
  }
  return(\@all_rects);

}

###########################
### Subroutine products
###
### Finds reaction products
### -Input arguments:
###  [0] - reaction
### -Returned values: ref to array with products
###########################


sub products {
  
  local($rx)=($_[0]);
  @all_prods=();

  @f_phyr=$cyc->get_slot_values($rx,"REACTION-DIRECTION");
  if($f_phyr[0]=~/LEFT-TO-RIGHT|REVERSIBLE/){
    @prods=();
    @prods=$cyc->get_slot_values($rx,"RIGHT");
    foreach $p (@prods){
      if(!($p=~/^PROTON$|^AMP$|^ATP$|^ADP$|^WATER$|^NAD.*$|^\|*Pi\|*$/i)){
	$p_nm=&IDS($p);
	push(@all_prods,$p_nm);
      }
    } 
  }
  if($f_phyr[0]=~/RIGHT-TO-LEFT|REVERSIBLE/){
    @prods=();
    @prods=$cyc->get_slot_values($rx,"LEFT");
    foreach $p (@prods){
      if(!($p=~/^PROTON$|^AMP$|^ATP$|^ADP$|^WATER$|^NAD.*$|^\|*Pi\|*$/i)){
	$p_nm=&IDS($p);
	push(@all_prods,$p_nm);
      }
    }
  }
  return(\@all_prods);
}
    


###########################
### IDS Function
###
### Turn metabolites' Ecocyc id's into common names
### -Input arguments:
###  [0] - ecocyc ID
### -Returned values: common name
###########################

sub IDS {
 
  local($name)=($_[0]);

  $q_name=quotemeta($name);
  open(SMS,$dictionary) || die "gu_compound_dictionary.txt file not found on dictionary. Please add.\n";
  while(<SMS>){
    if($_=~/^$q_name\t([^\t]+)\t/i){
      $name=$1;
      $name=&xml($name);
      return($name);
    }
  }
  if($name=~/^\|(.+)\|$/){
    return($1);
  }
  print"Small Molecule $name not found.\n";
  return($name);
}

###########################
### XML Function
###
### Remove XML code from molecule names
### -Input arguments:
###  [0] - line
### -Returned values: line in a single char value
###########################

sub xml{

  local($line)=$_[0];
  
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
  while($common=~/(.*)&omega;(.*)/i){
    $common=$1 . "omega" . $2;
  }
  return($line);
  
}

