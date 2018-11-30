#################################################################################
###### 01-02-18                                                            ######
###### GU Assembly From Groups with Enzymatic Regulation                   ######
###### Assemblies GUs (relational tables) for groups of genes in reg format #####
###### Changes over GU_Assembly:                                           ######
###### - Needs groups of genes                                             ######
###### - Ignores TFs                                                       ######
###### Argument usage:                                                     ######
###### [0] - Directory with genes in regulon file:                         ######
######       f1 - regulon name                                             ######
######       f2 - Genes separated by //                                    ######
###### [1] - Output directory. A new directory will be created per regulon ######
######       each containing 5 relational tables.                          ######
###### [2] - Library dir path. Must include at least: TF_active, TF_inactive, ###
######       gene_links, tu_TRN, tu_complex tu_links and                   ######
######       gu_compound_dictionary files                                  ######
###### History:                                                            ######
#################################################################################

use perlcyc;

###Open pathway tools, create new e.coli session.
#system("ptools -api -lisp");
$cyc = perlcyc -> new("ECOLI");

###Get args
$dir_in=$ARGV[0];
$dir_out=$ARGV[1];
$lib=$ARGV[2];
chomp($lib);

### Get library files
$tf_active=$lib . "TF_active.txt";
$tf_inactive=$lib . "TF_inactive.txt";
$gene_links=$lib . "gene_links.txt";
$dictionary=$lib . "gu_compound_dictionary.txt";
$tu_trn=$lib . "tu_TRN.txt";
$tu_links=$lib . "tu_links.txt";
$tu_complex=$lib . "tu_complex.txt";

##Check correct mode usage.
if(($dir_out) && ($dir_in) && ($lib)){
  open(DFL,$dir_in) || die "Cannot open DFL at $dir_in.\n";
  while(<DFL>){
    if($_=~/^\#/){
      next;
    }
    if($_=~/^([^\t]+)\t([^\t]+)$/){
      $name=$1;
      $genelist=$2;
      chomp($genelist);
      print"\n****$name****\n";

      ###Name output directory.
      $outdir=$dir_out . $name . "/";

      ###Create output directory.
      system("mkdir $outdir");
      $TF=$name;
      # Send genes and output file to sub
      &ecocyc($genelist,$outdir);
    }
  }
  close(DFL);
}else{
  die "Argument usage:\n[1] - Mode: must be \"d\" or \"s\"\nd--Directory mode. Build GU's for every regulon file in directory.\ns--Single mode. Build GU for a single regulon file.\n[2] - Read directory/file.\n[3] - Working directory.\n";
}






###########################
### Ecocyc Function
###
### - Looks for gene products, complexes and the reactions they catalyze.
### - Prints:
###     +Gene Products/Complexes
###     +Reactions
### Arguments: 
###        [0]> Gene vector ref 
###        [1]> Output dir
### Returns:
###        Termination value (1).
########################### 

sub ecocyc {

  local(%r_cpx, %r_prots);
  local($c_rxn, $c_cpx, @genes)=(1,1,());
  local($genelist,$d_out)=($_[0],$_[1]);

  undef %r_cpx;
  undef %r_prots;

  ###Create temporal output file.
  $temp=$d_out . "data_temp.txt";
  open(TMP,">$temp") || die "Cannot open TMP file at $temp\n";

  ##Get genes
  chop($genelist);
  chop($genelist); 		# Delete final diagonals on gene list
  @allgenelist=split(/\/\//,$genelist);
  foreach $gn (@allgenelist){
    $g=&IDS($gn,"g");
    if(!($g)){
      print"WARNING! $gn is not a gene and will not be considered.\n";
      next;
    }else{
      print"Gene $gn recognized as $g.\n";
      push(@genes,$g);
    }
  }
  close(GIN);
    
    
  foreach $g (@genes){
    $tu=&IDS($g,"gi");

    #print"started $tu.\n";
    @gprods=();
    @gprods=$cyc->all_products_of_gene($g);
    foreach $en (@gprods){
      #print"Enzyme= $en\n";
      #Ignore if product is a complex (will appear below from a monomer)
      #Needed for identifying the regulated monomer
      @mons=();
      @mons=$cyc->monomers_of_protein($en);
      $size=@mons;
      #Separate homomeric from heteromultimeric complexes. Ignore hetero-.
      if($size>1){
	next;
      }
      #Print enzyme; needed for translation reaction
      print TMP"$tu\t$en\n";
      #Check if product is part of a complex
      @cpxs=();
      @cpxs=$cyc->get_slot_values($en,"COMPONENT-OF");
      foreach $cp (@cpxs){
	@mons=();
	$all_mons="";
	@mons=$cyc->monomers_of_protein($cp);
	$size=@mons;
	if($size>1){
	  foreach $mn (@mons){
	    $all_mons=$all_mons . "$mn,";
	  }
	  print TMP"$tu\t$en\t$cp\t$all_mons\n";
	  $ref=&enzyme($cp);
	  @rxns=(@{$ref});
	  foreach $rx (@rxns){
	    print TMP"$tu\t$en\t$cp\t\t$rx\n";
	  }
	}else{
	  $ref=&enzyme($cp);
	  @rxns=(@{$ref});
	  foreach $rx (@rxns){
	    print TMP"$tu\t$en\t\t\t$rx\n";
	  }
	}
      }
      #Check for product reactions
      $ref=&enzyme($en);
      @rxns=(@{$ref});
      foreach $rx (@rxns){
	print TMP"$tu\t$en\t\t\t$rx\n";
      }
    }
  }
  close(TMP);

  ###Create output files.
  $f_objects=$d_out . "ob_temp.txt";
  $f_reactions=$d_out . "reactions.txt";
  $f_complexes=$d_out . "complexes.txt";
  $f_recprod=$d_out . "reactants_products_temp.txt";
  $f_modification=$d_out . "t_modification.txt";
  $uniq_objects=$d_out . "objects.txt";
  $uniq_mods=$d_out . "modification.txt";

  ##Open output files.
  open(OBJ,">$f_objects") || die "Cannot open OBJ file for $TF.\n";
  open(RXN,">$f_reactions") || die "Cannot open RXN file for $TF.\n";
  open(CPX,">$f_complexes") || die "Cannot open CPX file for $TF.\n";
  open(MOD,">$f_modification") || die "Cannot open MOD file for $TF.\n";
  
  ##Create reactants & products temporal file, print trnascription reactions and close (it will be appended inside the next loop).
  open(RPD,">$f_recprod") || die "Cannot create RPD file for $TF at $f_recprod.\n";


  # Print transcription reactions
  foreach $g (@genes){
    $tu=&IDS($g,"gi");
    
    print RXN"re$c_rxn\tTRANSCRIPTION\n";
    print RPD"re$c_rxn\t$tu\treactant\n";
    print RPD"re$c_rxn\t$tu\_mRNA\tproduct\n";
    print OBJ"$tu\tGENE\n$tu\_mRNA\tRNA\n";
    $c_rxn++;
  }

  close(RPD);

  undef %r_rxn;


  # Print all else reactions.
  open(TMP,$temp) || die "Cannot open TMP file at $temp\n";
  while(<TMP>){

    #Open RPD inside the loop to avoid conflicts in repeated reactions.
    open(RPD,">>$f_recprod") || die "Cannot open RPD file for $TF at $f_recprod.\n";

    ## Print translation reactions
    if($_=~/^([^\t]+)\t([^\t]+)$/){
      $tu=$1;
      $pr_id=$2;
      chomp($pr_id);
      
      #Get protein name
      $pr_nm=&IDS($pr_id,"p");

      #Ignore repeated proteins from the same TU
      $flag=0;
      $flag1=0;
      foreach $rp_tu (@{$r_prots{$pr_nm}}){
	if($rp_tu eq $tu){
	  $flag=1;
	}
	if($rp_tu eq "EFF"){
	  $flag1=1;
	}
      }
      if($flag){
	next;
      }

      push(@{$r_prots{$pr_nm}},$tu);
      print RXN"re$c_rxn\tTRANSLATION\n";
      #print"re$c_rxn\n";                                                                     ########### Uncomment this line for verbose
      print RPD"re$c_rxn\t$tu\_mRNA\treactant\n";
      print RPD"re$c_rxn\t$pr_nm\tproduct\n";
      #Ignore proteins already annotated in &effectors (avoid duplicated OBJ entries)
      if(!($flag1)){
	print OBJ"$pr_nm\tPROTEIN\t$pr_id\n";
      }
      $c_rxn++;
    }

    ## Print complexes
    if($_=~/^[^\t]+\t[^\t]+\t([^\t]+)\t([^\t]+)$/){
      $cpx_id=$1;
      $monos=$2;
      chomp($monos);

      #Get complex name
      $cpx_nm=&IDS($cpx_id,"cx");
      
      #Ignore repeated complexes
      $flag=0;
      foreach $rp_cp (keys %r_cpx){
	if($rp_cp eq $cpx_nm){
	  $flag=1;
	  last;
	}
      }
      if($flag){
	next;
      }
      $r_cpx{$cpx_nm}++;

      ## Print complex formation reaction
      print RXN"re$c_rxn\tSTATE_TRANSITION\tL2R\n";
      print OBJ"$cpx_nm\tCOMPLEX\t$cpx_id\n";
      @monos=split(/,/,$monos);
      foreach $mono_id (@monos){
	$mono_nm=&IDS($mono_id,"p");

	#Omit OBJ print of proteins printed in &effectors
	$flag1=0;
	foreach $rp_pr (@{$r_prots{$pr_nm}}){
	  if($rp_pr eq "EFF"){
	    $flag1=1;
	  }
	}
	if(!($flag1)){
	  print OBJ"$mono_nm\tPROTEIN\t$mono_id\n";
	}
	print RPD"re$c_rxn\t$mono_nm\treactant\n";
	
	## Print complex id
	print CPX"csa$c_cpx\t$cpx_nm\t$mono_nm\n"
      }
      print RPD"re$c_rxn\t$cpx_nm\tproduct\n";
      $c_rxn++;
      $c_cpx++;
    }
      
    ## Print reactions
    if($_=~/^[^\t]+\t([^\t]+)\t([^\t]*)\t\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)$/){
      $mono=$1;
      $cpx=$2;
      $rxn_id=$3;
      $rxn_dir=$4;
      $reactants=$5;
      $products=$6;
      chomp($products);

      #Identify catalyzer
      if($cpx){
	$cat_id=$cpx;
	$cat_nm=&IDS($cat_id,"cx");
      }else{
	$cat_id=$mono;
	$cat_nm=&IDS($cat_id,"p");
      }

      #Ignore repeated reactions from same catalyzer 
      $flag=0;
      foreach $rp_cat (@{$r_rxn{$rxn_id}}){
	if($rp_cat eq $cat_nm){
	  $flag=1;
	  last;
	}
      }
      if($flag){
	next;
      }
      push(@{$r_rxn{$rxn_id}},$cat_nm);  

      ###Print modification of repeated reactions
      $flag=0;
      undef %rxn_index;

      ##Look for repeated reactions in RPD
      #Get substrates and products of reaction
      #Get individual substrates & products
      @substrates=split(/,/,$reactants);
      @products=split(/,/,$products);

      foreach $s (@substrates){
	if($s=~/^(.+)_Ext$/){
	  $sm=$1;
	  $s_nm=&IDS($sm,"sm");
	  $s_nm=$s_nm . "_Ext" . "//" . $s;
	}else{
	  $s_nm=&IDS($s,"sm");
	  $s_nm=$s_nm . "//" . $s;
	}
	push(@{$rxn_index{"substrates"}},$s_nm);
      }

      foreach $p (@products){
	if($p=~/^(.+)_Ext$/){
	  $sm=$1;
	  $p_nm=&IDS($sm,"sm");
	  $p_nm=$p_nm . "_Ext" . "//" . $p;
	}else{
	  $p_nm=&IDS($p,"sm");
	  $p_nm=$p_nm . "//" . $p;
	}
	push(@{$rxn_index{"products"}},$p_nm);
      }

      close(RPD);

      ##Compare to RPD reactions by sm name
      #Get RPD substrates and products by reaction
      for($inx=1;$inx<$c_rxn;$inx++){
	@rpd_subs=();
	@rpd_prods=();
	$fhl="fhl" . "$inx";
	open($fhl,$f_recprod) || die "Cannot open RPD1 at $f_recprod\n";
	while(<$fhl>){
	  if($_=~/^re$inx\t([^\t]+)\treactant/){
	    push(@rpd_subs,$1);
	  }elsif($_=~/^re$inx\t([^\t]+)\tproduct/){
	    push(@rpd_prods,$1);
	  }
	}
	close($fhl);
	
	#Compare reaction/products length
	if((@rpd_subs == @substrates) && (@rpd_prods == @products)){

	  #Compare reactants
	  $subs=@substrates;
	  $c_subs=0;
	  foreach $s (@{$rxn_index{"substrates"}}){
	    foreach $r (@rpd_subs){
	      $q_r=quotemeta($r);
	      if($s=~/^$q_r\/\/.+/i){
		$c_subs++;
	      }
	    }
	  }
	  #Compare products
	  $prods=@products;
	  $c_prods=0;
	  foreach $p (@{$rxn_index{"products"}}){
	    foreach $r (@rpd_prods){
	      $q_r=quotemeta($r);
	      if($p=~/^$q_r\/\/.+/i){
		$c_prods++;
	      }
	    }
	  }
	  
	  #Compare intersections
	  if(($c_subs == $subs) && ($c_prods == $prods)){
	    #Print modification of reaction
	    print MOD"re$inx\tCATALYSIS\t$cat_nm\n";
	    $flag=1;
	    last;
	  }
	}
      }
	 
      #If reaction has not been printed
      if(!($flag)){
	#Print reaction
	print MOD"re$c_rxn\tCATALYSIS\t$cat_nm\n";
	#Identify transport reactions
	if(($reactants=~/_Ext/) || ($products=~/_Ext/)){
	  print RXN"re$c_rxn\tTRANSPORT\t$rxn_dir\t$rxn_id\n";
	}else{
	  print RXN"re$c_rxn\tSTATE_TRANSITION\t$rxn_dir\t$rxn_id\n";
	}
	
	##Print Reactants & Products
	open(RPD,">>$f_recprod") || die "Cannot open RPD file for $TF at $f_recprod.\n";
	foreach $s (@{$rxn_index{"substrates"}}){
	  if($s=~/^(.+)\/\/(.+)$/){
	    print OBJ"$1\tSIMPLE_MOLECULE\t$2\n";
	    print RPD"re$c_rxn\t$1\treactant\n";
	  }
	}
	foreach $p (@{$rxn_index{"products"}}){
	  if($p=~/^(.+)\/\/(.+)$/){
	    print OBJ"$1\tSIMPLE_MOLECULE\t$2\n";
	    print RPD"re$c_rxn\t$1\tproduct\n";
	  }
	}
	$c_rxn++;
	close(RPD);
      }
      

    }

  }
  close(TMP);
  close(OBJ);
  close(RXN);
  close(CPX);
  close(RPD);
  close(MOD);


  ### Sort reactants & products by length in RPD

  # Open OUT file
  $f_recprod2=$d_out . "reactants_products.txt";
  open(RPO,">$f_recprod2") || die "Cannot open RPO at $f_recprod2/\n";

  for($index=1;$index<=$c_rxn;$index++){

    #Avoids dead due to no $f_recprod when no reactions were found (no genes in the group were matched to their ID).
    if($c_rxn<=1){
      #print"Avoided formatting RPD.\n";
      $index=3;
      next;
    }

    # Open IN file
    open(RPD,$f_recprod) || die "Cannot open RPD file for $TF because $index.\n";
    
    #Get objects
    @reactants=();
    @products=();
    while(<RPD>){
      if($_=~/^re$index\t([^\t]+)\treactant$/){
	push(@reactants,$1);
      }elsif($_=~/^re$index\t([^\t]+)\tproduct$/){
	push(@products,$1);
      }
    }
    close(RPD);
    
    ## Print reactants
    #Get words length
    undef %length;
    foreach $w (@reactants){
      @chars=split(//,$w);
      $size=@chars;
      $length{$w}=$size;
    }
    
    #Get longest and print in decreasing order
    @rp_word=();
    $longest="";
    $l_size=0;
    foreach $k (keys %length){
      $flag=0;
      foreach $r (@rp_word){
	if($r eq $k){
	  $flag=1;
	}
      }
      if(!($flag)){
	push(@rp_word,$k);
	if($length{$k} > $l_size){
	  $l_size=$length{$k};
	  $longest=$k;
	}
      }
    }
    
    for($i=$l_size;$i>0;$i--){
      foreach $k (keys %length){
	if($length{$k} == $i){
	  print RPO"re$index\t$k\treactant\n";
	}
      }
    }

    ## Print products
    #Get words length
    undef %length;
    foreach $w (@products){
      @chars=split(//,$w);
      $size=@chars;
      $length{$w}=$size;
    }
    
    #Get longest and print in decreasing order
    @rp_word=();
    $longest="";
    $l_size=0;
    foreach $k (keys %length){
      $flag=0;
      foreach $r (@rp_word){
	if($r eq $k){
	  $flag=1;
	}
      }
      if(!($flag)){
	push(@rp_word,$k);
	if($length{$k} > $l_size){
	  $l_size=$length{$k};
	  $longest=$k;
	}
      }
    }
    
    for($i=$l_size;$i>0;$i--){
      foreach $k (keys %length){
	if($length{$k} == $i){
	  print RPO"re$index\t$k\tproduct\n";
	}
      }
    }
  }


  #Print TFs (if not printed already)
  
  ##Get individual TFs
  @tfs=split(/_/,$TF);
  
  foreach $stf (@tfs){
    
    open(OBJ,$f_objects) || die "Cannot open OBJ at $f_objects.\n";
    $flag=0;
    while(<OBJ>){
      if($_=~/^$stf\t/){
	$flag=1;
      }
    }
    close(OBJ);
    if(!($flag)){
      open(OBJ,">>$f_objects") || die "Cannot open OBJ at $f_objects.\n";
      print OBJ"$stf\tPROTEIN\n";
      close(OBJ);
    }
  }

  #Delete duplicated objects
  system("cat $f_objects | sort | uniq > $uniq_objects");
  #Remove temporal files - except modifications, it will be appended in &regenz
  system("rm $f_objects $f_recprod $temp");

  #Calls function to add enzymatic regulation
  $renz=&regenz($d_out);
  if($renz){
    #print"Exited enzymatic regulation successfully.\n";
  }else{
    print"WARNING! Did not exit enzymatic regulation successfully.\n";
  }

  #Delete duplicated modifications
  system("cat $f_modification | sort | uniq > $uniq_mods");
  #Remove temporal files.
  system("rm $f_modification");

  return(1);

}


###########################
### enzyme Function
###
### Look for enzyme reactions and annotate them
### -Input arguments:
###  [0] - enzyme name
### -Returned values: array with data_temp lines
###########################

sub enzyme{

  local($enz)=($_[0]);

  @v_out=();
  @rxns=();

  @rxns=$cyc->reactions_of_enzyme($enz);
  foreach $rx (@rxns){
    @f_phyr=();
    $all_subs="";
    $all_prods="";
    @f_phyr=$cyc->get_slot_values($rx,"REACTION-DIRECTION");
    #If no reaction direction
    if(!(@f_phyr)){
      $f_phyr[0]="REVERSIBLE";
      print"Reaction $rx without direction. Considered reversible.\n";
    }
    $f_trans=0;
    $f_trans=$cyc->transport_rxn($rx);
    if($f_trans){
      if($f_phyr[0]=~/LEFT-TO-RIGHT/){
	@subs=();
	@subs=$cyc->get_slot_values($rx,"LEFT");
	foreach $s (@subs){
	  if(!($s=~/^PROTON$|^AMP$|^ATP$|^ADP$|^WATER$|^NAD.*$|^\|*Pi\|*$/i)){
	    $all_subs=$all_subs . "$s\_Ext,";
	  }
	}
	@prods=();
	@prods=$cyc->get_slot_values($rx,"RIGHT");
	foreach $p (@prods){
	  if(!($p=~/^PROTON$|^AMP$|^ATP$|^ADP$|^WATER$|^NAD.*$|^\|*Pi\|*$/i)){
	    $all_prods=$all_prods . "$p,";
	  }
	}
	#Complete and save line
	$line="$rx\tL2R\t$all_subs\t$all_prods";
	push(@v_out,$line);
      }elsif($f_phyr[0]=~/RIGHT-TO-LEFT/){
	@subs=();
	@subs=$cyc->get_slot_values($rx,"RIGHT");
	foreach $s (@subs){
	  if(!($s=~/^PROTON$|^AMP$|^ATP$|^ADP$|^WATER$|^NAD.*$|^\|*Pi\|*$/i)){
	    $all_subs=$all_subs . "$s\_Ext,";
	  }
	}
	@prods=();
	@prods=$cyc->get_slot_values($rx,"LEFT");
	foreach $p (@prods){
	  if(!($p=~/^PROTON$|^AMP$|^ATP$|^ADP$|^WATER$|^NAD.*$|^\|*Pi\|*$/i)){
	    $all_prods=$all_prods . "$p,";
	  }
	}
	$line="$rx\tL2R\t$all_subs\t$all_prods";
	push(@v_out,$line);
      }elsif($f_phyr[0]=~/REVERSIBLE/){
	@subs=();
	@subs=$cyc->get_slot_values($rx,"LEFT");
	foreach $s (@subs){
	  if(!($s=~/^PROTON$|^AMP$|^ATP$|^ADP$|^WATER$|^NAD.*$|^\|*Pi\|*$/i)){
	    $all_subs=$all_subs . "$s\_Ext,";
	  }
	}
	@prods=();
	@prods=$cyc->get_slot_values($rx,"RIGHT");
	foreach $p (@prods){
	  if(!($p=~/^PROTON$|^AMP$|^ATP$|^ADP$|^WATER$|^NAD.*$|^\|*Pi\|*$/i)){
	    $all_prods=$all_prods . "$p,";
	  }
	}
	$line="$rx\tRVB\t$all_subs\t$all_prods";
	push(@v_out,$line);
      }
    }else{
      if($f_phyr[0]=~/LEFT-TO-RIGHT/){
	@subs=();
	@subs=$cyc->get_slot_values($rx,"LEFT");
	foreach $s (@subs){
	  if(!($s=~/^PROTON$|^AMP$|^ATP$|^ADP$|^WATER$|^NAD.*$/i)){
	    $all_subs=$all_subs . "$s,";
	  }
	}
	@prods=();
	@prods=$cyc->get_slot_values($rx,"RIGHT");
	foreach $p (@prods){
	  if(!($p=~/^PROTON$|^AMP$|^ATP$|^ADP$|^WATER$|^NAD.*$/i)){
	    $all_prods=$all_prods . "$p,";
	  }
	}
	$line="$rx\tL2R\t$all_subs\t$all_prods";
	push(@v_out,$line);
      }elsif($f_phyr[0]=~/RIGHT-TO-LEFT/){
	@subs=();
	@subs=$cyc->get_slot_values($rx,"RIGHT");
	foreach $s (@subs){
	  if(!($s=~/^PROTON$|^AMP$|^ATP$|^ADP$|^WATER$|^NAD.*$/i)){
	    $all_subs=$all_subs . "$s,";
	  }
	}
	@prods=();
	@prods=$cyc->get_slot_values($rx,"LEFT");
	foreach $p (@prods){
	  if(!($p=~/^PROTON$|^AMP$|^ATP$|^ADP$|^WATER$|^NAD.*$/i)){
	    $all_prods=$all_prods . "$p,";
	  }
	}
	$line="$rx\tL2R\t$all_subs\t$all_prods";
	push(@v_out,$line);
      }elsif($f_phyr[0]=~/REVERSIBLE/){
	@subs=();
	@subs=$cyc->get_slot_values($rx,"LEFT");
	foreach $s (@subs){
	  if(!($s=~/^PROTON$|^AMP$|^ATP$|^ADP$|^WATER$|^NAD.*$/i)){
	    $all_subs=$all_subs . "$s,";
	  }
	}
	@prods=();
	@prods=$cyc->get_slot_values($rx,"RIGHT");
	foreach $p (@prods){
	  if(!($p=~/^PROTON$|^AMP$|^ATP$|^ADP$|^WATER$|^NAD.*$/i)){
	    $all_prods=$all_prods . "$p,";
	  }
	}
	$line="$rx\tRVB\t$all_subs\t$all_prods";
	push(@v_out,$line);
      }
    }
  }

  return(\@v_out);
}

###########################
### regenz Function
###
### Finds out if any molecules in GU regulates an enzyme also present in the GU
### -Input arguments: none
### -Returned values: none - directly modifies modification.txt and prints a message on screen
###########################

sub regenz{
  local($d_out)=$_[0];
  local($f_rxn,$f_obj,$molid,@regs,$rg,$ent,@enzymes,@enties,@rxns,$enz,@effect,$qenz,$f_mod,$ename,$qename,$mrx,$molnm,$rx,$qrx)=("","","",(),"","",(),(),(),"",(),"","","","","","","","");

  #Open OBJ and RXN files. 
  $f_rxn=$d_out . "reactions.txt";
  $f_obj=$d_out . "objects.txt";
  $f_mod=$d_out . "t_modification.txt";

  #Find simple molecules present in the GU (omits external forms).
  open(OBJ,$f_obj) || die "Cannot open OBJ at $f_obj.\n";
  while(<OBJ>){
    if($_=~/^([^\t]+)\tSIMPLE_MOLECULE\t([^\t]+)/){
      $molid=$2;
      $molnm=$1;
      chomp($molid);
      if(($molnm=~/^phosphate$/) || ($molnm=~/^diphosphate$/) || ($molnm=~/^.+_Ext$/)){
	next;
      }
      #print"==$molnm == $molid\n";
      #Find proteins regulated by molecule
      @regs=$cyc->get_slot_values($molid,"REGULATES"); # Gets regulatory interactions ids
      foreach $rg (@regs){
	#print"$molid - $rg\n";
	@enties=$cyc->get_slot_values($rg,"REGULATED-ENTITY"); # Gets regulates entity. Name does not always coincides with IDs in OBJ.
	foreach $ent (@enties){
	  #print"$molid - $rg - $ent\n";
	  #@enzymes=$cyc->get_slot_values($ent,"ENZYME"); # Gets more common name of the regulated enzyme - to use, activate proceding code.
	  @rxns=$cyc->get_slot_values($ent,"REACTION"); # Gets the ID of the regulated reactions.
	  
	  ##
	  ## The following chunk annotates all reactions in the GU that are catalyzed by an enzyme regulated by $mol as regulated. If used, also omit comment in the last two lines
	  ##
	  #foreach $enz (@enzymes){ # Look for regulated enzymes
	  #  print"$molid - $rg - $ent - enzyme=$enz\n";
	  #  $qenz=quotemeta($enz);
	  #  open(OBZ,$f_obj) || die "Cannot open OBZ at $f_obj.\n";
	  #  while(<OBZ>){
	  #    if($_=~/^([^\t]+)\t[^\t]+\t$qenz$/){
	  #	$ename=$1;
	  #	$qename=quotemeta($ename);
		# Print modification in MOD file.
		# Get modification effect
	  #	@effect=$cyc->get_slot_values($rg,"MODE"); # Gets more common name of the regulated enzyme
		# Find modified reactions
	  #	open(RMD,$f_mod) || die "Cannot open RMD at $f_mod.\n";
	  #	while(<RMD>){
	  #	  if($_=~/^(re\d+)\tCATALYSIS\t$qename$/){
	  #	    $mrx=$1;
		    # Open MOD file to append
	  #	    open(MOX,">>$f_mod") || die "Cannot open MOX at $f_mod.\n";
	  #	    if($effect[0]=="-"){ # REVISAR LAS OPCIONES
	  #	      print MOX"$mrx\tINHIBITION\t$molnm\n";
	  #	    }elsif($effect[0]=="+"){
	  #	      print MOX"$mrx\tPHYSICAL_STIMULATIONt$molnm\n";
	  #	    }
	  #	    close(MOX);
	  #	  }
	  #	}
	  #    }
	   # }
	  #  close(OBZ);
	  #}

	  #Look for regulated reactions
	  foreach $rx (@rxns){
	    #print"$molid - $rg - $ent - reaction=$rx\n";
	    $qrx=quotemeta($rx);
	    open(RXR,$f_rxn) || die "Cannot open RXR at $f_rxn.\n";
	    while(<RXR>){
	      if($_=~/^(re\d+)\t[^\t]+\t[^\t]*\t$qrx$/){
		$mrx=$1;
		@effect=$cyc->get_slot_values($rg,"MODE"); # Gets more common name of the regulated enzyme
		#Open MOD file to append
		open(MOX,">>$f_mod") || die "Cannot open MOX at $f_mod.\n";
		if($effect[0]=="-"){ # Only options are - and +
		  print MOX"$mrx\tINHIBITION\t$molnm\t$rg\n";
		}elsif($effect[0]=="+"){
		  print MOX"$mrx\tPHYSICAL_STIMULATIONt$molnm\t$rg\n";
		}
		close(MOX);
	      }
	    }
	    close(RXR);
	  }
	}
      }
    }
  }
  close(OBJ);


  return(1);
}
		  

###########################
### IDS Function
###
### Turn Ecocyc id's into common names and vv.
### -Input arguments:
###  [0] - id/common name
###  [1] - type of entity (g=gene, p=protein, pn=protein to convert to gene, cx=complex, sm=simple molecule).
### -Returned values: the id or common name, opposite of input.
###########################

sub IDS {
 
  local($name,$type)=($_[0],$_[1]);

  ##Genes to ecocyc ID.
  if($type eq "g"){
    local($id);
    open(IDS, $gene_links) || die "Cannot open IDS at $gene_links.\n";
    while(<IDS>){
      if($_=~/^([^\t]+)\t[^\t]*\t[^\t]*\t[^\t]*\t[^\t]*\t[^\t]*\t$name$/i){
	$id=$1;
	close(IDS);
	return($id);
      }
    }
    close(IDS);
    open(IDS, $gene_links) || die "Cannot open IDS at $gene_links.\n";
    # Check if gene names are in ecocyc name (third column of gene names)
    while(<IDS>){
      if($_=~/^([^\t]+)\t[^\t]*\t[^\t]*\t$name\t/i){
	$id=$1;
	close(IDS);
	return($id);
      }
    }
    close(IDS);
    print"WARNING 2! Gene $name is not present in gene_links.\n";
    return(0);
  }

  ##Gene name from Ecocyc ID
  if($type eq "gi"){
    local($nm);
    open(IDS, $gene_links) || die "Cannot open IDS at $gene_links.\n";
    while(<IDS>){
      if($_=~/^$name\t[^\t]*\t[^\t]*\t[^\t]*\t[^\t]*\t[^\t]*\t([^\t]+)$/i){
	$nm=$1;
	chomp($nm);
	close(IDS);
	return($nm);
      }
    }
    close(IDS);
    print"Warning 6! $name was not found at $gene_links in &IDS type \"gi\".\n";
    return($name);
  }

  ##Proteins
  if($type eq "p"){

    #RNA
    if(($name=~/RNA/) && (!($name=~/MONOMER/))){
      #Declare local variables
      local($rna_id);     
      #If class name (between pipes)
      if($name=~/^\|(.+)\|$/){
	$name=$1;
	return($name);
      }else{
	#Get synonyms
	$rna_id=$cyc->get_slot_value($name, "SYNONYMS");
	#If no synonyms
	if(!($rna_id)){
	  #Get names
	  $rna_id=$cyc->get_slot_value($name, "NAMES");
	  #If no names
	  if(!($rna_id)){
	    return("RNA");
	    #If names
	  }else{
	    $rna_id=&xml($rna_id);
	    return("RNA:$rna_id");
	  }
	  #If synonyms
	}else{
	  $rna_id=&xml($rna_id);
	  return("RNA:$rna_id");
	}
      }
    }

    ###Proteins
    #Declare local variables
    local(@gn_id,$sz_id,@p_id,$pgi,$gn_nm,$p_id);
    #Get gene to find similar name (if only 1 gene)
    @gn_id=$cyc->genes_of_protein($name);
    $sz_id=@gn_id;
    if($sz_id == 1){
      $gn_nm=&IDS($gn_id[0],"gi");
      @p_id=$cyc->get_slot_values($name,"NAMES");
      foreach $pgi (@p_id){
	if($pgi=~/^$gn_nm$/i){
	  return($pgi);
	}
      }
    }
    $p_id=$cyc->get_slot_value($name,"SYNONYMS");
    $p_id=&xml($p_id);
    #If no synonyms
    if((!($p_id=~/^\w{3,7}$/)) || ($p_id=~/^B\d+/)){
      $p_id=&IDS($name,"pn");
      return($p_id);
    }else{
      return($p_id);
    }

  }

  ##Heteromeric Complexes
  if($type eq "cx"){
    #Declare local variables
    local(@mons,$m,$c_name)=((),"","");
    @mons=$cyc->monomers_of_protein($name);
    foreach $m (@mons){
      $m_nm=&IDS($m,"p");
      $c_name=$c_name . "-" . $m_nm;
    }
    if($c_name=~/^-(.+)/){
      $c_name=$1;
      return($c_name);
    }else{
      print"Warning 7! No complex name for $name.\n";
      return"NP";
    }
  }

  ##Small Molecules from Ecocyc IDs
  if($type eq "sm"){

    #Declare local variables
    local($q_name);
    $q_name=quotemeta($name);
    open(SMS,$dictionary) || die "Cannot open SMS file at $dictionary.\n";
    while(<SMS>){
      if($_=~/^$q_name\t([^\t]+)\t/i){
	$name=$1;
	$name=&xml($name);
	return($name);
      }
    }
    close(SMS);
    if($name=~/^\|(.+)\|$/){
      return($1);
    }
    print"Warning 8! Small Molecule $name not found in dictionary.\n";
    return($name);
  }

  ##Small Molecules from Regulon DB
  if($type eq "rDB"){
    #Declare local variables
    local($q_name);
    $q_name=quotemeta($name);
    open(SMS,$dictionary) || die "Cannot open SMS file at $dictionary.\n";
    while(<SMS>){
      if($_=~/^[^\t]+\t([^\t]+)\t.*\/\/$q_name\/\//i){
	$name=$1;
	$name=&xml($name);
	return($name);
      }
    }
    close(SMS);
    print"WARNING 3! $name from RegulonDB not found in dictionary.\n";
    $name=&xml($name);
    return($name);
  }

  ##Proteins short name (gene in form XxxX)
  if($type eq "pn"){
    #Declare local variables
    local($g_pt,@g_pt,$g_nm);
    @g_pt=$cyc->genes_of_protein($name);
    open(IDS, $gene_links) || die "Cannot open IDS at $gene_links.\n";
    while(<IDS>){
      if($_=~/^$g_pt[0]\t[^\t]*\t[^\t]*\t[^\t]*\t[^\t]*\t[^\t]*\t([^\t]+)/i){
	$g_nm=$1;
	$g_nm=ucfirst($g_nm);
	chomp($g_nm);
	return($g_nm);
      }
    }
    close(IDS);
    print"Warning 4! No gene name for $g_pt in $name.\n";
    return("NP");
  }
      

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
  while($line=~/(.*)&omega;(.*)/i){
    $line=$1 . "omega" . $2;
  }
  return($line);
  
}

###########################
### is_protein function
###
### Looks for a protein name in protein-links. Returns 0 or 1 if absent/present
### -Input arguments: protein name
### -Returned values: line in a single char value
###########################

sub is_protein {

  local($name,$q_name)=($_[0],"");

  open(PRT,$gene_links) || die "Cannot open PRT in file $gene_links at &is_protein\n";
  $q_name=quotemeta($name);
  while(<PRT>){
    if(($_=~/\t$q_name\t/i) || ($_=~/\t$q_name$/i)) {
      close(PRT);
      return(1);
    }
  }
  close(PRT);
  return(0);
}
