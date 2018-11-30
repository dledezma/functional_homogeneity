#################################################################################
###### 29-01-18                                                            ######
###### Connectivity for big datasets with enzymatic regulation             ######
###### Calculates connectivity omitting the search for short paths.        ######
###### Finds components from list of reactions through an adyacence matrix ######
###### Considers enzymes as "connected" when a molecule (in the GU)modifies #####
###### a rxn they catalyze                                                 ######
###### [0] - gu_assembly_vX.pl output directory                            ######
###### [1] - Output file                                                   ######
###### Output: Creates a file for each GU with all paths                   ######
###### [1]> GU name                                                        ######
###### [2]> Total enzymes                                                  ######
###### [2]> Total connected enzymes                                        ######
###### [3]> Connected enzymes                                              ######
###### [4]> Total components                                               ######
###### [5]> Components                                                     ######
###### [6]> Total Enzymatic reactions                                      ######
###### [7]> Connectivity                                                   ######
###### History:                                                            ######
###### Notes:                                                              ######
###### - Considers only reactant/product relationships in connected        ######
######   enzymes calculation but reactant/reactant & prod/prod for         ######
######   components. Does so by calculating adyacence matrix twice.        ######
######   Prints connectivity with only 1 decimal.                          ######
###### - Prints count of enzymatic reactions and transport reactions       ######
######  separately for filtering of GUs that only include transport        ######
#################################################################################

$time=localtime();
print"Begin: $time\n";

#Get args
$dir_path=$ARGV[0];
$out_file=$ARGV[1];
chomp($out_file);


#Get temp file name (for reversible reactions)
if($out_file=~/^(.+)\.txt/){
  $temp=$1 . "_temp.txt";
}else{
  die "Cannot create temp file, outfile is not .txt\n";
}

#Open outfile
open(OUT,">$out_file") || die "Cannot open OUT at $out_file.\n";
print OUT"# From $0 on $time\n# Input = $dir_path\n# GU name\tTotal enzymes\tEnzymes\tTotal connected enzymes\tConnected enzymes\tTotal components\tComponents\tTotal Enzymatic reactions\tTotal transport reactions\tConnectivity\n";

#Open dir 
opendir(DIR,$dir_path) || die "Cannot open DIR at $in_dir.\n";
@GU_files=readdir(DIR);
foreach $file (sort @GU_files){
  if(!($file=~/^\./)){

    print"\n**$file\n";
    $TF=$file;
    
    $time=localtime();
    print"Started: $time\n";

    #Get files paths
    $products=$dir_path . $file . "/reactants_products.txt";
    $objects=$dir_path . $file . "/objects.txt";
    $reactions=$dir_path . $file . "/reactions.txt";
    $modifications=$dir_path . $file . "/modification.txt";


    #### Get Total Catalysis Reactions ####

    #Reset variables
    $cat_rxns=0;
    @catrxns=();
    $tpt_rxns=0;
    @tptrxns=();

    #Open reactions files
    open(RXN,$reactions) || die "Cannot open RXN file at $reactions\n";
    while(<RXN>){
      if($_=~/^(re\d+)\tSTATE_TRANSITION\t[^\t]+\t/){  # Only enzymatic reactions are considered since comples formation will not have a tab after direction column.
	$rxn=$1;
	$cat_rxns++;
	push(@catrxns,$rxn);
      }elsif($_=~/^(re\d+)\tTRANSPORT\t[^\t]+\t/){
	$rxn=$1;
	$cat_rxns++;
	$tpt_rxns++;
	push(@tptrxns,$rxn);
	push(@catrxns,$rxn);
      }
    }
    close(RXN);

    #print"$cat_rxns catalysis reactions:";
    foreach $rxn (@catrxns){
      #print"$rxn,";
    }
    #print"\n";

    ########

    #### Get Total Enzymes ####

    #Initialize variables
    %enzymes=();

    #Find CATALYSIS modifications
    open(MOD,$modifications) || die "Cannot open MOD at $modifications.\n"; 
    while(<MOD>){
      if($_=~/^re\d+\tCATALYSIS\t([^\t]+)$/i){
	     $prot=$1;
	     chomp($prot);
	     $enzymes{$prot}++;
      }
    }
    close(MOD);

    @enzymes=keys %enzymes;
    $c_enzymes=@enzymes;
    #print"$c_enzymes enzymes.\n";

    ########

    #### Create temporal file with reversible reactions ####
    
    #Create temp file
    open(TMP,">$temp") || die "Cannot open TMP file at TMP.\n";

    #Create temp reactants_products file with duplicated entries for REV reactions.
    open(RXN,$reactions) || die "Cannot open RXN at $reactions.\n";
    while(<RXN>){
      if(($_=~/^(re\d+)\tSTATE_TRANSITION\tRVB\t/i) || ($_=~/^(re\d+)\tTRANSPORT\tRVB\t/i)){
	     $rex=$1;
	     open(RPD,$products) || die "Cannot open RPD at $products.\n";
	     while(<RPD>){
	       if($_=~/^$rex\t([^\t]+)\treactant$/i){
	         print TMP"$rex\t$1\treactant\n$rex\t$1\tproduct\n";
	       }
	       if($_=~/^$rex\t([^\t]+)\tproduct$/i){
	         print TMP"$rex\t$1\treactant\n$rex\t$1\tproduct\n";
	       }
	     }
	     close(RPD);
      }elsif(($_=~/^(re\d+)\tSTATE_TRANSITION\tL2R\t/) || ($_=~/^(re\d+)\tTRANSPORT\tL2R\t/)){
	     $rex=$1;
	     open(RPD,$products) || die "Cannot open RPD at $products.\n";
	     while(<RPD>){
	       if($_=~/^$rex\t[^\t]+\t[^\t]+$/i){
	         print TMP"$_";
	       }
	     }
	     close(RPD);
      }elsif(($_=~/^(re\d+)\tSUPER_REACTION\t/) || ($_=~/^(re\d+)\tSTATE_TRANSITION$/)){   #If reaction is 2ndry reaction or complex formation. Complex formation are needed to consider in connectivity single reactions that produce the effector. These contribute to connectivity but are omitted as heads of search for components (are only on one side of the matrix).
	     $rex=$1;
	     open(RPD,$products) || die "Cannot open RPD at $products.\n";
	     while(<RPD>){
	       if($_=~/^$rex\t[^\t]+\t[^\t]+$/i){
	         print TMP"$_";
	       }
	     }
	     close(RPD);
      }	
    }
    close(RXN);
    close(TMP);

    ########                                                                 ########
    #### Create adyacence matrix for Connected Enzymes                           ####
    #### Only considers reactant/product connections and enzymatic modifications ####
    #### An enzyme is connected whenever a catalyzed reaction is modified by a   ####
    #### molecule in the GU                                                      ####
    ########                                                                 ########

    #Reset matrixes
    %adyacence=();
    %components=(); #Will be curcial in next chunk of code (find components) but is declared here to start including rxns that are totally disconnected and therefore will not be in %adyacence, i.e. those for which $f_connected flag remains 0.
    %used=();  #Same as %components. Is used for independent reactions to mantain uniformity in results.
    $matrix_size=@catrxns;
    for($im=0;$im<$matrix_size;$im++){   #Only begins with catalysis reactions although supe reactions and complex formation reactions will be in TMP file and will be considered. Matrix will not be diagonal, but only half will be filled to optimize memory usage. Only 1's are being added instead of going through all the combinations of coordinates.
      $f_connected=0; #Flag to signal that a reaction in @catrxns is not connected to any other reaction in TMP, will not be included in %adyacence and should be cnsidered as an independent component.
      open(TMP,$temp) || die "Cannot open TMP file at $temp.\n";
      while(<TMP>){
	if($_=~/^$catrxns[$im]\t([^\t]+)\t/){ 
	  $reactant=$1;
	  $qreactant=quotemeta($reactant);
	  if($qreactant=~/^phosphate$/){
	    next;
	  }

	  # See if molecule is a modifier
	  open(MDX,$modifications) || die "Cannot open MDX at $modifications.\n";
	  while(<MDX>){
	    if($_=~/^(re\d+)\t[^\t]+\t$qreactant$/){ # Enzymatic modifications include the Ecocyc ID in 4th column.
	      $nrx=$1;
	      if($catrxns[$im] ne $nrx){
		$adyacence{$catrxns[$im]}{$nrx}=2; #Second index will be longer because it will include 2dry reactions and complex formation reactions.
		$f_connected=1;
	      }
	    }
	  }
	  close(MDX);

	  ### Look for other reactions that include de molecule irrespective of role.
	  open(TMP1,$temp) || die "Cannot open TMP1 file at $temp.\n";
	  while(<TMP1>){
	    if($_=~/^(re[^\t]+)\t$qreactant\t/){
	      $nrx=$1;
	      #print"$catrxns[$im] / $nrx / r-p\n";
	      if($catrxns[$im] ne $nrx){
		$adyacence{$catrxns[$im]}{$nrx}=1; #Second index will be longer because it will include 2dry reactions and complex formation reactions.
		$f_connected=1;
	      }
	    }
	  }
	}
      }
      close(TMP);

    }

    #Print adyacence matrix                                                           ##### Uncomment to print adyacence matrix
    #print"Matrix for connected enzymes:\n";
    #foreach $dn (sort keys %adyacence){
    #  print"$dn\t";
    #  foreach $nx (sort keys %{$adyacence{$dn}}){
    #	     print"$nx-$adyacence{$dn}{$nx}\t";
    #  }
    #  print"\n";
    #}

    ########
    
    #### Find components of connected enzymes adyacence matrix - NOT final component calculation####
    #### Only considers reactant/product connections ####

    foreach $head (sort keys %adyacence){
      #print"Head=$head\n";
      
      if($used{$head}){   #Omit reactions that are already in a component to optimize process.
	next;
      }
      @members=();
      push(@members,$head);
      $ref=&component_members(\@members,$head);
      @members=@{$ref};
      $new_component=1;
      foreach $mb (@members){
	if($used{$mb}){ #If a member of the new component (@members) is already in another component, they will be fused.
	  #Insert new members to @{$components{$head}} avoiding repetitions.
	  foreach $mb1 (@members){
	    $rep_flag=0;
	    foreach $mb2 (@{$components{$used{$mb}}}){  #$used{$mb} holds the key in %components that includes the component to which @members will be fused.
	      if($mb1 eq $mb2){
		$rep_flag=1;  #If new member is already present in %components, change flag.
		last;
	      }
	    }
	    if(!($rep_flag)){ #If new member is not yet present in %components, add it.
	      push(@{$components{$used{$mb}}},$mb1);
	      $used{$mb1}=$used{$mb}; # Save %components key in which the rxn resides and mark it as used at tha same time. It is inside "if" to avoid reassigning the key to members that were already present.
	    }
	  }
	  $new_component=0; #Marks that @members is not a new component and was fused to an already existing %components key.
	  last; #Breaks "foreach $mb (@members)". If one member was found used, it is enough to stop the search for another already used member. It will cause some components to overlap and not be fused. A last step of looking for overlaps once all components are in %components was added below to overcome this issue.
	}
      }

      if($new_component){ #If no member of @members is already part of a component in %components.
	push(@{$components{$head}},@members);  #New component key will be $head.
	foreach $mb (@members){
	  $used{$mb}=$head; # Save component name in which the rxn resides and mark it as used at tha same time.
	}
      }
    }

    #print components                                                                              ######### Uncomment to print components
    #foreach $cm (sort keys %components){
    #  print"=$cm\t";
    #  foreach $rx (@{$components{$cm}}){
    #print"$rx,";
    #  }
    #  print"\n";
    #}

    #Find final overlaps between individual components in %components.
    @heads=sort keys %components;
    $size_heads=@heads;
    for($in=0;$in<$size_heads;$in++){
      $head1=$heads[$in];
      for($gn=($in+1);$gn<$size_heads;$gn++){
	$head2=$heads[$gn];
	$f_overlap=0;
	foreach $mb1 (@{$components{$head1}}){
	  foreach $mb2 (@{$components{$head2}}){
	    if($mb1 eq $mb2){  #If a member is in both components.
	      #print"$mb1 overlap with $mb2, components of $head2 will be fused to $head1.\n";
	      $f_overlap=1;  #Turn on flag of overlap
	      last;  #break foreach
	    }
	  }
	  if($f_overlap){   #Break second foreach. One common member is enough to fuse components.
	    last;
	  }
	}
	if($f_overlap){
	  foreach $mba (@{$components{$head2}}){ #Component from $head2 will be added to $head1
	    $f_new=1;
	    foreach $mbb (@{$components{$head1}}){ #Check if member of @{$components{$head2}} is already in @{$components{$head1}}
	      if($mba eq $mbb){
		$f_new=0;
		last;
	      }
	    }
	    if($f_new){ #If @{$components{$head2}} member is not already in @{$components{$head1}}
	      push(@{$components{$head1}},$mba); #Add member from @{$components{$head2}}
	    }
	  }

	  #Delete $head2 from %components, since it is already in $head1 key. No array will be affected since for loops are on @heads. Note that it is impossible that $gn is a position lower than $in, which would mess up the indexes. When $in moves forward, it will find no problems.
	  delete $components{$head2};
	  #Delete $head2 from @heads to avoid spurious comparisons.
	  splice(@heads,$gn,1);
	  $size_heads--;  #Reduce @heads size because $head2 was deleted.
	  $gn=$in;  #Check all other heads again in case new elements find new intersections, $gn will be $in+1 when next iteration begins because of $gn++. Groups previous to $in+1 are not needed to be considered because if an overlap existed it would have been found when previous group was $in. 
	}
      }
    }

    #print components                                                           ######## Uncomment to print new components
    #print"New components:\n";
    #foreach $cm (sort keys %components){
    #  print"=$cm\t";
    #  foreach $rx (@{$components{$cm}}){#
	#print"$rx,";
    #  }
    #  print"\n";
    #}

    ########

    #### Get Total Connected Enzymes ####

    #Obtain connected enzymes by omitting components of 1 rxn and looking for catalyzers.
    %conn_enz=(); #Needed for TFs with no components.
    foreach $head (keys %components){
      $size=(@{$components{$head}});
      if($size==1){  #Omit components of 1 reaction since they are not connected reactions.
	next;
      }
      #%conn_enz=();
      foreach $rxn (@{$components{$head}}){
	open(MOD,$modifications) || die "Cannot open MOD at $modifications.\n";
	while(<MOD>){
	  if($_=~/^$rxn\tCATALYSIS\t([^\t]+)/){
	    $enzyme=$1;
	    chomp($enzyme);
	    $conn_enz{$enzyme}++;
	  }
	}
	close(MOD);
      }
    }
      
    #Save total connected enzymes
    @heads=keys %conn_enz;
    $total_conn_enz=@heads;
    #print"$total_conn_enz connected enzymes\n";

    ########
    
    #### Create adyacence matrix for components ####
    #### Considers all connections between reactions: reactant/reactant, product/product, reactant/product ####

    #Reset matrixes
    %adyacence=();
    %components=(); #Will be curcial in next chunk of code (find components) but is declared here to start including rxns that are totally disconnected and therefore will not be in %adyacence, i.e. those for which $f_connected flag remains 0.
    %used=();  #Same as %components. Is used for independent reactions to mantain uniformity in results.

    $matrix_size=@catrxns;
    for($im=0;$im<$matrix_size;$im++){   #Only begins with catalysis reactions although supe reactions and complex formation reactions will be in TMP file and will be considered. Matrix will not be diagonal, but only half will be filled to optimize memory usage. Only 1's are being added instead of going through all the combinations of coordinates.
      $f_connected=0; #Flag to signal that a reaction in @catrxns is not connected to any other reaction in TMP, will not be included in %adyacence and should be cnsidered as an independent component.
      open(TMP,$temp) || die "Cannot open TMP file at $temp.\n";
      while(<TMP>){
	if($_=~/^$catrxns[$im]\t([^\t]+)\t/){ 
	  $reactant=$1;
	  $qreactant=quotemeta($reactant);
	  if($qreactant=~/^phosphate$/){
	    next;
	  }
	  open(TMP1,$temp) || die "Cannot open TMP1 file at $temp.\n";
	  while(<TMP1>){
	    if($_=~/^(re[^\t]+)\t$qreactant\t/){
	      $nrx=$1;
	      #print"$catrxns[$im] / $nrx / r-p\n";
	      if($catrxns[$im] ne $nrx){
		$adyacence{$catrxns[$im]}{$nrx}=1; #Second index will be longer because it will include 2dry reactions and complex formation reactions.
		$f_connected=1;
	      }
	    }
	  }
	}
      }
      close(TMP);
      
      #Save reactions that have already been identified as individual components since they are not connected with any other reaction in TMP. They are no included in %adyacence to avoid increasing comparisons and execution time. They will only be printed with components.
      if(!($f_connected)){  
	push(@{$components{$catrxns[$im]}},$catrxns[$im]);
	$used{$catrxns[$im]}=$catrxns[$im];
      }

    }

    #Print adyacence matrix                                                       ########### Uncomment to print matrix
    #print"Matrix for components:\n";
    #foreach $dn (sort keys %adyacence){
    #  print"$dn\t";
    #  foreach $nx (sort keys %{$adyacence{$dn}}){
#	     print"$nx-$adyacence{$dn}{$nx}\t";
    #  }
    #  print"\n";
    #}
		

    ########
    
    #### Find components ####	  
    foreach $head (sort keys %adyacence){
      #print"Head=$head\n";
      
      if($used{$head}){   #Omit reactions that are already in a component to optimize process.
	next;
      }
      @members=();
      push(@members,$head);
      $ref=&component_members(\@members,$head);
      @members=@{$ref};
      $new_component=1;
      foreach $mb (@members){
	if($used{$mb}){ #If a member of the new component (@members) is already in another component, they will be fused.
	  #Insert new members to @{$components{$head}} avoiding repetitions.
	  foreach $mb1 (@members){
	    $rep_flag=0;
	    foreach $mb2 (@{$components{$used{$mb}}}){  #$used{$mb} holds the key in %components that includes the component to which @members will be fused.
	      if($mb1 eq $mb2){
		$rep_flag=1;  #If new member is already present in %components, change flag.
		last;
	      }
	    }
	    if(!($rep_flag)){ #If new member is not yet present in %components, add it.
	      push(@{$components{$used{$mb}}},$mb1);
	      $used{$mb1}=$used{$mb}; # Save %components key in which the rxn resides and mark it as used at tha same time. It is inside "if" to avoid reassigning the key to members that were already present.
	    }
	  }
	  $new_component=0; #Marks that @members is not a new component and was fused to an already existing %components key.
	  last; #Breaks "foreach $mb (@members)". If one member was found used, it is enough to stop the search for another already used member. It will cause some components to overlap and not be fused. A last step of looking for overlaps once all components are in %components was added below to overcome this issue.
	}
      }

      if($new_component){ #If no member of @members is already part of a component in %components.
	push(@{$components{$head}},@members);  #New component key will be $head.
	foreach $mb (@members){
	  $used{$mb}=$head; # Save component name in which the rxn resides and mark it as used at tha same time.
	}
      }
    }

    #print components                                                            ####### Uncomment to print components
    #foreach $cm (sort keys %components){
    #  print"=$cm\t";
    #  foreach $rx (@{$components{$cm}}){
#	     print"$rx,";
    #  }
    #  print"\n";
    #}

    #Find final overlaps between individual components in %components.
    @heads=sort keys %components;
    $size_heads=@heads;
    for($in=0;$in<$size_heads;$in++){
      $head1=$heads[$in];
      for($gn=($in+1);$gn<$size_heads;$gn++){
	$head2=$heads[$gn];
	$f_overlap=0;
	foreach $mb1 (@{$components{$head1}}){
	  foreach $mb2 (@{$components{$head2}}){
	    if($mb1 eq $mb2){  #If a member is in both components.
	      print"$mb1 overlap with $mb2, components of $head2 will be fused to $head1.\n";
	      $f_overlap=1;  #Turn on flag of overlap
	      last;  #break foreach
	    }
	  }
	  if($f_overlap){   #Break second foreach. One common member is enough to fuse components.
	    last;
	  }
	}
	if($f_overlap){
	  foreach $mba (@{$components{$head2}}){ #Component from $head2 will be added to $head1
	    $f_new=1;
	    foreach $mbb (@{$components{$head1}}){ #Check if member of @{$components{$head2}} is already in @{$components{$head1}}
	      if($mba eq $mbb){
		$f_new=0;
		last;
	      }
	    }
	    if($f_new){ #If @{$components{$head2}} member is not already in @{$components{$head1}}
	      push(@{$components{$head1}},$mba); #Add member from @{$components{$head2}}
	    }
	  }

	  #Delete $head2 from %components, since it is already in $head1 key. No array will be affected since for loops are on @heads. Note that it is impossible that $gn is a position lower than $in, which would mess up the indexes. When $in moves forward, it will find no problems.
	  delete $components{$head2};
	  #Delete $head2 from @heads to avoid spurious comparisons.
	  splice(@heads,$gn,1);
	  $size_heads--;  #Reduce @heads size because $head2 was deleted.
	  $gn=$in;  #Check all other heads again in case new elements find new intersections, $gn will be $in+1 when next iteration begins because of $gn++. Groups previous to $in+1 are not needed to be considered because if an overlap existed it would have been found when previous group was $in. 
	}
      }
    }

    #print components                                                               ###### Uncomment to print new components
    #print"New components:\n";
    #foreach $cm (sort keys %components){
    #  print"=$cm\t";
    #  foreach $rx (@{$components{$cm}}){
#	print"$rx,";
    #  }
    #  print"\n";
    #}

    #Save total number of componentes (path length = 0)
    @heads=keys %components;
    $total_components=@heads;



    ######################################################
    #### Use to print components < limit of reactions ####
    ######################################################
    #if($limit){  #If limit=0, entirely omit this step.
    #  @heads=sort keys %components;
    #  $size_heads=@heads;
    #  for($in=0;$in<$size_heads;$in++){
    #	$head1=$heads[$in];
    #	$path_len=(@{$components{$head1}});
    #	$path_len--; #Path length= reactions in component - 1. E.g. a path length of 1 involves 2 reactions
    #	if($path_len<$limit){
    #	  delete $components{$head1};
    #	}
    #	splice(@heads,$in,1);
    #	$size_heads--;
    #	$in--;
    #  }
    #}
    ######################################################

    ########

    #### Calculate Connectivity ####

    #Calculate connectivity = connected enzymes/(total enzymes - (components-1))
    $cal_components=$total_components - 1;
    $denominator=$c_enzymes + $cal_components;
    if($denominator){ 		# Avoid illegal division by zero.
      $connectivity=$total_conn_enz / $denominator;
    }else{
      $connectivity=0;
    }
    $rounded = sprintf "%.1f", $connectivity;

    ########

    #### Print Output ####    
    print OUT"$TF\t$c_enzymes\t";
    foreach $mb (keys %enzymes){
      print OUT"$mb,";
    }
    print OUT"\t$total_conn_enz\t";
    foreach $mb (keys %conn_enz){
      print OUT"$mb,";
    }
    print OUT"\t$total_components\t";
    foreach $mb (keys %components){
      print OUT"*";
      foreach $el (@{$components{$mb}}){
	print OUT"$el,";
      }
    }
    print OUT"*\t$cat_rxns\t$tpt_rxns\t$rounded\n";

    ########


  }
}

#
system("rm $temp");


###############################
### subroutine component_members
### Finds all members of a component by recursively identifying connections among reactions
### input: 
### [0] > array with component members
### [1] > newest member
### output: array with all component members
### NOTE: SRs and complex formation reactions will cause function to return since they are not in first index of %adyacence
  # New member is already in array    
###############################

sub component_members {

  local(@members)=(@{$_[0]});
  local($new,$ky,$ref,$mb,$cycle)=($_[1],"","","",0);
  #print"Evaluating $new\n";

  foreach $ky (sort keys %{$adyacence{$new}}){
    $cycle=0;
    foreach $mb (@members){
      if($mb eq $ky){ #Avoid infinite cycles
	     $cycle=1;
	     #print"Cycle @ $new: reaction $mb already present in \@members.\n"
      }
    }
    if($cycle){
      next;
    }

    if($adyacence{$new}{$ky}){
      #print"New member added - $new / $ky\n";
      push(@members,$ky);
      $ref=&component_members(\@members,$ky);
      @members=@{$ref};
    }
  }
  #print"Back from $new\n";
  return(\@members);
}

