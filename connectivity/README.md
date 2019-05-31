# Connectivity Analysis
Scripts, pipeline and results of connectivity analysis can be found here.

## connectivity_library ##
All Perl scripts and data sets used as input for connectivity calculation.

## connectivity_pipeline.pl ##
Perl script that executes the complete connectivity pipeline, including generation of random sets and calculation of random connectivity. Arguments received are:
    1) path to connectivity library*.
    2) path to original data sets directory*.
    3) path to output directory*.
    * Path should end in "/" for correct execution.
A Pathway Tools session must be active in -api -lisp mode.
Example of execution:
$ perl /Users/user/calculate_connectivity_wenzymes.pl /Users/user/connectivity/connectivity_library/ /Users/user/original_data_sets/data_sets/ /Users/user/connectivity/results/

## results ##
Contains results of connectivity pipeline using Pathway Tools v22 on original data sets included in this repository.

## general_regulons_noenz ##
Contains results of replicating the pipeline reported in [Ledezma-Tejeida et. al., 2017](https://www.frontiersin.org/articles/10.3389/fmicb.2017.01466/full) using version 22 of Pathway Tools. GENSOR Units included here does not contaon enzymatic regulation.

Network representations of GENSOR Units of general regulons are available at [RegulonDB](http://regulondb.ccg.unam.mx/central_panel_menu/integrated_views_and_tools/gensor_unit_groups).

###### Questions and comments please contact Daniela Ledezma @ dledezma@lcg.unam.mx 
