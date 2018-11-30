## go_analysis.pl ##
Perl script that executes the complete go analysis pipeline, including generation of random sets and their analysis. Arguments received are:
1) path to go_analysis library*.
2) path to original data sets directory*.
3) path to output directory*.
* Path should end in "/" for correct execution.
A Pathway Tools session must be active in -api -lisp mode.
Example of execution:
$ perl /Users/user/go_analysis.pl /Users/user/go_analysis/goanalysis_library/ /Users/user/original_data_sets/data_sets/ /Users/user/go_analsisy/results/


## goanalysis_library ##
All files necessary for GO analysis. File names should be preserved.

## results ##
Results from using pipeline on reported data sets.


###### Questions and comments please contact Daniela Ledezma @ dledezma@lcg.unam.mx 
