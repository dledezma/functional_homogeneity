# Coexpression Analysis

## Getting started

### Prerequisites

The programs and dependencies needed are listed. The versions used to run the analyses are in parentheses.

* R (version 3.4.1)

    The followings libraries for R:
    
    * ggplot2 (2.2.1)
    * gplots (3.0.1)
    * RColorBrewer (1.1-2)
    * parallel (3.4.1)

* Unix-like terminal interface to run bash commands and manage codes and data files

## Directory structure and descriptions

After downloading this repository the following files and directories should be found

```
.
├── Coexpresion_12.0.R                                       # Main code
├── Data                                                     # Directory with all the data needed
│   ├── Colombos                                            # Colombos database directory
│   │   └── colombos_ecoli_exprdata_20151029.txt             # Colombos database (version 3.0)
│   └── GeneSets                                             # Gene sets directory
│       ├── gos_bnum.txt                                     # The gene set files, including the positive control
  [...]   
│       └── random                                           # The directory with the random controls
│           ├── random_GOs                                   # Directories for each data set.
│           │   ├── ran_1.txt                                # There are 100 randomnizations of the gene set
  [...]
│           │   └── ran_100.txt
  [...]
│           └── random_simple
└── Results                                                  # All the results will be here after running the script
    ├── Plots                                                # Graphic results
    └── Tables                                               # Tabular results, these data can be used to plot the graphic results
```

Additionally, there are already the correlation tables (Results/Tables), there is no need to run the code for this.

## Pre-processing files

Through the terminal, before running the R code, the files must be pre-processed as shown next. Be sure the working directory is ```Coexpresion/Data/GeneSets```

```
# For the gene set files:

for i in *.txt
do
   grep -v "#" $i | sed 's$/$\t$g' > processed_$i
done

# For the random files:

for i in *
do
    cd $i

    for j in ran_*.txt
    do
        grep -v "#" $j | sed 's$/$\t$g' > processed_$j
    done

    echo "$i ... done"
    cd ..
done
```
## Running the code

To run the analyses be sure to be in the principal directory, i.e. the directory with the Data and Results directories and the Coexpresion_11.0.R file

```
nohup Rscript Coexpresion_11.0.R Data/GeneSets/ processed_general_regulons_bnum.txt processed_strict_regulons_bnum.txt processed_simple_regulons_bnum.txt processed_rcomplex_bnum.txt processed_simplex_regulons_bnum.txt processed_regulons_by_conformation_bnum.txt processed_regulons_by_conformation_effect_bnum.txt processed_gos_bnum.txt processed_pathways_bnum.txt processed_all_tus_bnum_uniq.txt processed_clusters_coexp_bnum.txt General Strict Simple Complex Simplex Conf Conf+Eff GOs Pathways TUs Clusters > resultsCoex.out 2>&1 &
```
### Important notes

The code uses several cores, to change this, the code should be modified. See internal documantation to know more about it.

Please note that COLOMBOS Compendia must be downloaded [here](http://colombos.net/cws_data/compendium_data/ecoli_compendium_data.zip) and placed in Data/Colombos. GitHub file size restrictions forbid us from keeping it here.

The code can be run with the next command, the output will be printed to screen, as it will be running in the foreground. Any interruption of the connection will halt the analyses. The previous way is recommended instead of the following.

```
Rscript Coexpresion_11.0.R Data/GeneSets/ processed_general_regulons_bnum.txt processed_strict_regulons_bnum.txt processed_simple_regulons_bnum.txt processed_rcomplex_bnum.txt processed_simplex_regulons_bnum.txt processed_regulons_by_conformation_bnum.txt processed_regulons_by_conformation_effect_bnum.txt processed_gos_bnum.txt processed_pathways_bnum.txt processed_all_tus_bnum_uniq.txt processed_clusters_coexp_bnum.txt General Strict Simple Complex Simplex Conf Conf+Eff GOs Pathways TUs Clusters
```

The progress of the code can be partially monitored by reading the "resultsCoex.out" file.
```
tail -f resultsCoex.out
```

## References

Moretto, M., Sonego, P., Dierckxsens, N., Brilli, M., Bianco, L., Ledezma-Tejeida, D., ... & Collado-Vides, J. (2015). COLOMBOS v3. 0: leveraging gene expression compendia for cross-species analyses. Nucleic acids research, 44(D1), D620-D623.

###### Questions and comments please contact Daniela Ledezma or Fernando Altamirano @ dledezma@lcg.unam.mx / feraltp@lcg.unam.mx
