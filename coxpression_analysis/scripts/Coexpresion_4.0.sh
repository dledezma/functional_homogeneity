
# Coexpresion 
# pre-processing of the data

# Downloading colombos data:
cd /home/feraltp/PGC/Coexpresion
wget https://www.dropbox.com/sh/tomhk946s8bqywl/AABez1HUaSD4oHc2EdPwv--Qa/colombos_ecoli_exprdata_20151029.txt?dl=0
mv colombos_ecoli_exprdata_20151029.txt\?dl\=0 colombos_ecoli_exprdata_20151029.txt
sed -i '1,6d' colombos_ecoli_exprdata_20151029.txt
grep -v "#" all_regulons_bnumber.txt | sed 's$/$\t$g' > processed_all_regulons_bnumber.txt


cd /home/feraltp/PGC/Coexpresion/Data/GeneSets
# With the input files in the directory ... :

for i in *.txt
do
   grep -v "#" $i | sed 's$/$\t$g' > processed_$i
done

# Renaming the random carpets according to the names to pass as arguments to the script
cd random

mv random_clusterscoexp random_Clusters
mv random_complex random_Complex
mv random_conformation random_Conformation
mv random_conformation_effect random_ConfEffect
mv random_general random_General
mv random_pathways random_Pathways
mv random_simple random_Simple
mv random_simplex random_Simplex
mv random_strict random_Strict
mv random_tus random_TUs

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

cd /home/feraltp/PGC/Coexpresion/

# Running the code:
nohup Rscript Coexpresion_11.0.R Data/GeneSets/ processed_general_regulons_bnum.txt processed_strict_regulons_bnum.txt processed_simple_regulons_bnum.txt processed_rcomplex_bnum.txt processed_simplex_regulons_bnum.txt processed_regulons_by_conformation_bnum.txt processed_regulons_by_conformation_effect_bnum.txt processed_gos_bnum.txt processed_pathways_bnum.txt processed_all_tus_bnum_uniq.txt processed_clusters_coexp_bnum.txt General Strict Simple Complex Simplex Conf Conf+Eff GOs Pathways TUs Clusters > resultsCoex11.out 2>&1 &

# To test the script into R, omit the line of the arguments and put instead one of the following lines:

# args <- c("Data/GeneSets/", "processed_general_regulons_bnum.txt", "processed_regulons_by_conformation_bnum.txt", "Regulon", "RegulonsByConformationBnum")

# args <- c("Data/GeneSets/", "processed_general_regulons_bnum.txt", "processed_strict_regulons_bnum.txt", "processed_simple_regulons_bnum.txt", "processed_rcomplex_bnum.txt", "processed_simplex_regulons_bnum.txt", "processed_regulons_by_conformation_bnum.txt", "processed_regulons_by_conformation_effect_bnum.txt", "processed_gos_bnum.txt", "processed_pathways_bnum.txt", "processed_all_tus_bnum_uniq.txt", "processed_clusters_coexp_bnum.txt", "General", "Strict", "Simple", "Complex", "Simplex", "Conformation", "ConfEffect", "GOs", "Pathways", "TUs", "Clusters")

# Downloading the results:
scp feraltp@tepeu.lcg.unam.mx:/home/feraltp/PGC/Coexpresion/Results/Tables/* .

# The pictures were converted to jpg using http://convert-my-image.com/PdfToJpg_Es (300 dpi's)