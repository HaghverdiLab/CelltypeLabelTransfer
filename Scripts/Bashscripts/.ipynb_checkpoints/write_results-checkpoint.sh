#!/bin/sh

path=/fast/AG_Haghverdi/Carla_Moelbert/Celltype_annotation
source /home/cmoelbe/bin/anaconda3/bin/activate r-environment 

path_scripts="$path/Scripts/"
resultFolder="$path/Results/Files/"
folder="$path/Data/Predictions_mosaic/"

for set in PBMC; do # Lung PBMC HumanMotorCortex Kidney
    umapFile=$resultFolder/results_curated_$set".csv"
    qsub  -N $set -l m_mem_free=16G -l h_rt=8:30:00 -o out_$set.txt  -e error_$set.txt \
    $path/Scripts/Function_calls/write_results.sh $folder $set $path_scripts $umapFile  \
    $path/Data/Fulldata/$set"_Query"/meta.csv
done