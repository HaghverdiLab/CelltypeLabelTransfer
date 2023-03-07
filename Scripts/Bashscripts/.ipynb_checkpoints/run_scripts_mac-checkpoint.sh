#!/bin/bash
# -cwd

path=/fast/AG_Haghverdi/Carla_Moelbert/Celltype_annotation
path_scripts="$path/Scripts"
path_data="$path/Data/Subsets_runtime/"
output="$path/Data/Predictions/"


for set in Filelist_length=$(ls -d $path_data/$set* | wc -l); do 
  echo $set
  test=$path/Data/Fulldata/PBMC_Query
  echo $test
  
  sh $path_scripts/Function_calls/itClust.sh $set $path_data $test $output/ItClust 
  
  sh $path_scripts/Function_calls/seurat.sh $set $path_data $test $output/Seurat
  
  sh $path_scripts/Function_calls/cellid.sh $set $path_data $test $output/CellID 
 
  sh $path_scripts/Function_calls/scn.sh $set $path_data $test $output/SCN
  
  sh $path_scripts/Function_calls/singler.sh $set $path_data $test $output/SingleR

done