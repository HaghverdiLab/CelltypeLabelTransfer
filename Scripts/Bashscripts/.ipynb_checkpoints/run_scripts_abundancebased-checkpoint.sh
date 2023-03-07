#!/bin/bash
# -cwd

path=/fast/AG_Haghverdi/Carla_Moelbert/Celltype_annotation
path_scripts="$path/Scripts"
path_data="$path/Data/Subsets_mosaic/"
output="$path/Data/Predictions_mosaic/"


for set in PBMC; do 
  echo $set
  test=$path/Data/Fulldata/$set"_Query"
  echo $test
  
  Filelist_length=$(ls -d $path_data/$set"Mosaic"* | wc -l)
  echo ${Filelist_length} 
 # qsub -cwd -N ItClust  -l m_mem_free=32G -l h_rt=0:45:00  -t 1-${Filelist_length}  \
 # $path_scripts/Function_calls/itClust.sh $set $path_data $test $output/ItClust 
  
 # qsub -cwd -N Seurat  -l m_mem_free=32G -l h_rt=0:45:00 -t 1-${Filelist_length} \
 # $path_scripts/Function_calls/seurat.sh $set $path_data $test $output/Seurat
  
 # qsub -cwd -N CellID  -l m_mem_free=64G -l h_rt=0:45:00  -t 1-${Filelist_length} \
 # $path_scripts/Function_calls/cellid.sh $set $path_data $test $output/CellID 
 
  qsub -cwd -N SCN     -l m_mem_free=64G -l h_rt=0:45:00  -t 1-${Filelist_length} \
  $path_scripts/Function_calls/scn.sh $set $path_data $test $output/SingleCellNet
  
 # qsub -cwd -N singleR -l m_mem_free=64G -l h_rt=0:45:00  -t 1-${Filelist_length} \
 # $path_scripts/Function_calls/singler.sh $set $path_data $test $output/SingleR

done