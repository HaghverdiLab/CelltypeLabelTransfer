#!/bin/bash
# -cwd

path=/fast/AG_Haghverdi/Carla_Moelbert/Celltype_annotation
path_scripts="$path/Scripts"
path_data="$path/Data/Subsets_bonemarrow/"
output="$path/Data/Predictions_bonemarrow/"
set=BoneMarrow
test=$path/Data/Fulldata/$set"_Query"

f=all

Filelist_length=$(ls -d $path_data/$set* | wc -l)
echo ${Filelist_length} 
  
#qsub -cwd -N Seurat  -l m_mem_free=32G -l h_rt=0:45:00 -l cpu=8 -t 1-${Filelist_length} \
#$path_scripts/Function_calls/seurat.sh $set $path_data $test $output/Seurat $f

#qsub -cwd -N CellID  -l m_mem_free=64G -l h_rt=0:45:00 -l cpu=8  -t 1-${Filelist_length} \
#$path_scripts/Function_calls/cellid.sh $set $path_data $test $output/CellID  $f

#qsub -cwd -N SCN     -l m_mem_free=64G -l h_rt=0:45:00 -l cpu=8  -t 1-${Filelist_length} \
#$path_scripts/Function_calls/scn.sh $set $path_data $test $output/SingleCellNet $f

#qsub -cwd -N singleR -l m_mem_free=64G -l h_rt=0:45:00 -l cpu=8  -t 1-${Filelist_length} \
#$path_scripts/Function_calls/singler.sh $set $path_data $test $output/SingleR $f

qsub -cwd -N ItClust  -l m_mem_free=64G -l h_rt=0:45:00 -l cpu=8  -t 1-${Filelist_length} \
$path_scripts/Function_calls/itClust.sh $set $path_data $test $output/ItClust 
