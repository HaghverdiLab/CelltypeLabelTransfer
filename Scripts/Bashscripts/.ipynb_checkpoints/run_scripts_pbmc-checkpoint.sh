#!/bin/bash
# -cwd

path=/fast/AG_Haghverdi/Carla_Moelbert/Celltype_annotation
path_scripts="$path/Scripts"
path_data="$path/Data/Subsets/"
output="$path/Data/Predictions/"


for set in PBMC; do #PBMC
  echo $set
  test=$path/Data/Fulldata/$set"_Query"
  echo $test
  
  Filelist_length=$(ls -d $path_data/$set"10x_3090"* | wc -l)


  for genes in  2000; do # 200 1000
     # qsub -cwd -N Seurat  -l m_mem_free=64G -l h_rt=0:15:00  -l cpu=8 -t 1-${Filelist_length} \
     # $path_scripts/Function_calls/seurat.sh $set $path_data $test $output/Seurat $genes

      qsub -cwd -N CellID  -l m_mem_free=64G -l h_rt=0:15:00  -l cpu=8  -t 1-2 \
      $path_scripts/Function_calls/cellid.sh $set $path_data $test $output/CellID $genes

     # qsub -cwd -N SCN   -l m_mem_free=64G -l h_rt=0:15:00 -l cpu=8  -t 1-${Filelist_length} \
     # $path_scripts/Function_calls/scn.sh $set $path_data $test $output/SingleCellNet $genes

     # qsub -cwd -N singleR -l m_mem_free=64G -l h_rt=0:15:00 -l cpu=8 -t 1-${Filelist_length} \
     # $path_scripts/Function_calls/singler.sh $set $path_data $test $output/SingleR $genes

done

#  echo ${Filelist_length} 
#  qsub -cwd -N ItClust  -l m_mem_free=64G -l h_rt=1:35:00 -l cpu=8   -t 30-${Filelist_length}  \
#  $path_scripts/Function_calls/itClust.sh $set $path_data $test $output/ItClust 
  
done
