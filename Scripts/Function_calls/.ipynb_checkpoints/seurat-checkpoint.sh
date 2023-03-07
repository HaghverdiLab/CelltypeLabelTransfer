#!/bin/sh
source  /home/cmoelbe/bin/anaconda3/bin/activate r-environment

Filelist=($(ls -d $2/$1* | xargs -n1  basename))
echo Test ${Filelist[${SGE_TASK_ID}-1]}
echo $5
start=`date +%s`
Rscript /fast/AG_Haghverdi/Carla_Moelbert/Celltype_annotation/Scripts/Methods/seurat.R $2/${Filelist[${SGE_TASK_ID}-1]} $3 $4/${Filelist[${SGE_TASK_ID}-1]} $5 
end=`date +%s`

echo Seurat  $5 ${Filelist[${SGE_TASK_ID}-1]} `expr $end - $start` >>  /fast/AG_Haghverdi/Carla_Moelbert/Celltype_annotation/Results/Files/runtime.txt

