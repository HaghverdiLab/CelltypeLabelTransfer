#!/bin/sh
source  /home/cmoelbe/bin/anaconda3/bin/activate cell_annotation
start=`date +%s`
Filelist=($(ls -d $2/$1"10x_3090"* | xargs -n1  basename))

Rscript /fast/AG_Haghverdi/Carla_Moelbert/Celltype_annotation/Scripts/Methods/cellid.R $2 ${Filelist[${SGE_TASK_ID}-1]}  $3 $4 $5

end=`date +%s`

echo CellID  $5 ${Filelist[${SGE_TASK_ID}-1]} `expr $end - $start` >>  /fast/AG_Haghverdi/Carla_Moelbert/Celltype_annotation/Results/Files/runtime.txt
