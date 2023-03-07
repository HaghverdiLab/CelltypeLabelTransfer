#!/bin/sh
source  /home/cmoelbe/bin/anaconda3/bin/activate r-environment
start=`date +%s`
Filelist=($(ls -d $2/$1* | xargs -n1  basename))

Rscript /fast/AG_Haghverdi/Carla_Moelbert/Celltype_annotation/Scripts/Methods/singler.R $2 ${Filelist[${SGE_TASK_ID}-1]} $3 $4 $5

end=`date +%s`

echo SingleR  $5 ${Filelist[${SGE_TASK_ID}-1]} `expr $end - $start` >>  /fast/AG_Haghverdi/Carla_Moelbert/Celltype_annotation/Results/Files/runtime.txt
