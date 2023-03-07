#!/bin/sh
source  /home/cmoelbe/bin/anaconda3/bin/activate r-environment

echo $1
echo $2/$1
Filelist=($(ls -d $2/$1* | xargs -n1  basename))

echo ${SGE_TASK_ID}
echo ${Filelist[${SGE_TASK_ID}-1]}
start=`date +%s`

Rscript /fast/AG_Haghverdi/Carla_Moelbert/Celltype_annotation/Scripts/Methods/scn.R $2 ${Filelist[${SGE_TASK_ID}-1]}  $3 $4 $5

end=`date +%s`

echo SCN  $5 ${Filelist[${SGE_TASK_ID}-1]} `expr $end - $start` >>  /fast/AG_Haghverdi/Carla_Moelbert/Celltype_annotation/Results/Files/runtime.txt


