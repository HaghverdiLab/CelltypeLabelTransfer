#!/bin/sh



Filelist=($(ls -d $2/$1* | xargs -n1  basename))

source  /home/cmoelbe/bin/anaconda3/bin/activate r-environment
Rscript /fast/AG_Haghverdi/Carla_Moelbert/Celltype_annotation/Scripts/Methods/prep_ItClust.R  ${Filelist[${SGE_TASK_ID}-1]} $2 $3 $4

source /home/cmoelbe/bin/anaconda3/bin/activate ItClust
echo "set" ${Filelist[${SGE_TASK_ID}-1]}
echo "-----------------------------------------------"
start=`date +%s`
python3 /fast/AG_Haghverdi/Carla_Moelbert/Celltype_annotation/Scripts/Methods/itClust.py $4/${Filelist[${SGE_TASK_ID}-1]} $4/${Filelist[${SGE_TASK_ID}-1]} $4/${Filelist[${SGE_TASK_ID}-1]}
end=`date +%s`

echo ItClust all ${Filelist[${SGE_TASK_ID}-1]} `expr $end - $start` >>  /fast/AG_Haghverdi/Carla_Moelbert/Celltype_annotation/Results/Files/runtime.txt

echo "-----------------------------------------------"
source  /home/cmoelbe/bin/anaconda3/bin/activate r-environment
Rscript /fast/AG_Haghverdi/Carla_Moelbert/Celltype_annotation/Scripts/itClust_transformResults.R $4/${Filelist[${SGE_TASK_ID}-1]}