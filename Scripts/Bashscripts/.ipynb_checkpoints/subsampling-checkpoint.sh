#!/bin/sh
source /home/cmoelbe/bin/anaconda3/bin/activate r-environment

path=/fast/AG_Haghverdi/Carla_Moelbert/Celltype_annotation
path_scripts=$path/Scripts
output="$path/Data/Subsets/"
data=$path/Data/Fulldata/


#qsub -N balanced -l m_mem_free=64G  -l h_rt=08:00:00 \
#$path_scripts/Function_calls/subsampling.sh $data/PBMC10x_Reference $output 20 \
#"100,250,500,1000,2000,3000" $path_scripts 
    
#qsub -N bootstrap -l m_mem_free=32G -l h_rt=03:00:00 \
#$path_scripts/Function_calls/subsampling.sh $data/PBMC10x_Reference $output 20 \
#"3090,2418,1373,1022,703,623,273,126,38" $path_scripts 
    
#qsub -N mosaic -l m_mem_free=32G -l h_rt=03:00:00 \
#$path_scripts/Function_calls/subsampling.sh $data/PBMCMosaic_Reference $output 20 \ #"2274,3763,551,4321,6611,209,890,818,102" $path_scripts 

qsub -N m -l m_mem_free=32G -l h_rt=03:00:00 \
$path_scripts/Function_calls/subsampling.sh $data/PBMCMosaic_Reference $output 20 "4321" $path_scripts 

#qsub -N mosaic -l m_mem_free=32G -l h_rt=03:00:00 \
#$path_scripts/Function_calls/subsampling.sh $data/PBMCMosaic_Reference $output 20 \ "890,818,102" $path_scripts 

#qsub -N subsets -l m_mem_free=8G -l h_rt=01:00:00 \
#$path_scripts/Function_calls/subsampling.sh $data/BoneMarrow_Reference $output 10 \
#"18,199,486,492,666,1454,1908,2681,3060,4019,4081,5438,9848,14707" $path_scripts 


 
                               
                                    