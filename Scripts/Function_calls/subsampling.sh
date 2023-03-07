#!/bin/sh
source /home/cmoelbe/bin/anaconda3/bin/activate r-environment
start=`date +%s`
echo $1 $2 $3 $4 
Rscript $5/subsampling.R $1 $2 $3 $4 

end=`date +%s`
