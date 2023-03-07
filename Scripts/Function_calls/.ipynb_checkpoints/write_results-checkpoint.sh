#!/bin/sh
source /home/cmoelbe/bin/anaconda3/bin/activate r-environment

#echo "---------------------- Seurat ----------------------"
#for tag in $1/Seurat/$2* ; do
#   Rscript $3/getUMAPfile.R Seurat $tag $4
#done

#echo "---------------------- SCN ----------------------"
#for tag in $1/SCN/$2* ; do
#    Rscript $3/getUMAPfile.R SCN $tag $4
#done

#echo "---------------------- CellID ----------------------"
#for tag in $1/CellID/$2* ; do
#   Rscript $3/getUMAPfile.R CellID $tag $4
#done

#echo "---------------------- ItClust ----------------------"
for tag in $1/ItClust/$2* ; do
     Rscript $3/itClust_transformResults.R $tag
#     Rscript $3/getUMAPfile.R ItClust $tag $4
done

#echo "---------------------- SingleR ----------------------"
#for tag in $1/SingleR/$2* ; do
#    Rscript $3/getUMAPfile.R SingleR $tag $4
#done



