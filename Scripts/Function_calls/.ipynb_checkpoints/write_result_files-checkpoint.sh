folder=/fast/AG_Haghverdi/Carla_Moelbert/Celltype_annotation/
folder_data=$folder"/Data/"
meta=$folder_data"/Fulldata/PBMC_Query/meta.csv"
name=PBMC10x
folder_results=$folder"/Results/Files"

#Rscript $folder/Notebooks/get_result_file.r $folder_data"/Predictions/" $meta $name $folder_results"/result_general.csv" #$folder_results"/results_long.csv"

Rscript $folder/Scripts/get_confidence_scores.r $folder_data"/Predictions/" $name 38,100,1000

#Rscript $folder/Notebooks/get_umap_file.r $name"_100_" PBMCMosaicBalanced_102_ $folder_data"/Predictions/" $meta $folder_results"/result_general.csv" $folder_results"/results_curated.csv" $folder_results"/umap_data.csv"



#Rscript $folder/Notebooks/get_result_file.r $folder_data"/Predictions_curated/" $meta $name $folder_results"/results_curated.csv" $folder_results"/results_curated_long.csv"

