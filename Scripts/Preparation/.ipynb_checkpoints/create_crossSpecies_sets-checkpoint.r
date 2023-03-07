setwd("/fast/AG_Haghverdi/Carla_Moelbert/Benchmarking_TrainingDataSelection")

source("Scripts/functions.R")
source("Scripts/preparation_functions.r")
library(ggplot2)
library(dplyr)

path_raw="Data/raw/CrossSpecies"
path_out ="Data/processed/CrossSpecies"
path_subsets="Data/subsets"
if(!dir.exists(path_out))dir.create(path_out)
#if(!dir.exists(paste0(path_out,"CrossSpecies"))) dir.create(paste0(path_out,"CrossSpecies"))
#print("..............")
# Read raw data 
#data_train <- t(as.matrix(read.table(paste(path_raw,"mouse/mouse_source_matrix_t.txt", sep="/"))))
#meta_train <- read.csv(paste(path_raw,"mouse/mouse_source_metadata.csv", sep="/"))
#meta_train <- getMetaFormat(meta_train)
#colnames(data_train) <- updateID(colnames(data_train))
#print("..............")
#folder_human <- paste(sep="/", path_raw,"human/human_kidney_use/" )
#matrixfile <- paste0(folder_human, "matrix.mtx")
#barcodefile <- paste0(folder_human, "barcodes.tsv")
#featuresfile <- paste0(folder_human, "genes.tsv")
#print("..............")
#data_test <- Matrix::readMM(file= matrixfile)
#colnames(data_test) <- as.character(read.table(barcodefile, header=F)[,1])
#rownames(data_test) <- as.character(read.table(featuresfile, header=F)[,1])
#genes <- rownames(data_test)
#print("..............")
#mart <- biomaRt::useMart("ensembl", dataset = "hsapiens_gene_ensembl",
#                         host = "www.ensembl.org" )
#gene_ids <- biomaRt::getBM(attributes=c("ensembl_gene_id","hgnc_symbol"),
 #                          mart = mart, values = genes)
#translations <- setNames(gene_ids$hgnc_symbol, gene_ids$ensembl_gene_id)
#translations <- translations[translations != ""]
#genes <- stringr::str_replace_all(rownames(data_test),
 #                                 stringr::fixed(translations))
#print("..............")
#rownames(data_test) <- genes
#data_test <- data_test[stringr::str_detect(rownames(data_test), "ENSG0" ,
 #                                    negate = TRUE),]
#colnames(data_test) <- updateID(colnames(data_test))
#meta_test <- getMetaFormat(read.csv(paste(sep="/", path_raw,"human/human_kidney_label_JH.csv" )))

#print("..............")
#genes <- Reduce(intersect, list(rownames(data_train), rownames(data_test)))                                    
#data_train <- as.matrix(data_train[rownames(data_train) %in% genes,])
#data_test <- as.matrix(data_test[rownames(data_test) %in% genes,])

#x <- max(table(meta_train$class_))
#min_n <- min(table(meta_train$class_))
#print("..............")
#print("Get Subset")
#get_subset(data_train,meta_train, paste(sep="/", path_out,"data_train.txt") ,
 #        paste(sep="/", path_out,"meta_train.txt"))
        
#get_subset(data_test,meta_test, paste(sep="/", path_out,"data_test.txt") ,
#         paste(sep="/", path_out,"meta_test.txt"))
print("Get random") 
get_random_sets(input=paste0(path_out), output=paste0(path_subsets), 20, "CrossSpecies", steps=c(228,100,500,1000,2000,3000))
                        
