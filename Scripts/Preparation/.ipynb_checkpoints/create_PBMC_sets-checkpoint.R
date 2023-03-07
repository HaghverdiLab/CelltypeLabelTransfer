#' #############################################################################
#' Prepare data for the cell type classification approaches 
#' #############################################################################
#library(dorothea)
library(tibble)
library(tidyr)
library(dplyr)
library(ggplot2)
library(readr)

setwd("/fast/AG_Haghverdi/Carla_Moelbert/Benchmarking_TrainingDataSelection/Data")
source("../Scripts/preparation_functions.r")
source("../Scripts/functions.R")
#' ################################# PBMC data #################################
#' Prepare the full set to have the files in the right format
#' #############################################################################

path_raw="raw/PBMC"
path_out ="processed"
path_subsets="subsets"

#fullData <- as.matrix(Matrix::readMM(paste(sep="/", path_raw, "counts.umi.txt.gz")))
#colnames(fullData) <-  read.csv(paste(sep="/", path_raw, "cells.umi.new.txt"), header=F)[,1]
#rownames(fullData)<- read.csv(paste(sep="/", path_raw, "genes.umi.txt"), header=F)[,1]
#rownames(fullData) <- unlist(lapply(rownames(fullData),
 #                                   function(id) unlist(stringr::str_split(id, "_"))[2]))
#fullData <- fullData[!duplicated(rownames(fullData)),]
#fullMeta <- read.csv(paste(sep="/", path_raw, "SCP424/metadata/meta.txt"), sep="\t")
#fullMeta <- fullMeta[-1,]
#print(unique(fullMeta$Method))
#fullMeta <- fullMeta[fullMeta$CellType != "Unassigned",]
#print("------------------------------------")                                    
#print(table(fullMeta$Method))
#meta_mono <- fullMeta[grepl( "inDrops", fullMeta$Method, fixed = TRUE),]
#print("Mono")
#print(nrow(meta_mono))                                    
#meta_test <- fullMeta[fullMeta$Method == "Drop-seq",]
#print("Test")
#print(nrow(meta_test))                                    
#meta_mosaic  <- fullMeta[fullMeta$Method != "Drop-seq",]
#print("Mosaic")
#print(nrow(meta_mosaic))
                                    
#print("------------------------------------")                                    
#meta_mono <- meta_mono[,c("NAME","CellType", "Method")]
#colnames(meta_mono) <- c("id", "class_", "tech")

#meta_mosaic <- meta_mosaic[,c("NAME","CellType", "Method")]
#colnames(meta_mosaic) <- c("id", "class_", "tech")

#meta_test <- meta_test[,c("NAME","CellType")]
#colnames(meta_test) <- c("id", "class_")

#print("------------------------------------")
                                 
#get_subset(fullData, meta_mono, paste(sep="/", path_out,"PBMC_mono", "data_train.txt") ,
#           paste(sep="/", path_out,"PBMC_mono", "meta_train.txt"))

                                    print("Mosaic")                                    
#get_subset(fullData, meta_mosaic, paste(sep="/", path_out,"PBMC_Mosaic", "data_train.txt") ,
#         paste(sep="/", path_out,"PBMC_Mosaic", "meta_train.txt"))

                                    
#print(head(meta_test))
#print(meta_test)
#get_subset(fullData,meta_test, paste(sep="/", path_out,"PBMC_mono", "data_test.txt") ,
#         paste(sep="/", path_out,"PBMC_mono", "meta_test.txt"))

                                    
#print("Get Subsets....")                                    
#x <- max(table(meta_mono$class_))
#min_n <- min(table(meta_mono$class_))                                   
get_random_sets(input=paste(sep="/",path_out, "PBMC_mono"),
                output=path_subsets, 20, "Mono", steps=c(42, 100,500,1000,2000,3000))

#x <- max(table(meta_mosaic$class_))
#min_n <- min(table(meta_mosaic$class_))
#print(table(meta_mosaic$class_))                                     
#print(paste(x, min_n))                                    
get_random_sets(input=paste(sep="/",path_out, "PBMC_Mosaic"),
              output=path_subsets, 20, "Mosaic", steps=c(136))# 100,500,1000,2000,3000