#!/usr/bin/env Rscript

visualizeItClust <- function(folder){
  clustering <- read.csv(paste(folder, "results/clustering_results.csv", sep="/"))
  #cells <- read.csv(paste(folder, "cells_test.csv", sep="/"), header=F)

  #clustering$id <- cells$V1
  ids <- read.csv(paste(folder, "results/celltype_assignment.txt", sep="/"), header=F)
  ids <- strsplit(as.character(ids$V1), " be ")
  ids <- as.data.frame(do.call(rbind,ids))
  ids$V1 <- gsub("Cluster ","",ids$V1)
  ids$V2 <- gsub(" cell","",ids$V2)

  truth <- read.csv(paste(folder, "meta_test.csv", sep="/"), sep=",")
  print("Truth")
  print(head(truth))
  print(".......................................")
  x <- as.data.frame(do.call(rbind, strsplit(ids$V1, " ")))
  data <- data.frame(ids$V2, x$V1)
  colnames(data) <- c("class_", "cluster")
  clustering$predicted_celltype <- data$class_[ match(clustering$cluster, data$cluster ) ]

  print("Clustering updated")
  print(head(clustering))
  clustering$cell_id <- gsub("-target","", clustering$cell_id)
  x <- merge(clustering, truth, by.x="cell_id", by.y = "id", ) # cell_id
  x$class_ <- gsub(" cell","", x$class_)
  
  print("final results")
  print(head(x))
  name <- paste(sep="/", folder, "results.txt")
  write.table(x, name, col.names = T, row.names = T, quote = F, sep="\t")
}

args = commandArgs(trailingOnly=TRUE)

folder <- args[1]
print(folder)
visualizeItClust(folder)
print("----------------------------------------------------------------------")
