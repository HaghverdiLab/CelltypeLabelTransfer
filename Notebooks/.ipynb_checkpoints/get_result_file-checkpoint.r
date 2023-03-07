source("../Scripts/functions.R")
source("../Scripts/visulizations.r")
library(ggplot2)
library(dplyr)
library(purrr)
library(Seurat)
library(viridis)
library("RColorBrewer")


getName <- function(folder, position){
  name <- unlist(stringr::str_split(folder,"/"))
  name <- name[length(name)]
  name <- unlist(stringr::str_split(name,"_"))[position]
  return(name)
}

getVector <- function(method, folder){
    if(method %in% c("Seurat", "SingleCellNet", "CellID", "SingleR")){
  if(stringr::str_detect(folder, "predict")) stop("Wrong file")
  name <- paste(sep="_", getName(folder,1), method ,getName(folder,2),
                stringr::str_replace(getName(folder,3), ".txt", "") )
  df <- read.csv(folder, sep="\t")
  df <- df[,c("id", "predicted")] #.match
  colnames(df) <- c("id", name)
    
} else if (method == "ItClust"){
  name <- paste(sep="_",  getName(folder,1), "ItClust",getName(folder,2),
                getName(folder,3))
  df <- read.csv(paste(folder, "results.txt", sep="/"), sep="\t")
  df <- df[,c("class_", "predicted_celltype", "cell_id")]
  colnames(df) <- c("class", "predicted", "id")
  df[,name] <- df$predicted 
  df <- df[, c("id",  name)]
}
    rownames(df) <- df$id
    return(df)
}

get_results_method <- function(resultfolder, id, method){
    print(paste("Start", method, "..."))
    files <- list.files(resultfolder, pattern=id, full.names = T)
    if (method == "Seurat") files <- files[stringr::str_detect(files,
                                                               "predict", negate=T)]
    data <- lapply(files, function(file) getVector(method, file))
    summary <- data %>% reduce(full_join, by = "id")
    return(summary)
}
                   
translate <- function(col){
  
    col[col == "B"] <- "B cell"
    col[col == "Cytotoxic T"] <- "Cytotoxic T cell"
    col[col == "CD4+ T"] <- "CD4+ T cell"
    col[col == "Dendritic"] <- "Dendritic cell"
    col[col == "Natural killer"] <- "Natural killer cell"
    col[col == "Plasmacytoid dendritic"] <- "Plasmacytoid dendritic cell"
    return(col)
} 
                   
adjust_names <- function(data, name="ItClust"){
    cols <- colnames(data)[stringr::str_detect(colnames(data), name)] 
    x<- do.call(cbind,lapply(cols, function(col) translate(data[,col])))
    colnames(x) <- cols  
    id <- data[,!(colnames(data) %in% cols)]

    data <- cbind(x,id)
    return(as.data.frame(data))
}
                             
transform_PBMC_results <- function(data_file, celltypes, methods, sizes,
                                   cols= c("id","nGene", "nUMI", "percent.mito",
                                           "Cluster", "class_", "Experiment", "Method")){
    data <- read.csv(data_file)
    x <- reshape2::melt(data,  id.vars =cols,value.name = "Prediction")
    
    x[c('Reference', 'Approach', "Size", "Set")] <- stringr::str_split_fixed(x$variable, '_', 4)
    x$Match <- x$Prediction == x$class_ 

    x$refSize <- sizes[ match(x$class_, celltypes ) ]
    x$Approach <- factor(x$Approach, levels=methods)
    x$class <- factor(x$class_, levels=celltypes)
    x <- x[,c("id", "Prediction", "Reference", "Approach", "Size", "Set",
              "class", "refSize", "Match")]
    x$Size <- as.numeric(x$Size)
    x <- x[!is.na(x$Match),]
    return(x)
}  
args = commandArgs(trailingOnly=TRUE)
folder <- args[1]
metafile <- args[2]
name <- args[3]                        
out <- args[4]
out_long <- args[5]
 
celltypes = c("Cytotoxic T cell", "CD4+ T cell", "CD14+ monocyte", "B cell",
              "Megakaryocyte", "Natural killer cell",
              "CD16+ monocyte", "Dendritic cell", "Plasmacytoid dendritic cell")
methods <- c("Seurat",  "SingleR","CellID", "SCN", "ItClust")
sizes <- c(3090, 2418, 1373, 1022, 703, 623, 273, 126, 38)
names(sizes) <- celltypes
                             
                             
seurat <- get_results_method(paste(sep="/",folder,"Seurat"), name, "Seurat")
scn <- get_results_method(paste(sep="/",folder,"SingleCellNet" ), name, "SingleCellNet")
itclust <- get_results_method(paste(sep="/",folder,"ItClust" ), name, "ItClust")
itclust <- adjust_names(itclust)
singleR <- get_results_method(paste(sep="/",folder,"SingleR" ), name, "SingleR")
cellid <- get_results_method(paste(sep="/",folder,"CellID" ), name, "CellID")
                             
data <- list(cellid, seurat, scn, singleR) %>% reduce(full_join, by = "id")#, itclust
meta <- read.csv(metafile)
data <- merge(data,meta)  
rownames(data) <- data$id

                             
                             
write.table(data, out, sep=",", col.names=T, row.names=T, quote=F, append=F)  
 


long <- transform_PBMC_results(out, celltypes, methods, sizes)  
                             
write.table(data, out_long, sep=",", col.names=T, row.names=T, quote=F, append=F)