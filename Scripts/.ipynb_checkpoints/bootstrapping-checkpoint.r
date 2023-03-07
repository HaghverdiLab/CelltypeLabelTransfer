library(ggplot2)
library(dplyr)
#setwd("/fast/AG_Haghverdi/Carla_Moelbert/Cell_annotation/Scripts/")
get_combinations <- function(files, size, seed){
    set.seed(seed)
    return(sample(files, size = size))
}

translate <- function(vector){
    vector <- stringr::str_replace_all(vector, "\\.\\.", "+ ")
    vector <- stringr::str_replace_all(vector, "\\.", " ")
    return(vector)
}

readSeurat <- function(file){
    seurat <- read.csv(file, sep="\t")
    seurat$prediction.score.max <- NULL
    seurat <- seurat[, stringr::str_detect(colnames(seurat), "prediction.score.")]
    colnames(seurat) <- stringr::str_replace(colnames(seurat), "prediction.score.", "")
    colnames(seurat) <- translate(colnames(seurat))
    seurat[is.na(seurat)] <- 0
    return(seurat)
}

readSCN <- function(file){
    scn <- read.csv(file, sep="\t")
    rownames(scn) <- scn$id 
    scn <- scn[, !colnames(scn) %in% c("id", "predicted", "nGene", "nUMI", "percent.mito",
                                       "Cluster", "class_", "Experiment", "Method",
                                       "prediction.match")]
    colnames(scn) <- translate(colnames(scn))
    scn[is.na(scn)] <- 0
   
    return(as.data.frame(scn))
}

readSingleR <- function(file){
    singler <- read.csv(file, sep="\t")
    rownames(singler) <- singler$id 
    
    singler <- singler[, stringr::str_detect(colnames(singler), "scores.")]
    singler <- singler[, !stringr::str_detect(colnames(singler), "tuning")]
    colnames(singler) <- stringr::str_replace(colnames(singler), "scores.", "")
    colnames(singler) <- translate(colnames(singler))
    singler[is.na(singler)] <- 0
    return(singler)
}

readItclust <- function(file){

    itclust <- read.csv(paste(sep="/", file, "results/clustering_prob.csv"))
    ids <- read.csv(paste(sep="/", file, "results/celltype_assignment.txt"), header=F)
    rownames(itclust) <- itclust$X
    itclust$X <- NULL

    ids[c("x", 'Cluster', 'x', "confidence", "x", "x", "type")] <-  stringr::str_split_fixed(ids$V1, ' ', 7)
    ids$type <- unlist(lapply(ids$type, function(type) substr(type, 1, nchar(type)-5)))
    ids$Cluster <- paste0("cluster", ids$Cluster)

    colnames(itclust) <- ids$type[ids$Cluster ==colnames(itclust)]
    itclust[is.na(itclust)] <- 0
    return(itclust)

}
                              
filter_set<- function(set, cutoff=NULL){
    set[is.na(set)] <- 0
    if(!is.null(cutoff)) set[set < cutoff]<- 0
    return(set)
} 
                              
summarize_binary_prediction <- function(set){
   x <- which.max(table(set))
   return(set[x]) 
}

get_binary_prediction <- function(set){
    predictions <-  apply(set, 1, which.max)
    predictions <- colnames(set)[predictions]
    names(predictions) <- rownames(set)
    
    return(predictions)
}

filter_data <- function(data, cutoff){
    
    data[data < cutoff] <- 0
    return(data)
}
get_confidence_df <- function(subsets, nr_sets){
    combs <- Reduce('+', subsets)
    combs <- combs / nr_sets                    
    prediction <- apply(combs, 1, which.max)
    prediction <- colnames(combs)[prediction] 
    return(prediction)
}
get_predicton <- function(sets, combinations, method, set, cutoff=0.25){
    nr_sets <- length(combinations)
    subsets <- sets[names(sets) %in% combinations]
    
    prediction <- sapply(subsets, function(set) get_binary_prediction(set)) 
    prediction_binary <- apply(prediction, 1, summarize_binary_prediction)
    
    df <- data.frame(id=rownames(combs), method=method, set=set,prediction_binary = prediction_binary)   
    return(df)
}
    
                         
print("Start....")
size=20
cutoff = 0.3
sequence <- seq(1,20,1)
                              
files <- list.files("/fast/AG_Haghverdi/Carla_Moelbert/Celltype_annotation/Data/Predictions/ItClust/",
                    pattern = "PBMCMosaicBalanced_102_")
                         
print(length(files))
print(head(files))
combinations <- lapply(sequence, function(seed) get_combinations(files, size, seed))

print("Get CellID...")
files_cellid <- lapply(files, function(file) list.files("/fast/AG_Haghverdi/Carla_Moelbert/Celltype_annotation/Data/Predictions/CellID/",
                                                     pattern = paste0(file, ".txt"),
                                                     full.names = T))

sets_cellid  <- lapply(files_cellid, function(file) readCellID(file))                     
names(sets_cellid)<- files
print("--------------------------------------------------------")
cellid_pred <- do.call(rbind,lapply(sequence, function(seq) get_predicton(sets_cellid,
                                                                       unlist(combinations[seq]),                                                                      "CellID", seq)))                         
                       
print("Get Seurat...")
files_seurat <- lapply(files, function(file) list.files("/fast/AG_Haghverdi/Carla_Moelbert/Celltype_annotation/Data/Predictions/Seurat/",pattern = paste0(file, "_predictions.txt"),
                                                        full.names = T))
              
sets_seurat  <- lapply(files_seurat, function(file) readSeurat(file))
names(sets_seurat)<- files
seurat_pred <- do.call(rbind,lapply(sequence, function(seq) get_predicton(sets_seurat,
                                                                       unlist(combinations[seq]),
                                                                       "Seurat", seq, cutoff)))
                                   
print("Get SCN...")                                    
files_scn <- lapply(files, function(file) list.files("/fast/AG_Haghverdi/Carla_Moelbert/Celltype_annotation/Data/Predictions/SCN/",
                                                     pattern = paste0(file, ".txt"),
                                                     full.names = T))
sets_scn  <- lapply(files_scn, function(file) readSCN(file))
names(sets_scn)<- files
scn_pred <- do.call(rbind,lapply(sequence, function(seq) get_predicton(sets_scn,
                                                                       unlist(combinations[seq]),
                                                                       "SCN", seq)))

print("Get SingleR")
files_singler <- lapply(files, function(file) list.files("/fast/AG_Haghverdi/Carla_Moelbert/Celltype_annotation/Data/Predictions/SingleR//",
                                                         pattern = paste0(file, ".txt"), 
                                                         full.names = T))
sets_singler  <- lapply(files_singler, function(file) readSingleR(file))
names(sets_singler )<- files
singler_pred <- do.call(rbind,lapply(sequence, function(seq) get_predicton(sets_singler,
                                                                       unlist(combinations[seq]),
                                                                       "SingleR", seq)))

print("Get ItClust...")
files_itclust <- paste0("/fast/AG_Haghverdi/Carla_Moelbert/Celltype_annotation/Data/Predictions/ItClust/", files)
sets_itclust  <- lapply(files_itclust, function(file) readItclust(file))
names(sets_itclust )<- files
itclust_pred <- do.call(rbind,lapply(sequence, function(seq) get_predicton(sets_itclust,
                                                                       unlist(combinations[seq]),
                                                                       "ItClust", seq)))
itclust_pred$id <- stringr::str_replace(itclust_pred$id, "-target", "")
                              
                                   
                                     
print("Get Reference...")
ref <- read.csv("/fast/AG_Haghverdi/Carla_Moelbert/Celltype_annotation/Data/Fulldata/PBMC_Query/meta.csv")
data <- read.csv("/fast/AG_Haghverdi/Carla_Moelbert/Celltype_annotation/Results/Files/results_PBMCMosaic.csv")
print(head(data))                                  
celltypes = c("Cytotoxic T cell", "CD4+ T cell", "CD14+ monocyte", "B cell", "Megakaryocyte",
              "Natural killer cell", "CD16+ monocyte", "Dendritic cell",
              "Plasmacytoid dendritic cell")
x <- reshape2::melt(data,  id.vars =c("id"),value.name = "Prediction")
x[c('Reference', 'method', "Size", "Set")] <- stringr::str_split_fixed(x$variable, '_', 4)

x$Prediction[x$Prediction =="Cytotoxic T"] <- "Cytotoxic T cell"
x$Prediction[x$Prediction == "B"] <- "B cell"
x$Prediction[x$Prediction == "CD4+ T"] <- "CD4+ T cell"
x$Prediction[x$Prediction == "Natural killer"] <- "Natural killer cell"
x$Prediction[x$Prediction == "Dendritic"] <- "Dendritic cell"
x$Prediction[x$Prediction == "Plasmacytoid dendritic"] <- "Plasmacytoid dendritic cell"
                                     
x <- merge(x, ref, by="id")                                     
x$match_full <- x$Prediction == x$class_ 
x$match_full[is.na(x$match_full)] <- FALSE                                     
                                     
query <- x[x$Size == 5503,]
print(head(query))
print(unique(query$method))   
                                     
print("Get files....")
predictions <- do.call(rbind, list(seurat_pred, scn_pred, singler_pred, itclust_pred))
print(unique(predictions$method))
                                     
data <- merge(query[,c("id", "method", "match_full", "class_")], predictions, by=c("id", "method"))
print(unique(data$method))
data$match_confidence <- data$class_ == data$prediction_confidence
data$match_binary <- data$class_ == data$prediction_binary
data$match_filtered <- data$class_ == data$prediction_filtered
data$match_cf <- data$class_ == data$prediction_confidenceFiltered

write.table(data,
            "/fast/AG_Haghverdi/Carla_Moelbert/Celltype_annotation/Results/Files/PBMC_bootstrap_mosaic_02.csv",
            col.names=T, row.names=F, quote=F,sep=",")
                                     
summary <- data %>% 
           dplyr::group_by(method, class_, set) %>% 
dplyr::summarize(accuracy_confidence= mean(match_confidence),
                 accuracy_full=mean(match_full),
                 accuracy_binary=mean(match_binary),
                 accuracy_BinaryFiltered=mean(match_filtered),
                 accuracy_confidenceFiltered= mean(match_cf)) 
print(min(summary$match_full))                                    
                                     
write.table(summary,
            "/fast/AG_Haghverdi/Carla_Moelbert/Celltype_annotation/Results/Files/PBMC_bootstrap_summary_mosaic_02.csv",
            col.names=T, row.names=F, quote=F,sep=",")