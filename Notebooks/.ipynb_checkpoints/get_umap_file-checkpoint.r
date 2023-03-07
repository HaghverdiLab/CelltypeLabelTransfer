source("/fast/AG_Haghverdi/Carla_Moelbert/Celltype_annotation/Notebooks/functions.r")
  
getVectors <- function(method, data){

    summary <- data[data$method == method,]  %>% 
           dplyr::group_by(id, class_, tech, full, curated) %>% 
           dplyr::summarize("mono" = mean(bootstrap_mono),
                            "mosaic" = mean(bootstrap_mosaic),
                            "individual" = mean(individual),
                            "mono_dif" = mean(mono_dif),
                            "mosaic_dif" = mean(mosaic_dif))
    colnames(summary) <- c("id", "class","tech", paste(sep="_", "full", method) ,
                           paste(sep="_", "curated", method),
                           paste(sep="_","mono", method),
                           paste(sep="_","mosaic", method),
                           paste(sep="_","individual", method),
                           paste(sep="_","mono_dif", method),
                           paste(sep="_","mosaic_dif", method))
    return(summary)
}

getVectors_set <- function(method, data){

    summary <- data[data$method == method,]  %>% 
           dplyr::group_by(set, class_, tech, full, curated) %>% 
           dplyr::summarize("mono" = mean(bootstrap_mono),
                            "mosaic" = mean(bootstrap_mosaic),
                            "individual" = mean(individual),
                            "mono_dif" = mean(mono_dif),
                            "mosaic_dif" = mean(mosaic_dif))
    colnames(summary) <- c("id", "class","tech", paste(sep="_", "full", method) ,
                           paste(sep="_", "curated", method),
                           paste(sep="_","mono", method),
                           paste(sep="_","mosaic", method),
                           paste(sep="_","individual", method),
                           paste(sep="_","mono_dif", method),
                           paste(sep="_","mosaic_dif", method))
    return(summary)
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

readCellID <- function(folder, id){
    file <- list.files(paste(sep="/",folder, "CellID"),pattern = paste0(id, ".txt"),full.names = T)
    data <- read.csv(file, sep="\t")
    data <- data[, c("id", "predicted")]
    colnames(data)[2] <- id
    return(data)
}  
size=20
cutoff = 0.3
sequence <- seq(1,20,1)
args = commandArgs(trailingOnly=TRUE)

name_mono <- args[1]
name_mosaic <- args[2]
path <- args[3]
ref_file <- args[4]
result_file <- args[5]
curated_file <- args[6]
out_umap <- args[7]
print("Start...")

mosaic <- get_bootstrap_df(name_mosaic, path)
colnames(mosaic)[4] <- "bootstrap_mosaic"

mono <- get_bootstrap_df(name_mono, path)
colnames(mono)[4] <- "bootstrap_mono"
bootstrap <- merge(mono, mosaic, by=c("id", "method", "set"))

ref <- read.csv(ref_file)
ref <- ref[,c("id", "class_")]

data <- merge(ref, bootstrap, by=c("id"), all=T)
data$match_mosaic <- data$class_ == data$bootstrap_mosaic
data$match_mosaic[data$match_mosaic == TRUE] <- 1
data$match_mono <- data$class_ == data$bootstrap_mono
data$match_mono[data$match_mono == TRUE] <- 1

summary <- data %>% 
           dplyr::group_by(method, class_, id) %>% 
dplyr::summarize(bootstrap_mono= mean(match_mono),
                 bootstrap_mosaic= mean(match_mosaic)) 


data <- read.csv(result_file)
full <- data[, stringr::str_detect(colnames(data), "3090")]
full$class_ <- data$class_
full$id <- data$id
full$tech <- data$Method
full <- reshape2::melt(full,id=c("class_", "id", "tech"), value.name = "full")
full[c('reference', 'method', "size", "set")] <- stringr::str_split_fixed(full$variable, '_', 4)
full <- full[, c("id", "class_", "method", "full", "tech")]
                                     
individual <- data[, stringr::str_detect(colnames(data), "_100_")]
individual$class_ <- data$class_
individual$id <- data$id
individual <- reshape2::melt(individual,id=c("class_", "id"), value.name = "individual")
individual[c('reference', 'method', "size", "set")] <- stringr::str_split_fixed(individual$variable, '_', 4)
individual <- individual[, c("id", "class_", "method", "individual")]


individual$match <- individual$class_ == individual$individual
individual$match[individual$match == TRUE] <- 1

ind<- individual %>% 
           dplyr::group_by(method, class_, id) %>% 
dplyr::summarize(individual= mean(match)) 

umap_data <- merge(full, ind, by=c("id", "class_", "method"))
umap_data <- merge(summary, umap_data, by=c("id", "class_", "method")) 
umap_data$mono_dif <- umap_data$bootstrap_mono - umap_data$individual
umap_data$mosaic_dif <- umap_data$bootstrap_mosaic - umap_data$bootstrap_mono    
                                     
data <- read.csv(curated_file)
full <- data[, stringr::str_detect(colnames(data), "3090")]
full$class_ <- data$class_
full$id <- data$id
full$tech <- data$Method
full <- reshape2::melt(full,id=c("class_", "id", "tech"), value.name = "full")
full[c('reference', 'method', "size", "set")] <- stringr::str_split_fixed(full$variable, '_', 4)
full <- full[, c("id", "class_", "method", "full", "tech")]

colnames(full)[4] <- "curated"
umap_data <- merge(umap_data, full, by=c("id", "class_","method", "tech"),all =TRUE)
head(umap_data)
                                     
methods <- c("Seurat",  "SingleR","CellID", "SCN", "ItClust")
df <- lapply(methods, function(method) getVectors(method, umap_data))
length(df)
lapply(df, nrow)
df <- Reduce(function(x, y) merge(x, y, by=c("id", "class", "tech")),df)
rownames(df)<- df$id
write.table(df, out_umap, col.names=T, row.names=T, quote=T, sep=",")            