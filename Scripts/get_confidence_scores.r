library(dplyr)

get_combinations <- function(files, size, seed){
    set.seed(seed)
    return(sample(files, size = size))
}

translate <- function(vector){
    vector <- stringr::str_replace_all(vector, "\\.\\.", "+ ")
    vector <- stringr::str_replace_all(vector, "\\.", " ")
    return(vector)
}

adjust_df <- function(data, tag, method){
    data$method = method
    x <-  stringr::str_split_fixed(tag, '_', 4)  
    data$reference <- x[1]
    data$size <- x[2]
    data$set <- x[3]
    if(x[4]=="")data$genes <-  0
    else data$genes <-  x[4]
    
    return(data)
}


readConfidence_Seurat <- function(folder, tag){
    file <- paste0(paste(sep="/",folder, tag), "_predictions.txt")
    
    seurat <- read.csv(file, sep="\t")
    seurat$id <- rownames(seurat)
    seurat <- seurat[, c("predicted.id", "prediction.score.max",
                                         "class_", "prediction.match", "id")]
    colnames(seurat) <- c("predicted", "score", "class", "match", "id")
    seurat <- adjust_df(seurat, tag, "Seurat")
    return(seurat)
}





readConfidence_SCN <- function(folder, tag){
    scn <- read.csv(paste0(paste(sep="/",folder, tag), ".txt"), sep="\t")
    rownames(scn) <- scn$id 
    scn$predicted[is.na(scn$predicted)]<- "rand"
    scn <- scn[!is.na(scn$class_),] 
    scn <- scn[!is.na(scn$prediction.match),]
    
 
    scn$predicted <- stringr::str_replace_all(scn$predicted, "\\+", "\\.")
    scn$predicted <- stringr::str_replace_all(scn$predicted, " ", "\\.")
    
    scn$score <- apply(scn, 1, function(row) row[names(row) == row["predicted"]])
    
    scn <-  scn[, c("predicted", "score","class_", "prediction.match", "id")]
    colnames(scn) <- c("predicted", "score", "class", "match", "id")
    scn <- adjust_df(scn, tag, "SCN")   
    return(as.data.frame(scn))
}

readConfidence_SingleR <- function(folder, tag){
    singler <- read.csv(paste0(paste(sep="/",folder, tag), ".txt"), sep="\t")
    rownames(singler) <- singler$id 
    
    singler$predicted <- stringr::str_replace_all(singler$predicted, "\\+", "\\.")
    singler$predicted <- stringr::str_replace_all(singler$predicted, " ", "\\.")
    singler$score <- apply(singler, 1, function(row) row[names(row) == paste(sep=".","scores",row["predicted"])])

    singler <-  singler[, c("predicted", "score","class_", "prediction.match", "id")]
    colnames(singler) <- c("predicted", "score", "class", "match", "id")
    singler <- adjust_df(singler, tag, "SingleR")   
                     
    return(singler)
}

readConfidence_Itclust <- function(folder, tag, meta){
    file <- paste(sep="/", folder, tag)
    itclust <- read.csv(paste(sep="/", file, "results/clustering_prob.csv"))
    ids <- read.csv(paste(sep="/", file, "results/celltype_assignment.txt"), header=F)
    rownames(itclust) <- itclust$X
    itclust$X <- NULL
    ids[c("x", 'Cluster', 'x', "confidence", "x", "x", "type")] <-  stringr::str_split_fixed(ids$V1,
                                                                                             ' ', 7)
    ids$type <- unlist(lapply(ids$type, function(type) substr(type, 1, nchar(type)-5)))
    ids$Cluster <- paste0("cluster", ids$Cluster)
    colnames(itclust) <- ids$type[ids$Cluster ==colnames(itclust)]
    colnames(itclust) <- translate(colnames(itclust))
    itclust$predicted <- unlist(apply(as.data.frame(itclust), 1,
                                      function(row) names(row)[which.max(row)]))
    itclust$id <- stringr::str_replace(rownames(itclust),"-target", "")
                                       
    itclust <- merge(itclust, meta)                          
    itclust$match <- itclust$class_ == itclust$predicted 

    itclust$score <- apply(itclust, 1, function(row) row[names(row) == row["predicted"]])

    itclust <-  itclust[, c("predicted", "score","class_", "match", "id")]
    colnames(itclust) <- c("predicted", "score", "class", "match", "id")
                                                             
    itclust <- adjust_df(itclust, tag, "ItClust")                                    
                            
    return(itclust)

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

get_predicton <- function(sets, combinations, method, set, cutoff=0.25){
    nr_sets <- length(combinations)
    subsets <- sets[names(sets) %in% combinations]
    if(method != "CellID") prediction <- sapply(subsets, function(set) get_binary_prediction(set)) 
    else { 
        prediction <- Reduce(function(x, y) merge(x, y, by="id", all=TRUE), subsets)
        rownames(prediction) <- prediction$id
        prediction <- prediction[colnames(prediction) != "id",]
                             }                       
    prediction_binary <- apply(prediction, 1, summarize_binary_prediction)

    df <- data.frame(id=rownames(prediction), method=method,
                     set=set, prediction_binary = prediction_binary)   
    return(df)
}
                             
get_bootstrap_df <- function(pattern, path,  size=20, sequence = seq(1,20,1)){
    files <- list.files(paste(sep="/", path, "ItClust"), pattern = pattern)
    combinations <- lapply(sequence, function(seed) get_combinations(files, size, seed))
      print("Get CellID...")
    sets_cellid <- lapply(files, function(file) readCellID(path, file))
    names(sets_cellid) <- files
    cellid_pred <- do.call(rbind,lapply(sequence, function(seq) get_predicton(sets_cellid,
                                                                       unlist(combinations[seq]),
                                                                       "CellID", seq, cutoff)))
                             
    print("Get SCN...")                                    
    files_scn <- lapply(files, function(file) list.files(paste(sep="/", path, "SingleCellNet"),
                                                     pattern = paste0(file, ".txt"),
                                                     full.names = T))
    sets_scn  <- lapply(files_scn, function(file) readSCN(file))
    names(sets_scn)<- files
    scn_pred <- do.call(rbind,lapply(sequence, function(seq) get_predicton(sets_scn,
                                                                       unlist(combinations[seq]),"SCN", seq)))                      
 
    print("Get Seurat...")
    files_seurat <- lapply(files, function(file) list.files(paste(sep="/", path, "Seurat"),
                                                            pattern = paste0(file, "_predictions.txt"),
                                                            full.names = T))
              
    sets_seurat  <- lapply(files_seurat, function(file) readSeurat(file))
    names(sets_seurat)<- files
    seurat_pred <- do.call(rbind,lapply(sequence, function(seq) get_predicton(sets_seurat,
                                                                       unlist(combinations[seq]),
                                                                       "Seurat", seq, cutoff)))
                               
    

      print("Get SingleR")
    files_singler <- lapply(files, function(file) list.files(paste(sep="/", path, "SingleR"),
                                                         pattern = paste0(file, ".txt"), 
                                                         full.names = T))
    sets_singler  <- lapply(files_singler, function(file) readSingleR(file))
    names(sets_singler )<- files
    singler_pred <- do.call(rbind,lapply(sequence, function(seq) get_predicton(sets_singler,
                                                                       unlist(combinations[seq]),
                                                                       "SingleR", seq)))
  print("Get ItClust...")
files_itclust <- paste(sep="/", path, "ItClust", files)
sets_itclust  <- lapply(files_itclust, function(file) readItclust(file))
names(sets_itclust )<- files
itclust_pred <- do.call(rbind,lapply(sequence, function(seq) get_predicton(sets_itclust,
                                                                       unlist(combinations[seq]),
                                                                       "ItClust", seq)))
itclust_pred$id <- stringr::str_replace(itclust_pred$id, "-target", "")
  predictions <- do.call(rbind, list(seurat_pred, scn_pred, singler_pred, itclust_pred, cellid_pred))
  return(predictions)                      
 
}
                                         
                                                                                                           
args = commandArgs(trailingOnly=TRUE)
folder <- args[1]
name <- args[2]
meta <- read.csv(args[3]) 

sizes <- c(38, 100, 1000)
nr_genes <- c(200,1000,2000)
                                     
x <- expand.grid(list(name, sizes, seq(1,20,1)))

tags_all <- list.files(paste(sep="/", folder, "ItClust"), pattern=name)

tags <- apply(x, 1, function(row) paste(collapse ="_", unlist(row), sep=""))
tags <- stringr::str_replace_all(tags, " ", "")

tags <- tags[tags %in% tags_all]

data <- do.call(rbind,lapply(tags, function(tag) readConfidence_Itclust(paste(sep="/", folder, "ItClust"), tag, meta)))
print(head(data))

    seurat <- do.call(rbind,lapply(tags, function(tag) readConfidence_Seurat(paste(sep="/", folder, "Seurat"), paste(sep="_",tag, 1000))))

    scn <- do.call(rbind,lapply(tags, function(tag) readConfidence_SCN(paste(sep="/", folder, "SingleCellNet"), paste(sep="_",tag,1000))))

    singleR <- do.call(rbind,lapply(tags, function(tag) readConfidence_SingleR(paste(sep="/", folder, "SingleR"), paste(sep="_",tag, 1000))))
                                    
    data <- rbind(data, seurat, scn, singleR) 
                                                                         
write.csv(data, args[4], col.names = T, row.names=F, sep=",", quote=F)