######## General Functions ######## 
 
#' Read the expression data matrix (fast)
#' Uses fread to increase the reading spead of the expression data
#' @param data file containing the expression matrix´
readData <- function(data){
  df <-  as.data.frame(data.table::fread(data,sep = ",", verbose = F))
  rownames(df) <- df$V1
  df[,1] <- NULL
  return(df)
}

## Highly variable gene selection ##
pearson_residuals <- function(counts, theta){
    counts_sum1 = rowSums(counts)
    counts_sum0 = colSums(counts)
    counts_sum  = sum(counts)

    #get residuals
    mu = (counts_sum1  %*% t(counts_sum0)) / counts_sum
    z = (counts - mu) / sqrt(mu + mu**2/theta)

    #clip to sqrt(n)
    n = ncol(counts)
    z[z >  sqrt(n)] = sqrt(n)
    z[z < -sqrt(n)] = -sqrt(n)
    return(as.matrix(z))
}

select_hvg <- function(data, hvgs){
  residuals = pearson_residuals(data,200)
  residual_var = matrixStats::rowVars(residuals)
  names(residual_var) <- rownames(residuals)
  residual_var <- names(residual_var[order(residual_var, decreasing = T)])
  data = data[head(residual_var, hvgs), ]
  return(data)
}


## Preprocessing  the data ##
prepareData <- function(data, meta, project, hvgs=200, features=NULL){
    counts <- as.data.frame(readData(data))
    metadata <- read.csv(meta)
    metadata <- metadata[metadata$id %in% colnames(counts),]
    counts <- counts[, colnames(counts) %in% metadata$id]
    rownames(metadata)<- metadata$id
    
    counts <- preprocessing(counts, features=features, hvgs=hvgs)
    obj <-  CreateSeuratObject(counts = counts, meta.data = metadata) 

  return(obj)
}

preprocessing <- function(data, features=NULL, hvgs=200){
    if(!is.null(features)){
      data <- data[rownames(data) %in% features, ]
    } else if (!is.null(hvgs)) data <- select_hvg(data, hvgs)
    data <- log2(data + 0.001)
    data <- as.matrix(data)
    data <- data/colSums(data)[col(data)]
    return(data)
}

normalize_l1<- function(x) {
  return ((x - min(x)) / (max(x) - min(x)))
}                        
   

getFileName <- function(folder, type, set){
  name <- paste(folder, paste(type, set, sep="_"), sep="/")
  return(paste(name, "csv", sep="."))
}

getFiles <- function(datafile, metafile=NULL,set=NULL, output=NULL, features=NULL, hvgs=500){
    data <- readData(datafile)
    metadata <- read.csv(metafile)
    metadata <- metadata[metadata$id %in% colnames(data),]
    metadata <- metadata[order(metadata$class_),]
    data <- data[, colnames(data) %in% metadata$id]
    rownames(metadata)<- metadata$id
    if(!is.null(features))  data <- data[rownames(data) %in% features, ]
    # We skip the preprocessing step since ItClust has internal processing
    #data <- preprocessing(data, features=features, hvgs=hvgs)
    data <- data[order(rownames(data)),]
    
    metadata <- metadata[metadata$id %in% colnames(data),]
    data <- data[, metadata$id]

    write.table(t(data), getFileName(output, "data", set), sep=",", quote=F,
              row.names = F, col.names = T)
    
    write.table(colnames(data), getFileName(output, "cells", set), sep=",", quote=F,
                          row.names = F, col.names = F)

    write.table(metadata, getFileName(output, "meta", set), sep=",", quote=F,
              row.names = F, col.names = T)
    return(data)
}

getExperiment <- function(datafile, metafile=NULL, features=NULL, hvgs=200){
    data <- readData(datafile)
    metadata <- read.csv(metafile)
    metadata <- metadata[metadata$id %in% colnames(data),]
    data <- data[,colnames(data) %in% metadata$id]
    rownames(metadata)<- metadata$id
    data <- preprocessing(data, features=features, hvgs=hvgs)
    data <- data[, order(colnames(data))]

    experiment <- SingleCellExperiment::SingleCellExperiment(list(logcounts=data))
    metadata <- metadata[metadata$id %in% colnames(data),]
    SingleCellExperiment::colLabels(experiment) <- metadata$class_[order(metadata$id)]

  return(experiment)
}


######## Get the quality measures ######## 
calculate_measures <-  function(data, type){
    tp <- length(data$predicted[data$predicted == type & data$class == type])
    fp <- length(data$predicted[data$predicted == type & data$class != type])
    fn <- length(data$predicted[data$predicted != type & data$class == type])
    tn <- length(data$predicted[data$predicted != type & data$class != type])
    precision <- tp / (tp + fp)
    if(is.na(precision)) precision <- 0
    recall <- tp / (tp + fn)
    f1 <- 2*(precision * recall) / (precision + recall)
    if(is.na(f1)) f1 <- 0
    accuracy <- (tp) / length(data$predicted[data$class == type])
    
   return(c(precision,recall,f1, accuracy))
}

get_measures <- function(data, type, ref, method, size, set){
    
    data <- data[data$Reference == ref & data$Approach == method & data$Size == size & data$Set == set,] #
    predictions <- calculate(data, type)
   
    return(data.frame("class"=type,"reference"=ref,"method"=method,"size"=size,"set"=set, 
                      "precision"=predictions[1],"recall"=predictions[2],
                      "f1"=predictions[3], "accuracy"=predictions[4]))
}

get_measures_set <- function(data, type, ref, method, set){
    data <- data[data$ref == ref & data$method == method & data$set == set,]
    predictions <- calculate(data, type)
    return(data.frame("class"=type,"reference"=ref,"method"=method,"set"=set,
                      "precision"=predictions[1],"recall"=predictions[2],
                      "f1"=predictions[3], "accuracy"=predictions[4]))
}

get_accuracy_umap <- function(data){
    data$correct <- data$class_ == data$predicted
    data$correct[data$correct == "True"] <- 1
    summary <- data %>% 
           dplyr::group_by(id, method, reference) %>% 
           dplyr::summarize( accuracy = sum(correct) / n())


}

get_overall_measures <- function(measures, celltypes, sizes){+
    measures$obs <- sizes[ match(measures$class, celltypes ) ]

    overall_measures <-  measures %>% dplyr::group_by( method) %>% 
           dplyr::summarize(f1 = sum(f1 * obs) / sum(obs) ,
                            precision = sum(precision*obs)/ sum(obs),
                            accuracy = sum(accuracy*obs)/ sum(obs))
  
    overall_measures$method <- factor(overall_measures$method, levels=methods)
    return(overall_measures)
}


################## Translating the cell names back after ItClust changed them ##################                 
adjust_names <- function(data, name, dataset="PBMC"){
  cols <- colnames(data)[stringr::str_detect(colnames(data), name)]
  data[,cols] <- do.call(rbind,lapply(cols, function(col) translate(data[,col], dataset)))
                                                          
  return(data)                                                             
}

translate <- function(col, dataset="PBMC"){
    if (dataset == "PBMC"){
        col[col == "B"] <- "B cell"
        col[col == "Cytotoxic T"] <- "Cytotoxic T cell"
        col[col == "CD4+ T"] <- "CD4+ T cell"
        col[col == "Dendritic"] <- "Dendritic cell"
        col[col == "Natural killer"] <- "Natural killer cell"
        col[col == "Plasmacytoid dendritic"] <- "Plasmacytoid dendritic cell"
        
    } else if (dataset =="BM"){
        col[col == "B"] <- "B cell"
        col[col == "Dendritic"] <- "Dendritic cell"
        col[col == 'Eosinophil-basophil-mast progenitors'] <- 'Eosinophil-basophil-mast cell progenitors'
        col[col == 'Mesenchymal'] <- 'Mesenchymal cell'
    } 
    return(col)
} 
                                      
                                      
#################################
                                      
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
                stringr::str_replace(getName(folder,3), ".txt", ""), getName(folder,4))
  df <- read.csv(folder, sep="\t")
  df <- df[,c("id", "predicted")] #.match
  colnames(df) <- c("id", name)
    
} else if (method == "ItClust"){
  name <- paste(sep="_",  getName(folder,1), "ItClust",getName(folder,2),
                getName(folder,3), 0)
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
                   
make_long <- function(data, celltypes, methods, sizes, ids =c("id", "class_")){
    long <- reshape2::melt(data, id.vars = ids, value.name = "predicted")
    long <- long %>% tidyr::separate(variable, c('Reference', 'Approach', "Size", "Set", "Genes"))

    return(long)
}                                    