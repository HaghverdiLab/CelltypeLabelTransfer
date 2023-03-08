################## Read the input files ##################

read_data_results <- function(file, path, method){
    data <- read.csv(paste(path,file, sep="/"), sep="\t")
    data$tag <- stringr::str_replace(file, ".txt", "")
    data$method <- method
    data <- data[, c("id", "predicted", "tag", "method")]
    return(data)
}

get_data <- function(path, method, pattern){
   files <- list.files(path, full.names = F, pattern= pattern)
   if(method =="Seurat") files <- files[stringr::str_detect(files, "pred",
                                                                          negate = T)]
   seurat <- do.call(rbind, lapply(files,
                                   function(file) read_data_results(file, path, method)))
   return(seurat)
 
}
                                  
getItClust <- function(file, path, method){
    data <- read.csv(paste(path,file, "results.txt", sep="/"), sep="\t")
    data$tag <- stringr::str_replace(file, ".txt", "")
    data$method <- method
    data <- data[, c("cell_id", "predicted_celltype", "tag", "method")]
    colnames(data) <- c("id", "predicted", "tag", "method")
    return(data)
}


                                
################## Get the summary files ##################                                  
get_summary <- function(folder, meta, outfolder, query, name, celltypes, methods, seqs, full,
                        maxsize=NULL, weighted=TRUE, pattern=""){
    print("Get predictions....")
    data <- get_predictions(folder, meta, methods, maxsize, weighted = weighted, pattern=pattern)
 
    print("Merge with query...")
    data <- merge(data, query, by="id")

    print("Get measures...")
    if(length(data$prediction[is.na(data$prediciton)]) > 0 ) data$prediction[is.na(data$prediciton)]<- "unassigned"
    
    data_set <- do.call(rbind, lapply(celltypes, function(type) 
            do.call(rbind, lapply(methods,   function(method) 
            do.call(rbind, lapply(seqs,function(set)get_measures_set(data, type, pattern,
                                                                             method, set)))))))
    print("Merge with full...")
    print(head(full))
    data_set <-  merge(data_set, full, by=c("class", "method"))
    print(head(data_set))
                                  
    print("Get measures id...")                         
    data_id <- get_accuracy_umap(data) 
    print("Write files...")
    if(!is.null(maxsize)) name <- paste0(name, maxsize) 
                                  
    if(weighted == FALSE) name <- paste0(name, "unweighted") 
    write.table(data_set, paste(sep="/", outfolder, paste0("summary_", name, ".csv" )),
                col.names=T, row.names=T, quote=T, sep=",")
    return(data_id)
}

                               

                   
   
