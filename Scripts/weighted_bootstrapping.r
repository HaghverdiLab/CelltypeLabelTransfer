library(dplyr)
source("/fast/AG_Haghverdi/Carla_Moelbert/Celltype_annotation/Scripts/functions.R")
#######################  Subsampling functions ###########################
write_random <- function(data, meta, step, set, name, output){
    out <- paste(sep="/", output, paste(sep="_",name, step, set))
    if(!dir.exists(out)) dir.create(out)
    meta_sub <- do.call(rbind,
                        lapply(unique(meta$class_),
                               function(type) get_random(meta, step, type,set)))
   
    write.table(meta_sub, paste(sep="/", out, "meta_train.txt"), sep=",",
              row.names = F, col.names = T, quote = F)
  
    data_sub <- data[,colnames(data) %in% meta_sub$id]                          
    write.table(data_sub, paste(sep="/", out, "data_train.txt"), sep=",",
              row.names = T, col.names = T, quote = F)  

}           

get_random <- function(data, nr ,type, seed){
  df <- data[data$class_ == type,]
  if(nrow(df)<=nr)return(df)  
  set.seed(seed)
  df <- df[sample(nrow(df), nr),]
  return(df)
} 
                               
write_full <- function(data, meta, name, output){
    out <- paste(sep="/", output,
                 paste(sep="_",name, max(table(meta$class_)), "0"))
    
    if(!dir.exists(out)) dir.create(out)
    
    write.table(meta, paste(sep="/", out, "meta_train.txt"), sep=",",
              row.names = F, col.names = T, quote = F)
    
    data_sub <- data[,colnames(data) %in% meta$id]
                               
    write.table(data_sub, paste(sep="/", out, "data_train.txt"), sep=",",
              row.names = T, col.names = T, quote = F)  
} 
                               
samples_wb <- function(datafile, metafile, name, nr_sets=1, output="",
                       maxvalue=NULL){
    if(output != "" & !dir.exists(output)) dir.create(output) 

    meta <- read.csv(metafile)
    data <- readData(datafile)
    steps <- table(meta$class_)
    if(!is.null(maxvalue)){
    steps <- steps[steps < maxvalue]
    steps <- append(steps, maxvalue)}
    
    x <- lapply(steps,
                function(step) lapply(seq(1,nr_sets,1),
                function(set) write_random(data, meta, step, set,name, output)))

}

#######################  Bootstrapping functions ###########################                               
weighted_bootstrap <- function(ref, data, weighted=TRUE){
    print("Start weighted boostrapping.....")
    n <- data.frame(table(ref$class_))

    data$score <- n$Freq[match( data$predicted, n$Var1)]
    if(weighted== TRUE)data$score <- 1 / (abs(data$score - as.integer(data$size)) + 1)
    data <- data[data$set != 0,]
    
    summary <- data %>% 
           dplyr::group_by(id, predicted, set, method, reference, genes) %>% 
           dplyr::summarize(score = sum(score)) 
    summary$predicted[is.na(summary$predicted)]<- "unassigned"
    
    df.wide <- tidyr::pivot_wider(summary, names_from = predicted, values_from = score)
    df.wide[is.na(df.wide)] <- 0

    names <-  colnames(df.wide[, 6:ncol(df.wide)])
    
    df.wide$predicted <- apply(df.wide[, 6:ncol(df.wide)], 1, which.max)
    df.wide$predicted <- names[df.wide$predicted]
    df <- df.wide[, c("id", "set", "method", "predicted", "reference", "genes")]

    return(df)
}
                                      
get_predictions <- function(path,ref,methods, maxsize=NULL, weighted=TRUE, pattern="", nr_genes=1000){
    data <- do.call(rbind,
                    lapply(methods[methods != "ItClust"],
                           function(method) get_data(paste(sep="/", path, method),
                                                    method, pattern)))
    if("ItClust" %in% methods){
        itclust_files <- list.files(paste(sep="/", path,"ItClust"), full.names = F)
        itclust_files <- itclust_files[stringr::str_detect(itclust_files, "_0", negate=T)]
        itclust <- do.call(rbind, 
                           lapply(itclust_files,
                                  function(file) getItClust(file,
                                                            paste(sep="/", 
                                                                  path,"ItClust"),
                                                            "ItClust")))
        itclust$predicted <-   translate(itclust$predicted) 
        data <- do.call(rbind,list(itclust, data))
    }
                                  
    data[c('reference', 'size', "set", "genes")] <- stringr::str_split_fixed(data$tag, '_', 4)
    data$size <- as.integer(data$size)
    data$genes[data$genes == ""] <- nr_genes                               

    if(!is.null(maxsize)) data <-  data[data$size <= maxsize,]  
                                  
    predictions <-  weighted_bootstrap(read.csv(ref), data, weighted)
    return(predictions)
}