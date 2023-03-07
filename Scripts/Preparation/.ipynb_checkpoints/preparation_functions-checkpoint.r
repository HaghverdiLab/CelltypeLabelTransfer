######## Data preparation - Construction of Training and Validation sets ######## 

#' Match expression matrix and meta data
#' 
#' Reduces the expression matrix to cells that are also included in the meta data
#' @param data expression matrix with cells in the columns
#' @param meta meta data for the expression matrix 
#' @param out_data file in which to save the resulting expression matrix
#' @param out_meta file in which to save the resulting meta data
get_subset <- function(data, meta, out_data, out_meta){
  
  data <- data[, colnames(data) %in% meta$id]
  print(paste(ncol(data), nrow(data)))
  
  write.table(data, out_data, row.names = T,
              col.names = T, quote = F, sep=",")
  print(nrow(meta))
  meta <- meta[meta$id %in% colnames(data),]
  write.table(meta, out_meta, row.names = F,
              col.names = T, quote = F, sep=",")
  print(nrow(meta))
}


                                                                  
                               
get_validation_set <- function(full, meta, train, output){
  data_train <- read.csv(train)
  meta <- meta[!(meta$id %in% data_train$id),]
  data <- full[,(colnames(full) %in% meta$id)]
  if(!dir.exists(output))dir.create(output)
  write.table(meta, paste(output, "meta_test.txt", sep="/"), col.names = T,
              row.names = T, sep=",", quote = F)
  write.table(as.matrix(data), paste(output, "data_test.txt", sep="/"),
              col.names = T, row.names = T, sep=",", quote = F)
}

addX <- function(name){
  if (startsWith(name, "4")) return(paste0("X",name))
  else return(name)
}
updateID <- function(idVector){
  idVector <- gsub("-", ".", idVector)
  idVector <- unlist(lapply(idVector, function(n) addX(n)))
  return(idVector)
}
getMetaFormat <- function(meta){
  meta <- meta[,c("celltype", "index")]
  colnames(meta) <- c("class_", "id")
  meta$id <- updateID(meta$id)
  return(meta)
}
