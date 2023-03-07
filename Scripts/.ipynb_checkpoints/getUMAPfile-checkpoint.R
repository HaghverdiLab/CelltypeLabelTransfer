getName <- function(folder, position){
  name <- unlist(stringr::str_split(folder,"/"))
  print(name)
  name <- name[length(name)]
  print(name)
  name <- unlist(stringr::str_split(name,"_"))[position]
  return(name)
}
print("Start....")
args = commandArgs(trailingOnly=TRUE)
method = args[1] 
folder = args[2]
output = args[3]
if(length(args) == 4) metafile = args[4]

if(method %in% c("Seurat", "SCN", "CellID", "SingleR")){
  if(stringr::str_detect(folder, "predict")) stop("Wrong file")
  print(method)
  name <- paste(sep="_", getName(folder,1), method ,getName(folder,2),
                stringr::str_replace(getName(folder,3), ".txt", "") )
  df <- read.csv(folder, sep="\t")
    print(head(df))
  df <- df[,c("id", "predicted")] #.match
  colnames(df) <- c("id", name)
    
} else if (method == "ItClust"){
  name <- paste(sep="_",  getName(folder,1), "ItClust",getName(folder,2),
                getName(folder,3))
  
  df <- read.csv(paste(folder, "results.txt", sep="/"), sep="\t")
  df <- df[,c("class_", "predicted_celltype", "cell_id")]
  print(head(df))
  colnames(df) <- c("class", "predicted", "id")
  df[,name] <- df$predicted # df$class == 
  df <- df[, c("id",  name)]
  
} else if (method == "MLP"){
  if(is.null(metafile)) return(NULL)
  name <- paste(sep="_",  getName(folder,1),"MLP", getName(folder,2),
                getName(folder,3))
  df <- read.csv(metafile)
  if(colnames(df)[1] == "id") colnames(df) <- c("id","class")
  if(colnames(df)[2] == "id") colnames(df) <- c("class","id")  
  pred <- read.csv(paste(folder, "predictions.csv", sep="/"), header=F)
  df[,name] <-   pred$V2 # pred$V1 ==
  df <- df[, c("id", name)]
}

if(file.exists(output)){
  print("Compare files")
  x <- as.data.frame(data.table::fread(output,sep = ",", verbose = F))
  print(colnames(df))
  x <- x[x$id %in% df$id,]
  df <- df[df$id %in% x$id,]  
  if(all(x$id == df$id)) x[,name] <- df[,name] 
  else x <-  merge(x, df, by=c("id"),  all=T)
  print("Write file...")
  data.table::fwrite(x, output, row.names = FALSE, col.names = TRUE, append = FALSE, quote=F, sep=",") 
} else write.table(df, output, append = F, row.names = F, col.names = T, quote = F,
                   sep=",")
print("done")
