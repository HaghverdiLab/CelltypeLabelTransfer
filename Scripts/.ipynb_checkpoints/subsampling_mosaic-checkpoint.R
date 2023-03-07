source("/fast/AG_Haghverdi/Carla_Moelbert/Celltype_annotation/Scripts/preparation_functions.r")
source("/fast/AG_Haghverdi/Carla_Moelbert/Celltype_annotation/Scripts/functions.R")

args = commandArgs(trailingOnly=TRUE)

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
  print(paste(nr, type, seed))
  df <- data[data$class_ == type,]
  if(nrow(df)<=nr)return(df)  
  set.seed(seed)
  df <- df[sample(nrow(df), nr),]
  return(df)
} 
                               
get_random_mosaic <- function(data, nr,type, seed, col, x){
  print(paste(nr, type, seed))
  df <- data[data$class_ == type & data[col] == x,]
  if(nrow(df)<=nr)return(df)  
  set.seed(seed)
  df <- df[sample(nrow(df), nr),]
  return(

write_full <- function(data, meta, name, output){
    
    
    out <- paste(sep="/", output, paste(sep="_",name, max(table(meta$class_)), "0"))
    if(!dir.exists(out)) dir.create(out)

    write.table(meta, paste(sep="/", out, "meta_train.txt"), sep=",",
              row.names = F, col.names = T, quote = F)
  
    data_sub <- data[,colnames(data) %in% meta$id]
                               
    write.table(data_sub, paste(sep="/", out, "data_train.txt"), sep=",",
              row.names = T, col.names = T, quote = F)  

} 

                               
input = args[1]
output = args[2]
nr = args[3]
steps = as.numeric(unlist(stringr::str_split(args[4], ",")))

name =  unlist(stringr::str_split(basename(input), "_"))[1]
name <- paste0(name, "2")
if(!dir.exists(output)) dir.create(output) 
meta <- read.csv(paste(sep="/", input, "meta.csv"))
data <- readData(paste(sep="/", input, "data.csv")) 

if(min(table(meta$class_)) > 10) steps <- append(steps, min(table(meta$class_)))
steps <- steps[steps <= max(table(meta$class_))]
print(paste("Steps:", steps))
print("----------------------------------")
 
print(data[1:3,1:5])
ids <- Reduce(intersect,list(meta$id, colnames(data)))
print(length(ids))
meta <- meta[meta$id  %in% ids,]    
   
         
full <- write_full(data, meta, name, output) 
sets <- seq(1,nr,1)
x <- lapply(steps,function(step) lapply(sets, function(set) write_random(data, meta, step, set,name, output)))

                                   
