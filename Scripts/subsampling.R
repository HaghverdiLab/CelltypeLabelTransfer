source("/fast/AG_Haghverdi/Carla_Moelbert/Celltype_annotation/Scripts/preparation_functions.r")
source("/fast/AG_Haghverdi/Carla_Moelbert/Celltype_annotation/Scripts/functions.R")
source("/fast/AG_Haghverdi/Carla_Moelbert/Celltype_annotation/Scripts/weighted_bootstrapping.R")

print("----------------------------------------------------")
args = commandArgs(trailingOnly=TRUE)
                               
input = args[1]
output = args[2]
nr = args[3]
steps = as.numeric(unlist(stringr::str_split(args[4], ",")))

name =  unlist(stringr::str_split(basename(input), "_"))[1]

if(!dir.exists(output)) dir.create(output) 
print(input)
meta <- read.csv(paste(sep="/", input, "meta.csv"))
print(head(meta))
data <- readData(paste(sep="/", input, "data.csv")) 

ids <- Reduce(intersect,list(meta$id, colnames(data)))

meta <- meta[meta$id  %in% ids,]    
   
         
full <- write_full(data, meta, name, output) 
sets <- seq(1,nr,1)
x <- lapply(steps,function(step) lapply(sets, function(set) write_random(data, meta, step, set,name, output)))

                                   
