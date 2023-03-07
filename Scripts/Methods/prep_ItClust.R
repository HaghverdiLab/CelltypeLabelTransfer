library("Seurat")
source("/fast/AG_Haghverdi/Carla_Moelbert/Celltype_annotation/Scripts/functions.R")

args = commandArgs(trailingOnly=TRUE)
print(args)

folder <- args[2]
tag <- args[1]
test <- args[3]
out <- args[4]
print("Starting...")
input <- paste(sep="/", folder, tag)

output <- paste(sep="/", out,tag)


if(!dir.exists(output))dir.create(output)

train <- getFiles(paste(sep="/", input, "data_train.txt"),
                     paste(sep="/", input, "meta_train.txt"),
                    "train", output)

test <- getFiles(paste(sep="/", test, "data.csv"),
                        paste(sep="/", test, "meta.csv"), "test", output,
                        features=rownames(train))


print("Done")