source("/fast/AG_Haghverdi/Carla_Moelbert/Celltype_annotation/Scripts/functions.R")

if (!require("Seurat")) install.packages("Seurat")
library("Seurat")


predictCelltype_seurat <- function(dataTrain, dataTest, metaTrain, metaTest,
                                   output, path, features=NULL, hvgs=200){
  
  print(dataTrain)
  train <- train <- prepareData(dataTrain, metaTrain,"train", features=features, hvgs=hvgs)
  
  print(dataTest)
  test <- prepareData(dataTest, metaTest, "test", features=rownames(GetAssayData(train)))
    
  print("Scale data")
  #train <- ScaleData(train, verbose = FALSE)
  print(train)
    
  print("Get anchors")
  target.anchors <- FindTransferAnchors(reference = train, query = test,   
                                        features=rownames(GetAssayData(train)),
                                        dims = 1:30, verbose = F)
  print("Transfer data")
  
  predictions <- TransferData(anchorset = target.anchors,
                              refdata = train$class_, dims = 1:30, verbose = F)
  predictions$class_<- test$class_
  print(head(predictions[c("predicted.id", "class_")]))
  predictions$prediction.match <- predictions$predicted.id == predictions$class_
  
  write.table(predictions, paste0( path,"_predictions.txt") ,sep="\t", quote=F)
  target.copy <- test
  target.copy <- AddMetaData(target.copy, metadata = predictions)
  target.copy$prediction.match <- target.copy$predicted.id == target.copy$class_
  
  #print("Write data")
    
  data <- data.frame(target.copy$id, target.copy$class_, target.copy$predicted.id,
                         target.copy$prediction.match)
  colnames(data) <- c("id", "class_", "predicted", "prediction.match")
  write.table(data, paste0(path,".txt") ,sep="\t", quote=F)
  print("Done")
}


#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
print(length(args))
input <- args[1]
test <- args[2]
hvgs=200
if(length(args)== 4){
    if(args[4] == "all"){
        hvgs = NULL
        features= NULL
    } else if(is.na(as.integer(args[4]))) features <- read.csv(args[4])$features
    else {
        hvgs= as.numeric(args[4])
        features= NULL
    }
} else features <- NULL 

output = args[3]
if (length(args)== 4) output <- paste(sep="_", output, hvgs) 
print(output)
predictCelltype_seurat(paste(sep="/", input, "data_train.txt"),
                       paste(sep="/", test, "data.csv"),
                       paste(sep="/", input, "meta_train.txt"),
                       paste(sep="/", test, "meta.csv"),
                       basename(input), output, features, hvgs)

