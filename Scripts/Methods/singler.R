source("/fast/AG_Haghverdi/Carla_Moelbert/Celltype_annotation/Scripts/functions.R")


args = commandArgs(trailingOnly=TRUE)
folder <- args[1]
tag <- args[2]
testfile <- args[3]
out <- args[4]
hvgs=200
if(length(args)== 5){
    if(args[5] == "all"){
        hvgs = NULL
        features= NULL
    }else if(is.na(as.integer(args[5]))) features <- read.csv(args[5])$features
    else {
        hvgs= as.numeric(args[5])
        features= NULL
    }
} else features <- NULL


print("Starting....")
input <- paste(sep="/", folder, tag)
output <- paste(sep="/", out,tag)

if (length(args)== 5)output <- paste(sep="_", output, hvgs)
print("Get training data....")
train <- getExperiment(paste(sep="/", input, "data_train.txt"),
                     paste(sep="/", input, "meta_train.txt"), features=features, hvgs=hvgs)

print("Get Query data....")
print(head(rownames(train)))
test <- getExperiment(paste(sep="/", testfile, "data.csv"),
                        paste(sep="/", testfile, "meta.csv"), 
                        features=rownames(train))

print("-------------------------------------------------------------------------")

print("Run single R")

predictions<- as.data.frame(SingleR::SingleR(test=test, ref=train,labels=train$label, de.method="wilcox"))



#print(head(predictions))
predictions$id <- rownames(predictions)




meta_test <- read.csv(paste(sep="/", testfile, "meta.csv"))
results <- merge(meta_test, predictions, by="id")


colnames(results)[colnames(results)== "labels"] <- "predicted"

results$prediction.match <- results$predicted == results$class_

print(table(results$prediction.match))
write.table(results, paste0(output,".txt") ,sep="\t", quote=F)
