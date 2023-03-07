library(singleCellNet)
source("/fast/AG_Haghverdi/Carla_Moelbert/Celltype_annotation/Scripts/functions.R")


args = commandArgs(trailingOnly=TRUE)
print(args)
folder <- args[1]
tag <- args[2]
testfolder <- args[3]
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

set.seed(100) #can be any random seed number

print("Get training data...")
data <- readData(paste(sep="/", input, "data_train.txt"))
data <- preprocessing(data, hvgs=hvgs, features=features)

meta <- read.csv(paste(sep="/", input, "meta_train.txt"))

print("Get test data...")
test <- readData(paste(sep="/", testfolder, "data.csv"))
test <- preprocessing(test, features=rownames(data))
meta_test <- read.csv(paste(sep="/", testfolder, "meta.csv")) 
print(head(meta_test))

print("1. ------------")
ncells <- as.numeric(stringr::str_split(tag,pattern="_", simplify = T)[2])

x<- table(meta$class_)
meta <- meta[meta$class_ %in% names(x)[x >= 3],]

stList = splitCommon(sampTab=meta, ncells=ncells, dLevel="class_")

meta = stList[[1]]
data = data[,meta$id]

print("2. -------------------")
class_info<-scn_train(stTrain = meta, expTrain = as.matrix(data), nTopGenes = 10, nRand = min(100,nrow(meta)),
                      nTrees = 1000,  nTopGenePairs = 25, dLevel = "class_", colName_samp = "id")

print("3. -------------------")
nqRand = 50
print(head(class_info[['cnProc']]))
predictions <- as.data.frame(t(scn_predict(class_info[['cnProc']], test, nrand=nqRand)))
print("4. -------------------")
predicted <- unlist(colnames(predictions)[apply(predictions,1,which.max)])


predictions$id = rownames(predictions)
predictions$predicted = predicted
print("5. -------------------")

results <- merge(predictions, meta_test, by="id", all = T)
print(head(results))
results$prediction.match <- results$predicted == results$class_
write.table(results, paste0(output,".txt") ,sep="\t", quote=F)
