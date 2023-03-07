#BiocManager::install("CelliD")
library("CelliD")
source("/fast/AG_Haghverdi/Carla_Moelbert/Celltype_annotation/Scripts/functions.R")

args = commandArgs(trailingOnly=TRUE)
folder <- args[1]
tag <- args[2]
testfile <- args[3]
out <- args[4]
hvgs=200
if(length(args)== 5){
    print(args[5])
    if(args[5] == "all"){
        features= NULL
    }else if(is.na(as.integer(args[5]))) features <- read.csv(args[5])$features
    else {
        hvgs= as.numeric(args[5])
        features= NULL
    }
} else features <- NULL 

input <- paste(sep="/", folder, tag)
output <- paste(sep="/", out,tag)
if (length(args)== 5)output <- paste(sep="_", output, hvgs)
print("Prepare training data")


print("Start....")
train <- prepareData(paste(sep="/", input, "data_train.txt"),
                     paste(sep="/", input, "meta_train.txt"),
                    "train", hvgs=5000, features=features)

print("1. ---------------")    
# Library-size normalization, log-transformation, and centering and scaling of gene expression values
train <- RunMCA(train)

print("2. ---------------")   
if (is.null(hvgs)) hvgs <- nrow(train)
print(hvgs)
train_cell_gs <- GetCellGeneSet(train, dims = 1:50, n.features = hvgs)
train_group_gs <- GetGroupGeneSet(train, dims = 1:50, n.features = hvgs, group.by = "class_")



print("Prepare test data")
test <- prepareData(paste(sep="/", testfile, "data.csv"),
                        paste(sep="/", testfile, "meta.csv"), "test",
                        features=rownames(GetAssayData(train)))

print(".........")
test <- RunMCA(test, nmcs = 50)

print("Cell type prediction")
predicted_cell_gs <- RunCellHGT(test, pathways = train_group_gs, dims = 1:50)
predicted_cell_gs_match <- rownames(predicted_cell_gs)[apply(predicted_cell_gs, 2, which.max)]


#predicted_cell_gs_prediction <- train$class_[predicted_cell_gs_match]

predicted_cell_gs_prediction_signif <- ifelse(apply(predicted_cell_gs, 2, max)>2, yes =predicted_cell_gs_match, "unassigned")
test$predicted <- predicted_cell_gs_prediction_signif

print("Save results")
test$prediction.match <- test$predicted == test$class_


data <- data.frame(test$id, test$class_, test$predicted,
                         test$prediction.match)
colnames(data) <- c("id","class_", "predicted", "prediction.match")
print(head(data))
write.table(data, paste0(output,".txt") ,sep="\t", quote=F)
print("Done")
