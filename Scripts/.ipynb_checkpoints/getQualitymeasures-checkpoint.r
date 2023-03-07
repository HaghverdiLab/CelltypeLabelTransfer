data <- read.csv("/fast/AG_Haghverdi/Carla_Moelbert/Celltype_annotation//Results/Files/results_PBMC10x_long.csv")

get_measures <- function(data, type, ref, method, size, set){
    data <- data[data$Reference == ref & data$Approach == method & data$Size == size & data$Set == set,] #
    tp <- length(data$Prediction[data$Prediction == type & data$class == type])
    fp <- length(data$Prediction[data$Prediction == type & data$class != type])
    fn <- length(data$Prediction[data$Prediction != type & data$class == type])
    tn <- length(data$Prediction[data$Prediction != type & data$class != type])
    precision <- tp / (tp + fp)
    recall <- tp / (tp + fn)
    f1 <- 2*(precision * recall) / (precision + recall)
    accuracy <- (tp) / length(data$Prediction[data$class == type])
    return(data.frame("class"=type,"reference"=ref,"method"=method,"size"=size,"set"=set,
                      "precision"=precision,"recall"=recall,"f1"=f1, "accuracy"=accuracy))
}

celltypes = c("Cytotoxic T cell", "CD4+ T cell", "CD14+ monocyte", "B cell", "Megakaryocyte", "Natural killer cell",
              "CD16+ monocyte", "Dendritic cell", "Plasmacytoid dendritic cell")
methods <- c("Seurat",  "SingleR","CellID", "SCN", "ItClust")

measures <- do.call(rbind, lapply(celltypes,         function(type) 
            do.call(rbind, lapply(methods,           function(method) 
            do.call(rbind, lapply(unique(data$Size), function(size)
            do.call(rbind, lapply(unique(data$Set),  function(set) get_measures(data, type, "PBMC10x", method, size,set)))))
                                                                                 
                                                                                
write.table(measures, "/fast/AG_Haghverdi/Carla_Moelbert/Celltype_annotation/Results/Files/results_PBMC10x_qualitymeasures.csv",
            row.names=F, col.names=F, quote=F, sep=",")  