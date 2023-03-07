addFormatting <- function(plot, ylab="", xlab="", legend="bottom", xtext="none",
                          celltypes=NULL, ylim = c(0,1), size= 8, keysize = 0.5){
    if(!is.null(celltypes)) plot <- plot + scale_fill_manual(values=celltypes) +
                                           scale_color_manual(values=celltypes)
    plot <- plot + 
            theme_bw() + ylab(ylab) + xlab(xlab) + 
            theme(axis.text=element_text(size=size), axis.title=element_text(size=size),
                 legend.position = legend, legend.key.size = unit(keysize, 'lines'),
                 legend.text=element_text(size = size),
                 strip.text = element_text(size = size), title=element_text(size = size)) +
            ylim(ylim) 
    
    if(xtext=="none"){
        plot <- plot + theme( axis.text.x = element_blank(),
                             axis.ticks.x = element_blank())
    } else if (xtext == "angle") {
        plot <- plot + theme(axis.text.x=element_text(angle=45, hjust = 1))
    }

    return(plot)
}

################################ UMAPS ################################

prepare_umap <- function(file,  meta_data, normalization_method="L1",
                         min_cells=3, hvgs=500, split=NULL){
    set.seed(0)
    data <- readData(file)
  
    #data <- data[, colnames(data) %in% rownames(meta_data)]
    meta_data <- meta_data[rownames(meta_data)%in% colnames(data),]
    data <- data[!(rowSums(data != 0) < min_cells),]

    meta_data <- meta_data[meta_data$id %in% colnames(data),]
    rownames(meta_data) <- meta_data$id
    print(hvgs)
    data <- select_hvg(data, hvgs)
    if(!is.null(split)){
      sets <- lapply(unique(meta_data[, split]), 
                     function(s) data[, colnames(data) %in% meta_data$id[meta_data[,split] == s]])
      data <- batchelor::fastMNN(sets)
      print("....")
      data <-  SummarizedExperiment::assay(data, "reconstructed")
    } else{

        # Cell wise (Result: Sum of each row = 1)
        if(normalization_method == "L1"){
            data <- apply(data, 2, normalize_l1)
        }else if(normalization_method == "L2") data <- apply(data, 2, normalize_l2)

    }

    print("Scaling...")
    # Gene wise (Variance of each column is 1, mean = 0) (Colum should be Gene)
    data <- t(scale(t(data), center = TRUE, scale = TRUE))
    
    print("create seurat object....")   
    pbmc <- Seurat::CreateSeuratObject(as.matrix(data), meta.data = meta_data)

    print("Set assay....")                 
    pbmc <- Seurat::SetAssayData(pbmc, assay = "RNA", slot = "scale.data",
                                 new.data = as.matrix(data))
    
    print("PCA....")
    pbmc <- Seurat::RunPCA(pbmc, features = rownames(data),verbose = F)
    print("Get UMAP...")
    pbmc <- Seurat::RunUMAP(pbmc, dims = 1:10, verbose = F)
    return(pbmc)
}


plot_umap <- function(data,groups, split=NULL,color="Spectral",
                      nrow=1, ncol=1, cells=NULL, return_plotlist=FALSE, legend_ncol=3, title=""){
    umaps <- lapply(groups, function(method) Seurat::DimPlot(data,reduction="umap", split.by=split, 
                                                     group.by=method, 
                                                     cols=color[names(color) %in% data@meta.data[,method]],
                                                     cells = cells,label=F)+
                    labs(title=title)+ 
                    theme(axis.text=element_blank(), axis.title=element_text(size=8), axis.ticks=element_blank(),
                          legend.text=element_text(size=8), legend.position="right", 
                          plot.title=element_text(size=8),
                          legend.key.size = unit(0, 'lines'))+
                    guides(color = guide_legend(ncol = legend_ncol, override.aes = list(size=3))) 
                   ) 
    print("....")              
    names(umaps) <- groups
    if(return_plotlist)return(umaps)
    #plot <- ggpubr::ggarrange(plotlist=umaps, common.legend = F, nrow=nrow, ncol=ncol)
    return(umaps)
}
                     

plot_features <- function(data,groups, split=NULL,nrow=1, ncol=1, cells=NULL, return_plotlist=FALSE,
                          legend_ncol=3, title=""){

    umaps <- lapply(groups, function(method) Seurat::FeaturePlot(data, reduction="umap", 
                                                                 split.by=split,
                                                                 features=method) +
                    labs(title=title, color= "Accuracy")+ 
                    theme(axis.text=element_blank(), axis.title=element_text(size=8), axis.ticks=element_blank(),
                          plot.title=element_text(size=8),
                          legend.text=element_text(size=8), legend.title=element_text(size=8),
                          legend.key.size = unit(0, 'lines'),
                         legend.justification = "center")+
                    guides(color = guide_legend(ncol = legend_ncol, override.aes = list(size=3)))+
                    viridis::scale_color_viridis(guide = "colourbar")
                   ) 
                   
    names(umaps) <- groups
    if(return_plotlist)return(umaps)
    #plot <- ggpubr::ggarrange(plotlist=umaps, common.legend = F, nrow=nrow, ncol=ncol)
    return(umaps)
}  
                     
plot_difference <- function(data,groups, split=NULL,nrow=1, ncol=1, cells=NULL, return_plotlist=FALSE,
                          legend_ncol=3, title=""){

    umaps <- lapply(groups, function(method) Seurat::FeaturePlot(data, reduction="umap", 
                                                                 split.by=split,
                                                                 features=method) +
                    labs(title=title, color= "Difference in Accuracy")+ 
                    theme(axis.text=element_blank(), axis.title=element_text(size=8), axis.ticks=element_blank(),
                          plot.title=element_text(size=8),
                          legend.text=element_text(size=8), legend.title=element_text(size=8),
                          legend.key.size = unit(0, 'lines'),
                         legend.justification = "center")+
                    guides(color = guide_legend(ncol = legend_ncol, override.aes = list(size=3)))+ lims(color=c(-1,1))+ 
                    scale_colour_gradient2(low ="red",mid = "grey",high = "green",midpoint = 0) 
                   ) 
                   
    names(umaps) <- groups
    if(return_plotlist)return(umaps)
    #plot <- ggpubr::ggarrange(plotlist=umaps, common.legend = F, nrow=nrow, ncol=ncol)
    return(umaps)
}                    
                      
                     
                     
################################ Lineplots ################################
                        


get_data_lineplot <- function(data, maxsize=3000, groupby = c("size", "method",
                                                              "reference", "class") ){
    print(head(data))
    summary <- data %>% 
           dplyr::group_by(across(all_of(groupby))) %>% 
           dplyr::summarize(mean_accuracy = mean(accuracy),
                            p25_accuracy = quantile(accuracy, probs = c(0.25)),
                            p75_accuracy = quantile(accuracy, probs = c(0.75)),
                            mean_precision = mean(precision),
                            mean_f1 = mean(f1),
                            p25_f1 = quantile(f1, probs = c(0.25)),
                            p75_f1 = quantile(f1, probs = c(0.75)),
                            mean_precision = mean(precision),
                            p25_precision = quantile(precision, probs = c(0.25)),
                            p75_precision = quantile(precision, probs = c(0.75))) 
    
    return(summary)
}

get_lineplot <- function(data, value, log=TRUE, breaks=c(38,100, 500, 1000, 3000) ){
    if(value == "Accuracy"){
        plot = ggplot(data, aes(size, mean_accuracy, group=method)) +
           geom_ribbon(aes(ymin = p25_accuracy, ymax = p75_accuracy),
                       alpha = 0.3)+
           geom_hline(aes(yintercept = full_accuracy), color="red",
                      linetype="dashed")+
           geom_errorbar(aes(ymin=p25_accuracy, ymax=p75_accuracy), width=.5) 
    } else if(value=="Precision") {
        plot = ggplot(data, aes(size, mean_precision, group=method)) +
               geom_ribbon(aes(ymin = p25_precision, ymax = p75_precision),
                       alpha = 0.3)+
           geom_hline(aes(yintercept = full_precision), color="red",
                      linetype="dashed")+
           geom_errorbar(aes(ymin=p25_precision, ymax=p75_precision), width=.5) 
    } else if (value == "F1"){
        plot = ggplot(data, aes(size, mean_f1, group=method)) +
               geom_ribbon(aes(ymin = p25_f1, ymax = p75_f1),
                       alpha = 0.3)+
           geom_hline(aes(yintercept = full_f1), color="red",
                      linetype="dashed")+
           geom_errorbar(aes(ymin=p25_f1, ymax=p75_f1), width=.5) 
    }
    
    plot <- plot + geom_point(size=0.5) + geom_line()+ 
            geom_vline(aes(xintercept=refSize), linetype="dotted",
                       color="black")+
            facet_grid(rows = vars(class), cols =vars(method), 
                       labeller = label_wrap_gen(width=10)) + 
            scale_y_continuous(breaks=c(0,0.5,1))
    
    plot <- addFormatting(plot, value, "Cells per cell type", xtext="angle", )
    
    if(log == TRUE) scale_x_continuous(trans='log2', breaks=breaks)
    else scale_x_continuous(breaks=breaks)
    return(plot)
}
                                      
                        
################################# Other Plots ###########################################

plot_confidence_scores <- function(data){
   plot <- ggplot(data, aes(as.factor(size), score, color=match)) + 
    geom_boxplot(alpha=0.3, width=0.5, position = position_dodge(0.5),  outlier.size=0.1)+
    facet_grid(cols=vars(method), rows=vars(class),
               labeller = label_wrap_gen(width=15))+
    theme_bw()+
    theme(legend.position = "bottom", legend.text=element_text(size = 8),
          legend.title=element_text(size = 8), legend.key.size = unit(0, 'lines'),
          axis.text=element_text(size = 8), axis.title=element_text(size = 8),
          strip.text = element_text(size = 8)) +
    xlab("Cells per cell type")+ylab("Confidence Scores")+
    labs(color="Predictions")+ylim(0,1)+
    scale_y_continuous(breaks=c(0,0.5,1))+
    stat_compare_means(aes(group = match), label = "p.signif", hide.ns = TRUE, label.y = 0.9,
                      method = "t.test",
                       symnum.args = list(cutpoints = c(0,0.05, 1), symbols = c("*", "ns")),
                       method.args = list(alternative = "two.sided",var.equal=FALSE)
                       )
    return(plot)
}
                        
get_violin_plot <- function(data, colors, celltypes, methods, title, value, color){
    data$class <- factor(data$class, levels=celltypes)
    data$method <- factor(data$method, levels=methods)
    
    if(value =="Accuracy"){
      plot <- ggplot(data, aes(class, accuracy, color=class, fill=class)) +
              geom_boxplot(alpha=0.3)+
              geom_point(aes(y=full_accuracy), color="black",
                         size= 0.5, show.legend = F)
    } else if (value == "Precision"){
      plot <- ggplot(data, aes(class, precision, color=class, fill=class)) + 
              geom_boxplot(alpha=0.3)+
              geom_point(aes(y=full_precision), color="black", size= 0.5,
                         show.legend = F)
    } else if (value == "F1"){
        plot <- ggplot(data, aes(class, f1, color=class, fill=class)) + 
              geom_boxplot(alpha=0.3)+
              geom_point(aes(y=full_f1), color="black", size= 0.5,
                         show.legend = F)
    }
    plot <- plot + 
            facet_wrap(facets=vars(method),labeller = label_wrap_gen(width=12),
                       nrow = 1) + labs(color="", fill="", title=title)
    plot <- addFormatting(plot, value, "Cell type", "bottom", "none", color)
    return(plot)
}
                     
get_plot_comparison <-  function(data, value, colors){
    if(value == "Accuracy"){
        plot <- ggplot(data, aes(class, accuracy, color=class, size=as.factor(genes), #shape=set, 
                                  group=genes))
    } else if (value == "Precision"){plot <-  ggplot(data, aes(class, precision,color=class,
                                                              size=as.factor(genes), #shape=set,
                                      group=genes))
    }else{plot <-  ggplot(data, aes(class, f1,color=class, size=as.factor(genes), #shape=set,
                                      group=genes))
    }
    
    plot <- plot + geom_point(alpha=0.5) + #size=2
            labs(color="Celltype", size= "Gene set used") + 
            facet_grid(cols = vars(method)) + scale_size_manual(values=c(2,4,6))
     plot <- addFormatting(plot, ylab=value, xlab="Cell types", legend="left", xtext="none",
                           celltypes= colors) 

}   