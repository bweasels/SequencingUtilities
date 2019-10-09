#thanks to Hannah Thel for writing this very useful plotting function!

heatmap_plot <- function(gene_panel, p.threshold = 0.1, cluster_genes = TRUE, centered = TRUE, scaled = FALSE, horizontalHeatmaps = FALSE){
  require(pheatmap)
  require(ggplot2)
  
  #read in files, convert to matrices
  KMU_data <- read.csv('/Dropbox (Partners HealthCare)/CTC-Seq/input/2019-08-12_KMUCounts_star_genes_erc.csv', row.names = 1)
  metadata_KMU <- read.csv('/Dropbox (Partners HealthCare)/CTC-Seq/input/2019-08-12_KMUMetadata.csv', row.names = 1)
  metadata_KMU$Sample.ID <- sub("^", "X", metadata_KMU$Sample.ID )
  metadata_KMU <- metadata_KMU[order(metadata_KMU$HCC.CLD),]
  KMU_data <- KMU_data[,match(metadata_KMU$Sample.ID,colnames(KMU_data))]
  
  MGH_data <- read.csv('/Dropbox (Partners HealthCare)/CTC-Seq/input/HCC CLD RawCounts copy.csv', row.names = 1)
  
  #remove Biopsy, HepG2, HD, HL from mgh data
  MGH_data <- MGH_data[,!grepl('Biopsy|HepG2|HD|HL', colnames(MGH_data))]

  #remove samples with count < 100,000
  #don't remove low q Flow sorted data because looses too many samples
  KMU_data <- KMU_data[,colSums(KMU_data)>100000]
  metadata_KMU <- metadata_KMU[match(colnames(KMU_data), metadata_KMU$Sample.ID),]
  MGH_data <- MGH_data[,colSums(MGH_data)>100000]

  #RPM normalize, log10 scale
  RPM_KMU <- apply(KMU_data, 2, function(x) x*1e6/sum(x))
  RPM_KMU <- log10(RPM_KMU+1)
  RPM_MGH <- apply(MGH_data, 2, function(x) x*1e6/sum(x))
  RPM_MGH <- log10(RPM_MGH+1)
  
  #make KMU selected for best samples as well
  liverPanel <- c('FGB', 'FGG', 'ALB', 'AFP', 'FABP1', 'GPC3', 'RBP4', 'RBP4', 'AHSG', 'APOH', 'TF')
  liverScore <- colSums(RPM_KMU[rownames(RPM_KMU) %in% liverPanel,])
  
  #remove all samples with a liver score < 0.5 (less than 3 RPM for the panel)
  RPM_KMU_Best <- RPM_KMU[,liverScore>0.5]
  metadata_KMU_best <- metadata_KMU[liverScore>0.5,]
  
  #If desired, center expression. Else just use gene panel to subset data
  RPM_KMU_Best <- t(scale(t(RPM_KMU_Best[match(gene_panel, rownames(RPM_KMU_Best)),]), center = centered, scale = F))
  RPM_MGH <- t(scale(t(RPM_MGH[match(gene_panel, rownames(RPM_MGH)),]), center = centered, scale = F))
  

  #make annotations
  annotation_MGH <- data.frame(row.names=colnames(RPM_MGH),
                           diagnosis=ifelse(grepl("CLD",colnames(RPM_MGH)),"CLD","HCC"))
  annotation_KMU_Best <- data.frame(row.names=metadata_KMU_best$Sample.ID,
                               diagnosis=ifelse(grepl("CLD",metadata_KMU_best$HCC.CLD),"CLD","HCC"))
  KMUBest_HCC_metadata <- metadata_KMU_best[grepl('HCC', metadata_KMU_best$HCC.CLD),]
  KMUBest_CLD_metadata <- metadata_KMU_best[grepl('CLD', metadata_KMU_best$HCC.CLD),]

  #Make sure that all the data doesn't have NAs
  RPM_MGH <- RPM_MGH[!is.na(RPM_MGH[,1]),]
  RPM_KMU_Best <- RPM_KMU_Best[!is.na(RPM_KMU_Best[,1]),]
  
  #plot the heatmaps
  if(horizontalHeatmaps){
    pheatmap(t(RPM_MGH), main = 'MGH Data', 
             annotation_row = annotation_MGH, 
             show_rownames = FALSE,
             cluster_cols = cluster_genes,
             cluster_rows = FALSE,
             gaps_row = c(length(grep("CLD", annotation_MGH$diagnosis))))
    pheatmap(t(RPM_KMU_Best), main = 'KMU Best Samples',
             annotation_row = annotation_KMU_Best,
             show_rownames = FALSE,
             cluster_rows = FALSE,
             cluster_cols = cluster_genes,
             gaps_row = c(length(grep('CLD', annotation_KMU_Best))))
  }else{
    pheatmap(RPM_MGH, main = 'MGH Data', 
             annotation_col = annotation_MGH, 
             show_colnames = FALSE,
             cluster_rows = cluster_genes,
             cluster_cols = FALSE,
             gaps_col = c(length(grep("CLD", annotation_MGH$diagnosis))))
    pheatmap(RPM_KMU_Best, main = 'KMU Best Samples',
             annotation_col = annotation_KMU_Best,
             show_colnames = FALSE,
             cluster_cols = FALSE,
             cluster_rows = cluster_genes,
             gaps_col = c(length(grep('CLD', annotation_KMU_Best))))
  }
  #for each gene, run T Test and compare HCC and CLD
  genes_to_plot <- c()
  for (i in 1:nrow(RPM_MGH)){
    gene_counts_MGH <- data.frame(t(unlist(RPM_MGH[i,])))
    gene_counts_KMU <- data.frame(t(unlist(RPM_KMU_Best[i,])))
    
    gene_MGH_HCC<- gene_counts_MGH[,grepl('HCC', colnames(gene_counts_MGH))]
    gene_MGH_CLD <- gene_counts_MGH[,grepl('CLD', colnames(gene_counts_MGH))]
    
    gene_KMU_HCC <-gene_counts_KMU[colnames(gene_counts_KMU) %in% KMUBest_HCC_metadata$Sample.ID]
    gene_KMU_CLD <-gene_counts_KMU[colnames(gene_counts_KMU) %in% KMUBest_CLD_metadata$Sample.ID]
    
    MGH_test <- t.test(gene_MGH_HCC, gene_MGH_CLD)
    KMU_test <- t.test(gene_KMU_HCC, gene_KMU_CLD)
    
    #since these are centered, constant expression genes will give a var of 0, producing a nan p.value
    if(!is.nan(MGH_test$p.value) & !is.nan(KMU_test$p.value)){
      if (MGH_test$p.value<p.threshold & KMU_test$p.value<p.threshold){
        genes_to_plot <-c(genes_to_plot,rownames(RPM_MGH)[i])
      }
    }
  }
  
  #make boxplots of genes w/ p vals under 0.1
  MGH_HCC_all <- RPM_MGH[,grepl("HCC",colnames(RPM_MGH))]
  MGH_CLD_all <- RPM_MGH[,grepl("CLD",colnames(RPM_MGH))]
  KMUBest_CLD_all <-RPM_KMU_Best[,colnames(RPM_KMU_Best) %in% KMUBest_CLD_metadata$Sample.ID]
  KMUBest_HCC_all <-RPM_KMU_Best[,colnames(RPM_KMU_Best) %in% KMUBest_HCC_metadata$Sample.ID]
  
  for (gene in genes_to_plot){
    MGH_HCC<- MGH_HCC_all[gene,]
    MGH_CLD <- MGH_CLD_all[gene,]
    KMUBest_HCC <- KMUBest_HCC_all[gene,]
    KMUBest_CLD <- KMUBest_CLD_all[gene,]
    
    x <- c('MGH_HCC', 'MGH_CLD', 'KMUBest_HCC', 'KMUBest_CLD')
    dataList <- lapply(x, get, envir=environment())
    names(dataList) <- x
    
    par(mar = c(6, 5, 3, 2))
    boxplot(dataList, main=gene, ylab='RPM',las=2)
    stripchart(dataList, 
               col=ifelse(grepl('HCC',names(dataList)), 'blue',
                          'red'),
               vertical =TRUE, method="jitter", add=TRUE)
  }
  return(genes_to_plot)
}
