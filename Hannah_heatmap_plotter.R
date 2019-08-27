heatmap_plot <- function(gene_panel, p.threshold = 0.1, cluster_genes = TRUE, centered = TRUE, scaled = FALSE){
  require(pheatmap)
  require(ggplot2)
  
  #read in files, convert to matrices
  KMU_data <- read.csv('/MGH/LiverProject/CTC-Seq/input/2019-08-12_KMUCounts_star_genes_erc.csv', row.names = 1)
  metadata_KMU <- read.csv('/MGH/LiverProject/CTC-Seq/input/2019-08-12_KMUMetadata.csv', row.names = 1)
  metadata_KMU$Sample.ID <- sub("^", "X", metadata_KMU$Sample.ID )
  metadata_KMU <- metadata_KMU[order(metadata_KMU$HCC.CLD),]
  KMU_data <- KMU_data[,match(metadata_KMU$Sample.ID,colnames(KMU_data))]
  
  MGH_data <- read.csv('/MGH/LiverProject/CTC-Seq/input/HCC CLD RawCounts copy.csv', row.names = 1)
  
  PBMC_data <- read.csv('/MGH/06-25-19/IFD RNA Seq/PBMCExp_countsTable.csv',row.names = 1)
  
  FlowSort_data <- read.csv('/MGH/LiverProject/CTC-Seq/feature/Cross Sample Feature Set/IFD_flowsort_RawCounts_tbl.csv', row.names = 1)

  #remove Biopsy, HepG2, HD, HL from mgh data
  MGH_data <- MGH_data[,!grepl('Biopsy|HepG2|HD|HL', colnames(MGH_data))]

  #remove samples with count < 100,000
  #don't remove low q Flow sorted data because looses too many samples
  KMU_data <- KMU_data[,colSums(KMU_data)>100000]
  metadata_KMU <- metadata_KMU[match(colnames(KMU_data), metadata_KMU$Sample.ID),]
  MGH_data <- MGH_data[,colSums(MGH_data)>100000]
  PBMC_data <- PBMC_data[,colSums(PBMC_data)>100000]

  #RPM normalize, log10 scale
  RPM_KMU <- apply(KMU_data, 2, function(x) x*1e6/sum(x))
  RPM_KMU <- log10(RPM_KMU+1)
  RPM_MGH <- apply(MGH_data, 2, function(x) x*1e6/sum(x))
  RPM_MGH <- log10(RPM_MGH+1)
  RPM_PBMC <- apply(PBMC_data, 2, function(x) x*1e6/sum(x))
  RPM_PBMC <- log10(RPM_PBMC+1)
  RPM_FlowSort <- apply(FlowSort_data, 2, function(x) x*1e6/sum(x))
  RPM_FlowSort <- log10(RPM_FlowSort+1)
  
  #If desired, center expression. Else just use gene panel to subset data
  RPM_KMU <- t(scale(t(RPM_KMU[match(gene_panel, rownames(RPM_KMU)),]), center = centered, scale = F))
  RPM_MGH <- t(scale(t(RPM_MGH[match(gene_panel, rownames(RPM_MGH)),]), center = centered, scale = F))
  RPM_PBMC <- t(scale(t(RPM_PBMC[match(gene_panel, rownames(RPM_PBMC)),]), center = centered, scale = F))
  RPM_FlowSort <- t(scale(t(RPM_FlowSort[match(gene_panel, rownames(RPM_FlowSort)),]), center = centered, scale = F))
  
  #breakout Flowsorted Data into HCC and CLD samples and similarly average the cell within them
  RPM_FlowSort_HCC <- RPM_FlowSort[,grep('HCC', colnames(RPM_FlowSort))]
  RPM_FlowSort_CLD <- RPM_FlowSort[,grep('CLD', colnames(RPM_FlowSort))]
  cellType <- c('B', 'G', 'M', 'N', 'C', 'H')
  RPM_FlowSort_Sub <- matrix(nrow = nrow(RPM_FlowSort), ncol = length(cellType), dimnames = list(rownames(RPM_FlowSort), c('B.Cell', 'Granulocytes', 'Monocytes', 'NK.Cells', 'Cytotoxic.T.Cells', 'Helper.T.Cells')))
  cellType <- c('B', 'G', 'M', 'N', 'C', 'H')
  for (i in 1:length(cellType)) {
    HCC_temp <- rowMeans(RPM_FlowSort_HCC[,grep(cellType[i], colnames(RPM_FlowSort_HCC))])
    CLD_temp <- rowMeans(RPM_FlowSort_CLD[,grep(cellType[i], colnames(RPM_FlowSort_CLD))])
    RPM_FlowSort_Sub[,i] <- HCC_temp - CLD_temp
  }
  
  #make annotations
  annotation_MGH <- data.frame(row.names=colnames(RPM_MGH),
                           diagnosis=ifelse(grepl("CLD",colnames(RPM_MGH)),"CLD","HCC"))
  annotation_PBMC <- data.frame(row.names=colnames(RPM_PBMC),
                           diagnosis=ifelse(grepl("CLD",colnames(RPM_PBMC)),"CLD","HCC"))
  annotation_KMU <- data.frame(row.names=metadata_KMU$Sample.ID,
                               diagnosis=ifelse(grepl("CLD",metadata_KMU$HCC.CLD),"CLD","HCC"))
  KMU_HCC_metadata <- metadata_KMU[grepl('HCC',metadata_KMU$HCC.CLD),]
  KMU_CLD_metadata <- metadata_KMU[grepl('CLD',metadata_KMU$HCC.CLD),]

  #Make sure that all the data doesn't have NAs
  RPM_MGH <- RPM_MGH[!is.na(RPM_MGH[,1]),]
  RPM_PBMC <- RPM_PBMC[!is.na(RPM_PBMC[,1]),]
  RPM_KMU <- RPM_KMU[!is.na(RPM_KMU[,1]),]
  RPM_FlowSort_Sub <- RPM_FlowSort_Sub[!is.na(RPM_FlowSort_Sub[,1]),]
  
  #plot the heatmaps
  pheatmap(RPM_MGH, main = 'MGH Data', 
           annotation_col = annotation_MGH, 
           show_colnames = FALSE,
           cluster_cols = FALSE,
           cluster_rows = cluster_genes,
           gaps_col = c(length(grep("CLD", annotation_MGH$diagnosis))))
  pheatmap(RPM_PBMC, main='PBMC Data', 
           annotation_col= annotation_PBMC, 
           show_colnames=FALSE,
           cluster_cols=FALSE,
           cluster_rows = cluster_genes,
           gaps_col = c(length(grep("CLD", annotation_PBMC$diagnosis))))
  pheatmap(RPM_KMU, main='KMU Data', 
           annotation_col= annotation_KMU, 
           show_colnames=FALSE,
           cluster_cols=FALSE,
           cluster_rows = cluster_genes,
           gaps_col = c(length(grep("CLD", annotation_KMU$diagnosis)))-1)
  pheatmap(RPM_FlowSort_Sub, main = 'Flow Sorted Data (HCC-CLD)',
           show_colnames = T,
           cluster_cols = F,
           cluster_rows = cluster_genes)

  #for each gene, run T Test and compare HCC and CLD
  genes_to_plot <- c()
  for (i in 1:nrow(RPM_MGH)){
    gene_counts_MGH <- data.frame(t(unlist(RPM_MGH[i,])))
    gene_counts_KMU <- data.frame(t(unlist(RPM_KMU[i,])))
    
    gene_MGH_HCC<- gene_counts_MGH[,grepl('HCC', colnames(gene_counts_MGH))]
    gene_MGH_CLD <- gene_counts_MGH[,grepl('CLD', colnames(gene_counts_MGH))]
    
    gene_KMU_HCC <-gene_counts_KMU[colnames(gene_counts_KMU) %in% KMU_HCC_metadata$Sample.ID]
    gene_KMU_CLD <-gene_counts_KMU[colnames(gene_counts_KMU) %in% KMU_CLD_metadata$Sample.ID]
    
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
  PBMC_HCC_all <- RPM_PBMC[,grepl("CTC",colnames(RPM_PBMC))]
  PBMC_CLD_all <- RPM_PBMC[,grepl("CLD",colnames(RPM_PBMC))]
  KMU_HCC_all <-RPM_KMU[,colnames(RPM_KMU) %in% KMU_HCC_metadata$Sample.ID]
  KMU_CLD_all <-RPM_KMU[,colnames(RPM_KMU) %in% KMU_CLD_metadata$Sample.ID]
  
  for (gene in genes_to_plot){
    MGH_HCC<- MGH_HCC_all[gene,]
    MGH_CLD <- MGH_CLD_all[gene,]
    PBMC_HCC <- PBMC_HCC_all[gene,]
    PBMC_CLD <- PBMC_CLD_all[gene,]
    KMU_HCC <- KMU_HCC_all[gene,]
    KMU_CLD <- KMU_CLD_all[gene,]
    
    x <- c('MGH_HCC', 'MGH_CLD', 'PBMC_HCC', 'PBMC_CLD', 'KMU_HCC', 'KMU_CLD')
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
