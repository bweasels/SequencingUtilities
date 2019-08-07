heatmap_plot <- function(gene_panel, p.threshold = 0.1, cluster_genes = TRUE){

  #read in files, convert to matrices
  KMU_data <- read.csv('/MGH/LiverProject/CTC-Seq/input/2019-07-23_KMUCounts_star_genes_erc.counts.csv')
  KMU_data_matrix <- as.matrix(KMU_data[,-1])
  rownames(KMU_data_matrix) <-KMU_data[,1]
  metadata_KMU <- read.csv('/MGH/LiverProject/CTC-Seq/input/2019-07-23_KMUMetadata.csv')
  metadata_KMU$Sample.ID <- sub("^", "X", metadata_KMU$Sample.ID )
  metadata_KMU<-metadata_KMU[order(metadata_KMU$HCC.CLD),]
  KMU_data_matrix<-KMU_data_matrix[,match(metadata_KMU$Sample.ID,colnames(KMU_data_matrix))]
  
  MGH_data <- read.csv('/MGH/LiverProject/CTC-Seq/input/HCC CLD RawCounts copy.csv')
  MGH_data_matrix <- as.matrix(MGH_data[,-1])
  rownames(MGH_data_matrix) <-MGH_data[,1]
  
  PBMC_data <- read.csv('/MGH/06-25-19/IFD RNA Seq/PBMCExp_countsTable.csv')
  PBMC_data_matrix <- as.matrix(PBMC_data[,-1])
  rownames(PBMC_data_matrix) <-PBMC_data[,1]

  #remove Biopsy, HepG2, HD, HL from mgh data
  MGH_data_matrix <- MGH_data_matrix[,!grepl('Biopsy|HepG2|HD|HL', colnames(MGH_data_matrix))]

  #remove samples with count < 100,000, then remove sum column
  KMU_sum <-addmargins(KMU_data_matrix)
  filter_KMU <- KMU_sum[,KMU_sum['Sum',]>100000]
  filter_KMU <- subset( filter_KMU, select = -Sum )
  MGH_sum <-addmargins(MGH_data_matrix)
  filter_MGH <- MGH_sum[,MGH_sum['Sum',]>100000]
  filter_MGH <- subset( filter_MGH, select = -Sum )
  PBMC_sum <-addmargins(PBMC_data_matrix)
  filter_PBMC <- PBMC_sum[,PBMC_sum['Sum',]>100000]
  filter_PBMC <- subset( filter_PBMC, select = -Sum )
  
  #RPM normalize, log10 scale
  RPM_KMU <- apply(filter_KMU, 2, function(x)x*1e6/sum(x))
  RPM_KMU <- log10(RPM_KMU+1)
  RPM_MGH <- apply(filter_MGH, 2, function(x)x*1e6/sum(x))
  RPM_MGH <- log10(RPM_MGH+1)
  RPM_PBMC <- apply(filter_PBMC, 2, function(x)x*1e6/sum(x))
  RPM_PBMC <- log10(RPM_PBMC+1)
  
  #Use gene panel to subset data, center expression
  RPM_KMU <- scale(RPM_KMU[match(gene_panel, rownames(RPM_KMU)),], center = T, scale = F)
  RPM_MGH <- scale(RPM_MGH[match(gene_panel, rownames(RPM_MGH)),], center = T, scale = F)
  RPM_PBMC <- scale(RPM_PBMC[match(gene_panel, rownames(RPM_PBMC)),], center = T, scale = F)

  #make annotations
  annotation_MGH <- data.frame(row.names=colnames(filter_MGH),
                           diagnosis=ifelse(grepl("CLD",colnames(filter_MGH)),"CLD","HCC"))
  annotation_PBMC <- data.frame(row.names=colnames(filter_PBMC),
                           diagnosis=ifelse(grepl("CLD",colnames(filter_PBMC)),"CLD","HCC"))
  annotation_KMU <- data.frame(row.names=metadata_KMU$Sample.ID,
                               diagnosis=ifelse(grepl("CLD",metadata_KMU$HCC.CLD),"CLD","HCC"))
  KMU_HCC_metadata <- metadata_KMU[grepl('HCC',metadata_KMU$HCC.CLD),]
  KMU_CLD_metadata <- metadata_KMU[grepl('CLD',metadata_KMU$HCC.CLD),]

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
    
    if (MGH_test$p.value<p.threshold && KMU_test$p.value<p.threshold){
      genes_to_plot <-c(genes_to_plot,rownames(RPM_MGH)[i])
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
}
