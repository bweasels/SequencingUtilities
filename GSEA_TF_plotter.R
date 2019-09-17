#Plots the transcription factor binding sites and tells if they are significant or not
#requires 
#library(TFEA.ChIP)

GSEA_TF_Plotter <- function(results, outputDir, timepoint){
  #get the data and make the filename for outputs
  fileName <- results@elementMetadata
  fileName <- grep('Condition',fileName$description, value = T)
  fileName <- gsub('.*(Condition.*)', '\\1', fileName[1])
  fileName <- gsub(' ', '.', fileName)
  results <- data.frame(results@listData, Genes = rownames(results))
  
  if(nrow(results) > 5){
    #prep data for the GSEA run & generate the plot of TFs
    data <- preprocessInputData(results)
    GSEA.results <- GSEA_run(data$Genes, data$log2FoldChange, get.RES = T)
    p <- plot_ES(GSEA.results, data$log2FoldChange)
    htmlwidgets::saveWidget(p, file = paste0(outputDir, timepoint,'hr.', fileName, '.TFPlots.html'))
    
    #Find the significant TFs (pval < 0.05) 
    enrichTable <- GSEA.results[['Enrichment.table']]
    enrichTable <- enrichTable[enrichTable$pval.ES<.05,]
    
    #create the list of unique elements & get their counts
    Tfactors <- unique(GSEA.results[['Enrichment.table']]$TF)
    counts <- vector(length = length(Tfactors))
    names(counts) <- Tfactors
    
    for (i in 1:length(Tfactors)){
      counts[i] <- sum(GSEA.results[['Enrichment.table']]$TF %in% Tfactors[i])
    }
    
    counts <- sort(counts, decreasing = T)
    counts <- counts[counts>5]
    
    counts <- data.frame(TF = factor(names(counts), levels = names(counts)), count = counts)
    
    pdf(paste0(outputDir, timepoint,'hr.', fileName, 'TFCounts.pdf'), width = 20)
    gplot <- ggplot(counts, aes(x = TF, y = count)) + geom_bar(stat='identity')
    print(gplot)
    dev.off()
  }
}