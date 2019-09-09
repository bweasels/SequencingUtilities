###########################
#File: GSEA_plotter.R
#Purpose: Plot out pathway enrichments using the clusterProfiler package in R
#Inputs: results: The Results dataframe from DESeq2::results(deseqObject, tidy = T) ***Before handing the results, filter out all genes with a p.adj == NA
#        gmt: .GMT File with the pathway weightings, which can be found at http://software.broadinstitute.org/gsea/msigdb/index.jsp
#        fdr: False discovery rate - default is 0.1, but can set lower or higher if needed
#
#Outputs: None, but does generate plots in the function, so have to capture them outside the function
#Required Packages: clusterProfiler, ggplot2, DESeq2
###########################

GSEA_plotter <- function(results, gmt, fdr = 0.1){
  #if there are more than 10 genes go 2 town
  require(org.Hs.eg.db)
  require(clusterProfiler)
  require(enrichplot)
  
  output <- NULL
  noGenes = F
  if (nrow(results) > 10){
    
    #make the gene names a decreasing ordered l2fc named vector
    results <- results[order(results$log2FoldChange, decreasing = T),]
    genes <- results$log2FoldChange
    names(genes) <- results$row
    
    #make the GSEA object
    GSEA_output <- GSEA(geneList = genes, TERM2GENE = gmt, pvalueCutoff = fdr)
    output <- GSEA_output@result
    output <- cbind(ID = output$ID, Enrichment.Score = output$enrichmentScore, Norm.Enrichment.Score = output$NES, pvalue = output$pvalue, p.adjust = output$p.adjust, qvalues = output$qvalues, max.rank = output$rank)
        
    #Plot!
    if(nrow(GSEA_output) > 0){
      print(heatplot(GSEA_output, foldChange = genes))
      print(ridgeplot(GSEA_output))
      for(i in 1:nrow(GSEA_output@result)){
        print(gseaplot(GSEA_output, geneSetID=i, title = GSEA_output@result$Description[i]))
      }
      #If there are no enriched GSEA pathways or no singificant results, print the no genes statement
    }else{noGenes = T}
  }else{noGenes = T}
  
  #Code to print the noGenes Statment
  if (noGenes){
    par(mar = c(0,0,0,0))
    plot(c(0,1), c(0,1), ann = F, bty = 'n', type = 'n', xaxt='n', yaxt='n')
    text(x = 0.5, y = 0.5, paste('There were no significant genes for', unique(gsub('^([A-Z]*)_.*', '\\1', gmt$ont)), '.'))
    par(mar = c(5, 4, 4, 2) + 0.1)
  }
  return(output)
}
