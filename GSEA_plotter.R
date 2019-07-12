GSEA_plotter <- function(results, gmt){
  #if there are more than 10 genes go 2 town
  require(org.Hs.eg.db)
  require(clusterProfiler)
  require(enrichplot)
  
  noGenes = F
  if (nrow(results) > 10){
    
    #make the gene names a decreasing ordered l2fc named vector
    results <- results[order(results$log2FoldChange, decreasing = T),]
    genes <- results$log2FoldChange
    names(genes) <- rownames(results)
    
    #make the GSEA object
    GSEA_output <- GSEA(geneList = genes, TERM2GENE = gmt, pvalueCutoff = 0.1)
    
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
}
