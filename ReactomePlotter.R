#plots any significant pathways enriched in the reactome database
#require('org.Hs.eg.db')
#require('ReactomePA')

reactomePlotter <- function(results, squeeze, timepoint){
  #clean up any NAs
  genes <- rownames(results)
  genes_entrez <- bitr(genes, fromType = 'SYMBOL', toType = 'ENTREZID', OrgDb = 'org.Hs.eg.db')
  genes_entrez <- genes_entrez$ENTREZID
  genes_entrez <- genes_entrez[!is.na(genes_entrez)]
  genes_enrich <- enrichPathway(genes_entrez, organism = 'human', pvalueCutoff = 0.05, readable = T)
  
  if(length(genes_enrich)!=0){
    print(dotplot(genes_enrich, title = paste('Enriched Pathays in',squeeze, 'at hour', timepoint)))
  }else{
    par(mar = c(0,0,0,0))
    plot(c(0,1), c(0,1), ann = F, bty = 'n', type = 'n', xaxt='n', yaxt='n')
    text(x = 0.5, y = 0.5, paste('There were no significant genes for ReactomePA.'))
    par(mar = c(5, 4, 4, 2) + 0.1)
  }  
}