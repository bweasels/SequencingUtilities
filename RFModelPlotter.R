RFModelPlotter <- function(RPM.df, metadata, model, modelName){
  #remove the leafs <- we only want the nodes which lead somewhere to make the dendrograph
  nodes <- modelInfo[!is.na(modelInfo$leftChild),]
  
  #make the model frame to give to graph
  model_frame <- data.frame(from = rep(nodes$nodeID, 2), to = c(nodes$leftChild,nodes$rightChild))
  model_frame <- model_frame[!is.na(model_frame$to),]
  
  #tack on filler to the nodes dataframe so that there are NAs in appropriate areas
  filler <- matrix(nrow = nrow(modelInfo)-nrow(nodes), ncol = ncol(modelInfo))
  colnames(filler) <- colnames(modelInfo)
  nodes <- rbind(nodes, filler)
  
  #make the graph
  graph <- graph_from_data_frame(model_frame)
  
  #match predictions to the order that graph gave to the nodes/leaves
  predictions <- as.character(modelInfo$prediction[as.numeric(V(graph)$name)+1])
  
  #feed graph the various information
  V(graph)$node_label <- nodes$splitvarName
  V(graph)$leaf_label <- predictions
  V(graph)$split <- as.character(round(nodes$splitval, digits = 3))
  
  #plot the dendrograph
  pdf(paste0(modelName,'_RFDecisionPlots.pdf'))
  plot <- ggraph(graph, 'dendrogram') + 
    theme_bw() +
    geom_edge_link() +
    geom_node_point() +
    geom_node_text(aes(label = node_label), na.rm = TRUE, repel = TRUE) +
    geom_node_label(aes(label = split), vjust = 2.5, na.rm = TRUE, fill = "white") +
    geom_node_label(aes(label = leaf_label, fill = leaf_label), na.rm = TRUE, 
                    repel = TRUE, colour = "white", fontface = "bold", show.legend = FALSE) +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(),
          panel.background = element_blank(),
          plot.background = element_rect(fill = "white"),
          panel.border = element_blank(),
          axis.line = element_blank(),
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          plot.title = element_text(size = 18)) + 
    ggtitle(paste0('Random Forest Decision Tree: ', modelName, '\nRight branch is less than, left is greater'))
  
  print(plot)
  
  #get the requsite info for the box plots
  genes <- modelInfo$splitvarName[!is.na(modelInfo$splitvarName)]
  splitVal <- modelInfo$splitval[!is.na(modelInfo$splitvarName)]
  nodes <- modelInfo$nodeID[!is.na(modelInfo$splitvarName)]
  
  #cut RPM.df down to the important data & append the diagnosis to it
  RPM.df <- RPM.df[rownames(RPM.df)%in%genes,]
  data <- data.frame(t(RPM.df), diagnosis = metadata$diagnosis)
  
  #plot box plots
  for (i in 1:length(genes)){
    g <- ggplot(data, aes_string(x='diagnosis', y = genes[i])) + geom_boxplot(outlier.shape = NA)
    g <- g + geom_jitter(position = position_jitter(width = 0.2, height = 0)) + geom_hline(yintercept = splitVal[i], color = 'red')
    g <- g + ggtitle(paste('Decision at Node:', nodes[i], '\nSplit Value at', splitVal[i]))
    print(g)
  }
  dev.off()
}