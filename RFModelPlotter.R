RFModelPlotter <- function(RPM.df, metadata, modelInfo, modelName){
  require(igraph)
  require(ggraph)

  date <- Sys.Date()
  
  #find out how many sample IDs are left for each node
  #since the first node is a special case, initialize it outside of the loop
  samples <- vector('list', length = nrow(modelInfo))
  samples[[1]] <- colnames(RPM.df)
  names(samples)[1] <- modelInfo$splitvarName[1]
  
  for(i in 1:nrow(modelInfo)){
    #get the gene used to decide positive/negative 
    gene <- modelInfo$splitvarName[i]
    if(!is.na(gene)){
      #get the necessary info to count the number of samples on each side of the split
      names(samples)[i] <- gene
      lNode <- modelInfo$leftChild[i]
      rNode <- modelInfo$rightChild[i]
      sep <- modelInfo$splitval[i]
      
      #Select for the splitting gene, and the samples left in this node
      tempDF <- RPM.df[gene==rownames(RPM.df),match(samples[[i]],colnames(RPM.df))]
      samples[[lNode+1]] <- names(tempDF[tempDF<sep])
      samples[[rNode+1]] <- names(tempDF[tempDF>sep])
    }
  }
  
  #add the nSamples to the modelInfo for the dendriograph
  modelInfo$nSamples <- lengths(samples)
  
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
  
  #match predictions & nSamples to the order that graph gave to the nodes/leaves
  predictions <- as.character(modelInfo$prediction[as.numeric(V(graph)$name)+1])
  nSamples <- as.character(modelInfo$nSamples[as.numeric(V(graph)$name)+1])
  
  #feed graph the various info
  V(graph)$node_label <- nodes$splitvarName
  V(graph)$leaf_label <- predictions
  V(graph)$nSamplesRemain <- nSamples 
  V(graph)$split <- as.character(round(nodes$splitval, digits = 3))
  
  #plot the dendrograph
  dimSize <- 7 + 0.1*nrow(modelInfo)
  pdf(paste0(date, '_', modelName,'_RFDecisionPlots.pdf'), width = dimSize+1, height = dimSize )
  plot <- ggraph(graph, 'dendrogram') + 
    theme_bw() +
    geom_edge_link() +
    geom_node_point() +
    geom_node_label(aes(filter = !is.na(node_label), label = paste0('Gene: ', node_label, '\nSplit Val: ', split, '\nSamples: ', nSamplesRemain))) +
    geom_node_label(aes(filter = is.na(node_label), label = paste0('Decision: ', leaf_label, '\nSamples: ', nSamplesRemain), fill = leaf_label), fontface = 'bold', show.legend = F) +
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
  #also cut the samples down to only the forking nodes (no leaves)
  RPM.df <- RPM.df[rownames(RPM.df)%in%genes,]
  data <- data.frame(t(RPM.df), diagnosis = metadata$diagnosis)
  samples <- samples[match(genes, names(samples))]
  
  #plot box plots
  for (i in 1:length(genes)){
    #select only the relevant samples
    subData <- data[match(samples[[i]],rownames(data)),]
    g <- ggplot(subData, aes_string(x='diagnosis', y = genes[i])) + geom_boxplot(outlier.shape = NA)
    g <- g + geom_jitter(position = position_jitter(width = 0.2, height = 0)) + geom_hline(yintercept = splitVal[i], color = 'red')
    g <- g + ggtitle(paste('Decision at Node:', nodes[i], '\nSplit Value at', splitVal[i]))
    print(g)
  }
  dev.off()
}