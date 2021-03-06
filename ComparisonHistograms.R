ComparisonHistogram <- function(dataset1, dataset2, labels, nBins = 50, title = 'Comparison Histogram', xTitle = NULL){
  require(ggplot2)
  data <- data.frame(dat = c(dataset1, dataset2), 
                     lab = c(rep(labels[1], length(dataset1)),
                                rep(labels[2], length(dataset2))))
  g <- ggplot(data, aes(x=dat))
  g <- g + geom_histogram(data=subset(data,lab == labels[1]), 
                          fill = 'red', 
                          alpha = 0.3,
                          bins = nBins)
  g <- g + geom_histogram(data=subset(data,lab == labels[2]), 
                          fill = 'blue', 
                          alpha = 0.3,
                          bins = nBins)
  g <- g + ggtitle(paste(title, '\n Red:', labels[1], '| Blue:', labels[2]))
  if(!is.null(xTitle)){
    g <- g + labs(x = xTitle)
  }
  plot(g)
}