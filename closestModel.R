#Utility to pull an model which is closest to the overall average
closestModel <- function(auc, aucHolder, nGenes){
  
  #Clean up the output
  data <- aucHolder[,grep('AUC', colnames(aucHolder))]
  colnames(data) <- c(10,20,50,100,200,500,1000,2000,5000,10000)
  rownames(data) <- gsub('10(.*)', '\\1', rownames(data))
  
  #Get the set of runs that we're interested in
  runsofInterest <- data[,match(nGenes, colnames(data))]
  minIndex <- 1
  minDifference <- runsofInterest[1]-auc
  
  #Loop through all aucs to find the one that is closest
  for(i in 1:length(runsofInterest)){
    if(abs(runsofInterest[i]-auc) < minDifference){
      minIndex <- i
      minDifference <- runsofInterest[i]-auc
    }
  }
  return(paste0(nGenes, rownames(data)[minIndex]))
}