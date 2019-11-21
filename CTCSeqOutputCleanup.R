CTCSeqOutputCleanup <- function(falseAUC, trueAUC){
  
  #Reset rownames to be compatable with quickly finding the folder
  rownames(falseAUC) <- gsub('10_(.*)', '\\1', rownames(falseAUC))
  rownames(trueAUC) <- gsub('10_(.*)', '\\1', rownames(trueAUC))
  
  #truncate falseAUCs to just AUCs
  falseAUC <- falseAUC[,grep('AUC', colnames(falseAUC))]
  trueAUC <- trueAUC[,grep('AUC', colnames(trueAUC))]
  
  #reset colnames for quick foldername paste
  colnames(falseAUC) <- c('10_', '20_', '50_', '100_', '200_', '500_', '1000_', '2000_', '5000_', '10000_')
  colnames(trueAUC) <- c('10_', '20_', '50_', '100_', '200_', '500_', '1000_', '2000_', '5000_', '10000_')
  
  #return it
  return(list(falseAUC, trueAUC))
}

extractBestModels <- function(falseAUC, trueAUC, randFalseAUC, randTrueAUC, threshold = 0.2){
  
  #compress the CTCSeq output into a neat dataframe
  out <- CTCSeqOutputCleanup(falseAUC, trueAUC)
  falseAUC <- out[[1]]
  trueAUC <- out[[2]]
  
  randOut <- CTCSeqOutputCleanup(randFalseAUC, randTrueAUC)
  randFalseAUC <- randOut[[1]]
  randTrueAUC <- randOut[[2]]
  
  #get the average reference AUC for each nGenes
  randFalseAUC <- colMeans(randFalseAUC)
  randTrueAUC <- colMeans(randTrueAUC)
  
  bestFalseModels <- c()
  bestTrueModels <- c()
  
  #Go through each column & get the threshold clearing models
  for(i in 1:ncol(falseAUC)){
    
    #take the models that supercede their threshold over random
    samps <- falseAUC[,i,drop=F]
    samps <- samps[samps>randFalseAUC[i]+threshold,,drop = F]
    bestFalse <- unlist(samps)
    if(length(bestFalse)>0){
      names(bestFalse) <- paste0(colnames(samps), rownames(samps))
      bestFalseModels <- c(bestFalseModels, bestFalse)
    }
    
    #Redo for True AUCs
    samps <- trueAUC[,i,drop=F]
    samps <- samps[samps>randTrueAUC[i]+threshold,,drop=F]
    bestTrue <- unlist(samps)
    if(length(bestTrue)>0){
      names(bestTrue) <- paste0(colnames(samps), rownames(samps))
      bestTrueModels <- c(bestTrueModels, bestTrue)
    }
  }
  
  return(list(bestFalseModels, bestTrueModels))
}