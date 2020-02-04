#################################
##FUNCTION TO CALCULATE PCA ERROR
#################################
PCAErrorMatrix <- function(rawData, PCs, selectedPCs='All'){
  
  #Set the defaults for ncol/selectedRows
  if('All'%in%selectedPCs){
    selectedPCs <- 1:ncol(PCs)
  }
  if(sum(rownames(PCs)!=rownames(rawData))>0){
    stop('PC Features do not align with raw data rownmaes')
  }
  
  #Trim data down to optioned length
  PCs <- PCs[,selectedPCs]
  
  #Intialize the Error data frame and the geneMeans for calculating euclidian distance
  errorDF <- matrix(nrow = ncol(rawData), ncol = ncol(PCs), dimnames = list(colnames(rawData), paste0('PC', seq(1, ncol(PCs)))))
  geneMeans <- rowMeans(rawData)
  
  #Iterate through each PC and find the error
  for(j in 1:ncol(PCs)){
    for(i in 1:ncol(rawData)){
      #Get the raw data to calculate the true euclidian distance from the origin
      #formula is sqrt((Val1-Mean1)^2+(Val2-Mean2)^2...(ValN-MeanN)^2)
      tempData <- rawData[,i]
      rawDist <- sqrt(sum((tempData-geneMeans)^2))
      
      #Multiply the Gene expression by the PC Rotation to get the transformed value
      #pcDist <- tempData*PCs[,j]
      
      #center the gene expression data, then multiply by PC rotation to get transformed value
      normDist <- tempData-geneMeans
      pcDist <- sqrt(sum((normDist*PCs[,j])^2))
      
      #Error is the difference between the optimized distance and PC Distance
      #errorDF[i,j] <- (rawDist-pcDist)^2
      #Output unsquared error to see if they are all positive (should be)
      errorDF[i,j] <- rawDist-pcDist
    }
  }
  return(errorDF)
}