#Structure
#Average the Actin value for each cell line
#Subtract the actin average from the primer CTs (Primer - Actin) <- this is Delta CT
####Break out into a matrix of control and condtions 
#Create a table of averaged per primer CTs from the control matrix
#Subtract the average primer CTs from the control matrix and the condition matrices (control or condition - average)
#2^ the table values
#Stick everything together into a matrix of SD and Means and primer and stick all together

QPCR_plotter <- function(filename, dox){
  
  data <- read.csv(file= filename, stringsAsFactors = F)
  data <- data.frame(Sample = data$Sample.Name, Target = data$Target.Name, CT = data$C.)
  
  targets <- unique(data$Target)
  samples <- unique(data$Sample)
  
  #gather the actin average
  actinAvg <- vector('numeric', length = length(samples))
  
  actinVals <- filter(data, Target == 'Actin' | Target == 'ACTIN')
  for (i in 1:length(samples)){
    temp <- filter(actinVals, Sample == samples[i])
    actinAvg[i] <- mean(temp$CT, na.rm = T)
  }
  
  #subtract the actin average from the other CTs
  data <- data[!(data$Target=='Actin'|data$Target=='ACTIN'),]
  for (i in 1:length(samples)){
    selection <- data$Sample==samples[i]
    data[selection,]$CT <- data[selection,]$CT - actinAvg[i]
  }
  
  #Break out a matrix of control and conditions
  #doxSamples <- data[grep('[^no] dox|with Dox', data$Sample),]
  #noDoxSamples <- data[grep('no Dox|no dox', data$Sample),]
  
  if(dox){
    #doxSamples <- data[grep('[^-]DOX', data$Sample),]
    #noDoxSamples <- data[grep('-DOX', data$Sample),]
    
    doxSamples <- data[grep('10', data$Sample),]
    noDoxSamples <- data[!grepl('10', data$Sample),]
    
    #Gather noDox averages
    noDoxSamps <- unique(noDoxSamples$Sample)
    doxSamps <- unique(doxSamples$Sample)
    targets <- unique(noDoxSamples$Target)
    doxAvgs <- matrix(ncol = length(noDoxSamps), nrow = length(targets), dimnames = list(targets, noDoxSamps))
    
    #Copy the dox averages to the matrix
    for (i in 1:length(targets)){
      temp <- filter(noDoxSamples, Target == targets[i])
      for (j in 1:length(noDoxSamps)){
        temp2 <- filter(temp, Sample == noDoxSamps[j])
        doxAvgs[i,j] <- mean(temp2$CT, na.rm = T)
      }
    }
    
    #make average and SD matrices
    averagesNoDox <- matrix(nrow = length(targets), ncol = length(noDoxSamps), dimnames = list(targets, noDoxSamps))
    averagesDox <- matrix(nrow = length(targets), ncol = length(doxSamps), dimnames = list(targets, doxSamps))
    
    SDNoDox <- matrix(nrow = length(targets), ncol = length(noDoxSamps), dimnames = list(targets, noDoxSamps))
    SDDox <- matrix(nrow = length(targets), ncol = length(doxSamps), dimnames = list(targets, doxSamps))
    
    samplesDoxOrder <- averagesDox
    samplesNoDoxOrder <- averagesNoDox
    targetsOrder <- averagesDox
    targetsNoDoxOrder <- averagesNoDox
    
    #Subtract the noDox avgs from noDox values
    for (i in 1:length(targets)){
      for (j in 1:length(noDoxSamps)) {
        noDoxCT <- noDoxSamples[noDoxSamples$Target==targets[i]&noDoxSamples$Sample==noDoxSamps[j],]$CT - doxAvgs[i,j]
        noDoxCT <- 2^noDoxCT
        averagesNoDox[i,j] <- mean(noDoxCT, na.rm = T)
        SDNoDox[i,j] <- sd(noDoxCT, na.rm = T)
        samplesNoDoxOrder[i,j] <- as.character(noDoxSamps[j])
        targetsNoDoxOrder[i,j] <- as.character(targets[i])
      }
    }
    
    for (i in 1:length(targets)){
      for (j in 1:length(doxSamps)){
        doxCT <- doxSamples[doxSamples$Target==targets[i]&doxSamples$Sample==doxSamps[j],]$CT - doxAvgs[i,floor(j/2+0.5)] ###TEMP FIX
        doxCT <- 2^doxCT
        averagesDox[i,j] <- mean(doxCT, na.rm = T)
        SDDox[i,j] <- sd(doxCT, na.rm = T)
        samplesDoxOrder[i,j] <- as.character(doxSamps[j])
        targetsOrder[i,j] <- as.character(targets[i])
      }
    }
    
    outputNoDox <- data.frame(Sample = as.character(as.vector(samplesNoDoxOrder)), Target = as.character(as.vector(targetsNoDoxOrder)), CT = as.vector(averagesNoDox), SD = as.vector(SDNoDox))
    outputDox <- data.frame(Sample = as.character(as.vector(samplesDoxOrder)), Target = as.character(as.vector(targetsOrder)), CT = as.vector(averagesDox), SD = as.vector(SDDox))
    
    output <- rbind(outputNoDox, outputDox)
    output$Sample <- factor(output$Sample, levels = unique(output$Sample[order(as.character(output$Sample))]))
  }else{ #for non Dox situations
    
    ###Change the grep depending on the normalization condition
    control <- data[grep('Untreat|untreat', data$Sample),]
    data <- data[!grepl('Untreat|untreat', data$Sample),]
    
    #get the list of targets
    targets <- unique(data$Target)
    controlSamples <- unique(control$Sample)
    dataSamples <- unique(data$Sample)
    
    controlAvgs <- vector('numeric', length = length(targets))
    
    #get the average for each target within the control group
    for (i in 1:length(controlAvgs)){
      temp <- control[control$Target == targets[i],]
      controlAvgs[i] <- mean(temp$CT, na.rm = T)
    }
    
    #make the matrices to hold the control data
    ctrlAverage <- matrix(nrow = length(targets), ncol = length(controlSamples), dimnames = list(targets, controlSamples))
    ctrlSd <- matrix(nrow = length(targets), ncol = length(controlSamples), dimnames = list(targets, controlSamples))
    ctrlSamplesOrder <- matrix(nrow = length(targets), ncol = length(controlSamples), dimnames = list(targets, controlSamples))
    ctrlTargetsOrder <- matrix(nrow = length(targets), ncol = length(controlSamples), dimnames = list(targets, controlSamples))
    
    #make the matrices to hold the other data
    dataAverage <- matrix(nrow = length(targets), ncol = length(dataSamples), dimnames = list(targets, dataSamples))
    dataSd <- matrix(nrow = length(targets), ncol = length(dataSamples), dimnames = list(targets, dataSamples))
    dataSamplesOrder <- matrix(nrow = length(targets), ncol = length(dataSamples), dimnames = list(targets, dataSamples))
    dataTargetsOrder <- matrix(nrow = length(targets), ncol = length(dataSamples), dimnames = list(targets, dataSamples))
    
    #assort the relevant control data
    for (i in 1:length(targets)){
      temp <- control[control$Target==targets[i],]
      for (j in 1:length(controlSamples)){
        #get the relevant CT values and subtract the control values from it
        temp_samp <- temp[temp$Sample==controlSamples[j],]
        temp_ct <- temp_samp$CT - controlAvgs[i]
        
        #raise everything to the second power
        temp_ct <- 2^temp_ct
        
        #put the data in their relevant dataobjects
        ctrlAverage[i,j] <- mean(temp_ct, na.rm = T)
        ctrlSd[i,j] <- sd(temp_samp$CT, na.rm = T)
        ctrlSamplesOrder[i,j] <- as.character(controlSamples[j])
        ctrlTargetsOrder[i,j] <- as.character(targets[i])
      }
    }
    
    #assort the relevant other data data
    for (i in 1:length(targets)){
      temp <- data[data$Target==targets[i],]
      for (j in 1:length(dataSamples)){
        #get the relevant CT values and subtract the control values from it
        temp_samp <- temp[temp$Sample==dataSamples[j],]
        temp_ct <- temp_samp$CT - controlAvgs[i]
        
        #raise everything to the second pwoer
        temp_ct <- 2^temp_ct
        dataAverage[i,j] <- mean(temp_ct, na.rm = T)
        dataSd[i,j] <- sd(temp_samp$CT, na.rm = T)
        dataSamplesOrder[i,j] <- as.character(dataSamples[j])
        dataTargetsOrder[i,j] <- as.character(targets[i])
      }
    }
    
    #Flatten everything down and stick it into a dataframe
    outputCtrl <- data.frame(Sample = as.character(as.vector(ctrlSamplesOrder)), 
                             Target = as.character(as.vector(ctrlTargetsOrder)), 
                             CT = as.vector(ctrlAverage), 
                             SD = as.vector(ctrlSd))
    outputData <- data.frame(Sample = as.character(as.vector(dataSamplesOrder)), 
                             Target = as.character(as.vector(dataTargetsOrder)), 
                             CT = as.vector(dataAverage), 
                             SD = as.vector(dataSd))
    
    
    #flatten them all into the dataframe for the output
    output <- rbind(outputCtrl, outputData)
    output$Sample <- factor(output$Sample, levels = unique(output$Sample[order(as.character(output$Sample))]))
  }
  
  #plot everything out
  pdf(paste0(gsub('(*.).csv', '\\1', filename), '.pdf'), width = 12)
  p <- ggplot(data=output, aes(x=Sample, y=CT, fill = Target), alpha = Sample) 
  p <- p + geom_bar(stat='identity', position ='dodge') 
  p <- p + scale_fill_hue(l=75, c=100) 
  p <- p + ggtitle(gsub('(*.).csv', '\\1', filename)) 
  p <- p + geom_errorbar(aes(ymin=CT-SD, ymax=CT+SD), width=0.2, position = position_dodge(.9)) 
  p <- p + guides(fill=guide_legend(title="Target"))
  print(p)
  dev.off()
  return(output)
}
#Stick everything together into a matrix of SD and Means




