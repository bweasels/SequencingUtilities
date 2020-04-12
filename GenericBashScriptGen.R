###This is a flexible function that will hopefully allow me to generate bash scripts easily###
genBashScript <- function(commands, filename){
  
  #Make sure our inputs are correct
  if(typeof(commands)!='character'){
    stop('Commands must be a character vector with the commands in order')
  }
  if(typeof(filename)!='character'){
    stop('Filename must be a character vector with the commands in order')
  }
  if(length(filename)>1){
    stop('Filename must be one element')
  }
  
  #Add .sh to the filename if it doesn't exist already
  if(!grepl('(sh$)|(bash$)', filename)){
    filename <- paste0(filename, '.sh')
  }
  
  sink(filename, append = F)
  cat('#!/bin/bash\n')
  
  for(i in 1:length(commands)){
    cat(paste0(commands[i], '\n'))
  }
  
  #write file
  sink()
}

genLsfScript <- function(commands, filename, jobID){
  
  #Make sure our inputs are correct
  if(typeof(commands)!='character'){
    stop('Commands must be a character vector with the commands in order')
  }
  if(typeof(filename)!='character'){
    stop('Filename must be a character vector with the commands in order')
  }
  if(length(filename)>1){
    stop('Filename must be one element')
  }
  
  #Add .lsf to the filename if it doesn't exist already
  if(!grepl('lsf$', filename)){
    filename <- paste0(filename, '.lsf')
  }
  
  #Make lsf header
  sink(filename, append = F)
  cat('#!/bin/bash\n')
  cat(paste0('#BSUB -J ', jobID), '\n')
  cat(paste0('#BSUB -o output/', jobID,'-%J.out'), '\n')
  cat(paste0('#BSUB -e output/', jobID,'-%J.err'), '\n\n')
  
  for(i in 1:length(commands)){
    cat(paste0(commands[i], '\n'))
  }
  
  sink()
}