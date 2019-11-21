setwd('/data/rama/labMembers/bw792/mechanomutationFastqs/pool19/output')

#Get the files and make the output folders
files <- list.files(recursive = T)
files <- files[grep('*.sortedByCoord.out.bam', files)]
fullPath <- '/data/rama/labMembers/bw792/mechanomutationFastqs/pool19/output/'

#get the sample ids
samples <- unique(gsub('(.*)Aligned.sortedByCoord.out.bam', '\\1', files))

for (i in 1:length(samples)){  
  #Make the requsite header
  sink(paste0(samples[i], '-SortByName.lsf'), append =FALSE)
  cat('#!/bin/bash\n')
  cat(paste0('#BSUB -J ', samples[i]), '\n')
  cat(paste0('#BSUB -o output/', samples[i],'-%J.out'), '\n')
  cat(paste0('#BSUB -e output/', samples[i],'-%J.err'), '\n\n')
  
  #Print the variable names
  print(paste0("file='", fullPath, files[i],"'\n"))
  print(paste0("sortedFile='", fullPath, samples[i], '.Aligned.sortedByName.out.bam'))
  print('module load samtools/1.4.1\n\n')
  print("samtools sort -n -o $sortedFile $file\n")
}
