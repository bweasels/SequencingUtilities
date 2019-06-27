##As advertised
##Based on the format: /data/rama/lawrence/software/Salmon/Salmon-0.8.2_linux_x86_64/bin/salmon quant -i /data/rama/lawrence/db/GRCh37/Sequence/Transcriptome/Homo_sapiens.GRCh37.75.cdna.all_index -l A -1 {read1.fastq} -2 {read2.fastq} -o {output_dir} -p {num_threads} --seqBias --gcBias
##Many thanks to Adam Langebucher for this gud shit

setwd('/PHShome/bw792/MassGeneralWork/FACS_NK_Cell_RNASeq/')
files <- list.files('.', recursive = T)

files <- files[grepl('*.fastq.gz', files)]

#define the samples, index and salmon call I'll use
samples <- gsub('fastq/NS91-(.*)_L00[0-9]/.*fastq.gz', '\\1', files)
samples <- unique(samples)
index <- '/PHShome/bw792/starIndex/Homo_sapiens.GRCh38_Index'
cmd <- '/data/rama/lawrence/software/Salmon/Salmon-0.8.2_linux_x86_64/bin/salmon'

for (i in 1:length(samples)){
  #make the individual output folder for each run
  #Salmon makes identically named files for each fastq, 
  #so you need to create individual output folders to retain sample separation
  outputDir <- paste0('output/',samples[i], '/')
  dir.create(outputDir)
  
  #for each run get a list of the R1 and R2 files
  tempFiles <- files[grep(samples[i], files)]
  tempR1 <- tempFiles[grep('R1', tempFiles)]
  tempR2 <- tempFiles[grep('R2', tempFiles)]
  
  #Make the requsite header
  sink(paste0(samples[i], '-align.lsf'), append =FALSE)
  cat('#!/bin/bash\n')
  cat(paste0('#BSUB -J ', samples[i]), '\n')
  cat(paste0('#BSUB -o output/', samples[i],'-%J.out'), '\n')
  cat(paste0('#BSUB -e output/', samples[i],'-%J.err'), '\n\n')
  
  #Print the body of the command
  # -l A tells Salmon to automatically detect what kind of fastq files I've handed to it 
  # the docs recommend it unless there's a specific reason to not
  cat(cmd, 'quant -i', index, '-l A -1 ')
  
  #add the R1 and R2 fastq files - multiple fastqs are supported, just seperate with a ' '
  for (j in 1:length(tempR1)){
    cat(tempR1[j])
    if (j != length(tempR1)){
      cat(' ')
    }
  }
  cat(' -2 ')
  for (j in 1:length(tempR2)){
    cat(tempR2[j])
    if (j != length(tempR2)){
      cat(' ')
    }
  }
  
  #Select the output sub folder for your sample
  cat(' -o', paste0('/PHShome/bw792/MassGeneralWork/FACS_NK_Cell_RNASeq/', outputDir))
  
  #-p tells the number of processors to use, and --seqBias and --gcBias corrects for biases
  cat(' -p 6 ')
  cat('--seqBias --gcBias')
  sink()
}

