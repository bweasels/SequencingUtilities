setwd('/Users/benku/Desktop/IFD WBC FastQ/')

files <- list.files(recursive = T)
files <- files[grep('*.fastq.gz', files)]
fullPath <- '/Users/benku/Desktop/IFD WBC FastQ/'
files <- paste0(fullPath, files)

dir.create('output_star')

samples <- unique(gsub('.*NS112-(.*)_S.*', '\\1', files))

for (i in 1:length(samples)){
  #make the individual output folder for each run
  outputDir <- paste0('output_star/',samples[i], '/')
  dir.create(outputDir)
  
  #for each run get a list of the R1 and R2 files
  tempFiles <- files[grep(samples[i], files)]
  tempR1 <- tempFiles[grep('.*_R1_001.fastq.gz', tempFiles)]
  tempR2 <- tempFiles[grep('.*_R2_001.fastq.gz', tempFiles)]
  
  #Make the requsite header
  sink(paste0(samples[i], '-Staralign.lsf'), append =FALSE)
  cat('#!/bin/bash\n')
  cat(paste0('#BSUB -J ', samples[i]), '\n')
  cat(paste0('#BSUB -o output/', samples[i],'-%J.out'), '\n')
  cat(paste0('#BSUB -e output/', samples[i],'-%J.err'), '\n\n')
  
  #Print the variable names
  cat(paste0("samp='", samples[i],"'\n"))
  cat("gFasta='/PHShome/bw792/starIndex/ercc_gfp_hg38_nochr.fa'\n")
  cat("sjdbGTF='/PHShome/bw792/StarIndex/Homo_sapiens.GRCh38.79.gtf'\n")
  cat("ercGTF='/PHShome/bw792/starIndex/hg38_ercc.gtf'\n")
  cat(paste0("countsOutput='./",samples[i],".txt'\n\n"))
  
  #Print the output filesystems and load the modules
  cat("bamDir='./output/'\n")
  cat("bamInput=$bamDir$samp'Aligned.sortedByCoord.out.bam'\n\n")
  
  cat('module load fastqc/0.11.2\n\n')
  cat('zcat ', tempR1, ' | fastqc stdin\n')
  cat('zcat ', tempR2, ' | fastqc stdin\n\n')
  
  cat('module load aryee/star-2.4.0h\n')
  #make the command for star aligner
  
  cat("STAR --outTmpDir './temp_'$samp",
      "--genomeDir ~/starIndex/",
      "--genomeFastaFiles $gFasta",
      "--sjdbGTFfile $sjdbGTF",
      "--outSAMtype BAM SortedByCoordinate",
      "--outFileNamePrefix $bamDir$samp",
      "--readFilesCommand zcat",
      "--readFilesIn ")
  for (j in 1:length(tempR1)){
    cat(tempR1[j])
    if (j != length(tempR1)){
      cat(',')
    }
  }
  cat(' ')
  for (j in 1:length(tempR2)){
    cat(tempR2[j])
    if (j != length(tempR2)){
      cat(',')
    }
  }
  
  #Step 2 generate the counts table
  cat('\n\n')
  cat('module load python/2.7.3\n')
  cat('module load htseq/0.6.1\n')
  cat('module load samtools-0.1.18\n\n')
  
  #Generate the counts table
  cat('samtools view $bamInput | htseq-count',
      '-r pos', 
      '-m intersection-strict',
      '--stranded=no',
      '-i gene_name - $ercGTF > $countsOutput')
  sink()
}
