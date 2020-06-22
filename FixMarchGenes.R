##############################################################
##Cmd: mkala001.R --vanilla
##Purpose: Fix September/March/December genes
##Output: The RNA seq table with dates fixed
##Notes: ugh this is such a stupid excersize - fix yo shit Excel
##############################################################

fixMarchGenes <- function(reads_table){
  #as a result, the rownames is actually the first column
  rownames <- reads_table[,1]
  
  #Finds the excel corrected dates in the format dd-Mmm (ie 14-Sep) and flips it
  rownames <- toupper(gsub("^(\\d+)-([A-Za-z0-9]+)$","\\2\\1", rownames))
  #Outputs in this format (SEP14)
  
  #convert to the proper gene names if inverted form is incorrect
  rownames <- gsub("^SEP([0-9]+)","SEPT\\1", rownames)
  rownames <- gsub("^MAR([0-9]+)","MARCH\\1", rownames)
  
  #get all the MARCH genes, because the first MARCH1 and MARCH2 are actually
  #MARC1 and MARC2
  marchGenes <- grep("^MARCH([0-9]+)$", rownames)
  
  #variables to save the second MARCH1 and MARCH2
  fMARC1 <- FALSE
  fMARC2 <- FALSE
  
  if (sum(grepl("^MARC[0-9]", rownames)) > 0) { ##IF for some reason MARC didn't become march, don't touch them
    fMARC1 <- T
    fMARC2 <- T
  }
  
  #go through the marchGenes and change the first MARCH1 & MARCH2 to MARC1 & MARC2 
  for (index in marchGenes){
    if ((rownames[index]=='MARCH1') & (!fMARC1)){
      rownames[index] <- 'MARC1'
      fMARC1 <- TRUE
    }
    if ((rownames[index]=='MARCH2') & (!fMARC2)){
      rownames[index] <- 'MARC2'
      fMARC2 <- TRUE
    }
  }
  
  #replace the default names with the correct names
  rownames(reads_table) <- rownames
  
  #remove the column with the names
  reads_table<- subset(reads_table, select =-c(X))
}
