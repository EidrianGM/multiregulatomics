# Load libraries
setwd("/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/")
library(openxlsx)
# library(ggVennDiagram)
# library(Biobase)
# library(limma)
# library(ggplot2)
# library(circlize)
# library(PCAtools)
# library("scatterplot3d")
# library(RColorBrewer)
# library(ggfortify)
# library(EnhancedVolcano)
# library(stringr)
# library(gdata)
# library(corrplot)
# library(ComplexHeatmap)
# library(dplyr)
# library(fgsea)
# library(httr)
# library(jsonlite)
# library(RCurl)
# library(plyr)
# library(reshape2)
# source("scripts/GC4libR.R")

## Helper function to save tables in TSV and/or Excel
saveTablesTsvExc <- function(DF,outdir,completeNdedup=T,excel=T,bycompleteFC=T,rownames=F){
  dir.create(outdir,F,T)
  outname <- deparse(substitute(DF))
  if (class(DF) == "data.frame"){
    if (nrow(DF) == 0){
      print("Empty DF")
    }else{
      outFile <- paste0(file.path(outdir,outname),".tsv")
      if (completeNdedup){
        #DF <- DF[complete.cases(DF),]
        if(bycompleteFC){
          fcORnetCol <- grep("log2|netChanges",colnames(DF))
          DF <- DF[!is.na(DF[,fcORnetCol]),]
          DF <- DF[!is.na(DF[order(DF[,fcORnetCol],decreasing = T),]),]
        }
        DF <- DF[!duplicated(DF),]
      }
      write.table(DF,outFile,quote = F,sep = '\t',col.names = T,row.names = rownames)
      if (excel){
        write.xlsx(DF,gsub(".tsv",".xlsx",outFile), asTable = FALSE, overwrite = TRUE,rowNames = rownames)
      }
    }
  } else{
    myvector <- unique(na.omit(DF))
    if (length(myvector) == 0){
      print("Empty Vector")
    }else{
      outFile <- paste0(file.path(outdir,outname),".txt")
      write(myvector,outFile,sep = "\n")
    }
  }
}

## Get Results of GC4 and simplify the output.
summaryGC4results <- function(GC4res,outfile){
  allQC <- c()
  allRes <- c()
  for (res in names(GC4res$quality_controls[[1]])){
    resX <- GC4res$quality_controls[[1]][[res]]
    annot <- gsub(".*-","",res)
    if (res == "notInDB"){
      allQC <- rbind(allQC,c("notInDB","","","","",
                             length(resX$invalidInput),
                             paste(resX$invalidInput,collapse = " "),
                             paste(resX$invalidUniverse,collapse = " "),
                             paste(resX$notMapped,collapse = " ")))
    }else{
      allQC <- rbind(allQC,c(annot,length(resX$annotated),length(resX$noAnnotated),
                             paste(resX$noAnnotated,collapse = " "),resX$universe,
                             "",
                             "",
                             "",
                             ""))
    }
  }
  allQC <- as.data.frame(allQC)
  colnames(allQC) <- c("annotation","annotated","noannotated","noannotatedlist","universe","NinvalidInput","invalidInput","invalidUniverse","notMapped")
  for (res in 1:length(GC4res$stats_tables)){
    annot <- gsub(".*-","",names(GC4res$stats_tables)[res])
    resX <- GC4res$stats_tables[[res]]
    resX$annotation <- annot
    allRes <- rbind(allRes,resX)
  }
  return(list(qc=allQC,enr=allRes))
}

completeNdedup <- function(DF){
  if (class(DF) == "data.frame"){
    return(DF[complete.cases(DF),])
  }else{
    return(unique(na.omit(DF)))
  }
}

getProteinsMatrixFromRAW <- function(omeFile){
  myDF <- read.delim(omeFile,sep = ',')
  rawDCol <- grep("Raw",colnames(myDF)) # Raw Normalized Values are not the ones!
  specDCol <- grep("Spectral",colnames(myDF))
  expresionMatrix <- myDF[rawDCol:(specDCol-1)] 
  samplesNames <- paste(gsub("_.*","",expresionMatrix[2,]),gsub(".*\\s","",expresionMatrix[1,]),sep = '_')
  colnames(expresionMatrix) <- samplesNames
  expresionMatrix <- expresionMatrix[3:nrow(expresionMatrix),]
  expresionMatrix[is.na(expresionMatrix)] <- 0
  expresionMatrix <- type.convert(expresionMatrix)
  expresionMatrixProteins <- aggregate(expresionMatrix, by=list(myDF$X.10[3:nrow(myDF)]), FUN = sum)
  colnames(expresionMatrixProteins)[1] <- 'proteins'
  return(expresionMatrixProteins)
}
