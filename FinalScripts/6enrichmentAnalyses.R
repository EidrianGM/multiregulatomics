source("scripts/GC4libR.R")
library(fgsea)
setwd("/home/eidriangm/Desktop/toDo/surrey/multiregulatomics")

################################################################################
############### 1. OVERREPRESENTATION ANALYSIS WITH GENECODIS4 #################
################################################################################
doGC4ORA <- function(dirWithGeneLists,orthologs=TRUE){
  GSAlists <- list.files(dirWithGeneLists,pattern = "*.txt",all.files = T,full.names = T)
  outdir <- gsub("GeneLists","EnrichmentResults",dirWithGeneLists)
  dir.create(outdir,recursive = T,showWarnings = F)
  for (GSAlist in GSAlists){
    siggenes <- read.delim(GSAlist,header = F)[,1]
    if (any(grepl(";",siggenes)) & orthologs){
      siggenes <- unique(unlist(strsplit(siggenes,";")))
    }
    enrFile <- gsub(".txt",".Enr.tsv",gsub("GeneLists","EnrichmentResults",GSAlist),fixed = T) # ORA results
    qcFile <- gsub(".txt",".QCEnr.tsv",gsub("GeneLists","EnrichmentResults",GSAlist),fixed = T) # Quality Control of ORA
    if (file.exists(qcFile)){next}
    resultado <- launchAnalysis(organism = "Saccharomyces cerevisiae",
                                inputType = "genes",
                                inputQuery = siggenes,
                                annotationsDBs = c("KEGG","GO_BP","GO_CC","GO_MF"),
                                inputCoannotation = "no",
                                universeScope = "annotated",
                                enrichmentStat = "hypergeom",
                                ReportName=gsub(".tsv|.txt","",basename(GSAlist)))
    
    results <- summaryGC4results(resultado)
    write.table(results[["enr"]],enrFile,quote = F,sep = '\t',col.names = T,row.names = F)
    write.table(results[["qc"]],qcFile,quote = F,sep = '\t',col.names = T,row.names = F)
  }
}

dirWithGeneLists <- "FinalData/GeneLists/ORA"
doGC4ORA(dirWithGeneLists,orthologs=TRUE)

dirWithGeneLists <- "FinalData/GeneLists/ORAseparatedUpnDw"
doGC4ORA(dirWithGeneLists,orthologs=TRUE)

################################################################################
################################# 2. GSEA  #####################################

################################################################################
################# 2.1. CREATION OF GSEA DATABASE #################################
################################################################################
annotsInfo <- read.delim("/home/eidriangm/Desktop/annotation_info_table.tsv")

# Download manually annotation files from GeneCodis 4
annotationsFilesFolder <- "data/databases/uniprotGC4"
annotFiles <- list.files(annotationsFilesFolder,pattern = "*.tsv",full.names = T)
allDBs <- c()
for (annotFile in annotFiles){
  DB <- read.delim(annotFile)
  DB$db <- unlist(strsplit(basename(annotFile),"-"))[1]
  allDBs <- rbind(allDBs,DB)
}

allDBs <- merge(allDBs,annotsInfo,all.x=T)
allDBs$term <- paste(allDBs$term,allDBs$annotation_id,allDBs$db,sep = " ")

gseaDB <- list()
maxGeneset <- 1
minGeneset <- 1
for (db in unique(allDBs$db)){
  gseaDB[[db]] <- list()
  allterms <- unique(allDBs$term[allDBs$db == db])
  for (term in allterms){
    gseaDB[[db]][[term]] <- allDBs$synonyms[allDBs$term == term]
    maxGeneset <- max(c(maxGeneset,length(gseaDB[[db]][[term]]))) # Necessary for the analysis 
    minGeneset <- min(c(minGeneset,length(gseaDB[[db]][[term]]))) # later to not filter any annotation
  }
}

################################################################################
################################################################################
library(tidyr)
library(fgsea)


dirname("a/b/c/asda.txt")

doGSEAs <- function(GSEAFile,gseaDB,minGeneset,maxGeneset, orthologs=T) {
  outfile <- gsub("GeneLists","EnrichmentResults",GSEAFile)
  dir.create(dirname(outfile),showWarnings = F)
  if (file.exists(outfile)){
    cat("Done")
    return()
  }
  GSEAdf <- read.delim(GSEAFile)
  if (any(grepl(";",GSEAdf[,1])) & orthologs){
    GSEAdf <- as.data.frame(separate_rows(GSEAdf,1,sep = ";"))
  }
  dupgeprots <- sum(duplicated(GSEAdf[,1]))
  if (dupgeprots > 0){
    cat("Gene/Prots Duplicated: ",dupgeprots, "... collapsing by FC average")
    GSEAdf <- aggregate(GSEAdf[,2], list(GSEAdf[,1]), FUN=mean) 
  }
  colnames(GSEAdf) <- c("GeneProt","log2FC")
  GSEAdf <- GSEAdf[order(GSEAdf$log2FC,decreasing = T),]
  
  inputRank <- as.array(GSEAdf$log2FC)
  names(inputRank) <- GSEAdf$GeneProt
  allGSEAres <- c()
  for (db in names(gseaDB)){
    dbName <- gsub("\\(|\\)","",db)
    gseaRes <- fgsea(pathways = gseaDB[[db]], stats = inputRank, minSize  = minGeneset, maxSize  = maxGeneset, eps=0)
    gseaRes$db <- db
    allGSEAres <- rbind(allGSEAres,gseaRes)
  }
  allGSEAres$leadingEdge <- unlist(lapply(allGSEAres$leadingEdge, function(x) return(paste(x, collapse = ", "))))
  allGSEAres <- allGSEAres[order(allGSEAres$padj,decreasing = F),]
  allGSEAres$pathway <- gsub("\\s*\\w*$","",allGSEAres$pathway)
  write.table(allGSEAres,outfile,quote = F,sep = '\t',col.names = T,row.names = F)
  return(allGSEAres)
}

SAdegsGSEAFile <- "FinalData/GeneLists/GSEA/SAdegsGSEA.tsv"
SAProteomeFAXgseaFile <- "FinalData/GeneLists/GSEA/SAProteomeFAXgsea.tsv"
SARBPomeFAXgseaFile <- "FinalData/GeneLists/GSEA/SARBPomeFAXgsea.tsv"
SAnetchangesFAXgseaFile <- "FinalData/GeneLists/GSEA/SAnetchangesFAXgsea.tsv"
SAdegsGSEAFile <- "FinalData/GeneLists/GSEA/SAdegsGSEA.tsv"
SAProteomeUVXgseaFile <- "FinalData/GeneLists/GSEA/SAProteomeUVXgsea.tsv"
SARBPomeUVXgseaFile <- "FinalData/GeneLists/GSEA/SARBPomeUVXgsea.tsv"
SAnetchangesUVXgseaFile <- "FinalData/GeneLists/GSEA/SAnetchangesUVXgsea.tsv"

SAdegsGSEA <- doGSEAs(SAdegsGSEAFile,gseaDB,minGeneset,maxGeneset)
SAProteomeFAXgsea <- doGSEAs(SAProteomeFAXgseaFile,gseaDB,minGeneset,maxGeneset)

SARBPomeFAXgsea <- doGSEAs(SARBPomeFAXgseaFile,gseaDB,minGeneset,maxGeneset)
SAnetchangesFAXgsea <- doGSEAs(SAnetchangesFAXgseaFile,gseaDB,minGeneset,maxGeneset)

SAdegsGSEA <- doGSEAs(SAdegsGSEAFile,gseaDB,minGeneset,maxGeneset)
SAProteomeUVXgsea <- doGSEAs(SAProteomeUVXgseaFile,gseaDB,minGeneset,maxGeneset)

SARBPomeUVXgsea <- doGSEAs(SARBPomeUVXgseaFile,gseaDB,minGeneset,maxGeneset)
SAnetchangesUVXgsea <- doGSEAs(SAnetchangesUVXgseaFile,gseaDB,minGeneset,maxGeneset)
