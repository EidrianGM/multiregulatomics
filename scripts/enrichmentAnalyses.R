library(httr)
library(jsonlite)
library(RCurl)
library(plyr)
source("GC4libR.R")
library(openxlsx)

GSAlists <- list.files("/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/data/geneLists/GSA",pattern = "*.txt",all.files = T,full.names = T)
for (GSAlist in GSAlists){
  siggenes <- read.delim(GSAlist,header = F)[,1]
  enrFile <- gsub(".txt",".Enr.tsv",gsub("geneLists","enrichmentRes",GSAlist),fixed = T)
  qcFile <- gsub(".txt",".QCEnr.tsv",gsub("geneLists","enrichmentRes",GSAlist),fixed = T)
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


####
library(reshape2)

gobpDB <- read.delim("/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/data/databases/go_bp-559292taxId_GeneCodis4.tsv")
gobpDB$db <- "(GO BP)"
goccDB <- read.delim("/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/data/databases/go_cc-559292taxId_GeneCodis4.tsv")
goccDB$db <- "(GO CC)"
gomfDB <- read.delim("/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/data/databases/go_mf-559292taxId_GeneCodis4.tsv")
gomfDB$db <- "(GO MF)"
keggDB <- read.delim("/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/data/databases/kegg-559292taxId_GeneCodis4.tsv")
keggDB$db <- "(KEGG)"
panthDB <- read.delim("/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/data/databases/panther-559292taxId_GeneCodis4.tsv")
panthDB$db <- "(Panther)"
allDBs <- rbind(gobpDB,goccDB,gomfDB,keggDB,panthDB)

annotsInfo <- read.delim("/home/eidriangm/Desktop/annotation_info_table.tsv")
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
    maxGeneset <- max(c(maxGeneset,length(gseaDB[[db]][[term]])))
    minGeneset <- min(c(minGeneset,length(gseaDB[[db]][[term]])))
    write(paste(c(term,"na",gseaDB[[db]][[term]]),collapse = "\t"),paste0("/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/data/databases/",gsub("\\(|\\)","",db),".gmt"),append = T)
  }
}

genesprots <- allDBs$synonyms[allDBs$term == term]
term <- "high-affinity zinc transmembrane transporter activity GO:0000006 (GO MF)"

gseaDB[[db]][[term]]

View(read.delim("/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/data/databases/all.gmt",col.names = F))

mRNABPomeUVfcSA <- read.delim("/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/geneLists/GSEA/mRNABPomeUVfcSA.tsv")
mRNABPomeUVfcSA <- mRNABPomeUVfcSA[,c("UniprotACC","mRNABPomeUV.log2ratio_PolyARNAUVwithSA")]
mRNABPomeUVfcSA <- aggregate(mRNABPomeUVfcSA$mRNABPomeUV.log2ratio_PolyARNAUVwithSA, list(mRNABPomeUVfcSA$UniprotACC), FUN=mean) 
colnames(mRNABPomeUVfcSA) <- c("UniprotACC","mRNABPomeUV.log2ratio_PolyARNAUVwithSA")
mRNABPomeUVfcSA <- mRNABPomeUVfcSA[!duplicated(mRNABPomeUVfcSA),]
mRNABPomeUVfcSA <- mRNABPomeUVfcSA[order(mRNABPomeUVfcSA$mRNABPomeUV.log2ratio_PolyARNAUVwithSA),]
ibtissamDir <- "/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/toIbtissam"
saveTablesTsvExc(mRNABPomeUVfcSA,ibtissamDir)

mRNABPomeUVfcSA <- read.delim("/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/geneLists/GSEA/mRNABPomeUVfcSA.tsv")
mRNABPomeUVfcSAgene <- mRNABPomeUVfcSA[,c("GeneName","mRNABPomeUV.log2ratio_PolyARNAUVwithSA")]
mRNABPomeUVfcSA <- aggregate(mRNABPomeUVfcSAgene$mRNABPomeUV.log2ratio_PolyARNAUVwithSA, list(mRNABPomeUVfcSAgene$GeneName), FUN=mean) 
colnames(mRNABPomeUVfcSAgene) <- c("GeneName","mRNABPomeUV.log2ratio_PolyARNAUVwithSA")
mRNABPomeUVfcSAgene <- mRNABPomeUVfcSAgene[!duplicated(mRNABPomeUVfcSAgene),]
mRNABPomeUVfcSAgene <- mRNABPomeUVfcSAgene[order(mRNABPomeUVfcSAgene$mRNABPomeUV.log2ratio_PolyARNAUVwithSA),]
ibtissamDir <- "/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/toIbtissam"
saveTablesTsvExc(mRNABPomeUVfcSAgene$GeneName,ibtissamDir)


mRNABPomeUVfcSAfull <- merge(mRNABPomeUVfcSA,yeastGenesProtMap,all.x = T)

mRNABPomeUVfcSAfull[,]

all(!is.na(mRNABPomeUVfcSAfull$ORFgene))
all(!is.na(mRNABPomeUVfcSAfull$GeneName))
all(!is.na(mRNABPomeUVfcSAfull$Gene.Names..synonym.))

gsea_mRNABPomeUVfcSA <- as.array(mRNABPomeUVfcSA$mRNABPomeUV.log2ratio_PolyARNAUVwithSA) 
names(gsea_mRNABPomeUVfcSA) <- mRNABPomeUVfcSA$UniprotACC



mRNABPomeUVfcSAfGSEAResBP <- fgsea(pathways = gseaDB[["(GO BP)"]], stats = gsea_mRNABPomeUVfcSA, minSize  = minGeneset, maxSize  = biggestGeneset)
mRNABPomeUVfcSAfGSEAResBP$db <- "GO BP"
mRNABPomeUVfcSAfGSEAResCC <- fgsea(pathways = gseaDB[["(GO CC)"]], stats = gsea_mRNABPomeUVfcSA, minSize  = minGeneset, maxSize  = biggestGeneset)
mRNABPomeUVfcSAfGSEAResCC$db <- "GO CC"
mRNABPomeUVfcSAfGSEAResMF <- fgsea(pathways = gseaDB[["(GO MF)"]], stats = gsea_mRNABPomeUVfcSA, minSize  = minGeneset, maxSize  = biggestGeneset)
mRNABPomeUVfcSAfGSEAResMF$db <- "GO MF"
mRNABPomeUVfcSAfGSEAResP <- fgsea(pathways = gseaDB[["(Panther)"]], stats = gsea_mRNABPomeUVfcSA, minSize  = minGeneset, maxSize  = biggestGeneset)
mRNABPomeUVfcSAfGSEAResP$db <- "Panther"
mRNABPomeUVfcSAfGSEAResK <- fgsea(pathways = gseaDB[["(KEGG)"]], stats = gsea_mRNABPomeUVfcSA, minSize  = minGeneset, maxSize  = biggestGeneset)
mRNABPomeUVfcSAfGSEAResK$db <- "KEGG"

mRNABPomeUVfcSAfGSEARes <- rbind(mRNABPomeUVfcSAfGSEAResBP,mRNABPomeUVfcSAfGSEAResCC,mRNABPomeUVfcSAfGSEAResMF,mRNABPomeUVfcSAfGSEAResP,mRNABPomeUVfcSAfGSEAResK)

mRNABPomeFAXfcSAOR <- read.delim("/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/geneLists/GSEA/mRNABPomeFAXfcSA.tsv")
mRNABPomeFAXfcSA <- mRNABPomeFAXfcSAOR[,c("UniprotACC","mRNABPomeFAX.log2ratio_PolyARNAFAXwithSA")]
mRNABPomeFAXfcSA <- aggregate(mRNABPomeFAXfcSA$mRNABPomeFAX.log2ratio_PolyARNAFAXwithSA, list(mRNABPomeFAXfcSA$UniprotACC), FUN=mean) 
colnames(mRNABPomeFAXfcSA) <- c("UniprotACC","mRNABPomeFAX.log2ratio_PolyARNAFAXwithSA")
mRNABPomeFAXfcSA <- mRNABPomeFAXfcSA[!duplicated(mRNABPomeFAXfcSA),]
mRNABPomeFAXfcSA <- mRNABPomeFAXfcSA[order(mRNABPomeFAXfcSA$mRNABPomeFAX.log2ratio_PolyARNAFAXwithSA),]

mRNABPomeFAXfcSAgeneName <- merge(mRNABPomeFAXfcSA,mRNABPomeFAXfcSAOR,all.x=T)

[,c("GeneName","mRNABPomeFAX.log2ratio_PolyARNAFAXwithSA")]

saveTablesTsvExc(mRNABPomeUVfcSA,ibtissamDir)
gsea_mRNABPomeFAXfcSA <- as.array(mRNABPomeFAXfcSA$mRNABPomeFAX.log2ratio_PolyARNAFAXwithSA) 
names(gsea_mRNABPomeFAXfcSA) <- mRNABPomeFAXfcSA$UniprotACC

mRNABPomeFAXfcSAfGSEAResBP <- fgsea(pathways = gseaDB[["(GO BP)"]], stats = gsea_mRNABPomeFAXfcSA, minSize  = minGeneset, maxSize  = biggestGeneset)
mRNABPomeFAXfcSAfGSEAResBP$db <- "GO BP"
mRNABPomeFAXfcSAfGSEAResCC <- fgsea(pathways = gseaDB[["(GO CC)"]], stats = gsea_mRNABPomeFAXfcSA, minSize  = minGeneset, maxSize  = biggestGeneset)
mRNABPomeFAXfcSAfGSEAResCC$db <- "GO CC"
mRNABPomeFAXfcSAfGSEAResMF <- fgsea(pathways = gseaDB[["(GO MF)"]], stats = gsea_mRNABPomeFAXfcSA, minSize  = minGeneset, maxSize  = biggestGeneset)
mRNABPomeFAXfcSAfGSEAResMF$db <- "GO MF"
mRNABPomeFAXfcSAfGSEAResP <- fgsea(pathways = gseaDB[["(Panther)"]], stats = gsea_mRNABPomeFAXfcSA, minSize  = minGeneset, maxSize  = biggestGeneset)
mRNABPomeFAXfcSAfGSEAResP$db <- "Panther"
mRNABPomeFAXfcSAfGSEAResK <- fgsea(pathways = gseaDB[["(KEGG)"]], stats = gsea_mRNABPomeFAXfcSA, minSize  = minGeneset, maxSize  = biggestGeneset)
mRNABPomeFAXfcSAfGSEAResK$db <- "KEGG"

mRNABPomeFAXfcSAfGSEARes <- rbind(mRNABPomeFAXfcSAfGSEAResBP,mRNABPomeFAXfcSAfGSEAResCC,mRNABPomeFAXfcSAfGSEAResMF,mRNABPomeFAXfcSAfGSEAResP,mRNABPomeFAXfcSAfGSEAResK)

mRNABPomeUVfcSAfGSEARes$leadingEdge <- unlist(lapply(mRNABPomeUVfcSAfGSEARes$leadingEdge, function(x) return(paste(x, collapse = ", "))))
mRNABPomeFAXfcSAfGSEARes$leadingEdge <- unlist(lapply(mRNABPomeFAXfcSAfGSEARes$leadingEdge, function(x) return(paste(x, collapse = ", "))))

mRNABPomeUVfcSAgseaFile <- "/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/enrichmentRes/GSEA/mRNABPomeUVfcSA.tsv"
mRNABPomeUVfcSAfGSEARes <- mRNABPomeUVfcSAfGSEARes[order(mRNABPomeUVfcSAfGSEARes$padj),]
write.table(mRNABPomeUVfcSAfGSEARes,mRNABPomeUVfcSAgseaFile,quote = F,sep = '\t',col.names = T,row.names = F)
mRNABPomeFAXfcSAgseaFile <- "/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/enrichmentRes/GSEA/mRNABPomeFAXfcSA.tsv"
mRNABPomeFAXfcSAfGSEARes <- mRNABPomeFAXfcSAfGSEARes[order(mRNABPomeFAXfcSAfGSEARes$padj),]
write.table(mRNABPomeFAXfcSAfGSEARes,mRNABPomeFAXfcSAgseaFile,quote = F,sep = '\t',col.names = T,row.names = F)


####



