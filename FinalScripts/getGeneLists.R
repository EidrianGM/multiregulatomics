wholeDFFile <- "/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/FinalData/allDataDF.tsv"
wholeDF <- read.delim(wholeDFFile,quote = "")

################################################################################
######### Generation of GeneLists for Enrichment Analyses ORA and GSEA #########
################################################################################
outdir <- "/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/FinalData/GeneLists"

IDs <- c("UniprotACC") # "GeneName","P<rotein.names") # rbpomeID <- IDs
DEGpvalCut <- 0.01
DEGFCCut <- 3
protsPvalCut <- 0.05
fcCut   <- 1

# ORA
getIDsList <- function(DF,logicalRows,ID){
  DF <- DF[logicalRows,ID]
  if (class(DF) == "character"){
    DF <- unique(DF) 
    print(paste(length(DF),"are significant"))
  }else{
    DF <- DF[!duplicated(DF),]
    print(paste(nrow(DF),"are significant")) 
  }
  return(DF)
} 

# 2- DEGS SA ORA 
SAdegsQSig <- wholeDF$DEGs.adj.pval.SA.YPD < DEGpvalCut
SAdegsUpFC <- wholeDF$DEGs.log2FC.SA.YPD > DEGFCCut
SAdegsDwFC <- wholeDF$DEGs.log2FC.SA.YPD < -DEGFCCut

SADEGsUpFCUniProts <- getIDsList(wholeDF,SAdegsQSig & SAdegsUpFC,IDs)
SADEGsDwFCUniProts <- getIDsList(wholeDF,SAdegsQSig & SAdegsDwFC,IDs)

outdir <- "/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/FinalData/GeneLists/ORAseparatedUpnDw"
saveTablesTsvExc(SADEGsUpFCUniProts,outdir)
saveTablesTsvExc(SADEGsDwFCUniProts,outdir)

outdir <- "/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/FinalData/GeneLists/ORA"
SADEGsUniProts <- c(SADEGsUpFCUniProts,SADEGsDwFCUniProts)
saveTablesTsvExc(SADEGsUniProts,outdir)

# 2- FAX SA ORA 
SAFAXprotQSig <- wholeDF$ProteomeFAX.qValue_FAXwithSA < protsPvalCut
SAFAXrbpmQsig <- wholeDF$RBPomeFAX.qValue_PolyARNAFAXwithSA < protsPvalCut
SAFAXprotUpFC <- wholeDF$ProteomeFAX.log2ratio_FAXwithSA > fcCut
SAFAXprotDwFC <- wholeDF$ProteomeFAX.log2ratio_FAXwithSA < -fcCut
SAFAXrbpmUpFC <- wholeDF$RBPomeFAX.log2ratio_PolyARNAFAXwithSA > fcCut
SAFAXrbpmDwFC <- wholeDF$RBPomeFAX.log2ratio_PolyARNAFAXwithSA < -fcCut
SAFAXupNC <- wholeDF$FAXnetchangesSA > fcCut
SAFAXdwNC <- wholeDF$FAXnetchangesSA < -fcCut

SAFAXprotUpFCUniProts <- getIDsList(wholeDF, SAFAXprotQSig & SAFAXprotUpFC, IDs)
SAFAXprotDwFCUniProts <- getIDsList(wholeDF, SAFAXprotQSig & SAFAXprotDwFC, IDs)
SAFAXrbpmUpFCUniProts <- getIDsList(wholeDF, SAFAXrbpmQsig & SAFAXrbpmUpFC, IDs)
SAFAXrbpmDwFCUniProts <- getIDsList(wholeDF, SAFAXrbpmQsig & SAFAXrbpmDwFC, IDs)
SAFAXnetcupNCUniProts <- getIDsList(wholeDF, SAFAXrbpmQsig & SAFAXupNC, IDs)
SAFAXnetcdwNCUniProts <- getIDsList(wholeDF, SAFAXrbpmQsig & SAFAXdwNC, IDs)

outdir <- "/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/FinalData/GeneLists/ORAseparatedUpnDw"
saveTablesTsvExc(SAFAXprotUpFCUniProts,outdir)
saveTablesTsvExc(SAFAXprotDwFCUniProts,outdir)
saveTablesTsvExc(SAFAXrbpmUpFCUniProts,outdir)
saveTablesTsvExc(SAFAXrbpmDwFCUniProts,outdir)
saveTablesTsvExc(SAFAXnetcupNCUniProts,outdir)
saveTablesTsvExc(SAFAXnetcdwNCUniProts,outdir)

outdir <- "/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/FinalData/GeneLists/ORA"
SAFAXprotUniProts <- c(SAFAXprotUpFCUniProts,SAFAXprotDwFCUniProts)
SAFAXrbpmFCUniProts <- c(SAFAXrbpmUpFCUniProts,SAFAXrbpmDwFCUniProts)
SAFAXnetcUniProts <- c(SAFAXnetcupNCUniProts,SAFAXnetcdwNCUniProts)
saveTablesTsvExc(SAFAXprotUniProts,outdir)
saveTablesTsvExc(SAFAXrbpmFCUniProts,outdir)
saveTablesTsvExc(SAFAXnetcUniProts,outdir)

# 2- UVX SA ORA
SAUVXprotQSig <- wholeDF$ProteomeUVX.qValue_Condition4 < protsPvalCut
SAUVXrbpmQsig <- wholeDF$RBPomeUVX.qValue_PolyARNAUVwithSA < protsPvalCut
SAUVXprotUpFC <- wholeDF$ProteomeUVX.log2ratio_Condition4 > fcCut
SAUVXprotDwFC <- wholeDF$ProteomeUVX.log2ratio_Condition4 < -fcCut
SAUVXrbpmUpFC <- wholeDF$RBPomeUVX.log2ratio_PolyARNAUVwithSA > fcCut
SAUVXrbpmDwFC <- wholeDF$RBPomeUVX.log2ratio_PolyARNAUVwithSA < -fcCut
SAUVXupNC <- wholeDF$UVXnetchangesSA > fcCut
SAUVXdwNC <- wholeDF$UVXnetchangesSA < -fcCut

SAUVXprotUpFCUniProts <- getIDsList(wholeDF,SAUVXprotQSig & SAUVXprotUpFC,IDs)
SAUVXprotDwFCUniProts <- getIDsList(wholeDF,SAUVXprotQSig & SAUVXprotDwFC,IDs)
SAUVXrbpmUpFCUniProts <- getIDsList(wholeDF,SAUVXrbpmQsig & SAUVXrbpmUpFC,IDs)
SAUVXrbpmDwFCUniProts <- getIDsList(wholeDF,SAUVXrbpmQsig & SAUVXrbpmDwFC,IDs)
SAUVXnetcupNCUniProts <- getIDsList(wholeDF,SAUVXrbpmQsig & SAUVXupNC,IDs)
SAUVXnetcdwNCUniProts <- getIDsList(wholeDF,SAUVXrbpmQsig & SAUVXdwNC,IDs)

outdir <- "/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/FinalData/GeneLists/ORAseparatedUpnDw"
saveTablesTsvExc(SAUVXprotUpFCUniProts,outdir)
saveTablesTsvExc(SAUVXprotDwFCUniProts,outdir)
saveTablesTsvExc(SAUVXrbpmUpFCUniProts,outdir)
saveTablesTsvExc(SAUVXrbpmDwFCUniProts,outdir)
saveTablesTsvExc(SAUVXnetcupNCUniProts,outdir)
saveTablesTsvExc(SAUVXnetcdwNCUniProts,outdir)

outdir <- "/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/FinalData/GeneLists/ORA"
SAUVXprotUniProts <- c(SAUVXprotUpFCUniProts,SAUVXprotDwFCUniProts)
SAUVXrbpmFCUniProts <- c(SAUVXrbpmUpFCUniProts,SAUVXrbpmDwFCUniProts)
SAUVXnetcUniProts <- c(SAUVXnetcupNCUniProts,SAUVXnetcdwNCUniProts)

saveTablesTsvExc(SAUVXprotUniProts,outdir)
saveTablesTsvExc(SAUVXrbpmFCUniProts,outdir)
saveTablesTsvExc(SAUVXnetcUniProts,outdir)


# GSEA 
UVXnetchangesSA

getGSEArankDF <- function(DF,dataype,crosslink,treatment,IDs,rankCol="log2"){
  #crosslink <- ifelse(dataype=="DEGs","",crosslink)
  rankcol <- colnames(DF)[grep(paste0(dataype,crosslink,"\\.",rankCol,".*",treatment,".*"),colnames(DF))]
  if (rankCol != "log2"){
    rankcol <- colnames(DF)[grep(paste0(crosslink,rankCol,treatment),colnames(DF))]
  }
  if (length(rankcol) == 0){
    treatment <- gsub("DRR","Condition2",gsub("H2O2|H202","Condition3",gsub("SA","Condition4",treatment)))
    rankcol <- colnames(DF)[grep(paste0(dataype,crosslink,"\\.log2ratio.*",treatment),colnames(DF))]
  }
  print(rankcol)
  DF <- DF[order(DF[,rankcol],decreasing = T),c(IDs,rankcol)]
  DF <- DF[!duplicated(DF),]
  DF <- DF[complete.cases(DF),]
  print(paste(nrow(DF),"to GSEA"))
  return(DF)
} 
outdir <- "/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/FinalData/GeneLists/GSEA"
DF <- wholeDF

dataype <- "DEGs"; crosslink <- ""; treatment <- "SA"
SAdegsGSEA <- getGSEArankDF(wholeDF,dataype,crosslink,treatment,IDs)
saveTablesTsvExc(SAdegsGSEA,outdir,completeNdedup = T,excel = F,bycompleteFC = F,rownames = F)

dataype <- "Proteome"; crosslink <- "FAX"; treatment <- "SA"
SAProteomeFAXgsea <- getGSEArankDF(wholeDF,dataype,crosslink,treatment,IDs)
saveTablesTsvExc(SAProteomeFAXgsea,outdir,completeNdedup = F,excel = F,bycompleteFC = F,rownames = F)

dataype <- "RBPome"; crosslink <- "FAX"; treatment <- "SA"
SARBPomeFAXgsea <- getGSEArankDF(wholeDF,dataype,crosslink,treatment,IDs)
saveTablesTsvExc(SARBPomeFAXgsea,outdir,completeNdedup = F,excel = F,bycompleteFC = F,rownames = F)

SAnetchangesFAXgsea <- getGSEArankDF(wholeDF,dataype,crosslink,treatment,IDs,"netchanges")
saveTablesTsvExc(SAnetchangesFAXgsea,outdir,completeNdedup = F,excel = F,bycompleteFC = F,rownames = F)

dataype <- "Proteome"; crosslink <- "UVX"; treatment <- "SA"
SAProteomeUVXgsea <- getGSEArankDF(wholeDF,dataype,crosslink,treatment,IDs)
saveTablesTsvExc(SAProteomeUVXgsea,outdir,completeNdedup = F,excel = F,bycompleteFC = F,rownames = F)

dataype <- "RBPome"; crosslink <- "UVX"; treatment <- "SA"
SARBPomeUVXgsea <- getGSEArankDF(wholeDF,dataype,crosslink,treatment,IDs)
saveTablesTsvExc(SARBPomeUVXgsea,outdir,completeNdedup = F,excel = F,bycompleteFC = F,rownames = F)

SAnetchangesUVXgsea <- getGSEArankDF(wholeDF,dataype,crosslink,treatment,IDs,"netchanges")
saveTablesTsvExc(SAnetchangesUVXgsea,outdir,completeNdedup = F,excel = F,bycompleteFC = F,rownames = F)

