wholeDFFile <- "/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/FinalData/allDataDF.tsv"
wholeDF <- read.delim(wholeDFFile,quote = "")

################################################################################
######### Generation of GeneLists for Enrichment Analyses ORA and GSEA #########
################################################################################
outdir <- "/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/FinalData/GeneLists"

DEGsIDs <- c("UniprotACC") # Is the one with most significant annotated in GC4
ProtIDs <- c("proteinName") # To also consider orthologues that will be splited for ORA

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

SADEGsUpFC <- getIDsList(wholeDF,which(SAdegsQSig & SAdegsUpFC),DEGsIDs)
SADEGsDwFC <- getIDsList(wholeDF,which(SAdegsQSig & SAdegsDwFC),DEGsIDs)

outdir <- "/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/FinalData/GeneLists/ORAseparatedUpnDw"
saveTablesTsvExc(SADEGsUpFC,outdir)
saveTablesTsvExc(SADEGsDwFC,outdir)

outdir <- "/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/FinalData/GeneLists/ORA"
SADEGs <- c(SADEGsUpFC,SADEGsDwFC)
saveTablesTsvExc(SADEGs,outdir)

# 2- FAX SA ORA 
SAFAXprotQSig <- wholeDF$ProteomeFAX.qValue_FAXwithSA < protsPvalCut
SAFAXrbpmQsig <- wholeDF$RBPomeFAX.qValue_PolyARNAFAXwithSA < protsPvalCut
SAFAXprotUpFC <- wholeDF$ProteomeFAX.log2ratio_FAXwithSA > fcCut
SAFAXprotDwFC <- wholeDF$ProteomeFAX.log2ratio_FAXwithSA < -fcCut
SAFAXrbpmUpFC <- wholeDF$RBPomeFAX.log2ratio_PolyARNAFAXwithSA > fcCut
SAFAXrbpmDwFC <- wholeDF$RBPomeFAX.log2ratio_PolyARNAFAXwithSA < -fcCut
SAFAXupNC <- wholeDF$FAXnetchangesSA > fcCut
SAFAXdwNC <- wholeDF$FAXnetchangesSA < -fcCut

SAFAXprotUpFC <- getIDsList(wholeDF, which(SAFAXprotQSig & SAFAXprotUpFC), ProtIDs)
SAFAXprotDwFC <- getIDsList(wholeDF, which(SAFAXprotQSig & SAFAXprotDwFC), ProtIDs)
SAFAXrbpmUpFC <- getIDsList(wholeDF, which(SAFAXrbpmQsig & SAFAXrbpmUpFC), ProtIDs)
SAFAXrbpmDwFC <- getIDsList(wholeDF, which(SAFAXrbpmQsig & SAFAXrbpmDwFC), ProtIDs)
SAFAXnetcupNC <- getIDsList(wholeDF, which(SAFAXrbpmQsig & SAFAXupNC), ProtIDs)
SAFAXnetcdwNC <- getIDsList(wholeDF, which(SAFAXrbpmQsig & SAFAXdwNC), ProtIDs)

outdir <- "/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/FinalData/GeneLists/ORAseparatedUpnDw"
saveTablesTsvExc(SAFAXprotUpFC,outdir)
saveTablesTsvExc(SAFAXprotDwFC,outdir)
saveTablesTsvExc(SAFAXrbpmUpFC,outdir)
saveTablesTsvExc(SAFAXrbpmDwFC,outdir)
saveTablesTsvExc(SAFAXnetcupNC,outdir)
saveTablesTsvExc(SAFAXnetcdwNC,outdir)

outdir <- "/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/FinalData/GeneLists/ORA"
SAFAXprot <- c(SAFAXprotUpFC,SAFAXprotDwFC)
SAFAXrbpmFC <- c(SAFAXrbpmUpFC,SAFAXrbpmDwFC)
SAFAXnetc <- c(SAFAXnetcupNC,SAFAXnetcdwNC)
saveTablesTsvExc(SAFAXprot,outdir)
saveTablesTsvExc(SAFAXrbpmFC,outdir)
saveTablesTsvExc(SAFAXnetc,outdir)

# 2- UVX SA ORA
SAUVXprotQSig <- wholeDF$ProteomeUVX.qValue_Condition4 < protsPvalCut
SAUVXrbpmQsig <- wholeDF$RBPomeUVX.qValue_PolyARNAUVwithSA < protsPvalCut
SAUVXprotUpFC <- wholeDF$ProteomeUVX.log2ratio_Condition4 > fcCut
SAUVXprotDwFC <- wholeDF$ProteomeUVX.log2ratio_Condition4 < -fcCut
SAUVXrbpmUpFC <- wholeDF$RBPomeUVX.log2ratio_PolyARNAUVwithSA > fcCut
SAUVXrbpmDwFC <- wholeDF$RBPomeUVX.log2ratio_PolyARNAUVwithSA < -fcCut
SAUVXupNC <- wholeDF$UVXnetchangesSA > fcCut
SAUVXdwNC <- wholeDF$UVXnetchangesSA < -fcCut

SAUVXprotUpFC <- getIDsList(wholeDF,which(SAUVXprotQSig & SAUVXprotUpFC),ProtIDs)
SAUVXprotDwFC <- getIDsList(wholeDF,which(SAUVXprotQSig & SAUVXprotDwFC),ProtIDs)
SAUVXrbpmUpFC <- getIDsList(wholeDF,which(SAUVXrbpmQsig & SAUVXrbpmUpFC),ProtIDs)
SAUVXrbpmDwFC <- getIDsList(wholeDF,which(SAUVXrbpmQsig & SAUVXrbpmDwFC),ProtIDs)
SAUVXnetcupNC <- getIDsList(wholeDF,which(SAUVXrbpmQsig & SAUVXupNC),ProtIDs)
SAUVXnetcdwNC <- getIDsList(wholeDF,which(SAUVXrbpmQsig & SAUVXdwNC),ProtIDs)

outdir <- "/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/FinalData/GeneLists/ORAseparatedUpnDw"
saveTablesTsvExc(SAUVXprotUpFC,outdir)
saveTablesTsvExc(SAUVXprotDwFC,outdir)
saveTablesTsvExc(SAUVXrbpmUpFC,outdir)
saveTablesTsvExc(SAUVXrbpmDwFC,outdir)
saveTablesTsvExc(SAUVXnetcupNC,outdir)
saveTablesTsvExc(SAUVXnetcdwNC,outdir)

outdir <- "/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/FinalData/GeneLists/ORA"
SAUVXprot <- c(SAUVXprotUpFC,SAUVXprotDwFC)
SAUVXrbpm <- c(SAUVXrbpmUpFC,SAUVXrbpmDwFC)
SAUVXnetc <- c(SAUVXnetcupNC,SAUVXnetcdwNC)

saveTablesTsvExc(SAUVXprot,outdir)
saveTablesTsvExc(SAUVXrbpm,outdir)
saveTablesTsvExc(SAUVXnetc,outdir)

# GSEA 

getGSEArankDF <- function(DF,dataype,crosslink,treatment,IDs,rankCol="log2"){
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
SAdegsGSEA <- getGSEArankDF(wholeDF,dataype,crosslink,treatment,DEGsIDs)
saveTablesTsvExc(SAdegsGSEA,outdir,completeNdedup = T,excel = F,bycompleteFC = F,rownames = F)

dataype <- "Proteome"; crosslink <- "FAX"; treatment <- "SA"
SAProteomeFAXgsea <- getGSEArankDF(wholeDF,dataype,crosslink,treatment,ProtIDs)
saveTablesTsvExc(SAProteomeFAXgsea,outdir,completeNdedup = F,excel = F,bycompleteFC = F,rownames = F)

dataype <- "RBPome"; crosslink <- "FAX"; treatment <- "SA"
SARBPomeFAXgsea <- getGSEArankDF(wholeDF,dataype,crosslink,treatment,ProtIDs)
saveTablesTsvExc(SARBPomeFAXgsea,outdir,completeNdedup = F,excel = F,bycompleteFC = F,rownames = F)

SAnetchangesFAXgsea <- getGSEArankDF(wholeDF,dataype,crosslink,treatment,ProtIDs,"netchanges")
saveTablesTsvExc(SAnetchangesFAXgsea,outdir,completeNdedup = F,excel = F,bycompleteFC = F,rownames = F)

dataype <- "Proteome"; crosslink <- "UVX"; treatment <- "SA"
SAProteomeUVXgsea <- getGSEArankDF(wholeDF,dataype,crosslink,treatment,ProtIDs)
saveTablesTsvExc(SAProteomeUVXgsea,outdir,completeNdedup = F,excel = F,bycompleteFC = F,rownames = F)

dataype <- "RBPome"; crosslink <- "UVX"; treatment <- "SA"
SARBPomeUVXgsea <- getGSEArankDF(wholeDF,dataype,crosslink,treatment,ProtIDs)
saveTablesTsvExc(SARBPomeUVXgsea,outdir,completeNdedup = F,excel = F,bycompleteFC = F,rownames = F)

SAnetchangesUVXgsea <- getGSEArankDF(wholeDF,dataype,crosslink,treatment,ProtIDs,"netchanges")
saveTablesTsvExc(SAnetchangesUVXgsea,outdir,completeNdedup = F,excel = F,bycompleteFC = F,rownames = F)


