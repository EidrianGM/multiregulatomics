wholeDFFile <- "FinalData/allDataDF.tsv"
wholeDF <- read.delim(wholeDFFile,quote = "")

FAXacceptedProts <- read.delim("FinalData/BackgroundRemoval/FAXAccBackGrAcceptedFC3.tsv",header = F)[,1]
UVXacceptedProts <- read.delim("FinalData/BackgroundRemoval/UVXAccBackGrAcceptedFC3.tsv",header = F)[,1]

################################################################################
######### Generation of GeneLists for Enrichment Analyses ORA and GSEA #########
################################################################################
outdir <- "FinalData/GeneLists"

DEGsIDs <- c("UniprotACC") # Is the one with most significant annotated in GC4
ProtIDs <- c("proteinName") # To also consider orthologues that will be splited for ORA

DEGpvalCut <- 0.01
DEGFCCut <- 3
protsPvalCut <- 0.05
fcCut   <- 1

# ORA
getIDsList <- function(DF,logicalRows,ID, keep=c()){
  DF <- DF[logicalRows,ID]
  if (class(DF) == "character"){
    DF <- unique(DF) 
    print(paste(length(DF),"are significant"))
    if (length(keep) != 0){
      cat("Removed by back ",sum(!DF %in% keep))
      DF <- DF[DF %in% keep]
      print(paste(length(DF),"are remaining"))
    }
  }else{
    DF <- DF[!duplicated(DF),]
    print(paste(nrow(DF),"are significant")) 
    if (length(keep) > 0){
      cat("Removed by back ",sum(!DF[,ID] %in% keep))
      DF <- DF[DF[,ID] %in% keep,]
      print(paste(nrow(DF)," are remaining")) 
    }
  }
  return(DF)
} 

# 2- DEGS SA ORA 
SAdegsQSig <- wholeDF$DEGs.adj.pval.SA.YPD < DEGpvalCut
SAdegsUpFC <- wholeDF$DEGs.log2FC.SA.YPD > DEGFCCut
SAdegsDwFC <- wholeDF$DEGs.log2FC.SA.YPD < -DEGFCCut
SADEGs <- getIDsList(wholeDF,which(SAdegsQSig),DEGsIDs) # DEGsIDs <- 'DEGs.target'
SADEGsUpFC <- getIDsList(wholeDF,which(SAdegsQSig & SAdegsUpFC),DEGsIDs)
SADEGsDwFC <- getIDsList(wholeDF,which(SAdegsQSig & SAdegsDwFC),DEGsIDs)

outdir <- "FinalData/GeneLists/ORAseparatedUpnDw"
saveTablesTsvExc(SADEGsUpFC,outdir)
saveTablesTsvExc(SADEGsDwFC,outdir)

outdir <- "FinalData/GeneLists/ORA"
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

SAfaxDEGsProt <- getIDsList(wholeDF,which(SAFAXprotQSig),ProtIDs)
SAfaxDEGsmRBP <- getIDsList(wholeDF,which(SAFAXrbpmQsig),ProtIDs)

SAFAXprotUpFC <- getIDsList(wholeDF, which(SAFAXprotQSig & SAFAXprotUpFC), ProtIDs)
SAFAXprotDwFC <- getIDsList(wholeDF, which(SAFAXprotQSig & SAFAXprotDwFC), ProtIDs)
SAFAXrbpmUpFC <- getIDsList(wholeDF, which(SAFAXrbpmQsig & SAFAXrbpmUpFC), ProtIDs, FAXacceptedProts)
SAFAXrbpmDwFC <- getIDsList(wholeDF, which(SAFAXrbpmQsig & SAFAXrbpmDwFC), ProtIDs, FAXacceptedProts)
SAFAXnetcupNC <- getIDsList(wholeDF, which(SAFAXrbpmQsig & SAFAXupNC), ProtIDs, FAXacceptedProts)
SAFAXnetcdwNC <- getIDsList(wholeDF, which(SAFAXrbpmQsig & SAFAXdwNC), ProtIDs, FAXacceptedProts)

outdir <- "FinalData/GeneLists/ORAseparatedUpnDw"
saveTablesTsvExc(SAFAXprotUpFC,outdir)
saveTablesTsvExc(SAFAXprotDwFC,outdir)
saveTablesTsvExc(SAFAXrbpmUpFC,outdir)
saveTablesTsvExc(SAFAXrbpmDwFC,outdir)
saveTablesTsvExc(SAFAXnetcupNC,outdir)
saveTablesTsvExc(SAFAXnetcdwNC,outdir)

outdir <- "FinalData/GeneLists/ORA"
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

SAuvxDEGsProt <- getIDsList(wholeDF,which(SAUVXprotQSig),ProtIDs)
SAuvxDEGsmRBP <- getIDsList(wholeDF,which(SAUVXrbpmQsig),ProtIDs)

SAUVXprotUpFC <- getIDsList(wholeDF,which(SAUVXprotQSig & SAUVXprotUpFC),ProtIDs)
SAUVXprotDwFC <- getIDsList(wholeDF,which(SAUVXprotQSig & SAUVXprotDwFC),ProtIDs)
SAUVXrbpmUpFC <- getIDsList(wholeDF,which(SAUVXrbpmQsig & SAUVXrbpmUpFC),ProtIDs, UVXacceptedProts) # 1 change
SAUVXrbpmDwFC <- getIDsList(wholeDF,which(SAUVXrbpmQsig & SAUVXrbpmDwFC),ProtIDs, UVXacceptedProts)
SAUVXnetcupNC <- getIDsList(wholeDF,which(SAUVXrbpmQsig & SAUVXupNC),ProtIDs, UVXacceptedProts) # 1 change
SAUVXnetcdwNC <- getIDsList(wholeDF,which(SAUVXrbpmQsig & SAUVXdwNC),ProtIDs, UVXacceptedProts)

outdir <- "FinalData/GeneLists/ORAseparatedUpnDw"
saveTablesTsvExc(SAUVXprotUpFC,outdir)
saveTablesTsvExc(SAUVXprotDwFC,outdir)
saveTablesTsvExc(SAUVXrbpmUpFC,outdir)
saveTablesTsvExc(SAUVXrbpmDwFC,outdir)
saveTablesTsvExc(SAUVXnetcupNC,outdir)
saveTablesTsvExc(SAUVXnetcdwNC,outdir)

outdir <- "FinalData/GeneLists/ORA"
SAUVXprot <- c(SAUVXprotUpFC,SAUVXprotDwFC)
SAUVXrbpm <- c(SAUVXrbpmUpFC,SAUVXrbpmDwFC)
SAUVXnetc <- c(SAUVXnetcupNC,SAUVXnetcdwNC)

saveTablesTsvExc(SAUVXprot,outdir)
saveTablesTsvExc(SAUVXrbpm,outdir)
saveTablesTsvExc(SAUVXnetc,outdir)

# GSEA 

getGSEArankDF <- function(DF,dataype,crosslink,treatment,IDs,rankCol="log2", keep=c()){
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
  if (length(keep) > 0){
    cat("Removed by back ",sum(!DF[,IDs] %in% keep))
    DF <- DF[DF[,IDs] %in% keep,]
    print(paste(nrow(DF)," are remaining")) 
  }
  return(DF)
} 
outdir <- "FinalData/GeneLists/GSEA"
DF <- wholeDF

DEGsIDs <- 'UniprotACC'
dataype <- "DEGs"; crosslink <- ""; treatment <- "SA"
SAdegsGSEA <- getGSEArankDF(wholeDF,dataype,crosslink,treatment,DEGsIDs)
saveTablesTsvExc(SAdegsGSEA,outdir,completeNdedup = T,excel = F,bycompleteFC = F,rownames = F)

dataype <- "Proteome"; crosslink <- "FAX"; treatment <- "SA"
SAProteomeFAXgsea <- getGSEArankDF(wholeDF,dataype,crosslink,treatment,ProtIDs)
saveTablesTsvExc(SAProteomeFAXgsea,outdir,completeNdedup = F,excel = F,bycompleteFC = F,rownames = F)

dataype <- "RBPome"; crosslink <- "FAX"; treatment <- "SA"
SARBPomeFAXgsea <- getGSEArankDF(wholeDF,dataype,crosslink,treatment,ProtIDs,keep = FAXacceptedProts) # 2 change
saveTablesTsvExc(SARBPomeFAXgsea,outdir,completeNdedup = F,excel = F,bycompleteFC = F,rownames = F)

SAnetchangesFAXgsea <- getGSEArankDF(wholeDF,dataype,crosslink,treatment,ProtIDs,"netchanges",FAXacceptedProts) # 2 change
saveTablesTsvExc(SAnetchangesFAXgsea,outdir,completeNdedup = F,excel = F,bycompleteFC = F,rownames = F)

dataype <- "Proteome"; crosslink <- "UVX"; treatment <- "SA"
SAProteomeUVXgsea <- getGSEArankDF(wholeDF,dataype,crosslink,treatment,ProtIDs)
saveTablesTsvExc(SAProteomeUVXgsea,outdir,completeNdedup = F,excel = F,bycompleteFC = F,rownames = F)

dataype <- "RBPome"; crosslink <- "UVX"; treatment <- "SA"
SARBPomeUVXgsea <- getGSEArankDF(wholeDF,dataype,crosslink,treatment,ProtIDs,keep = UVXacceptedProts) # 9 change
saveTablesTsvExc(SARBPomeUVXgsea,outdir,completeNdedup = F,excel = F,bycompleteFC = F,rownames = F)

SAnetchangesUVXgsea <- getGSEArankDF(wholeDF,dataype,crosslink,treatment,ProtIDs,"netchanges",UVXacceptedProts) # 9 change
saveTablesTsvExc(SAnetchangesUVXgsea,outdir,completeNdedup = F,excel = F,bycompleteFC = F,rownames = F)


outdir <- "FinalData/GeneLists/ORA"
files <- list.files(outdir,pattern = "*.txt",full.names = T)
outdir <- "FinalData/GeneLists/GSEA"
files <- c(files, list.files(outdir,pattern = "*.tsv",full.names = T))

stats <- c()
upreg <- "-"; downreg <- "-"
for (file in files){
  enrchTec <- gsub(".*/","",dirname(file)) 
  dataype <- gsub("\\..*","",basename(file))
  if (grepl('degs',dataype,ignore.case = T)) {
    fcCut <- 3; pvalCut <- 0.01
  }else{
    fcCut <- 1; pvalCut <- 0.05
  }
  if (enrchTec == "GSEA"){
    pvalCut <- "-"
    degs <- read.delim(file)  
    northologspep <- sum(grepl(';',degs[,1]))
    singleprots <- sum(!grepl(';',degs[,1]))
    upreg <- degs[which(degs[,2] > fcCut),1]
    downreg <- degs[which(degs[,2] < -fcCut),1]
    upregorthologspep <- sum(grepl(";",upreg))
    downregorthologspep <- sum(grepl(";",downreg))
    upreg <- length(upreg) - upregorthologspep
    downreg <- length(downreg) - downregorthologspep
    #northologspep <- upregorthologspep + downregorthologspep
    sigleuniprots <- nrow(degs) - northologspep
    singleprots == sigleuniprots
  }else{
    degs <- read.delim(file,header = F)  
    northologspep <- sum(grepl(";",degs[,1]))
    sigleuniprots <- nrow(degs) - northologspep
    upregorthologspep <- downregorthologspep <- '-'
  }
  total <- nrow(degs)
  stats <- rbind(stats,cbind(enrchTec,dataype,total,sigleuniprots,northologspep,
                             fcCut,pvalCut,upreg,downreg,
                             upregorthologspep,downregorthologspep))
}
stats <- as.data.frame(stats)
outdir <- "FinalData/GeneLists"
saveTablesTsvExc(stats,outdir,completeNdedup = F,excel = T,bycompleteFC = F,rownames = F)

View(stats)

wholeDFFile <- "FinalData/allDataDF.tsv"
wholeDF <- read.delim(wholeDFFile,quote = "")

sum(na.omit(wholeDF$DEGs.adj.pval.SA.YPD) < 0.01)
sum(na.omit(wholeDF$ProteomeFAX.qValue_FAXwithSA) < 0.05)
sum(na.omit(wholeDF$ProteomeUVX.qValue_Condition4) < 0.05)
sum(na.omit(wholeDF$RBPomeFAX.qValue_PolyARNAFAXwithSA) < 0.05)
sum(na.omit(wholeDF$RBPomeUVX.qValue_PolyARNAUVwithSA) < 0.05)

pvalsignificants <- c(sum(na.omit(wholeDF$DEGs.adj.pval.SA.YPD) < 0.01),
                      sum(na.omit(wholeDF$ProteomeFAX.qValue_FAXwithSA) < 0.05),
                      sum(na.omit(wholeDF$ProteomeUVX.qValue_Condition4) < 0.05),
                      sum(na.omit(wholeDF$RBPomeFAX.qValue_PolyARNAFAXwithSA) < 0.05),
                      sum(na.omit(wholeDF$RBPomeUVX.qValue_PolyARNAUVwithSA) < 0.05))

