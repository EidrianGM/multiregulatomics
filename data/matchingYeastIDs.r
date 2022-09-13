library(openxlsx)
library(rtracklayer)
library(plyr)
library(reshape2)

completeNdedup <- function(DF){
  if (class(DF) == "data.frame"){
    return(DF[complete.cases(DF),])
  }else{
    return(unique(na.omit(DF)))
  }
}

saveTablesTsvExc <- function(DF,outdir){
  outname <- deparse(substitute(DF))
  if (class(DF) == "data.frame"){
    if (nrow(DF) == 0){
      print("Empty DF")
    }else{
      outFile <- paste0(file.path(outdir,outname),".tsv")
      DF <- DF[!duplicated(DF),]
      write.table(DF,outFile,quote = F,sep = '\t',col.names = T,row.names = F)
      write.xlsx(DF,gsub(".tsv",".xlsx",outFile), asTable = FALSE, overwrite = TRUE)
    }
  } else{
    if (length(DF) == 0){
      print("Empty Vector")
    }else{
      outFile <- paste0(file.path(outdir,outname),".txt")
      write(DF,outFile,sep = "\n")
    }
  }
}

##################################
##### DATA and MAPPING FILES #####
##################################

wholeDFfile <- "/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/data/mapped/allDataMapped_together.tsv"
wholeDF <- read.delim(wholeDFfile,quote = "")
mappingFinalFile <- "/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/data/yeastReference/mappingFile.tsv"
yeastGenesProtMap <- read.delim(mappingFinalFile,quote = "")

#######################################
##### Protein Net Activities Data #####
#######################################

ibtissamDir <- "/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/data/mapped/forIbtissam"

# 1- RNABProteins FAX

netChangesFAXDTT <- wholeDF$mRNABPomeFAX.log2ratio_PolyARNAFAXwithDTT - wholeDF$proteomeFAX.log2ratio_Condition2
wholeDF$netChangesFAXDTT <- netChangesFAXDTT
netChangesFAXH2O2 <- wholeDF$mRNABPomeFAX.log2ratio_PolyARNAFAXwithH2O2 - wholeDF$proteomeFAX.log2ratio_Condition3
wholeDF$netChangesFAXH2O2 <- netChangesFAXH2O2
netChangesFAXSA <- wholeDF$mRNABPomeFAX.log2ratio_PolyARNAFAXwithSA - wholeDF$proteomeFAX.log2ratio_Condition4
wholeDF$netChangesFAXSA <- netChangesFAXSA

FAXcolsNetChanges <- c("mRNABPomeFAX.pValue_PolyARNAFAXwithDTT","proteomeFAX.pValue_Condition2",
                       "mRNABPomeFAX.qValue_PolyARNAFAXwithDTT","proteomeFAX.qValue_Condition2",
                       "mRNABPomeFAX.log2ratio_PolyARNAFAXwithDTT","proteomeFAX.log2ratio_Condition2", "netChangesFAXDTT",
                       
                       "mRNABPomeFAX.pValue_PolyARNAFAXwithH2O2","proteomeFAX.pValue_Condition3",
                       "mRNABPomeFAX.qValue_PolyARNAFAXwithH2O2","proteomeFAX.qValue_Condition3",
                       "mRNABPomeFAX.log2ratio_PolyARNAFAXwithH2O2","proteomeFAX.log2ratio_Condition3", "netChangesFAXH2O2",
                       
                       "mRNABPomeFAX.pValue_PolyARNAFAXwithSA","proteomeFAX.pValue_Condition4",
                       "mRNABPomeFAX.qValue_PolyARNAFAXwithSA","proteomeFAX.qValue_Condition4",
                       "mRNABPomeFAX.log2ratio_PolyARNAFAXwithSA","proteomeFAX.log2ratio_Condition4", "netChangesFAXSA")

FAXprotNrbpome <- wholeDF[(!is.na(netChangesFAXDTT) | !is.na(netChangesFAXH2O2) | !is.na(netChangesFAXSA)), c(colnames(yeastGenesProtMap),FAXcolsNetChanges)]
saveTablesTsvExc(FAXprotNrbpome,ibtissamDir)

# 2- RNABProteins UV

netChangesUVDTT <- wholeDF$mRNABPomeUV.log2ratio_PolyARNAUVwithDTT - wholeDF$proteomeUV.log2ratio_Condition2
wholeDF$netChangesUVDTT <- netChangesUVDTT
netChangesUVH2O2 <- wholeDF$mRNABPomeUV.log2ratio_PolyARNAUVwithH202 - wholeDF$proteomeUV.log2ratio_Condition3
wholeDF$netChangesUVH2O2 <- netChangesUVH2O2
netChangesUVSA <- wholeDF$mRNABPomeUV.log2ratio_PolyARNAUVwithSA - wholeDF$proteomeUV.log2ratio_Condition4
wholeDF$netChangesUVSA <- netChangesUVSA

UVprotNrbpomeCols <- c("mRNABPomeUV.pValue_PolyARNAUVwithDTT","proteomeUV.pValue_Condition2",
                       "mRNABPomeUV.qValue_PolyARNAUVwithDTT","proteomeUV.qValue_Condition2",
                       "mRNABPomeUV.log2ratio_PolyARNAUVwithDTT","proteomeUV.log2ratio_Condition2", "netChangesUVDTT",
                       
                       "mRNABPomeUV.pValue_PolyARNAUVwithH202","proteomeUV.pValue_Condition3",
                       "mRNABPomeUV.qValue_PolyARNAUVwithH202","proteomeUV.qValue_Condition3",
                       "mRNABPomeUV.log2ratio_PolyARNAUVwithH202","proteomeUV.log2ratio_Condition3", "netChangesUVH2O2",
                       
                       "mRNABPomeUV.pValue_PolyARNAUVwithSA","proteomeUV.pValue_Condition4",
                       "mRNABPomeUV.qValue_PolyARNAUVwithSA","proteomeUV.qValue_Condition4",
                       "mRNABPomeUV.log2ratio_PolyARNAUVwithSA","proteomeUV.log2ratio_Condition4", "netChangesUVSA")

UVprotNrbpome <- wholeDF[(!is.na(netChangesUVDTT) | !is.na(netChangesUVH2O2) | !is.na(netChangesUVSA)), c(colnames(yeastGenesProtMap),UVprotNrbpomeCols)]
saveTablesTsvExc(UVprotNrbpome,ibtissamDir)

FAXnUVprotNrbpome <- merge(UVprotNrbpome,FAXprotNrbpome)
saveTablesTsvExc(FAXnUVprotNrbpome,ibtissamDir)

##########################
##### FILTERING DATA #####
##########################
toPO4outdir <- "/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/data/mapped/toPaintomics"

# 1- Differentially Expressed Genes
DEGpvalCut <- 0.01
DEGfcCut <- 1

DEGsH202 <- wholeDF$GeneName[wholeDF$DEGs.adj.pval.H202.YPD < DEGpvalCut & abs(wholeDF$DEGs.log2FC.H202.YPD) > DEGfcCut]
DEGsDTT <- wholeDF$GeneName[wholeDF$DEGs.adj.pval.DTT.YPD < DEGpvalCut  & abs(wholeDF$DEGs.log2FC.DTT.YPD) > DEGfcCut]
DEGsSA <- wholeDF$GeneName[wholeDF$DEGs.adj.pval.SA.YPD < DEGpvalCut & abs(wholeDF$DEGs.log2FC.SA.YPD) > DEGfcCut]
DEGsH202 <- unique(na.omit(DEGsH202))
DEGsDTT <- unique(na.omit(DEGsDTT))
DEGsSA <- unique(na.omit(DEGsSA))
saveTablesTsvExc(DEGsH202,toPO4outdir)
saveTablesTsvExc(DEGsDTT,toPO4outdir)
saveTablesTsvExc(DEGsSA,toPO4outdir)

DEGsH202po4 <- wholeDF[order(wholeDF$DEGs.log2FC.H202.YPD,decreasing = T),c("GeneName","DEGs.log2FC.H202.YPD")]
DEGsDTTpo4 <- wholeDF[order(wholeDF$DEGs.log2FC.DTT.YPD,decreasing = T),c("GeneName","DEGs.log2FC.DTT.YPD")]
DEGsSApo4 <- wholeDF[order(wholeDF$DEGs.log2FC.SA.YPD,decreasing = T),c("GeneName","DEGs.log2FC.SA.YPD")]
DEGsH202po4 <- DEGsH202po4[complete.cases(DEGsH202po4),]
DEGsDTTpo4 <- DEGsDTTpo4[complete.cases(DEGsDTTpo4),]
DEGsSApo4 <- DEGsSApo4[complete.cases(DEGsSApo4),]

saveTablesTsvExc(DEGsDTTpo4,toPO4outdir)
saveTablesTsvExc(DEGsDTTpo4,toPO4outdir)
saveTablesTsvExc(DEGsSApo4,toPO4outdir)

# 2- Proteins and mRNABProteins
protsPvalCut <- 0.05
protsFcCut <- 1
netChangesCut <- 1

# 2.1- FAX
proteomeFAXfcDTT <- wholeDF[,c("UniprotACC","proteomeFAX.log2ratio_Condition2")]
proteomeFAXfcH2O2 <- wholeDF[,c("UniprotACC","proteomeFAX.log2ratio_Condition3")]
proteomeFAXfcSA <- wholeDF[,c("UniprotACC","proteomeFAX.log2ratio_Condition4")]
mRNABPomeFAXfcDTT <- wholeDF[,c("UniprotACC","mRNABPomeFAX.log2ratio_PolyARNAFAXwithDTT")]
mRNABPomeFAXfcH2O2 <- wholeDF[,c("UniprotACC","mRNABPomeFAX.log2ratio_PolyARNAFAXwithH2O2")]
mRNABPomeFAXfcSA <- wholeDF[,c("UniprotACC","mRNABPomeFAX.log2ratio_PolyARNAFAXwithSA")]
netChangesFAXDTT <- wholeDF[,c("UniprotACC","netChangesFAXDTT")]
netChangesFAXH2O2 <- wholeDF[,c("UniprotACC","netChangesFAXH2O2")]
netChangesFAXSA <- wholeDF[,c("UniprotACC","netChangesFAXSA")]

proteomeFAXfcDTT <- proteomeFAXfcDTT[complete.cases(proteomeFAXfcDTT),]
proteomeFAXfcH2O2 <- proteomeFAXfcH2O2[complete.cases(proteomeFAXfcH2O2),]
proteomeFAXfcSA <- proteomeFAXfcSA[complete.cases(proteomeFAXfcSA),]
mRNABPomeFAXfcDTT <- mRNABPomeFAXfcDTT[complete.cases(mRNABPomeFAXfcDTT),]
mRNABPomeFAXfcH2O2 <- mRNABPomeFAXfcH2O2[complete.cases(mRNABPomeFAXfcH2O2),]
mRNABPomeFAXfcSA <- mRNABPomeFAXfcSA[complete.cases(mRNABPomeFAXfcSA),]
netChangesFAXDTT <- netChangesFAXDTT[complete.cases(netChangesFAXDTT),]
netChangesFAXH2O2 <- netChangesUVH2O2[complete.cases(netChangesFAXH2O2),]
netChangesFAXSA <- netChangesUVSA[complete.cases(netChangesFAXSA),]

proteomeFAXfcDTT <- proteomeFAXfcDTT[order(proteomeFAXfcDTT[,2],decreasing = T),]
proteomeFAXfcH2O2 <- proteomeFAXfcH2O2[order(proteomeFAXfcH2O2[,2],decreasing = T),]
proteomeFAXfcSA <- proteomeFAXfcSA[order(proteomeFAXfcSA[,2],decreasing = T),]
mRNABPomeFAXfcDTT <- mRNABPomeFAXfcDTT[order(mRNABPomeFAXfcDTT[,2],decreasing = T),]
mRNABPomeFAXfcH2O2 <- mRNABPomeFAXfcH2O2[order(mRNABPomeFAXfcH2O2[,2],decreasing = T),]
mRNABPomeFAXfcSA <- mRNABPomeFAXfcSA[order(mRNABPomeFAXfcSA[,2],decreasing = T),]
netChangesFAXDTT <- netChangesFAXDTT[order(netChangesFAXDTT[,2],decreasing = T),]
netChangesFAXH2O2 <- netChangesUVH2O2[order(netChangesUVH2O2[,2],decreasing = T),]
netChangesFAXSA <- netChangesUVSA[order(netChangesUVSA[,2],decreasing = T),]

saveTablesTsvExc(proteomeFAXfcDTT,toPO4outdir)
saveTablesTsvExc(proteomeFAXfcH2O2,toPO4outdir)
saveTablesTsvExc(proteomeFAXfcSA,toPO4outdir)
saveTablesTsvExc(mRNABPomeFAXfcDTT,toPO4outdir)
saveTablesTsvExc(mRNABPomeFAXfcH2O2,toPO4outdir)
saveTablesTsvExc(mRNABPomeFAXfcSA,toPO4outdir)
saveTablesTsvExc(netChangesFAXDTT,toPO4outdir)
saveTablesTsvExc(netChangesFAXH2O2,toPO4outdir)
saveTablesTsvExc(netChangesFAXSA,toPO4outdir)

## 2.1.1 - FAX DTT Selection
bypValueProteomeFAXDTT <- wholeDF$proteomeFAX.pValue_Condition2 < protsPvalCut
byqValueProteomeFAXDTT <- wholeDF$proteomeFAX.qValue_Condition2 < protsPvalCut
byFCproteomeFAXDTT <- abs(wholeDF$proteomeFAX.log2ratio_Condition2) > protsFcCut
bypValuemRNABPomeFAXDTT <- wholeDF$mRNABPomeFAX.pValue_PolyARNAFAXwithDTT < protsPvalCut
byqValuemRNABPomeFAXDTT <- wholeDF$mRNABPomeFAX.qValue_PolyARNAFAXwithDTT < protsPvalCut
byFCmRNABPomeFAXDTT <- abs(wholeDF$mRNABPomeFAX.log2ratio_PolyARNAFAXwithDTT) > protsFcCut
byNetChangesFAXDTT <- abs(wholeDF$netChangesFAXDTT) > protsFcCut

pValFCproteomeFAXDTTselected <- na.omit(unique(wholeDF$UniprotACC[bypValueProteomeFAXDTT & byFCproteomeFAXDTT])); length(pValFCproteomeFAXDTTselected)
qValFCproteomeFAXDTTselected <- na.omit(unique(wholeDF$UniprotACC[byqValueProteomeFAXDTT & byFCproteomeFAXDTT])); length(qValFCproteomeFAXDTTselected)
pValFCmRNABPomeFAXDTTselected <- na.omit(unique(wholeDF$UniprotACC[bypValuemRNABPomeFAXDTT & byFCmRNABPomeFAXDTT])); length(pValFCmRNABPomeFAXDTTselected)
qValFCmRNABPomeFAXDTTselected <- na.omit(unique(wholeDF$UniprotACC[byqValuemRNABPomeFAXDTT & byFCmRNABPomeFAXDTT])); length(qValFCmRNABPomeFAXDTTselected)
saveTablesTsvExc(pValFCproteomeFAXDTTselected,toPO4outdir)
saveTablesTsvExc(qValFCproteomeFAXDTTselected,toPO4outdir)
saveTablesTsvExc(pValFCmRNABPomeFAXDTTselected,toPO4outdir)
saveTablesTsvExc(qValFCmRNABPomeFAXDTTselected,toPO4outdir)

netChangesFAXDTTprotNrbpomeSelPval <- na.omit(unique(wholeDF$UniprotACC[byNetChangesFAXDTT & bypValueProteomeFAXDTT & bypValuemRNABPomeFAXDTT])) # 0
netChangesFAXDTTprotNrbpomeSelQval <- na.omit(unique(wholeDF$UniprotACC[byNetChangesFAXDTT & byqValueProteomeFAXDTT & byqValuemRNABPomeFAXDTT])) # 0
netChangesFAXDTTSelFC <- na.omit(unique(wholeDF$UniprotACC[byNetChangesFAXDTT]))

saveTablesTsvExc(netChangesFAXDTTprotNrbpomeSelPval,toPO4outdir)
saveTablesTsvExc(netChangesFAXDTTprotNrbpomeSelQval,toPO4outdir)
saveTablesTsvExc(netChangesFAXDTTSelFC,toPO4outdir)

## 2.1.2 - FAX H2O2 Selection
bypValueProteomeFAXH2O2 <- wholeDF$proteomeFAX.pValue_Condition2 < protsPvalCut
byqValueProteomeFAXH2O2 <- wholeDF$proteomeFAX.qValue_Condition2 < protsPvalCut
byFCproteomeFAXH2O2 <- abs(wholeDF$proteomeFAX.log2ratio_Condition2) > protsFcCut
bypValuemRNABPomeFAXH2O2 <- wholeDF$mRNABPomeFAX.pValue_PolyARNAFAXwithH2O2 < protsPvalCut
byqValuemRNABPomeFAXH2O2 <- wholeDF$mRNABPomeFAX.qValue_PolyARNAFAXwithH2O2 < protsPvalCut
byFCmRNABPomeFAXH2O2 <- abs(wholeDF$mRNABPomeFAX.log2ratio_PolyARNAFAXwithH2O2) > protsFcCut
byNetChangesFAXH2O2 <- abs(wholeDF$netChangesFAXH2O2) > protsFcCut

pValFCproteomeFAXH2O2selected <- na.omit(unique(wholeDF$UniprotACC[bypValueProteomeFAXH2O2 & byFCproteomeFAXH2O2])); length(pValFCproteomeFAXH2O2selected)
qValFCproteomeFAXH2O2selected <- na.omit(unique(wholeDF$UniprotACC[byqValueProteomeFAXH2O2 & byFCproteomeFAXH2O2])); length(qValFCproteomeFAXH2O2selected)
pValFCmRNABPomeFAXH2O2selected <- na.omit(unique(wholeDF$UniprotACC[bypValuemRNABPomeFAXH2O2 & byFCmRNABPomeFAXH2O2])); length(pValFCmRNABPomeFAXH2O2selected)
qValFCmRNABPomeFAXH2O2selected <- na.omit(unique(wholeDF$UniprotACC[byqValuemRNABPomeFAXH2O2 & byFCmRNABPomeFAXH2O2])); length(qValFCmRNABPomeFAXH2O2selected)
saveTablesTsvExc(pValFCproteomeFAXH2O2selected,toPO4outdir)
saveTablesTsvExc(qValFCproteomeFAXH2O2selected,toPO4outdir)
saveTablesTsvExc(pValFCmRNABPomeFAXH2O2selected,toPO4outdir)
saveTablesTsvExc(qValFCmRNABPomeFAXH2O2selected,toPO4outdir)

netChangesFAXH2O2protNrbpomeSelPval <- na.omit(unique(wholeDF$UniprotACC[byNetChangesFAXH2O2 & bypValueProteomeFAXH2O2 & bypValuemRNABPomeFAXH2O2])) # 0
netChangesFAXH2O2protNrbpomeSelQval <- na.omit(unique(wholeDF$UniprotACC[byNetChangesFAXH2O2 & byqValueProteomeFAXH2O2 & byqValuemRNABPomeFAXH2O2])) # 0
netChangesFAXH2O2SelFC <- na.omit(unique(wholeDF$UniprotACC[byNetChangesFAXH2O2]))

saveTablesTsvExc(netChangesFAXH2O2protNrbpomeSelPval,toPO4outdir)
saveTablesTsvExc(netChangesFAXH2O2protNrbpomeSelQval,toPO4outdir)
saveTablesTsvExc(netChangesFAXH2O2SelFC,toPO4outdir)

## 2.1.3 - FAX SA Selection
bypValueProteomeFAXSA <- wholeDF$proteomeFAX.pValue_Condition2 < protsPvalCut
byqValueProteomeFAXSA <- wholeDF$proteomeFAX.qValue_Condition2 < protsPvalCut
byFCproteomeFAXSA <- abs(wholeDF$proteomeFAX.log2ratio_Condition2) > protsFcCut
bypValuemRNABPomeFAXSA <- wholeDF$mRNABPomeFAX.pValue_PolyARNAFAXwithSA < protsPvalCut
byqValuemRNABPomeFAXSA <- wholeDF$mRNABPomeFAX.qValue_PolyARNAFAXwithSA < protsPvalCut
byFCmRNABPomeFAXSA <- abs(wholeDF$mRNABPomeFAX.log2ratio_PolyARNAFAXwithSA) > protsFcCut
byNetChangesFAXSA <- abs(wholeDF$netChangesFAXSA) > protsFcCut

pValFCproteomeFAXSAselected <- na.omit(unique(wholeDF$UniprotACC[bypValueProteomeFAXSA & byFCproteomeFAXSA])); length(pValFCproteomeFAXSAselected)
qValFCproteomeFAXSAselected <- na.omit(unique(wholeDF$UniprotACC[byqValueProteomeFAXSA & byFCproteomeFAXSA])); length(qValFCproteomeFAXSAselected)
pValFCmRNABPomeFAXSAselected <- na.omit(unique(wholeDF$UniprotACC[bypValuemRNABPomeFAXSA & byFCmRNABPomeFAXSA])); length(pValFCmRNABPomeFAXSAselected)
qValFCmRNABPomeFAXSAselected <- na.omit(unique(wholeDF$UniprotACC[byqValuemRNABPomeFAXSA & byFCmRNABPomeFAXSA])); length(qValFCmRNABPomeFAXSAselected)
saveTablesTsvExc(pValFCproteomeFAXSAselected,toPO4outdir)
saveTablesTsvExc(qValFCproteomeFAXSAselected,toPO4outdir)
saveTablesTsvExc(pValFCmRNABPomeFAXSAselected,toPO4outdir)
saveTablesTsvExc(qValFCmRNABPomeFAXSAselected,toPO4outdir)

netChangesFAXSAprotNrbpomeSelPval <- na.omit(unique(wholeDF$UniprotACC[byNetChangesFAXSA & bypValueProteomeFAXSA & bypValuemRNABPomeFAXSA])) # 0
netChangesFAXSAprotNrbpomeSelQval <- na.omit(unique(wholeDF$UniprotACC[byNetChangesFAXSA & byqValueProteomeFAXSA & byqValuemRNABPomeFAXSA])) # 0
netChangesFAXSASelFC <- na.omit(unique(wholeDF$UniprotACC[byNetChangesFAXSA]))

saveTablesTsvExc(netChangesFAXSAprotNrbpomeSelPval,toPO4outdir)
saveTablesTsvExc(netChangesFAXSAprotNrbpomeSelQval,toPO4outdir)
saveTablesTsvExc(netChangesFAXSASelFC,toPO4outdir)

# 2.2 - UV
proteomeUVfcDTT <- wholeDF[,c("UniprotACC","proteomeUV.log2ratio_Condition2")]
proteomeUVfcH2O2 <- wholeDF[,c("UniprotACC","proteomeUV.log2ratio_Condition3")]
proteomeUVfcSA <- wholeDF[,c("UniprotACC","proteomeUV.log2ratio_Condition4")]
mRNABPomeUVfcDTT <- wholeDF[,c("UniprotACC","mRNABPomeUV.log2ratio_PolyARNAUVwithDTT")]
mRNABPomeUVfcH2O2 <- wholeDF[,c("UniprotACC","mRNABPomeUV.log2ratio_PolyARNAUVwithH202")]
mRNABPomeUVfcSA <- wholeDF[,c("UniprotACC","mRNABPomeUV.log2ratio_PolyARNAUVwithSA")]
netChangesUVDTT <- wholeDF[,c("UniprotACC","netChangesUVDTT")]
netChangesUVH2O2 <- wholeDF[,c("UniprotACC","netChangesUVH2O2")]
netChangesUVSA <- wholeDF[,c("UniprotACC","netChangesUVSA")]

proteomeUVfcDTT <- proteomeUVfcDTT[complete.cases(proteomeUVfcDTT),]
proteomeUVfcH2O2 <- proteomeUVfcH2O2[complete.cases(proteomeUVfcH2O2),]
proteomeUVfcSA <- proteomeUVfcSA[complete.cases(proteomeUVfcSA),]
mRNABPomeUVfcDTT <- mRNABPomeUVfcDTT[complete.cases(mRNABPomeUVfcDTT),]
mRNABPomeUVfcH2O2 <- mRNABPomeUVfcH2O2[complete.cases(mRNABPomeUVfcH2O2),]
mRNABPomeUVfcSA <- mRNABPomeUVfcSA[complete.cases(mRNABPomeUVfcSA),]
netChangesUVDTT <- netChangesUVDTT[complete.cases(netChangesUVDTT),]
netChangesUVH2O2 <- netChangesUVH2O2[complete.cases(netChangesUVH2O2),]
netChangesUVSA <- netChangesUVSA[complete.cases(netChangesUVSA),]

proteomeUVfcDTT <- proteomeUVfcDTT[order(proteomeUVfcDTT[,2],decreasing = T),]
proteomeUVfcH2O2 <- proteomeUVfcH2O2[order(proteomeUVfcH2O2[,2],decreasing = T),]
proteomeUVfcSA <- proteomeUVfcSA[order(proteomeUVfcSA[,2],decreasing = T),]
mRNABPomeUVfcDTT <- mRNABPomeUVfcDTT[order(mRNABPomeUVfcDTT[,2],decreasing = T),]
mRNABPomeUVfcH2O2 <- mRNABPomeUVfcH2O2[order(mRNABPomeUVfcH2O2[,2],decreasing = T),]
mRNABPomeUVfcSA <- mRNABPomeUVfcSA[order(mRNABPomeUVfcSA[,2],decreasing = T),]
netChangesUVDTT <- netChangesUVDTT[order(netChangesUVDTT[,2],decreasing = T),]
netChangesUVH2O2 <- netChangesUVH2O2[order(netChangesUVH2O2[,2],decreasing = T),]
netChangesUVSA <- netChangesUVSA[order(netChangesUVSA[,2],decreasing = T),]

saveTablesTsvExc(proteomeUVfcDTT,toPO4outdir)
saveTablesTsvExc(proteomeUVfcH2O2,toPO4outdir)
saveTablesTsvExc(proteomeUVfcSA,toPO4outdir)
saveTablesTsvExc(mRNABPomeUVfcDTT,toPO4outdir)
saveTablesTsvExc(mRNABPomeUVfcH2O2,toPO4outdir)
saveTablesTsvExc(mRNABPomeUVfcSA,toPO4outdir)
saveTablesTsvExc(netChangesUVDTT,toPO4outdir)
saveTablesTsvExc(netChangesUVH2O2,toPO4outdir)
saveTablesTsvExc(netChangesUVSA,toPO4outdir)

## 2.1- UV DTT Selection
bypValueProteomeUVDTT <- wholeDF$proteomeUV.pValue_Condition2 < protsPvalCut
byqValueProteomeUVDTT <- wholeDF$proteomeUV.qValue_Condition2 < protsPvalCut
byFCproteomeUVDTT <- abs(wholeDF$proteomeUV.log2ratio_Condition2) > protsFcCut
bypValuemRNABPomeUVDTT <- wholeDF$mRNABPomeUV.pValue_PolyARNAUVwithDTT < protsPvalCut
byqValuemRNABPomeUVDTT <- wholeDF$mRNABPomeUV.qValue_PolyARNAUVwithDTT < protsPvalCut
byFCmRNABPomeUVDTT <- abs(wholeDF$mRNABPomeUV.log2ratio_PolyARNAUVwithDTT) > protsFcCut
byNetChangesUVDTT <- abs(wholeDF$netChangesUVDTT) > protsFcCut

pValFCproteomeUVDTTselected <- na.omit(unique(wholeDF$UniprotACC[bypValueProteomeUVDTT & byFCproteomeUVDTT])); length(pValFCproteomeUVDTTselected)
qValFCproteomeUVDTTselected <- na.omit(unique(wholeDF$UniprotACC[byqValueProteomeUVDTT & byFCproteomeUVDTT])); length(qValFCproteomeUVDTTselected)
pValFCmRNABPomeUVDTTselected <- na.omit(unique(wholeDF$UniprotACC[bypValuemRNABPomeUVDTT & byFCmRNABPomeUVDTT])); length(pValFCmRNABPomeUVDTTselected)
qValFCmRNABPomeUVDTTselected <- na.omit(unique(wholeDF$UniprotACC[byqValuemRNABPomeUVDTT & byFCmRNABPomeUVDTT])); length(qValFCmRNABPomeUVDTTselected)
saveTablesTsvExc(pValFCproteomeUVDTTselected,toPO4outdir)
saveTablesTsvExc(qValFCproteomeUVDTTselected,toPO4outdir)
saveTablesTsvExc(pValFCmRNABPomeUVDTTselected,toPO4outdir)
saveTablesTsvExc(qValFCmRNABPomeUVDTTselected,toPO4outdir)

netChangesUVDTTprotNrbpomeSelPval <- na.omit(unique(wholeDF$UniprotACC[byNetChangesUVDTT & bypValueProteomeUVDTT & bypValuemRNABPomeUVDTT])) # 0
netChangesUVDTTprotNrbpomeSelQval <- na.omit(unique(wholeDF$UniprotACC[byNetChangesUVDTT & byqValueProteomeUVDTT & byqValuemRNABPomeUVDTT])) # 0
netChangesUVDTTSelFC <- na.omit(unique(wholeDF$UniprotACC[byNetChangesUVDTT]))

saveTablesTsvExc(netChangesUVDTTprotNrbpomeSelPval,toPO4outdir)
saveTablesTsvExc(netChangesUVDTTprotNrbpomeSelQval,toPO4outdir)
saveTablesTsvExc(netChangesUVDTTSelFC,toPO4outdir)

## 2.2- UV H2O2 Selection
bypValueProteomeUVH2O2 <- wholeDF$proteomeUV.pValue_Condition2 < protsPvalCut
byqValueProteomeUVH2O2 <- wholeDF$proteomeUV.qValue_Condition2 < protsPvalCut
byFCproteomeUVH2O2 <- abs(wholeDF$proteomeUV.log2ratio_Condition2) > protsFcCut
bypValuemRNABPomeUVH2O2 <- wholeDF$mRNABPomeUV.pValue_PolyARNAUVwithH2O2 < protsPvalCut
byqValuemRNABPomeUVH2O2 <- wholeDF$mRNABPomeUV.qValue_PolyARNAUVwithH2O2 < protsPvalCut
byFCmRNABPomeUVH2O2 <- abs(wholeDF$mRNABPomeUV.log2ratio_PolyARNAUVwithH202) > protsFcCut
byNetChangesUVH2O2 <- abs(wholeDF$netChangesUVH2O2) > protsFcCut

pValFCproteomeUVH2O2selected <- na.omit(unique(wholeDF$UniprotACC[bypValueProteomeUVH2O2 & byFCproteomeUVH2O2])); length(pValFCproteomeUVH2O2selected)
qValFCproteomeUVH2O2selected <- na.omit(unique(wholeDF$UniprotACC[byqValueProteomeUVH2O2 & byFCproteomeUVH2O2])); length(qValFCproteomeUVH2O2selected)
pValFCmRNABPomeUVH2O2selected <- na.omit(unique(wholeDF$UniprotACC[bypValuemRNABPomeUVH2O2 & byFCmRNABPomeUVH2O2])); length(pValFCmRNABPomeUVH2O2selected)
qValFCmRNABPomeUVH2O2selected <- na.omit(unique(wholeDF$UniprotACC[byqValuemRNABPomeUVH2O2 & byFCmRNABPomeUVH2O2])); length(qValFCmRNABPomeUVH2O2selected)
saveTablesTsvExc(pValFCproteomeUVH2O2selected,toPO4outdir)
saveTablesTsvExc(qValFCproteomeUVH2O2selected,toPO4outdir)
saveTablesTsvExc(pValFCmRNABPomeUVH2O2selected,toPO4outdir)
saveTablesTsvExc(qValFCmRNABPomeUVH2O2selected,toPO4outdir)

netChangesUVH2O2protNrbpomeSelPval <- na.omit(unique(wholeDF$UniprotACC[byNetChangesUVH2O2 & bypValueProteomeUVH2O2 & bypValuemRNABPomeUVH2O2])) # 0
netChangesUVH2O2protNrbpomeSelQval <- na.omit(unique(wholeDF$UniprotACC[byNetChangesUVH2O2 & byqValueProteomeUVH2O2 & byqValuemRNABPomeUVH2O2])) # 0
netChangesUVH2O2SelFC <- na.omit(unique(wholeDF$UniprotACC[byNetChangesUVH2O2]))

saveTablesTsvExc(netChangesUVH2O2protNrbpomeSelPval,toPO4outdir)
saveTablesTsvExc(netChangesUVH2O2protNrbpomeSelQval,toPO4outdir)
saveTablesTsvExc(netChangesUVH2O2SelFC,toPO4outdir)

## 2.3- UV SA Selection
bypValueProteomeUVSA <- wholeDF$proteomeUV.pValue_Condition2 < protsPvalCut
byqValueProteomeUVSA <- wholeDF$proteomeUV.qValue_Condition2 < protsPvalCut
byFCproteomeUVSA <- abs(wholeDF$proteomeUV.log2ratio_Condition2) > protsFcCut
bypValuemRNABPomeUVSA <- wholeDF$mRNABPomeUV.pValue_PolyARNAUVwithSA < protsPvalCut
byqValuemRNABPomeUVSA <- wholeDF$mRNABPomeUV.qValue_PolyARNAUVwithSA < protsPvalCut
byFCmRNABPomeUVSA <- abs(wholeDF$mRNABPomeUV.log2ratio_PolyARNAUVwithSA) > protsFcCut
byNetChangesUVSA <- abs(wholeDF$netChangesUVSA) > protsFcCut

pValFCproteomeUVSAselected <- na.omit(unique(wholeDF$UniprotACC[bypValueProteomeUVSA & byFCproteomeUVSA])); length(pValFCproteomeUVSAselected)
qValFCproteomeUVSAselected <- na.omit(unique(wholeDF$UniprotACC[byqValueProteomeUVSA & byFCproteomeUVSA])); length(qValFCproteomeUVSAselected)
pValFCmRNABPomeUVSAselected <- na.omit(unique(wholeDF$UniprotACC[bypValuemRNABPomeUVSA & byFCmRNABPomeUVSA])); length(pValFCmRNABPomeUVSAselected)
qValFCmRNABPomeUVSAselected <- na.omit(unique(wholeDF$UniprotACC[byqValuemRNABPomeUVSA & byFCmRNABPomeUVSA])); length(qValFCmRNABPomeUVSAselected)
saveTablesTsvExc(pValFCproteomeUVSAselected,toPO4outdir)
saveTablesTsvExc(qValFCproteomeUVSAselected,toPO4outdir)
saveTablesTsvExc(pValFCmRNABPomeUVSAselected,toPO4outdir)
saveTablesTsvExc(qValFCmRNABPomeUVSAselected,toPO4outdir)

netChangesUVSAprotNrbpomeSelPval <- na.omit(unique(wholeDF$UniprotACC[byNetChangesUVSA & bypValueProteomeUVSA & bypValuemRNABPomeUVSA])) # 0
netChangesUVSAprotNrbpomeSelQval <- na.omit(unique(wholeDF$UniprotACC[byNetChangesUVSA & byqValueProteomeUVSA & byqValuemRNABPomeUVSA])) # 0
netChangesUVSASelFC <- na.omit(unique(wholeDF$UniprotACC[byNetChangesUVSA]))

saveTablesTsvExc(netChangesUVSAprotNrbpomeSelPval,toPO4outdir)
saveTablesTsvExc(netChangesUVSAprotNrbpomeSelQval,toPO4outdir)
saveTablesTsvExc(netChangesUVSASelFC,toPO4outdir)

###########################
##### OUR DATA MAPPED #####
###########################

DEGs_file <- "/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/data/mapped/DEGsFullMapping.tsv"
mRNABPomeFAX_file <- "/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/data/mapped/mRNABPomeFAX.tsv"
mRNABPomeUV_file <- "/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/data/mapped/mRNABPomeUV.tsv"
proteomeFAX_file <- "/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/data/mapped/proteomeFAX.tsv"
proteomeNOX_file <- "/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/data/mapped/proteomeNOX.tsv"
proteomeUV_file <- "/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/data/mapped/proteomeUV.tsv"

DEGs_DF <- read.delim("/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/data/mapped/DEGsFullMapping.tsv",quote = "")
mRNABPomeFAX_DF <- read.delim("/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/data/mapped/mRNABPomeFAX.tsv",quote = "")
mRNABPomeUV_DF <- read.delim("/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/data/mapped/mRNABPomeUV.tsv",quote = "")
proteomeFAX_DF <- read.delim("/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/data/mapped/proteomeFAX.tsv",quote = "")
proteomeNOX_DF <- read.delim("/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/data/mapped/proteomeNOX.tsv",quote = "")
proteomeUV_DF <- read.delim("/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/data/mapped/proteomeUV.tsv",quote = "")

unique(c(colnames(DEGs_DF),colnames(mRNABPomeFAX_DF),colnames(mRNABPomeUV_DF),colnames(proteomeFAX_DF),colnames(proteomeNOX_DF),colnames(proteomeUV_DF)))

colnames(yeastGenesProtMap)
colnames(DEGs_DF)
which(colnames(DEGs_DF) %in% colnames(yeastGenesProtMap))
colnames(DEGs_DF)[which(!colnames(DEGs_DF) %in% colnames(yeastGenesProtMap))] <- paste0("DEGs.",colnames(DEGs_DF)[which(!colnames(DEGs_DF) %in% colnames(yeastGenesProtMap))])
colnames(proteomeFAX_DF)[which(!colnames(proteomeFAX_DF) %in% colnames(yeastGenesProtMap))] <- paste0("proteomeFAX.",colnames(proteomeFAX_DF)[which(!colnames(proteomeFAX_DF) %in% colnames(yeastGenesProtMap))])
colnames(proteomeNOX_DF)[which(!colnames(proteomeNOX_DF) %in% colnames(yeastGenesProtMap))] <- paste0("proteomeNOX.",colnames(proteomeNOX_DF)[which(!colnames(proteomeNOX_DF) %in% colnames(yeastGenesProtMap))])
colnames(proteomeUV_DF)[which(!colnames(proteomeUV_DF) %in% colnames(yeastGenesProtMap))] <- paste0("proteomeUV.",colnames(proteomeUV_DF)[which(!colnames(proteomeUV_DF) %in% colnames(yeastGenesProtMap))])
colnames(mRNABPomeFAX_DF)[which(!colnames(mRNABPomeFAX_DF) %in% colnames(yeastGenesProtMap))] <- paste0("mRNABPomeFAX.",colnames(mRNABPomeFAX_DF)[which(!colnames(mRNABPomeFAX_DF) %in% colnames(yeastGenesProtMap))])
colnames(mRNABPomeUV_DF)[which(!colnames(mRNABPomeUV_DF) %in% colnames(yeastGenesProtMap))] <- paste0("mRNABPomeUV.",colnames(mRNABPomeUV_DF)[which(!colnames(mRNABPomeUV_DF) %in% colnames(yeastGenesProtMap))])

wholeDataDataFrame <- merge(DEGs_DF,proteomeFAX_DF,all=TRUE)
wholeDataDataFrame <- merge(wholeDataDataFrame,proteomeNOX_DF,all=TRUE)
wholeDataDataFrame <- merge(wholeDataDataFrame,proteomeUV_DF,all=TRUE)
wholeDataDataFrame <- merge(wholeDataDataFrame,mRNABPomeFAX_DF,all=TRUE)
wholeDataDataFrame <- merge(wholeDataDataFrame,mRNABPomeUV_DF,all=TRUE)

outexcel <- "/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/data/mapped/allDataMapped_together.xlsx"
write.xlsx(wholeDataDataFrame,outexcel, asTable = FALSE, overwrite = TRUE)
write.table(wholeDataDataFrame,gsub('.xlsx','.tsv',outexcel),quote = F,sep = '\t',col.names = T,row.names = F)


columnsToremove <- c("Gene.Names..ordered.locus.","Gene.Names..ORF.","Gene.Names..primary.")

outexcel <- "/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/data/mapped/allDataMapped_separated.xlsx"
write.xlsx(list(DEGs=DEGs_DF,mRNABPomeFAX=mRNABPomeFAX_DF,mRNABPomeUV=mRNABPomeUV_DF,
                proteomeFAX=proteomeFAX_DF,proteomeNOX=proteomeNOX_DF,proteomeUV=proteomeUV_DF),
           outexcel, asTable = FALSE, overwrite = TRUE)

outexcel <- "/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/data/yeastReference/mappingFile.xlsx"
write.xlsx(yeastGenesProtMap,outexcel, asTable = FALSE, overwrite = TRUE)

####################
##### OUR DATA #####
####################

DEGsFile <- "/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/data/RNA-Seq data/result/DE gene testing statistics.csv"
proteomeFAXfile <- "/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/data/P445_Prot_DIA_iBAQ_noOutliers_nocomb_20220826_Final/P445_Prot_DIA_iBAQ_noOutliers_nocomb_20220826_FAX_Final/P445_Prot_DIA_iBAQ_noOutliers_nocomb_20220826_FAX_Final_PROTEIN.tsv"
proteomeNOXfile <- "/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/data/P445_Prot_DIA_iBAQ_noOutliers_nocomb_20220826_Final/P445_Prot_DIA_iBAQ_noOutliers_nocomb_20220826_NOX_Final/P445_Prot_DIA_iBAQ_noOutliers_nocomb_20220826_NOX_Final_PROTEIN.tsv"  
proteomeUVfile <- "/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/data/P445_Prot_DIA_iBAQ_noOutliers_nocomb_20220826_Final/P445_Prot_DIA_iBAQ_noOutliers_nocomb_20220826_UV_Final/P445_Prot_DIA_iBAQ_noOutliers_nocomb_20220826_UV_Final_PROTEIN.tsv"
mRNABPomeFAXfile <- "/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/data/P445_AG_PolyA_40_all_3_cond_20220611_Final/P445_AG_PolyA_40_20220617_gmin_FAX_p2_Final/P445_AG_PolyA_40_20220617_gmin_FAX_p2_PROTEIN.tsv"
mRNABPomeUVffile <- "/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/data/P445_AG_PolyA_40_all_3_cond_20220611_Final/P445_AG_PolyA_40_20220617_gmin_UV_p2_no511_Final/P445_AG_PolyA_40_20220617_gmin_UV_p2_no511_PROTEIN.tsv"
#? mRNABPomeNOXfile <- "/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/data/P445_AG_PolyA_40_all_3_cond_20220611_Final/P445_AG_PolyA_40_20220617_gmin_UV_nonorm_p2_NOX_Final/P445_AG_PolyA_40_20220617_gmin_UV_nonorm_p2_NOX_PROTEIN.tsv

DEGsDF <- read.delim(DEGsFile, sep = ",") 
proteomeFAXDF <- read.delim(proteomeFAXfile,quote = "")
proteomeNOXDF <- read.delim(proteomeNOXfile,quote = "")
proteomeUVDF <- read.delim(proteomeUVfile,quote = "")
mRNABPomeFAXDF <- read.delim(mRNABPomeFAXfile,quote = "")
mRNABPomeUVDF <- read.delim(mRNABPomeUVffile,quote = "")

###################
##### MAPPING #####
###################
# proteomeFAXDF
proteomeFAXDFMaped <- merge(proteomeFAXDF,yeastGenesProtMap,by.x="ac",by.y="UniprotACC")
proteomeFAXDFMaped$UniprotACC <- proteomeFAXDFMaped$ac
any(duplicated(proteomeFAXDFMaped))
proteomeFAXDFMaped <- proteomeFAXDFMaped[!duplicated(proteomeFAXDFMaped),]
proteomeFAXDFNotMapped <- proteomeFAXDF[!(proteomeFAXDF$ac %in% proteomeFAXDFMaped$ac),]

proteomeFAXDFMaped2 <- merge(proteomeFAXDFNotMapped,yeastGenesProtMap,by.x="geneName",by.y="GeneName")
proteomeFAXDFMaped2$GeneName <- proteomeFAXDFMaped2$geneName

proteomeFAXDFMaped_full <- rbind(proteomeFAXDFMaped,proteomeFAXDFMaped2)

proteomeFAXMapedFullFile <- '/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/data/mapped/proteomeFAX.tsv'
write.table(proteomeFAXDFMaped_full,proteomeFAXMapedFullFile,quote = F,sep = '\t',col.names = T,row.names = F)

# proteomeNOXDF
proteomeNOXDFMaped <- merge(proteomeNOXDF,yeastGenesProtMap,by.x="ac",by.y="UniprotACC")
proteomeNOXDFMaped$UniprotACC <- proteomeNOXDFMaped$ac
any(duplicated(proteomeNOXDFMaped))
proteomeNOXDFMaped <- proteomeNOXDFMaped[!duplicated(proteomeNOXDFMaped),]
proteomeNOXDFNotMapped <- proteomeNOXDF[!(proteomeNOXDF$ac %in% proteomeNOXDFMaped$ac),]

proteomeNOXDFMaped2 <- merge(proteomeNOXDFNotMapped,yeastGenesProtMap,by.x="geneName",by.y="GeneName")
proteomeNOXDFMaped2$GeneName <- proteomeNOXDFMaped2$geneName

proteomeNOXDFMaped_full <- rbind(proteomeNOXDFMaped,proteomeNOXDFMaped2)

proteomeNOXMapedFullFile <- '/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/data/mapped/proteomeNOX.tsv'
write.table(proteomeNOXDFMaped_full,proteomeNOXMapedFullFile,quote = F,sep = '\t',col.names = T,row.names = F)

# proteomeUVDF
proteomeUVDFMaped <- merge(proteomeUVDF,yeastGenesProtMap,by.x="ac",by.y="UniprotACC")
proteomeUVDFMaped$UniprotACC <- proteomeUVDFMaped$ac
any(duplicated(proteomeUVDFMaped))
proteomeUVDFMaped <- proteomeUVDFMaped[!duplicated(proteomeUVDFMaped),]
proteomeUVDFNotMapped <- proteomeUVDF[!(proteomeUVDF$ac %in% proteomeUVDFMaped$ac),]

proteomeUVDFMaped2 <- merge(proteomeUVDFNotMapped,yeastGenesProtMap,by.x="geneName",by.y="GeneName")
proteomeUVDFMaped2$GeneName <- proteomeUVDFMaped2$geneName

proteomeUVDFMaped_full <- rbind(proteomeUVDFMaped,proteomeUVDFMaped2)

proteomeUVMapedFullFile <- '/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/data/mapped/proteomeUV.tsv'
write.table(proteomeUVDFMaped_full,proteomeUVMapedFullFile,quote = F,sep = '\t',col.names = T,row.names = F)

# mRNABPomeFAXDF
mRNABPomeFAXDF$isoformsAC <- mRNABPomeFAXDF$ac
mRNABPomeFAXDF$ac <- gsub('-.*','',mRNABPomeFAXDF$ac)
mRNABPomeFAXDFMaped <- merge(mRNABPomeFAXDF,yeastGenesProtMap,by.x="ac",by.y="UniprotACC")
mRNABPomeFAXDFMaped$UniprotACC <- mRNABPomeFAXDFMaped$ac
mRNABPomeFAXDFMaped <- mRNABPomeFAXDFMaped[!duplicated(mRNABPomeFAXDFMaped),]
mRNABPomeFAXDFNotMapped <- mRNABPomeFAXDF[!(mRNABPomeFAXDF$ac %in% mRNABPomeFAXDFMaped$ac),]
nrow(mRNABPomeFAXDFNotMapped)
mRNABPomeFAXMapedFullFile <- '/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/data/mapped/mRNABPomeFAX.tsv'
write.table(mRNABPomeFAXDFMaped,mRNABPomeFAXMapedFullFile,quote = F,sep = '\t',col.names = T,row.names = F)

# mRNABPomeUVfDF
mRNABPomeUVDF$isoformsAC <- mRNABPomeUVDF$ac
mRNABPomeUVDF$ac <- gsub('-.*','',mRNABPomeUVDF$ac)
mRNABPomeUVDFMaped <- merge(mRNABPomeUVDF,yeastGenesProtMap,by.x="ac",by.y="UniprotACC")
mRNABPomeUVDFMaped$UniprotACC <- mRNABPomeUVDFMaped$ac
mRNABPomeUVDFMaped <- mRNABPomeUVDFMaped[!duplicated(mRNABPomeUVDFMaped),]
mRNABPomeUVDFNotMapped <- mRNABPomeUVDF[!(mRNABPomeUVDF$ac %in% mRNABPomeUVDFMaped$ac),]
nrow(mRNABPomeUVDFNotMapped)
mRNABPomeUVMapedFullFile <- '/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/data/mapped/mRNABPomeUV.tsv'
write.table(mRNABPomeUVDFMaped,mRNABPomeUVMapedFullFile,quote = F,sep = '\t',col.names = T,row.names = F)

########################
##### MAPPING DEGs #####
########################
DEGsDF <- read.delim(DEGsFile, sep = ",") 
DEGsDF <- reshape(DEGsDF, direction = "wide", idvar = "target", timevar = "contrast")
#colnames(DEGsDF) <- gsub("X","",gsub("\\.","",colnames(DEGsDF))); colnames(DEGsDF)

DEGsDFMaped <- merge(DEGsDF,yeastGenesProtMap,by.x="target",by.y="ORFgene")
DEGsDFMaped$ORFgene <- DEGsDFMaped$target
any(duplicated(DEGsDFMaped))
DEGsDFMaped <- DEGsDFMaped[!duplicated(DEGsDFMaped),]
DEGsNotMapped <- DEGsDF[!(DEGsDF$target %in% DEGsDFMaped$target),]

DEGsNotMapped$target <- gsub("_.*","",toupper(DEGsNotMapped$target))
DEGsDFMaped2 <- merge(DEGsNotMapped,yeastGenesProtMap,by.x="target",by.y="ORFgene")
DEGsDFMaped2$ORFgene <- DEGsDFMaped2$target
any(duplicated(DEGsDFMaped2))
DEGsDFMaped2 <- DEGsDFMaped2[!duplicated(DEGsDFMaped2),]
DEGsNotMapped2 <- DEGsNotMapped[!(DEGsNotMapped$target %in% DEGsDFMaped2$target),]

DEGsDFMaped3 <- merge(DEGsNotMapped2,yeastGenesProtMap,by.x="target",by.y="GeneName")
any(duplicated(DEGsDFMaped3))
DEGsDFMaped3$GeneName <- DEGsDFMaped3$target
DEGsDFMaped3 <- DEGsDFMaped3[!duplicated(DEGsDFMaped3),]
DEGsNotMapped3 <- DEGsNotMapped2[!(DEGsNotMapped2$target %in% DEGsDFMaped3$target),]
nrow(DEGsNotMapped3)
DEGsDFMaped4 <- cbind(DEGsNotMapped3,yeastGenesProtMap[grepl(paste(DEGsNotMapped3$target,collapse = '|'),yeastGenesProtMap$Alias),])
DEGsNotMapped4 <- DEGsNotMapped3[!(DEGsNotMapped3$target %in% DEGsDFMaped4$target),]
nrow(DEGsNotMapped4)

DEGsDFMapedFull <- rbind(DEGsDFMaped,DEGsDFMaped2,DEGsDFMaped3,DEGsDFMaped4)
any(duplicated(DEGsDFMapedFull))
colnames(DEGsDFMapedFull)

DEGsDFMapedFull <- DEGsDFMapedFull[c("target","adj.pval.H202-YPD","log2FC.H202-YPD","adj.pval.DTT-YPD","log2FC.DTT-YPD","adj.pval.SA-YPD","log2FC.SA-YPD",'seqnames','start','end','type','ORFgene','ORF','SGDID','GeneName','UniprotACC','UniprotName','Gene.Names','Protein.names','Gene.Names..ordered.locus.','Gene.Names..ORF.','Gene.Names..primary.','Gene.Names..synonym.','Alias')]

DEGsDFMapedFullFile <- '/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/data/mapped/DEGsFullMapping.tsv'
write.table(DEGsDFMapedFull,DEGsDFMapedFullFile,quote = F,sep = '\t',col.names = T,row.names = F)


#################################
##### CREATING MAPPING DATA #####
#################################

yeastdbGFFfile <- "/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/data/yeastReference/S288C_reference_genome_Current_Release/S288C_reference_genome_R64-3-1_20210421/NOFASTAsaccharomyces_cerevisiae_R64-3-1_20210421.gff"
yeastdbGFFdf <- read.table(yeastdbGFFfile,quote = "",sep = "\t",comment.char = "#",header = F,as.is = T,fill = T,encoding = "utf-8")
table(yeastdbGFFdf$V7)
yeastdbGFFdf$V7[yeastdbGFFdf$V7 == 0] <- "*"
yeastdbGFFdf$V7[yeastdbGFFdf$V7 == "."] <- "*"
table(yeastdbGFFdf$V7)

yeastdbGFFmyfile <- '/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/data/yeastReference/S288C_reference_genome_Current_Release/S288C_reference_genome_R64-3-1_20210421/myyeast.gff'
write.table(yeastdbGFFdf,yeastdbGFFmyfile,col.names = F,row.names = F, sep = "\t",quote = F)
yeastdbGFFdf2 <- read.delim(yeastdbGFFmyfile,quote = "",header = F,sep = "\t",comment.char = "#",stringsAsFactors = F)
table(yeastdbGFFdf2$V7)


yeastdbGFF <- as.data.frame(rtracklayer::import(yeastdbGFFmyfile))
paste(colnames(yeastdbGFF),collapse = "','")
allcolumns <- c('seqnames','start','end','width','strand','source','type','score','phase','ID','dbxref','Name','Note','display','curie','Parent','Ontology_term','orf_classification','protein_id','Alias','gene','transcript_id','conditions')
table(yeastdbGFF$score)
table(yeastdbGFF$type)

micron2gffile <- "/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/data/yeastReference/NOFASTAscerevisiae_2-micron.gff"
micron2gffDF <- as.data.frame(rtracklayer::import(micron2gffile))

yeastdbGFF <- rbind.fill(yeastdbGFF,micron2gffDF)

yeastdbGFF_genes <- yeastdbGFF[yeastdbGFF$type != "CDS",]
table(yeastdbGFF_genes$protein_id)
yeastdbGFF_genes$protein_id <- NULL

yeastdbGFF_proteins <- yeastdbGFF[yeastdbGFF$type == "CDS",]
yeastdbGFF_proteins$Name <- gsub("_CDS","",yeastdbGFF_proteins$Name)
yeastdbGFF_proteinsSub <- yeastdbGFF_proteins[c("Name","protein_id")]

geneNprots <- merge(yeastdbGFF_genes,yeastdbGFF_proteinsSub,by.x="ID",by.y="Name",all.x=T)
geneNprots$dbxref <- gsub("SGD:","",geneNprots$dbxref)
geneNprots$protein_id <- gsub("UniProtKB:","",geneNprots$protein_id)

interestingCols <- c('seqnames','start','end','type','ID','dbxref','protein_id','Alias','gene')
selectedMapping <- geneNprots[interestingCols]
selectedMapping

uniprotFile <- "/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/data/yeastReference/uniprot-compressed_true_download_true_fields_accession_2Cid_2Cprotei-2022.09.07-13.46.19.80.tsv"
uniprotDF <- read.delim(uniprotFile)
fullmerge <- merge(selectedMapping,uniprotDF,by.x="protein_id",by.y="Entry",all.x=T)
fullmerge$ORFgene <- gsub("-.*","",fullmerge$ID)

paste(colnames(fullmerge),collapse = "','")
fullmerge <- fullmerge[c('seqnames','start','end','type','ORFgene','ID','dbxref','gene','protein_id','Entry.Name','Gene.Names','Protein.names','Gene.Names..ordered.locus.','Gene.Names..ORF.','Gene.Names..primary.','Gene.Names..synonym.','Alias')]
fullmerge$Alias <- unlist(lapply(fullmerge$Alias,function(x) paste0(x,collapse = ' | ')))
colnames(fullmerge) <- c('seqnames','start','end','type','ORFgene','ORF','SGDID','GeneName','UniprotACC','UniprotName','Gene.Names','Protein.names','Gene.Names..ordered.locus.','Gene.Names..ORF.','Gene.Names..primary.','Gene.Names..synonym.','Alias')

mappingFinalFile <- "/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/data/yeastReference/mappingFile.tsv"
write.table(fullmerge,mappingFinalFile,quote = F,sep = '\t',col.names = T)

mappingFinalFile <- "/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/data/yeastReference/mappingFile.tsv"
yeastGenesProtMap <- read.delim(mappingFinalFile,quote = "")


# library(biomaRt)
# mymart <- useMart("ensembl")#; View(listFilters(humanmart)); View(listAttributes(humanmart))
# G_list <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol"),values=ensgWithoutMapping,mart= humanmart)

