##########################
##### FILTERING DATA #####
##########################
wholeDFFile <- "/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/FinalData/wholeDF.tsv"
wholeDF <- read.delim(wholeDFFile,quote = "")

outdir <- "/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/FinalData/GeneLists"

IDs <- c("UniprotACC","GeneName","Protein.names") # rbpomeID <- IDs
DEGpvalCut <- 0.01
DEGfcCut   <- 1

protsPvalCut <- 0.05
protsFcCut <- 1 

# 1- Differentially Expressed Genes

wholeDF$DEGs.log2FC.SA.YPD
wholeDF$DEGs.log2FC.SA.YPD
wholeDF$ProteomeFAX.log2ratio_FAXwithSA
wholeDF$ProteomeFAX.qValue_FAXwithSA
wholeDF$RBPomeFAX.log2ratio_PolyARNAFAXwithSA
wholeDF$RBPomeFAX.qValue_PolyARNAFAXwithSA
wholeDF$FAXnetchangesSA

DEGsH202 <- wholeDF[wholeDF$DEGs.adj.pval.H202.YPD < DEGpvalCut & abs(wholeDF$DEGs.log2FC.H202.YPD) > DEGfcCut,"GeneName"]
DEGsDTT <- wholeDF[wholeDF$DEGs.adj.pval.DTT.YPD < DEGpvalCut  & abs(wholeDF$DEGs.log2FC.DTT.YPD) > DEGfcCut,"GeneName"]
DEGsSA <- wholeDF[wholeDF$DEGs.adj.pval.SA.YPD < DEGpvalCut & abs(wholeDF$DEGs.log2FC.SA.YPD) > DEGfcCut,"GeneName"]
outdir <- "/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/data/geneLists/GSA"
saveTablesTsvExc(DEGsH202,outdir)
saveTablesTsvExc(DEGsDTT,outdir)
saveTablesTsvExc(DEGsSA,outdir)

DEGsH202fc <- wholeDF[order(wholeDF$DEGs.log2FC.H202.YPD,decreasing = T),c(IDs,"DEGs.log2FC.H202.YPD")]
DEGsDTTfc <- wholeDF[order(wholeDF$DEGs.log2FC.DTT.YPD,decreasing = T),c(IDs,"DEGs.log2FC.DTT.YPD")]
DEGsSAfc <- wholeDF[order(wholeDF$DEGs.log2FC.SA.YPD,decreasing = T),c(IDs,"DEGs.log2FC.SA.YPD")]
outdir <- "/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/data/geneLists/GSEA"
saveTablesTsvExc(DEGsH202fc,outdir,excel=F)
saveTablesTsvExc(DEGsDTTfc,outdir,excel=F)
saveTablesTsvExc(DEGsSAfc,outdir,excel=F)

# 2- Proteins and mRNAB Proteins FAX n UVW

########################
#### CUT-OFFs Prots ####
########################

protFAX.DTT.FCsig <- abs(wholeDF$proteomeFAX.log2ratio_Condition2) > protsFcCut
protFAX.H2O2.FCsig <- abs(wholeDF$proteomeFAX.log2ratio_Condition3) > protsFcCut
protFAX.SA.FCsig <- abs(wholeDF$proteomeFAX.log2ratio_Condition4) > protsFcCut
mrbpFAX.DTT.FCsig <- abs(wholeDF$mRNABPomeFAX.log2ratio_PolyARNAFAXwithDTT) > protsFcCut
mrbpFAX.H2O2.FCsig <- abs(wholeDF$mRNABPomeFAX.log2ratio_PolyARNAFAXwithH2O2) > protsFcCut
mrbpFAX.SA.FCsig <- abs(wholeDF$mRNABPomeFAX.log2ratio_PolyARNAFAXwithSA) > protsFcCut
protUV.DTT.FCsig <- abs(wholeDF$proteomeUV.log2ratio_Condition2) > protsFcCut
protUV.H2O2.FCsig <- abs(wholeDF$proteomeUV.log2ratio_Condition3) > protsFcCut
protUV.SA.FCsig <- abs(wholeDF$proteomeUV.log2ratio_Condition4) > protsFcCut
mrbpUV.DTT.FCsig <- abs(wholeDF$mRNABPomeUV.log2ratio_PolyARNAUVwithDTT) > protsFcCut
mrbpUV.H2O2.FCsig <- abs(wholeDF$mRNABPomeUV.log2ratio_PolyARNAUVwithH202) > protsFcCut
mrbpUV.SA.FCsig <- abs(wholeDF$mRNABPomeUV.log2ratio_PolyARNAUVwithSA) > protsFcCut

netFAX.DTT.FCsig <- abs(wholeDF$netChangesFAXDTT) > protsFcCut
netFAX.H2O2.FCsig <- abs(wholeDF$netChangesFAXH2O2) > protsFcCut
netFAX.SA.FCsig <- abs(wholeDF$netChangesFAXSA) > protsFcCut
netUV.DTT.FCsig <- abs(wholeDF$netChangesUVDTT) > protsFcCut
netUV.H2O2.FCsig <- abs(wholeDF$netChangesUVH2O) > protsFcCut
netUV.SA.FCsig <- abs(wholeDF$netChangesUVSA) > protsFcCut

protFAX.DTT.QValsig <- wholeDF$proteomeFAX.qValue_Condition2 < protsPvalCut
protFAX.H2O2.QValsig <- wholeDF$proteomeFAX.qValue_Condition3 < protsPvalCut
protFAX.SA.QValsig <- wholeDF$proteomeFAX.qValue_Condition4 < protsPvalCut
mrbpFAX.DTT.QValsig <- wholeDF$mRNABPomeFAX.qValue_PolyARNAFAXwithDTT < protsPvalCut
mrbpFAX.H2O2.QValsig <- wholeDF$mRNABPomeFAX.qValue_PolyARNAFAXwithH2O2 < protsPvalCut
mrbpFAX.SA.QValsig <- wholeDF$mRNABPomeFAX.qValue_PolyARNAFAXwithSA < protsPvalCut
protUV.DTT.QValsig <- wholeDF$proteomeUV.qValue_Condition2 < protsPvalCut
protUV.H2O2.QValsig <- wholeDF$proteomeUV.qValue_Condition3 < protsPvalCut
protUV.SA.QValsig <- wholeDF$proteomeUV.qValue_Condition4 < protsPvalCut
mrbpUV.DTT.QValsig <- wholeDF$mRNABPomeUV.qValue_PolyARNAUVwithDTT < protsPvalCut
mrbpUV.H2O2.QValsig <- wholeDF$mRNABPomeUV.qValue_PolyARNAUVwithH202 < protsPvalCut
mrbpUV.SA.QValsig <- wholeDF$mRNABPomeUV.qValue_PolyARNAUVwithSA < protsPvalCut

protFAX.DTT.PValsig <- wholeDF$proteomeFAX.pValue_Condition2 < protsPvalCut
protFAX.H2O2.PValsig <- wholeDF$proteomeFAX.pValue_Condition3 < protsPvalCut
protFAX.SA.PValsig <- wholeDF$proteomeFAX.pValue_Condition4 < protsPvalCut
mrbpFAX.DTT.PValsig <- wholeDF$mRNABPomeFAX.pValue_PolyARNAFAXwithDTT < protsPvalCut
mrbpFAX.H2O2.PValsig <- wholeDF$mRNABPomeFAX.pValue_PolyARNAFAXwithH2O2 < protsPvalCut
mrbpFAX.SA.PValsig <- wholeDF$mRNABPomeFAX.pValue_PolyARNAFAXwithSA < protsPvalCut
protUV.DTT.PValsig <- wholeDF$proteomeUV.pValue_Condition2 < protsPvalCut
protUV.H2O2.PValsig <- wholeDF$proteomeUV.pValue_Condition3 < protsPvalCut
protUV.SA.PValsig <- wholeDF$proteomeUV.pValue_Condition4 < protsPvalCut
mrbpUV.DTT.PValsig <- wholeDF$mRNABPomeUV.pValue_PolyARNAUVwithDTT < protsPvalCut
mrbpUV.H2O2.PValsig <- wholeDF$mRNABPomeUV.pValue_PolyARNAUVwithH202 < protsPvalCut
mrbpUV.SA.PValsig <- wholeDF$mRNABPomeUV.pValue_PolyARNAUVwithSA < protsPvalCut

protFAX.DTT.Psig <- unique(na.omit(wholeDF[protFAX.DTT.FCsig & protFAX.DTT.PValsig, "UniprotACC"]))
protFAX.H2O2.Psig <- unique(na.omit(wholeDF[protFAX.H2O2.FCsig & protFAX.H2O2.PValsig, "UniprotACC"]))
protFAX.SA.Psig <- unique(na.omit(wholeDF[protFAX.SA.FCsig & protFAX.SA.PValsig, "UniprotACC"]))
mrbpFAX.DTT.Psig <- unique(na.omit(wholeDF[mrbpFAX.DTT.FCsig & mrbpFAX.DTT.PValsig, "UniprotACC"]))
mrbpFAX.H2O2.Psig <- unique(na.omit(wholeDF[mrbpFAX.H2O2.FCsig & mrbpFAX.H2O2.PValsig, "UniprotACC"]))
mrbpFAX.SA.Psig <- unique(na.omit(wholeDF[mrbpFAX.SA.FCsig & mrbpFAX.SA.PValsig, "UniprotACC"]))
protUV.DTT.Psig <- unique(na.omit(wholeDF[protUV.DTT.FCsig & protUV.DTT.PValsig, "UniprotACC"]))
protUV.H2O2.Psig <- unique(na.omit(wholeDF[protUV.H2O2.FCsig & protUV.H2O2.PValsig, "UniprotACC"]))
protUV.SA.Psig <- unique(na.omit(wholeDF[protUV.SA.FCsig & protUV.SA.PValsig, "UniprotACC"]))
mrbpUV.DTT.Psig <- unique(na.omit(wholeDF[mrbpUV.DTT.FCsig & mrbpUV.DTT.PValsig, "UniprotACC"]))
mrbpUV.H2O2.Psig <- unique(na.omit(wholeDF[mrbpUV.H2O2.FCsig & mrbpUV.H2O2.PValsig, "UniprotACC"]))
mrbpUV.SA.Psig <- unique(na.omit(wholeDF[mrbpUV.SA.FCsig & mrbpUV.SA.PValsig, "UniprotACC"]))
netFAX.DTT.PsigPnR <- unique(na.omit(wholeDF[netFAX.DTT.FCsig & mrbpFAX.DTT.PValsig  & protFAX.DTT.PValsig, "UniprotACC"]))
netFAX.H2O2.PsigPnR <- unique(na.omit(wholeDF[netFAX.H2O2.FCsig & mrbpFAX.H2O2.PValsig & protFAX.H2O2.PValsig, "UniprotACC"]))
netFAX.SA.PsigPnR <- unique(na.omit(wholeDF[netFAX.SA.FCsig & mrbpFAX.SA.PValsig   & protFAX.SA.PValsig, "UniprotACC"]))
netUV.DTT.PsigPnR <- unique(na.omit(wholeDF[netUV.DTT.FCsig  & mrbpUV.DTT.PValsig   & protUV.DTT.PValsig, "UniprotACC"]))
netUV.H2O2.PsigPnR <- unique(na.omit(wholeDF[netUV.H2O2.FCsig  & mrbpUV.H2O2.PValsig  & protUV.H2O2.PValsig, "UniprotACC"]))
netUV.SA.PsigPnR <- unique(na.omit(wholeDF[netUV.SA.FCsig  & mrbpUV.SA.PValsig    & protUV.SA.PValsig, "UniprotACC"]))
netFAX.DTT.PsigP <- unique(na.omit(wholeDF[netFAX.DTT.FCsig & protFAX.DTT.PValsig, "UniprotACC"]))
netFAX.H2O2.PsigP <- unique(na.omit(wholeDF[netFAX.H2O2.FCsig & protFAX.H2O2.PValsig, "UniprotACC"]))
netFAX.SA.PsigP <- unique(na.omit(wholeDF[netFAX.SA.FCsig & protFAX.SA.PValsig, "UniprotACC"]))
netUV.DTT.PsigP <- unique(na.omit(wholeDF[netUV.DTT.FCsig  & protUV.DTT.PValsig, "UniprotACC"]))
netUV.H2O2.PsigP <- unique(na.omit(wholeDF[netUV.H2O2.FCsig  & protUV.H2O2.PValsig, "UniprotACC"]))
netUV.SA.PsigP <- unique(na.omit(wholeDF[netUV.SA.FCsig  & protUV.SA.PValsig, "UniprotACC"]))
netFAX.DTT.PsigR <- unique(na.omit(wholeDF[netFAX.DTT.FCsig & mrbpFAX.DTT.PValsig, "UniprotACC"]))
netFAX.H2O2.PsigR <- unique(na.omit(wholeDF[netFAX.H2O2.FCsig & mrbpFAX.H2O2.PValsig, "UniprotACC"]))
netFAX.SA.PsigR <- unique(na.omit(wholeDF[netFAX.SA.FCsig & mrbpFAX.SA.PValsig, "UniprotACC"]))
netUV.DTT.PsigR <- unique(na.omit(wholeDF[netUV.DTT.FCsig  & mrbpUV.DTT.PValsig, "UniprotACC"]))
netUV.H2O2.PsigR <- unique(na.omit(wholeDF[netUV.H2O2.FCsig  & mrbpUV.H2O2.PValsig, "UniprotACC"]))
netUV.SA.PsigR <- unique(na.omit(wholeDF[netUV.SA.FCsig  & mrbpUV.SA.PValsig, "UniprotACC"]))
protFAX.DTT.Qsig <- unique(na.omit(wholeDF[protFAX.DTT.FCsig & protFAX.DTT.QValsig, "UniprotACC"]))
protFAX.H2O2.Qsig <- unique(na.omit(wholeDF[protFAX.H2O2.FCsig & protFAX.H2O2.QValsig, "UniprotACC"]))
protFAX.SA.Qsig <- unique(na.omit(wholeDF[protFAX.SA.FCsig & protFAX.SA.QValsig, "UniprotACC"]))
mrbpFAX.DTT.Qsig <- unique(na.omit(wholeDF[mrbpFAX.DTT.FCsig & mrbpFAX.DTT.QValsig, "UniprotACC"]))
mrbpFAX.H2O2.Qsig <- unique(na.omit(wholeDF[mrbpFAX.H2O2.FCsig & mrbpFAX.H2O2.QValsig, "UniprotACC"]))
mrbpFAX.SA.Qsig <- unique(na.omit(wholeDF[mrbpFAX.SA.FCsig & mrbpFAX.SA.QValsig, "UniprotACC"]))
protUV.DTT.Qsig <- unique(na.omit(wholeDF[protUV.DTT.FCsig & protUV.DTT.QValsig, "UniprotACC"]))
protUV.H2O2.Qsig <- unique(na.omit(wholeDF[protUV.H2O2.FCsig & protUV.H2O2.QValsig, "UniprotACC"]))
protUV.SA.Qsig <- unique(na.omit(wholeDF[protUV.SA.FCsig & protUV.SA.QValsig, "UniprotACC"]))
mrbpUV.DTT.Qsig <- unique(na.omit(wholeDF[mrbpUV.DTT.FCsig & mrbpUV.DTT.QValsig, "UniprotACC"]))
mrbpUV.H2O2.Qsig <- unique(na.omit(wholeDF[mrbpUV.H2O2.FCsig & mrbpUV.H2O2.QValsig, "UniprotACC"]))
mrbpUV.SA.Qsig <- unique(na.omit(wholeDF[mrbpUV.SA.FCsig & mrbpUV.SA.QValsig, "UniprotACC"]))
netFAX.DTT.QsigPnR <- unique(na.omit(wholeDF[netFAX.DTT.FCsig & mrbpFAX.DTT.QValsig  & protFAX.DTT.QValsig, "UniprotACC"]))
netFAX.H2O2.QsigPnR <- unique(na.omit(wholeDF[netFAX.H2O2.FCsig & mrbpFAX.H2O2.QValsig & protFAX.H2O2.QValsig, "UniprotACC"]))
netFAX.SA.QsigPnR <- unique(na.omit(wholeDF[netFAX.SA.FCsig & mrbpFAX.SA.QValsig   & protFAX.SA.QValsig, "UniprotACC"]))
netUV.DTT.QsigPnR <- unique(na.omit(wholeDF[netUV.DTT.FCsig  & mrbpUV.DTT.QValsig   & protUV.DTT.QValsig, "UniprotACC"]))
netUV.H2O2.QsigPnR <- unique(na.omit(wholeDF[netUV.H2O2.FCsig  & mrbpUV.H2O2.QValsig  & protUV.H2O2.QValsig, "UniprotACC"]))
netUV.SA.QsigPnR <- unique(na.omit(wholeDF[netUV.SA.FCsig  & mrbpUV.SA.QValsig    & protUV.SA.QValsig, "UniprotACC"]))
netFAX.DTT.QsigP <- unique(na.omit(wholeDF[netFAX.DTT.FCsig & protFAX.DTT.QValsig, "UniprotACC"]))
netFAX.H2O2.QsigP <- unique(na.omit(wholeDF[netFAX.H2O2.FCsig & protFAX.H2O2.QValsig, "UniprotACC"]))
netFAX.SA.QsigP <- unique(na.omit(wholeDF[netFAX.SA.FCsig & protFAX.SA.QValsig, "UniprotACC"]))
netUV.DTT.QsigP <- unique(na.omit(wholeDF[netUV.DTT.FCsig  & protUV.DTT.QValsig, "UniprotACC"]))
netUV.H2O2.QsigP <- unique(na.omit(wholeDF[netUV.H2O2.FCsig  & protUV.H2O2.QValsig, "UniprotACC"]))
netUV.SA.QsigP <- unique(na.omit(wholeDF[netUV.SA.FCsig  & protUV.SA.QValsig, "UniprotACC"]))
netFAX.DTT.QsigR <- unique(na.omit(wholeDF[netFAX.DTT.FCsig & mrbpFAX.DTT.QValsig, "UniprotACC"]))
netFAX.H2O2.QsigR <- unique(na.omit(wholeDF[netFAX.H2O2.FCsig & mrbpFAX.H2O2.QValsig, "UniprotACC"]))
netFAX.SA.QsigR <- unique(na.omit(wholeDF[netFAX.SA.FCsig & mrbpFAX.SA.QValsig, "UniprotACC"]))
netUV.DTT.QsigR <- unique(na.omit(wholeDF[netUV.DTT.FCsig  & mrbpUV.DTT.QValsig, "UniprotACC"]))
netUV.H2O2.QsigR <- unique(na.omit(wholeDF[netUV.H2O2.FCsig  & mrbpUV.H2O2.QValsig, "UniprotACC"]))
netUV.SA.QsigR <- unique(na.omit(wholeDF[netUV.SA.FCsig  & mrbpUV.SA.QValsig, "UniprotACC"]))

outdir <- "/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/data/geneLists/GSA"
saveTablesTsvExc(protFAX.DTT.Psig,outdir,excel=F)
saveTablesTsvExc(protFAX.H2O2.Psig,outdir,excel=F)
saveTablesTsvExc(protFAX.SA.Psig,outdir,excel=F)
saveTablesTsvExc(mrbpFAX.DTT.Psig,outdir,excel=F)
saveTablesTsvExc(mrbpFAX.H2O2.Psig,outdir,excel=F)
saveTablesTsvExc(mrbpFAX.SA.Psig,outdir,excel=F)
saveTablesTsvExc(protUV.DTT.Psig,outdir,excel=F)
saveTablesTsvExc(protUV.H2O2.Psig,outdir,excel=F)
saveTablesTsvExc(protUV.SA.Psig,outdir,excel=F)
saveTablesTsvExc(mrbpUV.DTT.Psig,outdir,excel=F)
saveTablesTsvExc(mrbpUV.H2O2.Psig,outdir,excel=F)
saveTablesTsvExc(mrbpUV.SA.Psig,outdir,excel=F)

saveTablesTsvExc(netFAX.DTT.PsigPnR,outdir,excel=F)
saveTablesTsvExc(netFAX.H2O2.PsigPnR,outdir,excel=F)
saveTablesTsvExc(netFAX.SA.PsigPnR,outdir,excel=F)
saveTablesTsvExc(netUV.DTT.PsigPnR,outdir,excel=F)
saveTablesTsvExc(netUV.H2O2.PsigPnR,outdir,excel=F)
saveTablesTsvExc(netUV.SA.PsigPnR,outdir,excel=F)
saveTablesTsvExc(netFAX.DTT.PsigP,outdir,excel=F)
saveTablesTsvExc(netFAX.H2O2.PsigP,outdir,excel=F)
saveTablesTsvExc(netFAX.SA.PsigP,outdir,excel=F)
saveTablesTsvExc(netUV.DTT.PsigP,outdir,excel=F)
saveTablesTsvExc(netUV.H2O2.PsigP,outdir,excel=F)
saveTablesTsvExc(netUV.SA.PsigP,outdir,excel=F)
saveTablesTsvExc(netFAX.DTT.PsigR,outdir,excel=F)
saveTablesTsvExc(netFAX.H2O2.PsigR,outdir,excel=F)
saveTablesTsvExc(netFAX.SA.PsigR,outdir,excel=F)
saveTablesTsvExc(netUV.DTT.PsigR,outdir,excel=F)
saveTablesTsvExc(netUV.H2O2.PsigR,outdir,excel=F)
saveTablesTsvExc(netUV.SA.PsigR,outdir,excel=F)

saveTablesTsvExc(protFAX.DTT.Qsig,outdir,excel=F)
saveTablesTsvExc(protFAX.H2O2.Qsig,outdir,excel=F)
saveTablesTsvExc(protFAX.SA.Qsig,outdir,excel=F)
saveTablesTsvExc(mrbpFAX.DTT.Qsig,outdir,excel=F)
saveTablesTsvExc(mrbpFAX.H2O2.Qsig,outdir,excel=F)
saveTablesTsvExc(mrbpFAX.SA.Qsig,outdir,excel=F)
saveTablesTsvExc(protUV.DTT.Qsig,outdir,excel=F)
saveTablesTsvExc(protUV.H2O2.Qsig,outdir,excel=F)
saveTablesTsvExc(protUV.SA.Qsig,outdir,excel=F)
saveTablesTsvExc(mrbpUV.DTT.Qsig,outdir,excel=F)
saveTablesTsvExc(mrbpUV.H2O2.Qsig,outdir,excel=F)
saveTablesTsvExc(mrbpUV.SA.Qsig,outdir,excel=F)
saveTablesTsvExc(netFAX.DTT.QsigPnR,outdir,excel=F)
saveTablesTsvExc(netFAX.H2O2.QsigPnR,outdir,excel=F)
saveTablesTsvExc(netFAX.SA.QsigPnR,outdir,excel=F)
saveTablesTsvExc(netUV.DTT.QsigPnR,outdir,excel=F)
saveTablesTsvExc(netUV.H2O2.QsigPnR,outdir,excel=F)
saveTablesTsvExc(netUV.SA.QsigPnR,outdir,excel=F)
saveTablesTsvExc(netFAX.DTT.QsigP,outdir,excel=F)
saveTablesTsvExc(netFAX.H2O2.QsigP,outdir,excel=F)
saveTablesTsvExc(netFAX.SA.QsigP,outdir,excel=F)
saveTablesTsvExc(netUV.DTT.QsigP,outdir,excel=F)
saveTablesTsvExc(netUV.H2O2.QsigP,outdir,excel=F)
saveTablesTsvExc(netUV.SA.QsigP,outdir,excel=F)
saveTablesTsvExc(netFAX.DTT.QsigR,outdir,excel=F)
saveTablesTsvExc(netFAX.H2O2.QsigR,outdir,excel=F)
saveTablesTsvExc(netFAX.SA.QsigR,outdir,excel=F)
saveTablesTsvExc(netUV.DTT.QsigR,outdir,excel=F)
saveTablesTsvExc(netUV.H2O2.QsigR,outdir,excel=F)
saveTablesTsvExc(netUV.SA.QsigR,outdir,excel=F)

allSigDFs <- list(protFAX.DTT.Psig, protFAX.H2O2.Psig, protFAX.SA.Psig, mrbpFAX.DTT.Psig, mrbpFAX.H2O2.Psig, 
                  mrbpFAX.SA.Psig, protUV.DTT.Psig, protUV.H2O2.Psig, protUV.SA.Psig, mrbpUV.DTT.Psig, 
                  mrbpUV.H2O2.Psig, mrbpUV.SA.Psig, netFAX.DTT.PsigPnR, netFAX.H2O2.PsigPnR, netFAX.SA.PsigPnR, 
                  netUV.DTT.PsigPnR, netUV.H2O2.PsigPnR, netUV.SA.PsigPnR, netFAX.DTT.PsigP, netFAX.H2O2.PsigP, 
                  netFAX.SA.PsigP, netUV.DTT.PsigP, netUV.H2O2.PsigP, netUV.SA.PsigP, netFAX.DTT.PsigR, netFAX.H2O2.PsigR, 
                  netFAX.SA.PsigR, netUV.DTT.PsigR, netUV.H2O2.PsigR, netUV.SA.PsigR, protFAX.DTT.Qsig, protFAX.H2O2.Qsig, 
                  protFAX.SA.Qsig, mrbpFAX.DTT.Qsig, mrbpFAX.H2O2.Qsig, mrbpFAX.SA.Qsig, protUV.DTT.Qsig, protUV.H2O2.Qsig, 
                  protUV.SA.Qsig, mrbpUV.DTT.Qsig, mrbpUV.H2O2.Qsig, mrbpUV.SA.Qsig, netFAX.DTT.QsigPnR, netFAX.H2O2.QsigPnR, 
                  netFAX.SA.QsigPnR, netUV.DTT.QsigPnR, netUV.H2O2.QsigPnR, netUV.SA.QsigPnR, netFAX.DTT.QsigP, 
                  netFAX.H2O2.QsigP, netFAX.SA.QsigP, netUV.DTT.QsigP, netUV.H2O2.QsigP, netUV.SA.QsigP, netFAX.DTT.QsigR, 
                  netFAX.H2O2.QsigR, netFAX.SA.QsigR, netUV.DTT.QsigR, netUV.H2O2.QsigR, netUV.SA.QsigR)

countSigs <- unlist(lapply(allSigDFs, function(x){length(x)}))

countSigsNames <- deparse(substitute(c(protFAX.DTT.Psig, protFAX.H2O2.Psig, protFAX.SA.Psig, mrbpFAX.DTT.Psig, mrbpFAX.H2O2.Psig, 
                                       mrbpFAX.SA.Psig, protUV.DTT.Psig, protUV.H2O2.Psig, protUV.SA.Psig, mrbpUV.DTT.Psig, 
                                       mrbpUV.H2O2.Psig, mrbpUV.SA.Psig, netFAX.DTT.PsigPnR, netFAX.H2O2.PsigPnR, netFAX.SA.PsigPnR, 
                                       netUV.DTT.PsigPnR, netUV.H2O2.PsigPnR, netUV.SA.PsigPnR, netFAX.DTT.PsigP, netFAX.H2O2.PsigP, 
                                       netFAX.SA.PsigP, netUV.DTT.PsigP, netUV.H2O2.PsigP, netUV.SA.PsigP, netFAX.DTT.PsigR, netFAX.H2O2.PsigR, 
                                       netFAX.SA.PsigR, netUV.DTT.PsigR, netUV.H2O2.PsigR, netUV.SA.PsigR, protFAX.DTT.Qsig, protFAX.H2O2.Qsig, 
                                       protFAX.SA.Qsig, mrbpFAX.DTT.Qsig, mrbpFAX.H2O2.Qsig, mrbpFAX.SA.Qsig, protUV.DTT.Qsig, protUV.H2O2.Qsig, 
                                       protUV.SA.Qsig, mrbpUV.DTT.Qsig, mrbpUV.H2O2.Qsig, mrbpUV.SA.Qsig, netFAX.DTT.QsigPnR, netFAX.H2O2.QsigPnR, 
                                       netFAX.SA.QsigPnR, netUV.DTT.QsigPnR, netUV.H2O2.QsigPnR, netUV.SA.QsigPnR, netFAX.DTT.QsigP, 
                                       netFAX.H2O2.QsigP, netFAX.SA.QsigP, netUV.DTT.QsigP, netUV.H2O2.QsigP, netUV.SA.QsigP, netFAX.DTT.QsigR, 
                                       netFAX.H2O2.QsigR, netFAX.SA.QsigR, netUV.DTT.QsigR, netUV.H2O2.QsigR, netUV.SA.QsigR)))

countSigsNames <- gsub(" ","",unlist(strsplit(gsub(")","",gsub("c(","",countSigsNames,fixed = T),fixed = T),", ")))
protNrbpSigs <- as.data.frame(cbind(countSigsNames,countSigs))
colnames(protNrbpSigs) <- c("data","n.Psig.prots")
outdir <- "/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/data/geneLists"
saveTablesTsvExc(protNrbpSigs,outdir,excel=F,bycompleteFC = F)

#########################
#### SORTED FOR GSEA ####
#########################
proteomeFAXfcDTT <- wholeDF[,c(IDs,"proteomeFAX.log2ratio_Condition2")]
proteomeFAXfcH2O2 <- wholeDF[,c(IDs,"proteomeFAX.log2ratio_Condition3")]
proteomeFAXfcSA <- wholeDF[,c(IDs,"proteomeFAX.log2ratio_Condition4")]
mRNABPomeFAXfcDTT <- wholeDF[,c(IDs,"mRNABPomeFAX.log2ratio_PolyARNAFAXwithDTT")]
mRNABPomeFAXfcH2O2 <- wholeDF[,c(IDs,"mRNABPomeFAX.log2ratio_PolyARNAFAXwithH2O2")]
mRNABPomeFAXfcSA <- wholeDF[,c(IDs,"mRNABPomeFAX.log2ratio_PolyARNAFAXwithSA")]
netChangesFAXDTT <- wholeDF[,c(IDs,"netChangesFAXDTT")]
netChangesFAXH2O2 <- wholeDF[,c(IDs,"netChangesFAXH2O2")]
netChangesFAXSA <- wholeDF[,c(IDs,"netChangesFAXSA")]

proteomeUVfcDTT <- wholeDF[,c(IDs,"proteomeUV.log2ratio_Condition2")]
proteomeUVfcH2O2 <- wholeDF[,c(IDs,"proteomeUV.log2ratio_Condition3")]
proteomeUVfcSA <- wholeDF[,c(IDs,"proteomeUV.log2ratio_Condition4")]
mRNABPomeUVfcDTT <- wholeDF[,c(IDs,"mRNABPomeUV.log2ratio_PolyARNAUVwithDTT")]
mRNABPomeUVfcH2O2 <- wholeDF[,c(IDs,"mRNABPomeUV.log2ratio_PolyARNAUVwithH202")]
mRNABPomeUVfcSA <- wholeDF[,c(IDs,"mRNABPomeUV.log2ratio_PolyARNAUVwithSA")]
netChangesUVDTT <- wholeDF[,c(IDs,"netChangesUVDTT")]
netChangesUVH2O2 <- wholeDF[,c(IDs,"netChangesUVH2O")]
netChangesUVSA <- wholeDF[,c(IDs,"netChangesUVSA")]

outdir <- "/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/data/geneLists/GSEA"
saveTablesTsvExc(proteomeFAXfcDTT,outdir,excel=F)
saveTablesTsvExc(proteomeFAXfcH2O2,outdir,excel=F)
saveTablesTsvExc(proteomeFAXfcSA,outdir,excel=F)
saveTablesTsvExc(mRNABPomeFAXfcDTT,outdir,excel=F)
saveTablesTsvExc(mRNABPomeFAXfcH2O2,outdir,excel=F)
saveTablesTsvExc(mRNABPomeFAXfcSA,outdir,excel=F)
saveTablesTsvExc(netChangesFAXDTT,outdir,excel=F)
saveTablesTsvExc(netChangesFAXH2O2,outdir,excel=F)
saveTablesTsvExc(netChangesFAXSA,outdir,excel=F)
saveTablesTsvExc(proteomeUVfcDTT,outdir,excel=F)
saveTablesTsvExc(proteomeUVfcH2O2,outdir,excel=F)
saveTablesTsvExc(proteomeUVfcSA,outdir,excel=F)
saveTablesTsvExc(mRNABPomeUVfcDTT,outdir,excel=F)
saveTablesTsvExc(mRNABPomeUVfcH2O2,outdir,excel=F)
saveTablesTsvExc(mRNABPomeUVfcSA,outdir,excel=F)
saveTablesTsvExc(netChangesUVDTT,outdir,excel=F)
saveTablesTsvExc(netChangesUVH2O2,outdir,excel=F)
saveTablesTsvExc(netChangesUVSA,outdir,excel=F)



