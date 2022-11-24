ProteomeFAXfile <- "/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/data/raw/ProteomeFAX_np2_norm_gmin_no15-25-26-36-39/ProteomeFAX_np2_norm_gmin_no15-25-26-36-39_PROTEIN.tsv"
ProteomeUVXfile <- "/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/data/raw/ProteomeUVX_np2_norm_gmin_no4-6-23-31-32/ProteomeUVX_np2_norm_gmin_no4-6-23-31-32_PROTEIN.tsv"
RBPomeFAXfile <-  "/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/data/raw/RBPomeFAX_np2_norm_gmin_no12/RBPomeFAX_np2_norm_gmin_no12_PROTEIN.tsv"
RBPomeUVXfile <- "/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/data/raw/RBPomeUVX_np2_norm_gmin_no10-7-16-3/RBPomeUVX_np2_norm_gmin_no10-7-16-3_PROTEIN.tsv"

ProteomeFAX <- read.delim(ProteomeFAXfile)
ProteomeUVX <- read.delim(ProteomeUVXfile)
RBPomeFAX <- read.delim(RBPomeFAXfile)
RBPomeUVX <- read.delim(RBPomeUVXfile)
commonCols <- intersect(colnames(ProteomeFAX),colnames(ProteomeUVX))

colnames(ProteomeFAX) <- paste0("ProteomeFAX.",colnames(ProteomeFAX))
colnames(ProteomeUVX) <- paste0("ProteomeUVX.",colnames(ProteomeUVX))
colnames(RBPomeFAX) <- paste0("RBPomeFAX.",colnames(RBPomeFAX))
colnames(RBPomeUVX) <- paste0("RBPomeUVX.",colnames(RBPomeUVX))

colnames(ProteomeFAX)[1:length(commonCols)] <- commonCols
colnames(ProteomeUVX)[1:length(commonCols)] <- commonCols
colnames(RBPomeFAX)[1:length(commonCols)] <- commonCols
colnames(RBPomeUVX)[1:length(commonCols)] <- commonCols

# P05755 a paralogue that evidences the issue in raw data see repairOfRawProteome.R  
## Replace single orthologues annotation with both
ProteomeFAX$ac[grepl(";",ProteomeFAX$proteinName)] <- ProteomeFAX$proteinName[grepl(";",ProteomeFAX$proteinName)]
ProteomeUVX$ac[grepl(";",ProteomeUVX$proteinName)] <- ProteomeUVX$proteinName[grepl(";",ProteomeUVX$proteinName)]
RBPomeFAX$ac[grepl(";",RBPomeFAX$proteinName)] <- RBPomeFAX$proteinName[grepl(";",RBPomeFAX$proteinName)]
RBPomeUVX$ac[grepl(";",RBPomeUVX$proteinName)] <- RBPomeUVX$proteinName[grepl(";",RBPomeUVX$proteinName)]

ProteomeFAX$geneName[grepl(";",ProteomeFAX$proteinName)] <- ProteomeFAX$proteinName[grepl(";",ProteomeFAX$proteinName)]
ProteomeUVX$geneName[grepl(";",ProteomeUVX$proteinName)] <- ProteomeUVX$proteinName[grepl(";",ProteomeUVX$proteinName)]
RBPomeFAX$geneName[grepl(";",RBPomeFAX$proteinName)] <- RBPomeFAX$proteinName[grepl(";",RBPomeFAX$proteinName)]
RBPomeUVX$geneName[grepl(";",RBPomeUVX$proteinName)] <- RBPomeUVX$proteinName[grepl(";",RBPomeUVX$proteinName)]

# Merge Proteome and RBPome proteinName is the Accession in Raw and is the common column
FAXdf <- merge(ProteomeFAX,RBPomeFAX,all = T, by = "proteinName",suffixes = c(".ProteomeFAX",".RBPomeFAX")); nrow(FAXdf)
UVXdf <- merge(ProteomeUVX,RBPomeUVX,all = T, by = "proteinName",suffixes = c(".ProteomeUVX",".RBPomeUVX")); nrow(UVXdf)

FAXdf$FAXnetchangesDTT <- FAXdf$RBPomeFAX.log2ratio_PolyARNAFAXwithDTT - FAXdf$ProteomeFAX.log2ratio_FAXwithDTT
FAXdf$FAXnetchangesH2O2 <- FAXdf$RBPomeFAX.log2ratio_PolyARNAFAXwithH2O2 - FAXdf$ProteomeFAX.log2ratio_FAXwithH2O2
FAXdf$FAXnetchangesSA <- FAXdf$RBPomeFAX.log2ratio_PolyARNAFAXwithSA - FAXdf$ProteomeFAX.log2ratio_FAXwithSA

UVXdf$UVXnetchangesDTT <- UVXdf$RBPomeUVX.log2ratio_PolyARNAUVwithDTT - UVXdf$ProteomeUVX.log2ratio_Condition2
UVXdf$UVXnetchangesH2O2 <- UVXdf$RBPomeUVX.log2ratio_PolyARNAUVwithH202 - UVXdf$ProteomeUVX.log2ratio_Condition3
UVXdf$UVXnetchangesSA <- UVXdf$RBPomeUVX.log2ratio_PolyARNAUVwithSA - UVXdf$ProteomeUVX.log2ratio_Condition4

outdir <- "/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/FinalData"
saveTablesTsvExc(FAXdf,outdir,completeNdedup=F,excel=T,bycompleteFC=F,rownames=F)
saveTablesTsvExc(UVXdf,outdir,completeNdedup=F,excel=T,bycompleteFC=F,rownames=F)

