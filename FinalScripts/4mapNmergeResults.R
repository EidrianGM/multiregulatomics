library(openxlsx)
mappingFinalFile <- "yeastReference/mappingFile.xlsx"
yeastGenesProtMap <- read.xlsx(mappingFinalFile)
yeastGenesProtMap <- yeastGenesProtMap[!duplicated(yeastGenesProtMap),]

########################
##### MAPPING DEGs #####
########################

DEGsFile <- "data/RNA-Seq data/result/DE gene testing statistics.csv"
DEGsDF <- read.delim(DEGsFile, sep = ",") 

TPMsFile <- "data/RNA-Seq data/result/TPM_genes.csv"
TPMsDF <- read.delim(TPMsFile, sep = ",") 

DEGsDF <- reshape(DEGsDF, direction = "wide", idvar = "target", timevar = "contrast")

DEGsDFfull <- merge(TPMsDF,DEGsDF,by.x="X",by.y="target")
colnames(DEGsDFfull)[1] <- "target"

DEGsDFMaped <- merge(DEGsDFfull,yeastGenesProtMap,by.x="target",by.y="ORFgene")
DEGsDFMaped$ORFgene <- DEGsDFMaped$target
any(duplicated(DEGsDFMaped))
DEGsDFMaped <- DEGsDFMaped[!duplicated(DEGsDFMaped),]
DEGsNotMapped <- DEGsDFfull[!(DEGsDFfull$target %in% DEGsDFMaped$target),]

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

colnames(DEGsDFMapedFull)[colnames(DEGsDFMapedFull) %in% colnames(yeastGenesProtMap)]
colnames(DEGsDFMapedFull)[!colnames(DEGsDFMapedFull) %in% colnames(yeastGenesProtMap)] <- paste0("DEGs.",colnames(DEGsDFMapedFull)[!colnames(DEGsDFMapedFull) %in% colnames(yeastGenesProtMap)])

DEGsDFMapedFull <- DEGsDFMapedFull[,c(colnames(DEGsDFMapedFull)[colnames(DEGsDFMapedFull) %in% colnames(yeastGenesProtMap)], colnames(DEGsDFMapedFull)[!colnames(DEGsDFMapedFull) %in% colnames(yeastGenesProtMap)])]

DEGsDFMapedFullFile <- '/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/FinalData/DEGsFullMapping.tsv'
write.table(DEGsDFMapedFull,DEGsDFMapedFullFile,quote = F,sep = '\t',col.names = T,row.names = F)


################################################################################
########################## Merging Proteomes and RBPomes #######################
################################################################################

ProteomeFAXfile <- "/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/FinalData/Raw/correctedRAW/ProteomeFAX_np2_norm_gmin_no15-25-26-36-39/ProteomeFAX_np2_norm_gmin_no15-25-26-36-39_PROTEIN.tsv"
ProteomeUVXfile <- "/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/FinalData/Raw/correctedRAW/ProteomeUVX_np2_norm_gmin_no4-6-23-31-32/ProteomeUVX_np2_norm_gmin_no4-6-23-31-32_PROTEIN.tsv"
RBPomeFAXfile <-  "/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/FinalData/Raw/correctedRAW/RBPomeFAX_np2_norm_gmin_no12/RBPomeFAX_np2_norm_gmin_no12_PROTEIN.tsv"
RBPomeUVXfile <- "/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/FinalData/Raw/correctedRAW/RBPomeUVX_np2_norm_gmin_no10-7-16-3/RBPomeUVX_np2_norm_gmin_no10-7-16-3_PROTEIN.tsv"

ProteomeFAX <- read.delim(ProteomeFAXfile); nrow(ProteomeFAX)
ProteomeUVX <- read.delim(ProteomeUVXfile); nrow(ProteomeUVX)
RBPomeFAX <- read.delim(RBPomeFAXfile); nrow(RBPomeFAX)
RBPomeUVX <- read.delim(RBPomeUVXfile); nrow(RBPomeUVX)

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
## Replace orthologues acc and geneName annotation with both UniProts
ProteomeFAX$ac[grepl(";",ProteomeFAX$proteinName)] <- ProteomeFAX$proteinName[grepl(";",ProteomeFAX$proteinName)]
ProteomeUVX$ac[grepl(";",ProteomeUVX$proteinName)] <- ProteomeUVX$proteinName[grepl(";",ProteomeUVX$proteinName)]
RBPomeFAX$ac[grepl(";",RBPomeFAX$proteinName)] <- RBPomeFAX$proteinName[grepl(";",RBPomeFAX$proteinName)]
RBPomeUVX$ac[grepl(";",RBPomeUVX$proteinName)] <- RBPomeUVX$proteinName[grepl(";",RBPomeUVX$proteinName)]

ProteomeFAX$geneName[grepl(";",ProteomeFAX$proteinName)] <- ""
ProteomeUVX$geneName[grepl(";",ProteomeUVX$proteinName)] <- ""
RBPomeFAX$geneName[grepl(";",RBPomeFAX$proteinName)] <- ""
RBPomeUVX$geneName[grepl(";",RBPomeUVX$proteinName)] <- ""

# Merge Proteome and RBPome proteinName is the Accession in Raw and is the common column
FAXdf <- merge(ProteomeFAX,RBPomeFAX,all = T, by = "proteinName",suffixes = c(".ProteomeFAX",".RBPomeFAX")); nrow(ProteomeFAX); nrow(FAXdf)
UVXdf <- merge(ProteomeUVX,RBPomeUVX,all = T, by = "proteinName",suffixes = c(".ProteomeUVX",".RBPomeUVX")); nrow(ProteomeUVX); nrow(UVXdf)

toVennList <- list(ProteomeFAX = ProteomeFAX$proteinName, RBPomeFAX = RBPomeFAX$proteinName)
ggVennDiagram(toVennList, color = 2, lwd = 0.7) + scale_fill_gradient(low = "#F4FAFE", high = "#4981BF") + theme(legend.position = "none")

if (nrow(FAXdf) == length(unique(c(RBPomeFAX$proteinName, ProteomeFAX$proteinName)))){print("GOOD MERGE")}
if (nrow(UVXdf) == length(unique(c(RBPomeUVX$proteinName, ProteomeUVX$proteinName)))){print("GOOD MERGE")}

FAXdf$FAXnetchangesDTT <- FAXdf$RBPomeFAX.log2ratio_PolyARNAFAXwithDTT - FAXdf$ProteomeFAX.log2ratio_FAXwithDTT
FAXdf$FAXnetchangesH2O2 <- FAXdf$RBPomeFAX.log2ratio_PolyARNAFAXwithH2O2 - FAXdf$ProteomeFAX.log2ratio_FAXwithH2O2
FAXdf$FAXnetchangesSA <- FAXdf$RBPomeFAX.log2ratio_PolyARNAFAXwithSA - FAXdf$ProteomeFAX.log2ratio_FAXwithSA

UVXdf$UVXnetchangesDTT <- UVXdf$RBPomeUVX.log2ratio_PolyARNAUVwithDTT - UVXdf$ProteomeUVX.log2ratio_Condition2
UVXdf$UVXnetchangesH2O2 <- UVXdf$RBPomeUVX.log2ratio_PolyARNAUVwithH202 - UVXdf$ProteomeUVX.log2ratio_Condition3
UVXdf$UVXnetchangesSA <- UVXdf$RBPomeUVX.log2ratio_PolyARNAUVwithSA - UVXdf$ProteomeUVX.log2ratio_Condition4

outdir <- "/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/FinalData"
saveTablesTsvExc(FAXdf,outdir,completeNdedup=F,excel=T,bycompleteFC=F,rownames=F)
saveTablesTsvExc(UVXdf,outdir,completeNdedup=F,excel=T,bycompleteFC=F,rownames=F)

################################################################################
########################## MERGING FAX and UVX #################################
################################################################################

FAXnUVXdf <- merge(FAXdf,UVXdf,all = T, by = "proteinName"); nrow(FAXdf); nrow(FAXnUVXdf)

toVennList <- list(FAX = FAXdf$proteinName, UVX = UVXdf$proteinName)
ggVennDiagram(toVennList, color = 2, lwd = 0.7) + scale_fill_gradient(low = "#F4FAFE", high = "#4981BF") + theme(legend.position = "none")

outdir <- "/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/FinalData"
saveTablesTsvExc(FAXnUVXdf,outdir,completeNdedup=F,excel=T,bycompleteFC=F,rownames=F)

################################################################################
########################## MAPPING FAX and UVX #################################
################################################################################

# 1 UVX
FAXnUVXFile <- "/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/FinalData/FAXnUVXdf.tsv"
FAXnUVXdf <- read.delim(FAXnUVXFile)

FAXnUVXOrthologs <- FAXnUVXdf[grepl(";",FAXnUVXdf$proteinName),]; nrow(FAXnUVXOrthologs)
FAXnUVXSingProts <- FAXnUVXdf[!grepl(";",FAXnUVXdf$proteinName),]; nrow(FAXnUVXSingProts)

toVennList <- list(mapping = yeastGenesProtMap$UniprotACC, FAXnUVX = FAXnUVXSingProts$proteinName)
ggVennDiagram(toVennList, color = 2, lwd = 0.7) + scale_fill_gradient(low = "#F4FAFE", high = "#4981BF") + theme(legend.position = "none")

FAXnUVXmapped <- merge(yeastGenesProtMap,FAXnUVXSingProts,by.x = "UniprotACC", by.y = "proteinName")
FAXnUVXmapped$proteinName <- FAXnUVXmapped$UniprotACC
nrow(FAXnUVXSingProts); nrow(FAXnUVXmapped)

FAXnUVXmapped$UniprotACC[duplicated(FAXnUVXmapped$UniprotACC)] 
View(FAXnUVXmapped[FAXnUVXmapped$UniprotACC %in% FAXnUVXmapped$UniprotACC[duplicated(FAXnUVXmapped$UniprotACC)],])
# These proteins are located in two different chromosomes of Yeast and are behind a duplicity

FAXnUVXunmapped <- FAXnUVXSingProts[!(FAXnUVXSingProts$proteinName %in% yeastGenesProtMap$UniprotACC),]

toVennList <- list(mapping = yeastGenesProtMap$GeneName, FAXnUVX = FAXnUVXunmapped$geneName.ProteomeFAX)
ggVennDiagram(toVennList, color = 2, lwd = 0.7) + scale_fill_gradient(low = "#F4FAFE", high = "#4981BF") + theme(legend.position = "none")

FAXnUVXmapped2 <- merge(yeastGenesProtMap,FAXnUVXunmapped,by.x="GeneName",by.y="geneName.ProteomeFAX")
FAXnUVXmapped2$geneName.ProteomeFAX <- FAXnUVXmapped2$GeneName

FAXnUVXunmapped <- FAXnUVXunmapped[!(FAXnUVXunmapped$geneName.ProteomeFAX %in% yeastGenesProtMap$GeneName),]
nrow(FAXnUVXunmapped)

FAXnUVXfullMapped <- rbind.fill(FAXnUVXmapped,FAXnUVXmapped2,FAXnUVXOrthologs)
FAXnUVXfullMapped <- FAXnUVXfullMapped[,c(colnames(yeastGenesProtMap), colnames(FAXnUVXfullMapped)[!(colnames(FAXnUVXfullMapped) %in% colnames(yeastGenesProtMap))])]

outdir <- "/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/FinalData"
saveTablesTsvExc(FAXnUVXfullMapped,outdir,completeNdedup=F,excel=T,bycompleteFC=F,rownames=F)

################################################################################
###################### MERGIN Transcriptome and FAXnUVX ########################
################################################################################

toVennList <- list(DEGs = DEGsDFMapedFull$UniprotACC, FAXnUVX = FAXnUVXfullMapped$UniprotACC)
ggVennDiagram(toVennList, color = 2, lwd = 0.7) + scale_fill_gradient(low = "#F4FAFE", high = "#4981BF") + theme(legend.position = "none")
allDataDF <- merge(DEGsDFMapedFull,FAXnUVXfullMapped,by = "UniprotACC",all = T)
outdir <- "/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/FinalData"
saveTablesTsvExc(allDataDF,outdir,completeNdedup=F,excel=T,bycompleteFC=F,rownames=F)

checkingMerging <- cbind(allDataDF$proteinName[!is.na(allDataDF$proteinName)],allDataDF$UniprotACC[!is.na(allDataDF$proteinName)])
checkingMerging <- checkingMerging[complete.cases(checkingMerging),]
if (identical(checkingMerging[,1],checkingMerging[,2])){
  print("Merge correct")
}

### FINAL TABLES FOR IBTISSAM

FAXacceptedProts <- read.delim("FinalData/BackgroundRemoval/FAXAccBackGrAcceptedFC3.tsv",header = F)[,1]
UVXacceptedProts <- read.delim("FinalData/BackgroundRemoval/UVXAccBackGrAcceptedFC3.tsv",header = F)[,1]

FAXwoBKGR <- wholeDF[which((!is.na(wholeDF$RBPomeFAX.log2ratio_PolyARNAFAXwithSA)) & wholeDF$proteinName %in% FAXacceptedProts),]
length(FAXwoBKGR$proteinName); length(unique(FAXwoBKGR$proteinName))

# The proteins P10081 P02309 P32324 P02994 are located in two different chromosomal position this means an extra 4 rows due to mapping
FAXwoBKGR$proteinName[duplicated(FAXwoBKGR$proteinName)]

View(
  FAXwoBKGR[FAXwoBKGR$proteinName %in% FAXwoBKGR$proteinName[duplicated(FAXwoBKGR$proteinName)],]
)
UVXwoBKGR <- wholeDF[which((!is.na(wholeDF$RBPomeUVX.cv_PolyARNAUVwithSA)) & wholeDF$proteinName %in% UVXacceptedProts),]
length(UVXwoBKGR$proteinName); length(unique(UVXwoBKGR$proteinName))
View(
  UVXwoBKGR[UVXwoBKGR$proteinName %in% UVXwoBKGR$proteinName[duplicated(UVXwoBKGR$proteinName)],]
)
# Like in FAX the proteins P10081 P02309 are in two locations, so 2 extra rows.
UVXwoBKGR$proteinName[duplicated(UVXwoBKGR$proteinName)]

outdir <- 'FinalData/'
saveTablesTsvExc(FAXwoBKGR,outdir,completeNdedup=F,excel=T,bycompleteFC=F,rownames=F)
saveTablesTsvExc(UVXwoBKGR,outdir,completeNdedup=F,excel=T,bycompleteFC=F,rownames=F)


# FAXwoBKGRfile <- "FinalData/FAXwoBKGR.tsv"
# UVXwoBKGRfile <- "FinalData/UVXwoBKGR.tsv"
# FAXwoBKGR <- read.delim(FAXwoBKGRfile)
# UVXwoBKGR <- read.delim(UVXwoBKGRfile)

cat(paste(colnames(FAXwoBKGR),collapse = "\n"))


columnsSelected <-  c('seqnames','start','end','type','ORF','SGDID','GeneName','proteinName','UniprotACC',
                      'UniprotName','Gene.Names','Protein.names','Gene.Names..ordered.locus.',
                      'Gene.Names..ORF.','Gene.Names..primary.','Gene.Names..synonym.','Alias','ORFgene',
                      'DEGs.H202.brep1','DEGs.H202.brep2','DEGs.H202.brep3','DEGs.DTT.brep1','DEGs.DTT.brep2',
                      'DEGs.DTT.brep3','DEGs.SA.brep1','DEGs.SA.brep2','DEGs.SA.brep3','DEGs.YPD.brep1','DEGs.YPD.brep2',
                      'DEGs.YPD.brep3','DEGs.adj.pval.H202.YPD','DEGs.log2FC.H202.YPD','DEGs.adj.pval.DTT.YPD',
                      'DEGs.log2FC.DTT.YPD','DEGs.adj.pval.SA.YPD','DEGs.log2FC.SA.YPD','ProteomeFAX.X014_Pr40.DIA_DIA_D22',
                      'ProteomeFAX.X011_Pr40.DIA_DIA_D22','ProteomeFAX.X013_Pr40.DIA_DIA_D22',
                      'ProteomeFAX.X033_Pr40.DIA_DIA_D22','ProteomeFAX.X040_Pr40.DIA_DIA_D22',
                      'ProteomeFAX.X019_Pr40.DIA_DIA_D22','ProteomeFAX.X020_Pr40.DIA_DIA_D22',
                      'ProteomeFAX.X021_Pr40.DIA_DIA_D22','ProteomeFAX.log2ratio_FAXwithDTT','ProteomeFAX.log2ratio_FAXwithH2O2',
                      'ProteomeFAX.log2ratio_FAXwithSA','ProteomeFAX.pValue_FAXwithDTT','ProteomeFAX.pValue_FAXwithH2O2',
                      'ProteomeFAX.pValue_FAXwithSA','ProteomeFAX.qValue_FAXwithDTT','ProteomeFAX.qValue_FAXwithH2O2',
                      'ProteomeFAX.qValue_FAXwithSA','RBPomeFAX.X026_PolyA.40_APMS_F21','RBPomeFAX.X023_PolyA.40_APMS_F21',
                      'RBPomeFAX.X025_PolyA.40_APMS_F21','RBPomeFAX.X027_PolyA.40_APMS_F21','RBPomeFAX.X005_PolyA.40_APMS_F21',
                      'RBPomeFAX.X038_PolyA.40_APMS_F21','RBPomeFAX.X031_PolyA.40_APMS_F21','RBPomeFAX.X032_PolyA.40_APMS_F21',
                      'RBPomeFAX.X033_PolyA.40_APMS_F21','RBPomeFAX.log2ratio_PolyARNAFAXwithDTT',
                      'RBPomeFAX.log2ratio_PolyARNAFAXwithH2O2','RBPomeFAX.log2ratio_PolyARNAFAXwithSA',
                      'RBPomeFAX.pValue_PolyARNAFAXwithDTT','RBPomeFAX.pValue_PolyARNAFAXwithH2O2',
                      'RBPomeFAX.pValue_PolyARNAFAXwithSA','RBPomeFAX.qValue_PolyARNAFAXwithDTT','RBPomeFAX.qValue_PolyARNAFAXwithH2O2',
                      'RBPomeFAX.qValue_PolyARNAFAXwithSA','FAXnetchangesDTT','FAXnetchangesH2O2','FAXnetchangesSA',
                      'ProteomeUVX.X038_Pr40.DIA_DIA_D22','ProteomeUVX.X005_Pr40.DIA_DIA_D22','ProteomeUVX.X007_Pr40.DIA_DIA_D22',
                      'ProteomeUVX.X009_Pr40.DIA_DIA_D22','ProteomeUVX.X028_Pr40.DIA_DIA_D22','ProteomeUVX.X035_Pr40.DIA_DIA_D22',
                      'ProteomeUVX.X022_Pr40.DIA_DIA_D22','ProteomeUVX.X024_Pr40.DIA_DIA_D22','ProteomeUVX.log2ratio_Condition2',
                      'ProteomeUVX.log2ratio_Condition3','ProteomeUVX.log2ratio_Condition4','ProteomeUVX.pValue_Condition2',
                      'ProteomeUVX.pValue_Condition3','ProteomeUVX.pValue_Condition4','ProteomeUVX.qValue_Condition2',
                      'ProteomeUVX.qValue_Condition3','ProteomeUVX.qValue_Condition4','RBPomeUVX.X039_PolyA.40_APMS_F21',
                      'RBPomeUVX.X017_PolyA.40_APMS_F21','RBPomeUVX.X019_PolyA.40_APMS_F21','RBPomeUVX.X021_PolyA.40_APMS_F21',
                      'RBPomeUVX.X004_PolyA.40_APMS_F21','RBPomeUVX.X040_PolyA.40_APMS_F21','RBPomeUVX.X034_PolyA.40_APMS_F21',
                      'RBPomeUVX.X035_PolyA.40_APMS_F21','RBPomeUVX.X036_PolyA.40_APMS_F21','RBPomeUVX.log2ratio_PolyARNAUVwithDTT',
                      'RBPomeUVX.log2ratio_PolyARNAUVwithH202','RBPomeUVX.log2ratio_PolyARNAUVwithSA','RBPomeUVX.pValue_PolyARNAUVwithDTT',
                      'RBPomeUVX.pValue_PolyARNAUVwithH202','RBPomeUVX.pValue_PolyARNAUVwithSA','RBPomeUVX.qValue_PolyARNAUVwithDTT',
                      'RBPomeUVX.qValue_PolyARNAUVwithH202','RBPomeUVX.qValue_PolyARNAUVwithSA','UVXnetchangesDTT','UVXnetchangesH2O2','UVXnetchangesSA')


FAXwoBKGR <- FAXwoBKGR[,columnsSelected]
UVXwoBKGR <- UVXwoBKGR[,columnsSelected]
outdir <- "FinalData/TablesForPublication"
saveTablesTsvExc(FAXwoBKGR,outdir,completeNdedup=F,excel=T,bycompleteFC=F,rownames=F)
saveTablesTsvExc(UVXwoBKGR,outdir,completeNdedup=F,excel=T,bycompleteFC=F,rownames=F)


