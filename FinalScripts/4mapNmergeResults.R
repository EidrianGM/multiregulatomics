library(openxlsx)
library(ggVennDiagram)
library(ggplot2)

mappingFinalFile <- "yeastReference/mappingFile.xlsx"
yeastGenesProtMap <- read.xlsx(mappingFinalFile)
yeastGenesProtMap <- yeastGenesProtMap[!duplicated(yeastGenesProtMap),]

########################
##### MAPPING DEGs #####
########################

DEGsFile <- "data/RNA-Seq data/correction/RNAseqRepaired.tsv"
DEGsDF <- read.delim(DEGsFile) 
DEGsDF <- reshape(DEGsDF, direction = "wide", idvar = "target", timevar = "contrast")

TPMsFile <- "data/RNA-Seq data/result/counts_trans.csv"
TPMsDF <- read.delim(TPMsFile, sep = ",") 
colnames(TPMsDF)[2:ncol(TPMsDF)] <- paste0(colnames(TPMsDF)[2:ncol(TPMsDF)],".TPMcounts")
TPMsDF$X <- gsub("_.*","",TPMsDF$X)

DEGsDFfull <- merge(TPMsDF,DEGsDF,by.x="X",by.y="target")
DEGsDFfull$target <- DEGsDFfull$X

colnames(DEGsDFfull) <- paste0("DEGs.",colnames(DEGsDFfull))

whatwehaveMapping <- merge(DEGsDFfull,yeastGenesProtMap,by.x = "DEGs.target",by.y = "ORF")
whatwehaveMapping$ORF <- whatwehaveMapping$DEGs.target
whatwehaveNOTMapping <- DEGsDFfull[!DEGsDFfull$DEGs.target %in% yeastGenesProtMap$ORF,]

whatwehaveMapping2 <- merge(whatwehaveNOTMapping,yeastGenesProtMap,by.x = "DEGs.target",by.y = "Alias")
whatwehaveMapping2$Alias <- whatwehaveMapping2$DEGs.target

whatwehaveNOTMapping2 <- whatwehaveNOTMapping[!whatwehaveNOTMapping$DEGs.target %in% yeastGenesProtMap$Alias,]

whatwehaveMapping3 <- merge(whatwehaveNOTMapping2,yeastGenesProtMap,by.x = "DEGs.target",by.y = "GeneName")
whatwehaveMapping3$GeneName <- whatwehaveMapping3$DEGs.target
whatwehaveNOTMapping3 <- whatwehaveNOTMapping2[!whatwehaveNOTMapping2$DEGs.target %in% yeastGenesProtMap$GeneName,]

theAliases <- sapply(whatwehaveNOTMapping3$DEGs.target, function(x) {grep(paste0("\\b",x,"\\b"),yeastGenesProtMap$Alias,ignore.case = F)},simplify = "array")
mappedByAlias <- as.data.frame(cbind(DEGs.target=names(unlist(theAliases)),Alias=yeastGenesProtMap$Alias[unlist(theAliases)]))
whatwehaveMapping4 <- merge(mappedByAlias,yeastGenesProtMap,by = "Alias")
whatwehaveNOTMapping4 <- names(theAliases)[!names(theAliases) %in% names(unlist(theAliases))]

theAliases2 <- sapply(whatwehaveNOTMapping4, function(x) {grep(x,yeastGenesProtMap$Alias,ignore.case = F,fixed = T)},simplify = "array")
mappedByAlias2 <- as.data.frame(cbind(DEGs.target=names(unlist(theAliases2)),Alias=yeastGenesProtMap$Alias[unlist(theAliases2)]))
whatwehaveMapping5 <- merge(mappedByAlias2,yeastGenesProtMap,by = "Alias")
whatwehaveNOTMapping5 <- names(theAliases2)[!names(theAliases2) %in% names(unlist(theAliases2))]
length(whatwehaveNOTMapping5)

whatwehaveMapping4n5 <- merge(whatwehaveNOTMapping3,rbind(whatwehaveMapping4,whatwehaveMapping5),by = "DEGs.target") 

DEGsMapped <- rbind(whatwehaveMapping,whatwehaveMapping2,whatwehaveMapping3,whatwehaveMapping4n5)
table(DEGsMapped$type)
DEGsMapped$DEGs.X <- NULL
DEGsMappedFile <- 'FinalData/DEGsMapped.tsv'
write.table(DEGsMapped,DEGsMappedFile,quote = F,sep = '\t',col.names = T,row.names = F)

################################################################################
########################## Merging Proteomes and RBPomes #######################
################################################################################

ProteomeFAXfile <- "FinalData/Raw/correctedRAW/ProteomeFAX_np2_norm_gmin_no15-25-26-36-39/ProteomeFAX_np2_norm_gmin_no15-25-26-36-39_PROTEIN.tsv"
ProteomeUVXfile <- "FinalData/Raw/correctedRAW/ProteomeUVX_np2_norm_gmin_no4-6-23-31-32/ProteomeUVX_np2_norm_gmin_no4-6-23-31-32_PROTEIN.tsv"
RBPomeFAXfile <-  "FinalData/Raw/correctedRAW/RBPomeFAX_np2_norm_gmin_no12/RBPomeFAX_np2_norm_gmin_no12_PROTEIN.tsv"
RBPomeUVXfile <- "FinalData/Raw/correctedRAW/RBPomeUVX_np2_norm_gmin_no10-7-16-3/RBPomeUVX_np2_norm_gmin_no10-7-16-3_PROTEIN.tsv"

ProteomeFAX <- read.delim(ProteomeFAXfile); nrow(ProteomeFAX)
ProteomeUVX <- read.delim(ProteomeUVXfile); nrow(ProteomeUVX)
RBPomeFAX <- read.delim(RBPomeFAXfile); nrow(RBPomeFAX)
RBPomeUVX <- read.delim(RBPomeUVXfile); nrow(RBPomeUVX)

# protinfo <- RBPomeUVX$proteinName
# totalProts <- length(protinfo); totalProts
# nProtsOrth <- sum(grepl(';',unique(protinfo))); nProtsOrth
# nProts <- totalProts - nProtsOrth; nProts

commonCols <- colnames(ProteomeUVX)[1:8] #intersect(colnames(ProteomeFAX),colnames(ProteomeUVX))

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

outdir <- "FinalData"
saveTablesTsvExc(FAXdf,outdir,completeNdedup=F,excel=T,bycompleteFC=F,rownames=F)
saveTablesTsvExc(UVXdf,outdir,completeNdedup=F,excel=T,bycompleteFC=F,rownames=F)

################################################################################
########################## MERGING FAX and UVX #################################
################################################################################

FAXnUVXdf <- merge(FAXdf,UVXdf,all = T, by = "proteinName"); nrow(FAXdf); nrow(FAXnUVXdf)

toVennList <- list(FAX = FAXdf$proteinName, UVX = UVXdf$proteinName)
ggVennDiagram(toVennList, color = 2, lwd = 0.7) + scale_fill_gradient(low = "#F4FAFE", high = "#4981BF") + theme(legend.position = "none")

outdir <- "FinalData"
saveTablesTsvExc(FAXnUVXdf,outdir,completeNdedup=F,excel=T,bycompleteFC=F,rownames=F)

################################################################################
########################## MAPPING FAX and UVX #################################
################################################################################

FAXnUVXFile <- "FinalData/FAXnUVXdf.tsv"
FAXnUVXdf <- read.delim(FAXnUVXFile)

FAXnUVXSingProts <- FAXnUVXdf[!grepl(";",FAXnUVXdf$proteinName),]; nrow(FAXnUVXSingProts)

toVennList <- list(mapping = yeastGenesProtMap$UniprotACC, FAXnUVX = FAXnUVXSingProts$proteinName)
ggVennDiagram(toVennList, color = 2, lwd = 0.7) + scale_fill_gradient(low = "#F4FAFE", high = "#4981BF") + theme(legend.position = "none")

FAXnUVXmapped <- merge(yeastGenesProtMap,FAXnUVXSingProts,by.x = "UniprotACC", by.y = "proteinName")
FAXnUVXmapped$proteinName <- FAXnUVXmapped$UniprotACC
nrow(FAXnUVXSingProts); nrow(FAXnUVXmapped)

# FAXnUVXmapped$UniprotACC[duplicated(FAXnUVXmapped$UniprotACC)] 
# View(FAXnUVXmapped[FAXnUVXmapped$UniprotACC %in% FAXnUVXmapped$UniprotACC[duplicated(FAXnUVXmapped$UniprotACC)],])
# These proteins are located in two different chromosomes of Yeast and are behind a duplicity

FAXnUVXunmapped <- FAXnUVXSingProts[!(FAXnUVXSingProts$proteinName %in% yeastGenesProtMap$UniprotACC),]

toVennList <- list(mapping = yeastGenesProtMap$GeneName, FAXnUVX = FAXnUVXunmapped$geneName.ProteomeFAX)
ggVennDiagram(toVennList, color = 2, lwd = 0.7) + scale_fill_gradient(low = "#F4FAFE", high = "#4981BF") + theme(legend.position = "none")

FAXnUVXmapped2 <- merge(yeastGenesProtMap,FAXnUVXunmapped,by.x="GeneName",by.y="geneName.ProteomeFAX")
FAXnUVXmapped2$geneName.ProteomeFAX <- FAXnUVXmapped2$GeneName

FAXnUVXunmapped <- FAXnUVXunmapped[!(FAXnUVXunmapped$geneName.ProteomeFAX %in% yeastGenesProtMap$GeneName),]
nrow(FAXnUVXunmapped)

# FAXnUVXmapped$UniprotACC[duplicated(FAXnUVXmapped$UniprotACC)] 
# View(FAXnUVXmapped[FAXnUVXmapped$UniprotACC %in% FAXnUVXmapped$UniprotACC[duplicated(FAXnUVXmapped$UniprotACC)],])
# These proteins are located in two different chromosomes of Yeast and are behind a duplicity

library(tidyr)

FAXnUVXOrthologs <- FAXnUVXdf[grepl(";",FAXnUVXdf$proteinName),]; nrow(FAXnUVXOrthologs)
FAXnUVXOrthologs$orthologsIDs <- 1:nrow(FAXnUVXOrthologs)

sepFAXnUVXOrthologs <- separate_rows(FAXnUVXOrthologs, proteinName, sep = ";", convert = FALSE)
nrow(FAXnUVXOrthologs); nrow(sepFAXnUVXOrthologs)

toVennList <- list(mapping = yeastGenesProtMap$UniprotACC, FAXnUVX = sepFAXnUVXOrthologs$proteinName)
ggVennDiagram(toVennList, color = 2, lwd = 0.7) + scale_fill_gradient(low = "#F4FAFE", high = "#4981BF") + theme(legend.position = "none")

sepFAXnUVXOrthologsMapped <- merge(yeastGenesProtMap,sepFAXnUVXOrthologs,by.x = "UniprotACC", by.y = "proteinName")
sepFAXnUVXOrthologsMapped$proteinName <- sepFAXnUVXOrthologsMapped$UniprotACC
nrow(sepFAXnUVXOrthologs); nrow(sepFAXnUVXOrthologsMapped)
sepFAXnUVXOrthologsUnmapped <- sepFAXnUVXOrthologs[!(sepFAXnUVXOrthologs$proteinName %in% yeastGenesProtMap$UniprotACC),]

sepFAXnUVXOrthologsUnmapped$proteinName <- gsub("-2","",sepFAXnUVXOrthologsUnmapped$proteinName)
sepFAXnUVXOrthologsMapped2 <- merge(yeastGenesProtMap,sepFAXnUVXOrthologsUnmapped,by.x = "UniprotACC", by.y = "proteinName")

sepFAXnUVXOrthologsUnmapped <- sepFAXnUVXOrthologsMapped2[!(sepFAXnUVXOrthologsMapped2$proteinName %in% yeastGenesProtMap$UniprotACC),]
nrow(sepFAXnUVXOrthologsUnmapped)

sepFAXnUVXOrthologsMapped2$UniprotACC <- paste0(sepFAXnUVXOrthologsMapped2$UniprotACC,"-2")
sepFAXnUVXOrthologsMapped2$proteinName <- sepFAXnUVXOrthologsMapped2$UniprotACC

orthologsMapped <- rbind(sepFAXnUVXOrthologsMapped,sepFAXnUVXOrthologsMapped2)
nrow(orthologsMapped)
orthologsMapped <- orthologsMapped[!duplicated(orthologsMapped),]
nrow(orthologsMapped)

orthologsMappedMinInfo <- orthologsMapped[c("orthologsIDs", "proteinName", colnames(yeastGenesProtMap))]
nrow(orthologsMappedMinInfo); #View(orthologsMappedMinInfo)

require(dplyr)


orthologsMappedMinInfo[is.na(orthologsMappedMinInfo)] <- ""
orthologsCollapsed <- summarise_each(group_by(orthologsMappedMinInfo,orthologsIDs),funs(paste(., collapse = ";")))
orthologsCollapsed <- orthologsCollapsed[order(orthologsCollapsed$orthologsIDs,decreasing = F),]
nrow(orthologsCollapsed); #View(FAXnUVXOrthologs)

orthologsFullMap <- merge(orthologsCollapsed,FAXnUVXOrthologs,by="orthologsIDs")
orthologsFullMap$proteinName.x <- NULL 
orthologsFullMap$orthologsIDs <- NULL
colnames(orthologsFullMap)[17] <- gsub("\\.y","",colnames(orthologsFullMap)[17])

FAXnUVXfullMapped <- rbind(FAXnUVXmapped,FAXnUVXmapped2,orthologsFullMap)
FAXnUVXfullMapped <- FAXnUVXfullMapped[,c(colnames(yeastGenesProtMap), colnames(FAXnUVXfullMapped)[!(colnames(FAXnUVXfullMapped) %in% colnames(yeastGenesProtMap))])]

outdir <- "FinalData"
saveTablesTsvExc(FAXnUVXfullMapped,outdir,completeNdedup=F,excel=T,bycompleteFC=F,rownames=F)

################################################################################
###################### MERGIN Transcriptome and FAXnUVX ########################
################################################################################

FAXnUVXfullMapped <- read.delim("FinalData/FAXnUVXfullMapped.tsv",quote = "")
intersect(colnames(DEGsMapped),colnames(FAXnUVXfullMapped))

FAXnUVXfullMap <- FAXnUVXfullMapped[!grepl(";",FAXnUVXfullMapped$proteinName),]; nrow(FAXnUVXfullMap)
FAXnUVXfullMapOrtho <- FAXnUVXfullMapped[grepl(";",FAXnUVXfullMapped$proteinName),]; nrow(FAXnUVXfullMapOrtho)

toVennList <- list(DEGs = DEGsMapped$UniprotACC, FAXnUVX = FAXnUVXfullMap$UniprotACC)
ggVennDiagram(toVennList, color = 2, lwd = 0.7) + scale_fill_gradient(low = "#F4FAFE", high = "#4981BF") + theme(legend.position = "none")

toVennList <- list(DEGs = DEGsMapped$ORF, FAXnUVX = FAXnUVXfullMap$ORF)
ggVennDiagram(toVennList, color = 2, lwd = 0.7) + scale_fill_gradient(low = "#F4FAFE", high = "#4981BF") + theme(legend.position = "none")

DEGsnFAXnUVX <- merge(DEGsMapped,FAXnUVXfullMap,by="ORF",all = T)
checkMerge <- as.data.frame(cbind(DEGsnFAXnUVX$Gene.Names.x,DEGsnFAXnUVX$Gene.Names.y))
checkMerge <- checkMerge[complete.cases(checkMerge),]
if (!identical(checkMerge$V1,checkMerge$V2)){
  View(as.data.frame(cbind(checkMerge$V1,checkMerge$V2))[checkMerge$V1!=checkMerge$V2,])
}

for (xCol in colnames(DEGsnFAXnUVX)[grep("\\.x",colnames(DEGsnFAXnUVX))]){
  yCol <- gsub("\\.x",".y",xCol)
  DEGsnFAXnUVX[,xCol] <- coalesce(DEGsnFAXnUVX[,xCol],DEGsnFAXnUVX[,yCol])
  DEGsnFAXnUVX[,yCol] <- NULL
}
colnames(DEGsnFAXnUVX) <- gsub("\\.x","",colnames(DEGsnFAXnUVX))

library(plyr)
allDataDF <- rbind.fill(DEGsnFAXnUVX,FAXnUVXfullMapOrtho)

outdir <- "FinalData"
saveTablesTsvExc(allDataDF,outdir,completeNdedup=F,excel=T,bycompleteFC=F,rownames=F)

### FINAL TABLES FOR IBTISSAM

FAXacceptedProts <- read.delim("FinalData/BackgroundRemoval/FAXAccBackGrAcceptedFC3.tsv",header = F)[,1]
UVXacceptedProts <- read.delim("FinalData/BackgroundRemoval/UVXAccBackGrAcceptedFC3.tsv",header = F)[,1]
allDataDF <- read.delim("FinalData/allDataDF.tsv")

FAXmRBP <- allDataDF[which((!is.na(allDataDF$RBPomeFAX.log2ratio_PolyARNAFAXwithSA))),]
length(FAXmRBP$proteinName); length(unique(FAXmRBP$proteinName))

FAXwoBKGR <- FAXmRBP[which(FAXmRBP$proteinName %in% FAXacceptedProts),]
length(FAXwoBKGR$proteinName); length(unique(FAXwoBKGR$proteinName))

# The proteins P10081 P02309 P32324 P02994 are located in two different chromosomal position this means an extra 4 rows due to mapping
FAXwoBKGR$proteinName[duplicated(FAXwoBKGR$proteinName)]
View(FAXwoBKGR[FAXwoBKGR$proteinName %in% FAXwoBKGR$proteinName[duplicated(FAXwoBKGR$proteinName)],])

UVXwoBKGR <- allDataDF[which((!is.na(allDataDF$RBPomeUVX.log2ratio_PolyARNAUVwithSA)) & allDataDF$proteinName %in% UVXacceptedProts),]
length(UVXwoBKGR$proteinName); length(unique(UVXwoBKGR$proteinName))

UVXmRBP <- allDataDF[which((!is.na(allDataDF$RBPomeUVX.log2ratio_PolyARNAUVwithSA))),]
length(UVXmRBP$proteinName); length(unique(UVXmRBP$proteinName))

UVXwoBKGR <- UVXmRBP[which(UVXmRBP$proteinName %in% UVXacceptedProts),]
length(UVXwoBKGR$proteinName); length(unique(UVXwoBKGR$proteinName))

View(UVXwoBKGR[UVXwoBKGR$proteinName %in% UVXwoBKGR$proteinName[duplicated(UVXwoBKGR$proteinName)],])
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
                      'Gene.Names..ORF.','Gene.Names..primary.','Gene.Names..synonym.','Alias',
                      'DEGs.H202.brep1.rawCounts','DEGs.H202.brep2.rawCounts','DEGs.H202.brep3.rawCounts',
                      'DEGs.DTT.brep1.rawCounts','DEGs.DTT.brep2.rawCounts','DEGs.DTT.brep3.rawCounts',
                      'DEGs.SA.brep1.rawCounts','DEGs.SA.brep2.rawCounts','DEGs.SA.brep3.rawCounts',
                      'DEGs.YPD.brep1.rawCounts','DEGs.YPD.brep2.rawCounts','DEGs.YPD.brep3.rawCounts',
                      'DEGs.adj.pval.H202-YPD','DEGs.log2FC.H202-YPD','DEGs.adj.pval.DTT-YPD',
                      'DEGs.log2FC.DTT-YPD','DEGs.adj.pval.SA-YPD','DEGs.log2FC.SA-YPD','ProteomeFAX.X014_Pr40.DIA_DIA_D22',
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

columnsSelected[!columnsSelected %in% colnames(FAXwoBKGR)]
columnsSelected[!columnsSelected %in% colnames(UVXwoBKGR)]

FAXwoBKGR <- FAXwoBKGR[,columnsSelected]
UVXwoBKGR <- UVXwoBKGR[,columnsSelected]
outdir <- "FinalData/TablesForPublication"
saveTablesTsvExc(FAXwoBKGR,outdir,completeNdedup=F,excel=T,bycompleteFC=F,rownames=F)
saveTablesTsvExc(UVXwoBKGR,outdir,completeNdedup=F,excel=T,bycompleteFC=F,rownames=F)


