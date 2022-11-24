library(openxlsx)
mappingFinalFile <- "/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/yeastReference/mappingFile.xlsx"
yeastGenesProtMap <- read.xlsx(mappingFinalFile)

########################
##### MAPPING DEGs #####
########################

DEGsFile <- "/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/data/RNA-Seq data/result/DE gene testing statistics.csv"
DEGsDF <- read.delim(DEGsFile, sep = ",") 

TPMsFile <- "/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/data/RNA-Seq data/result/TPM_genes.csv"
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

DEGsDFMapedFull <- DEGsDFMapedFull[c("target","adj.pval.H202-YPD","log2FC.H202-YPD","adj.pval.DTT-YPD","log2FC.DTT-YPD","adj.pval.SA-YPD","log2FC.SA-YPD",'seqnames','start','end','type','ORFgene','ORF','SGDID','GeneName','UniprotACC','UniprotName','Gene.Names','Protein.names','Gene.Names..ordered.locus.','Gene.Names..ORF.','Gene.Names..primary.','Gene.Names..synonym.','Alias')]

DEGsDFMapedFullFile <- '/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/FinalData/DEGsFullMapping.tsv'
write.table(DEGsDFMapedFull,DEGsDFMapedFullFile,quote = F,sep = '\t',col.names = T,row.names = F)


################################################################################
############################## MAPPING Proteome and RBPome #####################
################################################################################

# 1 UVX
UVXdfFile <- "/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/FinalData/UVXdf.tsv"
UVXdf <- read.delim(UVXdfFile); nrow(UVXdf)
UVXdf$ac.ProteomeUVX[is.na(UVXdf$ac.ProteomeUVX)] <- UVXdf$ac.RBPomeUVX[is.na(UVXdf$ac.ProteomeUVX)]
UVXdf$geneName.ProteomeUVX[is.na(UVXdf$geneName.ProteomeUVX)] <- UVXdf$geneName.RBPomeUVX[is.na(UVXdf$geneName.ProteomeUVX)]

length(UVXdf$proteinName); length(unique(UVXdf$proteinName))

UVXdfOrthologs <- UVXdf[grepl(";",UVXdf$proteinName),]; nrow(UVXdfOrthologs)
UVXdfSingProts <- UVXdf[!grepl(";",UVXdf$proteinName),]; nrow(UVXdfSingProts)

UVXdfMaped <- merge(UVXdfSingProts,yeastGenesProtMap,by.x="proteinName",by.y="UniprotACC")
UVXdfMaped$UniprotACC <- UVXdfMaped$proteinName
unmapped <- !(UVXdfSingProts$proteinName %in% UVXdfMaped$proteinName); sum(unmapped)
UVXdfUnMap <- UVXdfSingProts[unmapped,]

UVXdfMaped2 <- merge(UVXdfUnMap,yeastGenesProtMap,by.x="ac.ProteomeUVX",by.y="UniprotACC")
UVXdfMaped2$UniprotACC <- UVXdfMaped2$ac.ProteomeUVX
unmapped <- !(UVXdfUnMap$ac.ProteomeUVX %in% UVXdfMaped2$ac.ProteomeUVX); sum(unmapped)
UVXdfUnMap <- UVXdfUnMap[unmapped,] # The Weird Case of REP2 ProteomeUVX (Bad Roll Up?)

UVXdfMaped3 <- merge(UVXdfUnMap,yeastGenesProtMap,by.x="geneName.ProteomeUVX",by.y="GeneName")
UVXdfMaped3$GeneName <- UVXdfMaped3$geneName.ProteomeUVX
unmapped <- !(UVXdfUnMap$geneName.ProteomeUVX %in% UVXdfMaped3$geneName.ProteomeUVX); sum(unmapped)

UVXdfMapedFull <- rbind.fill(UVXdfMaped,UVXdfMaped2,UVXdfMaped3,UVXdfOrthologs)
UVXdfMapedFull <- UVXdfMapedFull[,c(colnames(yeastGenesProtMap),colnames(UVXdf))]

all(UVXdf$proteinName %in% UVXdfMapedFull$proteinName) # Must be True

outdir <- "/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/FinalData"
saveTablesTsvExc(UVXdfMapedFull,outdir,completeNdedup = F,excel = T,bycompleteFC = F,rownames = F)


# 2 FAX
FAXdfFile <- "/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/FinalData/FAXdf.tsv"
FAXdf <- read.delim(FAXdfFile); nrow(FAXdf)

FAXdf$ac.ProteomeFAX[is.na(FAXdf$ac.ProteomeFAX)] <- FAXdf$ac.RBPomeFAX[is.na(FAXdf$ac.ProteomeFAX)]
FAXdf$geneName.ProteomeFAX[is.na(FAXdf$geneName.ProteomeFAX)] <- FAXdf$geneName.RBPomeFAX[is.na(FAXdf$geneName.ProteomeFAX)]

length(FAXdf$proteinName); length(unique(FAXdf$proteinName))

FAXdfOrthologs <- FAXdf[grepl(";",FAXdf$proteinName),]; nrow(FAXdfOrthologs)
FAXdfSingProts <- FAXdf[!grepl(";",FAXdf$proteinName),]; nrow(FAXdfSingProts)

FAXdfMaped <- merge(FAXdfSingProts,yeastGenesProtMap,by.x="proteinName",by.y="UniprotACC")
FAXdfMaped$UniprotACC <- FAXdfMaped$proteinName
unmapped <- !(FAXdfSingProts$proteinName %in% FAXdfMaped$proteinName); sum(unmapped)
FAXdfUnMap <- FAXdfSingProts[unmapped,]

FAXdfMaped2 <- merge(FAXdfUnMap,yeastGenesProtMap,by.x="ac.ProteomeFAX",by.y="UniprotACC")
FAXdfMaped2$UniprotACC <- FAXdfMaped2$ac.ProteomeFAX
unmapped <- !(FAXdfUnMap$ac.ProteomeFAX %in% FAXdfMaped2$ac.ProteomeFAX); sum(unmapped)
FAXdfUnMap <- FAXdfUnMap[unmapped,] # The Weird Case of REP2 ProteomeFAX (Bad Roll Up?)

FAXdfMaped3 <- merge(FAXdfUnMap,yeastGenesProtMap,by.x="geneName.ProteomeFAX",by.y="GeneName")
FAXdfMaped3$GeneName <- FAXdfMaped3$geneName.ProteomeFAX
unmapped <- !(FAXdfUnMap$geneName.ProteomeFAX %in% FAXdfMaped3$geneName.ProteomeFAX); sum(unmapped)

FAXdfMapedFull <- rbind.fill(FAXdfMaped,FAXdfMaped2,FAXdfMaped3,FAXdfOrthologs)
FAXdfMapedFull <- FAXdfMapedFull[,c(colnames(yeastGenesProtMap),colnames(FAXdf))]

all(FAXdf$proteinName %in% FAXdfMapedFull$proteinName) # Must be True

outdir <- "/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/FinalData"
saveTablesTsvExc(FAXdfMapedFull,outdir,completeNdedup = F,excel = T,bycompleteFC = F,rownames = F)

