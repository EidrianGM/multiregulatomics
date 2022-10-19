library(openxlsx)
library(rtracklayer)
library(plyr)
library(reshape2)

##################################
##### DATA and MAPPING FILES #####
##################################

wholeDFfile <- "/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/data/mapped/allDataMapped_together.tsv"
wholeDF <- read.delim(wholeDFfile,quote = "")
mappingFinalFile <- "/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/yeastReference/mappingFile.xlsx"
yeastGenesProtMap <- read.xlsx(mappingFinalFile)

#######################################
##### Protein Net Activities Data #####
#######################################

wholeDF$netChangesFAXDTT <- wholeDF$mRNABPomeFAX.log2ratio_PolyARNAFAXwithDTT - wholeDF$proteomeFAX.log2ratio_Condition2
wholeDF$netChangesFAXH2O2 <- wholeDF$mRNABPomeFAX.log2ratio_PolyARNAFAXwithH2O2 - wholeDF$proteomeFAX.log2ratio_Condition3
wholeDF$netChangesFAXSA <- wholeDF$mRNABPomeFAX.log2ratio_PolyARNAFAXwithSA - wholeDF$proteomeFAX.log2ratio_Condition4

wholeDF$netChangesUVDTT <- wholeDF$mRNABPomeUV.log2ratio_PolyARNAUVwithDTT - wholeDF$proteomeUV.log2ratio_Condition2
wholeDF$netChangesUVH2O <- wholeDF$mRNABPomeUV.log2ratio_PolyARNAUVwithH202 - wholeDF$proteomeUV.log2ratio_Condition3
wholeDF$netChangesUVSA <- wholeDF$mRNABPomeUV.log2ratio_PolyARNAUVwithSA - wholeDF$proteomeUV.log2ratio_Condition4

saveTablesTsvExc(wholeDF,"mapped",F)

# 1- RNABProteins FAX
ibtissamDir <- "/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/data/mapped/forIbtissam"
degscols <- colnames(wholeDF)[grep("DEGs",colnames(wholeDF))]

FAXcolsNetChanges <- c("mRNABPomeFAX.pValue_PolyARNAFAXwithDTT","proteomeFAX.pValue_Condition2",
                       "mRNABPomeFAX.qValue_PolyARNAFAXwithDTT","proteomeFAX.qValue_Condition2",
                       "mRNABPomeFAX.log2ratio_PolyARNAFAXwithDTT","proteomeFAX.log2ratio_Condition2", "netChangesFAXDTT",
                       
                       "mRNABPomeFAX.pValue_PolyARNAFAXwithH2O2","proteomeFAX.pValue_Condition3",
                       "mRNABPomeFAX.qValue_PolyARNAFAXwithH2O2","proteomeFAX.qValue_Condition3",
                       "mRNABPomeFAX.log2ratio_PolyARNAFAXwithH2O2","proteomeFAX.log2ratio_Condition3", "netChangesFAXH2O2",
                       
                       "mRNABPomeFAX.pValue_PolyARNAFAXwithSA","proteomeFAX.pValue_Condition4",
                       "mRNABPomeFAX.qValue_PolyARNAFAXwithSA","proteomeFAX.qValue_Condition4",
                       "mRNABPomeFAX.log2ratio_PolyARNAFAXwithSA","proteomeFAX.log2ratio_Condition4", "netChangesFAXSA")

FAXprotNrbpome <- wholeDF[, c(colnames(yeastGenesProtMap),degscols,FAXcolsNetChanges)]
saveTablesTsvExc(FAXprotNrbpome,ibtissamDir)

# 2- RNABProteins UV



UVprotNrbpomeCols <- c("mRNABPomeUV.pValue_PolyARNAUVwithDTT","proteomeUV.pValue_Condition2",
                       "mRNABPomeUV.qValue_PolyARNAUVwithDTT","proteomeUV.qValue_Condition2",
                       "mRNABPomeUV.log2ratio_PolyARNAUVwithDTT","proteomeUV.log2ratio_Condition2", "netChangesUVDTT",
                       
                       "mRNABPomeUV.pValue_PolyARNAUVwithH202","proteomeUV.pValue_Condition3",
                       "mRNABPomeUV.qValue_PolyARNAUVwithH202","proteomeUV.qValue_Condition3",
                       "mRNABPomeUV.log2ratio_PolyARNAUVwithH202","proteomeUV.log2ratio_Condition3", "netChangesUVH2O2",
                       
                       "mRNABPomeUV.pValue_PolyARNAUVwithSA","proteomeUV.pValue_Condition4",
                       "mRNABPomeUV.qValue_PolyARNAUVwithSA","proteomeUV.qValue_Condition4",
                       "mRNABPomeUV.log2ratio_PolyARNAUVwithSA","proteomeUV.log2ratio_Condition4", "netChangesUVSA")

UVprotNrbpome <- wholeDF[, c(colnames(yeastGenesProtMap),degscols,UVprotNrbpomeCols)]
saveTablesTsvExc(UVprotNrbpome,ibtissamDir)

#UVprotNrbpomeOnlyNetChanges <- wholeDF[(!is.na(netChangesUVDTT) | !is.na(netChangesUVH2O2) | !is.na(netChangesUVSA)), c(colnames(yeastGenesProtMap),UVprotNrbpomeCols)]
#saveTablesTsvExc(UVprotNrbpomeOnlyNetChanges,ibtissamDir)

saveTablesTsvExc(wholeDF,"/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/data/mapped")

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

