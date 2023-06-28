require(dplyr)
library(tidyr)
# % of RBPome detected as significant in SA UVX and FAX

noxRBPome <- read.delim("/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/data/raw/MQtoSQ.csv", sep = ",")
faxRBPome <- read.delim("/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/data/raw/MQtoSQ_FAX.csv", sep = ",")
uvxRBPome <- read.delim("/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/data/raw/MQtoSQ_UV.csv", sep = ",")
proteome <- read.delim("/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/data/raw/SNtoSQ-Prot40_SafeQuant_Input.csv", sep = ",")

length(proteome$X.10) - 2
length(unique(proteome$X.10)) - 2

length(noxRBPome$X.10) - 2
length(unique(noxRBPome$X.10)) - 2 

length(faxRBPome$X.10) - 2
length(unique(faxRBPome$X.10)) - 2

length(uvxRBPome$X.10) - 2
length(unique(uvxRBPome$X.10)) - 2




wholeDF <- read.delim("/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/FinalData/allDataDF.tsv",quote = "")
SARBPFAXsubdf <- wholeDF[,c("proteinName","RBPomeFAX.qValue_PolyARNAFAXwithSA")]

length(unique(wholeDF$proteinName[!is.na(wholeDF$RBPomeFAX.qValue_PolyARNAFAXwithSA)]))
length(unique(wholeDF$proteinName[!is.na(wholeDF$RBPomeUVX.qValue_PolyARNAUVwithSA)]))

SARBPUVXsubdf <- wholeDF[,c("proteinName","RBPomeUVX.qValue_PolyARNAUVwithSA")]


NOXdfFile <- "/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/FinalData/BackgroundRemoval/sqaProteinMean.tsv"
NOXdf <- read.delim(NOXdfFile)
nrow(NOXdf)
sum(abs(NOXdf$Condition2) > 1)
sum(abs(NOXdf$Condition3) > 1)
  

wholeDF$RBPomeUVX.qValue_PolyARNAUVwithSA


SARBPFAXsubdf <- SARBPFAXsubdf[complete.cases(SARBPFAXsubdf),]
SARBPFAXsubdf <- SARBPFAXsubdf[!duplicated(SARBPFAXsubdf),]

SARBPFAXsigProtsNparalgs <- SARBPFAXsubdf$RBPomeFAX.qValue_PolyARNAFAXwithSA < 0.05

paste("From",nrow(SARBPFAXsubdf),"prots and paralogs",
      round(sum(SARBPFAXsigProtsNparalgs) / nrow(SARBPFAXsubdf) * 100,2),"is significant:", sum(SARBPFAXsigProtsNparalgs))

SARBPUVXsubdf <- wholeDF[,c("proteinName","RBPomeUVX.qValue_PolyARNAUVwithSA")]
SARBPUVXsubdf <- SARBPUVXsubdf[complete.cases(SARBPUVXsubdf),]
SARBPUVXsubdf <- SARBPUVXsubdf[!duplicated(SARBPUVXsubdf),]
SARBPUVXsigProtsNparalgs <- SARBPUVXsubdf$RBPomeUVX.qValue_PolyARNAUVwithSA < 0.05

paste("From",nrow(SARBPUVXsubdf),"prots and paralogs",
      round(sum(SARBPUVXsigProtsNparalgs) / nrow(SARBPUVXsubdf) * 100,2),"is significant:", sum(SARBPUVXsigProtsNparalgs))


FAXdffile <- "/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/FinalData/FAXdf.tsv"
FAXdfMapedFullfile <- "/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/FinalData/FAXdfMapedFull.tsv"
FAXdfMapedFull <- read.delim(FAXdfMapedFullfile, quote = "")
FAXdf <- read.delim(FAXdffile,quote = "")

length(unique(wholeDF$proteinName[!is.na(wholeDF)]))

wholeDF$proteinName

wholeDFsigSA <- wholeDF[which(wholeDF$RBPomeFAX.qValue_PolyARNAFAXwithSA < 0.05),c("proteinName","RBPomeFAX.qValue_PolyARNAFAXwithSA")]
length(which(wholeDFsigSA$RBPomeFAX.qValue_PolyARNAFAXwithSA < 0.05))
wholeDFsigSA <- wholeDFsigSA[!duplicated(wholeDFsigSA),]
length(which(wholeDFsigSA$RBPomeFAX.qValue_PolyARNAFAXwithSA < 0.05))

wholeDF$UniprotACC

wholeDFsigSA <- wholeDF[which(wholeDF$RBPomeFAX.qValue_PolyARNAFAXwithSA < 0.05),]
wholeDFsigSA <- wholeDF[duplicated(wholeDFsigSA),]


wholeDFsigSA <- wholeDF[which(wholeDF$RBPomeFAX.qValue_PolyARNAFAXwithSA < 0.05),c("UniprotACC","RBPomeFAX.qValue_PolyARNAFAXwithSA")]
length(which(wholeDFsigSA$RBPomeFAX.qValue_PolyARNAFAXwithSA < 0.05))
wholeDFsigSA <- wholeDFsigSA[!duplicated(wholeDFsigSA),]
length(which(wholeDFsigSA$RBPomeFAX.qValue_PolyARNAFAXwithSA < 0.05))



wholeDFsigSA <- FAXdfMapedFull[which(FAXdfMapedFull$RBPomeFAX.qValue_PolyARNAFAXwithSA < 0.05),c("proteinName","RBPomeFAX.qValue_PolyARNAFAXwithSA")]
length(which(wholeDFsigSA$RBPomeFAX.qValue_PolyARNAFAXwithSA < 0.05))
wholeDFsigSA <- wholeDFsigSA[!duplicated(wholeDFsigSA),]
length(which(wholeDFsigSA$RBPomeFAX.qValue_PolyARNAFAXwithSA < 0.05))

wholeDFsigSA <- FAXdf[which(FAXdf$RBPomeFAX.qValue_PolyARNAFAXwithSA < 0.05),c("proteinName","RBPomeFAX.qValue_PolyARNAFAXwithSA")]
length(which(wholeDFsigSA$RBPomeFAX.qValue_PolyARNAFAXwithSA < 0.05))
wholeDFsigSA <- wholeDFsigSA[!duplicated(wholeDFsigSA),]
length(which(wholeDFsigSA$RBPomeFAX.qValue_PolyARNAFAXwithSA < 0.05))


length(which(wholeDF$RBPomeFAX.qValue_PolyARNAFAXwithSA < 0.05))
length(which(FAXdf$RBPomeFAX.qValue_PolyARNAFAXwithSA < 0.05))



FAXdfMapedFull$RBPomeFAX.qValue_PolyARNAFAXwithSA < 0.05


### Merged file of UVX and FAX mRBPome where 1 column allows to filter by condition dectetion combinations.

FAXrbpomeFile <- '/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/FinalData/Raw/correctedRAW/FAXmRBPomeRAW.csv'
UVXrbpomeFile <- '/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/FinalData/Raw/correctedRAW/UVXmRBPomeRAW.csv'
FAXrbpomeDF <- getProteinsMatrixFromRAW(FAXrbpomeFile)
UVXrbpomeDF <- getProteinsMatrixFromRAW(UVXrbpomeFile)

#FAXrbpomeDF[FAXrbpomeDF$proteins == 'O13535;P0C2I9;P0C2J1;P47100;Q03612;Q04214;Q04670;Q04711;Q12112;Q12273;Q12490;Q99231',]
#UVXrbpomeDF[UVXrbpomeDF$proteins == 'O13535;P0C2I9;P0C2J1;P47100;Q03612;Q04214;Q04670;Q04711;Q12112;Q12273;Q12490;Q99231',]

FAXnUVXdf <- merge(FAXrbpomeDF,UVXrbpomeDF,by = "proteins",all = T)

myclasses <- gsub('.*_','',colnames(FAXnUVXdf))[2:ncol(FAXnUVXdf)]
myclass <- myclasses[1]

filteringColumn <- rep('',nrow(FAXnUVXdf))
for (myclass in unique(myclasses)){
  mySubExprMatrix <- FAXnUVXdf[grep(myclass,colnames(FAXnUVXdf))]
  proteinsDetected <- which(rowSums(is.na(mySubExprMatrix)) != ncol(mySubExprMatrix))
  filteringColumn[proteinsDetected] <- paste(filteringColumn,myclass,sep = ',')[proteinsDetected]
}
filteringColumn <- trimws(filteringColumn,whitespace = ',')
filteringColumn[filteringColumn == paste(unique(myclasses), collapse = ',')] <- 'All'
rawFAXnUVXdf_withTreatmentDetectionFilter <- cbind(FAXnUVXdf,classPresentFilter = filteringColumn)
outdir <- 'FinalData/Raw/correctedRAW/'
saveTablesTsvExc(rawFAXnUVXdf_withTreatmentDetectionFilter,outdir,completeNdedup=F,excel=T,bycompleteFC=F,rownames=F)

mappingFinalFile <- "yeastReference/mappingFile.xlsx"
yeastGenesProtMap <- read.xlsx(mappingFinalFile)
yeastGenesProtMap[which(yeastGenesProtMap$GeneName == 'REP2'),]$UniprotACC <- 'P03872'

rawFAXnUVXdf_withTreatmentDetectionFilter <- read.xlsx('FinalData/Raw/correctedRAW/rawFAXnUVXdf_withTreatmentDetectionFilter.xlsx')

FAXnUVXdf <- rawFAXnUVXdf_withTreatmentDetectionFilter 
colnames(FAXnUVXdf)[1] <- 'proteinName'

FAXnUVXSingProts <- FAXnUVXdf[!grepl(";",FAXnUVXdf$proteinName),]; nrow(FAXnUVXSingProts)

FAXnUVXmapped <- merge(yeastGenesProtMap,FAXnUVXSingProts,by.x = "UniprotACC", by.y = "proteinName")
FAXnUVXmapped$proteinName <- FAXnUVXmapped$UniprotACC
nrow(FAXnUVXSingProts); nrow(FAXnUVXmapped)

FAXnUVXunmapped <- FAXnUVXSingProts[!(FAXnUVXSingProts$proteinName %in% yeastGenesProtMap$UniprotACC),]

FAXnUVXmapped$GeneName[which(is.na(FAXnUVXmapped$GeneName))] <- FAXnUVXmapped$ORF[which(is.na(FAXnUVXmapped$GeneName))]

FAXnUVXOrthologs <- FAXnUVXdf[grepl(";",FAXnUVXdf$proteinName),]; nrow(FAXnUVXOrthologs)
FAXnUVXOrthologs$orthologsIDs <- 1:nrow(FAXnUVXOrthologs)

sepFAXnUVXOrthologs <- separate_rows(FAXnUVXOrthologs, proteinName, sep = ";", convert = FALSE)
nrow(FAXnUVXOrthologs); nrow(sepFAXnUVXOrthologs)
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
orthologsMapped <- orthologsMapped[!duplicated(orthologsMapped),]

orthologsMappedMinInfo <- orthologsMapped[c("orthologsIDs", "proteinName", colnames(yeastGenesProtMap))]
sum(is.na(orthologsMappedMinInfo$GeneName))

orthologsMappedMinInfo$GeneName[which(is.na(orthologsMappedMinInfo$GeneName))] <- orthologsMappedMinInfo$ORF[which(is.na(orthologsMappedMinInfo$GeneName))]

nrow(orthologsMappedMinInfo); #View(orthologsMappedMinInfo)

orthologsMappedMinInfo[is.na(orthologsMappedMinInfo)] <- ""
orthologsCollapsed <- summarise_each(group_by(orthologsMappedMinInfo,orthologsIDs),funs(paste(., collapse = ";")))
orthologsCollapsed <- orthologsCollapsed[order(orthologsCollapsed$orthologsIDs,decreasing = F),]
nrow(orthologsCollapsed); #View(FAXnUVXOrthologs)

orthologsFullMap <- merge(orthologsCollapsed,FAXnUVXOrthologs,by="orthologsIDs")
orthologsFullMap$proteinName.x <- NULL 
orthologsFullMap$orthologsIDs <- NULL
colnames(orthologsFullMap)[17] <- gsub("\\.y","",colnames(orthologsFullMap)[17])

FAXnUVXfullMapped <- rbind(FAXnUVXmapped,orthologsFullMap)
FAXnUVXfullMapped <- FAXnUVXfullMapped[,c(colnames(yeastGenesProtMap), colnames(FAXnUVXfullMapped)[!(colnames(FAXnUVXfullMapped) %in% colnames(yeastGenesProtMap))])]

rawFAXnUVXdf_withTreatmentDetectionFilter <- FAXnUVXfullMapped
outdir <- 'FinalData/Raw/correctedRAW/'
rawFAXnUVXdf_withTreatmentDetectionFilter_Mapped <- rawFAXnUVXdf_withTreatmentDetectionFilter[!duplicated(rawFAXnUVXdf_withTreatmentDetectionFilter),]
any(is.na(rawFAXnUVXdf_withTreatmentDetectionFilter_Mapped$proteinName))

saveTablesTsvExc(rawFAXnUVXdf_withTreatmentDetectionFilter_Mapped,outdir,completeNdedup=F,excel=T,bycompleteFC=F,rownames=F)


