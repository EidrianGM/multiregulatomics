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



