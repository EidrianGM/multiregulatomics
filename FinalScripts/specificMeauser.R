# % of RBPome detected as significant in SA UVX and FAX
wholeDF <- read.delim("/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/FinalData/wholeDF.tsv",quote = "")

SARBPFAXsubdf <- wholeDF[,c("proteinName","RBPomeFAX.qValue_PolyARNAFAXwithSA")]
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
