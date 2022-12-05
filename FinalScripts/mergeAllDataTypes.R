
DEGfile <- "/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/FinalData/DEGsFullMapping.tsv"
UVXfile <- "/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/FinalData/UVXdfMapedFull.tsv"
FAXfile <- "/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/FinalData/FAXdfMapedFull.tsv"

DEGdf <- read.delim(DEGfile,quote = "")
UVXdf <- read.delim(UVXfile,quote = "")
FAXdf <- read.delim(FAXfile,quote = "")

colnames(DEGdf)[grep("target",colnames(DEGdf)):ncol(DEGdf)] <- paste0("DEGs.",colnames(DEGdf)[grep("target",colnames(DEGdf)):ncol(DEGdf)])

allDataDF <- merge(merge(DEGdf,FAXdf,all = T),UVXdf,all = T)


badFormatedProtNames <- allDataDF$proteinName[which(grepl("\\|",allDataDF$proteinName))] 
protNameCorrected <- c()
for (badFormatedProtName in badFormatedProtNames){
  badFormatedProtNameSep <- unlist(strsplit(badFormatedProtName,split = "\\|"))
  badFormatedProtNameStr <- paste(badFormatedProtNameSep[seq(2,length(badFormatedProtNameSep),2)],collapse = ";")
  protNameCorrected <- c(protNameCorrected,badFormatedProtNameStr)
}
allDataDF$proteinName[which(grepl("\\|",allDataDF$proteinName))] <- protNameCorrected

checkingMerging <- cbind(allDataDF$proteinName[!is.na(allDataDF$proteinName)],allDataDF$UniprotACC[!is.na(allDataDF$proteinName)])
checkingMerging <- checkingMerging[complete.cases(checkingMerging),]
if (identical(checkingMerging[,1],checkingMerging[,2])){
  print("Merge correct")
}

outdir <- "/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/FinalData"
saveTablesTsvExc(allDataDF,outdir,completeNdedup = F,excel = T,bycompleteFC = F,rownames = F)



## Short script cause all worked smoothly













