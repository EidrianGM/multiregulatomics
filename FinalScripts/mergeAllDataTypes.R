
DEGfile <- "/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/FinalData/DEGsFullMapping.tsv"
UVXfile <- "/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/FinalData/UVXdfMapedFull.tsv"
FAXfile <- "/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/FinalData/FAXdfMapedFull.tsv"

DEGdf <- read.delim(DEGfile,quote = "")
UVXdf <- read.delim(UVXfile,quote = "")
FAXdf <- read.delim(FAXfile,quote = "")

colnames(DEGdf)[grep("target",colnames(DEGdf)):ncol(DEGdf)] <- paste0("DEGs.",colnames(DEGdf)[grep("target",colnames(DEGdf)):ncol(DEGdf)])

allDataDF <- merge(merge(DEGdf,FAXdf,all = T),UVXdf,all = T)

# Keep Only UniProts in Orthologues
orthologuesmanyset <- allDataDF$proteinName[grep(";sp",wholeDF$proteinName)]
orthologuesmanyCorrected <- c()
for (orthologuesmany in orthologuesmanyset){
  orthologuesmanySep <- unlist(strsplit(orthologuesmany,split = "\\|"))
  orthologuesmanyTxt <- paste(orthologuesmanySep[seq(2,length(orthologuesmanySep),2)],collapse = ";")
  orthologuesmanyCorrected <- c(orthologuesmanyCorrected,orthologuesmanyTxt)
}
allDataDF$proteinName[grep(";sp",allDataDF$proteinName)] <- orthologuesmanyCorrected


outdir <- "/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/FinalData"
saveTablesTsvExc(allDataDF,outdir,completeNdedup = F,excel = T,bycompleteFC = F,rownames = F)

## Short script cause all worked smoothly













