library(ggVennDiagram)
library(ggplot2)
RBPomeUVXAdriFile <- "/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/data/raw/RBPomeUVX_np2_norm_gmin_no10-7-16-3/RBPomeUVX_np2_norm_gmin_no10-7-16-3_PROTEIN.tsv"
RBPomeUVXAdriDF <- read.delim(RBPomeUVXAdriFile,quote = '"',check.names = T)
colnames(RBPomeUVXAdriDF) <- gsub("\\.$","",gsub("X.","",colnames(RBPomeUVXAdriDF)))

RBPomeUVXOldFile <- "/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/data/old/P445_AG_PolyA_40_all_3_cond_20220611_Final/P445_AG_PolyA_40_20220617_gmin_UV_p2_no511_Final/P445_AG_PolyA_40_20220617_gmin_UV_p2_no511_PROTEIN.tsv"
RBPomeUVXOldDF <- read.delim(RBPomeUVXOldFile,quote = '"')

RBPomeUVXOldProts <- RBPomeUVXOldDF$proteinName[RBPomeUVXOldDF$qValue_PolyARNAUVwithSA < 0.05]
RBPomeUVXAdriProts <- RBPomeUVXAdriDF$proteinName[RBPomeUVXAdriDF$qValue_PolyARNAUVwithSA < 0.05]
toVennList <- list(RBPomeUVX_no10_7 = RBPomeUVXOldProts, RBPomeUVXAdriProts_no10_7_16_3 = RBPomeUVXAdriProts)

ggvennSA <- ggVennDiagram(toVennList, color = 2, lwd = 0.7) + 
  scale_fill_gradient(low = "#F4FAFE", high = "#4981BF") +
  theme(legend.position = "none")

nameOut <- "/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/vennUVXRBPomeOldvsNew.png"
ggsave(filename = nameOut, plot = ggvennSA, width = 30, height = 15, units = 'cm', dpi = 'print')

##########

backgroundRBPomeFile <- "/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/FinalData/BackgroundRemoval/sqaProteinMean.tsv"
backgroundRBPomeDF <- read.delim(backgroundRBPomeFile)

# UniprotsCoded <- c()
# for (badFormatedProtName in row.names(backgroundRBPomeDF)){
#   badFormatedProtNameSep <- unlist(strsplit(badFormatedProtName,split = "\\|"))
#   badFormatedProtNameStr <- paste(badFormatedProtNameSep[seq(2,length(badFormatedProtNameSep),2)],collapse = ";")
#   UniprotsCoded <- c(UniprotsCoded,badFormatedProtNameStr)
# }

FCcutOff <- 1
abovebackgrProtsFAX <- row.names(backgroundRBPomeDF)[which(abs(backgroundRBPomeDF$Condition2) > FCcutOff)]
abovebackgrProtsUVX <- row.names(backgroundRBPomeDF)[which(abs(backgroundRBPomeDF$Condition3) > FCcutOff)]
abovebackgrProtsFAXCorrected <- UniprotsCoded[which(abs(backgroundRBPomeDF$Condition2) > FCcutOff)]
abovebackgrProtsUVXCorrected <- UniprotsCoded[which(abs(backgroundRBPomeDF$Condition3) > FCcutOff)]

toVennList <- list(allProteins = UniprotsCoded, FAXacceptedBG = abovebackgrProtsFAXCorrected, UVCacceptedBG=abovebackgrProtsUVXCorrected)
ggvennSA <- ggVennDiagram(toVennList, color = 2, lwd = 0.7) + 
  scale_fill_gradient(low = "#F4FAFE", high = "#4981BF") +
  theme(legend.position = "none")
nameOut <- "/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/FinalData/BackgroundRemoval/backgroundStatus.tiff"
ggsave(filename = nameOut, plot = ggvennSA, width = 30, height = 15, units = 'cm', dpi = 'print')


overlapORFs <- read.xlsx("/home/eidriangm/Downloads/overlap analysis.xlsx")

backgroundORFs <- read.delim("/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/FinalData/BackgroundRemoval/backgroundHMAP.tsv")

backgroundORFsMapped <-  merge(yeastGenesProtMap,backgroundORFs,by.x = "UniprotACC", by.y = "row.names",all.y = T)
backgroundORFsMapped <- backgroundORFsMapped[,c("UniprotACC","ORF")]
backgroundORFsMapped$ORF[is.na(backgroundORFsMapped$ORF)] <- backgroundORFsMapped$UniprotACC[is.na(backgroundORFsMapped$ORF)]
backgroundORFsMapped <- backgroundORFsMapped[!duplicated(backgroundORFsMapped),]

sum(grepl(";", backgroundORFsMapped$UniprotACC))

toVennList <- list(Dataset.1 = overlapORFs$Dataset.1, Dataset.2 = overlapORFs$Dataset.2, Background=backgroundORFsMapped$ORF)
ggvennSA <- ggVennDiagram(toVennList, color = 2, lwd = 0.7) + 
  scale_fill_gradient(low = "#F4FAFE", high = "#4981BF") +
  theme(legend.position = "none")
nameOut <- "/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/FinalData/BackgroundRemoval/ORFsBackgroundStatus.tiff"
ggsave(filename = nameOut, plot = ggvennSA, width = 30, height = 15, units = 'cm', dpi = 'print')




