library(ggVennDiagram)
library(ggplot2)

## Number of proteins detected in NOX, UVX and FAX in Proteome, for each stress 
## condition at least one reading per sample of the same class
Proteomefile <- "FinalData/Raw/correctedRAW/newProteome.csv"
ProteomeSamplesInfo <- read.delim(Proteomefile,nrows = 3, sep = ",")
myDF <- ProteomeSamplesInfo
normDCol <- grep("Norm",colnames(myDF)) # Raw Normalized Values are not the ones!
rawDCol <- grep("Raw",colnames(myDF))
expresionMatrix <- myDF[normDCol:(rawDCol-1)] 
expresionMatrix[2,] <- gsub("_.*","",expresionMatrix[2,])
phenData <- as.data.frame(t(expresionMatrix[1:2,])); colnames(phenData) <- c("class","id"); 
phenData$class <- gsub(".*\\s","",phenData$class); conditionCol <- "class"
rownames(phenData) <- 1:nrow(phenData)

ProteomeDF <- read.delim(Proteomefile,skip = 3, sep = ",",header = F)
expresionMatrix <- ProteomeDF[c(11,normDCol:(rawDCol-1))] 
colnames(expresionMatrix) <- c('protein',phenData$id)
expresionMatrix <- expresionMatrix[3:nrow(expresionMatrix),]

expresionMatrix$protein

phenData$class <- gsub("no.*",'',phenData$class)

totalProteins <- unique(expresionMatrix$protein); length(totalProteins)
proteinsDetectionPerClass <- list()
for (myclass in unique(phenData$class)){
  mySubExprMatrix <- expresionMatrix[phenData$id[phenData$class == myclass]]
  proteinsDetected <- unique(expresionMatrix$protein[which(rowSums(is.na(mySubExprMatrix)) != ncol(mySubExprMatrix))])
  proteinsDetectionPerClass[[myclass]] <- proteinsDetected
}
View(proteinsDetectionPerClass)
unique(phenData$class)

toVennListControls <- list(NOX = proteinsDetectionPerClass[['NOX']], UVX = proteinsDetectionPerClass[['UV']], FAX=proteinsDetectionPerClass[['FAX']])
toVennListDTT <- list(NOX = proteinsDetectionPerClass[['NOXwithDTT']], UVX = proteinsDetectionPerClass[['UVwithDTT']], FAX=proteinsDetectionPerClass[['FAXwithDTT']])
toVennListH202 <- list(NOX = proteinsDetectionPerClass[['NOXwithH2O2']], UVX = proteinsDetectionPerClass[['UVwithH202']], FAX=proteinsDetectionPerClass[['FAXwithH2O2']])
toVennListSA <- list(NOX = proteinsDetectionPerClass[['NOXwithSA']], UVX = proteinsDetectionPerClass[['UVwithSA']], FAX=proteinsDetectionPerClass[['FAXwithSA']])

ggvennControls <- ggVennDiagram(toVennListControls, color = 2, lwd = 0.7) + 
            scale_fill_gradient(low = "#F4FAFE", high = "#4981BF") +
            theme(legend.position = "none")
nameOut <- "FinalData/QualityControlPlots/ProteomeProteinsDetectedVenns/Controls.png"
ggsave(filename = nameOut, plot = ggvennControls, width = 30, height = 15, units = 'cm', dpi = 'print')

ggvennDTT <- ggVennDiagram(toVennListDTT, color = 2, lwd = 0.7) + 
  scale_fill_gradient(low = "#F4FAFE", high = "#4981BF") +
  theme(legend.position = "none")
nameOut <- "FinalData/QualityControlPlots/ProteomeProteinsDetectedVenns/DTT.png"
ggsave(filename = nameOut, plot = ggvennDTT, width = 30, height = 15, units = 'cm', dpi = 'print')

ggvennH202 <- ggVennDiagram(toVennListH202, color = 2, lwd = 0.7) + 
  scale_fill_gradient(low = "#F4FAFE", high = "#4981BF") +
  theme(legend.position = "none")
nameOut <- "FinalData/QualityControlPlots/ProteomeProteinsDetectedVenns/H2O2.png"
ggsave(filename = nameOut, plot = ggvennH202, width = 30, height = 15, units = 'cm', dpi = 'print')

ggvennSA <- ggVennDiagram(toVennListSA, color = 2, lwd = 0.7) + 
  scale_fill_gradient(low = "#F4FAFE", high = "#4981BF") +
  theme(legend.position = "none")
nameOut <- "FinalData/QualityControlPlots/ProteomeProteinsDetectedVenns/SA.png"
ggsave(filename = nameOut, plot = ggvennSA, width = 30, height = 15, units = 'cm', dpi = 'print')


################# 


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

FAXfinalDFfile <- "FinalData/TablesForPublication/FAXwoBKGR.tsv"
UVXfinalDFfile <- "FinalData/TablesForPublication/UVXwoBKGR.tsv"
FAXfinalDF <- read.delim(FAXfinalDFfile)
UVXfinalDF <- read.delim(UVXfinalDFfile)

toVennList <- list(FAX = FAXfinalDF$proteinName, UVX = UVXfinalDF$proteinName)
ggvennSA <- ggVennDiagram(toVennList, color = 2, lwd = 0.7) + 
  scale_fill_gradient(low = "#F4FAFE", high = "#4981BF") +
  theme(legend.position = "none")
nameOut <- "FinalData/TablesForPublication/FAXnUVXwoBKGR.tiff"
ggsave(filename = nameOut, plot = ggvennSA, width = 10, height = 10, units = 'cm', dpi = 'print')

common <- intersect(FAXfinalDF$proteinName,UVXfinalDF$proteinName)
FAXesp <- setdiff(FAXfinalDF$proteinName,UVXfinalDF$proteinName)
UVXesp <- setdiff(UVXfinalDF$proteinName,FAXfinalDF$proteinName)

FAXnUVXdf <- merge(FAXfinalDF,UVXfinalDF,all = T)
vennStatus <- 1:nrow(FAXnUVXdf)
vennStatus[FAXnUVXdf$proteinName %in% common] <- "common"
vennStatus[FAXnUVXdf$proteinName %in% FAXesp] <- "FAX"
vennStatus[FAXnUVXdf$proteinName %in% UVXesp] <- "UVX"
FAXnUVXdf <- FAXnUVXdf[,colnames(FAXnUVXdf)[1:17]]
FAXnUVXdf$vennStatus <- vennStatus
outdir <- "FinalData/TablesForPublication"
FAXnUVXvennDF <- FAXnUVXdf
saveTablesTsvExc(FAXnUVXvennDF,outdir,completeNdedup=F,excel=T,bycompleteFC=F,rownames=F)

FAXnUVXvennDF$proteinName[FAXnUVXvennDF$proteinName %in% FAXnUVXvennDF$proteinName[duplicated(FAXnUVXvennDF$proteinName)]]








