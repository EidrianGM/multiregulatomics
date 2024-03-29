library(ggplot2)

#### Mean/Median SELECTION
medianNOXfaxFile <- "FinalData/Raw/correctedRAW/medianRBPomeNOX_FAX_np2_NOnorm/medianRBPomeNOX_FAX_np2_NOnorm_PROTEIN.tsv"
medianNOXuvxFile <- "FinalData/Raw/correctedRAW/medianRBPomeNOX_UVX_np2_NOnorm/medianRBPomeNOX_UVX_np2_NOnorm_PROTEIN.tsv"
meanNOXfaxFile <- "FinalData/Raw/correctedRAW/meanRBPomeNOX_FAX_np2_NOnorm/meanRBPomeNOX_FAX_np2_NOnorm_PROTEIN.tsv"
meanNOXuvxFile <- "FinalData/Raw/correctedRAW/meanRBPomeNOX_UVX_np2_NOnorm/meanRBPomeNOX_UVX_np2_NOnorm_PROTEIN.tsv"
outFile <- "FinalData/BackgroundRemoval/ByProtstatsBGremoval.tsv"
medianNOXfaxDF <- read.delim(medianNOXfaxFile);
medianNOXuvxDF <- read.delim(medianNOXuvxFile);
meanNOXfaxDF <- read.delim(meanNOXfaxFile);
meanNOXuvxDF <- read.delim(meanNOXuvxFile);

# medianNOXfaxFile <- "FinalData/Raw/correctedRAW/byPEPmedianRBPomeNOX_FAX_np2_NOnorm/byPEPmedianRBPomeNOX_FAX_np2_NOnorm_PROTEIN.tsv"
# medianNOXuvxFile <- "FinalData/Raw/correctedRAW/byPEPmedianRBPomeNOX_UVX_np2_NOnorm/byPEPmedianRBPomeNOX_UVX_np2_NOnorm_PROTEIN.tsv"
# meanNOXfaxFile <- "FinalData/Raw/correctedRAW/byPEPmeanRBPomeNOX_FAX_np2_NOnorm/byPEPmeanRBPomeNOX_FAX_np2_NOnorm_PROTEIN.tsv"
# meanNOXuvxFile <- "FinalData/Raw/correctedRAW/byPEPmeanRBPomeNOX_UVX_np2_NOnorm/byPEPmeanRBPomeNOX_UVX_np2_NOnorm_PROTEIN.tsv"
# outFile <- "FinalData/BackgroundRemoval/ByPEPstatsBGremoval.tsv"
# medianNOXfaxDF <- read.delim(medianNOXfaxFile);
# medianNOXuvxDF <- read.delim(medianNOXuvxFile);
# meanNOXfaxDF <- read.delim(meanNOXfaxFile);
# meanNOXuvxDF <- read.delim(meanNOXuvxFile);


###### mRBPomes UVX and FAX filtering

FAXfile <- "FinalData/Raw/correctedRAW/RBPomeFAX_np2_norm_gmin_no12/RBPomeFAX_np2_norm_gmin_no12_PROTEIN.tsv"
UVXfile <- "FinalData/Raw/correctedRAW/RBPomeUVX_np2_norm_gmin_no10-7-16-3/RBPomeUVX_np2_norm_gmin_no10-7-16-3_PROTEIN.tsv"

FAXdf <- read.delim(FAXfile) 
UVXdf <- read.delim(UVXfile)

if(all(FAXdf$proteinName %in% medianNOXfaxDF$proteinName)){
  print("Perfect all proteins in FAX data are in NOX")
}

if(all(UVXdf$proteinName %in% meanNOXuvxDF$proteinName)){
  print("Perfect all proteins in UVX data are in NOX")
} 

#### VennDiagrams to quantify background removal
library(ggVennDiagram)
library(ggplot2)

statsBGremoval <- c("CrossLink","FCmethod","FCcutOff","protsInNOX2crossLink","protsInCrosslink","acceptedProts","removedProts","notInResults")
protsInNOXfax <- length(medianNOXfaxDF$proteinName) 
protsInNOXuvx <- length(medianNOXuvxDF$proteinName)
protsInFAX <- length(FAXdf$proteinName)
protsInUVX <- length(UVXdf$proteinName)

FCcutOff <- 3
for (FCcutOff in 1:5){
  meanabovebackgrProtsFAX <- meanNOXfaxDF$proteinName[which(abs(meanNOXfaxDF$log2ratio_Condition2) > FCcutOff)]
  medianabovebackgrProtsFAX <- medianNOXfaxDF$proteinName[which(abs(medianNOXfaxDF$log2ratio_Condition2) > FCcutOff)]
  
  toVennList <- list(FAXproteins = FAXdf$proteinName, meanFAXacceptedBGprots = meanabovebackgrProtsFAX, medianFAXacceptedBGprots=medianabovebackgrProtsFAX)
  ggvennSA <- ggVennDiagram(toVennList, color = 2, lwd = 0.7) + scale_fill_gradient(low = "#F4FAFE", high = "#4981BF") + theme(legend.position = "none")
  nameOut <- paste0("FinalData/BackgroundRemoval/vennFAX_FC",FCcutOff,"bg.png")
  ggsave(filename = nameOut, plot = ggvennSA, width = 30, height = 15, units = 'cm', dpi = 'print')

  meanabovebackgrProtsUVX <- meanNOXuvxDF$proteinName[which(abs(meanNOXuvxDF$log2ratio_Condition2) > FCcutOff)]
  medianabovebackgrProtsUVX <- medianNOXuvxDF$proteinName[which(abs(medianNOXuvxDF$log2ratio_Condition2) > FCcutOff)]
  
  toVennList <- list(UVXproteins = UVXdf$proteinName, meanUVXacceptedBGprots = meanabovebackgrProtsUVX,medianFAXacceptedBGprots=medianabovebackgrProtsUVX)
  ggvennSA <- ggVennDiagram(toVennList, color = 2, lwd = 0.7) + scale_fill_gradient(low = "#F4FAFE", high = "#4981BF") + theme(legend.position = "none")
  nameOut <- paste0("FinalData/BackgroundRemoval/vennUVX_FC",FCcutOff,"bg.png")
  ggsave(filename = nameOut, plot = ggvennSA, width = 30, height = 15, units = 'cm', dpi = 'print')
  
  meanFAXaccepedProts <- length(intersect(FAXdf$proteinName,meanabovebackgrProtsFAX))
  meanFAXremovedProts <- length(setdiff(FAXdf$proteinName,meanabovebackgrProtsFAX))
  notinFAXProts <- length(setdiff(meanabovebackgrProtsFAX,FAXdf$proteinName))
  medianFAXaccepedProts <- length(intersect(FAXdf$proteinName,medianabovebackgrProtsFAX))
  medianFAXremovedProts <- length(setdiff(FAXdf$proteinName,medianabovebackgrProtsFAX))
  notinFAXProts <- length(setdiff(medianabovebackgrProtsFAX,FAXdf$proteinName))
  
  meanUVXaccepedProts <- length(intersect(UVXdf$proteinName,meanabovebackgrProtsUVX))
  meanUVXremovedProts <- length(setdiff(UVXdf$proteinName,meanabovebackgrProtsUVX))
  notinUVXProts <- length(setdiff(meanabovebackgrProtsUVX,UVXdf$proteinName))
  medianUVXaccepedProts <- length(intersect(UVXdf$proteinName,medianabovebackgrProtsUVX))
  medianUVXremovedProts <- length(setdiff(UVXdf$proteinName,medianabovebackgrProtsUVX))
  notinUVXProts <- length(setdiff(medianabovebackgrProtsUVX,UVXdf$proteinName))
  
  statsBGremoval <- rbind(statsBGremoval,
    cbind("FAX","mean",FCcutOff,protsInNOXfax,protsInFAX,meanFAXaccepedProts,meanFAXremovedProts,notinFAXProts),
    cbind("FAX","median",FCcutOff,protsInNOXfax,protsInFAX,medianFAXaccepedProts,medianFAXremovedProts,notinFAXProts),
    cbind("UVX","mean",FCcutOff,protsInNOXuvx,protsInUVX,meanUVXaccepedProts,meanUVXremovedProts,notinUVXProts),
    cbind("UVX","median",FCcutOff,protsInNOXuvx,protsInUVX,medianUVXaccepedProts,medianUVXremovedProts,notinUVXProts))
}

write.table(statsBGremoval,outFile,row.names = F,col.names = F,quote = F, sep = "\t")


## Final Selection # Stablish a FCcutoff and FC calc method to select the background
### MEAN as FC method and NOX raw filtered by proteins
FCcutOff <- 3

meanNOXfaxFile <- "FinalData/Raw/correctedRAW/meanRBPomeNOX_FAX_np2_NOnorm/meanRBPomeNOX_FAX_np2_NOnorm_PROTEIN.tsv"
meanNOXuvxFile <- "FinalData/Raw/correctedRAW/meanRBPomeNOX_UVX_np2_NOnorm/meanRBPomeNOX_UVX_np2_NOnorm_PROTEIN.tsv"
meanNOXfaxDF <- read.delim(meanNOXfaxFile);
meanNOXuvxDF <- read.delim(meanNOXuvxFile);

abovebackgrProtsFAX <- meanNOXfaxDF$proteinName[which(abs(meanNOXfaxDF$log2ratio_Condition2) > FCcutOff)]
abovebackgrProtsUVX <- meanNOXuvxDF$proteinName[which(abs(meanNOXuvxDF$log2ratio_Condition2) > FCcutOff)]

write(abovebackgrProtsFAX,"FinalData/BackgroundRemoval/FAXAccBackGrAcceptedFC3.tsv", sep = "\n")
write(abovebackgrProtsUVX,"FinalData/BackgroundRemoval/UVXAccBackGrAcceptedFC3.tsv", sep = "\n")

belowbackgrProtsFAX <- meanNOXfaxDF$proteinName[which(abs(meanNOXfaxDF$log2ratio_Condition2) < FCcutOff)]
belowbackgrProtsUVX <- meanNOXuvxDF$proteinName[which(abs(meanNOXuvxDF$log2ratio_Condition2) < FCcutOff)]

write(belowbackgrProtsFAX,"FinalData/BackgroundRemoval/FAXAccBackGrRemovedFC3.tsv", sep = "\n")
write(belowbackgrProtsUVX,"FinalData/BackgroundRemoval/UVXAccBackGrRemovedFC3.tsv", sep = "\n")

nrow(meanNOXfaxDF)
length(abovebackgrProtsFAX)
length(belowbackgrProtsFAX)

nrow(meanNOXuvxDF)
length(abovebackgrProtsUVX)
length(belowbackgrProtsUVX)

##### GENERATE HEATMAPS
protsInFAX <- length(FAXdf$proteinName)
protsInUVX <- length(UVXdf$proteinName)

selectedProtsFAX <- intersect(FAXdf$proteinName,abovebackgrProtsFAX)
selectedProtsUVX <- intersect(UVXdf$proteinName,abovebackgrProtsUVX)

length(selectedProtsFAX)
sum(grepl(';',selectedProtsFAX))
sum(!grepl(';',selectedProtsFAX))

length(selectedProtsUVX)
sum(grepl(';',selectedProtsUVX))
sum(!grepl(';',selectedProtsUVX))

toVennList <- list(FAXproteins = FAXdf$proteinName, meanFAXacceptedBGprots = abovebackgrProtsFAX)
ggvennSA <- ggVennDiagram(toVennList, color = 2, lwd = 0.7) + scale_fill_gradient(low = "#F4FAFE", high = "#4981BF") + theme(legend.position = "none")
nameOut <- paste0("FinalData/BackgroundRemoval/FAX_ACCEPTEDbackground_FC",FCcutOff,"_venn.png")
ggsave(filename = nameOut, plot = ggvennSA, width = 30, height = 15, units = 'cm', dpi = 'print')

toVennList <- list(UVXproteins = UVXdf$proteinName, meanUVXacceptedBGprots = abovebackgrProtsUVX)
ggvennSA <- ggVennDiagram(toVennList, color = 2, lwd = 0.7) + scale_fill_gradient(low = "#F4FAFE", high = "#4981BF") + theme(legend.position = "none")
nameOut <- paste0("FinalData/BackgroundRemoval/UVX_ACCEPTEDbackground_FC",FCcutOff,"_venn.png")
ggsave(filename = nameOut, plot = ggvennSA, width = 30, height = 15, units = 'cm', dpi = 'print')

nPeps <- length(peptidesInfo$Accession) - nPepsOrth
nProts <- length(unique(peptidesInfo$Accession)) - nProtsOrth


##### GENERATE HEATMAPS
meanNOXfaxFile <- "FinalData/Raw/correctedRAW/meanRBPomeNOX_FAX_np2_NOnorm/meanRBPomeNOX_FAX_np2_NOnorm_PROTEIN.tsv"
meanNOXuvxFile <- "FinalData/Raw/correctedRAW/meanRBPomeNOX_UVX_np2_NOnorm/meanRBPomeNOX_UVX_np2_NOnorm_PROTEIN.tsv"
meanNOXfaxDF <- read.delim(meanNOXfaxFile)
meanNOXuvxDF <- read.delim(meanNOXuvxFile)

samples2include <- c(colnames(meanNOXfaxDF),colnames(meanNOXuvxDF))[grep("X0",c(colnames(meanNOXfaxDF),colnames(meanNOXuvxDF)))]
samples2include <- unique(as.numeric(gsub("X0|_.*","",samples2include)))

NOXrawFile <- "FinalData/Raw/correctedRAW/NOXmRBPomeRAW.csv"
myDF <- read.delim(NOXrawFile,nrows = 3,sep = ",")
normDCol <- grep("Norm",colnames(myDF)) # Raw Normalized Values are not the ones!
rawDCol <- grep("Raw",colnames(myDF))
expresionMatrix <- myDF[normDCol:(rawDCol-1)] 
expresionMatrix[2,] <- gsub("_.*","",expresionMatrix[2,])
subPhenData <- as.data.frame(t(expresionMatrix[1:2,])); colnames(subPhenData) <- c("class","id"); 
subPhenData$class <- gsub(".*\\s","",subPhenData$class); conditionCol <- "class"
rownames(subPhenData) <- 1:nrow(subPhenData)
subPhenData$class <- gsub("no.*|with.*","",subPhenData$class)
subPhenData$id <- as.numeric(subPhenData$id)
subPhenData <- subPhenData[subPhenData$id %in% samples2include,]
subPhenData$class <- gsub("UV","UVX",subPhenData$class)
subPhenData$class <- factor(subPhenData$class,levels = c("NOX","FAX", "UVX"))

saveTablesTsvExc(subPhenData,"FinalData/BackgroundRemoval",completeNdedup=F,excel=T,bycompleteFC=F,rownames=T)

#### Filter and get expression for HeatMap
meanNOXfaxDF <- meanNOXfaxDF[order(meanNOXfaxDF$log2ratio_Condition2,decreasing = T),]
meanNOXuvxDF <- meanNOXuvxDF[order(meanNOXuvxDF$log2ratio_Condition2,decreasing = T),]

FCcutOff <- 3
abovebackgrProtsFAX <- meanNOXfaxDF$proteinName[which(abs(meanNOXfaxDF$log2ratio_Condition2) > FCcutOff)]
abovebackgrProtsUVX <- meanNOXuvxDF$proteinName[which(abs(meanNOXuvxDF$log2ratio_Condition2) > FCcutOff)]
abovebackgrProts <- unique(c(abovebackgrProtsFAX,abovebackgrProtsUVX))

rownames(meanNOXfaxDF) <- meanNOXfaxDF$proteinName
meanNOXfaxDF <- meanNOXfaxDF[,grep("^X",colnames(meanNOXfaxDF))]
colnames(meanNOXfaxDF) <- as.numeric(gsub("X|_.*","",colnames(meanNOXfaxDF)))

rownames(meanNOXuvxDF) <- meanNOXuvxDF$proteinName
meanNOXuvxDF <- meanNOXuvxDF[,grep("^X",colnames(meanNOXuvxDF))]
colnames(meanNOXuvxDF) <- as.numeric(gsub("X|_.*","",colnames(meanNOXuvxDF)))

FAXnUVXbackgroundDF <- merge(meanNOXfaxDF, meanNOXuvxDF,by="row.names",all=T,suffixes = c("","y"))
row.names(FAXnUVXbackgroundDF) <- FAXnUVXbackgroundDF$Row.names
FAXnUVXbackgroundDF$Row.names <- NULL

FAXnUVXbackgroundDF <- FAXnUVXbackgroundDF[abovebackgrProts,]

library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)

outdir <- "FinalData/BackgroundRemoval"
tagname <- "Background"

FAXnUVXbackgroundDF <- FAXnUVXbackgroundDF[colnames(FAXnUVXbackgroundDF)[grep("y",colnames(FAXnUVXbackgroundDF),invert = T)]]
FAXnUVXbackgroundDF[is.na(FAXnUVXbackgroundDF)] <- 0

toHMP <- log10(FAXnUVXbackgroundDF+1)
toHMP <- toHMP[,as.character(subPhenData$id)]
toHMP <- toHMP[rev(row.names(toHMP)),]
colours <- brewer.pal(brewer.pal.info["PuRd",]$maxcolors,"PuRd")
colourStep <- (max(toHMP,na.rm = T) / brewer.pal.info["PuRd",]$maxcolors)
breaksColours <- seq(0,max(toHMP,na.rm = T)-colourStep,colourStep)
col_fun <- colorRamp2(breaksColours, colours) # #574a85
myHMPplot <- Heatmap(as.matrix(toHMP),cluster_columns = F, cluster_rows = T, col=col_fun,
                     show_row_names = F,show_row_dend = F, row_km=6,
                     name = "Log10(LFQ)", #top_annotation = colannot,
                     column_names_rot = 0, column_split = subPhenData$class,
                     heatmap_legend_param = list(
                       legend_direction = "horizontal", 
                       legend_width = unit(5, "cm")
                     ), column_names_side = "top", row_gap = unit(10, "mm")
)
tiff(file.path(outdir,paste0(gsub("\ ", "",tagname),'HeatMap.tiff')), 1000, nrow(toHMP)*5, pointsize=12,res = "print")
draw(myHMPplot, heatmap_legend_side="bottom", annotation_legend_side = "left")
dev.off()

ht <- draw(myHMPplot)
htrowOrder <- row_order(ht)

protsCluster <- c()
for (cluster in names(htrowOrder)){
  protsCluster <- rbind(protsCluster,cbind(htrowOrder[[cluster]],cluster))
}
protsCluster <- as.data.frame(protsCluster)
colnames(protsCluster) <- c("idx","cluster")
protsCluster$proteins <- row.names(toHMP)[as.numeric(protsCluster$idx)]
protsCluster$nProtsInClust <- table(protsCluster$cluster)[as.numeric(protsCluster$cluster)]
protsCluster <- protsCluster[order(as.numeric(protsCluster$idx)),]

Background <- cbind(toHMP,protsCluster)
Background$idx <- NULL; Background$proteins <- NULL
saveTablesTsvExc(Background,"FinalData/BackgroundRemoval",completeNdedup=F,excel=T,bycompleteFC=F,rownames=T)

rowannot <- rowAnnotation(Cluster=protsCluster$nProtsInClust, show_annotation_name = F, show_legend = F)
myHMPplot <- Heatmap(as.matrix(toHMP),cluster_columns = F, cluster_rows = T, col=col_fun,
                     show_row_names = F,show_row_dend = F, row_split = as.character(protsCluster$nProtsInClust),
                     name = "Log10(LFQ)", #top_annotation = colannot, 
                     left_annotation = rowannot, row_title_rot = 0,
                     column_names_rot = 0, column_split = subPhenData$class,
                     heatmap_legend_param = list(
                       legend_direction = "horizontal", 
                       legend_width = unit(5, "cm")
                     ), column_names_side = "top", row_gap = unit(10, "mm")
)

tiff(file.path(outdir,paste0(gsub("\ ", "",tagname),'HeatMapClustered.tiff')), 1000, nrow(toHMP)*5, pointsize=12,res = "print")
draw(myHMPplot, heatmap_legend_side="bottom", annotation_legend_side = "left")
dev.off()

#### 3) Table with Non Normalised Values from RBPome detected and the background control.

NOXrbpomeFile <- '/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/FinalData/Raw/correctedRAW/NOXmRBPomeRAW.csv'
FAXrbpomeFile <- '/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/FinalData/Raw/correctedRAW/FAXmRBPomeRAW.csv'
UVXrbpomeFile <- '/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/FinalData/Raw/correctedRAW/UVXmRBPomeRAW.csv'

NOXrbpomeDF <- getProteinsMatrixFromRAW(NOXrbpomeFile)
FAXrbpomeDF <- getProteinsMatrixFromRAW(FAXrbpomeFile)
UVXrbpomeDF <- getProteinsMatrixFromRAW(UVXrbpomeFile)

colnames(NOXrbpomeDF) <- paste(colnames(NOXrbpomeDF),'NOXraw', sep = '_')
colnames(NOXrbpomeDF)[1] <- 'proteins'
colnames(FAXrbpomeDF) <- paste(colnames(FAXrbpomeDF),'FAXraw', sep = '_')
colnames(FAXrbpomeDF)[1] <- 'proteins'
colnames(UVXrbpomeDF) <- paste(colnames(UVXrbpomeDF),'UVXraw', sep = '_')
colnames(UVXrbpomeDF)[1] <- 'proteins'

FAXnUVXdf <- merge(FAXrbpomeDF,UVXrbpomeDF,by = "proteins",all = T)
NOXnFAXnUVXdf <- merge(NOXrbpomeDF,FAXnUVXdf,by = "proteins",all = T)

meanNOXfaxFile <- "FinalData/Raw/correctedRAW/meanRBPomeNOX_FAX_np2_NOnorm/meanRBPomeNOX_FAX_np2_NOnorm_PROTEIN.tsv"
meanNOXuvxFile <- "FinalData/Raw/correctedRAW/meanRBPomeNOX_UVX_np2_NOnorm/meanRBPomeNOX_UVX_np2_NOnorm_PROTEIN.tsv"
meanNOXfaxDF <- read.delim(meanNOXfaxFile)[c('proteinName','log2ratio_Condition2')]
colnames(meanNOXfaxDF) <- c('proteins','log2ratio_FAX')
meanNOXuvxDF <- read.delim(meanNOXuvxFile)[c('proteinName','log2ratio_Condition2')]
colnames(meanNOXuvxDF) <- c('proteins','log2ratio_UVX')

NOXnFAXnUVX_withFAXlog2 <- merge(NOXnFAXnUVXdf,meanNOXfaxDF,all = T)
NOXnFAXnUVX_withFAXnUVXlog2 <- merge(NOXnFAXnUVX_withFAXlog2,meanNOXuvxDF,all = T)

table(gsub('.x|.y','',colnames(NOXnFAXnUVX_withFAXnUVXlog2)))

outdir <- 'FinalData/BackgroundRemoval'
saveTablesTsvExc(NOXnFAXnUVX_withFAXnUVXlog2,outdir,completeNdedup=F,excel=T,bycompleteFC=F,rownames=F)


#### MAP IT
NOXnFAXnUVX_withFAXnUVXlog2 <- read.delim('FinalData/BackgroundRemoval/NOXnFAXnUVX_withFAXnUVXlog2.tsv',quote = '')

mappingFinalFile <- "yeastReference/mappingFile.xlsx"
yeastGenesProtMap <- read.xlsx(mappingFinalFile)
yeastGenesProtMap[which(yeastGenesProtMap$GeneName == 'REP2'),]$UniprotACC <- 'P03872'

FAXnUVXdf <- NOXnFAXnUVX_withFAXnUVXlog2 
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

NOXnFAXnUVX_withFAXnUVXlog2Mapped <- FAXnUVXfullMapped
outdir <- 'FinalData/BackgroundRemoval'
NOXnFAXnUVX_withFAXnUVXlog2Mapped <- NOXnFAXnUVX_withFAXnUVXlog2Mapped[!duplicated(NOXnFAXnUVX_withFAXnUVXlog2Mapped),]
saveTablesTsvExc(NOXnFAXnUVX_withFAXnUVXlog2Mapped,outdir,completeNdedup=F,excel=T,bycompleteFC=F,rownames=F)

