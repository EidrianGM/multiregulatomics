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

medianNOXfaxFile <- "FinalData/Raw/correctedRAW/byPEPmedianRBPomeNOX_FAX_np2_NOnorm/byPEPmedianRBPomeNOX_FAX_np2_NOnorm_PROTEIN.tsv"
medianNOXuvxFile <- "FinalData/Raw/correctedRAW/byPEPmedianRBPomeNOX_UVX_np2_NOnorm/byPEPmedianRBPomeNOX_UVX_np2_NOnorm_PROTEIN.tsv"
meanNOXfaxFile <- "FinalData/Raw/correctedRAW/byPEPmeanRBPomeNOX_FAX_np2_NOnorm/byPEPmeanRBPomeNOX_FAX_np2_NOnorm_PROTEIN.tsv"
meanNOXuvxFile <- "FinalData/Raw/correctedRAW/byPEPmeanRBPomeNOX_UVX_np2_NOnorm/byPEPmeanRBPomeNOX_UVX_np2_NOnorm_PROTEIN.tsv"
outFile <- "FinalData/BackgroundRemoval/ByPEPstatsBGremoval.tsv"
medianNOXfaxDF <- read.delim(medianNOXfaxFile);
medianNOXuvxDF <- read.delim(medianNOXuvxFile);
meanNOXfaxDF <- read.delim(meanNOXfaxFile);
meanNOXuvxDF <- read.delim(meanNOXuvxFile);


###### mRBPomes UVX and FAX filtering

FAXfile <- "FinalData/SafeQuantResults/RBPomeFAX_np2_norm_gmin_no12/RBPomeFAX_np2_norm_gmin_no12_PROTEIN.tsv"
UVXfile <- "FinalData/SafeQuantResults/RBPomeUVX_np2_norm_gmin_no10-7-16-3/RBPomeUVX_np2_norm_gmin_no10-7-16-3_PROTEIN.tsv"

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
length(unique(meanNOXfaxDF$proteinName))
length(abovebackgrProtsFAX)
length(belowbackgrProtsFAX)

nrow(meanNOXuvxDF)
length(unique(meanNOXuvxDF$proteinName))
length(abovebackgrProtsUVX)
length(belowbackgrProtsUVX)



##### GENERATE HEATMAPS

FAXrbpomeFile <- "FinalData/SafeQuantResults/RBPomeFAX_np2_norm_gmin_no12/RBPomeFAX_np2_norm_gmin_no12_PROTEIN.tsv"
UVXrbpomeFile <- "FinalData/SafeQuantResults/RBPomeUVX_np2_norm_gmin_no10-7-16-3/RBPomeUVX_np2_norm_gmin_no10-7-16-3_PROTEIN.tsv"

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

meanNOXfaxFile <- "FinalData/Raw/correctedRAW/meanRBPomeNOX_FAX_np2_NOnorm/meanRBPomeNOX_FAX_np2_NOnorm_PROTEIN.tsv"
meanNOXuvxFile <- "FinalData/Raw/correctedRAW/meanRBPomeNOX_UVX_np2_NOnorm/meanRBPomeNOX_UVX_np2_NOnorm_PROTEIN.tsv"
meanNOXfaxDF <- read.delim(meanNOXfaxFile)
meanNOXuvxDF <- read.delim(meanNOXuvxFile)

FCcutOff <- 3
meanNOXfaxDF <- meanNOXfaxDF[which(abs(meanNOXfaxDF$log2ratio_Condition2) > FCcutOff),]
meanNOXuvxDF <- meanNOXuvxDF[which(abs(meanNOXuvxDF$log2ratio_Condition2) > FCcutOff),]

meanNOXfaxDF <- meanNOXfaxDF[order(meanNOXfaxDF$log2ratio_Condition2,decreasing = T),]
meanNOXuvxDF <- meanNOXuvxDF[order(meanNOXuvxDF$log2ratio_Condition2,decreasing = T),]

rownames(meanNOXfaxDF) <- meanNOXfaxDF$proteinName
meanNOXfaxDF <- meanNOXfaxDF[abovebackgrProtsFAX,grep("^X",colnames(meanNOXfaxDF))]
colnames(meanNOXfaxDF) <- as.numeric(gsub("X|_.*","",colnames(meanNOXfaxDF)))

rownames(meanNOXuvxDF) <- meanNOXuvxDF$proteinName
meanNOXuvxDF <- meanNOXuvxDF[abovebackgrProtsUVX,grep("^X",colnames(meanNOXuvxDF))]
colnames(meanNOXuvxDF) <- as.numeric(gsub("X|_.*","",colnames(meanNOXuvxDF)))

library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)

outdir <- "FinalData/BackgroundRemoval"
tagname <- "Background"

FAXnUVXbackgroundDF <- merge(meanNOXfaxDF, meanNOXuvxDF,by="row.names",all=T,suffixes = c("","y"))
row.names(FAXnUVXbackgroundDF) <- FAXnUVXbackgroundDF$Row.names
FAXnUVXbackgroundDF$Row.names <- NULL
FAXnUVXbackgroundDF <- FAXnUVXbackgroundDF[colnames(FAXnUVXbackgroundDF)[grep("y",colnames(FAXnUVXbackgroundDF),invert = T)]]
#toHMP <- log10(FAXnUVXbackgroundDF+1)
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

backgroundHMAP <- cbind(toHMP,protsCluster)
backgroundHMAP$idx <- NULL; backgroundHMAP$proteins <- NULL
# outFile <- gsub(" ", "", file.path("FinalData/BackgroundRemoval",tagname), fixed = TRUE)
# write.table(backgroundHMAP,paste0(outFile,".tsv"),row.names = T,col.names = T,quote = F, sep = "\t")
saveTablesTsvExc(backgroundHMAP,"FinalData/BackgroundRemoval",completeNdedup=F,excel=T,bycompleteFC=F,rownames=T)

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

nrow(toHMP)

tiff(file.path(outdir,paste0(gsub("\ ", "",tagname),'HeatMapClustered.tiff')), 1000, nrow(toHMP)*5, pointsize=12,res = "print")
draw(myHMPplot, heatmap_legend_side="bottom", annotation_legend_side = "left")
dev.off()



# crosslink <- "FAX"
# toHMP <- log10(meanNOXfaxDF+1)
crosslink <- "UVX"
toHMP <- log10(meanNOXuvxDF+1)

tagname <- paste("Background",crosslink)
subPhenData2 <- subPhenData[subPhenData$class %in% c("NOX",crosslink),]
subPhenData2$class <- factor(subPhenData2$class,levels = c("NOX",crosslink))
toHMP <- toHMP[,as.character(subPhenData2$id)]

colours <- brewer.pal(brewer.pal.info["PuRd",]$maxcolors,"PuRd")
colourStep <- (max(toHMP,na.rm = T) / brewer.pal.info["PuRd",]$maxcolors)
breaksColours <- seq(0,max(toHMP,na.rm = T)-colourStep,colourStep)
col_fun <- colorRamp2(breaksColours, colours) # #574a85
myHMPplot <- Heatmap(as.matrix(toHMP),cluster_columns = F, cluster_rows = T, col=col_fun,
                     show_row_names = F,show_row_dend = F, row_km=8,
                     name = "Log10(LFQ)", #top_annotation = colannot,
                     column_names_rot = 0, column_split = subPhenData2$class,
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

backgroundHMAP <- cbind(toHMP,protsCluster)
backgroundHMAP$idx <- NULL; backgroundHMAP$proteins <- NULL
outFile <- gsub(" ", "", file.path("FinalData/BackgroundRemoval",tagname), fixed = TRUE)
write.table(backgroundHMAP,paste0(outFile,".tsv"),row.names = F,col.names = F,quote = F, sep = "\t")

rowannot <- rowAnnotation(Cluster=protsCluster$nProtsInClust, show_annotation_name = F, show_legend = F)

myHMPplot <- Heatmap(as.matrix(toHMP),cluster_columns = F, cluster_rows = T, col=col_fun,
                     show_row_names = F,show_row_dend = F, row_split = as.character(protsCluster$nProtsInClust),
                     name = "Log10(LFQ)", #top_annotation = colannot, 
                     left_annotation = rowannot, row_title_rot = 0,
                     column_names_rot = 0, column_split = subPhenData2$class,
                     heatmap_legend_param = list(
                       legend_direction = "horizontal", 
                       legend_width = unit(5, "cm")
                     ), column_names_side = "top", row_gap = unit(10, "mm")
)

tiff(file.path(outdir,paste0(gsub("\ ", "",tagname),'HeatMapClustered.tiff')), 1000, nrow(toHMP)*5, pointsize=12,res = "print")
draw(myHMPplot, heatmap_legend_side="bottom", annotation_legend_side = "left")
dev.off()

