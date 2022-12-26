library(paletteer)
library(ggplot2)

getPlotsRawData <- function(inFile,tagname){
  myDF <- read.delim(inFile,sep = ",")
  outdir <- "/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/adrian/rawViz"
  normDCol <- grep("Norm",colnames(myDF))
  rawDCol <- grep("Raw",colnames(myDF))
  expresionMatrix <- myDF[normDCol:(rawDCol-1)]
  someVisualization(expresionMatrix,outdir,paste0("norm",tagname))
  #rawData <- myDF[rawDCol:(rawDCol+ncol(normData)-1)]
  #someVisualization(rawData,outdir,paste0("raw",tagname))
}

someVisualization <- function(expresionMatrix,outdir,tagname,subPhenData=NA){
  
  if (class(subPhenData) == "logical"){
    expresionMatrix[2,] <- gsub("_.*","",expresionMatrix[2,])
    subPhenData <- as.data.frame(t(expresionMatrix[1:2,])); colnames(subPhenData) <- c("class","id"); rownames(subPhenData) <- subPhenData[,2]
    colnames(expresionMatrix) <- expresionMatrix[2,]; expresionMatrix <- expresionMatrix[3:nrow(expresionMatrix),]; expresionMatrix <- type.convert(expresionMatrix)
    subPhenData$class <- gsub(".*\\s","",subPhenData$class)
  }
  
  conditionCol <- "class"
  classes <- levels(as.factor(subPhenData[,conditionCol]))
  colors <- paletteer_dynamic("cartography::multi.pal", length(classes))[1:length(classes)]; names(colors) <- classes
  allColors <- colors[match(subPhenData[,conditionCol],names(colors))]
  colorsMatch <- as.list(allColors); names(colorsMatch) <- subPhenData[,conditionCol]
  
  totalSignal <- apply(expresionMatrix, 2, sum,na.rm=T)
  totalSignal <- totalSignal / max(totalSignal) * 100  
  totalSignal <- cbind(subPhenData,signalpercent=totalSignal)
  totalSignal$class <- factor(totalSignal$class, levels = unique(totalSignal$class))
  totalSignal$id <- factor(totalSignal$id, levels = unique(totalSignal$id))
  myPlot <- ggplot(totalSignal, aes(x=id, y=signalpercent, fill=class)) +
    geom_bar(stat="identity")+theme_classic(base_size = 10)+scale_fill_manual(values = colorsMatch) + 
    geom_hline(yintercept = 50)
  ggsave(file.path(outdir,paste0(tagname,'_SignalPlot.tiff')),myPlot,
         width = length(subPhenData$id) * 100,device="tiff",units = "px",dpi=300)
    
    
  toBoxPlot <- melt(log2(expresionMatrix))
  toBoxPlot$Var2 <- factor(as.character(toBoxPlot$Var2),levels=unique(toBoxPlot$Var2))
  toBoxPlot$class <- subPhenData$class[match(toBoxPlot$Var2, as.character(as.numeric(subPhenData$id)))]
  myBox <- ggplot(toBoxPlot, aes(x=Var2, y=value, fill=class)) + geom_boxplot() +
    scale_fill_manual(values = unique(allColors),aesthetics = "fill") +
    theme_classic() + theme(legend.position="right",plot.title = element_text(size=11)) +
    ggtitle(gsub("_"," ",tagname)) + xlab("")
  ggsave(file.path(outdir,paste0(tagname,'_BoxPlot.tiff')),myBox,device="tiff",units = "px",dpi=300)
  
  sum(SamplesMissingInGenes == ncol(expresionMatrix))
  
  if (any(is.na(expresionMatrix))){
    SamplesMissingInGenes <- apply(is.na(expresionMatrix),1,sum)
    GenesMissingInSamples <- round(apply(is.na(expresionMatrix),2,sum) / nrow(expresionMatrix) * 100,2)
    column_ha = HeatmapAnnotation(GenesMissing = anno_barplot(GenesMissingInSamples,add_numbers = T,
                                                              numbers_rot=0,numbers_offset = unit(1, "mm"),
                                                              numbers_gp=gpar(fontsize = 6)),
                                  treatment = subPhenData[,conditionCol],col = colorsMatch)
    row_ha = rowAnnotation(SamplesMissing = anno_barplot(SamplesMissingInGenes[1:15000]))
    toHMP <- t(scale(t(expresionMatrix)))[1:15000,]
    col_fun <- colorRamp2(c(min(toHMP,na.rm = T),0,max(toHMP,na.rm = T)), c("blue","white","red"))
    tiff(file.path(outdir,paste0(tagname,'_NAs_HeatMap.tiff')), ncol(toHMP) * 120, 1500, pointsize=5, res = 300)
    print(
      Heatmap(toHMP, cluster_columns=F, cluster_rows=F, col=col_fun,
              top_annotation=column_ha, right_annotation=row_ha,
              show_row_names=F, show_row_dend=F, na_col="grey",name = paste0(tagname,'_NAs'),
              show_heatmap_legend=F,  column_title_gp=gpar(fontsize = 6),
              column_names_gp = gpar(fontsize = 6),
              row_names_gp = gpar(fontsize = 6))
      #, width = ncol(toHMP)*unit(10, "mm")
      #, height = nrow(toHMAP)*unit(5, "mm")) 
    )
    dev.off()
  }
  
  expresionMatrix <- na.omit(expresionMatrix)
  
  pcares <- prcomp(t(expresionMatrix))
  tiff(file.path(outdir,paste0(tagname,'_PCA2D.tiff')), 1500, 1500, pointsize=12, res = 300)
  print(ggplot2::autoplot(pcares,data = data.frame(t(expresionMatrix),class=subPhenData[,conditionCol]),
                          colour="class", label=T, label.size=2,label.repel=T,main = gsub("_"," ",tagname)) + 
          scale_color_manual(values = unique(as.vector(colors)))+ theme_classic() + theme(legend.title=element_blank())+
          guides(color = guide_legend(override.aes = aes(label = ""))))
  dev.off()
  
  # pcaVarGroup <- round((pcares$sdev)^2 / sum(pcares$sdev^2) *100, 2)
  # tiff(file.path(outdir,paste0(tagname,'_PCA3D.tiff')), 1200, 900, pointsize=20, res = 300)
  # pca3d <- scatterplot3d(pcares$x[,1], pcares$x[,3], pcares$x[,2],color = as.vector(allColors),
  #                        main = tagname, xlab = paste0("PC 1 (", pcaVarGroup[1], " %)"),
  #                        zlab = paste0("PC 2 (", pcaVarGroup[2], " %)"), ylab = paste0("PC 3 (", pcaVarGroup[3], " %)"),
  #                        grid=T, box = T, pch = 20, cex.symbols = 2.5, angle = 40,type = "h")
  # zz.coords <- pca3d$xyz.convert(pcares$x[,1], pcares$x[,3], pcares$x[,2])
  # legend("topright", pch=20, legend = classes, col = unique(as.vector(allColors)),inset = 0,y.intersp =0.8)
  # text(zz.coords$x, zz.coords$y, labels = row.names(pcares$x),cex = 0.5, pos = 4,col = as.vector(allColors))
  # dev.off()
  
  toHMP <- t(scale(t(expresionMatrix)))
  col_fun <- colorRamp2(c(min(toHMP,na.rm = T),0,max(toHMP,na.rm = T)), c("blue","white","red"))
  column_ha <- HeatmapAnnotation(treatment = subPhenData[,conditionCol],col = colorsMatch)
  tiff(file.path(outdir,paste0(tagname,'_HeatMap.tiff')), 1200, 900, pointsize=20, res = 300)
  print(Heatmap(toHMP,cluster_columns = F, cluster_rows = T,col=col_fun,
                top_annotation = column_ha,show_row_names = F,show_row_dend = F,na_col = "grey"),
        column_title_gp=gpar(fontsize = 6),
        column_names_gp = gpar(fontsize = 6),
        row_names_gp = gpar(fontsize = 6))
  dev.off()
  
  corMat <- cor(expresionMatrix,method = "pearson")
  
  col_fun <- colorRamp2(c(min(toHMP,na.rm = T),0,max(toHMP,na.rm = T)), c("blue","white","red"))
  tiff(file.path(outdir,paste0(tagname,'_CorrPlot.tiff')), ncol(toHMP) * 80, ncol(toHMP) * 80, pointsize=10, res = 300)
  corrplot(corMat, method = 'circle', type = 'upper', order = 'original', tl.cex=1, title=tagname,
           tl.col = allColors,tl.srt = 45,cl.ratio = 0.2, number.cex = 0.6,
           mar = c(5, 1, 2, 1),tl.pos = 'lt',addCoef.col = "grey",
           col=rev(COL2('RdBu', 200)))
  print(legend("bottom", legend=classes,fill=colors, horiz=T, cex=0.5,y = 0, bty = "n"))
  dev.off()
}


################################################################################
############################ Vizs of RAW #######################################
################################################################################

#### UV RBPome
inFile <- "/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/data/raw/MQtoSQ_UV.csv"
tagname <- "UV-RBPome"
getPlotsRawData(inFile,tagname)

#### FAX RBPome
inFile <- "/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/data/raw/MQtoSQ_FAX.csv"
tagname <- "FAX-RBPome"
getPlotsRawData(inFile,tagname)

#### NOX RBPome
inFile <- "/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/data/raw/MQtoSQ.csv"
tagname <- "NOX-RBPome"
getPlotsRawData(inFile,tagname)

#### Proteome
nFile <- "/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/data/raw/SNtoSQ-Prot40_SafeQuant_Input.csv"
tagname <- "Proteome"
getPlotsRawData(inFile,tagname)

################################################################################
############################ Imputation ########################################
################################################################################

#### UV RBPome
inFile <- "/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/data/raw/MQtoSQ_UV.csv"
myDF <- read.csv(inFile)

normDCol <- grep("Norm",colnames(myDF))
rawDCol <- grep("Raw",colnames(myDF))
expresionMatrix <- myDF[normDCol:(rawDCol-1)]
expresionMatrix[2,] <- gsub("_.*","",expresionMatrix[2,])
subPhenData <- as.data.frame(t(expresionMatrix[1:2,])); colnames(subPhenData) <- c("class","id"); rownames(subPhenData) <- subPhenData[,2]
colnames(expresionMatrix) <- expresionMatrix[2,]; expresionMatrix <- expresionMatrix[3:nrow(expresionMatrix),]; expresionMatrix <- type.convert(expresionMatrix)
subPhenData$class <- gsub(".*\\s","",subPhenData$class)

eset <- sqNormalize(eset, method=method)

if("global" %in% method){
  globalNormFactors <- getGlobalNormFactors(esetNorm,method=method)
  ### add normalization factors to ExpressionSet
  pData(esetNorm)$globalNormFactors <- globalNormFactors
  esetNorm <- globalNormalize(esetNorm,globalNormFactors)
}
tagname <- "UV-RBPome"


#### FAX RBPome
inFile <- "/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/data/raw/MQtoSQ_FAX.csv"
tagname <- "FAX-RBPome"
getPlotsRawData(inFile,tagname)

#### NOX RBPome
inFile <- "/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/data/raw/MQtoSQ.csv"
tagname <- "NOX-RBPome"
getPlotsRawData(inFile,tagname)

#### Proteome
nFile <- "/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/data/raw/SNtoSQ-Prot40_SafeQuant_Input.csv"
tagname <- "Proteome"
getPlotsRawData(inFile,tagname)

################################################################################
############################ Vizs of Imputed ###################################
################################################################################


################################################################################
############################ What's DONE? ######################################
################################################################################

compareRawNDone <- function(doneFile,rawFile,crosslink=NA){
  doneDF <- read.delim(doneFile,nrows = 3)
  mysamplesDone <- colnames(doneDF)[grep("^X",colnames(doneDF))]
  mysamplesDone <- unique(gsub("X0|_.*","",mysamplesDone))
  
  rawDF <- read.csv(rawFile,nrows = 3)
  if (!is.na(crosslink)){
    rawDF <- rawDF[,grep(crosslink,rawDF[1,])]
    mysamples <- rawDF[2,]
  }else{
    mysamples <- rawDF[2,][grep("PolyA",rawDF[2,])]
  }
  mysamples <- unique(gsub("^0|_.*","",mysamples))
  if (length(mysamples) == length(mysamplesDone) &  all(mysamples %in% mysamplesDone)){
    cat("All Okey")
  }else{
    cat("Samples Discarded: ")
    cat(paste0(setdiff(mysamples,mysamplesDone),collapse = " "))
  }
}

### NOX RBPome
rawFile <- "/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/data/raw/MQtoSQ.csv"
doneFile <- "/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/data/old/P445_AG_PolyA_40_all_3_cond_20220611_Final/P445_AG_PolyA_40_20220617_gmin_UV_nonorm_p2_NOX_Final/P445_AG_PolyA_40_20220617_gmin_UV_nonorm_p2_NOX_PROTEIN.tsv"
compareRawNDone(doneFile,rawFile)

doneFile <- "/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/data/new/PolyA_RIC_nonorm_20220928/P445_AG_PolyA_40_20220617_gmin_UV_nonorm_p2_NOX_Final/P445_AG_PolyA_40_20220617_gmin_UV_nonorm_p2_NOX_PROTEIN.tsv"
compareRawNDone(doneFile,rawFile)

### UVX RBPome
rawFile <- "/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/data/raw/MQtoSQ_UV.csv"
doneFile <- "/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/data/old/P445_AG_PolyA_40_all_3_cond_20220611_Final/P445_AG_PolyA_40_20220617_gmin_UV_p2_no511_Final/P445_AG_PolyA_40_20220617_gmin_UV_p2_no511_PROTEIN.tsv"
compareRawNDone(doneFile,rawFile)

doneFile <- "/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/data/new/PolyA_RIC_nonorm_20220928/P445_AG_PolyA_40_20220617_gmin_UV_p2_nonorm_no511_Final/P445_AG_PolyA_40_20220617_gmin_UV_p2_nonorm_no511_Final_PROTEIN.tsv"
compareRawNDone(doneFile,rawFile)
### FAX RBPome
rawFile <- "/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/data/raw/MQtoSQ_FAX.csv"
doneFile <- "/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/data/old/P445_AG_PolyA_40_all_3_cond_20220611_Final/P445_AG_PolyA_40_20220617_gmin_FAX_p2_Final/P445_AG_PolyA_40_20220617_gmin_FAX_p2_PROTEIN.tsv"
compareRawNDone(doneFile,rawFile)

doneFile <- "/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/data/new/PolyA_RIC_nonorm_20220928/P445_AG_PolyA_40_20220617_gmin_FAX_p2_nonorm_no12_Final/P445_AG_PolyA_40_20220617_gmin_FAX_p2_nonorm_no12_Final_PROTEIN.tsv"
compareRawNDone(doneFile,rawFile)

### Proteome
rawFile <- "/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/data/raw/SNtoSQ-Prot40_SafeQuant_Input.csv"

### NOX Proteome
crosslink <- "NOX"

doneFile <- "/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/data/new/P445_Prot_DIA_iBAQ_nocomb_20220826_all_NOX_Final/P445_Prot_DIA_iBAQ_noOutliers_nocomb_20220826_NOX_Final/P445_Prot_DIA_iBAQ_noOutliers_nocomb_20220826_NOX_Final_PROTEIN.tsv"
compareRawNDone(doneFile,rawFile,crosslink)

doneFile <- "/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/data/new/P445_Prot_DIA_iBAQ_nocomb_20220826_all_NOX_Final/P445_Prot_DIA_iBAQ_nocomb_20220826_NOX_Final2/P445_Prot_DIA_iBAQ_nocomb_20220826_NOX_Final2_PROTEIN.tsv"
compareRawNDone(doneFile,rawFile,crosslink)

doneFile <- "/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/data/old/P445_Prot_DIA_iBAQ_noOutliers_nocomb_20220826_Final/P445_Prot_DIA_iBAQ_noOutliers_nocomb_20220826_NOX_Final/P445_Prot_DIA_iBAQ_noOutliers_nocomb_20220826_NOX_Final_PROTEIN.tsv"
compareRawNDone(doneFile,rawFile,crosslink)

### UVX Proteome
crosslink <- "UV"
doneFile <- "/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/data/old/P445_Prot_DIA_iBAQ_noOutliers_nocomb_20220826_Final/P445_Prot_DIA_iBAQ_noOutliers_nocomb_20220826_UV_Final/P445_Prot_DIA_iBAQ_noOutliers_nocomb_20220826_UV_Final_PROTEIN.tsv"
compareRawNDone(doneFile,rawFile,crosslink)

### FAX Proteome
crosslink <- "FAX"
doneFile <- "/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/data/old/P445_Prot_DIA_iBAQ_noOutliers_nocomb_20220826_Final/P445_Prot_DIA_iBAQ_noOutliers_nocomb_20220826_FAX_Final/P445_Prot_DIA_iBAQ_noOutliers_nocomb_20220826_FAX_Final_PROTEIN.tsv"
compareRawNDone(doneFile,rawFile,crosslink)


################################################################################
############################ Comparing NOX vs FAX and UVX ######################
################################################################################

inFile <- "/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/data/raw/MQtoSQ.csv"
myDF <- read.csv(inFile)

normDCol <- grep("Norm",colnames(myDF))
rawDCol <- grep("Raw",colnames(myDF))
expresionMatrix <- myDF[normDCol:(rawDCol-1)]
expresionMatrix[2,] <- gsub("_.*","",expresionMatrix[2,])
samplesInNOXraw <- as.character(expresionMatrix[2,])
subPhenDataNOX <- as.data.frame(t(expresionMatrix[1:2,])); colnames(subPhenData) <- c("class","id"); rownames(subPhenData) <- subPhenData[,2]


UVXresultsFile <- "/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/FinalData/SafeQuantResults/RBPomeUVX_np2_norm_gmin_no10-7-16-3/RBPomeUVX_np2_norm_gmin_no10-7-16-3_PROTEIN.tsv"
UVXresultsDF <- read.delim(UVXresultsFile)
UVXresultsSamples <- colnames(UVXresultsDF)[grep("^X",colnames(UVXresultsDF))]
UVXresultsSamples <- gsub("^X","",gsub("_.*","",UVXresultsSamples))


FAXresultsFile <- "/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/FinalData/SafeQuantResults/RBPomeFAX_np2_norm_gmin_no12/RBPomeFAX_np2_norm_gmin_no12_PROTEIN.tsv"
FAXresultsDF <- read.delim(FAXresultsFile)
FAXresultsSamples <- colnames(FAXresultsDF)[grep("^X",colnames(FAXresultsDF))]
FAXresultsSamples <- gsub("^X","",gsub("_.*","",FAXresultsSamples))

samplesInResults <- c(UVXresultsSamples,FAXresultsSamples)
all(samplesInResults %in% samplesInNOXraw) ## TRUE

cat(intersect(samplesInNOXraw,FAXresultsSamples),sep = "\n")
cat(intersect(samplesInNOXraw,UVXresultsSamples),sep = "\n")

cat(setdiff(samplesInNOXraw,samplesInResults),sep = "\n")
cat(samplesInNOXraw,sep = "\n")

## Samples in RBPome
rawRBPomeFAXfile <- "/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/FinalData/Raw/MQtoSQ_FAX.csv"
rawRBPomeFAXSamplesInfo <- read.delim(rawRBPomeFAXfile,nrows = 3, sep = ",")

rawRBPomeUVXfile <- "/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/FinalData/Raw/MQtoSQ_UV.csv"
rawRBPomeUVXSamplesInfo <- read.delim(rawRBPomeUVXfile,nrows = 3, sep = ",")

rawRBPomeNOXfile <- "/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/FinalData/Raw/MQtoSQ.csv"
rawRBPomeNOXSamplesInfo <- read.delim(rawRBPomeNOXfile,nrows = 3, sep = ",")

phenoDataFAX <- as.data.frame(t(rawRBPomeFAXSamplesInfo[1:2,grep("Norm",colnames(rawRBPomeFAXSamplesInfo)):(grep("Raw",colnames(rawRBPomeFAXSamplesInfo))-1)]))
phenoDataUVX <- as.data.frame(t(rawRBPomeUVXSamplesInfo[1:2,grep("Norm",colnames(rawRBPomeUVXSamplesInfo)):(grep("Raw",colnames(rawRBPomeUVXSamplesInfo))-1)]))
phenoDataNOX <- as.data.frame(t(rawRBPomeNOXSamplesInfo[1:2,grep("Norm",colnames(rawRBPomeNOXSamplesInfo)):(grep("Raw",colnames(rawRBPomeNOXSamplesInfo))-1)]))
colnames(phenoDataFAX) <- c("classFAX","sample")
colnames(phenoDataUVX) <- c("classUVX","sample")
colnames(phenoDataNOX) <- c("classNOX","sample")

phenoDataFAX$classFAX <- gsub(".*\\s","",phenoDataFAX$class)
phenoDataUVX$classUVX <- gsub(".*\\s","",phenoDataUVX$class)
phenoDataNOX$classNOX <- gsub(".*\\s","",phenoDataNOX$class)

allPhenoD <- merge(phenoDataNOX,phenoDataUVX,by = "sample",all = T)
allPhenoD <- merge(allPhenoD,phenoDataFAX,by = "sample",all = T)

###### 
proteomeFile <- "/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/FinalData/Raw/correctedRAW/newProteome.csv"
proteomeDF <- read.delim(proteomeFile, sep = ",")
length(proteomeDF[,grep("Acc",proteomeDF[2,])]) - 2
length(unique(proteomeDF[,grep("Acc",proteomeDF[2,])])) - 2

rawRBPomeFAXfile <- "/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/FinalData/Raw/correctedRAW/FAXmRBPomeRAW.csv"
rawRBPomeFAXDF <- read.delim(rawRBPomeFAXfile, sep = ",")
length(rawRBPomeFAXDF[,grep("Acc",rawRBPomeFAXDF[2,])]) - 2
length(unique(rawRBPomeFAXDF[,grep("Acc",rawRBPomeFAXDF[2,])])) - 2

rawRBPomeUVXfile <- "/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/FinalData/Raw/correctedRAW/UVXmRBPomeRAW.csv"
rawRBPomeUVXDF <- read.delim(rawRBPomeUVXfile, sep = ",")
length(rawRBPomeUVXDF[,grep("Acc",rawRBPomeUVXDF[2,])]) - 2
length(unique(rawRBPomeUVXDF[,grep("Acc",rawRBPomeUVXDF[2,])])) - 2

rawRBPomeNOXfile <- "/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/FinalData/Raw/correctedRAW/NOXmRBPomeRAW.csv"
rawRBPomeNOXDF <- read.delim(rawRBPomeNOXfile, sep = ",")
length(rawRBPomeNOXDF[,grep("Acc",rawRBPomeNOXDF[2,])]) - 2
length(unique(rawRBPomeNOXDF[,grep("Acc",rawRBPomeNOXDF[2,])])) - 2

