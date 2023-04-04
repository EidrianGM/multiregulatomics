library(paletteer)
library(ggplot2)
library(corrplot)
library(reshape2)
library(circlize)
library(ggfortify)
library(ComplexHeatmap)

setwd('/home/eidrian/Desktop/realServices/multiregulatomics')

getPlotsRawData <- function(inFile,tagname){
  myDF <- read.delim(inFile,sep = ",")
  outdir <- "FinalData/QualityControlPlots"
  normDCol <- grep("Norm",colnames(myDF))
  rawDCol <- grep("Raw",colnames(myDF))
  peptidesInfo <- myDF[1:(normDCol-1)]
  colnames(peptidesInfo) <- peptidesInfo[2,]
  peptidesInfo <- peptidesInfo[3:nrow(peptidesInfo),]
  expresionMatrix <- myDF[normDCol:(rawDCol-1)]
  someVisualization(expresionMatrix,outdir,paste0("norm",tagname))
  #rawData <- myDF[rawDCol:(rawDCol+ncol(normData)-1)]
  #someVisualization(rawData,outdir,paste0("raw",tagname))
}

someVisualization <- function(inFile,tagname,subPhenData=NA,peptidesInfo=NA){
  extension <- '.png'
  mydevice <- png
  
  outdir <- file.path("FinalData/QualityControlPlots",tagname)
  dir.create(outdir,showWarnings = F)
  logFile <- file.path(outdir,paste0(tagname,'.log'))
  
  if (class(inFile) == 'data.frame'){
    expresionMatrix <- inFile
  }
  
  if (class(subPhenData) == "logical"){
    myDF <- read.delim(inFile,sep = ",")
    normDCol <- grep("Norm",colnames(myDF))
    rawDCol <- grep("Raw",colnames(myDF))
    peptidesInfo <- myDF[1:(normDCol-1)]
    colnames(peptidesInfo) <- peptidesInfo[2,]
    peptidesInfo <- peptidesInfo[3:nrow(peptidesInfo),]
    expresionMatrix <- myDF[normDCol:(rawDCol-1)]
    expresionMatrix[2,] <- gsub("_.*","",expresionMatrix[2,])
    subPhenData <- as.data.frame(t(expresionMatrix[1:2,])); colnames(subPhenData) <- c("class","id"); rownames(subPhenData) <- subPhenData[,2]
    subPhenData$id <- as.character(as.numeric(subPhenData$id))
    colnames(expresionMatrix) <- subPhenData$id; 
    expresionMatrix <- expresionMatrix[3:nrow(expresionMatrix),]; 
    expresionMatrix <- type.convert(expresionMatrix)
    subPhenData$class <- gsub(".*\\s","",subPhenData$class)
  }
  
  if (tagname == "Proteome"){
    crosslinks <- unique(gsub("with.*|no.*",'',subPhenData$class))
    for (crosslink in crosslinks){
      subsubPhenData <- subPhenData[grep(crosslink,subPhenData$class),]
      subsetExpMat <- expresionMatrix[subsubPhenData$id]
      newtag <- paste0(crosslink,tagname)
      someVisualization(subsetExpMat,newtag,subsubPhenData,peptidesInfo)
    }
    return()
  }
  
  if (!is.na(peptidesInfo)){
    nPepsOrth <- sum(grepl(';',peptidesInfo$Accession))
    nProtsOrth <- sum(grepl(';',unique(peptidesInfo$Accession)))
    
    nPeps <- length(peptidesInfo$Accession) - nPepsOrth
    nProts <- length(unique(peptidesInfo$Accession)) - nProtsOrth
    
    write(paste(c('nPeps','nPepsOrth','nProts','nProtsOrth'),collapse = '\t'),logFile,append = F)
    write(paste(c(nPeps,nPepsOrth,nProts,nProtsOrth),collapse = '\t'),logFile,append = T)
  }
  
  
  nSampls <- ncol(expresionMatrix)
  write(paste0('nSampls: ',nSampls),logFile,append = T)
  write(paste(paste(names(table(subPhenData$class)),table(subPhenData$class),sep=': '),collapse = '\t'),logFile,append = T)

  conditionCol <- "class"
  classes <- levels(as.factor(subPhenData[,conditionCol]))
  colors <- paletteer_dynamic("cartography::multi.pal", length(classes))[1:length(classes)]; names(colors) <- classes
  allColors <- colors[match(subPhenData[,conditionCol],names(colors))]
  colorsMatch <- as.vector(allColors); 
  names(colorsMatch) <- subPhenData[,conditionCol]
  #colorsMatch <- list(colorsMatch)
  
  totalSignal <- apply(expresionMatrix, 2, sum,na.rm=T)
  #totalSignal <- totalSignal / max(totalSignal) * 100  
  totalSignal <- cbind(subPhenData,summedSignal=totalSignal)
  totalSignal$class <- factor(totalSignal$class, levels = unique(totalSignal$class))
  totalSignal$id <- factor(totalSignal$id, levels = unique(totalSignal$id))
  myPlot <- ggplot(totalSignal, aes(x=id, y=summedSignal, fill=class)) +
    geom_bar(stat="identity")+theme_classic(base_size = 10)+scale_fill_manual(values = colorsMatch) + 
    geom_hline(yintercept = 50)
  ggsave(file.path(outdir,paste0(tagname,'_SignalPlot',extension)),myPlot,
         width = length(subPhenData$id) * 100, device=gsub('\\.','',extension),units = "px",dpi=300)
  
  toBoxPlot <- melt(log2(expresionMatrix+1))
  toBoxPlot$variable <- factor(as.character(toBoxPlot$variable),levels=unique(toBoxPlot$variable))
  if (tagname == "Transcriptome"){
    toBoxPlot$class <- subPhenData$class[match(toBoxPlot$variable, subPhenData$id)] # 
  }else{
    toBoxPlot$class <- subPhenData$class[match(toBoxPlot$variable, as.character(as.numeric(subPhenData$id)))]
  }
  myBox <- ggplot(toBoxPlot, aes(x=variable, y=value, fill=class)) + geom_boxplot() +
    scale_fill_manual(values = unique(allColors),aesthetics = "fill") +
    theme_classic() + theme(legend.position="right",plot.title = element_text(size=11)) +
    ggtitle(gsub("_"," ",tagname)) + xlab("")
  ggsave(file.path(outdir,paste0(tagname,'_BoxPlot',extension)),myBox,device=gsub('\\.','',extension),units = "px",dpi=300)
  
  if (any(is.na(expresionMatrix))){
    
    totalElements <- ncol(expresionMatrix) * nrow(expresionMatrix)
    missingElements <- sum(is.na(expresionMatrix))
    totalPercentMissing <- missingElements / totalElements * 100
    
    SamplesMissingInGenes <- apply(is.na(expresionMatrix),1,sum)
    GenesMissingInSamples <- round(apply(is.na(expresionMatrix),2,sum) / nrow(expresionMatrix) * 100,2)
    
    write(paste(c("totalElements: ",totalElements,"% missing: ",totalPercentMissing) ,collapse = '\t'),logFile,append = T)
    
    write('Missing Peps Values % Stats Per Sample:',logFile,append = T)
    write(paste(subPhenData$id,collapse = '\t'),logFile,append = T)
    write(paste(subPhenData$class,collapse = '\t'),logFile,append = T)
    write(paste(GenesMissingInSamples,collapse = '\t'),logFile,append = T)
    
    write('% of peptides with X samples missing:', logFile, append = T)
    write(paste(round(table(SamplesMissingInGenes) / nrow(expresionMatrix) * 100,digits = 2),collapse = '\t'),logFile,append = T)
    write(paste(names(table(SamplesMissingInGenes)),collapse = '\t'),logFile,append = T)
    
    column_ha = HeatmapAnnotation(GenesMissing = anno_barplot(GenesMissingInSamples,add_numbers = T,
                                                              numbers_rot=0,numbers_offset = unit(1, "mm"),
                                                              numbers_gp=gpar(fontsize = 6)),
                                  treatment = subPhenData[,conditionCol],col = list(class=colorsMatch))
    row_ha = rowAnnotation(SamplesMissing = anno_barplot(SamplesMissingInGenes))
    toHMP <- t(scale(t(expresionMatrix)))
    col_fun <- colorRamp2(c(min(toHMP,na.rm = T),0,max(toHMP,na.rm = T)), c("blue","white","red"))
    mydevice(file.path(outdir,paste0(tagname,'_NAs_HeatMap',extension)), ncol(toHMP) * 120, 1500, pointsize=5, res = 300)
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
  mydevice(file.path(outdir,paste0(tagname,'_PCA2D',extension)), 1500, 1500, pointsize=12, res = 300)
  print(ggplot2::autoplot(pcares,data = data.frame(t(expresionMatrix),class=subPhenData[,conditionCol]),
                          colour="class", label=T, label.size=2,label.repel=T,main = gsub("_"," ",tagname)) + 
          scale_color_manual(values = unique(as.vector(colors)))+ theme_classic() + theme(legend.title=element_blank())+
          guides(color = guide_legend(override.aes = aes(label = ""))))
  dev.off()
  
  # pcaVarGroup <- round((pcares$sdev)^2 / sum(pcares$sdev^2) *100, 2)
  # mydevice(file.path(outdir,paste0(tagname,'_PCA3D',extension)), 1200, 900, pointsize=20, res = 300)
  # pca3d <- scatterplot3d(pcares$x[,1], pcares$x[,3], pcares$x[,2],color = as.vector(allColors),
  #                        main = tagname, xlab = paste0("PC 1 (", pcaVarGroup[1], " %)"),
  #                        zlab = paste0("PC 2 (", pcaVarGroup[2], " %)"), ylab = paste0("PC 3 (", pcaVarGroup[3], " %)"),
  #                        grid=T, box = T, pch = 20, cex.symbols = 2.5, angle = 40,type = "h")
  # zz.coords <- pca3d$xyz.convert(pcares$x[,1], pcares$x[,3], pcares$x[,2])
  # legend("topright", pch=20, legend = classes, col = unique(as.vector(allColors)),inset = 0,y.intersp =0.8)
  # text(zz.coords$x, zz.coords$y, labels = row.names(pcares$x),cex = 0.5, pos = 4,col = as.vector(allColors))
  # dev.off()

  toHMP <- t(scale(t(expresionMatrix)))
  toHMP[is.na(toHMP)] <- 0
  col_fun <- colorRamp2(c(min(toHMP,na.rm = T),0,max(toHMP,na.rm = T)), c("blue","white","red"))
  column_ha <- HeatmapAnnotation(treatment = subPhenData[,conditionCol],col = list(class=colorsMatch))
  mydevice(file.path(outdir,paste0(tagname,'_HeatMap',extension)), 1200, 900, pointsize=20, res = 300)
  print(
    Heatmap(toHMP,cluster_columns = F, cluster_rows = T,col=col_fun,
                top_annotation = column_ha,show_row_names = F,show_row_dend = F,na_col = "grey",
        column_title_gp=gpar(fontsize = 6),
        column_names_gp = gpar(fontsize = 6),
        row_names_gp = gpar(fontsize = 6))
    )
  dev.off()
  
  corMat <- cor(expresionMatrix,method = "pearson")
  
  col_fun <- colorRamp2(c(min(toHMP,na.rm = T),0,max(toHMP,na.rm = T)), c("blue","white","red"))
  mydevice(file.path(outdir,paste0(tagname,'_CorrPlot',extension)), ncol(toHMP) * 100, ncol(toHMP) * 100, pointsize=8, res = 300)
  corrplot(corMat, method = 'circle', type = 'upper', order = 'original', tl.cex=1, title=tagname,
           tl.col = allColors,tl.srt = 45,cl.ratio = 0.2, number.cex = 0.5,
           mar = c(5, 1, 2, 1),tl.pos = 'lt',addCoef.col = "grey",
           col=rev(COL2('RdBu', 200)))
  print(legend("bottom", legend=classes,fill=colors, horiz=T, cex=0.5,y = 0, bty = "n"))
  dev.off()
}


################################################################################
############################ Vizs of RAW #######################################
################################################################################

inFile <- "/home/eidrian/Desktop/realServices/multiregulatomics/FinalData/RNAseqCorrected/counts_trans.csv"
inFile <- read.csv(inFile)
rownames(inFile) <- inFile$X; inFile$X <- NULL
subPhenData <- as.data.frame(cbind(id=colnames(expression),class=gsub('\\..*','',colnames(expression))))
subPhenData <- subPhenData[order(gsub("YPD","c.YPD",subPhenData$class)),]
inFile <- inFile[subPhenData$id]
subPhenData$id <- gsub("brep","",subPhenData$id)
colnames(inFile) <- subPhenData$id; 
rownames(subPhenData) <- subPhenData$id
tagname <- "Transcriptome"
someVisualization(inFile,tagname,subPhenData,peptidesInfo=NA)

length(row.names(inFile))


#### UV RBPome
inFile <- "FinalData/Raw/correctedRAW/UVXmRBPomeRAW.csv"
tagname <- "UV-mRBPome"
someVisualization(inFile,tagname)

#### FAX RBPome
inFile <- "FinalData/Raw/correctedRAW/FAXmRBPomeRAW.csv"
tagname <- "FAX-RBPome"
someVisualization(inFile,tagname)

#### NOX RBPome
inFile <- "FinalData/Raw/correctedRAW/NOXmRBPomeRAW.csv"
tagname <- "NOX-RBPome"
someVisualization(inFile,tagname)

#### Proteome
inFile <- "FinalData/Raw/correctedRAW/newProteome.csv"
tagname <- "Proteome"
someVisualization(inFile,tagname)




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
rawFile <- "data/raw/MQtoSQ.csv"
doneFile <- "data/old/P445_AG_PolyA_40_all_3_cond_20220611_Final/P445_AG_PolyA_40_20220617_gmin_UV_nonorm_p2_NOX_Final/P445_AG_PolyA_40_20220617_gmin_UV_nonorm_p2_NOX_PROTEIN.tsv"
compareRawNDone(doneFile,rawFile)

doneFile <- "data/new/PolyA_RIC_nonorm_20220928/P445_AG_PolyA_40_20220617_gmin_UV_nonorm_p2_NOX_Final/P445_AG_PolyA_40_20220617_gmin_UV_nonorm_p2_NOX_PROTEIN.tsv"
compareRawNDone(doneFile,rawFile)

### UVX RBPome
rawFile <- "data/raw/MQtoSQ_UV.csv"
doneFile <- "data/old/P445_AG_PolyA_40_all_3_cond_20220611_Final/P445_AG_PolyA_40_20220617_gmin_UV_p2_no511_Final/P445_AG_PolyA_40_20220617_gmin_UV_p2_no511_PROTEIN.tsv"
compareRawNDone(doneFile,rawFile)

doneFile <- "data/new/PolyA_RIC_nonorm_20220928/P445_AG_PolyA_40_20220617_gmin_UV_p2_nonorm_no511_Final/P445_AG_PolyA_40_20220617_gmin_UV_p2_nonorm_no511_Final_PROTEIN.tsv"
compareRawNDone(doneFile,rawFile)
### FAX RBPome
rawFile <- "data/raw/MQtoSQ_FAX.csv"
doneFile <- "data/old/P445_AG_PolyA_40_all_3_cond_20220611_Final/P445_AG_PolyA_40_20220617_gmin_FAX_p2_Final/P445_AG_PolyA_40_20220617_gmin_FAX_p2_PROTEIN.tsv"
compareRawNDone(doneFile,rawFile)

doneFile <- "data/new/PolyA_RIC_nonorm_20220928/P445_AG_PolyA_40_20220617_gmin_FAX_p2_nonorm_no12_Final/P445_AG_PolyA_40_20220617_gmin_FAX_p2_nonorm_no12_Final_PROTEIN.tsv"
compareRawNDone(doneFile,rawFile)

### Proteome
rawFile <- "data/raw/SNtoSQ-Prot40_SafeQuant_Input.csv"

### NOX Proteome
crosslink <- "NOX"

doneFile <- "data/new/P445_Prot_DIA_iBAQ_nocomb_20220826_all_NOX_Final/P445_Prot_DIA_iBAQ_noOutliers_nocomb_20220826_NOX_Final/P445_Prot_DIA_iBAQ_noOutliers_nocomb_20220826_NOX_Final_PROTEIN.tsv"
compareRawNDone(doneFile,rawFile,crosslink)

doneFile <- "data/new/P445_Prot_DIA_iBAQ_nocomb_20220826_all_NOX_Final/P445_Prot_DIA_iBAQ_nocomb_20220826_NOX_Final2/P445_Prot_DIA_iBAQ_nocomb_20220826_NOX_Final2_PROTEIN.tsv"
compareRawNDone(doneFile,rawFile,crosslink)

doneFile <- "data/old/P445_Prot_DIA_iBAQ_noOutliers_nocomb_20220826_Final/P445_Prot_DIA_iBAQ_noOutliers_nocomb_20220826_NOX_Final/P445_Prot_DIA_iBAQ_noOutliers_nocomb_20220826_NOX_Final_PROTEIN.tsv"
compareRawNDone(doneFile,rawFile,crosslink)

### UVX Proteome
crosslink <- "UV"
doneFile <- "data/old/P445_Prot_DIA_iBAQ_noOutliers_nocomb_20220826_Final/P445_Prot_DIA_iBAQ_noOutliers_nocomb_20220826_UV_Final/P445_Prot_DIA_iBAQ_noOutliers_nocomb_20220826_UV_Final_PROTEIN.tsv"
compareRawNDone(doneFile,rawFile,crosslink)

### FAX Proteome
crosslink <- "FAX"
doneFile <- "data/old/P445_Prot_DIA_iBAQ_noOutliers_nocomb_20220826_Final/P445_Prot_DIA_iBAQ_noOutliers_nocomb_20220826_FAX_Final/P445_Prot_DIA_iBAQ_noOutliers_nocomb_20220826_FAX_Final_PROTEIN.tsv"
compareRawNDone(doneFile,rawFile,crosslink)


################################################################################
############################ Comparing NOX vs FAX and UVX ######################
################################################################################

inFile <- "data/raw/MQtoSQ.csv"
myDF <- read.csv(inFile)

normDCol <- grep("Norm",colnames(myDF))
rawDCol <- grep("Raw",colnames(myDF))
expresionMatrix <- myDF[normDCol:(rawDCol-1)]
expresionMatrix[2,] <- gsub("_.*","",expresionMatrix[2,])
samplesInNOXraw <- as.character(expresionMatrix[2,])
subPhenDataNOX <- as.data.frame(t(expresionMatrix[1:2,])); colnames(subPhenData) <- c("class","id"); rownames(subPhenData) <- subPhenData[,2]


UVXresultsFile <- "FinalData/SafeQuantResults/RBPomeUVX_np2_norm_gmin_no10-7-16-3/RBPomeUVX_np2_norm_gmin_no10-7-16-3_PROTEIN.tsv"
UVXresultsDF <- read.delim(UVXresultsFile)
UVXresultsSamples <- colnames(UVXresultsDF)[grep("^X",colnames(UVXresultsDF))]
UVXresultsSamples <- gsub("^X","",gsub("_.*","",UVXresultsSamples))


FAXresultsFile <- "FinalData/SafeQuantResults/RBPomeFAX_np2_norm_gmin_no12/RBPomeFAX_np2_norm_gmin_no12_PROTEIN.tsv"
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
rawRBPomeFAXfile <- "FinalData/Raw/MQtoSQ_FAX.csv"
rawRBPomeFAXSamplesInfo <- read.delim(rawRBPomeFAXfile,nrows = 3, sep = ",")

rawRBPomeUVXfile <- "FinalData/Raw/MQtoSQ_UV.csv"
rawRBPomeUVXSamplesInfo <- read.delim(rawRBPomeUVXfile,nrows = 3, sep = ",")

rawRBPomeNOXfile <- "FinalData/Raw/MQtoSQ.csv"
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
proteomeFile <- "FinalData/Raw/correctedRAW/newProteome.csv"
proteomeDF <- read.delim(proteomeFile, sep = ",")
length(proteomeDF[,grep("Acc",proteomeDF[2,])]) - 2
length(unique(proteomeDF[,grep("Acc",proteomeDF[2,])])) - 2

rawRBPomeFAXfile <- "FinalData/Raw/correctedRAW/FAXmRBPomeRAW.csv"
rawRBPomeFAXDF <- read.delim(rawRBPomeFAXfile, sep = ",")
length(rawRBPomeFAXDF[,grep("Acc",rawRBPomeFAXDF[2,])]) - 2
length(unique(rawRBPomeFAXDF[,grep("Acc",rawRBPomeFAXDF[2,])])) - 2

rawRBPomeUVXfile <- "FinalData/Raw/correctedRAW/UVXmRBPomeRAW.csv"
rawRBPomeUVXDF <- read.delim(rawRBPomeUVXfile, sep = ",")
length(rawRBPomeUVXDF[,grep("Acc",rawRBPomeUVXDF[2,])]) - 2
length(unique(rawRBPomeUVXDF[,grep("Acc",rawRBPomeUVXDF[2,])])) - 2

rawRBPomeNOXfile <- "FinalData/Raw/correctedRAW/NOXmRBPomeRAW.csv"
rawRBPomeNOXDF <- read.delim(rawRBPomeNOXfile, sep = ",")
length(rawRBPomeNOXDF[,grep("Acc",rawRBPomeNOXDF[2,])]) - 2
length(unique(rawRBPomeNOXDF[,grep("Acc",rawRBPomeNOXDF[2,])])) - 2

