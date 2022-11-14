# Load libraries
library(ggVennDiagram)
library(Biobase)
library(limma)
library(ggplot2)
library(circlize)
library(paletteer)
library(PCAtools)
library("scatterplot3d")
library(RColorBrewer)
library(ggfortify)
library(openxlsx)
library(EnhancedVolcano)
library(scatterplot3d)
library(stringr)
library(gdata)
library(corrplot)
library(ComplexHeatmap)
library(dplyr)
library(fgsea)
library(httr)
library(jsonlite)
library(RCurl)
library(plyr)
source("scripts/GC4libR.R")
library(reshape2)

summaryGC4results <- function(GC4res,outfile){
  allQC <- c()
  allRes <- c()
  for (res in names(GC4res$quality_controls[[1]])){
    resX <- GC4res$quality_controls[[1]][[res]]
    annot <- gsub(".*-","",res)
    if (res == "notInDB"){
      allQC <- rbind(allQC,c("notInDB","","","","",
                             length(resX$invalidInput),
                             paste(resX$invalidInput,collapse = " "),
                             paste(resX$invalidUniverse,collapse = " "),
                             paste(resX$notMapped,collapse = " ")))
    }else{
      allQC <- rbind(allQC,c(annot,length(resX$annotated),length(resX$noAnnotated),
                             paste(resX$noAnnotated,collapse = " "),resX$universe,
                             "",
                             "",
                             "",
                             ""))
    }
  }
  allQC <- as.data.frame(allQC)
  colnames(allQC) <- c("annotation","annotated","noannotated","noannotatedlist","universe","NinvalidInput","invalidInput","invalidUniverse","notMapped")
  for (res in 1:length(GC4res$stats_tables)){
    annot <- gsub(".*-","",names(GC4res$stats_tables)[res])
    resX <- GC4res$stats_tables[[res]]
    resX$annotation <- annot
    allRes <- rbind(allRes,resX)
  }
  return(list(qc=allQC,enr=allRes))
}

completeNdedup <- function(DF){
  if (class(DF) == "data.frame"){
    return(DF[complete.cases(DF),])
  }else{
    return(unique(na.omit(DF)))
  }
}

saveTablesTsvExc <- function(DF,outdir,completeNdedup=T,excel=T,bycompleteFC=T,rownames=F){
  outname <- deparse(substitute(DF))
  if (class(DF) == "data.frame"){
    if (nrow(DF) == 0){
      print("Empty DF")
    }else{
      outFile <- paste0(file.path(outdir,outname),".tsv")
      if (completeNdedup){
        #DF <- DF[complete.cases(DF),]
        if(bycompleteFC){
          fcORnetCol <- grep("log2|netChanges",colnames(DF))
          DF <- DF[!is.na(DF[,fcORnetCol]),]
          DF <- DF[!is.na(DF[order(DF[,fcORnetCol],decreasing = T),]),]
        }
        DF <- DF[!duplicated(DF),]
      }
      write.table(DF,outFile,quote = F,sep = '\t',col.names = T,row.names = rownames)
      if (excel){
        write.xlsx(DF,gsub(".tsv",".xlsx",outFile), asTable = FALSE, overwrite = TRUE,rowNames = rownames)
      }
    }
  } else{
    if (length(DF) == 0){
      print("Empty Vector")
    }else{
      outFile <- paste0(file.path(outdir,outname),".txt")
      write(unique(na.omit(DF)),outFile,sep = "\n")
    }
  }
}

transformFCto1scale <- function(log2fcarray){
  log2fcarray <- log2fcarray + abs(min(log2fcarray))
  log2fcarray <- (log2fcarray - min(log2fcarray)) / (max(log2fcarray) - min(log2fcarray))
  return(log2fcarray)
}

createHMtop10ORA <- function(enrichFile){
  plotEnrFile <- gsub(".tsv",".png",gsub("enrichmentRes","visualizations/enrichment",enrichFile))
  GC4res <- read.delim(enrichFile)
  GC4res <- GC4res[order(GC4res$pval_adj,decreasing = F),]
  GC4res$pval_adj <- -log10(GC4res$pval_adj)
  res2HM <- c()
  for (annot in unique(GC4res$annotation)){
    subGC4res <- GC4res[GC4res$annotation == annot,][1:10,]
    res2HM <- rbind(res2HM,subGC4res)
  }
  #rownames(res2HM) <- paste0(res2HM$description,rownames(res2HM))
  toHMAP <- res2HM[,"pval_adj",drop=F]
  col_funGenF <- colorRamp2(c(min(res2HM$genes_found),max(res2HM$genes_found)), c("cadetblue1","navy"))
  col_funRelE <- colorRamp2(c(min(res2HM$relative_enrichment),max(res2HM$relative_enrichment)), c("coral","brown3"))
  col_funFDR <- colorRamp2(c(min(toHMAP$pval_adj),max(toHMAP$pval_adj)), c("gold","maroon3"))
  ha = rowAnnotation(
    "Genes Sig." = anno_simple(res2HM$genes_found, col = col_funGenF),
    "Rel. Enr." = anno_simple(res2HM$relative_enrichment, col = col_funRelE),
    labels = anno_text(res2HM$description, which = "row"),
    show_annotation_name=c("Genes Sig." = F,
                           "Rel. Enr." = F),
    width = max(grobWidth(textGrob(res2HM$description))))
  myHMplot <- Heatmap(toHMAP,cluster_columns = F, cluster_rows = F,col=col_funFDR, show_row_dend = F,
                      right_annotation = ha, show_column_names = F, show_row_names = F,
                      row_split = gsub("_"," ",res2HM$annotation), row_gap = unit(5, "mm"),
                      name="-Log10(FDR)", show_heatmap_legend = T,
                      width = ncol(toHMAP)*unit(5, "mm"), height = nrow(toHMAP)*unit(5, "mm"))
  lgd_FDR = Legend(title = "-Log10(FDR)", col_fun = col_funFDR)
  lgd_GenF = Legend(title = "Genes Sig.", col_fun = col_funGenF)
  lgd_pRelE = Legend(title = "Rel. Enr.", col_fun = col_funRelE)
  png(plotEnrFile, 1200, 900, units = "px", pointsize=12, res = "print")
  draw(myHMplot, heatmap_legend_side = "left", merge_legend = TRUE, 
       annotation_legend_list = list(lgd_GenF, lgd_pRelE),annotation_legend_side = "left")
  dev.off()
}

createHMtop10GSEA <- function(enrichFile){
  plotEnrFile <- gsub(".tsv",".png",gsub("enrichmentRes","visualizations/enrichment",enrichFile))
  GC4res <- read.delim(enrichFile)
  GC4res <- GC4res[order(GC4res$padj,decreasing = F),]
  GC4res$padj <- -log10(GC4res$padj)
  res2HM <- c()
  for (db in unique(GC4res$db)){
    subGC4res <- GC4res[GC4res$db == db,][1:10,]
    res2HM <- rbind(res2HM,subGC4res)
  }
  #rownames(res2HM) <- paste0(res2HM$description,rownames(res2HM))
  toHMAP <- res2HM[,"padj",drop=F]
  col_funGenF <- colorRamp2(c(min(res2HM$size),max(res2HM$size)), c("cadetblue1","navy"))
  col_funRelE <- colorRamp2(c(min(res2HM$NES),max(res2HM$NES)), c("coral","brown3"))
  col_funFDR <- colorRamp2(c(min(toHMAP$padj),max(toHMAP$padj)), c("gold","maroon3"))
  
  ha = rowAnnotation(
    "Leading Edge Genes" = anno_simple(res2HM$size, col = col_funGenF),
    "NES" = anno_simple(res2HM$NES, col = col_funRelE),
    labels = anno_text(res2HM$pathway, which = "row"),
    show_annotation_name=c("Leading Edge Genes" = F,
                           "NES" = F),
    width = max(grobWidth(textGrob(res2HM$pathway))))
  myHMplot <- Heatmap(toHMAP,cluster_columns = F, cluster_rows = F,col=col_funFDR, show_row_dend = F,
                      right_annotation = ha, show_column_names = F, show_row_names = F,
                      row_split = gsub("_"," ",res2HM$db), row_gap = unit(5, "mm"),
                      name="-Log10(FDR)", show_heatmap_legend = T,
                      width = ncol(toHMAP)*unit(5, "mm"), height = nrow(toHMAP)*unit(5, "mm"))
  lgd_FDR = Legend(title = "-Log10(FDR)", col_fun = col_funFDR)
  lgd_GenF = Legend(title = "Leading Edge Genes", col_fun = col_funGenF)
  lgd_pRelE = Legend(title = "NES", col_fun = col_funRelE)
  png(plotEnrFile, 1200, 900, units = "px", pointsize=12, res = "print")
  draw(myHMplot, heatmap_legend_side = "left", merge_legend = TRUE, 
       annotation_legend_list = list(lgd_GenF, lgd_pRelE),annotation_legend_side = "left")
  dev.off()
}

basicAnalysis <- function(dataFile,ourPhenoData,outdir,tagname,mytreatment,mappingDF,control="no",conditionCol="Treatment",pCutoff = 0.05,FCcutoff = 1){
  if (class(dataFile) == "data.frame"){
    myDF <- dataFile
  }else{
    myDF <- read.delim(dataFile)
  }
  
  dir.create(outdir,F)
  
  expresionMatrix <- myDF[,grep("X[0-9]+",colnames(myDF))]

  runnumbers <- gsub("X0+","",unlist(lapply(strsplit(colnames(expresionMatrix),"_"), function(x) return(x[[1]]))))
  expresionMatrix <- expresionMatrix[,colnames(expresionMatrix)[match(ourPhenoData$Run.number,runnumbers)]]
  runnumbers <- gsub("X0+","",unlist(lapply(strsplit(colnames(expresionMatrix),"_"), function(x) return(x[[1]]))))
  subPhenData <- as.data.frame(cbind(cols=colnames(expresionMatrix),"Run.number"=runnumbers))
  subPhenData <- merge(subPhenData,ourPhenoData,by="Run.number")
  subPhenData <- subPhenData[order(subPhenData[,conditionCol]),]
  row.names(subPhenData) <- subPhenData$cols
  
  expresionMatrix <- expresionMatrix[subPhenData$cols]
  
  myDF$ac <- gsub("-.*","",myDF$ac)
  expresionMatrix <- cbind(acc=myDF$ac,expresionMatrix)
  expresionMatrix <- as.data.frame(summarise_all(group_by(expresionMatrix,acc),mean))
  
  expresionMatrixMapped <- merge(mappingDF,expresionMatrix,by.x = "UniprotACC",by.y = "acc")
  saveTablesTsvExc(expresionMatrixMapped,outdir,completeNdedup=F,excel=T,bycompleteFC=F,rownames = F)
  saveTablesTsvExc(subPhenData,outdir,completeNdedup=F,excel=T,bycompleteFC=F,rownames = T)
  
  row.names(expresionMatrix) <- expresionMatrix$acc
  expresionMatrix$acc <- NULL
  ###################
  ## Visualiztions ##
  ###################
  
  ## To detect sample outliers (Outlier= 1*5 Inter Quartile Range)
  classes <- levels(as.factor(subPhenData[,conditionCol]))
  colors <- paletteer_dynamic("cartography::multi.pal", length(classes))[1:length(classes)]; names(colors) <- classes
  allColors <- colors[match(subPhenData[,conditionCol],names(colors))]
  
  png(file.path(outdir,paste0(tagname,'_BoxPlot.png')), 1000, 1000, pointsize = 20,res = "print")
  par(mar=c(17,2,1,1))
  print(boxplot(log2(expresionMatrix),main=paste("Log2 Expr",'mrbpome',crosslink),
                at=c(1:ncol(expresionMatrix)),las=2,horizontal=F,
                notch=F, bty='L',col=allColors,range=100,height=1500, units='px'))
  print(legend("bottom", legend=classes,fill=colors, horiz=TRUE, cex=1,y = -10))
  dev.off()
  
  pcares <- prcomp(t(expresionMatrix))
  png(file.path(outdir,paste0(tagname,'_PCA2D.png')), 600, 500, pointsize=50, res = 100)
  print(ggplot2::autoplot(pcares,data = data.frame(t(expresionMatrix),class=subPhenData[,conditionCol]),
                          colour="class", label=T, label.size=2,label.repel=T,main = tagname) + 
          scale_color_manual(values = unique(as.vector(colors)))+ theme_classic())
  dev.off()
  
  pcaVarGroup <- round((pcares$sdev)^2 / sum(pcares$sdev^2) *100, 2)
  png(file.path(outdir,paste0(tagname,'_PCA3D.png')), 1200, 900, pointsize=12,res = "print")
  pca3d <- scatterplot3d(pcares$x[,1], pcares$x[,3], pcares$x[,2],color = as.vector(allColors),
                         main = tagname, xlab = paste0("PC 1 (", pcaVarGroup[1], " %)"),
                         zlab = paste0("PC 2 (", pcaVarGroup[2], " %)"), ylab = paste0("PC 3 (", pcaVarGroup[3], " %)"),
                         grid=T, box = T, pch = 20, cex.symbols = 2.5, angle = 40,type = "h")
  zz.coords <- pca3d$xyz.convert(pcares$x[,1], pcares$x[,3], pcares$x[,2])
  legend("topright", pch=20, legend = classes, col = unique(as.vector(allColors)),inset = 0,y.intersp =0.8)
  text(zz.coords$x, zz.coords$y, labels = row.names(pcares$x),cex = 0.5, pos = 4,col = as.vector(allColors))
  dev.off()
  
  corMat <- cor(expresionMatrix,method = "pearson")
  png(file.path(outdir,paste0(tagname,'_CorrPlot.png')), 1200, 900, pointsize=12,res = "print")
  corrplot(corMat, method = 'circle', order = 'original',tl.cex=0.5,title=tagname,
           mar = c(0, 1, 2, 1),tl.pos = 'lt',addCoef.col = "grey",number.cex = 0.5,
           col=rev(COL2('RdBu', 200)))
  dev.off()
  
  toHMP <- t(scale(t(expresionMatrix)))
  toHMP[is.na(toHMP)] <- 0
  col_fun <- colorRamp2(c(min(toHMP,na.rm = T),0,max(toHMP,na.rm = T)), c("blue","white","red"))
  column_ha <- HeatmapAnnotation(treatment = subPhenData[,conditionCol])
  png(file.path(outdir,paste0(tagname,'_HeatMap.png')), 1200, 900, pointsize=12,res = "print")
  print(Heatmap(toHMP,cluster_columns = F, cluster_rows = T,col=col_fun,top_annotation = column_ha,show_row_names = F,show_row_dend = F))
  dev.off()
  
  ### PREPARE DATA 
  subPhenData[,conditionCol] <- as.factor(subPhenData[,conditionCol])
  phenoD <- AnnotatedDataFrame(data=subPhenData[,conditionCol,drop=FALSE])
  featD <- AnnotatedDataFrame(data=data.frame(row.names=row.names(expresionMatrix)))
  eset <- ExpressionSet(assayData = as.matrix(expresionMatrix),phenoData=phenoD,featureData = featD)
  exprs(eset) <- log2(exprs(eset))
  
  diffExpr <- data.frame(row.names=featureNames(eset))
  #selCol <- unlist(pData(eset)[,conditionCol] %in% c(control,mytreatment))
  fit <- eBayes(lmFit(eset, model.matrix(~factor(pData(eset)[,conditionCol]), data=pData(eset))))
  diffExpr$pVal <- fit$p.value[,2]
  diffExpr$FDR <- p.adjust(diffExpr$pVal,method="BH")
  
  eCtrl <- subset(exprs(eset),select=which(pData(eset)[,conditionCol] == control))
  eCase <- subset(exprs(eset),select=which(pData(eset)[,conditionCol] == mytreatment))
  ratios <- apply(eCase,1,"median", na.rm=T) - apply(eCtrl,1,"median", na.rm=T)
  diffExpr$log2FC <- ratios
  
  volcano <- EnhancedVolcano(diffExpr, x='log2FC', y='FDR', lab=row.names(diffExpr),
                             title = tagname,subtitle = '',
                             pCutoff = pCutoff,FCcutoff = FCcutoff,
                             legendLabels = c('NS', expression(Log[2]~FC),
                                              'FDR', expression(p~and~log[2]~FC)))
  ggsave(file.path(outdir,paste0(tagname,'_VulcanoPlot.png')), plot = volcano, 
         device = 'png', scale = 1, dpi = "print",width = 25,height = 25,units = 'cm')
  
  colnames(diffExpr) <- paste0(colnames(diffExpr),".",mytreatment)
  diffExprMapped <- merge(mappingDF,diffExpr,by.x = "UniprotACC",by.y = "row.names",all.y = T)
  diffExprMapped <- diffExprMapped[!duplicated(diffExprMapped),]
  diffExprMapped <- diffExprMapped[order(diffExprMapped[,paste0("FDR",".",mytreatment)]),]
  saveTablesTsvExc(diffExprMapped,outdir,completeNdedup=F,excel=T,bycompleteFC=F,rownames = F)
  
  fullMerged <- merge(expresionMatrixMapped,diffExprMapped)
  fullMerged <- fullMerged[!duplicated(fullMerged),]
  saveTablesTsvExc(fullMerged,outdir,completeNdedup=F,excel=T,bycompleteFC=F,rownames = F)
}

getPlotsRawData <- function(inFile,tagname){
  myDF <- read.delim(inFile,sep = ",")
  outdir <- "/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/adrian/rawViz"
  normDCol <- grep("Norm",colnames(myDF))
  rawDCol <- grep("Raw",colnames(myDF))
  normData <- myDF[normDCol:(rawDCol-1)]
  someVisualization(normData,outdir,paste0("norm",tagname))
  #rawData <- myDF[rawDCol:(rawDCol+ncol(normData)-1)]
  #someVisualization(rawData,outdir,paste0("raw",tagname))
}

expresionMatrix <- normData

someVisualization <- function(expresionMatrix,outdir,tagname,subPhenData=NA){
  logfile <- file.path(outdir,paste0(tagname,".log"))
  
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
  
  toBoxPlot <- melt(log2(expresionMatrix))
  toBoxPlot$Var2 <- factor(as.character(toBoxPlot$Var2),levels=unique(toBoxPlot$Var2))
  toBoxPlot$class <- subPhenData$class[match(toBoxPlot$Var2, as.character(as.numeric(subPhenData$id)))]
  myBox <- ggplot(toBoxPlot, aes(x=Var2, y=value, fill=class)) + geom_boxplot() +
        scale_fill_manual(values = unique(allColors),aesthetics = "fill") +
        theme_classic() + theme(legend.position="right",plot.title = element_text(size=11)) +
        ggtitle(gsub("_"," ",tagname)) + xlab("")
  ggsave(file.path(outdir,paste0(tagname,'_BoxPlot.tiff')),myBox,device="tiff",units = "px",dpi=300)
  
  if (any(is.na(expresionMatrix))){
    SamplesMissingInGenes <- apply(is.na(expresionMatrix),1,sum)
    GenesMissingInSamples <- round(apply(is.na(expresionMatrix),2,sum) / nrow(expresionMatrix) * 100,2)
    
    column_ha = HeatmapAnnotation(GenesMissing = anno_barplot(GenesMissingInSamples,add_numbers = T,
                                                              numbers_rot=0,numbers_offset = unit(1, "mm"),
                                                              numbers_gp=gpar(fontsize = 6)),
                                  treatment = subPhenData[,conditionCol],col = colorsMatch)
    row_ha = rowAnnotation(SamplesMissing = anno_barplot(SamplesMissingInGenes))
    toHMP <- t(scale(t(expresionMatrix)))
    col_fun <- colorRamp2(c(min(toHMP,na.rm = T),0,max(toHMP,na.rm = T)), c("blue","white","red"))
    tiff(file.path(outdir,paste0(tagname,'_NAs_HeatMap.tiff')), 2000, 1200, pointsize=5, res = 300)
    print(
    Heatmap(toHMP, cluster_columns=F, cluster_rows=F, col=col_fun,
            top_annotation=column_ha, right_annotation=row_ha,
            show_row_names=F, show_row_dend=F, na_col="grey",name = paste0(tagname,'_NAs'),
            show_heatmap_legend=F,  column_title_gp=gpar(fontsize = 6),
            column_names_gp = gpar(fontsize = 6),
            row_names_gp = gpar(fontsize = 6))
    )
    dev.off()
  }

  expresionMatrix <- na.omit(expresionMatrix)

  write(paste0("Quality Control plots calculated with peptides: ",nrow(expresionMatrix)),logfile,append = T)
  
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
  tiff(file.path(outdir,paste0(tagname,'_CorrPlot.tiff')), 1500, 1500, pointsize=10, res = 300)
  corrplot(corMat, method = 'circle', type = 'upper', order = 'original', tl.cex=1, title=tagname,
            tl.col = allColors,tl.srt = 45,cl.ratio = 0.2, number.cex = 0.7,
            mar = c(5, 1, 2, 1),tl.pos = 'lt',addCoef.col = "grey",
            col=rev(COL2('RdBu', 200)))
  print(legend("bottom", legend=classes,fill=colors, horiz=T, cex=0.5,y = 0, bty = "n"))
  dev.off()
}


getPlotsOfConditions <- function(DF,outdir,tagname,pvalcutoff,FCcutoff){
  FCcol <- grep("FC",colnames(DF))
  PVcol <- grep("pval|FDR",colnames(DF))
  
  upReg <-  DF[,FCcol] > FCcutoff
  dwReg <-  DF[,FCcol] < -FCcutoff
  posReg <- DF[,FCcol] > 0
  negReg <- DF[,FCcol] < 0
  
  signt <- DF[,PVcol] < pvalcutoff
  
  upSigDEG <- DF[signt & upReg,]
  dwSigDEG <- DF[signt & dwReg,]
  posSigDEG <- DF[signt & posReg,]
  negSigDEG <- DF[signt & negReg,]
  
  outdir <- paste0(outdir,"/EA/",tagname,"/ORA/"); dir.create(outdir,recursive = T,showWarnings = F)
  obtainGC4res(upSigDEG$GeneName,outdir)
  obtainGC4res(dwSigDEG$GeneName,outdir)
  obtainGC4res(posSigDEG$GeneName,outdir)
  obtainGC4res(negSigDEG$GeneName,outdir)
  enrichFiles <- list.files(path = outdir,pattern = "*\\ORA.tsv",full.names = T,)
  for (enrichFile in enrichFiles){
    createHMtop10ORA(enrichFile)
  }
  
  ### GSEA
  toGSEA <- DF[,c("GeneName",FCcol)]
  outdir <- paste0(outdir,"/EA/",tagname,"/GSEA/"); dir.create(outdir,recursive = T,showWarnings = F)
  doGSEAs(toGSEA,gseaDB,outdir,tagname,minGeneset,maxGeneset)
  
  enrichFiles <- list.files(path = enrichResultsDir,pattern = "*\\GSEA.tsv",full.names = T,)
  for (enrichFile in enrichFiles){
    createHMtop10GSEA(enrichFile)
  }
 
}


doGSEAs <- function(toGSEA,gseaDB,outdir,tagname,minGeneset,maxGeneset) {
  FCcol <- grep("FC",colnames(toGSEA))
  toGSEA <- toGSEA[complete.cases(toGSEA),]
  toGSEA <- aggregate(toGSEA[,FCcol], list(toGSEA$GeneName), FUN=mean) 
  colnames(toGSEA) <- c("GeneName","log2FC")
  toGSEA <- toGSEA[order(toGSEA$log2FC,decreasing = T),]
  outfile <- paste0(outdir,"RankedInputGeneaNames.tsv")
  write.table(toGSEA,outfile,quote = F,sep = '\t',col.names = T,row.names = F)
  inputRank <- as.array(toGSEA$log2FC)
  names(inputRank) <- toGSEA$GeneName
  allGSEAres <- c()
  for (db in names(gseaDB)){
    dbName <- gsub("\\(|\\)","",db)
    gseaRes <- fgsea(pathways = gseaDB[[db]], stats = inputRank, minSize  = minGeneset, maxSize  = maxGeneset)
    gseaRes$db <- db
    allGSEAres <- rbind(allGSEAres,gseaRes)
  }
  allGSEAres$leadingEdge <- unlist(lapply(allGSEAres$leadingEdge, function(x) return(paste(x, collapse = ", "))))
  outfile <- paste0(outdir,tagname,"GSEA.tsv")
  allGSEAres <- allGSEAres[order(allGSEAres$padj,decreasing = F),]
  allGSEAres$pathway <- gsub("\\s*\\w*$","",allGSEAres$pathway)
  write.table(allGSEAres,outfile,quote = F,sep = '\t',col.names = T,row.names = F)
  return(allGSEAres)
}

obtainGC4res <- function(myinputlist,outdir){
  reportName <- gsub("\\$","",deparse(substitute(myinputlist)))
  inputGenes <- na.omit(unique(myinputlist))
  GC4res <- launchAnalysis(organism = "Saccharomyces cerevisiae",
                           inputType = "genes",
                           inputQuery = inputGenes,
                           annotationsDBs = c("KEGG","GO_BP","GO_CC","GO_MF"),
                           inputCoannotation = "no",
                           universeScope = "annotated",
                           enrichmentStat = "hypergeom",
                           ReportName=reportName)
  GC4resSumd <- summaryGC4results(GC4res)
  dir.create(outdir,showWarnings = F, recursive = T)
  outfile <- paste0(outdir,reportName,length(inputGenes),"genelist.txt")
  write(inputGenes,outfile,sep = "\n")
  outfile <- paste0(outdir,reportName,"ORA.tsv")
  write.table(GC4resSumd[["enr"]],outfile,quote = F,sep = '\t',col.names = T,row.names = F)
  outfile <- paste0(outdir,reportName,"ORA_QC.tsv")
  write.table(GC4resSumd[["qc"]],outfile,quote = F,sep = '\t',col.names = T,row.names = F)
}

