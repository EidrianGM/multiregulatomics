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

summaryGC4results <- function(resultsDFgc4,outfile){
  allQC <- data.frame(annotation=character(),
                      annotated=integer(),
                      noannotated=integer(),
                      noannotatedlist=character(),
                      universe=integer(),stringsAsFactors=FALSE)
  allRes <- c()
  for (res in names(resultsDFgc4$quality_controls[[1]])){
    annotation <- gsub(".*-","",names(resultsDFgc4$stats_tables)[res])
    resX <- resultsDFgc4$quality_controls[[1]][[res]]
    if (res == "notInDB"){
      allQC <- rbind.fill(allQC,data.frame(invalidInput=paste(resX$invalidInput,collapse = " "),
                                           invalidUniverse=paste(resX$invalidUniverse,collapse = " "),
                                           notMapped=paste(resX$notMapped,collapse = " "), 
                                           annotation="notInDB"))
    }else{
      allQC <- rbind(allQC,data.frame(annotation=annotation,
                                      annotated=length(resX$annotated),
                                      noannotated=length(resX$noAnnotated),
                                      noannotatedlist=paste(resX$noAnnotated,collapse = " "),
                                      universe=resX$universe))
    }
  }
  
  for (res in 1:length(resultsDFgc4$stats_tables)){
    annotation <- gsub(".*-","",names(resultsDFgc4$stats_tables)[res])
    resX <- resultsDFgc4$stats_tables[[res]]
    resX$annotation <- annotation
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

createHMtop10 <- function(enrResFile){
  plotEnrFile <- gsub(".tsv",".png",gsub("enrichmentRes","visualizations/enrichment",enrResFile))
  GC4res <- read.delim(enrResFile)
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

basicAnalysis <- function(dataFile,ourPhenoData,outdir,mytreatment,mappingDF,control="no"){
  dir.create(outdir,F)
  myDF <- read.delim(dataFile)
  
  expresionMatrix <- myDF[,grep("X[0-9]+",colnames(myDF))]

  runnumbers <- gsub("X0+","",unlist(lapply(strsplit(colnames(expresionMatrix),"_"), function(x) return(x[[1]]))))
  expresionMatrix <- expresionMatrix[,colnames(expresionMatrix)[match(ourPhenoData$Run.number,runnumbers)]]
  runnumbers <- gsub("X0+","",unlist(lapply(strsplit(colnames(expresionMatrix),"_"), function(x) return(x[[1]]))))
  subPhenData <- as.data.frame(cbind(cols=colnames(expresionMatrix),"Run.number"=runnumbers))
  subPhenData <- merge(subPhenData,ourPhenoData,by="Run.number")
  subPhenData <- subPhenData[order(subPhenData$Treatment),]
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
  classes <- levels(as.factor(subPhenData$Treatment))
  colors <- brewer.pal(n = length(classes), name = 'Set1')[1:length(classes)]; names(colors) <- classes
  allColors <- colors[match(subPhenData$Treatment,names(colors))]
  
  png(file.path(outdir,paste0(datatype,crosslink,'_BoxPlot.png')), 1000, 1000, pointsize = 20,res = "print")
  par(mar=c(17,2,1,1))
  print(boxplot(log2(expresionMatrix),main=paste("Log2 Expr",'mrbpome',crosslink),
                at=c(1:ncol(expresionMatrix)),las=2,horizontal=F,
                notch=F, bty='L',col=allColors,range=100,height=1500, units='px'))
  print(legend("bottom", legend=classes,fill=colors, horiz=TRUE, cex=1,y = -10))
  dev.off()
  
  pcares <- prcomp(t(expresionMatrix))
  png(file.path(outdir,paste0(datatype,crosslink,'_PCA2D.png')), 600, 500, pointsize=50, res = 100)
  print(ggplot2::autoplot(pcares,data = data.frame(t(expresionMatrix),class=subPhenData$Treatment),
                          colour="class", label=T, label.size=2,label.repel=T,main = paste0(datatype,crosslink)) + 
          scale_color_manual(values = unique(as.vector(colors)))+ theme_classic())
  dev.off()
  
  pcaVarGroup <- round((pcares$sdev)^2 / sum(pcares$sdev^2) *100, 2)
  png(file.path(outdir,paste0(datatype,crosslink,'_PCA3D.png')), 1200, 900, pointsize=12,res = 150)
  pca3d <- scatterplot3d(pcares$x[,1], pcares$x[,3], pcares$x[,2],color = as.vector(allColors),
                         main = paste0(datatype,crosslink), xlab = paste0("PC 1 (", pcaVarGroup[1], " %)"),
                         zlab = paste0("PC 2 (", pcaVarGroup[2], " %)"), ylab = paste0("PC 3 (", pcaVarGroup[3], " %)"),
                         grid=T, box = T, pch = 20, cex.symbols = 2.5, angle = 40,type = "h")
  zz.coords <- pca3d$xyz.convert(pcares$x[,1], pcares$x[,3], pcares$x[,2])
  legend("topright", pch=20, legend = classes, col = unique(as.vector(allColors)),inset = 0,y.intersp =0.8)
  text(zz.coords$x, zz.coords$y, labels = row.names(pcares$x),cex = 0.5, pos = 4,col = as.vector(allColors))
  dev.off()
  
  corMat <- cor(expresionMatrix,method = "pearson")
  png(file.path(outdir,paste0(datatype,crosslink,'_CorrPlot.png')), 1200, 900, pointsize=12,res = 150)
  corrplot(corMat, method = 'circle', order = 'original',tl.cex=0.5,title=paste0(datatype,crosslink),
           mar = c(0, 1, 2, 1),tl.pos = 'lt',addCoef.col = "grey",number.cex = 0.5,
           col=rev(COL2('RdBu', 200)))
  dev.off()
  
  toHMP <- t(scale(t(expresionMatrix)))
  toHMP[is.na(toHMP)] <- 0
  col_fun <- colorRamp2(c(min(toHMP,na.rm = T),0,max(toHMP,na.rm = T)), c("blue","white","red"))
  column_ha <- HeatmapAnnotation(treatment = subPhenData$Treatment)
  png(file.path(outdir,paste0(datatype,crosslink,'_HeatMap.png')), 1200, 900, pointsize=12,res = 150)
  print(Heatmap(toHMP,cluster_columns = F, cluster_rows = T,col=col_fun,top_annotation = column_ha,show_row_names = F,show_row_dend = F))
  dev.off()
  
  ### PREPARE DATA 
  subPhenData$Treatment <- as.factor(subPhenData$Treatment)
  phenoD <- AnnotatedDataFrame(data=subPhenData[,"Treatment",drop=FALSE])
  
  featD <- AnnotatedDataFrame(data=data.frame(row.names=row.names(expresionMatrix)))
  eset <- ExpressionSet(assayData = as.matrix(expresionMatrix),phenoData=phenoD,featureData = featD)
  exprs(eset) <- log2(exprs(eset))
  
  diffExpr <- data.frame(row.names=featureNames(eset))
  selCol <- unlist(pData(eset)$Treatment %in% c(control,mytreatment))
  fit <- eBayes(lmFit(eset, model.matrix(~factor(eset$Treatment), data=pData(eset))))
  diffExpr$pVal <- fit$p.value[,2]
  diffExpr$FDR <- p.adjust(diffExpr$pVal,method="BH")
  
  eCtrl <- subset(exprs(eset),select=which(pData(eset)$Treatment == control))
  eCase <- subset(exprs(eset),select=which(pData(eset)$Treatment == mytreatment))
  ratios <- apply(eCase,1,"median", na.rm=T) - apply(eCtrl,1,"median", na.rm=T)
  diffExpr$log2FC <- ratios
  
  volcano <- EnhancedVolcano(diffExpr, x='log2FC', y='FDR', lab=row.names(diffExpr),
                             title = paste0(datatype,crosslink),subtitle = '',
                             pCutoff = 0.05,FCcutoff = 1,
                             legendLabels = c('NS', expression(Log[2]~FC),
                                              'FDR', expression(p~and~log[2]~FC)))
  ggsave(file.path(outdir,paste0(datatype,crosslink,'_VulcanoPlot.png')), plot = volcano, 
         device = 'png', scale = 1, dpi = "print",width = 25,height = 25,units = 'cm')
  
  colnames(diffExpr) <- paste0(colnames(diffExpr),".",mytreatment)
  
  diffExprMapped <- merge(mappingDF,diffExpr,by.x = "UniprotACC",by.y = "row.names",all.y = T)
  diffExprMapped <- diffExprMapped[order(diffExprMapped[,paste0("FDR",".",mytreatment)]),]
  saveTablesTsvExc(diffExprMapped,outdir,completeNdedup=F,excel=T,bycompleteFC=F,rownames = F)
}
