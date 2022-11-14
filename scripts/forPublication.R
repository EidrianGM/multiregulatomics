wholeDFfile <- "/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/data/mapping/wholeDF.tsv"
wholeDF <- read.delim(wholeDFfile,quote = "")

newFAXSAFull <- read.xlsx("/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/adrian/newFAXSAFull.xlsx")
newUVXSAFull <- read.xlsx("/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/adrian/newUVXSAFull.xlsx")

newFAXSAFull$GeneName

wholeDF2 <- merge(wholeDF,newFAXSAFull)
wholeDF2 <- merge(wholeDF2,newUVXSAFull)

saveTablesTsvExc(wholeDF2,outdir = "/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/adrian",completeNdedup = F,excel = T,bycompleteFC = F,rownames = F)

################################################################################
# RNA SEQ Data
################################################################################
basedir <- "/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/adrian/forPublication"
datatype <- "Transcriptome"

#### Quality Control Plots of Expression
RNAseqMatrixFile <- "/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/data/RNA-Seq data/result/counts_genes.csv"
RNAseqMatrix <- read.delim(RNAseqMatrixFile,sep = ",")
expresionMatrix <- RNAseqMatrix[2:ncol(RNAseqMatrix)]
subPhenData <- data.frame(class=gsub("\\..*","",colnames(expresionMatrix)))

tagname <- data
outdir <- file.path(basedir,tagname); dir.create(outdir,recursive = T,showWarnings = F)
#someVisualization(expresionMatrix,subPhenData,outdir,tagname)

##### HEATMAP HIHGLIGHTS 
wholeDF2 <- read.delim("/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/adrian/wholeDF2.tsv",quote = "")

infoColumns <- c("DEGs.log2FC.DTT.YPD","DEGs.log2FC.H202.YPD","DEGs.log2FC.SA.YPD",
                 "proteomeFAX.log2ratio_Condition2","mRNABPomeFAX.log2ratio_PolyARNAFAXwithDTT","netChangesFAXDTT",
                 "proteomeFAX.log2ratio_Condition3","mRNABPomeFAX.log2ratio_PolyARNAFAXwithH2O2","netChangesFAXH2O2",
                 "newPROT.FAX.SA.log2FC.SA","mRNABPomeFAX.log2ratio_PolyARNAFAXwithSA","netChangesFAXSA",
                 "proteomeUV.log2ratio_Condition2","mRNABPomeUV.log2ratio_PolyARNAUVwithDTT","netChangesUVDTT",
                 "proteomeUV.log2ratio_Condition3","mRNABPomeUV.log2ratio_PolyARNAUVwithH202","netChangesUVH2O",
                 "proteomeUV.log2ratio_Condition4","newRBP.UVX.SA.log2FC.SA","netChangesUVSA")

pvalsCols <- c("DEGs.adj.pval.DTT.YPD","DEGs.adj.pval.H202.YPD","DEGs.adj.pval.SA.YPD",
               "proteomeFAX.qValue_Condition2","mRNABPomeFAX.qValue_PolyARNAFAXwithDTT","netChangesFAXDTT",
               "proteomeFAX.qValue_Condition3","mRNABPomeFAX.qValue_PolyARNAFAXwithH2O2","netChangesFAXH2O2",
               "newPROT.FAX.SA.FDR.SA","mRNABPomeFAX.qValue_PolyARNAFAXwithSA","netChangesFAXSA",
               "proteomeUV.qValue_Condition2","mRNABPomeUV.qValue_PolyARNAUVwithDTT","netChangesUVDTT",
               "proteomeUV.qValue_Condition3","mRNABPomeUV.qValue_PolyARNAUVwithH202","netChangesUVH2O",
               "proteomeUV.qValue_Condition4","newRBP.UVX.SA.FDR.SA","netChangesUVSA")

phenoData <- data.frame(datatype=c(rep("Transcriptome",3),rep(c("Proteome","RBPome","NetChages"),2*3)),
                        treatment=c(c("DTT","H2O2","SA"),rep(c(rep("DTT",3),rep("H2O2",3),rep("SA",3)),2)),
                        crosslink=c(rep("DEGs",3),rep("FAX",3*3),rep("UVX",3*3)))

colannot <- HeatmapAnnotation(datatype=phenoData$datatype,
                              crosslink=phenoData$crosslink,
                              treatment=phenoData$treatment)

######### Custom Genes 
# top2DiffGensSA + sphingo2glycerolipid + KATgenes + acetyltransferases
#top2DiffGensSA <- wholeDF2$GeneName[order(wholeDF2$DEGs.log2FC.SA.YPD,decreasing = T)][1:2]
#sphingo2glycerolipid <- c("DPL1","HFD1","FAA1","FAA4","FAA2","FAA3","FAT1","PCS60","ACS1","ACS2") # PCS60 == FAT2
#KATgenes <- wholeDF2$GeneName[grep("KAT",wholeDF2$Alias)]; cat(paste(KATgenes,collapse = ", "))
acetyltransferases <- wholeDF2$GeneName[grep("acetyltransferase",wholeDF2$Alias)]; cat(paste(acetyltransferases,collapse = ", "))
deacetylases <- wholeDF2$GeneName[grep("deacetylase",wholeDF2$Alias)]; cat(paste(deacetylases,collapse = ", "))
#highlightGenes <- na.omit(unique(c(top2DiffGensSA,sphingo2glycerolipid,KATgenes,acetyltransferases)))
highlightGenes <- na.omit(unique(c(acetyltransferases,deacetylases)))

pValstoHMP <- wholeDF2[wholeDF2$GeneName %in% highlightGenes,c("GeneName",pvalsCols)]
rownames(pValstoHMP) <- pValstoHMP$GeneName
pValstoHMP <- pValstoHMP[highlightGenes,]
pValstoHMP[,colnames(pValstoHMP)[grep("net",colnames(pValstoHMP))]] <- 1
pValstoHMP$GeneName <- NULL

toHMP <- wholeDF2[wholeDF2$GeneName %in% highlightGenes,c("GeneName","Protein.names",infoColumns)]
rownames(toHMP) <- paste(toHMP$GeneName,gsub(" \\(.*","",toHMP$Protein.names))
toHMP$GeneName <- NULL
toHMP$Protein.names <- NULL
toHMP <- toHMP[highlightGenes,]

acetylases <- highlightGenes %in% acetyltransferases 
acetylases[acetylases == T] <- "Acetylases"
acetylases[acetylases == F] <- "Deacetylases"

leftannot <- rowAnnotation(acetylases = acetylases)

outdir <- "/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/adrian/highlightedGenes"
tagname <- "HighlightedGenes"
createHMP(toHMP,pValstoHMP,phenoData$crosslink,colannot,outdir,tagname,acetylases,leftannot)

#### DEGS SA
SADEGs <- wholeDF2[(wholeDF2$DEGs.adj.pval.SA.YPD < 0.01 & abs(wholeDF2$DEGs.log2FC.SA.YPD) > 3),]
SADEGs <- unique(na.omit(SADEGs$GeneName)); length(SADEGs)

highlightGenes <- SADEGs
pValstoHMP <- wholeDF2[wholeDF2$GeneName %in% highlightGenes,c("GeneName",pvalsCols)]
if (any(duplicated(pValstoHMP$GeneName))){
  rownames(pValstoHMP) <- paste(pValstoHMP$GeneName,1:nrow(pValstoHMP))
}else{
  rownames(pValstoHMP) <- pValstoHMP$GeneName
}
pValstoHMP[,colnames(pValstoHMP)[grep("net",colnames(pValstoHMP))]] <- 1
pValstoHMP$GeneName <- NULL

toHMP <- wholeDF2[wholeDF2$GeneName %in% highlightGenes,c("GeneName","Protein.names",infoColumns)]
if (any(duplicated(toHMP$GeneName))){
  rownames(toHMP) <- paste(toHMP$GeneName,gsub(" \\(.*","",toHMP$Protein.names),1:nrow(toHMP))
}else{
  rownames(toHMP) <- toHMP$GeneName
}
toHMP$GeneName <- NULL
toHMP$Protein.names <- NULL

outdir <- "/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/adrian/highlightedGenes"
tagname <- "AllSADegs"
colsplit <- phenoData$crosslink
createHMP(toHMP,pValstoHMP,colsplit,colannot,outdir,tagname)

#### Whole RBPomes UVX + FAX
RBPomesGenes <- wholeDF2[!is.na(wholeDF2$newRBP.UVX.SA.log2FC.SA) & !is.na(wholeDF2$mRNABPomeFAX.log2ratio_PolyARNAFAXwithSA),]
RBPomesGenes <- unique(na.omit(RBPomesGenes$GeneName)); length(RBPomesGenes)
highlightGenes <- RBPomesGenes
pValstoHMP <- wholeDF2[wholeDF2$GeneName %in% highlightGenes,c("GeneName",pvalsCols)]
if (any(duplicated(pValstoHMP$GeneName))){
  rownames(pValstoHMP) <- paste(pValstoHMP$GeneName,1:nrow(pValstoHMP))
}else{
  rownames(pValstoHMP) <- pValstoHMP$GeneName
}
pValstoHMP[,colnames(pValstoHMP)[grep("net",colnames(pValstoHMP))]] <- 1
pValstoHMP$GeneName <- NULL

toHMP <- wholeDF2[wholeDF2$GeneName %in% highlightGenes,c("GeneName","Protein.names",infoColumns)]
if (any(duplicated(toHMP$GeneName))){
  rownames(toHMP) <- paste(toHMP$GeneName,gsub(" \\(.*","",toHMP$Protein.names),1:nrow(toHMP))
}else{
  rownames(toHMP) <- toHMP$GeneName
}
toHMP$GeneName <- NULL
toHMP$Protein.names <- NULL

outdir <- "/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/adrian/highlightedGenes"
tagname <- "RBPomes"
colsplit <- phenoData$crosslink
createHMP(toHMP,pValstoHMP,colsplit,colannot,outdir,tagname)

#### NUA4 Targets
NUA4targets <- read.delim("NUA4targets.txt",col.names = F)[,1]; length(NUA4targets)
all(NUA4targets %in% wholeDF2$ORF)
NUA4targets <- wholeDF2[wholeDF2$ORF %in% NUA4targets,]
NUA4targets <- unique(na.omit(NUA4targets$UniprotName)); length(NUA4targets)
highlightGenes <- NUA4targets
pValstoHMP <- wholeDF2[wholeDF2$UniprotName %in% highlightGenes,c("UniprotName",pvalsCols)]
if (any(duplicated(pValstoHMP$UniprotName))){
  rownames(pValstoHMP) <- paste(pValstoHMP$UniprotName,1:nrow(pValstoHMP))
}else{
  rownames(pValstoHMP) <- pValstoHMP$UniprotName
}
pValstoHMP[,colnames(pValstoHMP)[grep("net",colnames(pValstoHMP))]] <- 1
pValstoHMP$UniprotName <- NULL

toHMP <- wholeDF2[wholeDF2$UniprotName %in% highlightGenes,c("UniprotName","Protein.names",infoColumns)]
if (any(duplicated(toHMP$UniprotName))){
  rownames(toHMP) <- paste(toHMP$UniprotName,gsub(" \\(.*","",toHMP$Protein.names),1:nrow(toHMP))
}else{
  rownames(toHMP) <- toHMP$UniprotName
}
toHMP$UniprotName <- NULL
toHMP$Protein.names <- NULL

outdir <- "/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/adrian/highlightedGenes"
tagname <- "NUA4targets"
colsplit <- phenoData$crosslink
createHMP(toHMP,pValstoHMP,colsplit,colannot,outdir,tagname)

####
rowsplit

createHMP <- function(toHMP,pValstoHMP,colsplit,colannot,outdir,tagname,rowsplit=NA,leftannot=NA){
  # row_split = rowsplit, left_annotation = leftannot,  
  col_fun <- colorRamp2(c(min(toHMP,na.rm = T),0,max(toHMP,na.rm = T)), c("blue","white","red"))
  myHMPplot <- Heatmap(toHMP,cluster_columns = F, cluster_rows = T,col=col_fun,show_row_names = T,show_row_dend = F,
                       column_split = colsplit, name = "Log2FC",
                       top_annotation = colannot, cluster_row_slices = T,
                       row_names_gp = gpar(fontsize = 9),column_names_gp = gpar(fontsize = 9),
                       width = ncol(toHMP)*unit(3, "mm"), height = nrow(toHMP)*unit(3, "mm"),
                       cell_fun = function(j, i, x, y, w, h, f) {
                         gb = textGrob("*")
                         gb_w = convertWidth(grobWidth(gb), "mm")
                         gb_h = convertHeight(grobHeight(gb), "mm")
                         if (!is.na(pValstoHMP[i, j])){
                           if(pValstoHMP[i, j] < 0.01) {
                             grid.text("**", x,  y - gb_h*0.5 + gb_w*0.4)
                           } else if(pValstoHMP[i, j] < 0.05) {
                             grid.text("*", x,  y - gb_h*0.5 + gb_w*0.4)
                           }
                         }
                       }
  )
  png(file.path(outdir,paste0(tagname,'_HeatMap.png')), 1000, nrow(toHMP)*15, pointsize=12,res = "print")
  draw(myHMPplot, heatmap_legend_side="left", annotation_legend_side = "left")
  dev.off()
}


vulcanoCols <- c("DEGs.adj.pval.SA.YPD","DEGs.log2FC.SA.YPD")
getVulcano(wholeDF2,vulcanoCols,highlightGenes,outdir,"TranscriptomeSA",pvalcutoff=0.01,FCcutoff=1)
vulcanoCols <- c("newPROT.FAX.SA.FDR.SA","newPROT.FAX.SA.log2FC.SA")
getVulcano(wholeDF2,vulcanoCols,highlightGenes,outdir,"ProteomeFAXSA",pvalcutoff=0.05,FCcutoff=1)
vulcanoCols <- c("mRNABPomeFAX.qValue_PolyARNAFAXwithSA","mRNABPomeFAX.log2ratio_PolyARNAFAXwithSA")
getVulcano(wholeDF2,vulcanoCols,highlightGenes,outdir,"RBPomeFAXSA",pvalcutoff=0.05,FCcutoff=1)
vulcanoCols <- c("proteomeUV.qValue_Condition4","proteomeUV.log2ratio_Condition4")
getVulcano(wholeDF2,vulcanoCols,highlightGenes,outdir,"ProteomeUVXSA",pvalcutoff=0.05,FCcutoff=1)
vulcanoCols <- c("newRBP.UVX.SA.FDR.SA","newRBP.UVX.SA.log2FC.SA")
getVulcano(wholeDF2,vulcanoCols,highlightGenes,outdir,"RBPomeUVXSA",pvalcutoff=0.05,FCcutoff=1)


getVulcano <- function(wholeDF2,vulcanoCols,highlightGenes,outdir,tagname,pvalcutoff=0.05,FCcutoff=1){
  diffExpr <- wholeDF2[,c("GeneName",vulcanoCols)]
  diffExpr <- diffExpr[complete.cases(diffExpr),]
  colnames(diffExpr) <- c("GeneName","FDR","log2FC")
  volcano <- EnhancedVolcano(diffExpr, x='log2FC', y='FDR', lab=diffExpr$GeneName,
                             title = tagname,subtitle = '',
                             pCutoff = pvalcutoff,FCcutoff = FCcutoff,
                             legendLabels = c('log2FC', expression(Log[2]~FC),
                                              'FDR', expression(FDR~and~log[2]~FC)))
  #selectLab = highlightGenes,drawConnectors = TRUE,widthConnectors = 1.0,colConnectors = 'black',
  ggsave(file.path(outdir,paste0(tagname,'_VulcanoPlot.png')), plot = volcano, 
         device = 'png', scale = 1, dpi = "print",width = 25,height = 25,units = 'cm')
}


#### Differential Expression Plots
pvalcutoff <- 0.01; FCcutoff <- 1
wholeDF <- read.xlsx("/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/data/mapping/wholeDF.xlsx")
treatments <- c("SA","H2O2","DTT")

for (treatment in treatments){
  tagname <- paste(data,treatment,sep = "_")
  treatmentCols <- colnames(wholeDF)[grep(paste0("DEGs.*",treatment),colnames(wholeDF))]
  
  diffExpr <- wholeDF[,c("GeneName",treatmentCols)]
  colnames(diffExpr) <- c("GeneName","FDR","log2FC")
  volcano <- EnhancedVolcano(diffExpr, x='log2FC', y='FDR', lab=diffExpr$GeneName,
                             title = tagname,subtitle = '',
                             pCutoff = pvalcutoff,FCcutoff = FCcutoff,
                             legendLabels = c('log2FC', expression(Log[2]~FC),
                                              'FDR', expression(FDR~and~log[2]~FC)))
  ggsave(file.path(outdir,paste0(tagname,'VulcanoPlot.png')), plot = volcano, 
         device = 'png', scale = 1, dpi = "print",width = 25,height = 25,units = 'cm')
  
  getPlotsOfConditions(diffExpr,outdir,tagname,pvalcutoff,FCcutoff)  
}

#### Differential Expression Plots

