library(ComplexHeatmap)
library(pheatmap)
library(WGCNA)
library(openxlsx)
library(pheatmap)
library(circlize)
library(stringr)

setwd('/home/eidrian/Desktop/realServices/multiregulatomics')

### Read Data
FAXacceptedProts <- read.delim("FinalData/BackgroundRemoval/FAXAccBackGrAcceptedFC3.tsv",header = F)[,1]
UVXacceptedProts <- read.delim("FinalData/BackgroundRemoval/UVXAccBackGrAcceptedFC3.tsv",header = F)[,1]

wholeDFfile <- "FinalData/allDataDF.tsv"
wholeDF <- read.delim(wholeDFfile,quote = "")
wholeDF$UniprotName <- gsub("_YEAST","",wholeDF$UniprotName)

samplesNames <- colnames(wholeDF)[grep('.X0',colnames(wholeDF))]
runids <- as.numeric(sapply(strsplit(samplesNames,'X0|_'), function(x) return(x[2])))
datatypes <- sapply(strsplit(samplesNames,'ome'), function(x) return(x[1]))
datatypes <- paste0(gsub('RBP','mRBP',datatypes),'ome')
crosslinks <- paste0(sapply(strsplit(samplesNames,'ome|X'), function(x) return(x[2])),'X')

samplesNamesPhenoD <- as.data.frame(cbind(runid=runids,samplesName=samplesNames,datatype=datatypes,
                                          crosslink=crosslinks,uniqueid=paste0(datatypes,crosslinks,runids)))


phenoData <- read.xlsx('data/phenoData.xlsx')
phenoData <- phenoData[!(phenoData$Crosslinking == 'NOX'),]
phenoData <- phenoData[c('Data','Treatment','Crosslinking','Run.number','Sample.ID')]
phenoData$Data <- gsub(' ','',phenoData$Data)
phenoData$Treatment <- gsub(' ','',phenoData$Treatment)
phenoData$Crosslinking <- gsub(' ','',phenoData$Crosslinking)
phenoData <- phenoData[!duplicated(phenoData),]
phenoData$uniqueid <- paste0(phenoData$Data,phenoData$Crosslinking,phenoData$Run.number)

totalInfoPhenoD <- merge(phenoData,samplesNamesPhenoD,by='uniqueid')
nrow(samplesNamesPhenoD) == nrow(totalInfoPhenoD)
totalInfoPhenoD <- totalInfoPhenoD[c('samplesName',colnames(phenoData))]

degsCols <- colnames(wholeDF)[grep('brep',colnames(wholeDF))]
repnum <- sapply(strsplit(degsCols,'\\.|brep'), function(x) return(x[4]))
degsTreat <- sapply(strsplit(degsCols,'\\.'), function(x) return(x[2]))
degsPhenoD <- cbind(samplesName=degsCols,Data='Transcriptome',Treatment=degsTreat,Crosslinking='-',Run.number='-',Sample.ID=repnum,uniqueid='-')

totalInfoPhenoD <- rbind(degsPhenoD,totalInfoPhenoD)
write.table(totalInfoPhenoD,file = 'FinalData/colNamesPhenoData.tsv',sep = '\t',quote = F,row.names = F)
View(totalInfoPhenoD)

crosslink <- 'FAX'
treatment <- 'SA'
if(crosslink == 'FAX'){
  subDF <- wholeDF[which(wholeDF$proteinName %in% FAXacceptedProts),] 
}else{
  subDF <- wholeDF[which(wholeDF$proteinName %in% UVXacceptedProts),] 
}

totalInfoPhenoD$Sample.ID <- gsub('.*\\(|\\)','',totalInfoPhenoD$Sample.ID)

samplesSelection <- totalInfoPhenoD$Treatment == treatment & (totalInfoPhenoD$Crosslinking == crosslink | totalInfoPhenoD$Data == 'Transcriptome')
subDF <- subDF[c('DEGs.target','UniprotName',totalInfoPhenoD$samplesName[samplesSelection])]
subDF <- subDF[complete.cases(subDF),]
row.names(subDF) <- paste(subDF$DEGs.target,subDF$UniprotName,sep = '_')

subDF$DEGs.target <- NULL
subDF$UniprotName <- NULL

mycolNames <- paste(crosslink, treatment,
                    totalInfoPhenoD$Data[samplesSelection],
                    totalInfoPhenoD$Sample.ID[samplesSelection], sep = '_')

colnames(subDF) <- mycolNames

do.wgcna<-function(data,moduleMinSize=5,tomtype='signed'){
  data<-t(data)
  power=seq(1,50,by=1)  ## Calculamos poder para WGCNA
  RpowerTable=pickSoftThreshold(data, powerVector=power)[[2]] 
  
  plot(RpowerTable[,1], -sign(RpowerTable[,3])*RpowerTable[,2],xlab="
     Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed
     R^2",type="n")
  text(RpowerTable[,1], -sign(RpowerTable[,3])*RpowerTable[,2],
       labels=power,cex=1,col="red") 
  abline(h=0.95,col="red")
  jamon <- RpowerTable$SFT.R.sq
  jamon <- which(max(jamon[jamon < 0.9]) == jamon)

  power <- readline(prompt=paste0("Enter power (best=",jamon,"):"))
  # try selected the most closest to 90 but do not surpass it
  power<-as.integer(power)
  
  ## Calculo de modulos
  net = blockwiseModules(data, power = power, corType = "bicor",
                         networkType=tomtype,TOMType=tomtype,
                         minModuleSize = moduleMinSize,maxBlockSize = 200,
                         reassignThreshold = 0, mergeCutHeight = 0.25,
                         numericLabels = TRUE, pamRespectsDendro = FALSE,
                         saveTOMs = F,
                         saveTOMFileBase = "netTOM",
                         verbose = 3)
  return(net)
}

datExpr <- t(subDF)
power = seq(1, 70, by=1)  ## Calculamos poder para WGCNA
RpowerTable=pickSoftThreshold(datExpr, powerVector=power)[[2]] 
plot(RpowerTable[,1], -sign(RpowerTable[,3])*RpowerTable[,2],xlab="
     Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed
     R^2",type="n")
text(RpowerTable[,1], -sign(RpowerTable[,3])*RpowerTable[,2],
     labels=power,cex=1,col="red") 
abline(h=0.90,col="red")
jamon <- RpowerTable$SFT.R.sq
jamon <- which(max(jamon[jamon < 0.9]) == jamon)
softPower <- min(jamon,30)

doWGCAanalysis <- function(datExpr,tomtype,outdir,outName='',softPower=30,minModuleSize=5,maxBlockSize=200){
  net <- blockwiseModules(datExpr, power = softPower, corType = "bicor",
                         networkType=tomtype,TOMType=tomtype,
                         minModuleSize = minModuleSize, maxBlockSize = maxBlockSize,
                         reassignThreshold = 0, mergeCutHeight = 0.25,
                         numericLabels = TRUE, pamRespectsDendro = FALSE,
                         saveTOMs = F,
                         saveTOMFileBase = "netTOM",
                         verbose = 0)
  table(labels2colors(net$colors))
  
  MEList = moduleEigengenes(datExpr, colors = labels2colors(net$colors))
  MEs = MEList$eigengenes
  MEDiss = 1-cor(MEs);
  METree = hclust(as.dist(MEDiss), method = "average");
  
  outName <- file.path(outdir,paste0(outName,tomtype,'Modules'))
  png(paste0(outName,'wgcnaDendrogram.png'), 1200, 600, pointsize=20)
  print(plot(METree, main = "Clustering of module eigengenes", xlab = "", sub = ""))
  dev.off()
  
  resDF <- merge(t(datExpr),data.frame(cluster=net$colors),by = 'row.names')
  resDF <- cbind(resDF, colors=labels2colors(resDF$cluster))
  colnames(resDF)[0] <- 'genes'
  write.xlsx(resDF,paste0(outName,'wgcnaDF.xlsx'), asTable = FALSE, overwrite = TRUE,rowNames = F)

  toHMP <- resDF[order(resDF$cluster),]
  toAnnot <- toHMP$colors
  names(toAnnot) <- toHMP$colors
  toHMP <- toHMP[,grep('_',colnames(resDF))]
  
  heatmapannot <- rowAnnotation(cluster=toAnnot,col=list(cluster=toAnnot))
  toHMP <- t(scale(t(toHMP)))
  col_fun <- colorRamp2(c(min(toHMP),0,max(toHMP)), c("blue", "white","red"))
  png(paste0(outName,'wgcnaHeatmap.png'), 1200, 1200, pointsize = 20)
  print(Heatmap(toHMP,cluster_columns = F, cluster_rows = F,col=col_fun, 
          right_annotation=heatmapannot,row_labels=rep('',nrow(toHMP))))
  dev.off()
}


outdir <- 'clustering'
softPower <- 30
minModuleSize <- 5

tomtype <- 'signed'
doWGCAanalysis(datExpr,tomtype,outdir,softPower=30,minModuleSize=5,maxBlockSize=200)

tomtype <- 'unsigned'
doWGCAanalysis(datExpr,tomtype,outdir,softPower=30,minModuleSize=5,maxBlockSize=200)

tomtype <- 'signed'
doWGCAanalysis(log(datExpr),tomtype,outdir,outName = 'log', softPower=30,minModuleSize=5,maxBlockSize=200)

tomtype <- 'unsigned'
doWGCAanalysis(log(datExpr),tomtype,outdir,outName = 'log', softPower=30,minModuleSize=5,maxBlockSize=200)

#### Concrete HMP - COBOS SELECTION
outdir <- '/data/realServices/ej_cobos/results/wgcna/medula_again'
toHMPfile <- '/data/realServices/ej_cobos/to_Heatmap.xlsx'
### DRG
drgExpr <- read.xlsx(toHMPfile,sheet = 1)
drgExpr$Colors.joint

toHMP <- drgExpr[,grep('_',colnames(drgExpr))]
toHMP <- t(scale(t(toHMP)))
col_fun <- colorRamp2(c(min(toHMP),0,max(toHMP)), c("blue", "white","red"))
hmoutfile <- file.path(outdir,paste0('drgWGCNAHeatmap.png'))
png(hmoutfile, 1200, 1200, pointsize = 20)
Heatmap(toHMP,cluster_columns = F, cluster_rows = F,col=col_fun, 
        row_labels=rep('',nrow(toHMP)),split = str_to_title(drgExpr$Colors.joint),
        gap = unit(5, "mm"),row_title_gp = gpar(fontsize = 10),
        row_title_rot = switch("left", "left" = 0))
dev.off()
svg(gsub('png','svg',hmoutfile), 10, 10, pointsize = 20)
Heatmap(toHMP,cluster_columns = F, cluster_rows = F,col=col_fun, 
        row_labels=rep('',nrow(toHMP)),split = str_to_title(drgExpr$Colors.joint),
        gap = unit(3, "mm"),row_title_gp = gpar(fontsize = 10), column_names_gp = gpar(fontsize = 10), 
        row_title_rot = switch("left", "left" = 0),heatmap_legend_param = list(title = ''))
dev.off()
toEig <- t(as.matrix(drgExpr[,grep('_',colnames(drgExpr))]))
MEList = moduleEigengenes(toEig, colors = drgExpr$Colors.joint)
MEs = MEList$eigengenes
colnames(MEs) <- str_to_title(gsub('ME','',colnames(MEs)))
WriteXLS(MEs, file.path(outdir,'drgEigengenes.xlsx'), row.names = T,col.names = T)

outdir <- '/data/realServices/ej_cobos/reorderWGCNA/'
MEs <- read.xlsx('/data/realServices/ej_cobos/results/wgcna/medula_again/drgEigengenes.xlsx',rowNames = T)
hmoutfile <- file.path(outdir,'NEW_drgWGCNAEigCor.png')
#paste(colnames(MEs),collapse="','")

MEs <- MEs[c('Black+Blue+Purple','Turquoise+Yellow','Red','Brown+Green+Magenta+Pink')]

png(hmoutfile, 700, 700, pointsize = 20)
plotEigengeneNetworks(MEs,"", marHeatmap = c(11,11,5,5), cex.lab = 0.8,
                      xLabelsAngle= 90,cex.adjacency = 0.8,plotDendrograms = F,
                      signed = FALSE, heatmapColors=blueWhiteRed(50))
dev.off()
svg(gsub('png','svg',hmoutfile), 8, 8, pointsize = 20)
plotEigengeneNetworks(MEs,"", marHeatmap = c(11,11,5,5), cex.lab = 0.8,
                      xLabelsAngle= 90,cex.adjacency = 0.8,plotDendrograms = F,
                      signed = FALSE,heatmapColors=blueWhiteRed(50))
dev.off()

corME = abs(cor(MEs, use = "p"))
write.xlsx(as.data.frame(corME), file.path(outdir,'drgEigengenesPearsonCorMatrix.xlsx'), rowNames = T,colNames = T)

### MEDULA
medExpr <- read_excel(toHMPfile,sheet = 2,col_names = T,trim_ws = T)
medExpr <- medExpr[medExpr$...23 != "grey",]

toHMP <- medExpr[,grep('_',colnames(medExpr))]
toHMP <- t(scale(t(toHMP)))
col_fun <- colorRamp2(c(min(toHMP),0,max(toHMP)), c("blue", "white","red"))
hmoutfile <- file.path(outdir,paste0('medWGCNAHeatmap.png'))
png(hmoutfile, 1200, 1200, pointsize = 20)
Heatmap(toHMP,cluster_columns = F, cluster_rows = F,col=col_fun, 
        row_labels=rep('',nrow(toHMP)),split = str_to_title(medExpr$...23),
        gap = unit(5, "mm"),row_title_gp = gpar(fontsize = 10),
        row_title_rot = switch("left", "left" = 0))
dev.off()
svg(gsub('png','svg',hmoutfile), 10, 10, pointsize = 20)
Heatmap(toHMP,cluster_columns = F, cluster_rows = F,col=col_fun, 
        row_labels=rep('',nrow(toHMP)),split = str_to_title(medExpr$...23),
        gap = unit(3, "mm"),row_title_gp = gpar(fontsize = 10), column_names_gp = gpar(fontsize = 10), 
        row_title_rot = switch("left", "left" = 0),heatmap_legend_param = list(title = ''))
dev.off()

toEig <- t(as.matrix(medExpr[,grep('_',colnames(medExpr))]))
MEList = moduleEigengenes(toEig, colors = medExpr$...23)
MEs = MEList$eigengenes
colnames(MEs) <- str_to_title(gsub('ME','',colnames(MEs)))

WriteXLS(MEs, file.path(outdir,'corEigengenes.xlsx'), row.names = T,col.names = T)

colnames(MEs) <- str_to_title(gsub('ME','',colnames(MEs)))
WriteXLS(MEs, file.path(outdir,'medEigengenes.xlsx'), row.names = T,col.names = T)

outdir <- '/data/realServices/ej_cobos/reorderWGCNA/'
MEs <- read.xlsx('/data/realServices/ej_cobos/reorderWGCNA/medEigengenes.xlsx',rowNames = T)
hmoutfile <- file.path(outdir,'NEW_medWGCNAEigCor.png')
MEs <- MEs[c('Turquoise+Red+Blue+Green','Yellow','Brown')]
#paste(colnames(MEs),collapse="','")
png(hmoutfile, 700, 700, pointsize = 20)
plotEigengeneNetworks(MEs,"", marHeatmap = c(11,11,5,5), cex.lab = 0.8,
                      xLabelsAngle= 90,cex.adjacency = 0.8,plotDendrograms = F,
                      signed = FALSE,heatmapColors=blueWhiteRed(50))
dev.off()
svg(gsub('png','svg',hmoutfile), 8, 8, pointsize = 20)
plotEigengeneNetworks(MEs,"", marHeatmap = c(11,11,5,5), cex.lab = 0.8,
                      xLabelsAngle= 90,cex.adjacency = 0.8,plotDendrograms = F,
                      signed = FALSE,heatmapColors=blueWhiteRed(50))
dev.off()

corME = abs(cor(MEs, use = "p"))
write.xlsx(as.data.frame(corME), file.path(outdir,'medEigengenesPearsonCorMatrix.xlsx'), rowNames = T,colNames = T)

plotEigengeneNetworks(MEs,"", marHeatmap = c(11,11,5,5), cex.lab = 0.8,
                      xLabelsAngle= 90,cex.adjacency = 0.8,plotDendrograms = F,
                      signed = FALSE,heatmapColors=blueWhiteRed(50),printAdjacency = T)

