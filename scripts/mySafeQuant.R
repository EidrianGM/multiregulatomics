setwd("/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/")
## 1 Control de Calidad
## 2 se eliminan las muestras y se repiten pasos de QC
## Imputar sobre los valores crudos ya que tenemos casos y controles juntos 
## Imputar con los valores sin colapsar
## Imputar con todas las clases juntas siempre que se hayan medido con el mismo aparato
## Se imputa sobre los normalizados solo si viene de distintos estudios

inFile <- "/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/data/raw/MQtoSQ_UV.csv"
tagname <- "RBPomeUVX"
getRawQC(inFile,tagname)

tagname <- "RBPomeUVXFiltered"
samples2remove <- c(10,07,16,03)
getRawQC(inFile,tagname,samples2remove)

inFile <- "/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/data/raw/MQtoSQ_FAX.csv"
tagname <- "RBPomeFAX"
getRawQC(inFile,tagname)

inFile <- "/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/data/raw/MQtoSQ.csv"
tagname <- "NOX_RBPome"
getRawQC(inFile,tagname)

inFile <- "/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/data/raw/SNtoSQ-Prot40_SafeQuant_Input.csv"
tagname <- "ProteomeFAX"
crosslink <- "FAX"
getRawQC(inFile,tagname,"FAX")
tagname <- "ProteomeUVX"
getRawQC(inFile,tagname,"UV")
tagname <- "ProteomeNOX"
getRawQC(inFile,tagname,"NOX")


myDF <- read.csv(inFile,allowEscapes=T, check.names=F)

getRawQC <- function(inFile,tagname,crosslink=NA,samples2remove=NA){
  myDF <- read.csv(inFile,allowEscapes=T, check.names=F)
  
  outdir <- "/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/trial"
  outdir <- file.path(outdir,tagname);  dir.create(outdir,showWarnings = F)
  logfile <- file.path(outdir,paste0(tagname,".log")); write(paste("##",tagname),logfile,append = F)
  
  normDCol <- grep("Norm",colnames(myDF)) # Raw Normalized Values are not the ones!
  rawDCol <- grep("Raw",colnames(myDF))
  spcDCol <- grep("Spec",colnames(myDF))
  normRaw <- myDF[normDCol:(rawDCol-1)] 
  expresionMatrix <- myDF[rawDCol:(spcDCol-1)] 
  expresionMatrix[2,] <- gsub("_.*","",expresionMatrix[2,])
  subPhenData <- as.data.frame(t(expresionMatrix[1:2,])); colnames(subPhenData) <- c("class","id"); rownames(subPhenData) <- subPhenData[,2]
  colnames(expresionMatrix) <- expresionMatrix[2,]; expresionMatrix <- expresionMatrix[3:nrow(expresionMatrix),]; expresionMatrix <- type.convert(expresionMatrix)
  colnames(normRaw) <- normRaw[2,]; normRaw <- normRaw[3:nrow(normRaw),]; normRaw <- type.convert(normRaw)
  subPhenData$class <- gsub(".*\\s","",subPhenData$class); conditionCol <- "class"
  subPhenData[,conditionCol] <- as.factor(subPhenData[,conditionCol])
  metaData <- myDF[,1:(normDCol-1)];  colnames(metaData) <- metaData[2,]; metaData <- metaData[3:nrow(metaData),]

  if (class(samples2remove) != "logical"){
    subPhenData <- subPhenData[!(as.character(as.numeric(subPhenData$id)) %in% as.character(samples2remove)),]
    expresionMatrix <- expresionMatrix[subPhenData$id]
  }

  if (class(crosslink) != "logical"){
    samplesSelected <- subPhenData$id[grep(crosslink,subPhenData$class)]
    subPhenData <- subPhenData[grep(crosslink,subPhenData$class),]
    expresionMatrix <- expresionMatrix[samplesSelected]
  }

  write(paste("Number of Peptides Measures:",length(expresionMatrix)),logfile,append = T)
  write(paste("Number of Missing Peptides Measures:",sum(is.na(expresionMatrix)), 
              round(sum(is.na(expresionMatrix)) / (nrow(expresionMatrix) * ncol(expresionMatrix)) * 100,2)
  ),logfile,append = T)
  
  write(paste0("Number of Samples: ",ncol(expresionMatrix)),logfile,append = T)
  write.table(t(subPhenData),logfile,append = T,sep = " ",quote = F,col.names = F, row.names = F)
  write.table(t(table(subPhenData$class)),logfile,append = T,sep = " ",quote = F,col.names = T, row.names = F)

  nPeptides <- length(metaData$Accession)
  nParalogsPep <- length(metaData$Accession[grep(";",metaData$Accession)])
  uniqueProtsNParalogues <- length(unique(metaData$Accession))
  uniqueParalogues <- length(unique(metaData$Accession[grep(";",metaData$Accession)]))
  uniqueProts <- uniqueProtsNParalogues - uniqueParalogues
  write(paste("Number of Peptides:",nPeptides),logfile,append = T)
  write(paste("Number of Unique Protein Peptides:",nPeptides - nParalogsPep),logfile,append = T)
  write(paste("Number of Paralogues Peptides:",nParalogsPep),logfile,append = T)
  write(paste("Number of Unique Proteins and Paralogues:",uniqueProtsNParalogues),logfile,append = T)
  write(paste("Number of Unique Proteins:",uniqueProtsNParalogues - uniqueParalogues),logfile,append = T)
  write(paste("Number of Unique Paralogues:",uniqueParalogues),logfile,append = T)
  
  minPepsPerProt <- 2
  pepPerProt <- table(unique(data.frame(metaData$Accession,metaData$Sequence))[,1])
  write(paste0("Filtered proteins with < 2 peptides: ",sum(pepPerProt < minPepsPerProt)),logfile,append = T)
  pepPerProt <- pepPerProt[metaData$Accession] 
  protsWithMinPeps <- pepPerProt > minPepsPerProt  
  write(paste0("Filtered Peptides: ",sum(!protsWithMinPeps)),logfile,append = T)
  
  SamplesMissingPeptides <- apply(is.na(expresionMatrix),1,sum)
  PeptidesMissingSamples <- round(apply(is.na(expresionMatrix),2,sum) / nrow(expresionMatrix) * 100,2)
  peptidesWithMeasure <- SamplesMissingPeptides != ncol(expresionMatrix)
  
  write("% of Missing Peptides Per Sample:",logfile,append = T)
  write.table(t(PeptidesMissingSamples),logfile,append = T,sep = " ",quote = F,col.names = T, row.names = F)
  write("Missing Samples Per Peptide:",logfile,append = T)
  SamplesMissingPeptidesTable <- t(as.data.frame(table(SamplesMissingPeptides)))
  row.names(SamplesMissingPeptidesTable) <- c("N.Samples","Peptides.Missing")
  write.table(SamplesMissingPeptidesTable,logfile,append = T,sep = " ",quote = F,col.names = F,row.names = T)
  write(paste("Removing Peptides Without Any Measure:", sum(!peptidesWithMeasure)),logfile,append = T)
  
  expresionMatrixFiltered <- expresionMatrix[(protsWithMinPeps | minPepsPerProt),]
  metaDataFiltered <- metaData[(protsWithMinPeps | minPepsPerProt),]
  normRawFiltered <- normRaw[(protsWithMinPeps | minPepsPerProt),]
  
  dim(expresionMatrix)
  sum(protsWithMinPeps)
  
  length(unique(metaData$Accession))
  length(unique(metaDataFiltered$Accession))
  
  length()

  featD <- AnnotatedDataFrame(data=metaData)
  phenoD <- AnnotatedDataFrame(data=subPhenData)
  eset <- ExpressionSet(assayData = as.matrix(expresionMatrix),phenoData=phenoD,featureData = featD)
  
  write.table(subPhenData,file.path(outdir,"phenoData.tsv"),append = T,sep = "\t",quote = F,col.names = T, row.names = F)
  write.table(cbind(metaData,expresionMatrix),file.path(outdir,"RawData.tsv"),append = T,sep = "\t",quote = F,col.names = T, row.names = F)
  write.table(cbind(metaData,normRaw),file.path(outdir,"NormRawData.tsv"),append = T,sep = "\t",quote = F,col.names = T, row.names = F)
  
  
  eD <- exprs(eset); method <- "median"
  rawDataIdx <- apply(eD,2, FUN=method, na.rm=T)
  globalNormFactors <- as.numeric(rawDataIdx[1]) / as.numeric(rawDataIdx)
  pData(eset)$globalNormFactors <- globalNormFactors
  normData <- sweep(exprs(eset), 2, globalNormFactors, FUN="*")
  colnames(normData) <- sampleNames(eset)
  esetNorm <- ExpressionSet(assayData=normData,phenoData=phenoD,featureData=featD)
  expresionMatrix <- normData
  write.table(cbind(metaData,normData),file.path(outdir,"NormData.tsv"),append = F,sep = "\t",quote = F,col.names = T, row.names = F)
  
  someVisualization(expresionMatrix,outdir,paste0("Norm_Peptides_",tagname),subPhenData)
}
gminImputed <- sqImpute(esetNorm,method = "gmin")
gminImputed
View(as.data.frame(exprs(gminImputed)))
fenD(gminImputed)
collaExpMat <- as.data.frame(summarise_all(group_by(as.data.frame(exprs(gminImputed)),fData(gminImputed)$Accession),sum))
# someVisualization(exprs(gminImputed),outdir,paste0("Imputed_Uncollapsed_",tagname),subPhenData)
# 
# collaImpExpMat <- as.data.frame(summarise_all(group_by(as.data.frame(exprs(gminImputed)),metaData$Accession),sum))
# collaImpExpMat$`metaData$Accession` <- NULL
# 
# someVisualization(collaImpExpMat,outdir,paste0("Imputed_Collapsed",tagname),subPhenData)
# 

removeOutliers <- function(x, na.rm = TRUE, ...){
  qnt <- quantile(x, probs=c(.25, .75), na.rm = na.rm, ...)
  H <- 1.5 * IQR(x, na.rm = na.rm)
  y <- x
  y[x < (qnt[1] - H)] <- NA
  y[x > (qnt[2] + H)] <- NA
  
  return(y)
}