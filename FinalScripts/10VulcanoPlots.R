library(EnhancedVolcano)
source("FinalScripts/functionOmics.R")

wholeDFFile <- "FinalData/allDataDF.tsv"
wholeDF <- read.delim(wholeDFFile,quote = "")

FAXacceptedProts <- read.delim("FinalData/BackgroundRemoval/FAXAccBackGrAcceptedFC3.tsv",header = F)[,1]
UVXacceptedProts <- read.delim("FinalData/BackgroundRemoval/UVXAccBackGrAcceptedFC3.tsv",header = F)[,1]

#wholeDF <- wholeDF[wholeDF$proteinName %in% c(FAXacceptedProts,UVXacceptedProts),]

fcColumns <- c("DEGs.log2FC.H202.YPD",                  "DEGs.log2FC.DTT.YPD",                     "DEGs.log2FC.SA.YPD",                     
               "ProteomeFAX.log2ratio_FAXwithDTT",      "RBPomeFAX.log2ratio_PolyARNAFAXwithDTT",  "FAXnetchangesDTT",                         
               "ProteomeFAX.log2ratio_FAXwithH2O2",     "RBPomeFAX.log2ratio_PolyARNAFAXwithH2O2", "FAXnetchangesH2O2", 
               "ProteomeFAX.log2ratio_FAXwithSA",       "RBPomeFAX.log2ratio_PolyARNAFAXwithSA",   "FAXnetchangesSA",                        
               "ProteomeUVX.log2ratio_Condition2",      "RBPomeUVX.log2ratio_PolyARNAUVwithDTT",   "UVXnetchangesDTT",
               "ProteomeUVX.log2ratio_Condition3",      "RBPomeUVX.log2ratio_PolyARNAUVwithH202",  "UVXnetchangesH2O2",
               "ProteomeUVX.log2ratio_Condition4",      "RBPomeUVX.log2ratio_PolyARNAUVwithSA",    "UVXnetchangesSA")

pvalsCols <- c("DEGs.adj.pval.H202.YPD",         "DEGs.adj.pval.DTT.YPD",               "DEGs.adj.pval.SA.YPD",                
               "ProteomeFAX.qValue_FAXwithDTT",  "RBPomeFAX.qValue_PolyARNAFAXwithDTT", "RBPomeFAX.qValue_PolyARNAFAXwithDTT",                                   
               "ProteomeFAX.qValue_FAXwithH2O2", "RBPomeFAX.qValue_PolyARNAFAXwithH2O2", "RBPomeFAX.qValue_PolyARNAFAXwithH2O2",         
               "ProteomeFAX.qValue_FAXwithSA",   "RBPomeFAX.qValue_PolyARNAFAXwithSA", "RBPomeFAX.qValue_PolyARNAFAXwithSA",                                  
               "ProteomeUVX.qValue_Condition2",  "RBPomeUVX.qValue_PolyARNAUVwithDTT", "RBPomeUVX.qValue_PolyARNAUVwithDTT",        
               "ProteomeUVX.qValue_Condition3",  "RBPomeUVX.qValue_PolyARNAUVwithH202", "RBPomeUVX.qValue_PolyARNAUVwithH202",        
               "ProteomeUVX.qValue_Condition4",  "RBPomeUVX.qValue_PolyARNAUVwithSA", "RBPomeUVX.qValue_PolyARNAUVwithSA") 

phenoData <- data.frame(datatype=c(rep("Transcriptome",3),rep(c("Proteome","RBPome","NetChanges"),2*3)),
                        treatment=c(c("DTT","H2O2","SA"),rep(c(rep("DTT",3),rep("H2O2",3),rep("SA",3)),2)),
                        crosslink=c(rep("DEGs",3),rep("FAX",3*3),rep("UVX",3*3)))

subDF <- completeNdedup(wholeDF[,c(fcColumns[1],pvalsCols[1])])
colnames(subDF) <- c("log2FoldChange","FDR")

limsData <- c()

SELpvalsCols <- pvalsCols[grep("SA|Condition4",fcColumns)]
SELfcColumns <- fcColumns[grep("SA|Condition4",fcColumns)]
subphenoData <- phenoData[phenoData$treatment == 'SA',]

mysubDF <- wholeDF[which(wholeDF$proteinName %in% UVXacceptedProts),c('proteinName',SELfcColumns[idx],SELpvalsCols[idx])]

length(mysubDF$proteinName)
length(unique(mysubDF$proteinName))

for (idx in 1:length(SELfcColumns)){
  outName <- paste0(rev(subphenoData[idx,]),collapse = "")
  outdir <- "FinalData/EnrichmentResults/Visualizations/VulcanoPlots"
  if (grepl("Transc",outName)){
    pvalCutOff <- 0.01
    FCcutoff <- 3
  }else{
    pvalCutOff <- 0.05
    FCcutoff <- 1
  }
  if (grepl("FAX", outName)){
    subDF <- completeNdedup(wholeDF[which(wholeDF$proteinName %in% FAXacceptedProts),c(SELfcColumns[idx],SELpvalsCols[idx])])
  }else if(grepl("UVX", outName)){
    subDF <- completeNdedup(wholeDF[which(wholeDF$proteinName %in% UVXacceptedProts),c(SELfcColumns[idx],SELpvalsCols[idx])])
  }else{
    subDF <- completeNdedup(wholeDF[,c(SELfcColumns[idx],SELpvalsCols[idx])])
  }
  if (grepl("NetChanges", outName)){
    xLabel <- "NetChanges"
    yLabel <- "mRBPFDR"
  }else{
    xLabel <- "log2FoldChange"
    yLabel <- "FDR"
  }
  colnames(subDF) <- c(xLabel,yLabel)
  if (grepl("SARBPome", outName)){
    lim <- 11
  }else{
    lim <- ceiling(max(abs(c(max(subDF[,xLabel]),min(subDF[,xLabel])))))
  }
  outName <- paste(outName,yLabel,pvalCutOff,xLabel,FCcutoff)
  outFile <- file.path(outdir,gsub(' ','',paste0(outName,'_vulcanoplot')))
  limsData <- rbind(limsData,c(outName,lim))
  volcano <- EnhancedVolcano(subDF,x = xLabel,y = yLabel, lab = outName, 
                             title = outName, subtitle = '',
                             FCcutoff = FCcutoff, pCutoff=pvalCutOff, selectLab = F,
                             legendLabels = c("No Cutoff", xLabel, yLabel, paste0(xLabel,"&",yLabel)),
                             legendLabSize = 10,  legendIconSize = 3,
                             xlab = xLabel,ylab = yLabel,
                             xlim = c(-lim,lim), axisLabSize = 10,  titleLabSize = 12,
                             subtitleLabSize = 10, captionLabSize = 10,pointSize = 1, labSize = 10) 
  ggsave(paste0(outFile,".tiff"), plot = volcano, device = 'tiff', scale = 1, dpi = 300,width = 17,height = 15,units = 'cm')
}


