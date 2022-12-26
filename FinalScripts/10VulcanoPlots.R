library(EnhancedVolcano)
source("FinalScripts/functionOmics.R")


wholeDFFile <- "FinalData/allDataDF.tsv"
wholeDF <- read.delim(wholeDFFile,quote = "")

FAXacceptedProts <- read.delim("FinalData/BackgroundRemoval/FAXAccBackGrAcceptedFC3.tsv",header = F)[,1]
UVXacceptedProts <- read.delim("FinalData/BackgroundRemoval/UVXAccBackGrAcceptedFC3.tsv",header = F)[,1]

wholeDF <- wholeDF[wholeDF$proteinName %in% c(FAXacceptedProts,UVXacceptedProts),]


fcColumns <- c("DEGs.log2FC.H202.YPD",                  "DEGs.log2FC.DTT.YPD",                     "DEGs.log2FC.SA.YPD",                     
               "ProteomeFAX.log2ratio_FAXwithDTT",      "RBPomeFAX.log2ratio_PolyARNAFAXwithDTT",  "FAXnetchangesDTT",                         
               "ProteomeFAX.log2ratio_FAXwithH2O2",     "RBPomeFAX.log2ratio_PolyARNAFAXwithH2O2", "FAXnetchangesH2O2", 
               "ProteomeFAX.log2ratio_FAXwithSA",       "RBPomeFAX.log2ratio_PolyARNAFAXwithSA",   "FAXnetchangesSA",                        
               "ProteomeUVX.log2ratio_Condition2",      "RBPomeUVX.log2ratio_PolyARNAUVwithDTT",   "UVXnetchangesDTT",
               "ProteomeUVX.log2ratio_Condition3",      "RBPomeUVX.log2ratio_PolyARNAUVwithH202",  "UVXnetchangesH2O2",
               "ProteomeUVX.log2ratio_Condition4",      "RBPomeUVX.log2ratio_PolyARNAUVwithSA",    "UVXnetchangesSA")

pvalsCols <- c("DEGs.adj.pval.H202.YPD",         "DEGs.adj.pval.DTT.YPD",               "DEGs.adj.pval.SA.YPD",                
               "ProteomeFAX.qValue_FAXwithDTT",  "RBPomeFAX.qValue_PolyARNAFAXwithDTT", "FAXnetchangesDTT",                                   
               "ProteomeFAX.qValue_FAXwithH2O2", "RBPomeFAX.qValue_PolyARNAFAXwithH2O2", "FAXnetchangesH2O2",         
               "ProteomeFAX.qValue_FAXwithSA",   "RBPomeFAX.qValue_PolyARNAFAXwithSA", "FAXnetchangesSA",                                  
               "ProteomeUVX.qValue_Condition2",  "RBPomeUVX.qValue_PolyARNAUVwithDTT", "UVXnetchangesDTT",        
               "ProteomeUVX.qValue_Condition3",  "RBPomeUVX.qValue_PolyARNAUVwithH202", "UVXnetchangesH2O2",        
               "ProteomeUVX.qValue_Condition4",  "RBPomeUVX.qValue_PolyARNAUVwithSA", "UVXnetchangesSA") 

phenoData <- data.frame(datatype=c(rep("Transcriptome",3),rep(c("Proteome","RBPome","NetChages"),2*3)),
                        treatment=c(c("DTT","H2O2","SA"),rep(c(rep("DTT",3),rep("H2O2",3),rep("SA",3)),2)),
                        crosslink=c(rep("DEGs",3),rep("FAX",3*3),rep("UVX",3*3)))

subDF <- completeNdedup(wholeDF[,c(fcColumns[1],pvalsCols[1])])
colnames(subDF) <- c("log2FoldChange","FDR")

for (idx in 1:length(fcColumns)){
  outName <- paste0(rev(phenoData[idx,]),collapse = "")
  outdir <- "FinalData/EnrichmentResults/Visualizations/VulcanoPlots"
  if (grepl("Transc",outName)){
    pvalCutOff <- 0.01
    FCcutoff <- 3
  }else{
    pvalCutOff <- 0.05
    FCcutoff <- 1
  }
  outName <- paste0(outName,"FDR",pvalCutOff,"FC",FCcutoff)
  outFile <- file.path(outdir,paste0(outName,'_vulcanoplot'))
  subDF <- completeNdedup(wholeDF[,c(fcColumns[idx],pvalsCols[idx])])
  colnames(subDF) <- c("log2FoldChange","FDR")
  volcano <- EnhancedVolcano(subDF,x = 'log2FoldChange',y = "FDR",
                             lab = outName, title = outName,subtitle = '',
                             FCcutoff = FCcutoff, pCutoff=pvalCutOff, selectLab = F,
                             legendLabels = c('FC', expression(Log[2]~FC),
                                              "FDR", expression(FDR~and~log[2]~FC)))
  ggsave(paste0(outFile,".tiff"), plot = volcano, device = 'tiff', scale = 1, dpi = 300,width = 25,height = 25,units = 'cm')
}


