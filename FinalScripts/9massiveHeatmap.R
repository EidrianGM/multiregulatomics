library(ComplexHeatmap)
library(circlize)

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

### Read Data

wholeDFfile <- "/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/FinalData/allDataDF.tsv"
wholeDF <- read.delim(wholeDFfile,quote = "")

FAXacceptedProts <- read.delim("FinalData/BackgroundRemoval/FAXAccBackGrAcceptedFC3.tsv",header = F)[,1]
UVXacceptedProts <- read.delim("FinalData/BackgroundRemoval/UVXAccBackGrAcceptedFC3.tsv",header = F)[,1]

wholeDF <- wholeDF[wholeDF$proteinName %in% c(FAXacceptedProts,UVXacceptedProts),]

outdir <- "/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/FinalData/HeatMaps"


##### HEATMAP HIHGLIGHTS 
colnames(wholeDF)[grep("log|net",colnames(wholeDF))]
infoColumns <- c("DEGs.log2FC.H202.YPD",                  "DEGs.log2FC.DTT.YPD",                     "DEGs.log2FC.SA.YPD",                     
                 "ProteomeFAX.log2ratio_FAXwithDTT",      "RBPomeFAX.log2ratio_PolyARNAFAXwithDTT",  "FAXnetchangesDTT",                         
                 "ProteomeFAX.log2ratio_FAXwithH2O2",     "RBPomeFAX.log2ratio_PolyARNAFAXwithH2O2", "FAXnetchangesH2O2", 
                 "ProteomeFAX.log2ratio_FAXwithSA",       "RBPomeFAX.log2ratio_PolyARNAFAXwithSA",   "FAXnetchangesSA",                        
                 "ProteomeUVX.log2ratio_Condition2",      "RBPomeUVX.log2ratio_PolyARNAUVwithDTT",   "UVXnetchangesDTT",
                 "ProteomeUVX.log2ratio_Condition3",      "RBPomeUVX.log2ratio_PolyARNAUVwithH202",  "UVXnetchangesH2O2",
                 "ProteomeUVX.log2ratio_Condition4",      "RBPomeUVX.log2ratio_PolyARNAUVwithSA",    "UVXnetchangesSA")

colnames(wholeDF)[grep("qVal|adj",colnames(wholeDF))]
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

colannot <- HeatmapAnnotation(datatype=phenoData$datatype,
                              crosslink=phenoData$crosslink,
                              treatment=phenoData$treatment)

#### Whole RBPomes UVX + FAX
highlightGenes <- which(wholeDF$RBPomeFAX.qValue_PolyARNAFAXwithSA < 0.05 | wholeDF$RBPomeUVX.qValue_PolyARNAUVwithSA < 0.05)

length(highlightGenes)

pValstoHMP <- wholeDF[highlightGenes,c("proteinName",pvalsCols)]
toHMP <- wholeDF[highlightGenes,c("proteinName",infoColumns)]
pValstoHMP <- pValstoHMP[!duplicated(pValstoHMP),]
toHMP <- toHMP[!duplicated(toHMP),]

pValstoHMP[,colnames(pValstoHMP)[grep("netc",colnames(pValstoHMP))]] <- 1
#rownames(pValstoHMP) <- paste0(pValstoHMP$proteinName,1:nrow(pValstoHMP))
rownames(pValstoHMP) <- pValstoHMP$proteinName
pValstoHMP$proteinName <- NULL

#rownames(toHMP) <- paste0(toHMP$proteinName,1:nrow(toHMP))
rownames(toHMP) <- toHMP$proteinName
toHMP$proteinName <- NULL

colsplit <- phenoData$crosslink
tagname <- "SARBPomes"

## Higlight all Fatty Acid Genes

#createHMP(toHMP,pValstoHMP,colsplit,colannot,outdir,tagname)
col_fun <- colorRamp2(c(min(toHMP,na.rm = T),0,max(toHMP,na.rm = T)), c("blue","white","red"))
myHMPplot <- Heatmap(toHMP,cluster_columns = F, cluster_rows = F,col=col_fun,show_row_names = T,show_row_dend = F,
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

tiff(file.path(outdir,paste0(tagname,'_HeatMap.tiff')), 1000, nrow(toHMP)*10, pointsize=12,res = "print")
draw(myHMPplot, heatmap_legend_side="left", annotation_legend_side = "left")
dev.off()

pdf(file.path(outdir,paste0(tagname,'_HeatMap.pdf')), height = nrow(toHMP) / 10 + 10)
draw(myHMPplot, heatmap_legend_side="left", annotation_legend_side = "left")
dev.off()


