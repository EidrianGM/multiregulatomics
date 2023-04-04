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
wholeDF$UniprotName <- gsub("_YEAST","",wholeDF$UniprotName)

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

#### Whole RBPomes UVX + FAX
#highlightGenes <- which(wholeDF$RBPomeFAX.qValue_PolyARNAFAXwithSA < 0.05 | wholeDF$RBPomeUVX.qValue_PolyARNAUVwithSA < 0.05)

pvalcutoff <- 0.05
# c("FAX", "UVX", 'both')
crosslinks <- c('both')
fccutoff <- 1
for (crosslink in crosslinks){
  for (fccutoff in 1:3){
    if (crosslink == "UVX"){
      highlightGenes <- which(wholeDF$RBPomeUVX.qValue_PolyARNAUVwithSA < pvalcutoff & abs(wholeDF$RBPomeUVX.log2ratio_PolyARNAUVwithSA) > fccutoff)
      infoColumnsSel <- infoColumns[grep("UVX|DEGs",infoColumns)]
      pvalsColsSel <- pvalsCols[grep("UVX|DEGs",pvalsCols)]
      subphenoData <- phenoData[grep("UVX|DEGs",phenoData$crosslink),]
      tagname <- "UVXSARBPomes" 
    } else if (crosslink == "FAX"){
      highlightGenes <- which(wholeDF$RBPomeFAX.qValue_PolyARNAFAXwithSA < pvalcutoff & abs(wholeDF$RBPomeFAX.log2ratio_PolyARNAFAXwithSA) > fccutoff) 
      infoColumnsSel <- infoColumns[grep("FAX|DEGs",infoColumns)]
      pvalsColsSel <- pvalsCols[grep("FAX|DEGs",pvalsCols)]
      subphenoData <- phenoData[grep("FAX|DEGs",phenoData$crosslink),]
      tagname <- "FAXSARBPomes"
    }else{
      highlightGenes <- which((wholeDF$RBPomeFAX.qValue_PolyARNAFAXwithSA < pvalcutoff & abs(wholeDF$RBPomeFAX.log2ratio_PolyARNAFAXwithSA) > fccutoff) 
                              |(wholeDF$RBPomeUVX.qValue_PolyARNAUVwithSA < pvalcutoff & abs(wholeDF$RBPomeUVX.log2ratio_PolyARNAUVwithSA) > fccutoff))
      infoColumnsSel <- infoColumns
      pvalsColsSel <- pvalsCols
      subphenoData <- phenoData
      tagname <- "SARBPomes"
    }
    
    colannot <- HeatmapAnnotation(datatype=subphenoData$datatype,
                                  treatment=subphenoData$treatment) #crosslink=phenoData$crosslink,
    colsplit <- subphenoData$crosslink
    ngenes <- length(highlightGenes)
    nomenclature <- "UniprotName"
    pValstoHMP <- wholeDF[highlightGenes,c(nomenclature,pvalsColsSel)]
    toHMP <- wholeDF[highlightGenes,c(nomenclature,infoColumnsSel)]
    pValstoHMP <- pValstoHMP[!duplicated(pValstoHMP),]
    toHMP <- toHMP[!duplicated(toHMP),]
    
    pValstoHMP[,colnames(pValstoHMP)[grep("netc",colnames(pValstoHMP))]] <- 1
    #rownames(pValstoHMP) <- paste0(pValstoHMP$proteinName,1:nrow(pValstoHMP))
    rownames(pValstoHMP) <- rownames(toHMP) <- pValstoHMP[,nomenclature]
    pValstoHMP[,nomenclature] <- toHMP[,nomenclature] <- NULL
    
    #rownames(toHMP) <- paste0(toHMP$proteinName,1:nrow(toHMP))
    #rownames(toHMP) <- toHMP$proteinName
    #toHMP$proteinName <- NULL
    toHMP[is.na(toHMP)] <- 0
    
    ####
    # clustering
    ####
    logfile <- '/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/FinalData/HeatMaps/log.log'
    res.dist <- dist(toHMP,method = "euclidean")
    res.hc <- hclust(d = res.dist, method = "average")

    optClustMethods <- list()
    optNClusts <- c()
    indexes <- c("kl", "ch", "hartigan", "ccc", "scott", "marriot", "trcovw", "tracew", "friedman", "rubin", "cindex", "db", "silhouette", "duda", "pseudot2", "beale", "ratkowsky", "ball", "ptbiserial", "gap", "frey", "mcclain", "gamma", "gplus", "tau", "dunn", "hubert", "sdindex", "dindex", "sdbw")
    toremove <- c('gplus','gamma','tau')
    indexes <- setdiff(indexes,toremove)
    for (idx in indexes){
      print(idx)
      try(
        optClustMethods[[idx]] <- NbClust(data = toHMP, diss = NULL, 
                                          distance = "euclidean", 
                                          min.nc = 5, max.nc = 20,
                                          method = "average", alphaBeale = 0.1, 
                                          index = idx),
      )
      if (idx %in% names(optClustMethods)){
        if (!(idx %in% c("hubert", "dindex"))){
          resultbest <- optClustMethods[[idx]]$Best.nc[1]
          names(resultbest) <- idx
          optNClusts <- c(optNClusts,resultbest)
        }
      }
    }
    optNClusts
    optNClustsRes <- table(optNClusts)
    optNClustsRes <- optNClustsRes[2:(length(optNClustsRes)-1)]
    if (!any(optNClustsRes > 1)){
      optNClustsRes <- table(optNClusts)[1:(length(optNClustsRes)-1)]
    }
    optNClustsRes <- sort(optNClustsRes,decreasing = T)
    optNClust <- as.numeric(names(optNClustsRes)[1:3])
    
    # HEEEERE HERE
    
    
    #createHMP(toHMP,pValstoHMP,colsplit,colannot,outdir,tagname)
    col_fun <- colorRamp2(c(min(toHMP,na.rm = T),-1,0,1,max(toHMP,na.rm = T)), c("#332288","#88CCEE","white","#CC6677",'#661100'))
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
    
  
    
    tiff(file.path(outdir,paste0('topDEGs',ngenes,tagname,'pVal',pvalcutoff,'FC',fccutoff,'_HeatMap.tiff')), 1000, nrow(toHMP) * 11 + 200, pointsize=12,res = "print")
    print(draw(myHMPplot, heatmap_legend_side="left", annotation_legend_side = "left"))
    dev.off()
    
    pdf(file.path(outdir,paste0('topDEGs',ngenes,tagname,'pVal',pvalcutoff,'FC',fccutoff,'_HeatMap.pdf')), height = nrow(toHMP) / 10 + 10)
    print(draw(myHMPplot, heatmap_legend_side="left", annotation_legend_side = "left"))
    dev.off()
  }
}


