library(dplyr)
library(reshape2)
library(ggplot2)
library(stringr)
range01 <- function(x){(x-min(x))/(max(x)-min(x))}
swr = function(string, nwrap=50) {paste(strwrap(string, width=nwrap), collapse="\n")}; swr = Vectorize(swr)
source("FinalScripts/functionOmics.R")
################################################################################
####################### 1. SELECT GC4 ENRICHMENT FILES #########################
################################################################################
#treatments <- c("SA","DTT","H2O2")
#crosslinks <- c("FAX","UVX")

ORAfilesDir <- "/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/FinalData/EnrichmentResults/ORAseparatedUpnDw"

treatment <- "SA"; crosslink <- "UVX"
ORAfiles <- list.files(ORAfilesDir,pattern = paste0(treatment,".*(DEGs|",crosslink,").*\\.Enr.tsv"),full.names = T)
getORAplots(ORAfiles)

treatment <- "SA"; crosslink <- "FAX"
ORAfiles <- list.files(ORAfilesDir,pattern = paste0(treatment,".*(DEGs|",crosslink,").*\\.Enr.tsv"),full.names = T)
getORAplots(ORAfiles)

ORAfilesDir <- "/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/FinalData/EnrichmentResults/ORA"

treatment <- "SA"; crosslink <- "UVX"
ORAfiles <- list.files(ORAfilesDir,pattern = paste0(treatment,".*(DEGs|",crosslink,").*\\.Enr.tsv"),full.names = T)
getORAplots(ORAfiles)

treatment <- "SA"; crosslink <- "FAX"
ORAfiles <- list.files(ORAfilesDir,pattern = paste0(treatment,".*(DEGs|",crosslink,").*\\.Enr.tsv"),full.names = T)
getORAplots(ORAfiles)


getORAplots <- function(ORAfiles){
  allORAresults <- c()
  for (ORAfile in ORAfiles){
    ORAdf <- read.delim(ORAfile)
    if (grepl("DEGs",ORAfile)){
      dataType <- "Transc."
    }else if (grepl("rbpm",ORAfile)){
      dataType <- "RBP"
    }else if (grepl("prot",ORAfile)){
      dataType <- "Prot."
    }else if (grepl("netc",ORAfile)){
      dataType <- "Net C."
    }
    else{
      dataType <- "Unknown"
    }
    if (grepl("up",basename(ORAfile),ignore.case = T)){
      # Manually determine the order of the data types
      datatypeOrder <- c("Transc.\nUp","Transc.\nDown","Prot.\nUp","Prot.\nDown","RBP\nUp","RBP\nDown","Net C.\nUp","Net C.\nDown") 
      regul <- "\nUp"
    }else if(grepl("dw",basename(ORAfile),ignore.case = T)){
      datatypeOrder <- c("Transc.\nUp","Transc.\nDown","Prot.\nUp","Prot.\nDown","RBP\nUp","RBP\nDown","Net C.\nUp","Net C.\nDown") 
      regul <- "\nDown"
    }else{
      datatypeOrder <- c("Transc.","Prot.","RBP","Net C.") 
      regul <- ""
    }
    ORAdf$datatype <- paste0(dataType,regul)
    allORAresults <- rbind(allORAresults,ORAdf)
  }
  allORAresults <- allORAresults[allORAresults$pval_adj < 0.05,]
  table(allORAresults$datatype); length(unique(allORAresults$annotation_id))
  
  masterEAdf <- type.convert(as.data.frame(allORAresults),as.is=T)
  masterEAdf$annotation <- gsub("_"," ",masterEAdf$annotation)
  masterEAdf$pval_adj <- -log10(masterEAdf$pval_adj)
  masterEAdf <- masterEAdf[order(masterEAdf$pval_adj,decreasing = T),]
  masterEAdf$fulldesc <- paste(masterEAdf$description,masterEAdf$annotation_id)
  
  masterEAdf$datatype <- factor(masterEAdf$datatype,levels = datatypeOrder)
  masterEAdf$fulldesc <- factor(masterEAdf$fulldesc,levels = rev(unique(masterEAdf$fulldesc)))
  
  # Change Uniprots for ORFs or GeneNames
  mappingDF <- read.delim("/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/yeastReference/mappingFile.tsv",quote = "")
  setOforfs <- c()
  for (setOfuniProts in masterEAdf$genes){
    uniprots <- unlist(str_split(setOfuniProts,", "))
    orfs <- unique(mappingDF$ORF[mappingDF$UniprotACC %in% uniprots])
    setOforfs <- c(setOforfs,paste(sort(orfs),collapse = ", "))
    #geneNames <- unique(mappingDF$GeneName[mappingDF$UniprotACC %in% uniprots])
    if (!identical(length(unique(uniprots)),length(orfs))){ #,length(geneNames)
      print("WERID")
      print(c(length(unique(uniprots)),length(orfs)))#,length(geneNames)))
    }
  }
  
  masterEAdf$orfs <- setOforfs
  ################################################################################
  ########################## 2. BUBBLE PLOT BY ANNOTTION DB ######################
  ################################################################################
  # SIZE = pval_adj & color = relative_enrichment
  outdir <- dirname(gsub("ORA","Visualizations/ORA",ORAfiles[1]))
  outfile <- file.path(outdir,paste(treatment,crosslink,sep = "_"))
  
  myplot <- ggplot(masterEAdf, aes(x = datatype, y = fulldesc, size = relative_enrichment, color = pval_adj)) + 
    geom_point(alpha = 0.7) +
    scale_size(range = c(3, 8), name = 'Rel. Enr.') +
    scale_color_gradientn(colours = c('maroon3', 'gold', 'turquoise3'), name = '-log10(FDR)') +
    ylab('') +
    xlab('Data Type') +
    facet_grid(rows=vars(annotation), scales = "free_y", space = "free_y") +
    theme(axis.text.y = element_text(size = 10, face = 'bold'),
          strip.text = element_text(hjust = 0.5, size = 12, face = 'bold'))
  
  outdivice <- "tiff"
  ggsave(paste0(outfile,".",outdivice), plot = myplot, device = outdivice, scale = 1, dpi = "print", 
         width = 4000, height = length(unique(masterEAdf$fulldesc)) * 40, units = 'px')

  ################################################################################
  ################# 3. Cluster Annotations by Genes sharing ######################
  ################################################################################
  # https://www.datanovia.com/en/lessons/agglomerative-hierarchical-clustering/
  # https://stats.stackexchange.com/questions/3685/where-to-cut-a-dendrogram
  myNames <- unique(unlist(strsplit(masterEAdf$genes, split=", ")))
  genesSharing <- setNames(data.frame(lapply(myNames, function(i) as.integer(grepl(i, masterEAdf$genes)))), myNames)
  res.dist <- dist(genesSharing, method = "binary")
  
  oldbestclustcorrScore <- 0
  logfile <- paste0(outfile,".log")
  write("Best hclustMethod:",logfile,append = F)
  hclustMethods <- c("ward.D","ward.D2","single","complete","average","mcquitty","median","centroid")
  for(hclustMethod in hclustMethods){
    res.hc <- hclust(d = res.dist, method = hclustMethod)
    res.coph <- cophenetic(res.hc)
    clustcorrScore <- cor(res.dist, res.coph)
    write(paste("hclustMethod:",hclustMethod,"clustcorrScore:", clustcorrScore),logfile,append = T)
    cat(clustcorrScore,oldbestclustcorrScore)
    if (clustcorrScore > oldbestclustcorrScore){
      bestHclustMethod <- hclustMethod
      oldbestclustcorrScore <- clustcorrScore
      bestRes.hc <- res.hc
    }
  }
  write(paste("Best Hclust Method:",bestHclustMethod),logfile,append = T)
  res.hc <- hclust(d = res.dist, method = bestHclustMethod)
  
  masterEAdf2plot_cluster <- masterEAdf
  masterEAdf2plot_cluster$datatype <- factor(masterEAdf2plot_cluster$datatype,levels = rev(unique(masterEAdf2plot_cluster$datatype)))
  
  for (k in 10:20){
    grp <- cutree(bestRes.hc, k = k)
    masterEAdf2plot_cluster$cluster <- grp
    masterEAdf2plot_cluster$datatype <- factor(masterEAdf2plot_cluster$datatype,levels = datatypeOrder)
    
    myplot <- ggplot(masterEAdf2plot_cluster, aes(x = datatype, y = fulldesc, size = relative_enrichment, color = pval_adj)) + 
      geom_point(alpha = 0.7) +
      scale_size(range = c(3, 8), name = 'Rel. Enr.') +
      scale_color_gradientn(colours = c('maroon3', 'gold', 'turquoise3'), name = '-log10(FDR)') +
      ylab('') +
      xlab('Data Type') +
      facet_grid(rows=vars(cluster), scales = "free_y", space = "free_y") +
      theme(axis.text.y = element_text(size = 10, face = 'bold'),
            strip.text = element_text(hjust = 0.5, size = 12, face = 'bold'))
    
    ggsave(paste0(outfile,k,"Clusters.",outdivice), plot = myplot, device = outdivice, scale = 1, dpi = "print", 
           width = 4000, height = length(unique(masterEAdf2plot_cluster$fulldesc)) * 40, units = 'px')
    
    masterEAdf2plot_clusterGenes <- masterEAdf2plot_cluster
    masterEAdf2plot_clusterGenes$clusterGenes <- 1:nrow(masterEAdf2plot_clusterGenes)
    for (k in unique(masterEAdf2plot_cluster$cluster)){
      genesOfCluster <- str_split(paste(masterEAdf2plot_cluster$orfs[masterEAdf2plot_cluster$cluster == k],collapse = ", "), ", ", simplify = T)
      masterEAdf2plot_clusterGenes$clusterGenes[masterEAdf2plot_cluster$cluster == k] <- paste(unique(as.vector(genesOfCluster)), collapse = ", ")
    }
    
    masterEAdf2plot_clusterGenes$clusterGenes <- swr(masterEAdf2plot_clusterGenes$clusterGenes,50)
    masterEAdf2plot_clusterGenes$datatype <- factor(masterEAdf2plot_clusterGenes$datatype,levels = datatypeOrder)
    
    myplot <- ggplot(masterEAdf2plot_clusterGenes, aes(x = datatype, y = fulldesc, size = relative_enrichment, color = pval_adj)) + 
      geom_point(alpha = 0.7) +
      scale_size(range = c(3, 8), name = 'Rel. Enr.') +
      scale_color_gradientn(colours = c('maroon3', 'gold', 'turquoise3'), name = '-log10(FDR)') +
      ylab('') +
      xlab('Data Type') +
      facet_grid(rows=vars(clusterGenes), scales = "free_y", space = "free_y") +
      theme(axis.text.y = element_text(size = 8, face = 'bold'),
            axis.text.x = element_text(size = 8, face = 'bold'),
            strip.text = element_text(hjust = 0.5, size = 8, face = 'bold'),
            strip.text.y.right = element_text(angle = 0, size = 4))
    
    ggsave(paste0(outfile,k,"Clusters_EnrGenes.",outdivice), plot = myplot, device = outdivice, scale = 1, dpi = "print", 
           width = 4000, height = length(unique(masterEAdf2plot_clusterGenes$fulldesc)) * 40, units = 'px')
  }
} 
  
