library(NbClust)
library(dplyr)
library(reshape2)
library(ggplot2)
library(stringr)
range01 <- function(x){(x-min(x))/(max(x)-min(x))}
swr = function(string, nwrap=50) {paste(strwrap(string, width=nwrap), collapse="\n")}; swr = Vectorize(swr)
source("FinalScripts/functionOmics.R")

FAAgenes <- list(YOR317W="FAA1",YER015W="FAA2",YIL009W="FAA3",YMR246W="FAA4")
FAAunipr <- list(P30624="FAA1",P39518="FAA2",P39002="FAA3",P47912="FAA4")

getFAAs <- function(genes){
  FAA <- FAAunipr[which(names(FAAunipr) %in% unlist(strsplit(genes,split = ', ')))]
  if (length(FAA) > 0)
    return(gsub('FAA','*',paste0(FAA,collapse = ",")))
  else
    return('')
}

getORAplots <- function(ORAfiles){
  allORAresults <- c()
  for (ORAfile in ORAfiles){
    ORAdf <- read.delim(ORAfile)
    if (grepl("DEGs",ORAfile)){
      dataType <- "Transc."
    }else if (grepl("rbpm",ORAfile)){
      dataType <- "mRBP"
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
      datatypeOrder <- c("Transc.\nUp","Transc.\nDown","Prot.\nUp","Prot.\nDown","mRBP\nUp","mRBP\nDown","Net C.\nUp","Net C.\nDown") 
      regul <- "\nUp"
    }else if(grepl("dw",basename(ORAfile),ignore.case = T)){
      datatypeOrder <- c("Transc.\nUp","Transc.\nDown","Prot.\nUp","Prot.\nDown","mRBP\nUp","mRBP\nDown","Net C.\nUp","Net C.\nDown") 
      regul <- "\nDown"
    }else{
      datatypeOrder <- c("Transc.","Prot.","mRBP","Net C.") 
      regul <- ""
    }
    ORAdf$datatype <- paste0(dataType,regul)
    allORAresults <- rbind(allORAresults,ORAdf)
  }
  #allORAresults <- allORAresults[allORAresults$pval_adj < 0.05,]
  interestingAnnots <- unique(allORAresults$annotation_id[which(allORAresults$pval_adj < 0.05 & allORAresults$datatype %in% c("mRBP","Net C."))])
  allORAresults <- allORAresults[which(allORAresults$annotation_id %in% interestingAnnots),]
  
  table(allORAresults$datatype); length(unique(allORAresults$annotation_id))
  
  masterEAdf <- type.convert(as.data.frame(allORAresults),as.is=T)
  masterEAdf$annotation <- gsub("_"," ",masterEAdf$annotation)
  masterEAdf$pval_adj <- -log10(masterEAdf$pval_adj)
  masterEAdf <- masterEAdf[order(masterEAdf$pval_adj,decreasing = T),]
  FAASloc <- unlist(lapply(masterEAdf$genes, FUN=getFAAs))
  masterEAdf$fulldesc <- paste(masterEAdf$description,masterEAdf$annotation_id,FAASloc)
  
  masterEAdf$datatype <- factor(masterEAdf$datatype,levels = datatypeOrder)
  masterEAdf$fulldesc <- factor(masterEAdf$fulldesc,levels = rev(unique(masterEAdf$fulldesc)))
  
  # Change Uniprots for ORFs or GeneNames
  mappingDF <- read.delim("yeastReference/mappingFile.tsv",quote = "")
  mappingDF$UniprotName <- gsub("_YEAST","",mappingDF$UniprotName)
  setOforfs <- c()
  for (setOfuniProts in masterEAdf$genes){
    uniprots <- unlist(str_split(setOfuniProts,", "))
    orfs <- unique(mappingDF$UniprotName[mappingDF$UniprotACC %in% uniprots])
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
  # SIZE = pval_adj, color = relative_enrichment
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
         width = 4000, height = length(unique(masterEAdf$fulldesc)) * 50, units = 'px', limitsize = FALSE)
  
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
    if (clustcorrScore > oldbestclustcorrScore){
      bestHclustMethod <- hclustMethod
      oldbestclustcorrScore <- clustcorrScore
      bestRes.hc <- res.hc
    }
  }
  write(paste("Best Hclust Method:",bestHclustMethod),logfile,append = T)
  res.hc <- hclust(d = res.dist, method = bestHclustMethod)
  
  optClustMethods <- list()
  optNClusts <- c()
  indexes <- c("kl", "ch", "hartigan", "ccc", "scott", "marriot", "trcovw", "tracew", "friedman", "rubin", "cindex", "db", "silhouette", "duda", "pseudot2", "beale", "ratkowsky", "ball", "ptbiserial", "gap", "frey", "mcclain", "gamma", "gplus", "tau", "dunn", "hubert", "sdindex", "dindex", "sdbw")
  #indexes <- c('kl','ch','hartigan','cindex','db','silhouette','duda','pseudot2','beale','ratkowsky','ball','ptbiserial','gap','frey','mcclain','tau','dunn','hubert','sdindex','dindex','sdbw')
  toremove <- c('gplus','gamma','tau')
  indexes <- setdiff(indexes,toremove)
  #c("hartigan","ball")
  #indexes <- setdiff(indexes,names(optClustMethods[[idx]]))
  for (idx in indexes){
    print(idx)
    try(
      optClustMethods[[idx]] <- NbClust(data = genesSharing, diss = NULL, 
                                        distance = "binary", 
                                        min.nc = 5,#round(nrow(genesSharing) * 0.05),
                                        max.nc = 20,#round(nrow(genesSharing) * 0.2), 
                                        method = bestHclustMethod, alphaBeale = 0.1, 
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
  
  masterEAdf2plot_cluster <- masterEAdf
  masterEAdf2plot_cluster$datatype <- factor(masterEAdf2plot_cluster$datatype,levels = rev(unique(masterEAdf2plot_cluster$datatype)))
  print(optNClust)
  for (k in optNClust){
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
           width = 4000, height = length(unique(masterEAdf2plot_cluster$fulldesc)) * 55, units = 'px', limitsize = FALSE)
    
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
           width = 4000, height = length(unique(masterEAdf2plot_clusterGenes$fulldesc)) * 55, units = 'px', limitsize = FALSE)
  }
} 

################################################################################
####################### 1. SELECT GC4 ENRICHMENT FILES #########################
################################################################################
#treatments <- c("SA","DTT","H2O2")
#crosslinks <- c("FAX","UVX")

# ORAfilesDir <- "FinalData/EnrichmentResults/ORAseparatedUpnDw"
# treatment <- "SA"; crosslink <- "UVX"
# ORAfiles <- list.files(ORAfilesDir,pattern = paste0(treatment,".*(DEGs|",crosslink,").*\\.Enr.tsv"),full.names = T)
# getORAplots(ORAfiles)
# 
# treatment <- "SA"; crosslink <- "FAX"
# ORAfiles <- list.files(ORAfilesDir,pattern = paste0(treatment,".*(DEGs|",crosslink,").*\\.Enr.tsv"),full.names = T)
# getORAplots(ORAfiles)

ORAfilesDir <- "FinalData/EnrichmentResults/ORA"
treatment <- "SA"; crosslink <- "UVX"
ORAfiles <- list.files(ORAfilesDir,pattern = paste0(treatment,".*(DEGs|",crosslink,").*\\.Enr.tsv"),full.names = T)
getORAplots(ORAfiles)

treatment <- "SA"; crosslink <- "FAX"
ORAfiles <- list.files(ORAfilesDir,pattern = paste0(treatment,".*(DEGs|",crosslink,").*\\.Enr.tsv"),full.names = T)
getORAplots(ORAfiles)
