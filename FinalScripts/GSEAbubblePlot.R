## Normalize Legend for NES and Rel .Enr

library(dplyr)
library(reshape2)
library(ggplot2)

#treatments <- c("SA","DTT","H2O2")
#crosslinks <- c("FAX","UVX")

GSEAfilesDir <- "/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/FinalData/EnrichmentResults/GSEA"
treatment <- "SA"; crosslink <- "UVX"

GSEAfiles <- list.files(GSEAfilesDir,pattern = paste0(treatment,".*(degs|",crosslink,").*.tsv"),full.names = T)
allGSEAresults <- c()
for (GSEAfile in GSEAfiles){
  if (grepl("degs",GSEAfile)){
    dataType <- "Transc."
  }else if (grepl("RBPome",GSEAfile)){
    dataType <- "RBP"
  }else if (grepl("Prot",GSEAfile)){
    dataType <- "Prot."
  }else if (grepl("netc",GSEAfile)){
    dataType <- "Net C."
  }
  else{
    dataType <- "Unknown"
  }
  GSEAdf <- read.delim(GSEAfile)
  GSEAdf$datatype <- dataType
  allGSEAresults <- rbind(allGSEAresults,GSEAdf)
}

allGSEAresults <- allGSEAresults[allGSEAresults$padj < 0.05,]
allGSEAresults <- allGSEAresults[!(allGSEAresults$db == "panther" | is.na(allGSEAresults$db)), ]

table(allGSEAresults$datatype); length(unique(allGSEAresults$annotation_id))

masterEAdf <- type.convert(as.data.frame(allGSEAresults),as.is=T)
masterEAdf$db <- toupper(gsub("_"," ",masterEAdf$db))
#masterEAdf <- masterEAdf[order(masterEAdf$padj,decreasing = F),]
masterEAdf$padj <- -log10(masterEAdf$padj)
masterEAdf <- masterEAdf[order(masterEAdf$padj,decreasing = T),]

# Manually determine the order of the data types
datatypeOrder <- c("Transc.","Prot.","RBP","Net C.") 
masterEAdf$datatype <- factor(masterEAdf$datatype,levels = datatypeOrder)
masterEAdf$pathway <- factor(masterEAdf$pathway,levels = rev(unique(masterEAdf$pathway)))

# Change Uniprots for ORFs or GeneNames
mappingDF <- read.delim("/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/yeastReference/mappingFile.tsv",quote = "")
setOforfs <- c()
#nleadingGenes <- c()
for (setOfuniProts in masterEAdf$leadingEdge){
  uniprots <- unique(unlist(str_split(setOfuniProts,", ")))
  #nleadingGenes <- c(nleadingGenes, length(uniprots))
  orfs <- unique(mappingDF$ORF[mappingDF$UniprotACC %in% uniprots])
  setOforfs <- c(setOforfs,paste(sort(orfs),collapse = ", "))
  #geneNames <- unique(mappingDF$GeneName[mappingDF$UniprotACC %in% uniprots])
  if (!identical(length(unique(uniprots)),length(orfs))){ #,length(geneNames)
    print("WERID")
    print(c(length(unique(uniprots)),length(orfs)))#,length(geneNames)))
  }
}
masterEAdf$orfs <- setOforfs
masterEAdf$absNES <- abs(masterEAdf$NES)
#masterEAdf$nleadingGenes <- nleadingGenes

################################################################################
########################## 2. BUBBLE PLOT BY ANNOTTION DB ######################
################################################################################
# SIZE = pval_adj & color = relative_enrichment
outdir <- dirname(gsub("GSEA/","Visualizations/GSEA/",GSEAfiles[1]))
outfile <- file.path(outdir,paste(treatment,crosslink,sep = "_"))

# regulation <- c(1:nrow(masterEAdf))
# regulation[masterEAdf$NES > 0] <- "Up Regulated"
# regulation[masterEAdf$NES < 0] <- "Down Regulated"
# masterEAdf$regulation <- regulation
#absNESscale <- round(c(min(masterEAdf$absNES),max(masterEAdf$absNES)/2,max(masterEAdf$absNES)),2)
#scale_size_manual(values = absNESscale, breaks = absNESscale, labels=absNESscale), name = '|NES|') +
#scale_size(name = '|NES|') +
  
myplot <- ggplot(masterEAdf, aes(x = datatype, y = pathway, alpha = padj, size = absNES, colour = NES)) + 
  geom_point() +
  scale_size(name = '|NES|') +
  scale_alpha(range = c(0.5,1), name = '-log10(FDR)') +
  scale_colour_gradientn(colours = c('dodgerblue3','dodgerblue3','gold3','firebrick3','firebrick3'), name = 'NES') +
  ylab('') +
  xlab('Data Type') +
  facet_grid(rows=vars(db), scales = "free_y", space = "free_y") +
  theme(axis.text.y = element_text(size = 8, face = 'bold'),
        strip.text = element_text(hjust = 0.5, size = 12, face = 'bold'))

ggsave(paste0(outfile,".pdf"), plot = myplot, device = 'pdf', scale = 1, dpi = "print", 
       width = 4000, height = length(unique(masterEAdf$pathway)) * 40, units = 'px')

################################################################################
################# 3. Cluster Annotations by Genes sharing ######################
################################################################################
# https://www.datanovia.com/en/lessons/agglomerative-hierarchical-clustering/
# https://stats.stackexchange.com/questions/3685/where-to-cut-a-dendrogram

myNames <- unique(unlist(strsplit(masterEAdf$leadingEdge, split=", ")))
genesSharing <- setNames(data.frame(lapply(myNames, function(i) as.integer(grepl(i, masterEAdf$leadingEdge)))), myNames)
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
  
  myplot <- ggplot(masterEAdf2plot_cluster, aes(x = datatype, y = pathway, alpha = padj, size = absNES, colour = NES)) + 
    geom_point() +
    scale_size(name = '|NES|') +
    scale_alpha(range = c(0.5,1), name = '-log10(FDR)') +
    scale_colour_gradientn(colours = c('dodgerblue3','dodgerblue3','gold3','firebrick3','firebrick3'), name = 'NES') +
    ylab('') +
    xlab('Data Type') +
    facet_grid(rows=vars(cluster), scales = "free_y", space = "free_y") +
    theme(axis.text.y = element_text(size = 8, face = 'bold'),
          strip.text = element_text(hjust = 0.5, size = 12, face = 'bold'))
  
  ggsave(paste0(outfile,k,"Clusters.pdf"), plot = myplot, device = 'pdf', scale = 1, dpi = "print", 
         width = 4000, height = length(unique(masterEAdf2plot_cluster$pathway)) * 40, units = 'px')
  
  masterEAdf2plot_clusterGenes <- masterEAdf2plot_cluster
  masterEAdf2plot_clusterGenes$clusterGenes <- 1:nrow(masterEAdf2plot_clusterGenes)
  for (k in unique(masterEAdf2plot_cluster$cluster)){
    genesOfCluster <- str_split(paste(masterEAdf2plot_cluster$leadingEdge[masterEAdf2plot_cluster$cluster == k],collapse = ", "), ", ", simplify = T)
    masterEAdf2plot_clusterGenes$clusterGenes[masterEAdf2plot_cluster$cluster == k] <- paste(unique(as.vector(genesOfCluster)), collapse = ", ")
  }
  
  masterEAdf2plot_clusterGenes$clusterGenes <- swr(masterEAdf2plot_clusterGenes$clusterGenes,50)
  masterEAdf2plot_clusterGenes$datatype <- factor(masterEAdf2plot_clusterGenes$datatype,levels = datatypeOrder)
    
  myplot <- ggplot(masterEAdf2plot_clusterGenes, aes(x = datatype, y = pathway, alpha = padj, size = absNES, colour = NES)) + 
    geom_point() +
    scale_size(name = '|NES|') +
    scale_alpha(range = c(0.5,1), name = '-log10(FDR)') +
    scale_colour_gradientn(colours = c('dodgerblue3','dodgerblue3','gold3','firebrick3','firebrick3'), name = 'NES') +
    ylab('') +
    xlab('Data Type') +
    facet_grid(rows=vars(clusterGenes), scales = "free_y", space = "free_y") +
    theme(axis.text.y = element_text(size = 8, face = 'bold'),
          axis.text.x = element_text(size = 8, face = 'bold'),
          strip.text = element_text(hjust = 0.5, size = 8, face = 'bold'),
          strip.text.y.right = element_text(angle = 0, size = 4))
  
  ggsave(paste0(outfile,k,"Clusters_EnrGenes",".pdf"), plot = myplot, device = 'pdf', scale = 1, dpi = "print", 
         width = 4000, height = length(unique(masterEAdf2plot_clusterGenes$pathway)) * 40, units = 'px')
}




