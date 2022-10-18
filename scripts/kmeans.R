outdir <- "/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/data/clustering"

library(openxlsx)
UVSA <- read.xlsx("/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/data/clustering/UVX and FAX mRBPome Adrian.xlsx",sheet = 1)

FAXSA <- read.xlsx("/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/data/clustering/UVX and FAX mRBPome Adrian.xlsx",sheet = 2)
DF <- FAXSA
DF$GeneName <- NULL
DF <- DF[,c("mRNABPome","netChange","proteome","DEGs")]
DF <- DF[complete.cases(DF),]
DF <- DF[!duplicated(DF),]
rownames(DF) <- paste(DF$UniprotACC,row.names(DF),sep = "_")
DF$UniprotACC <- NULL

k.max <- 10
wss <- sapply(1:k.max, function(k){kmeans(DF, k, nstart=50,iter.max = 15 )$tot.withinss})
plot(1:k.max, wss,
     type="b", pch = 19, frame = FALSE, 
     xlab="Number of clusters K",
     ylab="Total within-clusters sum of squares")

DFcluster <- kmeans(DF,centers = 6)
DF$cluster <- DFcluster$cluster

library(ComplexHeatmap)
library(circlize)
toHMAP <- DF[order(DF$cluster),]
toHMAP$cluster <- NULL
maxcolor <- round(max(abs(c(min(toHMAP),max(toHMAP)))))
col_fun <- colorRamp2(c(-maxcolor,0,maxcolor), c("yellow","white","dodgerblue4"))
ha <- rowAnnotation(df=DF[,"cluster",drop=F],col = list("cluster"=colorRamp2(c(1:6), rep("black", 6))))

Heatmap(toHMAP, cluster_columns = F, cluster_rows = T,col=col_fun,show_row_names = F,show_row_dend = F,
        row_split = DF$cluster,row_gap = unit(5, "mm"),right_annotation = ha)


png(file.path(outdir,'HeatMap_log2FCsFAX.png'), 1200, 900, pointsize=12,res = 150)
dev.off()
