


source("scripts/functionsOmics.R")
setwd("/home/eidriangm/Desktop/toDo/surrey/multiregulatomics")
################
##### DATA #####
################


wholeDFfile <- "/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/data/mapped/allDataMapped_together.tsv"
wholeDF <- read.delim(wholeDFfile,quote = "")

phenoDatafile <- "/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/data/phenodata.tsv" 
phenoData <- read.delim(phenoDatafile,check.names = F)

protcols <- colnames(wholeDF)[grep("prot.*X[0-9]+",colnames(wholeDF))]; length(protcols)
phenoData$protdatnames <- protcols[match(phenoData$proteomeNumber,as.numeric(gsub("X0+","",str_extract(protcols,"X0[0-9]+"))))]
sum(!is.na(phenoData$protdatnames))

mrbpmcols <- colnames(wholeDF)[grep("mRNABPome.*X[0-9]+",colnames(wholeDF))]; length(mrbpmcols)
phenoData$mrbpomedatnames <- mrbpmcols[match(phenoData$mrbpomeNumber,as.numeric(gsub("X0+","",str_extract(mrbpmcols,"X0[0-9]+"))))]
sum(!is.na(phenoData$mrbpomedatnames))

phenoData <- phenoData[,c('mrbpomedatnames','mrbpomeNumber','mrbpomeCrosslinker','mrbpomeTreatment','mrbpomeBatch',
                          'protdatnames','proteomeNumber','proteomeCrosslinker','proteomeTreatment','proteomeBatch')]

outdir <- "/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/data" 
saveTablesTsvExc(phenoData,outdir)
phenoDatafile <- "/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/data/phenoData.tsv" 
phenoData <- read.delim(phenoDatafile,check.names = F)

outdir <- "/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/data/visualizations"

crosslinkers <- c("UV","FAX")
treatments <- c("no","H202","SA","DTT")
treatment <- treatments[1]

###############
#### PCAs #####
###############

#samples2pca <- phenoData$mrbpomedatnames[phenoData$mrbpomeCrosslinker == crosslinker & phenoData$mrbpomeTreatment == treatment]

#####################
#### mrbopme UV #####
#####################

crosslinker <- "UV"
samples2pca <- phenoData$mrbpomedatnames[phenoData$mrbpomeCrosslinker == crosslinker]
samples2pca <- na.omit(samples2pca)

subphenoData <- phenoData[phenoData$mrbpomedatnames %in% samples2pca,]
sampleClas <- na.omit(paste(subphenoData$mrbpomeCrosslinker,subphenoData$mrbpomeTreatment,subphenoData$mrbpomeBatch,sep="_"))

toPCA <- wholeDF[,samples2pca] # cuant values
toPCA <- toPCA[complete.cases(toPCA),]
length(unique(sampleClas))
pcares <- prcomp(t(toPCA))
sampleClas
factornum <- as.numeric(as.factor(sampleClas))
colors <- paletteer_d("ggsci::category20_d3")[factornum]

png(file.path(outdir,paste0('mrbpome',crosslinker,'_2D.png')), 1500, 1200, pointsize=10,res = 300)
ggplot2::autoplot(pcares,data = data.frame(t(toPCA),class=sampleClas),
                  colour="class", label=T, label.size=1,label.repel=T,main = paste0('mrbpome ',crosslinker)) + 
  scale_color_manual(values = unique(as.vector(colors)))+ theme_classic() 
dev.off()

pcaVarGroup <- round((pcares$sdev)^2 / sum(pcares$sdev^2) *100, 2)
png(file.path(outdir,paste0('mrbpome',crosslinker,'_3D.png')), 1200, 900, pointsize=12,res = 150)
pca3d <- scatterplot3d(pcares$x[,1], pcares$x[,3], pcares$x[,2],color = as.vector(colors),
                       main = paste0('mrbpome_',crosslinker), xlab = paste0("PC 1 (", pcaVarGroup[1], " %)"),
                       zlab = paste0("PC 2 (", pcaVarGroup[2], " %)"), ylab = paste0("PC 3 (", pcaVarGroup[3], " %)"),
                       grid=T, box = T, pch = 20, cex.symbols = 2.5, angle = 40,type = "h")
zz.coords <- pca3d$xyz.convert(pcares$x[,1], pcares$x[,3], pcares$x[,2])
legend("topright", pch=20, legend = unique(sampleClas), col = unique(as.vector(colors)),inset = 0,y.intersp =0.8)
text(zz.coords$x, zz.coords$y, labels = row.names(pcares$x),cex = 0.5, pos = 4,col = as.vector(colors))
dev.off()

####################
#### mrbopme FAX ###
####################

crosslinker <- "FAX"
samples2pca <- phenoData$mrbpomedatnames[phenoData$mrbpomeCrosslinker == crosslinker]
samples2pca <- na.omit(samples2pca)

subphenoData <- phenoData[phenoData$mrbpomedatnames %in% samples2pca,]
sampleClas <- na.omit(paste(subphenoData$mrbpomeCrosslinker,subphenoData$mrbpomeTreatment,subphenoData$mrbpomeBatch,sep="_"))

toPCA <- wholeDF[,samples2pca] # cuant values
toPCA <- toPCA[complete.cases(toPCA),]
pcares <- prcomp(t(toPCA))
factornum <- as.numeric(as.factor(sampleClas))
colors <- paletteer_d("ggsci::category20_d3")[factornum]

png(file.path(outdir,paste0('mrbpome',crosslinker,'_2D.png')), 1500, 1200, pointsize=10,res = 300)
ggplot2::autoplot(pcares,data = data.frame(t(toPCA),class=sampleClas),
                  colour="class", label=T, label.size=1,label.repel=T,main = paste0('mrbpome ',crosslinker)) + 
  scale_color_manual(values = unique(as.vector(colors)))+ theme_classic() 
dev.off()

pcaVarGroup <- round((pcares$sdev)^2 / sum(pcares$sdev^2) *100, 2)
png(file.path(outdir,paste0('mrbpome',crosslinker,'_3D.png')), 1200, 900, pointsize=12,res = 150)
pca3d <- scatterplot3d(pcares$x[,1], pcares$x[,3], pcares$x[,2],color = as.vector(colors),
                       main = paste0('mrbpome_',crosslinker), xlab = paste0("PC 1 (", pcaVarGroup[1], " %)"),
                       zlab = paste0("PC 2 (", pcaVarGroup[2], " %)"), ylab = paste0("PC 3 (", pcaVarGroup[3], " %)"),
                       grid=T, box = T, pch = 20, cex.symbols = 2.5, angle = 40,type = "h")
zz.coords <- pca3d$xyz.convert(pcares$x[,1], pcares$x[,3], pcares$x[,2])
legend("topright", pch=20, legend = unique(sampleClas), col = unique(as.vector(colors)),inset = 0,y.intersp =0.8)
text(zz.coords$x, zz.coords$y, labels = row.names(pcares$x),cex = 0.5, pos = 4,col = as.vector(colors))
dev.off()

#####################
#### proteome FAX ###
#####################

crosslinker <- "FAX"
samples2pca <- phenoData$protdatnames[phenoData$proteomeCrosslinker == crosslinker]
samples2pca <- na.omit(samples2pca)

subphenoData <- phenoData[phenoData$protdatnames %in% samples2pca,]
sampleClas <- na.omit(paste(subphenoData$proteomeCrosslinker,subphenoData$proteomeTreatment,subphenoData$proteomeBatch,sep="_"))

toPCA <- wholeDF[,samples2pca] # cuant values
toPCA <- toPCA[complete.cases(toPCA),]
pcares <- prcomp(t(toPCA))

factornum <- as.numeric(as.factor(sampleClas))
colors <- paletteer_d("ggsci::category20_d3")[factornum]

png(file.path(outdir,paste0('proteome',crosslinker,'_2D.png')), 1500, 1200, pointsize=10,res = 300)
ggplot2::autoplot(pcares,data = data.frame(t(toPCA),class=sampleClas),
                  colour="class", label=T, label.size=0.999,label.repel=T,main = paste0('proteome ',crosslinker)) + 
  scale_color_manual(values = unique(as.vector(colors)))+ theme_classic() 
dev.off()

pcaVarGroup <- round((pcares$sdev)^2 / sum(pcares$sdev^2) *100, 2)
png(file.path(outdir,paste0('proteome',crosslinker,'_3D.png')), 1200, 900, pointsize=12,res = 150)
pca3d <- scatterplot3d(pcares$x[,1], pcares$x[,3], pcares$x[,2],color = as.vector(colors),
                       main = paste0('proteome_',crosslinker), xlab = paste0("PC 1 (", pcaVarGroup[1], " %)"),
                       zlab = paste0("PC 2 (", pcaVarGroup[2], " %)"), ylab = paste0("PC 3 (", pcaVarGroup[3], " %)"),
                       grid=T, box = T, pch = 20, cex.symbols = 2.5, angle = 40,type = "h")
zz.coords <- pca3d$xyz.convert(pcares$x[,1], pcares$x[,3], pcares$x[,2])
legend("topright", pch=20, legend = unique(sampleClas), col = unique(as.vector(colors)),inset = 0,y.intersp =0.8)
text(zz.coords$x, zz.coords$y, labels = row.names(pcares$x),cex = 0.5, pos = 4,col = as.vector(colors))
dev.off()

#####################
#### proteome UV ###
#####################

crosslinker <- "UV"
samples2pca <- phenoData$protdatnames[phenoData$proteomeCrosslinker == crosslinker]
samples2pca <- na.omit(samples2pca)

subphenoData <- phenoData[phenoData$protdatnames %in% samples2pca,]
sampleClas <- na.omit(paste(subphenoData$proteomeCrosslinker,subphenoData$proteomeTreatment,subphenoData$proteomeBatch,sep="_"))

toPCA <- wholeDF[,samples2pca] # cuant values
toPCA <- toPCA[complete.cases(toPCA),]
pcares <- prcomp(t(toPCA))

factornum <- as.numeric(as.factor(sampleClas))
colors <- paletteer_d("ggsci::category20_d3")[factornum]

png(file.path(outdir,paste0('proteome',crosslinker,'_2D.png')), 1500, 1200, pointsize=10,res = 300)
ggplot2::autoplot(pcares,data = data.frame(t(toPCA),class=sampleClas),
                  colour="class", label=T, label.size=0.999,label.repel=T,main = paste0('proteome ',crosslinker)) + 
  scale_color_manual(values = unique(as.vector(colors)))+ theme_classic() 
dev.off()

pcaVarGroup <- round((pcares$sdev)^2 / sum(pcares$sdev^2) *100, 2)
png(file.path(outdir,paste0('proteome',crosslinker,'_3D.png')), 1200, 900, pointsize=12,res = 150)
pca3d <- scatterplot3d(pcares$x[,1], pcares$x[,3], pcares$x[,2],color = as.vector(colors),
                       main = paste0('proteome_',crosslinker), xlab = paste0("PC 1 (", pcaVarGroup[1], " %)"),
                       zlab = paste0("PC 2 (", pcaVarGroup[2], " %)"), ylab = paste0("PC 3 (", pcaVarGroup[3], " %)"),
                       grid=T, box = T, pch = 20, cex.symbols = 2.5, angle = 40,type = "h")
zz.coords <- pca3d$xyz.convert(pcares$x[,1], pcares$x[,3], pcares$x[,2])
legend("topright", pch=20, legend = unique(sampleClas), col = unique(as.vector(colors)),inset = 0,y.intersp =0.8)
text(zz.coords$x, zz.coords$y, labels = row.names(pcares$x),cex = 0.5, pos = 4,col = as.vector(colors))
dev.off()

############################
#### CORRELATION PLOTS #####
############################

#### log2FCs FAX #####
log2FCsFAX <- cbind(wholeDF$proteomeFAX.log2ratio_Condition2,
                    wholeDF$proteomeFAX.log2ratio_Condition3,
                    wholeDF$proteomeFAX.log2ratio_Condition4,
                    wholeDF$mRNABPomeFAX.log2ratio_PolyARNAFAXwithDTT,
                    wholeDF$mRNABPomeFAX.log2ratio_PolyARNAFAXwithH2O2,
                    wholeDF$mRNABPomeFAX.log2ratio_PolyARNAFAXwithSA,
                    wholeDF$mRNABPomeFAX.log2ratio_PolyARNAFAXwithDTT - wholeDF$proteomeFAX.log2ratio_Condition2,
                    wholeDF$mRNABPomeFAX.log2ratio_PolyARNAFAXwithH2O2 - wholeDF$proteomeFAX.log2ratio_Condition3,
                    wholeDF$mRNABPomeFAX.log2ratio_PolyARNAFAXwithSA - wholeDF$proteomeFAX.log2ratio_Condition4,
                    wholeDF$DEGs.log2FC.DTT.YPD, wholeDF$DEGs.log2FC.H202.YPD, wholeDF$DEGs.log2FC.SA.YPD)

corrcolnames <- c("prot.DTT",
                  "prot.H2O2", 
                  "prot.SA",
                  "RBP.DTT", 
                  "RBP.H2O2", 
                  "RBP.SA", 
                  "NetCh.DTT", 
                  "NetCh.H2O2",
                  "NetCh.SA", 
                  "DEGs.DTT","DEGs.H2O2","DEGs.SA")

colnames(log2FCsFAX) <- corrcolnames
log2FCsFAX <- as.data.frame(log2FCsFAX)
log2FCsFAX <- log2FCsFAX[complete.cases(log2FCsFAX),]
log2FCsFAX <- log2FCsFAX[!duplicated(log2FCsFAX),]
row.names(log2FCsFAX) <- log2FCsFAX$ORF

corMat <- cor(log2FCsFAX,method = "pearson")
png(file.path(outdir,paste0('log2FCnNetCHsFAX.png')), 1200, 900, pointsize=12,res = 150)
corrplot(corMat, method = 'circle', order = 'hclust',tl.cex=0.5,title="log2FCnNetCHsFAX",
         mar = c(0, 1, 2, 1),tl.pos = 'lt',addCoef.col = "grey",number.cex = 0.5,
         col=rev(COL2('RdBu', 200)))
dev.off()


#### log2FCs UV #####
log2FCsUV <- cbind(wholeDF$proteomeUV.log2ratio_Condition2,
                   wholeDF$proteomeUV.log2ratio_Condition3,
                   wholeDF$proteomeUV.log2ratio_Condition4,
                   wholeDF$mRNABPomeUV.log2ratio_PolyARNAUVwithDTT,
                   wholeDF$mRNABPomeUV.log2ratio_PolyARNAUVwithH202,
                   wholeDF$mRNABPomeUV.log2ratio_PolyARNAUVwithSA,
                   wholeDF$mRNABPomeUV.log2ratio_PolyARNAUVwithDTT - wholeDF$proteomeUV.log2ratio_Condition2,
                   wholeDF$mRNABPomeUV.log2ratio_PolyARNAUVwithH202 - wholeDF$proteomeUV.log2ratio_Condition3,
                   wholeDF$mRNABPomeUV.log2ratio_PolyARNAUVwithSA - wholeDF$proteomeUV.log2ratio_Condition4,
                   wholeDF$DEGs.log2FC.DTT.YPD, wholeDF$DEGs.log2FC.H202.YPD, wholeDF$DEGs.log2FC.SA.YPD)

corrcolnames <- c("prot.DTT",
                  "prot.H2O2", 
                  "prot.SA",
                  "RBP.DTT", 
                  "RBP.H2O2", 
                  "RBP.SA", 
                  "NetCh.DTT", 
                  "NetCh.H2O2",
                  "NetCh.SA", 
                  "DEGs.DTT","DEGs.H2O2","DEGs.SA")

colnames(log2FCsUV) <- corrcolnames
log2FCsUV <- as.data.frame(log2FCsUV)
log2FCsUV <- log2FCsUV[complete.cases(log2FCsUV),]
log2FCsUV <- log2FCsUV[!duplicated(log2FCsUV),]
corMat <- cor(log2FCsUV,method = "pearson")
png(file.path(outdir,paste0('log2FCnNetCHsUV.png')), 1200, 900, pointsize=12,res = 150)
corrplot(corMat, method = 'circle', order = 'hclust',tl.cex=0.5,title="log2FCnNetCHsFAX",
         tl.pos = 'lt',addCoef.col = "grey",number.cex = 0.5,mar = c(0, 1, 2, 1),
         col=rev(COL2('RdBu', 200)))
dev.off()

###################
#### HeatMaps #####
###################

################### log2FCsFAX
library(ComplexHeatmap)
toHMAP <- log2FCsFAX
col_fun <- colorRamp2(c(min(toHMAP),0,max(toHMAP)), c("blue","white","red"))
column_ha <- HeatmapAnnotation(data=gsub("\\..*","",colnames(log2FCsFAX)),
                               treatment = gsub(".*\\.","",colnames(log2FCsFAX)))
png(file.path(outdir,'HeatMap_log2FCsFAX.png'), 1200, 900, pointsize=12,res = 150)
Heatmap(toHMAP,cluster_columns = T, cluster_rows = T,col=col_fun,top_annotation = column_ha,show_row_names = F,show_row_dend = F)
dev.off()

################### log2FCsUV
toHMAP <- log2FCsUV
col_fun <- colorRamp2(c(min(toHMAP),0,max(toHMAP)), c("blue","white","red"))
column_ha <- HeatmapAnnotation(data=gsub("\\..*","",colnames(log2FCsFAX)),
                               treatment = gsub(".*\\.","",colnames(log2FCsFAX)))
png(file.path(outdir,'HeatMap_log2FCsUV.png'), 1200, 900, pointsize=12,res = 150)
Heatmap(toHMAP,cluster_columns = T, cluster_rows = T,col=col_fun,top_annotation = column_ha,
        show_row_names = F,show_row_dend = F)
dev.off()

################### log2FCs SA UV + FAX

log2FCsSA <- cbind(wholeDF$DEGs.log2FC.SA.YPD,
                   wholeDF$proteomeFAX.log2ratio_Condition4,
                   wholeDF$mRNABPomeFAX.log2ratio_PolyARNAFAXwithSA,
                   wholeDF$mRNABPomeFAX.log2ratio_PolyARNAFAXwithSA - wholeDF$proteomeFAX.log2ratio_Condition4,
                   wholeDF$proteomeUV.log2ratio_Condition4,
                   wholeDF$mRNABPomeUV.log2ratio_PolyARNAUVwithSA,
                   wholeDF$mRNABPomeUV.log2ratio_PolyARNAUVwithSA - wholeDF$proteomeUV.log2ratio_Condition4)
log2FCsSA <- as.data.frame(log2FCsSA)

mycolnames <- c("DEGs.SA",
                "FAX.prot.SA",
                "FAX.RBP.SA", 
                "FAX.NetCh.SA",
                "UV.prot.SA",
                "UV.RBP.SA", 
                "UV.NetCh.SA")
colnames(log2FCsSA) <- mycolnames
log2FCsSA <- log2FCsSA[complete.cases(log2FCsSA),]
log2FCsSA <- log2FCsSA[!duplicated(log2FCsSA),]

toHMAP <- log2FCsSA
col_fun <- colorRamp2(c(min(toHMAP),0,max(toHMAP)), c("blue","white","red"))
column_ha <- HeatmapAnnotation(data=gsub("\\..*","",colnames(log2FCsSA)),
                               treatment = gsub(".*\\.","",colnames(log2FCsSA)))
png(file.path(outdir,'HeatMap_log2FCsSA.png'), 1200, 900, pointsize=12,res = 150)
Heatmap(toHMAP,cluster_columns = T, cluster_rows = T,col=col_fun,
        top_annotation = column_ha,show_row_names = F,show_row_dend = F)
dev.off()

corMat <- cor(log2FCsSA,method = "pearson")
png(file.path(outdir,paste0('corr_log2FCsSA.png')), 1200, 900, pointsize=12,res = 150)
corrplot(corMat, method = 'circle', order = 'hclust',tl.cex=0.5,title="log2FCnNetCHsFAX",
         mar = c(0, 1, 2, 1),tl.pos = 'lt',addCoef.col = "grey",number.cex = 0.5,
         col=rev(COL2('RdBu', 200)))
dev.off()

###########################
#### Enrichment Plots #####
###########################
netFAX.DTT.PsigP.Enr.tsv

enrichResultsDir <- "/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/data/enrichmentRes/GSA"
enrichFiles <- list.files(path = enrichResultsDir,pattern = "*\\.Enr.tsv",full.names = T,)
for (enrichFile in enrichFiles){
  createHMtop10(enrichFile)
}


toIbtissam <- read.xlsx("/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/data/ibtissamPlotting.xlsx")
toIbtissam$data <- factor(toIbtissam$data,levels = c("mRBPome","Proteome","net changes ","DGE"))
toIbtissam$fc <- round(toIbtissam$fc,digits = 2)
ggplot(toIbtissam, aes(x = factor(experiment), fill = factor(gene), y = fc,label=fc)) +
  geom_dotplot(binaxis = "y", stackdir = "center", position = "dodge") +
  facet_grid(~ data) + ylab("Log2FC(SA/Control)") + xlab("") + 
  geom_text(size = 3,position = position_identity(),vjust = "inward")

#nudge_x=0.1,nudge_y=0.1)

# order mRBPome, Proteome, NetC, DGE

#######################
#### Venn Diagram #####
#######################
library(ggVennDiagram);library(ggplot2);

mRNABPomeFAXfcSA <- read.delim("/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/data/geneLists/GSEA/mRNABPomeFAXfcSA.tsv",quote = "")
mRNABPomeUVfcSA <- read.delim("/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/data/geneLists/GSEA/mRNABPomeUVfcSA.tsv",quote = "")

commonmRNABPomeUVFAXSA <- unique(intersect(na.omit(mRNABPomeFAXfcSA$UniprotACC), na.omit(mRNABPomeUVfcSA$UniprotACC))); length(commonmRNABPomeUVFAXSA)
uniqmRNABPomeFAXSA <- unique(setdiff(na.omit(mRNABPomeFAXfcSA$UniprotACC), na.omit(mRNABPomeUVfcSA$UniprotACC))); length(uniqmRNABPomeFAXSA)
uniqmRNABPomeUVSA <- unique(setdiff(na.omit(mRNABPomeUVfcSA$UniprotACC),na.omit(mRNABPomeFAXfcSA$UniprotACC))); length(uniqmRNABPomeUVSA)

toVennList <- list(mRNABPomeFAXfcSA = na.omit(mRNABPomeFAXfcSA$UniprotACC), mRNABPomeUVfcSA = na.omit(mRNABPomeUVfcSA$UniprotACC))

ggvennSA <- ggVennDiagram(toVennList, color = 2, lwd = 0.7) + 
  scale_fill_gradient(low = "#F4FAFE", high = "#4981BF") +
  theme(legend.position = "none")

nameOut <- "/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/data/visualizations/vennFaxUVmRBPome.png"
ggsave(filename = nameOut, plot = ggvennSA, width = 30, height = 15, units = 'cm', dpi = 'print')

##################################
#### Differential Expression #####
##################################

#### Create Matrix for Comparatives #####
ddsFromMatrix <- DESeqDataSetFromMatrix(countData = rsem_expcounts,
                                        colData = phenoData,
                                        design = ~ subtype)

normFactors_1 <- calcNormFactors(rsem_expcounts,method = "TMM")
normFactors <- normFactors_1*(colSums(rsem_expcounts)/mean(colSums(rsem_expcounts)))
names(normFactors) <- colnames(rsem_expcounts)
sizeFactors(ddsFromMatrix) <- normFactors ## Assign TMM normalization factors
dds <- DESeq(ddsFromMatrix)   # replacing outlier value with estimated value as predicted by distrubution using "trimmed mean"

pvalCut <- 0.05
fcCut <- 1
#comparatives <- list(c('Sano','Naive'),c('Sano','CVdet'),c('Sano','CVindet'),c('Naive','CVdet'),c('Naive','CVindet'),c('CVdet','CVindet'))
comparatives <- list(c('HC','Det'),c('HC','Undet'),c('Det','Undet'))

for (ncomparative in c(1:length(comparatives))){
  cond1 <- comparatives[[ncomparative]][1]
  cond2 <- comparatives[[ncomparative]][2]
  classes <- c(cond1,cond2)
  outName <- paste(cond1, paste(cond2,collapse='-'),sep = '_VS_')
  outdir2 <- file.path(outdir,outName)
  dir.create(outdir2,recursive = T,showWarnings = F)
  myContrast <- c("subtype", classes)
  res <- results(dds, contrast = myContrast, pAdjustMethod = "fdr", independentFiltering = T) #### ATTENTION!! Alpha is the FDR cutoff
  res <- as.data.frame(res)
  res <- res[!is.na(res$padj),]
  samplesCompared <- phenoData$samplesName[phenoData$subtype %in% c(cond1,cond2)]
  subphenoData <- phenoData[phenoData$samplesName %in% samplesCompared,]
  TMMexpCompared <- TMMexp[,samplesCompared]
  mergedRes <- merge(res,TMMexpCompared,by.x='row.names',by.y='row.names',all.x = T)
  mergedRes <- mergedRes[order(mergedRes$padj,mergedRes$pvalue),]
  colnames(mergedRes)[1] <- 'miRNA'
  allGenes <- file.path(outdir2,paste0('no_cutoffs_',outName,'.xlsx'))
  write.xlsx(mergedRes, allGenes, rowNames = F,colNames = T)
  
  if (all(mergedRes$pvalue > 0.05)){
    next
  }else if(all(mergedRes$padj > 0.05)){
    sortpvalby <- 'pvalue'
  }else{
    sortpvalby <- 'padj'
  }
  
  outFile <- file.path(outdir2,paste0(outName,'_vulcanoplot.png'))
  volcano <- EnhancedVolcano(mergedRes,x = 'log2FoldChange',y = sortpvalby,
                             lab = mergedRes$miRNA,title = outName,subtitle = '',
                             pCutoff = 0.05,FCcutoff = 1,
                             legendLabels = c('NS', expression(Log[2]~FC),
                                              sortpvalby, expression(p~and~log[2]~FC)))
  ggsave(outFile, plot = volcano, device = 'png', scale = 1, dpi = 300,width = 25,height = 25,units = 'cm')
  DExp <- mergedRes[,sortpvalby] < pvalCut & abs(mergedRes$log2FoldChange) > fcCut
  
  if (!any(DExp)){
    next
  }
  
  resDEGs <- mergedRes[DExp,]
  DEGsFile <- file.path(outdir2,paste0(sortpvalby,pvalCut,'fc',fcCut,'_',outName,'.xlsx'))
  write.xlsx(resDEGs,DEGsFile, rowNames = F,colNames = T)
  
  # 3 create DEGs HeatMap
  toHMP <- t(scale(t(resDEGs[,samplesCompared])))
  row.names(toHMP) <- resDEGs$miRNA
  col_fun <- colorRamp2(c(min(toHMP),0,max(toHMP)), c("blue", "white","red"))
  column_ha <- HeatmapAnnotation(type = subphenoData$subtype,annotation_name_gp= gpar(fontsize = 10))
  outFile <- file.path(outdir2,paste0(sortpvalby,pvalCut,'fc',fcCut,'_',outName,'_HeatMap.png'))
  mywidth <- ncol(toHMP)
  myheight <- nrow(toHMP)
  png(outFile, width = ncol(toHMP)*8,height = myheight*5+50, units = 'mm', pointsize=30,res = 90)
  print(
    Heatmap(toHMP,cluster_rows = F, name = outName, cluster_columns = T,
            show_row_dend = F, show_column_dend = F, show_row_names = T,
            column_names_rot = 65,column_names_side = "top",
            column_names_gp = gpar(fontsize = 10),col=col_fun,top_annotation = column_ha,
            width = mywidth*unit(5, "mm") ,height = myheight*unit(5, "mm") )
    #,width = unit(8, "cm"), height = unit(16, "cm"))#,heatmap_width = unit(10, "cm"), heatmap_height = unit(20, "cm"))
  )
  dev.off()
}

