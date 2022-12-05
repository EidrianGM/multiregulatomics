
###### Analysis of NOX mRBPome for background Removal

# The idea is to obtain the a log2FC of the raw data. So NO normalization and NO imputation and NO samples removal.

rawRBPomeNOXfile <- "/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/FinalData/Raw/MQtoSQ.csv"
rawRBPomeNOXSamplesInfo <- read.delim(rawRBPomeNOXfile, sep = ",")
myDF <- rawRBPomeNOXSamplesInfo

normDCol <- grep("Norm",colnames(myDF)) # Raw Normalized Values are not the ones!
rawDCol <- grep("Raw",colnames(myDF))
spcDCol <- grep("Spec",colnames(myDF))
expresionMatrix <- myDF[rawDCol:(spcDCol-1)] 
metadata <-   myDF[1:(normDCol-1)]
colnames(metadata) <- metadata[2,]
expresionMatrix[2,] <- gsub("_.*","",expresionMatrix[2,])
subPhenData <- as.data.frame(t(expresionMatrix[1:2,])); colnames(subPhenData) <- c("class","id"); 
subPhenData$class <- gsub(".*\\s","",subPhenData$class); conditionCol <- "class"
rownames(subPhenData) <- 1:nrow(subPhenData)
metadata <- metadata[3:nrow(metadata),]
expresionMatrix <- expresionMatrix[3:nrow(expresionMatrix),]
colnames(expresionMatrix) <- subPhenData$id
subPhenData$class <- gsub("with.*|no.*","",subPhenData$class)

resultsDF <- read.delim("/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/data/new/PolyA_RIC_nonorm_20220928/P445_AG_PolyA_40_20220617_gmin_UV_nonorm_p2_NOX_Final/P445_AG_PolyA_40_20220617_gmin_UV_nonorm_p2_NOX_PROTEIN.tsv")

all(resultsDF$proteinName %in% metadata$Accession) # TRUE == All proteins in results are in raw
expresionMatrix <- type.convert(expresionMatrix)

filteredByAlex <- which(metadata$Accession %in% resultsDF$proteinName) ## to keep the same proteins as Alex
metadata <- metadata[filteredByAlex,]
expresionMatrix <- expresionMatrix[filteredByAlex,]

rawDataIdx <- apply(expresionMatrix,2, FUN="median", na.rm=T)
globalNormFactors <- as.numeric(rawDataIdx[1]) / as.numeric(rawDataIdx)
normData <- sweep(expresionMatrix, 2, globalNormFactors, FUN="*")

expresionMatrix <- normData
protExprMatrix <- as.data.frame(summarise_all(group_by(expresionMatrix,metadata$Accession),sum))
row.names(protExprMatrix) <- protExprMatrix$`metadata$Accession`
protExprMatrix$`metadata$Accession` <- NULL

eCtrl <- protExprMatrix[,subPhenData$class == "NOX"]
eFAX <- protExprMatrix[,subPhenData$class == "FAX"]
eUVX <- protExprMatrix[,subPhenData$class == "UV"]

method <- "median"
log2medianNOX <- apply(log2(eCtrl),1,method, na.rm=T)
log2medianFAX <- apply(log2(eFAX),1,method, na.rm=T)
log2medianUVX <- apply(log2(eUVX),1,method, na.rm=T)
log2ratiosFAX <- medianFAX - medianNOX
log2ratiosUVX <- medianUVX - medianNOX

library(openxlsx)
outdir <- "/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/FinalData/BackgroundRemoval"
backgroundRBPomeDF <- cbind(protExprMatrix,log2medianNOX,log2medianFAX,log2medianUVX,log2ratiosFAX,log2ratiosUVX)
backgroundRBPomeDF <- cbind(proteinName=rownames(backgroundRBPomeDF),backgroundRBPomeDF)
saveTablesTsvExc(backgroundRBPomeDF,outdir,completeNdedup = F,excel = T,bycompleteFC = F,rownames = F)