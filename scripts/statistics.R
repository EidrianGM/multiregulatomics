# Evaluation and Improvement of Quantification Accuracy in Isobaric Mass Tag-Based Protein Quantification Experiments
# https://pubs.acs.org/doi/10.1021/acs.jproteome.6b00066

## Shall we remove this one? 
# P445_Prot_DIA_iBAQ_noOutliers_nocomb_20220826_Final	P445_Prot_DIA_iBAQ_noOutliers_nocomb_20220826_UV_Final	UVX stress-induced proteome  	Proteome 	UVX	SA	iBAQ 22	22	UV+ sodium acetate (1)	Yes 	4	8/26/2022	old

projectBaseDir <- "/home/eidriangm/Desktop/toDo/surrey/multiregulatomics"
setwd(projectBaseDir)

mappingFile <- "/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/data/mapping/wholeDF.xlsx"
mappingDF <- read.xlsx(mappingFile)
mappingDF <- mappingDF[c(1:17)]
mappingDF <- mappingDF[!duplicated(mappingDF),]

basePhenoData <- read.xlsx("data/phenoData.xlsx")
basePhenoData <- as.data.frame(apply(basePhenoData,2,function(x)gsub('\\s+', '',x)))

normalised <- "Yes";  mytreatment <- "SA"; control <- "no"
################################################################################
########################   RNABPOME   ##########################################
################################################################################
datatype <- "mRBPome"

#### UVX
crosslink <- "UVX"

mrbpomeUVXphenoD <- basePhenoData[basePhenoData$Treatment %in% c(control,mytreatment) & 
                                    basePhenoData$normalised == normalised & 
                                    basePhenoData$Data == datatype & 
                                    basePhenoData$Crosslinking == crosslink,]

subfolderpath <- unique(paste(projectBaseDir,"data",mrbpomeUVXphenoD$folderdate,mrbpomeUVXphenoD$Folder.name,mrbpomeUVXphenoD$Subfolder.name,sep = "/"))
dataFile <- list.files(subfolderpath,pattern = "*PROTEIN.tsv",full.names = T)

ourPhenoData <- mrbpomeUVXphenoD[!mrbpomeUVXphenoD$Run.number %in% c(16,3),]
outdir <- paste0("/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/adrian/",datatype,crosslink,mytreatment,"no16-3")
tagname <- paste0(datatype,crosslink,mytreatment,"no16-3")
basicAnalysis(dataFile,ourPhenoData,outdir,tagname,mytreatment,mappingDF)


#### FAX
crosslink <- "FAX"

ourPhenoData <- basePhenoData[basePhenoData$Treatment %in% c(control,mytreatment) & 
                                basePhenoData$normalised == normalised & 
                                basePhenoData$Data == datatype & 
                                basePhenoData$Crosslinking == crosslink,]

subfolderpath <- unique(paste(projectBaseDir,ourPhenoData$folderdate,ourPhenoData$Folder.name,ourPhenoData$Subfolder.name,sep = "/"))
dataFile <- list.files(subfolderpath,pattern = "*PROTEIN.tsv",full.names = T)

outdir <- paste0("/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/adrian/",datatype,crosslink,mytreatment)
basicAnalysis(dataFile,ourPhenoData,outdir,mytreatment,mappingDF)


################################################################################
########################   PROTEOME   ##########################################
################################################################################
datatype <- "Proteome"

#### FAX
crosslink <- "FAX"

ourPhenoData <- basePhenoData[basePhenoData$Treatment %in% c(control,mytreatment) & 
                                basePhenoData$normalised == normalised & 
                                basePhenoData$Data == datatype & 
                                basePhenoData$Crosslinking == crosslink,]

subfolderpath <- unique(paste(projectBaseDir,ourPhenoData$folderdate,ourPhenoData$Folder.name,ourPhenoData$Subfolder.name,sep = "/"))
dataFile <- list.files(subfolderpath,pattern = "*PROTEIN.tsv",full.names = T)

ourPhenoData <- ourPhenoData[!ourPhenoData$Run.number %in% c(25,36),]
outdir <- paste0("/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/adrian/",datatype,crosslink,mytreatment,"no25-36")
basicAnalysis(dataFile,ourPhenoData,outdir,mytreatment,mappingDF)

#### FAX
crosslink <- "UVX"

ourPhenoData <- basePhenoData[basePhenoData$Treatment %in% c(control,mytreatment) & 
                                basePhenoData$normalised == normalised & 
                                basePhenoData$Data == datatype & 
                                basePhenoData$Crosslinking == crosslink,]

subfolderpath <- unique(paste(projectBaseDir,ourPhenoData$folderdate,ourPhenoData$Folder.name,ourPhenoData$Subfolder.name,sep = "/"))
dataFile <- list.files(subfolderpath,pattern = "*PROTEIN.tsv",full.names = T)

outdir <- paste0("/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/adrian/",datatype,crosslink,mytreatment)
basicAnalysis(dataFile,ourPhenoData,outdir,mytreatment,mappingDF)

################################################################################
########################   BACKGROUND   ########################################
################################################################################

mappingFile <- "/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/yeastReference/mappingFile_8.9.22.xlsx"
mappingDF <- read.xlsx(mappingFile)

normalised <- "No"; datatype <- "mRBPome"; control <- "no"; mytreatment <- "SA"; crosslink <- "NOX";

selectedPhenoD <- basePhenoData[basePhenoData$Treatment %in% c(control,mytreatment) & 
                                  basePhenoData$normalised == normalised & 
                                  basePhenoData$Data == datatype & 
                                  basePhenoData$Crosslinking == crosslink &
                                  basePhenoData$folderdate == "old",]
subfolderpath <- unique(paste(projectBaseDir,selectedPhenoD$folderdate,selectedPhenoD$Folder.name,selectedPhenoD$Subfolder.name,sep = "/"))
subfolderpath <- subfolderpath[c(1,2,5)]
dataFile <- list.files(subfolderpath,pattern = "*PROTEIN.tsv",full.names = T)
NOXmRBPDF <- read.delim(dataFile)

FAXmRBPfile <- "/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/data/old/P445_AG_PolyA_40_all_3_cond_20220611_Final/P445_AG_PolyA_40_20220617_gmin_FAX_p2_Final/P445_AG_PolyA_40_20220617_gmin_FAX_p2_PROTEIN.tsv"
FAXmRBPDF <- read.delim(FAXmRBPfile)
length(unique(FAXmRBPDF$ac))
UVXmRBPfile <- "/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/data/old/P445_AG_PolyA_40_all_3_cond_20220611_Final/P445_AG_PolyA_40_20220617_gmin_UV_p2_no511_Final/P445_AG_PolyA_40_20220617_gmin_UV_p2_no511_PROTEIN.tsv"
UVXmRBPDF <- read.delim(UVXmRBPfile)
length(unique(UVXmRBPDF$ac))

all(FAXmRBPDF$proteinName %in% NOXmRBPDF$proteinName)
all(UVXmRBPDF$proteinName %in% NOXmRBPDF$proteinName)

FAXbackgDF <- NOXmRBPDF[NOXmRBPDF$proteinName %in% FAXmRBPDF$proteinName,]
UVXbackgDF <- NOXmRBPDF[NOXmRBPDF$proteinName %in% UVXmRBPDF$proteinName,]

outdir <- "/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/data/background"
saveTablesTsvExc(FAXbackgDF,outdir,completeNdedup=F,excel=F,bycompleteFC=F,rownames = F)
saveTablesTsvExc(UVXbackgDF,outdir,completeNdedup=F,excel=F,bycompleteFC=F,rownames = F)

FAXbackgPhenoD <- selectedPhenoD[selectedPhenoD$Condition %in% c("1","2"),]
UVXbackgPhenoD <- selectedPhenoD[selectedPhenoD$Condition %in% c("1","3"),]
saveTablesTsvExc(FAXbackgPhenoD,outdir,completeNdedup=F,excel=F,bycompleteFC=F,rownames = F)
saveTablesTsvExc(UVXbackgPhenoD,outdir,completeNdedup=F,excel=F,bycompleteFC=F,rownames = F)


ncol(FAXbackgDF)

conditionCol <- "Condition"
pCutoff <- 0.05
FCcutoff <- 3
control <- "1"
mytreatment <- "2"

outdir <- "/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/adrian/background/FAX"
tagname <- "mRBPome_NOX_FAX"
basicAnalysis(FAXbackgDF,FAXbackgPhenoD,outdir,tagname,"2",mappingDF,"1",conditionCol,pCutoff,FCcutoff)

outdir <- "/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/adrian/background/UVX"
tagname <- "mRBPome_NOX_UVX"
basicAnalysis(UVXbackgDF,UVXbackgPhenoD,outdir,tagname,"3",mappingDF,"1",conditionCol,pCutoff,FCcutoff)




################################################################################
################################################################################ 

selectedPhenoD <- selectedPhenoD[]


myDFfax <- read.delim(dataFile[grep("FAX_p2",dataFile)])
myDFuvx <- read.delim(dataFile[grep("UV_p2",dataFile)])
myDFnox <- read.delim(dataFile[grep("p2_NOX",dataFile)])

colnames(myDFfax) <- paste0(colnames(myDFfax),".fax")
colnames(myDFuvx) <- paste0(colnames(myDFuvx),".uvx")
colnames(myDFnox) <- paste0(colnames(myDFnox),".nox")

allDFs <- merge(myDFnox,myDFuvx,by.x = "proteinName.nox", by.y = "proteinName.uvx")
allDFs <- merge(allDFs,myDFfax,by.x = "proteinName.nox", by.y = "proteinName.fax")

# allDFs <- merge(myDFnox,myDFuvx,by.x = "proteinName.nox", by.y = "proteinName.uvx",all = T)
# allDFs <- merge(allDFs,myDFfax,by.x = "proteinName.nox", by.y = "proteinName.fax",all = T)

expresionMatrix <- allDFs[,grep("X[0-9]+",colnames(allDFs))]
ncol(expresionMatrix)

selectedPhenoD

proteinOriginalnfo <- grep(paste0(gsub("mRBPome","mRNABPome",datatype),gsub("UVX","UV",crosslink),".allAccessions"),colnames(mappingDF))
if (length(proteinOriginalnfo) == 0){
  proteinOriginalnfo <- grep(paste0(tolower(datatype),gsub("UVX","UV",crosslink),".allAccessions"),colnames(mappingDF))
} 
subMapDF <- mappingDF[c(1:17,proteinOriginalnfo)]
subMapDF <- subMapDF[!is.na(subMapDF[,ncol(subMapDF)]),]

ourPhenoData <- selectedPhenoD
outdir <- paste0("/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/adrian/",datatype,crosslink,mytreatment)
basicAnalysis(myDF,ourPhenoData,outdir,mytreatment,subMapDF)

ourPhenoData <- selectedPhenoD[!selectedPhenoD$Run.number %in% c(25,36),]
outdir <- paste0("/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/adrian/",datatype,crosslink,mytreatment,"no25-36")
basicAnalysis(myDF,ourPhenoData,outdir,mytreatment,subMapDF)



################################################################################
########################   POSTANALYSIS   ######################################  
################################################################################

protNOX <- read.delim("/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/data/new/P445_Prot_DIA_iBAQ_nocomb_20220826_all_NOX_Final/P445_Prot_DIA_iBAQ_nocomb_20220826_NOX_Final2/P445_Prot_DIA_iBAQ_nocomb_20220826_NOX_Final2_PROTEIN.tsv")
row.names(protNOX) <- protNOX$proteinName
expresionMatrix <- protNOX[,grep("X[0-9]+",colnames(protNOX))]
runnumbers <- gsub("X0+","",unlist(lapply(strsplit(colnames(expresionMatrix),"_"), function(x) return(x[[1]]))))
subPhenData <- as.data.frame(cbind(cols=colnames(expresionMatrix),"Run.number"=runnumbers))
ourPhenoData <- ourPhenoData[ourPhenoData$Subfolder.name == "P445_Prot_DIA_iBAQ_nocomb_20220826_NOX_Final2",]
subPhenData <- merge(subPhenData,ourPhenoData,by="Run.number")
subPhenData <- subPhenData[order(subPhenData$Treatment),]
expresionMatrix <- expresionMatrix[,match(subPhenData$cols,colnames(expresionMatrix))]
row.names(subPhenData) <- subPhenData$cols

### PREPARE DATA 
subPhenData$Treatment <- as.factor(subPhenData$Treatment)
phenoD <- AnnotatedDataFrame(data=subPhenData[,"Treatment",drop=FALSE])
featD <- AnnotatedDataFrame(data=protNOX[,"proteinName",drop = FALSE])
eset <- ExpressionSet(assayData = as.matrix(expresionMatrix),phenoData=phenoD,featureData = featD)
exprs(eset) <- log2(exprs(eset))

## NOT PAIRWISE
design <- model.matrix(~0+Treatment, data=pData(eset))
colnames(design) <- gsub("^Treatment","",colnames(design))
fit <- lmFit(eset,design)
contrastMatrix <- makeContrasts(contrasts=paste(caseConditions,"-", controlCondition),levels=design)
colnames(contrastMatrix) <- caseConditions
fitContrasts <- eBayes(contrasts.fit(fit,contrastMatrix)) # fit contrasts coefficients
pvalues <- data.frame(fitContrasts$p.value[,caseConditions])
names(pvalues) <- caseConditions


## PAIRWISE
controlCondition <- "no"
caseConditions <- c("DTT","H202","SA")

pairWpvalues <- data.frame(row.names=featureNames(eset))
for(cC in caseConditions){
  ### at least one replicate of one condition requires
  selCol <- unlist(pData(eset)$Treatment %in% c(control,mytreatment))
  if(sum(selCol) > 2 ){
    esetPair <- eset[,selCol]
    # add subject term to allow for paired t-statistic
    fit <-eBayes(lmFit(esetPair, model.matrix(~factor(esetPair$Treatment), data=pData(esetPair))))
    p <- fit$p.value[,2]
  }else{
    p <- rep(NA,nrow(eset))
  }
  pairWpvalues <- cbind(pairWpvalues,p)
  colnames(pairWpvalues)[ncol(pairWpvalues)] <- paste0("pVal",cC)
  # pairW <- cbind(pairW,p.adjust(p,method="BH"))
  # colnames(pairW)[ncol(pairW)] <- paste0("pValAdj",cC)
}

ratios <- data.frame(row.names=featureNames(eset))
eCtrl <- subset(exprs(eset),select=which(pData(eset)$Treatment == controlCondition))
for(cC in caseConditions){
  eCase <- subset(exprs(eset),select=which(pData(eset)$Treatment == cC))
  ratios <- cbind(ratios,apply(log2(eCase),1,"median", na.rm=T) - apply(log2(eCtrl),1,"median", na.rm=T))
}
names(ratios) <- caseConditions
resultsDF <- merge(pvalues,ratios,by="row.names")

View(resultsDF)


ENO1control <- c(10,50,30)
ENO1n2control <- c(100,70,84)
ENO1case <- c(6,12,9)

log2(mean(ENO1control)) - log2(mean(ENO1case))
log2(mean(ENO1control+ENO1n2control)) - log2(mean(ENO1case+ENO1n2control))


- 
  
  
  ENO1sa <- c(10,50,30)
ENO1n2sa <- c(100,70,84)
ENO2sa <- c(6,12,9)



## moderated t-test
# do not perform stat test for filtered out features

out$pValue <- out$ratio
out$pValue[match(rownames(eset) , rownames(out$pValue)),] = NA
out$qValue <- 	out$pValue 
out$FPValue <-  rep(NA,nrow(out$pValue))
out$FQValue <-  rep(NA,nrow(out$pValue))

## we need at least two conditions
if(length(unique(pData(eset)$condition)) > 1){
  
  out$pValue[match(rownames(eset)[sel] , rownames(out$pValue)),] <- getAllEBayes(eset[sel,],adjust=F,method=method,...)	
  #out$pValue[rownames(eset)[sel], ] <- getAllEBayes(eset[sel, ], adjust = F, method = method)
  
  # apply ratio cut-off for adjustment. 
  adjustFilter <- data.frame((abs(out$ratio) < log2(fcThrs)))
  out$qValue[match(rownames(eset)[sel] , rownames(out$qValue)) ,] <- getAllEBayes(eset[sel,],adjust=T,method=method, adjustFilter=subset(adjustFilter,subset=sel),... )
  #out$qValue[rownames(eset)[sel], ] <- getAllEBayes(eset[sel, ], adjust = T, method = method, adjustFilter = subset(adjustFilter, subset = sel))
  
  out$FPValue[sel] = getFTestPValue(eset[sel,],...)
  
  # apply ratio cut-off for adjustment. Keep if one cond meets ratio cut cut-off critera 
  adjustFilterF = rowSums(!adjustFilter) > 0
  out$FQValue[adjustFilterF & sel] = p.adjust(out$FPValue[adjustFilterF & sel],method="BH")
  
}

#colnames(mappingDF)[grep("proteomeUV",colnames(mappingDF))]

uvxSA <- read.xlsx("/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/adrian/mRBPomeUVXSA/diffExprMapped.xlsx")
no3918uvxSA <- read.xlsx("/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/adrian/mRBPomeUVXSAno39-18/diffExprMapped.xlsx")
no163uvxSA <- read.xlsx("/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/adrian/mRBPomeUVXSAno16-3/diffExprMapped.xlsx")
no393uvxSA <- read.xlsx("/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/adrian/mRBPomeUVXSAno39-3/diffExprMapped.xlsx")

uvxSAgenes <- uvxSA$mRNABPomeUV.allAccessions[uvxSA$FDR.SA < 0.05]
no3918uvxSAgenes <- no3918uvxSA$mRNABPomeUV.allAccessions[no3918uvxSA$FDR.SA < 0.05]
no163uvxSAgenes <- no163uvxSA$mRNABPomeUV.allAccessions[no163uvxSA$FDR.SA < 0.05]
no393uvxSAgenes <- no393uvxSA$mRNABPomeUV.allAccessions[no393uvxSA$FDR.SA < 0.05]
toVennList <- list(uvxSA = uvxSAgenes, no3918 = no3918uvxSAgenes, no163 = no163uvxSAgenes,no393=no393uvxSAgenes)

ggvennSA <- ggVennDiagram(toVennList, color = 2, lwd = 0.7) + 
  scale_fill_gradient(low = "#F4FAFE", high = "#4981BF") +
  theme(legend.position = "none")
nameOut <- "/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/adrian/mrbpomeUVXSAx4_FDR05_venn.png"
ggsave(filename = nameOut, plot = ggvennSA, width = 30, height = 15, units = 'cm', dpi = 'print')

protFAXSA <- read.xlsx("/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/adrian/ProteomeFAXSA/diffExprMapped.xlsx")
no2536protFAXSA <- read.xlsx("/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/adrian/ProteomeFAXSAno25-36/diffExprMapped.xlsx")

protFAXSAgenes <- protFAXSA$proteomeFAX.allAccessions[protFAXSA$FDR.SA < 0.05]
no2536protFAXSAgenes <- no2536protFAXSA$proteomeFAX.allAccessions[no2536protFAXSA$FDR.SA < 0.05]
toVennList <- list(FAXSA = protFAXSAgenes, no2536 = no2536protFAXSAgenes)

ggvennSA <- ggVennDiagram(toVennList, color = 2, lwd = 0.7) + 
  scale_fill_gradient(low = "#F4FAFE", high = "#4981BF") +
  theme(legend.position = "none")
nameOut <- "/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/adrian/protFAXSA_FDR05_venn.png"
ggsave(filename = nameOut, plot = ggvennSA, width = 30, height = 15, units = 'cm', dpi = 'print')


############ MAPPING and Merging

mappingFile <- "/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/data/mapping/wholeDF.xlsx"
mappingDF <- read.xlsx(mappingFile)

extraCols <- c("DEGs.adj.pval.SA.YPD","DEGs.log2FC.SA.YPD")
mappingDF <- mappingDF[c(1:17, which(colnames(mappingDF) %in% extraCols))]


######
######  UVX

newRBPomeUVXSA <- read.xlsx("/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/adrian/mRBPomeUVXSAno16-3/diffExprMapped.xlsx")
colnames(newRBPomeUVXSA)[(ncol(newRBPomeUVXSA)-2):ncol(newRBPomeUVXSA)]  <- paste0("newRBP.UVX.SA.",colnames(newRBPomeUVXSA)[(ncol(newRBPomeUVXSA)-2):ncol(newRBPomeUVXSA)])
newProteomeUVXSA <- read.xlsx("/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/adrian/ProteomeUVXSA/diffExprMapped.xlsx")
colnames(newProteomeUVXSA)[(ncol(newProteomeUVXSA)-2):ncol(newProteomeUVXSA)]  <- paste0("newPROT.UVX.SA.",colnames(newProteomeUVXSA)[(ncol(newProteomeUVXSA)-2):ncol(newProteomeUVXSA)])
newUVXSA <- merge(newProteomeUVXSA,newRBPomeUVXSA,all=T)
newUVXSAFull <- merge(mappingDF,newUVXSA,all=T)

newUVXSAFull <- newUVXSAFull[rowSums(is.na(newUVXSAFull[(ncol(newUVXSAFull)-7):ncol(newUVXSAFull)])) != length((ncol(newUVXSAFull)-7):ncol(newUVXSAFull)),]
newUVXSAFull <- newUVXSAFull[!duplicated(newUVXSAFull),]

newUVXSAFull$newUVXSA.netchanges <- newUVXSAFull$newRBP.UVX.SA.log2FC.SA - newUVXSAFull$newPROT.UVX.SA.log2FC.SA

outdir <- "/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/adrian"
saveTablesTsvExc(newUVXSAFull,outdir,completeNdedup=F,excel=T,bycompleteFC=F,rownames = F)

length(unique(newUVXSAFull$UniprotACC[newUVXSAFull$newPROT.UVX.SA.FDR.SA <= 0.05])) - 1 
length(unique(newUVXSAFull$UniprotACC[newUVXSAFull$newRBP.UVX.SA.FDR.SA <= 0.05])) - 1 

toHMP <- newUVXSAFull[,c("DEGs.log2FC.SA.YPD","newPROT.UVX.SA.log2FC.SA","newRBP.UVX.SA.log2FC.SA","newUVXSA.netchanges")]
toHMP <- toHMP[complete.cases(toHMP),]
colnames(toHMP) <- c("DEGs","Prot","RBP","NetChanges")
col_fun <- colorRamp2(c(min(toHMP,na.rm = T),0,max(toHMP,na.rm = T)), c("blue","white","red"))
png(file.path(outdir,'newUVXSAFull_HeatMap.png'), 1200, 900, pointsize=12,res = 150)
Heatmap(toHMP,cluster_columns = F, cluster_rows = T,col=col_fun,show_row_names = F,show_row_dend = F,
        column_title = "UVX SA FCs", name = "FCs", column_names_rot = 0)
dev.off()

######
######  FAX

newProtFAXSA <- read.xlsx("/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/adrian/ProteomeFAXSAno25-36/diffExprMapped.xlsx")
colnames(newProtFAXSA)[(ncol(newProtFAXSA)-2):ncol(newProtFAXSA)]  <- paste0("newPROT.FAX.SA.",colnames(newProtFAXSA)[(ncol(newProtFAXSA)-2):ncol(newProtFAXSA)])
newRBPomeFAXSA <- read.xlsx("/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/adrian/mRBPomeFAXSA/diffExprMapped.xlsx")
colnames(newRBPomeFAXSA)[(ncol(newRBPomeFAXSA)-2):ncol(newRBPomeFAXSA)]  <- paste0("newRBP.FAX.SA.",colnames(newRBPomeFAXSA)[(ncol(newRBPomeFAXSA)-2):ncol(newRBPomeFAXSA)])
newFAXSA <- merge(newProtFAXSA,newRBPomeFAXSA,all=T)
newFAXSAFull <- merge(mappingDF,newFAXSA,all=T)

newFAXSAFull <- newFAXSAFull[rowSums(is.na(newFAXSAFull[(ncol(newFAXSAFull)-7):ncol(newFAXSAFull)])) != length((ncol(newFAXSAFull)-7):ncol(newFAXSAFull)),]
newFAXSAFull <- newFAXSAFull[!duplicated(newFAXSAFull),]
newFAXSAFull$newFAXSA.netchanges <- newFAXSAFull$newRBP.FAX.SA.log2FC.SA - newFAXSAFull$newPROT.FAX.SA.log2FC.SA

outdir <- "/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/adrian"
saveTablesTsvExc(newFAXSAFull,outdir,completeNdedup=F,excel=T,bycompleteFC=F,rownames = F)

newFAXSAFull$Un


toHMP <- newFAXSAFull[,c("DEGs.log2FC.SA.YPD","newPROT.FAX.SA.log2FC.SA","newRBP.FAX.SA.log2FC.SA","newFAXSA.netchanges")]
toHMP <- toHMP[complete.cases(toHMP),]
colnames(toHMP) <- c("DEGs","Prot","RBP","NetChanges")
col_fun <- colorRamp2(c(min(toHMP,na.rm = T),0,max(toHMP,na.rm = T)), c("blue","white","red"))
png(file.path(outdir,'newFAXSAFull_HeatMap.png'), 1200, 900, pointsize=12,res = 150)
Heatmap(toHMP,cluster_columns = F, cluster_rows = T,col=col_fun,show_row_names = F,show_row_dend = F,
        column_title = "FAX SA FCs", name = "FCs", column_title_rot = 0,column_names_rot = 0)
dev.off()


length(unique(newFAXSAFull$UniprotACC[newFAXSAFull$newPROT.FAX.SA.FDR.SA <= 0.05])) - 1 
length(unique(newFAXSAFull$UniprotACC[newFAXSAFull$newRBP.FAX.SA.FDR.SA <= 0.05])) - 1 


sum(na.omit(newFAXSAFull$newRBP.FAX.SA.FDR.SA <= 0.05))
sum(na.omit(newFAXSAFull$newPROT.FAX.SA.FDR.SA <= 0.05))

##########
########## TO ESCHER
##########

outdir <- "/home/eidriangm/Desktop/toDo/surrey/multilayerpathwaysviz/data/pathway/toEscher"

read.xlsx("/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/adrian/newUVXSAFull.xlsx")
myData <- fromJSON(file="/home/eidriangm/Desktop/toDo/surrey/multilayerpathwaysviz/data/pathway/genes_iMM904.json")
ORFgenesinMM904 <- lapply(myData$results, FUN = function(x) return(x$bigg_id))
ORFgenesinMM904 <- unlist(ORFgenesinMM904)

NamegenesinMM904 <- lapply(myData$results, FUN = function(x) return(x$name))
NamegenesinMM904 <- unlist(NamegenesinMM904)


newFAXSAFull <- read.xlsx("/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/adrian/newFAXSAFull.xlsx")
mydataORF <- newFAXSAFull[,grep("ORF|log2|netchanges",colnames(newFAXSAFull))]

mydataORF$ORFgene <- NULL
mydataORF$Gene.Names..ORF. <- NULL
newFAXSAFull$ORF <- gsub("-","_",newFAXSAFull$ORF)
mydataORF <- mydataORF[!duplicated(mydataORF),]
mapORFsNotInOurData <- ORFgenesinMM904[!ORFgenesinMM904 %in% newFAXSAFull$ORF]
mydataORF <- mydataORF[mydataORF$ORF %in% ORFgenesinMM904,]
any(duplicated(mydataORF$ORFgene))
colnames(mydataORF) <- c("","Transcriptome","Proteome","mRNABPome","NetChangs")

outFile <- file.path(outdir,"Tr.ORFsEscherFAXSAFull.csv")
write.table(mydataORF[,c(1,2)],outFile,sep = ",",na = "0",row.names = F,quote = FALSE)

outFile <- file.path(outdir,"Prot.ORFsEscherFAXSAFull.csv")
write.table(mydataORF[,c(1,3)],outFile,sep = ",",na = "0",row.names = F,quote = FALSE)

outFile <- file.path(outdir,"mRBP.ORFsEscherFAXSAFull.csv")
write.table(mydataORF[,c(1,4)],outFile,sep = ",",na = "0",row.names = F,quote = FALSE)

outFile <- file.path(outdir,"netCH.ORFsEscherFAXSAFull.csv")
write.table(mydataORF[,c(1,5)],outFile,sep = ",",na = "0",row.names = F,quote = FALSE)


mydataName <- newFAXSAFull[,grep("GeneName|log2",colnames(newFAXSAFull))]
mydataName <- mydataName[!duplicated(mydataName),]
mapNamessNotInOurData <- NamegenesinMM904[!NamegenesinMM904 %in% newFAXSAFull$GeneName]
mydataName <- mydataName[mydataName$GeneName %in% NamegenesinMM904,]
any(duplicated(mydataName$GeneName))
outFile <- file.path(outdir,"NamesEscherFAXSAFull.csv")
colnames(mydataName) <- c("","Transcriptome","Proteome","mRNABPome")
write.table(mydataName,outFile,sep = ",",na = "0",row.names = F,quote = FALSE)



mydataName <- mydataName[!duplicated(mydataName),]


write.table(mydata,outFile,sep = ",",na = "0",row.names = F,quote = FALSE)
newFAXSAFull$ORFgene


all(genesinMM904 %in% newFAXSAFull$ORFgene)


NamegenesinMM904[!NamegenesinMM904 %in% newFAXSAFull$GeneName]
