################################################################################
############################ Vizs of RAW #######################################
################################################################################

#### UV RBPome
inFile <- "/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/data/raw/MQtoSQ_UV.csv"
tagname <- "UV-RBPome"
getPlotsRawData(inFile,tagname)

#### FAX RBPome
inFile <- "/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/data/raw/MQtoSQ_FAX.csv"
tagname <- "FAX-RBPome"
getPlotsRawData(inFile,tagname)

#### NOX RBPome
inFile <- "/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/data/raw/MQtoSQ.csv"
tagname <- "NOX-RBPome"
getPlotsRawData(inFile,tagname)

#### Proteome
nFile <- "/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/data/raw/SNtoSQ-Prot40_SafeQuant_Input.csv"
tagname <- "Proteome"
getPlotsRawData(inFile,tagname)

################################################################################
############################ Imputation ########################################
################################################################################

#### UV RBPome
inFile <- "/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/data/raw/MQtoSQ_UV.csv"
myDF <- read.csv(inFile)

normDCol <- grep("Norm",colnames(myDF))
rawDCol <- grep("Raw",colnames(myDF))
expresionMatrix <- myDF[normDCol:(rawDCol-1)]
expresionMatrix[2,] <- gsub("_.*","",expresionMatrix[2,])
subPhenData <- as.data.frame(t(expresionMatrix[1:2,])); colnames(subPhenData) <- c("class","id"); rownames(subPhenData) <- subPhenData[,2]
colnames(expresionMatrix) <- expresionMatrix[2,]; expresionMatrix <- expresionMatrix[3:nrow(expresionMatrix),]; expresionMatrix <- type.convert(expresionMatrix)
subPhenData$class <- gsub(".*\\s","",subPhenData$class)

eset <- sqNormalize(eset, method=method)

if("global" %in% method){
  globalNormFactors <- getGlobalNormFactors(esetNorm,method=method)
  ### add normalization factors to ExpressionSet
  pData(esetNorm)$globalNormFactors <- globalNormFactors
  esetNorm <- globalNormalize(esetNorm,globalNormFactors)
}
tagname <- "UV-RBPome"


#### FAX RBPome
inFile <- "/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/data/raw/MQtoSQ_FAX.csv"
tagname <- "FAX-RBPome"
getPlotsRawData(inFile,tagname)

#### NOX RBPome
inFile <- "/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/data/raw/MQtoSQ.csv"
tagname <- "NOX-RBPome"
getPlotsRawData(inFile,tagname)

#### Proteome
nFile <- "/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/data/raw/SNtoSQ-Prot40_SafeQuant_Input.csv"
tagname <- "Proteome"
getPlotsRawData(inFile,tagname)

################################################################################
############################ Vizs of Imputed ###################################
################################################################################


################################################################################
############################ What's DONE? ######################################
################################################################################

compareRawNDone <- function(doneFile,rawFile,crosslink=NA){
  doneDF <- read.delim(doneFile,nrows = 3)
  mysamplesDone <- colnames(doneDF)[grep("^X",colnames(doneDF))]
  mysamplesDone <- unique(gsub("X0|_.*","",mysamplesDone))
  
  rawDF <- read.csv(rawFile,nrows = 3)
  if (!is.na(crosslink)){
    rawDF <- rawDF[,grep(crosslink,rawDF[1,])]
    mysamples <- rawDF[2,]
  }else{
    mysamples <- rawDF[2,][grep("PolyA",rawDF[2,])]
  }
  mysamples <- unique(gsub("^0|_.*","",mysamples))
  if (length(mysamples) == length(mysamplesDone) &  all(mysamples %in% mysamplesDone)){
    cat("All Okey")
  }else{
    cat("Samples Discarded: ")
    cat(paste0(setdiff(mysamples,mysamplesDone),collapse = " "))
  }
}

### NOX RBPome
rawFile <- "/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/data/raw/MQtoSQ.csv"
doneFile <- "/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/data/old/P445_AG_PolyA_40_all_3_cond_20220611_Final/P445_AG_PolyA_40_20220617_gmin_UV_nonorm_p2_NOX_Final/P445_AG_PolyA_40_20220617_gmin_UV_nonorm_p2_NOX_PROTEIN.tsv"
compareRawNDone(doneFile,rawFile)

doneFile <- "/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/data/new/PolyA_RIC_nonorm_20220928/P445_AG_PolyA_40_20220617_gmin_UV_nonorm_p2_NOX_Final/P445_AG_PolyA_40_20220617_gmin_UV_nonorm_p2_NOX_PROTEIN.tsv"
compareRawNDone(doneFile,rawFile)

### UVX RBPome
rawFile <- "/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/data/raw/MQtoSQ_UV.csv"
doneFile <- "/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/data/old/P445_AG_PolyA_40_all_3_cond_20220611_Final/P445_AG_PolyA_40_20220617_gmin_UV_p2_no511_Final/P445_AG_PolyA_40_20220617_gmin_UV_p2_no511_PROTEIN.tsv"
compareRawNDone(doneFile,rawFile)

doneFile <- "/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/data/new/PolyA_RIC_nonorm_20220928/P445_AG_PolyA_40_20220617_gmin_UV_p2_nonorm_no511_Final/P445_AG_PolyA_40_20220617_gmin_UV_p2_nonorm_no511_Final_PROTEIN.tsv"
compareRawNDone(doneFile,rawFile)
### FAX RBPome
rawFile <- "/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/data/raw/MQtoSQ_FAX.csv"
doneFile <- "/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/data/old/P445_AG_PolyA_40_all_3_cond_20220611_Final/P445_AG_PolyA_40_20220617_gmin_FAX_p2_Final/P445_AG_PolyA_40_20220617_gmin_FAX_p2_PROTEIN.tsv"
compareRawNDone(doneFile,rawFile)

doneFile <- "/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/data/new/PolyA_RIC_nonorm_20220928/P445_AG_PolyA_40_20220617_gmin_FAX_p2_nonorm_no12_Final/P445_AG_PolyA_40_20220617_gmin_FAX_p2_nonorm_no12_Final_PROTEIN.tsv"
compareRawNDone(doneFile,rawFile)

### Proteome
rawFile <- "/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/data/raw/SNtoSQ-Prot40_SafeQuant_Input.csv"

### NOX Proteome
crosslink <- "NOX"

doneFile <- "/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/data/new/P445_Prot_DIA_iBAQ_nocomb_20220826_all_NOX_Final/P445_Prot_DIA_iBAQ_noOutliers_nocomb_20220826_NOX_Final/P445_Prot_DIA_iBAQ_noOutliers_nocomb_20220826_NOX_Final_PROTEIN.tsv"
compareRawNDone(doneFile,rawFile,crosslink)

doneFile <- "/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/data/new/P445_Prot_DIA_iBAQ_nocomb_20220826_all_NOX_Final/P445_Prot_DIA_iBAQ_nocomb_20220826_NOX_Final2/P445_Prot_DIA_iBAQ_nocomb_20220826_NOX_Final2_PROTEIN.tsv"
compareRawNDone(doneFile,rawFile,crosslink)

doneFile <- "/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/data/old/P445_Prot_DIA_iBAQ_noOutliers_nocomb_20220826_Final/P445_Prot_DIA_iBAQ_noOutliers_nocomb_20220826_NOX_Final/P445_Prot_DIA_iBAQ_noOutliers_nocomb_20220826_NOX_Final_PROTEIN.tsv"
compareRawNDone(doneFile,rawFile,crosslink)

### UVX Proteome
crosslink <- "UV"
doneFile <- "/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/data/old/P445_Prot_DIA_iBAQ_noOutliers_nocomb_20220826_Final/P445_Prot_DIA_iBAQ_noOutliers_nocomb_20220826_UV_Final/P445_Prot_DIA_iBAQ_noOutliers_nocomb_20220826_UV_Final_PROTEIN.tsv"
compareRawNDone(doneFile,rawFile,crosslink)

### FAX Proteome
crosslink <- "FAX"
doneFile <- "/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/data/old/P445_Prot_DIA_iBAQ_noOutliers_nocomb_20220826_Final/P445_Prot_DIA_iBAQ_noOutliers_nocomb_20220826_FAX_Final/P445_Prot_DIA_iBAQ_noOutliers_nocomb_20220826_FAX_Final_PROTEIN.tsv"
compareRawNDone(doneFile,rawFile,crosslink)


