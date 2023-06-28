### Install SafeQuant
# Download the latest release https://github.com/eahrne/SafeQuant
# Use is locally 

commandExec <- "scripts/SafeQuant-2.3.4/exec/safeQuant.R"
# In the commandExec file change the sourceDirOSX to "scripts/SafeQuant-2.3.4/R/"

# Commands Instructions
# FinalScripts/commandsGuide.txt

###### Analysis of RBPome UVX 
## removing samples: 

rawRBPomeUVXfile <- "FinalData/Raw/correctedRAW/UVXmRBPomeRAW.csv"
rawRBPomeUVXSamplesInfo <- read.delim(rawRBPomeUVXfile,nrows = 3, sep = ",")
myDF <- rawRBPomeUVXSamplesInfo

normDCol <- grep("Norm",colnames(myDF)) # Raw Normalized Values are not the ones!
rawDCol <- grep("Raw",colnames(myDF))
expresionMatrix <- myDF[normDCol:(rawDCol-1)] 
expresionMatrix[2,] <- gsub("_.*","",expresionMatrix[2,])
subPhenData <- as.data.frame(t(expresionMatrix[1:2,])); colnames(subPhenData) <- c("class","id"); 
subPhenData$class <- gsub(".*\\s","",subPhenData$class); conditionCol <- "class"

sammples2remove <- c(10,7,16,3) 
subPhenData$filter <- as.numeric(subPhenData$id) %in% sammples2remove
rownames(subPhenData) <- 1:nrow(subPhenData)
expDesignTag <- c()
for (class in unique(subPhenData$class)){
  expDesignTag <- c(expDesignTag,paste(rownames(subPhenData)[subPhenData$class == class & !subPhenData$filter],collapse=","))
}
expDesignTag <- paste0('"',paste(expDesignTag,collapse = ":"),'"')
outfilesLabel <- "RBPomeUVX_np2_norm_gmin_no10-7-16-3"
cmd <- paste("Rscript",commandExec,"-i",rawRBPomeUVXfile,"-l",outfilesLabel,"--EX ",expDesignTag," --FN 2 --SM gmin")
cat(cmd)
system(cmd)

###### Analysis of RBPome FAX 
## removing samples: 

rawRBPomeFAXfile <- "FinalData/Raw/correctedRAW/FAXmRBPomeRAW.csv"
rawRBPomeFAXSamplesInfo <- read.delim(rawRBPomeFAXfile,nrows = 3, sep = ",")
myDF <- rawRBPomeFAXSamplesInfo

normDCol <- grep("Norm",colnames(myDF)) # Raw Normalized Values are not the ones!
rawDCol <- grep("Raw",colnames(myDF))
expresionMatrix <- myDF[normDCol:(rawDCol-1)] 
expresionMatrix[2,] <- gsub("_.*","",expresionMatrix[2,])
subPhenData <- as.data.frame(t(expresionMatrix[1:2,])); colnames(subPhenData) <- c("class","id"); 
subPhenData$class <- gsub(".*\\s","",subPhenData$class); conditionCol <- "class"

sammples2remove <- c(12) 
subPhenData$filter <- as.numeric(subPhenData$id) %in% sammples2remove
rownames(subPhenData) <- 1:nrow(subPhenData)
expDesignTag <- c()
for (class in unique(subPhenData$class)){
  expDesignTag <- c(expDesignTag,paste(rownames(subPhenData)[subPhenData$class == class & !subPhenData$filter],collapse=","))
}
expDesignTag <- paste0('"',paste(expDesignTag,collapse = ":"),'"')
outfilesLabel <- "RBPomeFAX_np2_norm_gmin_no12"
cmd <- paste("Rscript",commandExec,"-i",rawRBPomeFAXfile,"-l",outfilesLabel,"--EX ",expDesignTag," --FN 2 --SM gmin")
cat(cmd)
system(cmd)

###### Analysis of Proteome  
Proteomefile <- "FinalData/Raw/correctedRAW/newProteome.csv"
ProteomeSamplesInfo <- read.delim(Proteomefile,nrows = 3, sep = ",")
myDF <- ProteomeSamplesInfo

normDCol <- grep("Norm",colnames(myDF)) # Raw Normalized Values are not the ones!
rawDCol <- grep("Raw",colnames(myDF))
expresionMatrix <- myDF[normDCol:(rawDCol-1)] 
expresionMatrix[2,] <- gsub("_.*","",expresionMatrix[2,])
phenData <- as.data.frame(t(expresionMatrix[1:2,])); colnames(phenData) <- c("class","id"); 
phenData$class <- gsub(".*\\s","",phenData$class); conditionCol <- "class"
rownames(phenData) <- 1:nrow(phenData)

# NOX
crosslink <- "NOX"
phenData$class[grepl(paste0(crosslink,"no"),phenData$class)] <- paste0(crosslink,"no")
phenData$filter <- as.numeric(phenData$id) %in% sammples2remove | !(grepl(crosslink,phenData$class))
subPhenData <- phenData[!phenData$filter,]
expDesignTag <- c()
for (class in unique(subPhenData$class)){
  expDesignTag <- c(expDesignTag,paste(rownames(subPhenData)[subPhenData$class == class],collapse=","))
}
expDesignTag <- paste0('"',paste(expDesignTag,collapse = ":"),'"')
outfilesLabel <- "ProteomeNOX_np1_norm_knngmin"
cmd <- paste("Rscript",commandExec,"-i",Proteomefile,"-l",outfilesLabel,"--EX ",expDesignTag,"--SM knn")
cat(cmd)
system(cmd)

# UVX
crosslink <- "UV"
phenData$class[grepl(paste0(crosslink,"no"),phenData$class)] <- paste0(crosslink,"no")

sammples2remove <- c(4,6,31,32,23) 
phenData$filter <- as.numeric(phenData$id) %in% sammples2remove | !(grepl(crosslink,phenData$class))
subPhenData <- phenData[!phenData$filter,]
expDesignTag <- c()
for (class in unique(subPhenData$class)){
  expDesignTag <- c(expDesignTag,paste(rownames(subPhenData)[subPhenData$class == class],collapse=","))
}
expDesignTag <- paste0('"',paste(expDesignTag,collapse = ":"),'"')
outfilesLabel <- "ProteomeUVX_np1_norm_knngmin_no4-6-23-31-32"
cmd <- paste("Rscript",commandExec,"-i",Proteomefile,"-l",outfilesLabel,"--EX ",expDesignTag,"--SM knn")
cat(cmd)
system(cmd)

# FAX 
crosslink <- "FAX"
phenData$class[grepl(paste0(crosslink,"no"),phenData$class)] <- paste0(crosslink,"no")

sammples2remove <- c(15,25,26,36,39) 
phenData$filter <- as.numeric(phenData$id) %in% sammples2remove | !(grepl(crosslink,phenData$class))
subPhenData <- phenData[!phenData$filter,]

expDesignTag <- c()
for (class in unique(subPhenData$class)){
  expDesignTag <- c(expDesignTag,paste(rownames(subPhenData)[subPhenData$class == class],collapse=","))
}
expDesignTag <- paste0('"',paste(expDesignTag,collapse = ":"),'"')
outfilesLabel <- "ProteomeFAX_np1_norm_knngmin_no15-25-26-36-39"
cmd <- paste("Rscript",commandExec,"-i",Proteomefile,"-l",outfilesLabel,"--EX ",expDesignTag,"--SM knn")
cat(cmd)
system(cmd)

###### BACKGROUND REMOVAL
# Obtain log2FC of the proteins with NOX raw data == NO normalization,  NO imputation and NO samples removal
# SafeQuant script has been edited by Adrian Garcia Moreno to finish after log2FC calculation and to
# use the mean approach to calculate the log2FC (by default is the median)

# The calculation of background must be only done with the peptides available in FAX and UVX
# Since these are different, the initial NOX raw data must be divided in two
# NOX_FAX and # NOX_UVX

rawRBPomeNOXfile <- "FinalData/Raw/correctedRAW/NOXmRBPomeRAW.csv"
rawNOXdf <- read.delim(rawRBPomeNOXfile,skip = 2, sep = ",")
rawFAXdf <- read.delim("FinalData/Raw/correctedRAW/FAXmRBPomeRAW.csv",skip = 2, sep = ",")
rawUVXdf <- read.delim("FinalData/Raw/correctedRAW/UVXmRBPomeRAW.csv",skip = 2, sep = ",")



# ############### FILTERING OF NOX BY PEPTIDES ??? ------> Weird issue
# # UNCOMMON PEPTIDES REMOVAL
# rawNOXuniquePeptides <- paste0(rawNOXdf$Sequence,rawNOXdf$Modifications,rawNOXdf$Accession)
# rawFAXuniquePeptides <- paste0(rawFAXdf$Sequence,rawFAXdf$Modifications,rawFAXdf$Accession)
# rawUVXuniquePeptides <- paste0(rawUVXdf$Sequence,rawUVXdf$Modifications,rawUVXdf$Accession)
# 
# toVennList <- list(NOXpeps = rawNOXuniquePeptides, FAXpeps = rawFAXuniquePeptides, UVXpeps=rawUVXuniquePeptides)
# ggvennSA <- ggVennDiagram(toVennList, color = 2, lwd = 0.7) + scale_fill_gradient(low = "#F4FAFE", high = "#4981BF") + theme(legend.position = "none")
# nameOut <- "FinalData/BackgroundRemoval/vennRawRBPomesPeptides.png"
# ggsave(filename = nameOut, plot = ggvennSA, width = 30, height = 15, units = 'cm', dpi = 'print')
# 
# all(rawFAXuniquePeptides %in% rawNOXuniquePeptides) ## TRUE, All peptides of FAX are in NOX == Perfect
# all(rawUVXuniquePeptides %in% rawNOXuniquePeptides) ## TRUE, All peptides of UVX are in NOX == Perfect
# 
# View(rawNOXdf[duplicated(rawNOXuniquePeptides),])
# 
# length(rawNOXuniquePeptides); length(unique(rawNOXuniquePeptides))
# length(rawFAXuniquePeptides); length(unique(rawFAXuniquePeptides))
# length(rawUVXuniquePeptides); length(unique(rawUVXuniquePeptides))
# ## There are duplications becasue of columns like rawNOXdf$Deconvoluted.charges
# 
# ## The number of peptides to be selected in NOX should be the same available in FAX and UVX ---> BUT IT IS NOT
# sum(rawNOXuniquePeptides %in% rawFAXuniquePeptides) == length(rawFAXuniquePeptides) ## FALSE = some peptides in FAX are identified more than once in NOX !?
# sum(rawNOXuniquePeptides %in% rawUVXuniquePeptides) == length(rawUVXuniquePeptides) ## FALSE = "" but in UVX
# 
# rawNOXdf2FAX <- rawNOXdf[rawNOXuniquePeptides %in% rawFAXuniquePeptides,]
# rawNOXdf2UVX <- rawNOXdf[rawNOXuniquePeptides %in% rawUVXuniquePeptides,]
# 
# outFile <- "FinalData/Raw/correctedRAW/bypepsNOX_FAXmRBPomeRAW.csv"
# write.table(rawRBPomeNOXSamplesInfo,outFile,sep = ",",append = F,quote = T,row.names = F,col.names = F)
# write.table(rawNOXdf2FAX,outFile,sep = ",",append = T,quote = T,row.names = F,col.names = F)
# 
# outFile <- "FinalData/Raw/correctedRAW/bypepsNOX_UVXmRBPomeRAW.csv"
# write.table(rawRBPomeNOXSamplesInfo,outFile,sep = ",",append = F,quote = T,row.names = F,col.names = F)
# write.table(rawNOXdf2UVX,outFile,sep = ",",append = T,quote = T,row.names = F,col.names = F)
###############

################################ Filtering of NOX for background removal by protein level
# UNCOMMON PROTEINS REMOVAL

all(rawFAXdf$Accession %in% rawNOXdf$Accession) ## TRUE, All proteins of FAX are in NOX == Perfect
all(rawUVXdf$Accession %in% rawNOXdf$Accession) ## TRUE, All proteins of UVX are in NOX == Perfect

rawNOXdf2FAX <- rawNOXdf[rawNOXdf$Accession %in% rawFAXdf$Accession,]
rawNOXdf2UVX <- rawNOXdf[rawNOXdf$Accession %in% rawUVXdf$Accession,]

toVennList <- list(NOXprots = rawNOXdf$Accession, FAXprots = rawFAXdf$Accession, UVXprots=rawUVXdf$Accession)
ggvennSA <- ggVennDiagram(toVennList, color = 2, lwd = 0.7) + scale_fill_gradient(low = "#F4FAFE", high = "#4981BF") + theme(legend.position = "none")
nameOut <- "FinalData/BackgroundRemoval/vennRawRBPomesProteins.png"
ggsave(filename = nameOut, plot = ggvennSA, width = 30, height = 15, units = 'cm', dpi = 'print')

toVennList <- list(RBPomeUVX_no10_7 = RBPomeUVXOldProts, RBPomeUVXAdriProts_no10_7_16_3 = RBPomeUVXAdriProts)
ggvennSA <- ggVennDiagram(toVennList, color = 2, lwd = 0.7) + scale_fill_gradient(low = "#F4FAFE", high = "#4981BF") + theme(legend.position = "none")

nameOut <- "/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/vennUVXRBPomeOldvsNew.png"
ggsave(filename = nameOut, plot = ggvennSA, width = 30, height = 15, units = 'cm', dpi = 'print')

outFile <- "FinalData/Raw/correctedRAW/NOX_FAXmRBPomeRAW.csv"
write.table(rawRBPomeNOXSamplesInfo,outFile,sep = ",",append = F,quote = T,row.names = F,col.names = F)
write.table(rawNOXdf2FAX,outFile,sep = ",",append = T,quote = T,row.names = F,col.names = F)

outFile <- "FinalData/Raw/correctedRAW/NOX_UVXmRBPomeRAW.csv"
write.table(rawRBPomeNOXSamplesInfo,outFile,sep = ",",append = F,quote = T,row.names = F,col.names = F)
write.table(rawNOXdf2UVX,outFile,sep = ",",append = T,quote = T,row.names = F,col.names = F)

###################################
rawRBPomeNOXSamplesInfo <- read.delim(rawRBPomeNOXfile,nrows = 3, sep = ",",header = F)

rawRBPomeNOXSamplesInfoRow <- rawRBPomeNOXSamplesInfo[2,]
noxsamplesselected <- unique(as.numeric(gsub("_.*","",rawRBPomeNOXSamplesInfo[3,][grep("^0",rawRBPomeNOXSamplesInfo[3,])])))

rawRBPomeNOXSamplesInfo[3,][grep("^0",rawRBPomeNOXSamplesInfo[3,])]

infoData <- cbind(as.character(rawRBPomeNOXSamplesInfo[2,][grep("^0",rawRBPomeNOXSamplesInfo[3,])]),
           noxsamplesselected)
infoData <- infoData[!duplicated(infoData),]
infoData <- infoData[order(as.numeric(infoData[,2])),]
View(infoData)

setdiff(noxsamplesselected,samples2include)
setdiff(samples2include,noxsamplesselected)

FAXrbpomeFile <- "FinalData/Raw/correctedRAW/RBPomeFAX_np2_norm_gmin_no12/RBPomeFAX_np2_norm_gmin_no12_PROTEIN.tsv"
UVXrbpomeFile <- "FinalData/Raw/correctedRAW/RBPomeUVX_np2_norm_gmin_no10-7-16-3/RBPomeUVX_np2_norm_gmin_no10-7-16-3_PROTEIN.tsv"
FAXsamplesInfo <- read.delim(FAXrbpomeFile,nrows = 1,header = F)
FAXsamples <- gsub("_.*","",FAXsamplesInfo[grep("^0",FAXsamplesInfo)])
UVXsamplesInfo <- read.delim(UVXrbpomeFile,nrows = 1,header = F)
UVXsamples <- gsub("_.*","",UVXsamplesInfo[grep("^0",UVXsamplesInfo)])

FAXsamplesInfo <- read.delim(FAXrbpomeFile,nrows = 1,header = F)
FAXsamples <- gsub("_.*","",FAXsamplesInfo[grep("^0",FAXsamplesInfo)])

UVXsamplesInfo <- read.delim(UVXrbpomeFile,nrows = 1,header = F)
UVXsamples <- gsub("_.*","",UVXsamplesInfo[grep("^0",UVXsamplesInfo)])


samples2include <- as.numeric(c(FAXsamples,UVXsamples))
samples2remove <- c(1,6,2,9)

# #### USING THE UNCOMMON PEPTIDES REMOVAL IN NOX

# crosslinks <- c("FAX","UVX")
# for (crosslink in crosslinks){
#   rawRBPomeNOXfile <- paste0("FinalData/Raw/correctedRAW/bypepsNOX_",crosslink,"mRBPomeRAW.csv")
#   rawRBPomeNOXSamplesInfo <- read.delim(rawRBPomeNOXfile,nrows = 2, sep = ",",header = T)
#   
#   normDCol <- grep("Norm",colnames(rawRBPomeNOXSamplesInfo)) # Raw Normalized Values are not the ones!
#   rawDCol <- grep("Raw",colnames(rawRBPomeNOXSamplesInfo))
#   expresionMatrix <- rawRBPomeNOXSamplesInfo[normDCol:(rawDCol-1)] 
#   expresionMatrix[2,] <- gsub("_.*","",expresionMatrix[2,])
#   subPhenData <- as.data.frame(t(expresionMatrix[1:2,])); colnames(subPhenData) <- c("class","id"); 
#   subPhenData$class <- gsub(".*\\s","",subPhenData$class); conditionCol <- "class"
#   subPhenData$class <- gsub("UV","UVX",subPhenData$class)
#   rownames(subPhenData) <- 1:nrow(subPhenData)
#   #subPhenData <- subPhenData[grep("SA",subPhenData$class),]
#   subPhenData$class <- gsub("no.*|with.*","",subPhenData$class)
#   subPhenData$id <- as.numeric(subPhenData$id)
#   subPhenData <- subPhenData[(subPhenData$class == crosslink & subPhenData$id %in% samples2include) | (subPhenData$class == "NOX" & !(subPhenData$id %in% samples2remove)),]
#   expDesignTag <- c()
#   for (class in c("NOX",crosslink)){
#     expDesignTag <- c(expDesignTag,paste(rownames(subPhenData)[subPhenData$class == class],collapse=","))
#   }
#   expDesignTag <- paste0('"',paste(expDesignTag,collapse = ":"),'"')
#   ### MODIFY getRatios function in ExpressionAnalysis.R in SafeQuant to use mean or median
#   outfilesLabel <- paste0("byPEPmeanRBPomeNOX_",crosslink,"_np2_NOnorm")
#   cmd <- paste("Rscript scripts/SafeQuant-2.3.4/exec/safeQuant.R -i",rawRBPomeNOXfile,"-l",outfilesLabel,"--EX ",expDesignTag," --FN 2 --SM no --SR")
#   cat(cmd)
#   system(cmd)
# }

#### USING THE UNCOMMON PROTEINS REMOVAL IN NOX
crosslinks <- c("FAX","UVX")
crosslink <- crosslinks[2]
for (crosslink in crosslinks){
  rawRBPomeNOXfile <- paste0("FinalData/Raw/correctedRAW/NOX_",crosslink,"mRBPomeRAW.csv")
  rawRBPomeNOXSamplesInfo <- read.delim(rawRBPomeNOXfile,nrows = 2, sep = ",",header = T)
  
  normDCol <- grep("Norm",colnames(rawRBPomeNOXSamplesInfo)) # Raw Normalized Values are not the ones!
  rawDCol <- grep("Raw",colnames(rawRBPomeNOXSamplesInfo))
  expresionMatrix <- rawRBPomeNOXSamplesInfo[normDCol:(rawDCol-1)] 
  expresionMatrix[2,] <- gsub("_.*","",expresionMatrix[2,])
  subPhenData <- as.data.frame(t(expresionMatrix[1:2,])); colnames(subPhenData) <- c("class","id"); 
  subPhenData$class <- gsub(".*\\s","",subPhenData$class); conditionCol <- "class"
  subPhenData$class <- gsub("UV","UVX",subPhenData$class)
  rownames(subPhenData) <- 1:nrow(subPhenData)
  #subPhenData <- subPhenData[grep("SA",subPhenData$class),]
  subPhenData$class <- gsub("no.*|with.*","",subPhenData$class)
  subPhenData$id <- as.numeric(subPhenData$id)
  subPhenData <- subPhenData[(subPhenData$class == crosslink & subPhenData$id %in% samples2include) | (subPhenData$class == "NOX" & !(subPhenData$id %in% samples2remove)),]
  expDesignTag <- c()
  for (class in c("NOX",crosslink)){
    expDesignTag <- c(expDesignTag,paste(rownames(subPhenData)[subPhenData$class == class],collapse=","))
  }
  expDesignTag <- paste0('"',paste(expDesignTag,collapse = ":"),'"')
  ### MODIFY getRatios function in ExpressionAnalysis.R in SafeQuant to use mean or median
  outfilesLabel <- paste0("meanRBPomeNOX_",crosslink,"_np2_NOnorm")
  cmd <- paste("Rscript scripts/SafeQuant-2.3.4/exec/safeQuant.R -i",rawRBPomeNOXfile,"-l",outfilesLabel,"--EX ",expDesignTag," --FN 2 --SM no --SR")
  cat(cmd)
  system(cmd)
}

