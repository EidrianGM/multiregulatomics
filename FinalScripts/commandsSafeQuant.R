### Install SafeQuant
# Download the latest release https://github.com/eahrne/SafeQuant
# Use is locally 

commandExec <- "/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/scripts/SafeQuant-2.3.4/exec/safeQuant.R"
# In the commandExec file change the sourceDirOSX to "/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/scripts/SafeQuant-2.3.4/R/"

# Commands Instructions
# /home/eidriangm/Desktop/toDo/surrey/multiregulatomics/FinalScripts/commandsGuide.txt

###### Analysis of RBPome UVX 
## removing samples: 

rawRBPomeUVXfile <- "/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/data/raw/MQtoSQ_UV.csv"
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

rawRBPomeFAXfile <- "/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/data/raw/MQtoSQ_FAX.csv"
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
# ProteomeUVXfile <- "/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/data/raw/SNtoSQ-Prot40_SafeQuant_Input.csv"
Proteomefile <- "/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/data/raw/newProteome.csv"
ProteomeSamplesInfo <- read.delim(Proteomefile,nrows = 3, sep = ",")
myDF <- ProteomeSamplesInfo

normDCol <- grep("Norm",colnames(myDF)) # Raw Normalized Values are not the ones!
rawDCol <- grep("Raw",colnames(myDF))
expresionMatrix <- myDF[normDCol:(rawDCol-1)] 
expresionMatrix[2,] <- gsub("_.*","",expresionMatrix[2,])
phenData <- as.data.frame(t(expresionMatrix[1:2,])); colnames(phenData) <- c("class","id"); 
phenData$class <- gsub(".*\\s","",phenData$class); conditionCol <- "class"
rownames(phenData) <- 1:nrow(phenData)

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
outfilesLabel <- "ProteomeUVX_np2_norm_gmin_no4-6-23-31-32"
cmd <- paste("Rscript",commandExec,"-i",Proteomefile,"-l",outfilesLabel,"--EX ",expDesignTag,"--SM knn")
cat(cmd)
system(cmd)

# FAX 
crosslink <- "FAX"
phenData$class[grepl(paste0(crosslink,"no"),phenData$class)] <- paste0(crosslink,"no")

sammples2remove <- c(15,25,26,36,39) 
phenData$filter <- as.numeric(phenData$id) %in% sammples2remove | !(grepl(crosslink,phenData$class))
subPhenData <- phenData[!phenData$filter,]

subPhenData <- phenData[!phenData$filter,]
expDesignTag <- c()
for (class in unique(subPhenData$class)){
  expDesignTag <- c(expDesignTag,paste(rownames(subPhenData)[subPhenData$class == class],collapse=","))
}
expDesignTag <- paste0('"',paste(expDesignTag,collapse = ":"),'"')
outfilesLabel <- "ProteomeFAX_np2_norm_gmin_no15-25-26-36-39"
cmd <- paste("Rscript",commandExec,"-i",Proteomefile,"-l",outfilesLabel,"--EX ",expDesignTag,"--SM knn")
cat(cmd)
system(cmd)


###### BACKGROUND REMOVAL
## removing samples: 

rawRBPomeNOXfile <- "/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/FinalData/Raw/MQtoSQ.csv"
rawRBPomeNOXSamplesInfo <- read.delim(rawRBPomeNOXfile,nrows = 3, sep = ",")
myDF <- rawRBPomeNOXSamplesInfo

normDCol <- grep("Norm",colnames(myDF)) # Raw Normalized Values are not the ones!
rawDCol <- grep("Raw",colnames(myDF))
expresionMatrix <- myDF[normDCol:(rawDCol-1)] 
expresionMatrix[2,] <- gsub("_.*","",expresionMatrix[2,])
subPhenData <- as.data.frame(t(expresionMatrix[1:2,])); colnames(subPhenData) <- c("class","id"); 
subPhenData$class <- gsub(".*\\s","",subPhenData$class); conditionCol <- "class"
subPhenData$class <- gsub("no.*|with.*","",subPhenData$class)
rownames(subPhenData) <- 1:nrow(subPhenData)
expDesignTag <- c()
for (class in c("NOX","FAX","UV")){
  expDesignTag <- c(expDesignTag,paste(rownames(subPhenData)[subPhenData$class == class],collapse=","))
}
expDesignTag <- paste0('"',paste(expDesignTag,collapse = ":"),'"')
outfilesLabel <- "RBPomeNOX_np2_norm"
cmd <- paste("Rscript",commandExec,"-i",rawRBPomeNOXfile,"-l",outfilesLabel,"--EX ",expDesignTag," --FN 2 --SM no")
cat(cmd)
system(cmd)

