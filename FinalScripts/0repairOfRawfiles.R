
proteomeFile <- "data/raw/AlexRAW/SNtoSQ-Prot40_SafeQuant_Input.csv"
NOXRBPomeFile <- "data/raw/AlexRAW/MQtoSQ.csv"
FAXRBPomeFile <- "data/raw/AlexRAW/MQtoSQ_FAX.csv"
UVXRBPomeFile <- "data/raw/AlexRAW/MQtoSQ_UV.csv"

proteomeDF <- read.delim(proteomeFile,skip = 2,sep = ",")
NOXRBPomeDF <- read.delim(NOXRBPomeFile,skip = 2,sep = ",")
FAXRBPomeDF <- read.delim(FAXRBPomeFile,skip = 2,sep = ",")
UVXRBPomeDF <- read.delim(UVXRBPomeFile,skip = 2,sep = ",")

keyCols <- c("Sequence","Accession","Description")

rbpomeInfo <- rbind(NOXRBPomeDF[,keyCols],
                    FAXRBPomeDF[,keyCols],
                    UVXRBPomeDF[,keyCols])
rbpomeInfo <- rbpomeInfo[!duplicated(rbpomeInfo),]

# P05755 a paralogue that evidences the issue in raw data

# The mRBPome sequences that are paralogues are tagged with the proteins included 
# while the Proteome only notes one of them
# This makes us decide that the all the sequences common between Proteome an 
# mRBPome will be tagged following the mRBPome Accession column. 

newProteomeDF <- merge(proteomeDF,rbpomeInfo,by="Sequence", all.x = T)
newProteomeDF$Accession.y[is.na(newProteomeDF$Accession.y)] <- newProteomeDF$Accession.x[is.na(newProteomeDF$Accession.y)]
newProteomeDF$Description.y[is.na(newProteomeDF$Description.y)] <- newProteomeDF$Description.x[is.na(newProteomeDF$Description.y)]
newProteomeDF$Accession.x <- newProteomeDF$Accession.y
newProteomeDF$Description.x <- newProteomeDF$Description.y
newProteomeDF$Accession.y <- NULL
newProteomeDF$Description.y <- NULL

colnames(newProteomeDF) <- gsub(".x$","",colnames(newProteomeDF))
newProteomeDF <- newProteomeDF[colnames(proteomeDF)]

correctAccessionCol <- function(DF){
  UniprotsCoded <- c()
  for (badFormatedProtName in DF$Accession){
    badFormatedProtNameSep <- unlist(strsplit(badFormatedProtName,split = "\\|"))
    if (length(badFormatedProtNameSep) > 1){
      badFormatedProtNameSep <- paste(sort(badFormatedProtNameSep[seq(2,length(badFormatedProtNameSep),2)]),collapse = ";")
    }
    UniprotsCoded <- c(UniprotsCoded,badFormatedProtNameSep)
  }
  DF$Accession <- UniprotsCoded
  return(DF)
}

newProteomeDF <- correctAccessionCol(newProteomeDF)
oldproteomeDF <- read.delim(proteomeFile,nrows = 3,sep = ",",header = F)
outFile <- "FinalData/Raw/correctedRAW/newProteome.csv"
write.table(oldproteomeDF,outFile,sep = ",",append = F,quote = T,row.names = F,col.names = F)
write.table(newProteomeDF,outFile,sep = ",",append = T,quote = T,row.names = F,col.names = F)

# The column "Accession" in the mRBPomes are gonna be also corrected to avoid 
# further issues in SafeQuant peptides RollUP

newNOXRBPomeDF <- correctAccessionCol(NOXRBPomeDF)
newFAXRBPomeDF <- correctAccessionCol(FAXRBPomeDF)
newUVXRBPomeDF <- correctAccessionCol(UVXRBPomeDF)

oldNOXRBPomeDF <- read.delim(NOXRBPomeFile, sep = ",", header = F)
oldFAXRBPomeDF <- read.delim(FAXRBPomeFile, sep = ",", header = F)
oldUVXRBPomeDF <- read.delim(UVXRBPomeFile, sep = ",", header = F)

outFile <- "FinalData/Raw/correctedRAW/NOXmRBPomeRAW.csv"
write.table(oldNOXRBPomeDF[1:3,],outFile,sep = ",",append = F,quote = T,row.names = F,col.names = F)
write.table(newNOXRBPomeDF,outFile,sep = ",",append = T,quote = T,row.names = F,col.names = F)

outFile <- "FinalData/Raw/correctedRAW/FAXmRBPomeRAW.csv"
write.table(oldFAXRBPomeDF[1:3,],outFile,sep = ",",append = F,quote = T,row.names = F,col.names = F)
write.table(newFAXRBPomeDF,outFile,sep = ",",append = T,quote = T,row.names = F,col.names = F)

outFile <- "FinalData/Raw/correctedRAW/UVXmRBPomeRAW.csv"
write.table(oldUVXRBPomeDF[1:3,],outFile,sep = ",",append = F,quote = T,row.names = F,col.names = F)
write.table(newUVXRBPomeDF,outFile,sep = ",",append = T,quote = T,row.names = F,col.names = F)

##### Checking of Repairments

oldproteomeFile <- "data/raw/AlexRAW/SNtoSQ-Prot40_SafeQuant_Input.csv"
oldNOXRBPomeFile <- "data/raw/AlexRAW/MQtoSQ.csv"
oldFAXRBPomeFile <- "data/raw/AlexRAW/MQtoSQ_FAX.csv"
oldUVXRBPomeFile <- "data/raw/AlexRAW/MQtoSQ_UV.csv"
oldproteomeDF <- read.delim(oldproteomeFile,skip = 2,sep = ",")
oldNOXRBPomeDF <- read.delim(oldNOXRBPomeFile,skip = 2,sep = ",")
oldFAXRBPomeDF <- read.delim(oldFAXRBPomeFile,skip = 2,sep = ",")
oldUVXRBPomeDF <- read.delim(oldUVXRBPomeFile,skip = 2,sep = ",")


newproteomeFile <- "FinalData/Raw/correctedRAW/newProteome.csv"
newNOXRBPomeFile <- "FinalData/Raw/correctedRAW/NOXmRBPomeRAW.csv"
newFAXRBPomeFile <- "FinalData/Raw/correctedRAW/FAXmRBPomeRAW.csv"
newUVXRBPomeFile <- "FinalData/Raw/correctedRAW/UVXmRBPomeRAW.csv"
newproteomeDF <- read.delim(newproteomeFile,skip = 2,sep = ",")
newNOXRBPomeDF <- read.delim(newNOXRBPomeFile,skip = 2,sep = ",")
newFAXRBPomeDF <- read.delim(newFAXRBPomeFile,skip = 2,sep = ",")
newUVXRBPomeDF <- read.delim(newUVXRBPomeFile,skip = 2,sep = ",")


length(oldproteomeDF$Accession); length(unique(oldproteomeDF$Accession))
length(newproteomeDF$Accession); length(unique(newproteomeDF$Accession)) # Many paralogues rescued by the sequence merge

length(oldNOXRBPomeDF$Accession); length(unique(oldNOXRBPomeDF$Accession))
length(newNOXRBPomeDF$Accession); length(unique(newNOXRBPomeDF$Accession)) # ok

length(oldFAXRBPomeDF$Accession); length(unique(oldFAXRBPomeDF$Accession))
length(newFAXRBPomeDF$Accession); length(unique(newFAXRBPomeDF$Accession)) # ok

length(oldUVXRBPomeDF$Accession); length(unique(oldUVXRBPomeDF$Accession))
length(newUVXRBPomeDF$Accession); length(unique(newUVXRBPomeDF$Accession)) # ok

