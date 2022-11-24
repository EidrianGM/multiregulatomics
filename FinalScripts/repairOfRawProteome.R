
proteomeFile <- "/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/data/raw/SNtoSQ-Prot40_SafeQuant_Input.csv"
NOXRBPomeFile <- "/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/data/raw/MQtoSQ.csv"
FAXRBPomeFile <- "/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/data/raw/MQtoSQ_FAX.csv"
UVXRBPomeFile <- "/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/data/raw/MQtoSQ_UV.csv"

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
# properMapping <- merge(proteomeDF[,keyCols],rbpomeInfo,by="Sequence", all.x = T)
# properMapping$Accession.y[is.na(properMapping$Accession.y)] <- properMapping$Accession.x[is.na(properMapping$Accession.y)]
# properMapping$Description.y[is.na(properMapping$Description.y)] <- properMapping$Description.x[is.na(properMapping$Description.y)]
# properMapping$Accession.x <- properMapping$Accession.y
# properMapping$Description.x <- properMapping$Description.y

newProteomeDF <- merge(proteomeDF,rbpomeInfo,by="Sequence", all.x = T)
newProteomeDF$Accession.y[is.na(newProteomeDF$Accession.y)] <- newProteomeDF$Accession.x[is.na(newProteomeDF$Accession.y)]
newProteomeDF$Description.y[is.na(newProteomeDF$Description.y)] <- newProteomeDF$Description.x[is.na(newProteomeDF$Description.y)]
newProteomeDF$Accession.x <- newProteomeDF$Accession.y
newProteomeDF$Description.x <- newProteomeDF$Description.y
newProteomeDF$Accession.y <- NULL
newProteomeDF$Description.y <- NULL

colnames(newProteomeDF) <- gsub(".x$","",colnames(newProteomeDF))
newProteomeDF <- newProteomeDF[colnames(proteomeDF)]

oldproteomeDF <- read.delim(proteomeFile,sep = ",",header = F)

outFile <- "/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/data/raw/newProteome.csv"
write.table(oldproteomeDF[1:3,],outFile,sep = ",",append = F,quote = T,row.names = F,col.names = F)
write.table(newProteomeDF,outFile,sep = ",",append = T,quote = T,row.names = F,col.names = F)

