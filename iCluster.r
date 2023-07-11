library(iClusterPlus)
library(openxlsx)

### Read Data
FAXacceptedProts <- read.delim("FinalData/BackgroundRemoval/FAXAccBackGrAcceptedFC3.tsv",header = F)[,1]
UVXacceptedProts <- read.delim("FinalData/BackgroundRemoval/UVXAccBackGrAcceptedFC3.tsv",header = F)[,1]

wholeDFfile <- "FinalData/allDataDF.tsv"
wholeDF <- read.delim(wholeDFfile,quote = "")
wholeDF$UniprotName <- gsub("_YEAST","",wholeDF$UniprotName)

samplesNames <- colnames(wholeDF)[grep('.X0',colnames(wholeDF))]
runids <- as.numeric(sapply(strsplit(samplesNames,'X0|_'), function(x) return(x[2])))
datatypes <- sapply(strsplit(samplesNames,'ome'), function(x) return(x[1]))
datatypes <- paste0(gsub('RBP','mRBP',datatypes),'ome')
crosslinks <- paste0(sapply(strsplit(samplesNames,'ome|X'), function(x) return(x[2])),'X')

samplesNamesPhenoD <- as.data.frame(cbind(runid=runids,samplesName=samplesNames,datatype=datatypes,
                                          crosslink=crosslinks,uniqueid=paste0(datatypes,crosslinks,runids)))

phenoData <- read.xlsx('data/phenoData.xlsx')
phenoData <- phenoData[!(phenoData$Crosslinking == 'NOX'),]
phenoData <- phenoData[c('Data','Treatment','Crosslinking','Run.number','Sample.ID')]
phenoData$Data <- gsub(' ','',phenoData$Data)
phenoData$Treatment <- gsub(' ','',phenoData$Treatment)
phenoData$Crosslinking <- gsub(' ','',phenoData$Crosslinking)
phenoData <- phenoData[!duplicated(phenoData),]
phenoData$uniqueid <- paste0(phenoData$Data,phenoData$Crosslinking,phenoData$Run.number)

totalInfoPhenoD <- merge(phenoData,samplesNamesPhenoD,by='uniqueid')
nrow(samplesNamesPhenoD) == nrow(totalInfoPhenoD)
totalInfoPhenoD <- totalInfoPhenoD[c('samplesName',colnames(phenoData))]

degsCols <- colnames(wholeDF)[grep('brep',colnames(wholeDF))]
repnum <- sapply(strsplit(degsCols,'\\.|brep'), function(x) return(x[4]))
degsTreat <- sapply(strsplit(degsCols,'\\.'), function(x) return(x[2]))
degsPhenoD <- cbind(samplesName=degsCols,Data='Transcriptome',Treatment=degsTreat,Crosslinking='-',Run.number='-',Sample.ID=repnum,uniqueid='-')

totalInfoPhenoD <- rbind(degsPhenoD,totalInfoPhenoD)
write.table(totalInfoPhenoD,file = 'FinalData/colNamesPhenoData.tsv',sep = '\t',quote = F,row.names = F)
View(totalInfoPhenoD)

######## PHENODATA TO iCLUSTER expected
iClusterPhenoD <- totalInfoPhenoD
iClusterPhenoD$Sample.ID <- gsub('.*\\(|\\)','',iClusterPhenoD$Sample.ID)

#totalInfoPhenoD$Sample.ID[totalInfoPhenoD$Sample.ID == "UV control (without any treatment) (3)"] == c("UV control (without any treatment) (3)", "UV control (without any treatment) (3)", "UV control (without any treatment) (4)", "UV control (without any treatment) (3)")

crosslink <- 'FAX'
treatment <- 'SA'
if(crosslink == 'FAX'){
  subDF <- wholeDF[which(wholeDF$proteinName %in% FAXacceptedProts),] 
}else{
  subDF <- wholeDF[which(wholeDF$proteinName %in% UVXacceptedProts),] 
}
samplesSelection <- iClusterPhenoD$Treatment == treatment & 
                    (iClusterPhenoD$Crosslinking == crosslink | iClusterPhenoD$Data == 'Transcriptome')

subDF <- subDF[c('DEGs.target','UniprotName',iClusterPhenoD$samplesName[samplesSelection])]
subDF <- subDF[complete.cases(subDF),]
View(subDF)

dfsSelected <- list()
for (datatype in unique(iClusterPhenoD$Data)){
  datatypeSelection <- iClusterPhenoD$Data == datatype & samplesSelection
  mydatatypeDF <- subDF[iClusterPhenoD$samplesName[datatypeSelection]]
  #mydatatypeDF <- subDF[c('UniprotName',iClusterPhenoD$samplesName[datatypeSelection])]
  if (datatype == 'Transcriptome'){
    mydatatypeDF <- aggregate(mydatatypeDF,by = list(subDF$UniprotName),FUN = sum)
  }else{
    mydatatypeDF <- aggregate(mydatatypeDF,by = list(subDF$UniprotName),FUN = mean)
  }
  mycolNames <- paste(crosslink,treatment,iClusterPhenoD$Sample.ID[datatypeSelection],sep = '_')
  row.names(mydatatypeDF) <- mydatatypeDF$Group.1
  mydatatypeDF$Group.1 <- NULL
  colnames(mydatatypeDF) <- mycolNames
  dfsSelected[[datatype]] <- as.matrix(mydatatypeDF)
}
names(dfsSelected)

identical(colnames(dfsSelected[[1]]),colnames(dfsSelected[[2]]))
identical(colnames(dfsSelected[[1]]),colnames(dfsSelected[[3]]))

identical(rownames(dfsSelected[[1]]),rownames(dfsSelected[[2]]))
identical(rownames(dfsSelected[[1]]),rownames(dfsSelected[[3]]))

df1 <- as.matrix(dfsSelected[[1]])
df2 <- as.matrix(dfsSelected[[2]])
df3 <- as.matrix(dfsSelected[[3]])

bayfit <- tune.iClusterBayes(cpus=6,dt1=df1,dt2=df2,dt3=df3,
                            type=c("gaussian","gaussian","gaussian"),K=4:9,n.burnin=1000,
                            n.draw=1200,prior.gamma=c(0.5,0.5,0.5),sdev=0.05,thin=3)
allBIC = NULL
devratio = NULL
nK = length(bayfit$fit)
for(i in 1:nK){
  allBIC = c(allBIC,bayfit$fit[[i]]$BIC)
  devratio = c(devratio,bayfit$fit[[i]]$dev.ratio)
}

par(mar=c(4.0,4.0,0.5,0.5),mfrow=c(1,2))
plot(1:nK, allBIC,type="b",xlab="k",ylab="BIC",pch=c(1,1,19,1,1,1))
plot(1:nK,devratio,type="b",xlab="k",ylab="Deviance ratio",pch=c(1,1,19,1,1,1))

row.names(mydatatypeDF)

View(subDF[subDF$UniprotName %in% subDF$UniprotName[duplicated(subDF$UniprotName)],])

View(mydatatypeDF[mydatatypeDF$UniprotName %in% mydatatypeDF$UniprotName[duplicated(mydatatypeDF$UniprotName)],])

any(duplicated(mydatatypeDF))


colnames(selectDF) <- c('UniprotName',iClusterPhenoD$Sample.ID[samplesSelection])

colnames(selectDF)


iClusterPhenoD$Sample.ID[samplesSelection]

View(iClusterPhenoD)

proteome_InfoPhenoD <- totalInfoPhenoD[totalInfoPhenoD$Data == 'mRBPome',]
rbpome_InfoPhenoD <- totalInfoPhenoD[totalInfoPhenoD$Data == 'Proteome',]

table(totalInfoPhenoD$Sample.ID,totalInfoPhenoD$Treatment)
table(totalInfoPhenoD$Sample.ID)


outdir <- "/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/FinalData/HeatMaps"

