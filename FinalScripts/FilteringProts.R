


### Remove Proteins With < 2 peptides measures and < 2 samples per condition
proteinUUIDs <- metadata[,c("Sequence", "Accession")]
proteinUUIDs <- proteinUUIDs[!duplicated(proteinUUIDs),]
table(proteinUUIDs)
# subExpDFsel <- which(metadata$Sequence == proteinUUIDs[1,1] & metadata$Accession == proteinUUIDs[1,2])

filteredPeptides <- c()
for (classSel in unique(subPhenData$class)){
  Accs2Remove <- c()
  for (acc2study in unique(proteinUUIDs$Accession)){
    subExpDFsel <- which(metadata$Accession == acc2study)
    samplesOfClass <- subPhenData$id[subPhenData$class == classSel]
    subExpDF <- expresionMatrix[subExpDFsel,samplesOfClass]
    peptidesAv <- nrow(subExpDF)
    if (peptidesAv == 1){
      # Protein has only one peptide sequence ==> Remove
      Accs2Remove <- c(Accs2Remove,acc2study)
    }
    samplesMeasures <- apply(subExpDF,2, FUN = function(x){sum(!is.na(x))})
    peptidesInSamples <- apply(subExpDF,1, FUN = function(x){sum(!is.na(x))})
    if (sum(samplesMeasures > 0) < 2){
      # Less that 2 samples with measures 
      Accs2Remove <- c(Accs2Remove,acc2study)
    }
    if (sum(peptidesInSamples > 0) <= 1){
      # Only one peptide have measures ==> Remove
      Accs2Remove <- c(Accs2Remove,acc2study)
    }
  }
  filteredPeptides <- rbind(filteredPeptides, cbind(classSel, Accs2Remove))
}
filteredPeptides <- as.data.frame(unique(filteredPeptides))

outfile <- "/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/FinalData/BackgroundRemoval/NOXmRBPomeBackGroundFilter.tsv"
write.table(filteredPeptides,file = outfile, sep = "\t", quote = F, row.names = F)

for (myclass in filteredPeptides$classSel){
  proteinsOfClass2remove <- filteredPeptides$Accs2Remove[filteredPeptides$classSel == myclass]
  samplesOfClass <- subPhenData$id[subPhenData$class == myclass]
  data2remove <- which(metadata$Accession %in% proteinsOfClass2remove) 
  expresionMatrix[data2remove,samplesOfClass] <- NA
}

minPepsPerProt <- 2
pepPerProt <- table(unique(data.frame(metadata$Accession,metadata$Sequence))[,1])
pepPerProt <- pepPerProt[metadata$Accession] 
protsWithMinPeps <- which(pepPerProt > minPepsPerProt)  

expresionMatrixFiltered <- expresionMatrix[(protsWithMinPeps | minPepsPerProt),]


