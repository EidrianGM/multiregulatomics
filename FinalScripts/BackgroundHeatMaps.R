

backgroundRBPomeFile <- "/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/FinalData/BackgroundRemoval/backgroundRBPomeDF.tsv"
backgroundRBPomeDF <- read.delim(backgroundRBPomeFile)

FCcutOff <- 3
FAXnoise <- which(abs(backgroundRBPomeDF$log2ratiosFAX) > FCcutOff); length(FAXnoise)
UVXnoise <- which(abs(backgroundRBPomeDF$log2ratiosUVX) > FCcutOff); length(UVXnoise)

FCcutOff <- 2
FAXnoise <- which(abs(backgroundRBPomeDF$log2ratiosFAX) > FCcutOff); length(FAXnoise)
UVXnoise <-which(abs(backgroundRBPomeDF$log2ratiosUVX) > FCcutOff); length(UVXnoise)

FCcutOff <- 1
FAXnoise <- which(abs(backgroundRBPomeDF$log2ratiosFAX) > FCcutOff); length(FAXnoise)
UVXnoise <-which(abs(backgroundRBPomeDF$log2ratiosUVX) > FCcutOff); length(UVXnoise)

