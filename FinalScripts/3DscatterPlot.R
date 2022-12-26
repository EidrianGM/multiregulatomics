

wholeDFFile <- "/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/FinalData/allDataDF.tsv"
wholeDF <- read.delim(wholeDFFile,quote = "")

FAXSAsig <- which(wholeDF$RBPomeFAX.qValue_PolyARNAFAXwithSA < 0.05)
UVXSAsig <- which(wholeDF$RBPomeUVX.qValue_PolyARNAUVwithSA < 0.05)
length(unique(c(FAXSAsig,UVXSAsig)))

length(FAXSAsig)

FAXtranscSAfc <- wholeDF$DEGs.log2FC.SA.YPD[FAXSAsig]
FAXproteoSAfc <- wholeDF$ProteomeFAX.log2ratio_FAXwithSA[FAXSAsig]
FAXrbpomeSAfc <- wholeDF$RBPomeFAX.log2ratio_PolyARNAFAXwithSA[FAXSAsig]

#FAXSAfcs <- rbind(cbind(FAXtranscSAfc,"Trancriptome"),cbind(FAXproteoSAfc,"Proteome"),cbind(FAXrbpomeSAfc,"mRBPome"))
#colnames(FAXSAfcs) <- c("log2FCs","Datatype")

FAXSAfcs <- cbind(FAXtranscSAfc,FAXproteoSAfc,FAXrbpomeSAfc)


library("scatterplot3d") # load
data(iris)
head(iris)


colors <- c("#999999", "#E69F00", "#56B4E9")
shapes <- c(16, 17, 18) 
outdir <- "/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/FinalData"
tiff(file.path(outdir,paste0('FAX_FDR05SA_FCs3D.tiff')), 2000, 2000, pointsize=12,res = 300)
scatterplot3d(FAXSAfcs, pch = shapes[1], color=colors[1],
              xlab = "Transcriptome",
              ylab = "Proteome",
              zlab = "mRBPome",
              grid=TRUE, box=FALSE)
dev.off()


# legend("bottom", legend = levels(FAXSAfcs$Datatype)[1],
#        col =  colors, 
#        pch = shapes, 
#        inset = -0.25, xpd = TRUE, horiz = TRUE)


wholeDF$DEGs.log2FC.SA.YPD
wholeDF$RBPomeUVX.log2ratio_PolyARNAUVwithSA
wholeDF$ProteomeUVX.log2ratio_Condition4
