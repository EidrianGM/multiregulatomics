library(ggplot2)

wholeDFFile <- "/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/FinalData/allDataDF.tsv"
wholeDF <- read.delim(wholeDFFile,quote = "")

FAXSAsig <- which(wholeDF$RBPomeFAX.qValue_PolyARNAFAXwithSA < 0.05)
UVXSAsig <- which(wholeDF$RBPomeUVX.qValue_PolyARNAUVwithSA < 0.05)

FAXtranscSAfc <- wholeDF$DEGs.log2FC.SA.YPD[FAXSAsig]
FAXproteoSAfc <- wholeDF$ProteomeFAX.log2ratio_FAXwithSA[FAXSAsig]
FAXrbpomeSAfc <- wholeDF$RBPomeFAX.log2ratio_PolyARNAFAXwithSA[FAXSAsig]

FAXSAfcs <- as.data.frame(rbind(cbind(FAXtranscSAfc,"Trancriptome"),cbind(FAXproteoSAfc,"Proteome"),cbind(FAXrbpomeSAfc,"mRBPome")))
colnames(FAXSAfcs) <- c("log2FCs","Datatype")

#FAXSAfcs <- cbind(FAXtranscSAfc,FAXproteoSAfc,FAXrbpomeSAfc)

FAXSAfcs$log2FCs
FAXSAfcs$protidx <- 1:nrow(FAXSAfcs)

myplot <- ggplot(masterEAdf2plot_clusterGenes, aes(x = datatype, y = pathway, alpha = padj, size = absNES, colour = NES)) + 
  geom_point() +
  scale_size(name = '|NES|') +
  scale_alpha(range = c(0.5,1), name = '-log10(FDR)') +
  scale_color_discrete(name = 'Datatype') +
  ylab('') +
  xlab('Data Type') +
  facet_grid(rows=vars(clusterGenes), scales = "free_y", space = "free_y") +
  theme(axis.text.y = element_text(size = 8, face = 'bold'),
        axis.text.x = element_text(size = 8, face = 'bold'),
        strip.text = element_text(hjust = 0.5, size = 8, face = 'bold'),
        strip.text.y.right = element_text(angle = 0, size = 4))

FAXtranscSAfc
FAXproteoSAfc
FAXrbpomeSAfc

signif(cor(FAXtranscSAfc,FAXproteoSAfc,method = "pearson"),3)
signif(coef(lm(FAXtranscSAfc ~ FAXproteoSAfc))[1],3)
signif(coef(lm(FAXtranscSAfc ~ FAXproteoSAfc))[2],3)

signif(cor(FAXtranscSAfc,FAXrbpomeSAfc,method = "pearson"),3)
signif(coef(lm(RBPomeUVXAdriProtsFCs ~ RBPomeUVXOldProtsFCs))[1],3)
signif(coef(lm(RBPomeUVXAdriProtsFCs ~ RBPomeUVXOldProtsFCs))[2],3)

signif(cor(FAXproteoSAfc,FAXrbpomeSAfc,method = "pearson"),3)
signif(coef(lm(RBPomeUVXAdriProtsFCs ~ RBPomeUVXOldProtsFCs))[1],3)
signif(coef(lm(RBPomeUVXAdriProtsFCs ~ RBPomeUVXOldProtsFCs))[2],3)


corScatterSA <- ggplot(data = FAXSAfcs, mapping = aes(x = protidx, y = log2FCs, colour=Datatype)) +
  geom_point(shape = 21, fill = '#0f993d', color = 'white', size = 3) +
  geom_smooth(method = "lm",se = F) +
  scale_color_discrete(name = 'Datatype') +
  annotate("text", x = -1, y = 4, col = "black", size = 3,
           label = paste("Pearson r = ", signif(cor(RBPomeUVXOldProtsFCs,RBPomeUVXAdriProtsFCs,method = "pearson"),3))) +
  annotate("text", x = 2.5, y = -2,  col = "red", size = 3,
           label = paste("y = ", signif(coef(lm(RBPomeUVXAdriProtsFCs ~ RBPomeUVXOldProtsFCs))[1],3), 
                         signif(coef(lm(RBPomeUVXAdriProtsFCs ~ RBPomeUVXOldProtsFCs))[2],2), "x"))

nameOut <- "/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/corrScaterUVXRBPomeOldvsNew.png"
ggsave(filename = nameOut, plot = corScatterSA, width = 30, height = 15, units = 'cm', dpi = 'print')