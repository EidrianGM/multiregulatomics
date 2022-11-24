library(ggVennDiagram)
library(ggplot2)
RBPomeUVXAdriFile <- "/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/data/raw/RBPomeUVX_np2_norm_gmin_no10-7-16-3/RBPomeUVX_np2_norm_gmin_no10-7-16-3_PROTEIN.tsv"
RBPomeUVXAdriDF <- read.delim(RBPomeUVXAdriFile,quote = '"',check.names = T)
colnames(RBPomeUVXAdriDF) <- gsub("\\.$","",gsub("X.","",colnames(RBPomeUVXAdriDF)))

RBPomeUVXOldFile <- "/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/data/old/P445_AG_PolyA_40_all_3_cond_20220611_Final/P445_AG_PolyA_40_20220617_gmin_UV_p2_no511_Final/P445_AG_PolyA_40_20220617_gmin_UV_p2_no511_PROTEIN.tsv"
RBPomeUVXOldDF <- read.delim(RBPomeUVXOldFile,quote = '"')

RBPomeUVXOldProts <- RBPomeUVXOldDF$proteinName[RBPomeUVXOldDF$qValue_PolyARNAUVwithSA < 0.05]
RBPomeUVXAdriProts <- RBPomeUVXAdriDF$proteinName[RBPomeUVXAdriDF$qValue_PolyARNAUVwithSA < 0.05]
toVennList <- list(RBPomeUVX_no10_7 = RBPomeUVXOldProts, RBPomeUVXAdriProts_no10_7_16_3 = RBPomeUVXAdriProts)

ggvennSA <- ggVennDiagram(toVennList, color = 2, lwd = 0.7) + 
  scale_fill_gradient(low = "#F4FAFE", high = "#4981BF") +
  theme(legend.position = "none")

nameOut <- "/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/vennUVXRBPomeOldvsNew.png"
ggsave(filename = nameOut, plot = ggvennSA, width = 30, height = 15, units = 'cm', dpi = 'print')


commonProts <- intersect(RBPomeUVXOldProts,RBPomeUVXAdriProts)
RBPomeUVXOldProtsFCs <- RBPomeUVXOldDF$log2ratio_PolyARNAUVwithSA[RBPomeUVXOldDF$proteinName %in% commonProts]
RBPomeUVXAdriProtsFCs <- RBPomeUVXAdriDF$log2ratio_PolyARNAUVwithSA[RBPomeUVXAdriDF$proteinName %in% commonProts]


myDataDF <- data.frame(oldFCs=RBPomeUVXOldProtsFCs,newFCs=RBPomeUVXAdriProtsFCs)
corScatterSA <- ggplot(data = myDataDF, mapping = aes(x = oldFCs, y = newFCs)) +
  geom_point(shape = 21, fill = '#0f993d', color = 'white', size = 3) +
  geom_smooth(method = "lm",se = F)+
  annotate("text", x = -1, y = 4, col = "black", size = 3,
           label = paste("Pearson r = ", signif(cor(RBPomeUVXOldProtsFCs,RBPomeUVXAdriProtsFCs,method = "pearson"),3))) +
  annotate("text", x = 2.5, y = -2,  col = "red", size = 3,
           label = paste("y = ", signif(coef(lm(RBPomeUVXAdriProtsFCs ~ RBPomeUVXOldProtsFCs))[1],3), 
                         signif(coef(lm(RBPomeUVXAdriProtsFCs ~ RBPomeUVXOldProtsFCs))[2],2), "x"))
nameOut <- "/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/corrScaterUVXRBPomeOldvsNew.png"
ggsave(filename = nameOut, plot = corScatterSA, width = 30, height = 15, units = 'cm', dpi = 'print')


commonProts


