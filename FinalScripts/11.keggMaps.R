library(tidyr)

getKEGGresults <- function(ORAfiles){
  allORAresults <- c()
  for (ORAfile in ORAfiles){
    ORAdf <- read.delim(ORAfile)
    ORAdf <- ORAdf[ORAdf$annotation == 'KEGG',]
    if (grepl("DEGs",ORAfile)){
      dataType <- "Transc."
    }else if (grepl("rbpm",ORAfile)){
      dataType <- "mRBP"
    }else if (grepl("prot",ORAfile)){
      dataType <- "Prot."
    }else if (grepl("netc",ORAfile)){
      dataType <- "Net C."
    }
    else{
      dataType <- "Unknown"
    }
    if (grepl("up",basename(ORAfile),ignore.case = T)){
      # Manually determine the order of the data types
      datatypeOrder <- c("Transc.\nUp","Transc.\nDown","Prot.\nUp","Prot.\nDown","mRBP\nUp","mRBP\nDown","Net C.\nUp","Net C.\nDown") 
      regul <- "\nUp"
    }else if(grepl("dw",basename(ORAfile),ignore.case = T)){
      datatypeOrder <- c("Transc.\nUp","Transc.\nDown","Prot.\nUp","Prot.\nDown","mRBP\nUp","mRBP\nDown","Net C.\nUp","Net C.\nDown") 
      regul <- "\nDown"
    }else{
      datatypeOrder <- c("Transc.","Prot.","mRBP","Net C.") 
      regul <- ""
    }
    ORAdf$datatype <- paste0(dataType,regul)
    allORAresults <- rbind(allORAresults,ORAdf)
  }
  return(allORAresults)
}


### GET FOLD CHANGES
wholeDFfile <- "/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/FinalData/allDataDF.tsv"
wholeDF <- read.delim(wholeDFfile,quote = "")
FAXacceptedProts <- read.delim("FinalData/BackgroundRemoval/FAXAccBackGrAcceptedFC3.tsv",header = F)[,1]
UVXacceptedProts <- read.delim("FinalData/BackgroundRemoval/UVXAccBackGrAcceptedFC3.tsv",header = F)[,1]
wholeDF <- wholeDF[wholeDF$proteinName %in% c(FAXacceptedProts,UVXacceptedProts),]
row.names(wholeDF) <- wholeDF$ORF 
outdir <- "/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/FinalData/HeatMaps"
FAXinfoColumns <- c("DEGs.log2FC.SA.YPD", "ProteomeFAX.log2ratio_FAXwithSA",
                    "RBPomeFAX.log2ratio_PolyARNAFAXwithSA",   "FAXnetchangesSA")

UVXinfoColumns <- c("DEGs.log2FC.SA.YPD", "ProteomeUVX.log2ratio_Condition4",
                  "RBPomeUVX.log2ratio_PolyARNAUVwithSA","UVXnetchangesSA")

keggdir <- "/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/FinalData/EnrichmentResults/Visualizations/KEGGmaps/db"

treatment <- "SA"; crosslink <- "UVX"
ORAfilesDir <- "/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/FinalData/EnrichmentResults/ORA"
ORAfiles <- list.files(ORAfilesDir,pattern = paste0(treatment,".*(DEGs|",crosslink,").*\\.Enr.tsv"),full.names = T)
KEGGres <- getKEGGresults(ORAfiles)
pathwaysEnriched <- unique(KEGGres$annotation_id)

geneInfo <- wholeDF[,UVXinfoColumns]

# geneInfo['ORF'] <- row.names(geneInfo)
# geneInfo <- as.data.frame(separate_rows(geneInfo,'ORF',sep = ";"))
# row.names(geneInfo) <- geneInfo$ORF 

outsuffix <- paste0(crosslink,treatment)
outdir <- paste0("/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/FinalData/EnrichmentResults/Visualizations/KEGGmaps/",crosslink)
setwd(outdir)
for (pathway in pathwaysEnriched){
  pathview(gene.data = geneInfo, pathway.id = pathway, species = "sce", out.suffix = outsuffix,
           kegg.dir = keggdir, gene.idtype="ORF", kegg.native = T, map.symbol = T, 
           same.layer=FALSE, map.null=T, min.nnodes=1, multi.state = T)
}

treatment <- "SA"; crosslink <- "FAX"
ORAfilesDir <- "/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/FinalData/EnrichmentResults/ORA"
ORAfiles <- list.files(ORAfilesDir,pattern = paste0(treatment,".*(DEGs|",crosslink,").*\\.Enr.tsv"),full.names = T)
KEGGres <- getKEGGresults(ORAfiles)
pathwaysEnriched <- unique(KEGGres$annotation_id)

geneInfo <- wholeDF[,FAXinfoColumns]

outsuffix <- paste0(crosslink,treatment)
outdir <- paste0("/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/FinalData/EnrichmentResults/Visualizations/KEGGmaps/",crosslink)
setwd(outdir)
for (pathway in pathwaysEnriched){
  pathview(gene.data = geneInfo, pathway.id = pathway, species = "sce", out.suffix = outsuffix,
           kegg.dir = keggdir, gene.idtype="ORF", kegg.native = T, map.symbol = T, 
           same.layer=FALSE, map.null=T, min.nnodes=1, multi.state = T)
}






