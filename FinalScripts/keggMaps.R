library(pathview)


ourdata <- wholeDF[,c("proteinName","ORF",infoColumns)] # infoColumns ## FROM 9massiveheatmap

ORAfile <- '/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/FinalData/EnrichmentResults/ORA/SAFAXrbpmFC.Enr.tsv'
ORAres <- read.delim(ORAfile)
pathway <- ORAres$annotation_id[ORAres$annotation == 'KEGG'][1]

mappingFinalFile <- "yeastReference/mappingFile.xlsx"
yeastGenesProtMap <- read.xlsx(mappingFinalFile)
yeastGenesProtMap <- yeastGenesProtMap[!duplicated(yeastGenesProtMap),]
ncbiMapping <- read.delim('/home/eidriangm/Downloads/proteins_15_22535.csv',sep = ',')

ourdata$RBPomeFAX.log2ratio_PolyARNAFAXwithSA

geneInfo <- ourdata[,c('ORF','DEGs.log2FC.SA.YPD','ProteomeFAX.log2ratio_FAXwithSA','RBPomeFAX.log2ratio_PolyARNAFAXwithSA','FAXnetchangesSA')]
rownames(geneInfo) <- geneInfo$ORF; geneInfo$ORF <- NULL 

mapadata <- pathview(gene.data = geneInfo, pathway.id = pathway, species = "sce", out.suffix = "SCE_SYMBOL_MNULL",
                     gene.idtype="ORF",kegg.native = T,map.symbol = T, same.layer=FALSE, map.null=T,
                     min.nnodes=1,multi.state = T)


#geneInfo <- geneInfo[which(!is.na(geneInfo$DEGs.log2FC.SA.YPD) | !is.na(geneInfo$RBPomeFAX.log2ratio_PolyARNAFAXwithSA)),]
#geneInfo <- geneInfo[complete.cases(geneInfo),]
View(geneInfo)
#geneInfoNCBI <- merge(geneInfo,ncbiMapping,by.x='ORF',by.y='Locus.tag')
#IDsMapping <- cbind(ORF=geneInfoNCBI$ORF,ENTREZID=geneInfoNCBI$GeneID)

geneInfo2Map <- geneInfo$DEGs.log2FC.SA.YPD
names(geneInfo2Map) <- geneInfo$ORF

pathview(gene.data = geneInfo2Map, pathway.id = pathway, species = "sce", out.suffix = "SCE_DEFAULT",
         gene.idtype="ORF",kegg.native = T,map.symbol = F, same.layer=T)

pathview(gene.data = geneInfo2Map, pathway.id = pathway, species = "sce", out.suffix = "SCE_SYMBOL",
         gene.idtype="ORF",kegg.native = T,map.symbol = T, same.layer=FALSE)

gse16873.d[, 1:3]

mapadata <- pathview(gene.data = geneInfo, pathway.id = pathway, species = "sce", out.suffix = "SCE_SYMBOL_MNULL",
            gene.idtype="ORF",kegg.native = T,map.symbol = T, same.layer=FALSE, map.null=T,
            min.nnodes=1,multi.state = T)

View(mapadata$plot.data.gene)
mapadata$
#  TOO UGLY and ONLY PDF but WE REMOVE ENZYMES OF OTHER ORGs 
# pathview(gene.data = geneInfo2Map, pathway.id = pathway, species = "sce", out.suffix = "SCE_SYMBOL_MNULL_ALT",
#          gene.idtype="ORF",kegg.native = F,map.symbol = T, same.layer=FALSE, map.null=T)
pathview(gene.data = geneInfo2Map, pathway.id = pathway, species = "sce", out.suffix = "SCE_SYMBOL_ALTEXP",
         gene.idtype="ORF",kegg.native = F,map.symbol = T, same.layer=FALSE, map.null=T,expand.node=T)




geneInfo2Map



id.map = id.map.ensprot


colnames(id.map.ensprot)



## MAPPING
id.map.ensprot <- id2eg(ids = names(gene.ensprot), category = gene.idtype.list[4], org = "Hs")



require(org.Hs.eg.db)
gse16873.t <- apply(gse16873.d, 1, function(x) t.test(x, alternative = "two.sided")$p.value)


sel.genes <- names(gse16873.t)[gse16873.t < 0.1]

sel.cpds <- names(sim.cpd.data)[abs(sim.cpd.data) > 0.5]

pv.out <- pathview(gene.data = sel.genes, cpd.data = sel.cpds, id.map=id.map.ensprot,
                   pathway.id = demo.paths$sel.paths[i], species = "hsa", out.suffix = "sel.genes.sel.cpd", 
                   keys.align = "y", kegg.native = T, key.pos = demo.paths$kpos1[i], 
                   limit = list(gene = 5, cpd = 2), bins = list(gene = 5, cpd = 2), 
                   na.col = "gray", discrete = list(gene = T, cpd = T))

pv.out <- pathview(gene.data = sel.genes, cpd.data = sim.cpd.data, 
                   pathway.id = demo.paths$sel.paths[i], species = "hsa", out.suffix = "sel.genes.cpd",
                   keys.align = "y", kegg.native = T, key.pos = demo.paths$kpos1[i], 
                   limit = list(gene = 5, cpd = 1), bins = list(gene = 5, cpd = 10), 
                   na.col = "gray", discrete = list(gene = T, cpd = F))