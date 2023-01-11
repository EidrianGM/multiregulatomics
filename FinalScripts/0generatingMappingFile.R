library(plyr)
#################################
##### CREATING MAPPING DATA #####
#################################

yeastdbGFFfile <- "yeastReference/S288C_reference_genome_Current_Release/S288C_reference_genome_R64-3-1_20210421/NOFASTAsaccharomyces_cerevisiae_R64-3-1_20210421.gff"
yeastdbGFFdf <- read.table(yeastdbGFFfile,quote = "",sep = "\t",comment.char = "#",header = F,as.is = T,fill = T,encoding = "utf-8")
table(yeastdbGFFdf$V7)
yeastdbGFFdf$V7[yeastdbGFFdf$V7 == 0] <- "*"
yeastdbGFFdf$V7[yeastdbGFFdf$V7 == "."] <- "*"
table(yeastdbGFFdf$V7)

yeastdbGFFmyfile <- 'yeastReference/S288C_reference_genome_Current_Release/S288C_reference_genome_R64-3-1_20210421/myyeast.gff'
write.table(yeastdbGFFdf,yeastdbGFFmyfile,col.names = F,row.names = F, sep = "\t",quote = F)
yeastdbGFFdf2 <- read.delim(yeastdbGFFmyfile,quote = "",header = F,sep = "\t",comment.char = "#",stringsAsFactors = F)
table(yeastdbGFFdf2$V7)

yeastdbGFF <- as.data.frame(rtracklayer::import(yeastdbGFFmyfile))
paste(colnames(yeastdbGFF),collapse = "','")
allcolumns <- c('seqnames','start','end','width','strand','source','type','score','phase','ID','dbxref','Name','Note','display','curie','Parent','Ontology_term','orf_classification','protein_id','Alias','gene','transcript_id','conditions')
table(yeastdbGFF$score)
table(yeastdbGFF$type)

micron2gffile <- "yeastReference/NOFASTAscerevisiae_2-micron.gff"
micron2gffDF <- as.data.frame(rtracklayer::import(micron2gffile))

yeastdbGFF <- rbind.fill(yeastdbGFF,micron2gffDF)

yeastdbGFF_genes <- yeastdbGFF[yeastdbGFF$type != "CDS",]
table(yeastdbGFF_genes$protein_id)
yeastdbGFF_genes$protein_id <- NULL

yeastdbGFF_proteins <- yeastdbGFF[yeastdbGFF$type == "CDS",]
yeastdbGFF_proteins$Name <- gsub("_CDS","",yeastdbGFF_proteins$Name)
yeastdbGFF_proteinsSub <- yeastdbGFF_proteins[c("Name","protein_id")]

geneNprots <- merge(yeastdbGFF_genes,yeastdbGFF_proteinsSub,by.x="ID",by.y="Name",all.x=T)
geneNprots$dbxref <- gsub("SGD:","",geneNprots$dbxref)
geneNprots$protein_id <- gsub("UniProtKB:","",geneNprots$protein_id)

interestingCols <- c('seqnames','start','end','type','ID','dbxref','protein_id','Alias','gene')
selectedMapping <- geneNprots[interestingCols]
selectedMapping

uniprotFile <- "yeastReference/uniprot-compressed_true_download_true_fields_accession_2Cid_2Cprotei-2022.09.07-13.46.19.80.tsv"
uniprotDF <- read.delim(uniprotFile)
fullmerge <- merge(selectedMapping,uniprotDF,by.x="protein_id",by.y="Entry",all.x=T)

paste(colnames(fullmerge),collapse = "','")
fullmerge <- fullmerge[c('seqnames','start','end','type','ID','dbxref','gene','protein_id','Entry.Name','Gene.Names','Protein.names','Gene.Names..ordered.locus.','Gene.Names..ORF.','Gene.Names..primary.','Gene.Names..synonym.','Alias')]
fullmerge$Alias <- unlist(lapply(fullmerge$Alias,function(x) paste0(x,collapse = ' | ')))
colnames(fullmerge) <- c('seqnames','start','end','type','ORF','SGDID','GeneName','UniprotACC','UniprotName','Gene.Names','Protein.names','Gene.Names..ordered.locus.','Gene.Names..ORF.','Gene.Names..primary.','Gene.Names..synonym.','Alias')

# Add manually from UniProt the info for P0CX68
mappingFile <- rbind(fullmerge,c("","","","transposable_element_gene","","L1823","TY1A-LR1","P0CX68","YL11A_YEAST","TY1A-LR1; YLRCTy1-1 GAG","Transposon Ty1-LR1 Gag polyprotein","YLR035C-B","L1823","TY1A-LR1","YLRCTy1-1 GAG","D3DM68; Q12231"))

outdir <- "yeastReference"
saveTablesTsvExc(mappingFile,outdir,completeNdedup=F,excel=T,bycompleteFC=F,rownames=F)

