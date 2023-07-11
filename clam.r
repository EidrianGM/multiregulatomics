### CLAM
# http://sgd-archive.yeastgenome.org/curation/literature/

yeastInteractionFile <- '/home/eidrian/Desktop/realServices/CLAM/DB/interaction_data.tab'
yeastInteractionDF <- read.delim(yeastInteractionFile,sep = '\t',header = F)

# 1) Feature Name (Bait) (Required)       	- The feature name of the gene used as the bait
# 2) Standard Gene Name (Bait) (Optional) 	- The standard gene name of the gene used as the bait
# 3) Feature Name (Hit) (Required)        	- The feature name of the gene that interacts with the bait
# 4) Standard Gene Name (Hit) (Optional)  	- The standard gene name of the gene that interacts with the bait
# 5) Experiment Type (Required)   		- A description of the experimental used to identify the interaction
# 6) Genetic or Physical Interaction (Required)   - Indicates whether the experimental method is a genetic or physical interaction
# 7) Source (Required)    			- Lists the database source for the interaction
# 8) Manually curated or High-throughput (Required)	- Lists whether the interaction was manually curated from a publication or added as part of a high-throughput dataset
# 9) Notes (Optional)     			- Free text field that contains additional information about the interaction
# 10) Phenotype (Optional)        		- Contains the phenotype of the interaction
# 11) Reference (Required)        		- Lists the identifiers for the reference as an SGDID (SGD_REF:) or a PubMed ID (PMID:)
# 12) Citation (Required) 			- Lists the citation for the reference
yeastInteractionDF$V9

library(tidyverse)
library(rvest)

responseGET <- read_html(urlrequest)
tables <- html_nodes(responseGET, "table")
head(tables)
mytable <- tables[[5]]
myDF <- html_table(mytable, fill = F)
View(myDF)



responseGET$doc

View(responseGET)

?html_table

tables <- html_table(responseGET, fill = TRUE)

View(tables[[5]])


mytable2 <- mytable[3:nrow(mytable),]

datainvector <- as.character(mytable)
head(datainvector)

myDF <- as.data.frame(matrix(as.character(mytable),ncol = 8))



View(as.data.frame(mytable))




urlrequest <- 'http://yeastract-plus.org/yeastract/scerevisiae/view.php?existing=regulation&proteinname=YJL206C'

library(httr)
responseGET <- GET(urlrequest)
write(content(responseGET, "text"),'response.txt')
responseGET$content

myDF <- yeastInteractionDF[c('V2','V4')]
myDF <- myDF[!duplicated(myDF),]
colnames(myDF) <- c('Protein','Protein')
myDF$Input <- '1'
outfolder <- '/home/eidrian/Desktop/realServices/CLAM/DB/processed'
outFile <- file.path(outfolder,'ppi.tsv')
write.table(myDF,file = outFile,row.names = F,quote = F, sep = '\t')


yeastTFsInteractionFile <- '/home/eidrian/Desktop/realServices/CLAM/DB/yeastract_tfs-targets.csv'
yeastTFsInteractionDF <- read.delim(yeastTFsInteractionFile,sep = ';',header = F)
myDF <- yeastInteractionDF[c('V2','V3')]
myDF <- myDF[!duplicated(myDF),]
colnames(myDF) <- c('Protein','Protein')
myDF$Input <- '1'
outfolder <- '/home/eidrian/Desktop/realServices/CLAM/DB/processed'
outFile <- file.path(outfolder,'ppi.tsv')
write.table(myDF,file = outFile,row.names = F,quote = F, sep = '\t')




table(yeastInteractionDF$V10)


View(yeastInteractionDF)
