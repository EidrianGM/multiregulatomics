library(dplyr)
library(forcats)
library(ggplot2)
library(gridExtra)
library(grid)
library(ggplot2)
library(lattice)

# Define function
createPlotFunction <- function(CMDfile, CMDHfile, CMDLfile, nameOut){
  CMD <- CMD_results
  CMDH <- CMDH_results
  CMDL <- CMDL_results
  
  ## Add comparison column
  CMD$Comparison <- rep('CvsMD', nrow(CMD))
  CMDH$Comparison <- rep('CvsMDH', nrow(CMDH))
  CMDL$Comparison <- rep('CvsMDL', nrow(CMDL))
  
  # all comparisons
  all_results <- rbind(CMD, CMDH, CMDL)
  
  # Order factor Comparision
  all_results$Comparison = factor(all_results$Comparison, levels=c('CvsMD','CvsMDH','CvsMDL'))
  # Order by pvalue
  all_results <- all_results[order(all_results$p.adjust_over),]
  all_results <- mutate(all_results, term = fct_reorder(term, desc(p.adjust_over)))
  
  # Set names
  comparison_label <- c(CvsMD = 'C vs MD', CvsMDH = 'C vs MDH', CvsMDL = 'C vs MDL')
  
  # Create plot
  plot_results <- ggplot(all_results, aes(x = count, y = term, size = count, color = p.adjust_over)) +
    geom_point(alpha = 1) +
    scale_size(range = c(3, 8), name = "Count") +
    scale_color_gradientn(colours = c('maroon3', 'gold', 'turquoise3'), name = 'p.adjust') +
    ylab('') +
    xlab('Gene ratio') +
    facet_grid(. ~ Comparison, labeller = labeller(Comparision = comparison_label)) +
    theme(axis.text.y = element_text(size = 10, face = 'bold'),
          strip.text = element_text(hjust = 0.5, size = 12, face = 'bold'))
  plot_results
  
  ggsave(filename = paste0(nameOut, '.png'), plot = plot_results, width = 30, height = 15, units = 'cm', dpi = 'print')
}

enrRes <- "/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/data/enrichmentRes/GSA/DEGsSA.Enr.tsv"
GC4res <- read.delim(enrRes)
GC4res <- GC4res[order(GC4res$pval_adj,decreasing = F),]
GC4res$pval_adj <- -log10(GC4res$pval_adj)
plotList <- list()
for (annot in unique(GC4res$annotation)){
  subGC4res <- GC4res[GC4res$annotation == annot,][1:10,]
  plot_results <- ggplot(subGC4res, aes(x = relative_enrichment, y = description, size = pval_adj, color = genes_found)) +
                  geom_point(alpha = 1) + ggtitle(annot) + 
                  scale_size(range = c(3, 8), name = "-Log10(FDR)") +
                  scale_color_gradientn(colours = c('maroon3', 'gold', 'turquoise3'), name = 'Significant Genes') +
                  ylab('') + xlab('Relative Enrichment') + 
                  theme(axis.text.y = element_text(size = 6, face = 'bold'),
                        strip.text = element_text(hjust = 0.5, size = 12, face = 'bold'))
  plotList[[annot]] <- plot_results
}


ggsave(plotEnrFile,grid.arrange(grobs=plotList, nrow = 2),width = 60,height = 20,units = "cm",dpi = "print",limitsize = F)


enrResFile <- "/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/data/enrichmentRes/GSA/DEGsSA.Enr.tsv"





GC4res$annotation <- as.factor(GC4res$annotation)
GC4res <- GC4res[order(GC4res$pval_adj),]
GC4all_results <- mutate(GC4res, description = fct_reorder(description, desc(pval_adj)))
## Add comparison column
# CMD$Comparison <- rep('CvsMD', nrow(CMD))
# CMDH$Comparison <- rep('CvsMDH', nrow(CMDH))
# CMDL$Comparison <- rep('CvsMDL', nrow(CMDL))

# all comparisons
#all_results <- rbind(CMD, CMDH, CMDL)

# Order factor Comparision
#all_results$Comparison = factor(all_results$Comparison, levels=c('CvsMD','CvsMDH','CvsMDL'))
# Order by pvalue
#all_results <- mutate(CMD_results, term = fct_reorder(term, desc(p.adjust_over)))

# Set names
#comparison_label <- c(CvsMD = 'C vs MD', CvsMDH = 'C vs MDH', CvsMDL = 'C vs MDL')

# Create plot
GC4all_results$relative_enrichment <- as.numeric(GC4all_results$relative_enrichment)
GC4all_results$pval_adj <- as.numeric(GC4all_results$pval_adj)
plot_results <- ggplot(GC4all_results, aes(x = relative_enrichment, y = description, size = relative_enrichment, color = pval_adj)) +
  geom_point(alpha = 1) +
  scale_size(range = c(3, 8), name = "Rel. Enr.") +
  scale_color_gradientn(colours = c('maroon3', 'gold', 'turquoise3'), name = 'p.adjust') +
  ylab('') +
  xlab('Gene ratio') +
  facet_grid(. ~ annotation) +
  theme(axis.text.y = element_text(size = 10, face = 'bold'),
        strip.text = element_text(hjust = 0.5, size = 12, face = 'bold'))
plot_results

ggsave(filename = paste0(nameOut, '.png'), plot = plot_results, width = 30, height = 15, units = 'cm', dpi = 'print')



savePath <- '/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/data'
setwd(savePath)

# Results
dirPath <- '/home/alba/metiloma/results_allCpGs/'
## GO
CMD_results <- read.delim("CvMD_goseq_GO_formula_sign_over_genes.txt", header = T, sep = '\t', stringsAsFactors = F)
CMDH_results <- read.delim('CvMDH_goseq_GO_formula_sign_over_genes.txt', header = T, sep = '\t', stringsAsFactors = F)
CMDL_results <- read.delim('CvMDL_goseq_GO_formula_sign_over_genes.txt', header = T, sep = '\t', stringsAsFactors = F)
createPlotFunction(CMDfile = CMD_results, CMDHfile = CMDH_results, CMDLfile = CMDL_results, nameOut = '/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/data/results_goseq_GO_formula_plot')


## GO - hypo
CMD_results <- read.delim(file = file.path(dirPath, 'CMD/GO/CMD_goseq_formula_hypo/CvMD_hypo_goseq_GO_formula_sign_over_genes.txt'), header = T, sep = '\t', stringsAsFactors = F)
CMDH_results <- read.delim(file = file.path(dirPath, 'CMDH/GO/CMDH_goseq_formula_hypo/CvMDH_hypo_goseq_GO_formula_sign_over_genes.txt'), header = T, sep = '\t', stringsAsFactors = F)
CMDL_results <- read.delim(file = file.path(dirPath, 'CMDL/GO/CMDL_goseq_formula_hypo/CvMDL_hypo_goseq_GO_formula_sign_over_genes.txt'), header = T, sep = '\t', stringsAsFactors = F)
CMD_results <- CMD_results[1:10,]
CMDH_results <- CMDH_results[1:10,]
createPlotFunction(CMDfile = CMD_results, CMDHfile = CMDH_results, CMDLfile = CMDL_results, nameOut = 'results_goseq_GO_hypo_formula_plot_10')


## GO - hypo - CMD: p-adjust = 0.05-0.06
CMD_results_005006 <- read.delim(file = file.path(dirPath, 'CMD/GO/CMD_goseq_formula_hypo/CvMD_hypo_goseq_GO_formula_sign_005_006_over_genes.txt'), header = T, sep = '\t', stringsAsFactors = F)
CMDH_results <- read.delim(file = file.path(dirPath, 'CMDH/GO/CMDH_goseq_formula_hypo/CvMDH_hypo_goseq_GO_formula_sign_over_genes.txt'), header = T, sep = '\t', stringsAsFactors = F)
CMDL_results <- read.delim(file = file.path(dirPath, 'CMDL/GO/CMDL_goseq_formula_hypo/CvMDL_hypo_goseq_GO_formula_sign_over_genes.txt'), header = T, sep = '\t', stringsAsFactors = F)
CMDH_results_filt <- CMDH_results[CMDH_results$category %in% CMD_results_005006$category,]
CMDL_results_filt <- CMDL_results[CMDL_results$category %in% CMD_results_005006$category,]
createPlotFunction(CMDfile = CMD_results_005006, CMDHfile = CMDH_results_filt, CMDLfile = CMDL_results_filt, nameOut = 'results_goseq_GO_hypo_formula_plot_005006')


## KEGG
CMD_results <- read.delim(file = file.path(dirPath, 'CMD/KEGG/CMD_goseq_formula_KEGG/CvMD_goseq_KEGG_formula_sign_over_genes.txt'), header = T, sep = '\t', stringsAsFactors = F)
CMDH_results <- read.delim(file = file.path(dirPath, 'CMDH/KEGG/CMDH_goseq_formula_KEGG/CvMDH_goseq_KEGG_formula_sign_over_genes.txt'), header = T, sep = '\t', stringsAsFactors = F)
CMDL_results <- read.delim(file = file.path(dirPath, 'CMDL/KEGG/CMDL_goseq_formula_KEGG/CvMDL_goseq_KEGG_formula_sign_over_genes.txt'), header = T, sep = '\t', stringsAsFactors = F)
createPlotFunction(CMDfile = CMD_results, CMDHfile = CMDH_results, CMDLfile = CMDL_results, nameOut = 'results_goseq_KEGG_formula_plot')


## KEGG - hypo
CMD_results <- read.delim(file = file.path(dirPath, 'CMD/KEGG/CMD_goseq_formula_KEGG_hypo/CvMD_hypo_goseq_KEGG_formula_sign_over_genes.txt'), header = T, sep = '\t', stringsAsFactors = F)
CMDH_results <- read.delim(file = file.path(dirPath, 'CMDH/KEGG/CMDH_goseq_formula_KEGG_hypo/CvMDH_hypo_goseq_KEGG_formula_sign_over_genes.txt'), header = T, sep = '\t', stringsAsFactors = F)
CMDL_results <- read.delim(file = file.path(dirPath, 'CMDL/KEGG/CMDL_goseq_formula_KEGG_hypo/CvMDL_hypo_goseq_KEGG_formula_sign_over_genes.txt'), header = T, sep = '\t', stringsAsFactors = F)
createPlotFunction(CMDfile = CMD_results, CMDHfile = CMDH_results, CMDLfile = CMDL_results, nameOut = 'results_goseq_KEGG_hypo_formula_plot')




