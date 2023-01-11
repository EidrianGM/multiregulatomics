library(limma)
library(edgeR)

oldRes <- read.delim("/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/data/RNA-Seq data/result/DE genes Adrian.csv",sep = ",")
# transcriptsFile <- "data/RNA-Seq data/result/TPM_genes.csv" ## NOT THE SAME close  with voom but not with voomweighted
# transcriptsFile <- "data/RNA-Seq data/result/counts_genes.csv" ## NOT THE SAME not close nor with voom nor voomweighted
transcriptsFile <- "data/RNA-Seq data/result/counts_trans.csv"
exprMatrix <- read.delim(transcriptsFile, sep = ",")
nrow(exprMatrix)
length(unique(exprMatrix$X))
length(unique(gsub("_.*","",exprMatrix$X)))
exprMatrix$X <- gsub("_.*","",exprMatrix$X)
row.names(exprMatrix) <- exprMatrix$X
exprMatrix$X <- NULL


class <- gsub("\\..*","",colnames(exprMatrix))
pval_adj_method <- "BH"
norm_method <- 'TMM'
design <- condition2design(class)
contrast <- c("H202-YPD","DTT-YPD","SA-YPD")

dge <- DGEList(counts = exprMatrix,
               group  = class)
genes_dge <- suppressWarnings(calcNormFactors(dge,method = norm_method))
# myResultsNOnorm <- limma.pipeline(dge,contrast,0.5,design)
myResults <- limma.pipeline(genes_dge,contrast,0.5,design)

myResults$DE.pval


RNAseqRepaired <- myResults$DE.stat
outdir <- "/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/data/RNA-Seq data/correction"
saveTablesTsvExc(RNAseqRepaired,outdir,completeNdedup=F,excel=T,bycompleteFC=F,rownames=F)

View(myResults)


ibtiResfile <- "/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/data/RNA-Seq data/result/DE genes Adrian.csv"
ibtiResDF <- read.delim(ibtiResfile,sep = ",")
ibtiResDF <- ibtiResDF[which(ibtiResDF$target %in% row.names(myResults$DE.pval)),]

myResultsFCs <- as.data.frame(myResults$DE.lfc)
myResultsPvals <- as.data.frame(myResults$DE.pval)

for (mycontrast in contrast){
  ibSAYPDdf <- ibtiResDF[ibtiResDF$contrast == mycontrast,]
  row.names(ibSAYPDdf) <- ibtiResDF$target[ibtiResDF$contrast == mycontrast]
  
  comparableGens <- intersect(row.names(ibSAYPDdf),row.names(myResultsFCs))
  
  ibSAYPDdf <- ibSAYPDdf[comparableGens,]
  myFCs <- myResultsFCs[comparableGens, mycontrast]

  x <- myFCs
  y <- ibSAYPDdf$log2FC
  toCorPlotDF <- as.data.frame(cbind(x=x,y=y))
  
  corScatterSA <- ggplot(data = toCorPlotDF, mapping = aes(x = x, y = y)) +
    geom_point(shape = 21, fill = '#0f993d', color = 'white', size = 3) +
    geom_smooth(method = "lm",se = F) +
    scale_color_discrete(name = 'Datatype') +
    annotate("text", x = -1, y = 4, col = "black", size = 3,
             label = paste("Pearson r = ", signif(cor(x,y,method = "pearson"),3))) +
    annotate("text", x = 2.5, y = -2,  col = "red", size = 3,
             label = paste("y = ", signif(coef(lm(x ~ y))[1],3), 
                           signif(coef(lm(x ~ y))[2],2), "x"))
  
  nameOut <- paste0("data/RNA-Seq data/correction/",mycontrast,"FCs_corrplot.png")
  ggsave(filename = nameOut, plot = corScatterSA, width = 30, height = 15, units = 'cm', dpi = 'print')
  
    myPvals <- myResultsPvals[comparableGens, mycontrast]
  x <- -log10(myPvals)
  y <- -log10(ibSAYPDdf$adj.pval)
  toCorPlotDF <- as.data.frame(cbind(x=x,y=y))
  
  corScatterSA <- ggplot(data = toCorPlotDF, mapping = aes(x = x, y = y)) +
    geom_point(shape = 21, fill = '#0f993d', color = 'white', size = 3) +
    geom_smooth(method = "lm",se = F) +
    scale_color_discrete(name = 'Datatype') +
    annotate("text", x = -1, y = 4, col = "black", size = 3,
             label = paste("Pearson r = ", signif(cor(x,y,method = "pearson"),3))) +
    annotate("text", x = 2.5, y = -2,  col = "red", size = 3,
             label = paste("y = ", signif(coef(lm(x ~ y))[1],3), 
                           signif(coef(lm(x ~ y))[2],2), "x"))
  
  nameOut <- paste0("data/RNA-Seq data/correctionPlots/",mycontrast,"-log10Pvals_corrplot.png")
  ggsave(filename = nameOut, plot = corScatterSA, width = 30, height = 15, units = 'cm', dpi = 'print')
}

load("/home/eidriangm/Desktop/toDo/surrey/multiregulatomics/data/RNA-Seq data/data/genes_dge.RData")
# # genes_dge <- suppressWarnings(calcNormFactors(dge,method = norm_method))
# class <- gsub("\\..*","",colnames(genes_dge))
# designDGE <- condition2design(class)

# myresultsExact <- limma.pipeline(genes_dge,contrast,0.5,designDGE)
# myresultsWeigh <- limma.pipeline(genes_dge,contrast,0.5,designDGE,voomWeights = T)
# myresultsglmQL <- edgeR.pipeline(dge = genes_dge,design = design,deltaPS = NULL,contrast = contrast,diffAS = F,method = 'glmQL',adjust.method = pval_adj_method)
# myresultsglm <- edgeR.pipeline(dge = genes_dge,design = design,deltaPS = NULL,contrast = contrast,diffAS = F,method = 'glm',adjust.method = pval_adj_method)





cor(mySAYPDdf,ibSAYPDdf$log2FC)

View(mySAYPDdf)

myresultsExact$DE.pval["YGR088W",]

myresultsWeigh$DE.pval["YGR088W",]
myresultsglmQL$DE.pval["YGR088W",]
myresultsglm$DE.pval["YGR088W",]


myresultsExact <- limma.pipeline(genes_dge,contrast,0.5,designDGE,voomWeights = T)

myresultsExact <- limma.pipeline(genes_dge,contrast,0.5,designDGE,voomWeights = T)

#myresults <- limma.pipeline(dge,contrast,0.5,design)


# myresultsW <- limma.pipeline(dge,contrast,0.5,design,voomWeights = T)

myresults$DE.pval["snR31",]
myresults$DE.lfc["snR31",]

myresultsW$DE.pval["snR31",]
myresultsW$DE.lfc["snR31",]

limma.pipeline <- function(dge,
                           contrast,
                           span=0.5,
                           design,
                           deltaPS=NULL,
                           diffAS=F,
                           adjust.method='BH',
                           block=NULL,
                           voomWeights=F,...){
  start.time <- Sys.time()
  results <- list()
  if(is.null(dge$genes)){
    diffAS <- F
  } else {
    if(max(table(dge$genes$GENEID)) < 2)
      diffAS <- F
  }
  
  ##########################################################
  ##--limma voom
  cat('> Limma-voon to estimate mean-vriance trend ...','\n')
  if(voomWeights){
    voom.object<-voomWithQualityWeights(dge,design,plot=F,span = span)
  } else {
    voom.object<-voom(dge,design,plot=F,span = span,...)
  }
  
  results$voom.object <- voom.object
  targets <- rownames(voom.object$E)
  
  ##########################################################
  ##--Fit block
  if(!is.null(block)){
    cat('> Fit block information','\n')
    corfit <- duplicateCorrelation(voom.object,design,block=block)
    correlation <- corfit$consensus
    results$corfit <- corfit
  } else {
    correlation <- NULL
  }
  results$block <- block
  results$correlation <- correlation
  
  ##########################################################
  ##--Fit a basic linear model
  cat('> Fit a basic linear model ...','\n')
  fit.lmFit <- lmFit(voom.object, design,block = block,correlation = correlation)
  results$fit.lmFit <- fit.lmFit
  
  ##########################################################
  ##--Fit a basic linear model
  cat('> Fit the contrast model ...','\n')
  contrast.matrix <- makeContrasts(contrasts = contrast, levels=design)
  print(paste0('Contrast groups: ',paste0(contrast,collapse = '; ')))
  fit.contrast<-contrasts.fit(fit.lmFit, contrast.matrix)
  results$fit.contrast<-fit.contrast
  
  ##########################################################
  ##--Fit a eBayes model
  cat('> Fit a eBayes model ...','\n')
  fit.eBayes<-eBayes(fit.contrast)
  results$fit.eBayes<-fit.eBayes
  
  ##########################################################
  ##--Testing statistics for each contrast group
  cat('> Testing for each contrast group ...','\n')
  DE.pval.list <- lapply(contrast,function(i){
    x <- topTable(fit.eBayes,
                  coef=i,
                  adjust.method =adjust.method,
                  number = Inf)
    x <- x[targets,]
    x
  })
  names(DE.pval.list) <- contrast
  results$DE.pval.list<-DE.pval.list
  
  ###---DE pval and lfc
  DE.pval <- do.call(cbind,lapply(DE.pval.list,FUN = function(x) x$adj.P.Val))
  DE.lfc <- do.call(cbind,lapply(DE.pval.list,FUN = function(x) x$logFC))
  rownames(DE.pval) <- rownames(DE.lfc) <- targets
  colnames(DE.pval) <- colnames(DE.lfc) <- contrast
  
  DE.stat <- summaryStat(x = DE.pval,y = DE.lfc,
                         target = rownames(DE.pval),
                         contrast = contrast,
                         stat.type = c('adj.pval','log2FC'))
  
  # DE.stat <- cbind(DE.pval,DE.lfc)
  # colnames(DE.stat) <- c(paste0('pval:',contrast),paste0('lfc:',contrast))
  results$DE.pval<-DE.pval
  results$DE.lfc<-DE.lfc
  results$DE.stat<-DE.stat
  
  ##########################################################
  ##--Testing statistics for across all contrast groups
  cat('> Testing across all contrast groups ...','\n')
  DE.stat.overalltest <- topTable(fit.eBayes,number = Inf,
                                coef = contrast,adjust.method =adjust.method )
  # DE.stat.overalltest <- DE.stat.overalltest[targets,]
  col.idx <- gsub('-','.',contrast)
  col.idx <- grep(paste0(col.idx,collapse = '|'),colnames(DE.stat.overalltest))
  DE.stat.overalltest <- data.frame(target=rownames(DE.stat.overalltest),
                                    contrast='overall',DE.stat.overalltest,row.names = NULL)
  colnames(DE.stat.overalltest)[col.idx] <- gsub('[.]','-',colnames(DE.stat.overalltest)[col.idx])
  results$DE.stat.overalltest<-DE.stat.overalltest
  
  if(diffAS){
    cat('> Fit a splicing model ...','\n')
    if(is.null(deltaPS))
      stop('Please provide deltaPS for DAS analysis...')
    fit.splice<-diffSplice(fit.contrast, geneid = 'GENEID')
    results$fit.splice<-fit.splice
    
    # ##########################################################
    ##---DTU transcripts
    ##DTU transcript pval list
    genes.idx <- unique(fit.splice$genes$GENEID)
    trans.idx <- unique(fit.splice$genes$TXNAME)
    
    DTU.pval.list<-lapply(contrast,function(i){
      y<-topSplice(fit.splice, coef=i, test="t", number=Inf, FDR=10000)
      rownames(y)<-y$TXNAME
      z <- y[trans.idx,]
      z
    })
    names(DTU.pval.list) <- contrast
    
    ##DTU transcript pvals
    DTU.pval <- lapply(contrast,function(i){
      x <- DTU.pval.list[[i]]
      y <- data.frame(pval=x$FDR)
      rownames(y) <- rownames(x)
      y
    })
    DTU.pval <- do.call(cbind,DTU.pval)
    colnames(DTU.pval) <- contrast
    
    ##DTU transcript deltaPS
    DTU.deltaPS <- deltaPS[rownames(DTU.pval),,drop=F]
    
    DTU.stat <- summaryStat (x = DTU.pval,y = DTU.deltaPS,
                             target = rownames(DTU.pval),
                             contrast = contrast,
                             stat.type = c('adj.pval','deltaPS'))
    
    # DTU.stat <- cbind(DTU.pval,DTU.deltaPS)
    # colnames(DTU.stat) <- c(paste0('pval:',contrast),paste0('deltaPS:',contrast))
    
    results$DTU.pval.list<-DTU.pval.list
    results$DTU.pval<-DTU.pval
    results$DTU.deltaPS<-DTU.deltaPS
    results$DTU.stat<-DTU.stat
    
    # ##########################################################
    ##---DAS genes
    ##---max deltaPS
    maxdeltaPS <- by(deltaPS[fit.splice$genes$TXNAME,,drop=F],
                     INDICES = fit.splice$genes$GENEID,
                     function(x){
                       apply(x,2,function(i) i[abs(i)==max(abs(i))][1])
                     },simplify = F)
    maxdeltaPS <- do.call(rbind,maxdeltaPS)
    maxdeltaPS <- maxdeltaPS[genes.idx,,drop=F]
    results$maxdeltaPS<-maxdeltaPS
    
    ##---F test
    DAS.pval.F.list<-lapply(contrast,function(i){
      y<-topSplice(fit.splice, coef=i, test="F", number=Inf, FDR=10000)
      rownames(y)<-y$GENEID
      y[genes.idx,]
    })
    names(DAS.pval.F.list) <- contrast
    
    DAS.pval.F <- lapply(contrast,function(i){
      x <- DAS.pval.F.list[[i]]
      y <- data.frame(pval=x$FDR)
      rownames(y) <- rownames(x)
      y
    })
    
    DAS.pval.F <- do.call(cbind,DAS.pval.F)
    DAS.pval.F <- DAS.pval.F[genes.idx,,drop=F]
    colnames(DAS.pval.F) <- contrast
    
    DAS.F.stat <- summaryStat (x = DAS.pval.F,y = maxdeltaPS,
                               target = rownames(DAS.pval.F),
                               contrast = contrast,
                               stat.type = c('adj.pval','maxdeltaPS'))
    
    # DAS.F.stat <- cbind(DAS.pval.F,maxdeltaPS)
    # colnames(DAS.F.stat) <- c(paste0('pval:',contrast),paste0('MaxdeltaPS:',contrast))
    #
    results$DAS.pval.F.list<-DAS.pval.F.list
    results$DAS.pval.F<-DAS.pval.F
    results$DAS.F.stat<-DAS.F.stat
    
    ##---simes test
    DAS.pval.simes.list<-lapply(contrast,function(i){
      y<-topSplice(fit.splice, coef=i, test="simes", number=Inf, FDR=10000)
      rownames(y)<-y$GENEID
      y[genes.idx,]
    })
    names(DAS.pval.simes.list) <- contrast
    
    DAS.pval.simes <- lapply(contrast,function(i){
      x <- DAS.pval.simes.list[[i]]
      y <- data.frame(pval=x$FDR)
      rownames(y) <- rownames(x)
      y
    })
    
    DAS.pval.simes <- do.call(cbind,DAS.pval.simes)
    DAS.pval.simes <- DAS.pval.simes[genes.idx,,drop=F]
    colnames(DAS.pval.simes) <- contrast
    
    DAS.simes.stat <- summaryStat (x = DAS.pval.simes,y = maxdeltaPS,
                                   target = rownames(DAS.pval.simes),
                                   contrast = contrast,
                                   stat.type = c('adj.pval','maxdeltaPS'))
    #
    # DAS.simes.stat <- cbind(DAS.pval.simes,maxdeltaPS)
    # colnames(DAS.simes.stat) <- c(paste0('pval:',contrast),paste0('MaxdeltaPS:',contrast))
    
    results$DAS.pval.simes.list<-DAS.pval.simes.list
    results$DAS.pval.simes<-DAS.pval.simes
    results$DAS.simes.stat<-DAS.simes.stat
  }
  end.time <- Sys.time()
  time.taken <- end.time - start.time
  cat(paste0('Time for analysis: ',round(time.taken,3)))
  cat('\n','> Done!!! ','\n')
  return(results)
}

condition2design <- function(condition,batch.effect=NULL){
  design.data <- data.frame(condition=condition)
  if(!is.null(batch.effect)){
    colnames(batch.effect) <- paste0('batch',1:ncol(batch.effect))
    design.data <- data.frame(design.data,batch.effect)
  }
  design.fomula <- as.formula(paste0('~0+',paste0(colnames(design.data),collapse = '+')))
  design <- model.matrix(design.fomula, data=design.data)
  idx <- grep('condition',colnames(design))
  design.idx <- colnames(design)
  design.idx <- substr(design.idx,nchar('condition')+1,nchar(design.idx))
  colnames(design)[idx] <- design.idx[idx]
  design
}

summaryStat  <- function(x,y,target,
                         contrast = NULL,
                         stat.type = c('adj.pval','lfc'),
                         srot.by = stat.type[1]){
  x <- x[target,,drop=F]
  y <- y[target,,drop=F]
  if(is.null(contrast))
    contrast <- colnames(x)
  x <- x[,contrast,drop=F]
  y <- y[,contrast,drop=F]
  
  stat <- lapply(contrast,function(i){
    z <- data.frame(targets =target,contrast=i, x[,i],y[,i],row.names = NULL)
    colnames(z) <- c('target','contrast',stat.type)
    z <- z[order(z[,srot.by]),]
  })
  names(stat) <- contrast
  stat <- do.call(rbind,stat)
  rownames(stat) <- NULL
  stat
}

# span <- 0.5
# adjust.method <- 'BH'
# voom.object <- voomWithQualityWeights(dge,design,plot=F,span = span)
# #voom.object <- voom(dge,design,plot=F,span = span)
# targets <- rownames(voom.object$E)
# fit.lmFit <- lmFit(voom.object, design)
# ##########################################################
# ##--Fit a basic linear model
# cat('> Fit the contrast model ...','\n')
# contrast.matrix <- makeContrasts(contrasts = contrast, levels=design)
# print(paste0('Contrast groups: ',paste0(contrast,collapse = '; ')))
# fit.contrast <- contrasts.fit(fit.lmFit, contrast.matrix)
# 
# ##########################################################
# ##--Fit a eBayes model
# cat('> Fit a eBayes model ...','\n')
# fit.eBayes <- eBayes(fit.contrast)
# 
# ##########################################################
# ##--Testing statistics for each contrast group
# cat('> Testing for each contrast group ...','\n')
# DE.pval.list <- lapply(contrast,function(i){
#   x <- topTable(fit.eBayes,
#                 coef=i,
#                 adjust.method = adjust.method,
#                 number = Inf)
#   x <- x[targets,]
#   x
# })
# names(DE.pval.list) <- contrast
# 
# ###---DE pval and lfc
# DE.pval <- do.call(cbind,lapply(DE.pval.list,FUN = function(x) x$adj.P.Val))
# DE.lfc <- do.call(cbind,lapply(DE.pval.list,FUN = function(x) x$logFC))
# rownames(DE.pval) <- rownames(DE.lfc) <- targets
# colnames(DE.pval) <- colnames(DE.lfc) <- contrast
# 
# DE.pval.list$`SA-YPD`["snR31",]


