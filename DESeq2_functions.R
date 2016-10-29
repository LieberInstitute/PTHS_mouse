# Find Differential Expression using DESeq2 controlling for batch effects with SVA
DESeq2 = function(countData,colData,design, sva = TRUE,parallel=TRUE){
  suppressPackageStartupMessages(require(DESeq2))
  
  #############################
  # initiate DESeqDataSet object
  cat('Creating DESeq object.\n')
  geneDds <- DESeqDataSetFromMatrix(countData = countData, colData = colData, design = design)
  
  ############################################################
  # estimate sample size factors, catch data with sparse counts
  cat('Estimating sample size factors.\n')
  geneDds = tryCatch(estimateSizeFactors(geneDds),
           error = function(err){
           print(paste("MY_ERROR:  ",err))
           countData = countData+1 #add pseudocount
           geneDds <- DESeqDataSetFromMatrix(countData = countData, colData = colData, design = design)
           geneDds <- estimateSizeFactors(geneDds)
           return(geneDds)
           })
  
  ##############################
  # filter lowly expressed genes
  cat('Filtering genes.\n')
  geneDds = geneDds[rowMeans(counts(geneDds))>1,]
  
  ###############################################
  # if need to adjust for unmodeled batch effects
  if(sva){
    suppressPackageStartupMessages(require(sva))
    ##################################
    # statistical model and null model
    mod <- model.matrix(design, colData)
    mod0 <- model.matrix(~ 1, colData)
    
    ###########################
    # extract normalized counts
    geneDat <- counts(geneDds, normalized=TRUE)
    geneDat <- geneDat[rowMeans(geneDat) > 1,]
    
    ########################
    # find surrogate variables
    cat('Estimating surrogate variables.\n')
    svaGene <- svaseq(geneDat, mod, mod0)
    cat('\n')
    
    #################################
    # add surrogate variables to model
    cat('Adding surrogate variables.\n')
    svaGene$sv = data.frame(svaGene$sv)
    colnames(svaGene$sv) = c(paste0('SV',seq(svaGene$n.sv)))
    for (x in colnames(svaGene$sv)) colData(geneDds)[[x]] <- svaGene$sv[,x]
    design(geneDds) = as.formula(paste('~',paste(c(as.character(design)[2],colnames(svaGene$sv)), collapse= "+")))
  }
  
  ###########################
  # find DEGs and fold change
  if(parallel){
    suppressPackageStartupMessages(require(parallel))
    suppressPackageStartupMessages(require(BiocParallel))
    register(MulticoreParam(detectCores()))
  }
  cat('Finding differential expression.\n')
  geneDds <- DESeq(geneDds,parallel=parallel)
  cat('Finished.\n')
  return(geneDds)
}