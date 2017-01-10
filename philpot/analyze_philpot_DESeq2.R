### philpot mice differential expression analysis 
# qsub -V -l mf=200G,h_vmem=220G,h_stack=256M -cwd -b y R CMD BATCH analyze_philpot_DESeq2.R
source('../DESeq2_functions.R') #work-horse of differential expression

##############################################
# load phenotype data and RPKM expression data
load('./rdas/pheno.rda',envir = dat<-new.env())
load('/dcl01/lieber/ajaffe/Brady/philpot/rawCounts_philpot_OCT20_n58.rda')
pd = cbind(dat$pd,pd)
all(seq(nrow(pd))==unlist(sapply(pd$SAMPLE_ID,grep,pd$'BGI file name'))) #samples line up
rownames(pd) = pd$SAMPLE_ID

##############
# split by age
indList = split(seq(nrow(pd)), pd$Age)

##############################
# create and run DESeq objects
geneDds <- lapply(indList,function(i) {
  tmp = pd[i,]
  tmp$Genotype = droplevels(tmp$Genotype)
  DESeq2(countData = geneCounts[,i], colData = tmp, 
         design = ~Genotype,sva = TRUE,parallel=TRUE)})

############################################
# get DE results, and fold-change PTHS v. WT
resGene <- lapply(c('P1','Adult'),function(a){ 
  g = geneDds[[a]]
  lines = unique(colData(g)$Line)
  tmp =lapply(lines,function(l) results(g,contrast = c('Genotype',l,'WT'), alpha=0.05))
  names(tmp) = paste0(a,'.',lines)
  tmp})
resGene = unlist(resGene)

sapply(resGene,function(g) sum(g$padj < 0.05, na.rm=TRUE))

outGeneList <- lapply(resGene,function(g) {
  outGene <- as.data.frame(g)
  outGene = outGene[order(outGene$padj,outGene$pvalue),]
  outGene = cbind(outGene,geneMap[rownames(outGene),])
})

sigGeneList = lapply(outGeneList,function(g){
  sigGene = g[which(g$pvalue<.01),]
})

#######################################
# save all the differential expressions
library(WriteXLS)
WriteXLS(sigGeneList, ExcelFileName = 'tables/philpot_DE_table_ages.xls')
save(pd, outGeneList, file="rdas/philpot_DE_by_age_objects_DESeq2.rda")
save(geneDds,file = '/dcl01/lieber/ajaffe/Brady/philpot/philpot_by_age_DESeq2_svaAdj.rda')

