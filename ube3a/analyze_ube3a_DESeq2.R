### ube3a mice differential expression analysis 
# qsub -V -l mf=200G,h_vmem=220G,h_stack=256M -cwd -b y R CMD BATCH analyze_ube3a_DESeq2.R
source('../DESeq2_functions.R') #work-horse of differential expression

##############################################
# load phenotype data and RPKM expression data
load('./rdas/pheno.rda',envir = dat<-new.env())
load('/dcl01/lieber/ajaffe/Brady/ube3a/rawCounts_ube3a_nov17_n8.rda')
pd = cbind(dat$pd,pd)
all(seq(nrow(pd))==unlist(sapply(pd$SAMPLE_ID,grep,pd$'BGI file name'))) #samples line up
rownames(pd) = pd$SAMPLE_ID

#####################
# split by age groups
indList = split(seq(nrow(pd)), pd$Age)
names(indList)

##############################
# create and run DESeq objects
geneDds <- lapply(indList,function(i) {
  DESeq2(countData = geneCounts[,i], colData = pd[i,], 
         design = ~Case,sva = FALSE,parallel=TRUE)})

############################################################
# get DE results, and fold-change of each het mouse genotype
resGene <-lapply(geneDds,function(g) results(g,contrast = c('Case','Mutant','Control'), alpha=0.05))

sapply(resGene,function(g) sum(g$padj < 0.05, na.rm=TRUE))
# P1 Adult 
# 18    30

outGeneList <- lapply(resGene,function(g) {
  outGene <- as.data.frame(g)
  outGene = outGene[order(outGene$padj,outGene$pvalue),]
  outGene = cbind(outGene,geneMap[rownames(outGene),])
})

sigGeneList = lapply(outGeneList,function(g){
  sigGene = g[which(g$pvalue<.01),]
})

pdf('plots/DESeq2_MA_plots_ube3a_ages.pdf')
plotMA(resGene[['P1']], main="Ube3a P14 MA plot", ylim=c(-1,1))
plotMA(resGene[['Adult']], main="Ube3a Adult MA plot", ylim=c(-1,1))
dev.off()

#######################################
# save all the differential expressions
library(WriteXLS)
WriteXLS(sigGeneList,ExcelFileName = 'tables/ube3a_DE_ages_table.xls')
save(sigGeneList,outGeneList, file="rdas/ube3a_DE_ages_objects_DESeq2.rda")
save(geneDds,resGene, file = '/dcl01/lieber/ajaffe/Brady/ube3a/ube3a_ages_DESeq2_svaAdj.rda')

