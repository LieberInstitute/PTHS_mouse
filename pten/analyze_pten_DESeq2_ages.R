### PTEN mouse differential expression analysis using DESeq2
# qsub -V -l mf=100G,h_vmem=150G,h_stack=256M -cwd -b y R CMD BATCH analyze_pten_DESeq2.R
source('../DESeq2_functions.R') #work-horse of differential expression

##############################################
# load phenotype data and RPKM expression data
load('./rdas/pheno.rda',envir = dat<-new.env())
load('/dcl01/lieber/ajaffe/Brady/mouseRNAseq/TilotPTEN/rawCounts_pten_OCT20_n18.rda')
pd = cbind(dat$pd,pd)
all.equal(pd$FileID,pd$SAMPLE_ID) #samples line up
rownames(pd) = pd$SAMPLE_ID
pd$Age = factor(pd$Age,levels = levels(pd$Age),labels = c('P14','Adult'))

#####################
# split by age groups
indList = split(seq(nrow(pd)), pd$Age)
names(indList)

##############################
# create and run DESeq objects
geneDds <- lapply(indList,function(i) {
  DESeq2(countData = geneCounts[,i], colData = pd[i,], 
         design = ~Genotype,sva = TRUE,parallel=TRUE)})

############################################
# get DE results, and fold-change PTHS v. WT
resGene <- lapply(geneDds,function(g) results(g,contrast = c('Genotype','Pten m3m4/m3m4','Pten wt/wt'), alpha=0.05))

sapply(resGene,function(g) sum(g$padj < 0.05, na.rm=TRUE))
# P14 Adult 
# 205  2818

outGeneList <- lapply(resGene,function(g) {
  outGene <- as.data.frame(g)
  outGene = outGene[order(outGene$padj,outGene$pvalue),]
  outGene = cbind(outGene,geneMap[rownames(outGene),])
})

sigGeneList = lapply(outGeneList,function(g){
  sigGene = g[which(g$pvalue<.01),]
})

pdf('plots/DESeq2_MA_plots_pten.pdf')
plotMA(resGene[['P14']], main="Gene DEG MA plot at P14", ylim=c(-1,1))
plotMA(resGene[['Adult']], main="Gene DEG MA plot at Adult", ylim=c(-2,2))
dev.off()

#################
# save everything
library(WriteXLS)
WriteXLS(sigGeneList, ExcelFileName = 'tables/pten_DE_ages_table_DESeq2.xls',row.names=T)
save(sigGeneList,outGeneList,file = 'rdas/pten_ages_DE_objects_DESeq2.rda')
save(geneDds,resGene, file = '/dcl01/lieber/ajaffe/Brady/mouseRNAseq/TilotPTEN/pten_ages_DESeq2_svaAdj.rda')

