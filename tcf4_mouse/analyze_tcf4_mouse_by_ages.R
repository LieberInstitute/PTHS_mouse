### Brady's pths mouse differential expression analysis using DESeq2
# qsub -V -l mf=200G,h_vmem=250G,h_stack=256M -cwd -b y R CMD BATCH analyze_tcf4_mouse_DESeq2.R
source('../DESeq2_functions.R') #work-horse of differential expression

##############################################
# load phenotype data and RPKM expression data
load('./rdas/pheno.rda',envir = dat<-new.env())
load('/dcl01/lieber/ajaffe/Brady/mouseRNAseq/rawCounts_tcf4_mouse_OCT20_n36.rda')
pd = cbind(dat$pd,pd)
all.equal(pd$FileID,pd$SAMPLE_ID) #samples line up

##########################
# conform to DESeq2 data types
jCounts = as.matrix(as.data.frame(jCounts))
jIndex=which(jMap$code != "Novel")
jCounts = jCounts[jIndex,]
jMap = jMap[jIndex]
colnames(jCounts) = pd$SAMPLE_ID
rownames(pd) = pd$SAMPLE_ID

##########################
# remove outliers mouse 31
pd = pd[-c(31),]
geneCounts = geneCounts[,-c(31)]
exonCounts = exonCounts[,-c(31)]
jCounts = jCounts[,-c(31)]

#####################
# split by age groups
indList = split(seq(nrow(pd)), pd$Age)

##############################
# create and run DESeq objects
geneDds <- lapply(indList,function(i) {
  DESeq2(countData = geneCounts[,i], colData = pd[i,], 
         design = ~Genotype,sva = TRUE,parallel=TRUE)})

############################################
# get DE results, and fold-change PTHS v. WT
resGene <- lapply(geneDds,function(g) results(g,contrast = c('Genotype','HT','WT'), alpha=0.05))

sapply(resGene,function(g) sum(g$padj < 0.05, na.rm=TRUE))

outGeneList <- lapply(resGene,function(g) {
  outGene <- as.data.frame(g)
  outGene = outGene[order(outGene$padj,outGene$pvalue),]
  outGene = cbind(outGene,geneMap[rownames(outGene),])
})

sigGeneList = lapply(outGeneList,function(g){
  sigGene = g[which(g$pvalue<.01),]
})

pdf('plots/DESeq2_MA_plots_tcf4_ages.pdf')
plotMA(resGene[[1]], main="P1 Gene MA plot", ylim=c(-2,2))
plotMA(resGene[[2]], main="P21 Gene MA plot", ylim=c(-2,2))
plotMA(resGene[[3]], main="Adult Gene MA plot", ylim=c(-2,2))
dev.off()

#################
# save everything
library(WriteXLS)
WriteXLS(sigGeneList, ExcelFileName = 'tables/mouse_tcf4_ages_DE_table_DESeq2.xls',row.names=T)
save(outGene,file = 'rdas/mouse_tcf4_ages_DE_objects_DESeq2.rda')
save(geneDds, file = '/dcl01/lieber/ajaffe/Brady/mouseRNAseq/mouse_tcf4_ages_DESeq2_svaAdj.rda')

