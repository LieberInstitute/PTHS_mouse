### parkinson mice differential expression analysis 
# qsub -V -l mf=200G,h_vmem=220G,h_stack=256M -cwd -b y R CMD BATCH analyze_parkinson_DESeq2.R
source('../DESeq2_functions.R') #work-horse of differential expression

##############################################
# load phenotype data and RPKM expression data
load('./rdas/pheno.rda',envir = dat<-new.env())
load('/dcl01/lieber/ajaffe/Brady/parkinson/rawCounts_parkinson_nov28_n11.rda')
pd = cbind(dat$pd,pd)
all(seq(nrow(pd))==unlist(sapply(pd$SAMPLE_ID,grep,pd$SampleID))) #samples line up
rownames(pd) = pd$SampleID

##############################
# create and run DESeq objects
geneDds <- DESeq2(countData = geneCounts, colData = pd, design = ~Condition,sva = TRUE,parallel=TRUE)

############################################################
# get DE results, and fold-change of each het mouse genotype
resGene <- results(geneDds,contrast = c('Condition','Lesioned','Unlesioned'),alpha=0.05) 

sum(resGene$padj < 0.05, na.rm=TRUE)

outGene <- as.data.frame(resGene[order(resGene$padj,resGene$pvalue),])
outGene = cbind(outGene,geneMap[rownames(outGene),])
sigGene = outGene[which(outGene$pvalue<.01),]

pdf('plots/DESeq2_MA_plots_parkinson.pdf')
plotMA(resGene, main="Gene MA plot", ylim=c(-1,1))
dev.off()

#######################################
# save all the differential expressions
library(WriteXLS)
WriteXLS(list(Gene = sigGene, phenotype = pd), ExcelFileName = 'tables/parkinson_DE_table.xls')
save(pd, outGene, file="rdas/parkinson_DE_objects_DESeq2.rda")
save(geneDds, file = '/dcl01/lieber/ajaffe/Brady/parkinson/parkinson_DESeq2_svaAdj.rda')

