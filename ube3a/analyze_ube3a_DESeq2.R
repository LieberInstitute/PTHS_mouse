### ube3a mice differential expression analysis 
# qsub -V -l mf=200G,h_vmem=220G,h_stack=256M -cwd -b y R CMD BATCH analyze_ube3a_DESeq2.R
source('../DESeq2_functions.R') #work-horse of differential expression

##############################################
# load phenotype data and RPKM expression data
load('./rdas/pheno.rda',envir = dat<-new.env())
load('/dcl01/lieber/ajaffe/Brady/ube3a/rawCounts_ube3a_nov17_n8.rda')
pd = cbind(dat$pd,pd)
all(seq(nrow(pd))==unlist(sapply(pd$SAMPLE_ID,grep,pd$'BGI file name'))) #samples line up
rownames(pd) = pd$'File rename/Unique identifier'

##########################
# conform to DESeq2 data types
jCounts = as.matrix(as.data.frame(jCounts))
jIndex=which(jMap$code != "Novel")
jCounts = jCounts[jIndex,]
jMap = jMap[jIndex]
colnames(jCounts) = pd$SAMPLE_ID
rownames(pd) = pd$SAMPLE_ID

##############################
# create and run DESeq objects
geneDds <- DESeq2(countData = geneCounts, colData = pd, design = ~Genotype+Age,sva = TRUE,parallel=TRUE)
exonDds <- DESeq2(countData = exonCounts, colData = pd, design = ~Genotype+Age,sva = TRUE,parallel=TRUE)
jxnDds <- DESeq2(countData = jCounts, colData = pd, design = ~Genotype+Age,sva = TRUE,parallel=TRUE)

############################################################
# get DE results, and fold-change of each het mouse genotype
resGene <- results(geneDds,contrast = c('Genotype','Ube3am-/p+','Ube3a+/+'),alpha=0.05) 
resExon <- results(exonDds,contrast = c('Genotype','Ube3am-/p+','Ube3a+/+'),alpha=0.05) 
resJxn <- results(jxnDds,contrast = c('Genotype','Ube3am-/p+','Ube3a+/+'),alpha=0.05) 

sum(resGene$padj < 0.05, na.rm=TRUE)
sum(resExon$padj < 0.05, na.rm=TRUE)
sum(resJxn$padj < 0.05, na.rm=TRUE)

outGene <- as.data.frame(resGene[order(resGene$padj,resGene$pvalue),])
outGene = cbind(outGene,geneMap[rownames(outGene),])
sigGene = outGene[which(outGene$pvalue<.01),]

outExon <- as.data.frame(resExon[order(resExon$padj,resExon$pvalue),])
outExon = cbind(outExon,exonMap[rownames(outExon),])
sigExon = outExon[which(outExon$pvalue<.01),]

outJxn <- as.data.frame(resJxn[order(resJxn$padj,resJxn$pvalue),])
outJxn = cbind(outJxn, as.data.frame(jMap)[rownames(outJxn),])
sigJxn = outJxn[which(outJxn$pvalue<.01),]

pdf('plots/DESeq2_MA_plots_ube3a.pdf')
plotMA(resGene, main="Gene MA plot", ylim=c(-4,4))
plotMA(resExon, main="Exon MA plot",  ylim=c(-4,4))
plotMA(resJxn, main="Junction MA plot",  ylim=c(-4,4))
dev.off()


#######################################
# save all the differential expressions
library(WriteXLS)
WriteXLS(list(Gene = sigGene,Exon= sigExon, Junction = sigJxn, phenotype = pd),
         ExcelFileName = 'tables/ube3a_DE_table.xls')
save(pd, outGene, outExon, outJxn, file="rdas/ube3a_DE_objects_DESeq2.rda")
save(geneDds,exonDds,jxnDds, file = '/dcl01/lieber/ajaffe/Brady/ube3a/ube3a_DESeq2_svaAdj.rda')

