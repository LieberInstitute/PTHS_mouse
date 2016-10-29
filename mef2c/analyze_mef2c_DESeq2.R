### mef2c mouse differential expression analysis 
# qsub -V -l mf=10G,h_vmem=15G,h_stack=256M -cwd -b y R CMD BATCH analyze_mef2c.R
source('../DESeq2_functions.R') #work-horse of differential expression

##############################################
# load phenotype data and RPKM expression data
load('./rdas/pheno.rda',envir = dat<-new.env())
load('/dcl01/lieber/ajaffe/Brady/mef2c/rawCounts_MEF2C_OCT27_n6.rda')
pd = cbind(dat$pd,pd)
all.equal(pd$SampleID,pd$SAMPLE_ID) #samples line up

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
geneDds <- DESeq2(countData = geneCounts, colData = pd, design = ~Genotype,sva = TRUE,parallel=TRUE)
rm(geneCounts); gc()
exonDds <- DESeq2(countData = exonCounts, colData = pd, design = ~Genotype,sva = TRUE,parallel=TRUE)
rm(exonCounts); gc()
jxnDds <- DESeq2(countData = jCounts, colData = pd, design = ~Genotype,sva = TRUE,parallel=TRUE)
rm(jCounts); gc()

################################################################
# get DE results, and fold-change homozygous mutant v. wild-type
resGene <- results(geneDds,contrast = c('Genotype','Mef2c cKO (Mef2cfl/fl; Emx1-Cre)','control (Mef2cfl/fl)'),alpha=0.05) 
resExon <- results(exonDds,contrast = c('Genotype','Mef2c cKO (Mef2cfl/fl; Emx1-Cre)','control (Mef2cfl/fl)'),alpha=0.05) 
resJxn <- results(jxnDds,contrast = c('Genotype','Mef2c cKO (Mef2cfl/fl; Emx1-Cre)','control (Mef2cfl/fl)'),alpha=0.05) 

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

pdf('plots/DESeq2_MA_plots_mef2c.pdf')
plotMA(resGene, main="Gene MA plot", ylim=c(-4,4))
plotMA(resExon, main="Exon MA plot",  ylim=c(-4,4))
plotMA(resJxn, main="Junction MA plot",  ylim=c(-4,4))
dev.off()

#################
# save everything
library(WriteXLS)
WriteXLS(list(Gene = sigGene,Exon = sigExon,Junction = sigJxn), ExcelFileName = 'tables/mef2c_DE_table_DESeq2.xls',row.names=T)
save(outGene,outExon,outJxn,file = 'rdas/mef2c_DE_objects_DESeq2.rda')
save(geneDds,exonDds,jxnDds, file = '/dcl01/lieber/ajaffe/Brady/mef2c/mef2c_DESeq2_svaAdj.rda')
