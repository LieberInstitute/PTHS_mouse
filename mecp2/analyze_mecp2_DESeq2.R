### mecp2 mouse differential expression analysis 
# qsub -V -l mf=10G,h_vmem=15G,h_stack=256M -cwd -b y R CMD BATCH analyze_mecp2.R
source('../DESeq2_functions.R') #work-horse of differential expression

##############################################
# load phenotype data and RPKM expression data
load('./rdas/pheno.rda',envir = dat<-new.env())
load('/dcl01/lieber/ajaffe/Brady/mouseRNAseq/Rett/rawCounts_mecp2_OCT20_n6.rda')
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
geneDds <- DESeq2(countData = geneCounts, colData = pd, design = ~Genotype,sva = FALSE,parallel=TRUE)
exonDds <- DESeq2(countData = exonCounts, colData = pd, design = ~Genotype,sva = FALSE,parallel=TRUE)
jxnDds <- DESeq2(countData = jCounts, colData = pd, design = ~Genotype,sva = FALSE,parallel=TRUE)

################################################################
# get DE results, and fold-change homozygous mutant v. wild-type
resGene <- results(geneDds,contrast = c('Genotype','MeCP2 KO','wild type'),alpha=0.05) 
resExon <- results(exonDds,contrast = c('Genotype','MeCP2 KO','wild type'),alpha=0.05) 
resJxn <- results(jxnDds,contrast = c('Genotype','MeCP2 KO','wild type'),alpha=0.05) 

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

pdf('plots/DESeq2_MA_plots_mecp2.pdf')
plotMA(resGene, main="Gene MA plot", ylim=c(-.75,.75))
plotMA(resExon, main="Exon MA plot", ylim=c(-1,1))
plotMA(resJxn, main="Junction MA plot", ylim=c(-2,2))
dev.off()

#################
# save everything
library(WriteXLS)
WriteXLS(list(Gene = sigGene,Exon = sigExon,Junction = sigJxn), ExcelFileName = 'tables/mecp2_DE_table_DESeq2.xls',row.names=T)
save(outGene,outExon,outJxn,file = 'rdas/mecp2_DE_objects_DESeq2.rda')
save(geneDds,exonDds,jxnDds, file = '/dcl01/lieber/ajaffe/Brady/mouseRNAseq/Rett/mecp2_DESeq2.rda')
