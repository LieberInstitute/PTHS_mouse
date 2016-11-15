## mega analysis pooling Maher, Philpot, and Sweatt RNAseq datas using DESeq2
source('../DESeq2_functions.R') #work-horse of differential expression

##############################################
# load phenotype data and RPKM expression data
load('/dcl01/lieber/ajaffe/Brady/mouseRNAseq/rawCounts_mega_dataset_nov14_n110.rda')
pd$Age[pd$Age=='P1'] = 'p1'
pd$Age=droplevels(pd$Age)

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

##############################
# create and run DESeq objects
geneDds <- DESeq2(countData = geneCounts, colData = pd, design = ~Genotype+Age+Line+totalAssignedGene,sva = TRUE,parallel=TRUE)
rm(geneCounts); gc()

exonDds <- DESeq2(countData = exonCounts, colData = pd, design = ~Genotype+Age+Line+totalAssignedGene,sva = TRUE,parallel=TRUE)
rm(exonCounts); gc()

jxnDds <- DESeq2(countData = jCounts, colData = pd, design = ~Genotype+Age+Line+totalAssignedGene,sva = TRUE,parallel=TRUE)
rm(jCounts); gc()

############################################
# get DE results, and fold-change PTHS v. WT
resGene <- results(geneDds,contrast = c('Genotype','HT','WT'),alpha=0.05) 
resExon <- results(exonDds,contrast = c('Genotype','HT','WT'),alpha=0.05) 
resJxn <- results(jxnDds,contrast = c('Genotype','HT','WT'),alpha=0.05) 

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

pdf('plots/DESeq2_MA_plots_mega_dataset.pdf')
plotMA(resGene, main="Gene MA plot", ylim=c(-.75,.75))
plotMA(resExon, main="Exon MA plot", ylim=c(-.75,.75))
plotMA(resJxn, main="Junction MA plot", ylim=c(-.75,.75))
dev.off()

#################
# save everything
library(WriteXLS)
WriteXLS(list(Gene = sigGene,Exon = sigExon,Junction = sigJxn), ExcelFileName = 'tables/stable5_mega_dataset_DE_table_DESeq2.xls',row.names=T)
save(outGene,outExon,outJxn,file = 'rdas/mega_dataset_DE_objects_DESeq2.rda')
save(geneDds,exonDds,jxnDds, file = '/dcl01/lieber/ajaffe/Brady/mouseRNAseq/mega_dataset_DESeq2_svaAdj.rda')

#################
# make some plots
cleaningY = function(y, mod, P=ncol(mod)) {
  Hat=solve(t(mod)%*%mod)%*%t(mod)
  beta=(Hat%*%t(y))
  cleany=y-t(as.matrix(mod[,-c(1:P)])%*%beta[-c(1:P),])
  return(cleany)
}
library(beeswarm)
yGene = log(counts(geneDds,normalized = T)+1)
modGene = model.matrix(design(geneDds),data = colData(geneDds))
geneClean = 2^cleaningY(yGene[rownames(sigGene),],mod = modGene,P=2)-1

pdf('plots/mega_dataset_deg_plots.pdf',width = 9)
par(mar = c(7,4,2,2))
for (i in rownames(sigGene)[1:100]){
  boxplot(geneClean[i,]~Genotype*Line,data = colData(geneDds),
          main = sigGene[i,c('Symbol')],las = 2,
          ylab = 'Normalized Counts | SVs')
  beeswarm(geneClean[i,]~Genotype*Line,data = colData(geneDds),
             pch = 21,pwbg = factor(colData(geneDds)$Genotype),add = T)
  }
dev.off()