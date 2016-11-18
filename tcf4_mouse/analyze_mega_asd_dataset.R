## mega analysis pooling Tcf4, Mef2c, MeCP2, Pten, Ube3a  RNAseq datas using DESeq2
source('../DESeq2_functions.R') #work-horse of differential expression

##############################################
# load phenotype data and RPKM expression data
load('/dcl01/lieber/ajaffe/Brady/mouseRNAseq/rawCounts_mega_dataset_asd_mice_nov17_n141.rda')
colnames(geneCounts) = pd$SAMPLE_ID
rownames(pd)= pd$SAMPLE_ID
pd$Age = as.character(pd$Age)

##############################
# create and run DESeq objects
geneDds <- DESeq2(countData = geneCounts, colData = pd, design = ~Genotype+Age+Line+totalAssignedGene,sva = TRUE,parallel=TRUE)

############################################
# get DE results, and fold-change PTHS v. WT
resGene <- results(geneDds,contrast = c('Genotype','MT','WT'),alpha=0.05) 

sum(resGene$padj < 0.05, na.rm=TRUE)

outGene <- as.data.frame(resGene[order(resGene$padj,resGene$pvalue),])
outGene = cbind(outGene,geneMap[rownames(outGene),])
sigGene = outGene[which(outGene$pvalue<.01),]

pdf('plots/DESeq2_MA_plots_mega_asd_mice.pdf')
plotMA(resGene, main="Gene MA plot", ylim=c(-.75,.75))
dev.off()

#################
# save everything
library(WriteXLS)
WriteXLS(list(Gene = sigGene,phenotype = pd), ExcelFileName = 'tables/stable8_mega_asd_mice_DE_table_DESeq2.xls',row.names=T)
save(outGene,file = 'rdas/mega_asd_mice_DE_objects_DESeq2.rda')
save(geneDds, file = '/dcl01/lieber/ajaffe/Brady/mouseRNAseq/mega_asd_mice_DESeq2_svaAdj.rda')

#################
# make some plots
cleaningY = function(y, mod, P=ncol(mod)) {
  Hat=solve(t(mod)%*%mod)%*%t(mod)
  beta=(Hat%*%t(y))
  cleany=y-t(as.matrix(mod[,-c(1:P)])%*%beta[-c(1:P),])
  return(cleany)
}
library(beeswarm)
library(RColorBrewer)

yGene = log(counts(geneDds,normalized = T)+1)
pd = colData(geneDds)
pd$Genotype = factor(pd$Genotype,levels = c('WT',"MT"))
pd$Line = factor(pd$Line,levels = c('Maher','Act','Nest', 'Del', 'R579W', 'Sweatt','Mecp2','Mef2c', 'Pten', 'Ube3a'))
modGene = model.matrix(design(geneDds),data = pd)
geneClean = 2^cleaningY(yGene[rownames(sigGene),],mod = modGene,P=2)-1


n = length(unique(colData(geneDds)$Line))
col = brewer.pal(n,'Paired')
pdf('plots/mega_asd_mice_deg_plots.pdf',width = 9)
par(mar = c(7,4,2,2))
for (i in rownames(sigGene)[1:100]){
  boxplot(geneClean[i,]~Genotype*Line,data = pd,
          main = sigGene[i,c('Symbol')],las = 2,
          ylab = 'Normalized Counts | SVs',col = rep(col,each = 2))
  beeswarm(geneClean[i,]~Genotype*Line,data = pd,
             pch = 21,pwbg = as.numeric(pd$Genotype)-1,add = T)
  }
dev.off()
