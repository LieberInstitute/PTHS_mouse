# fold-change differential expression in adult chd8 mice v. mega_tcf4 mice
source('../DESeq2_functions.R') #work-horse of differential expression
library(DESeq2)
library(jaffelab)
library(ggplot2)
###############
# load the data
load('/dcl01/lieber/ajaffe/Brady/chd8/rawCounts_chd8_feb7_n12.rda')
load('/dcl01/lieber/ajaffe/Brady/chd8/rpkmCounts_chd8_feb7_n12.rda')
load('rdas/pheno_chd8.rda',envir = dat <-new.env())
pd = data.frame(dat$pd,metrics[dat$pd$Sample_s,],row.names = dat$pd$Sample_s)
rownames(geneCounts) = ss(rownames(geneCounts),'\\.')
rownames(geneRpkm) = ss(rownames(geneRpkm),'\\.')
rownames(geneMap) = ss(rownames(geneMap),'\\.')

##########################
# find logFC in adult mice
#ind = which(pd$Age=='Adult')
#pd = pd[ind,]
#geneCounts = geneCounts[,ind]

######################
# exploratory analysis
gExpr = log2(geneRpkm+1)
gExpr = gExpr[rowMeans(geneRpkm)>0.1,]

pca = prcomp(t(gExpr))
ggplot(data = cbind(pca$x,pd),aes(x = PC1,y = PC2,colour = Genotype,pch = Age))+geom_point()
ggplot(data = cbind(pca$x,pd),aes(x = Age,y = PC1,colour = Genotype))+geom_point()
ggplot(data = cbind(pca$x,pd),aes(x = Age,y = PC2,colour = Genotype))+geom_point()


#########################
# DE analysis by genotype
geneDds <- DESeq2(countData = geneCounts, colData = pd, 
                  design = ~Genotype+Age,sva=T,parallel=T)
resGene <-results(geneDds,contrast = c('Genotype','CHD8hetero','WildType'),alpha = 0.05)
sum(resGene$padj < 0.05, na.rm=TRUE)

outGene <- as.data.frame(resGene)
outGene = outGene[order(outGene$padj,outGene$pvalue),]
outGene = cbind(outGene,geneMap[rownames(outGene),])
sigGene = outGene[which(outGene$pvalue<.01),]

#################
# save everything
library(WriteXLS)
WriteXLS(sigGene, ExcelFileName = 'tables/chd8_DE_table_DESeq2.xls',row.names=T)
save(sigGene,outGene,file = 'rdas/chd8_DE_objects_DESeq2.rda')
save(geneDds,resGene, file = '/dcl01/lieber/ajaffe/Brady/chd8/chd8_DESeq2_svaAdj.rda')

#########################
# compare with mega_tcf4 ages DEGs
load('../tcf4_mouse/rdas/mega_tcf4_ages_DE_objects_DESeq2.rda')
lfc = log2(nCounts[,pd$Genotype=='CHD8hetero']/nCounts[,pd$Genotype=='WildType'])
lapply(sigGeneList,function(g){
  tmp =outGene[rownames(g),]
  ind = tmp$pvalue < 0.05
  (t1 = with(g[ind,], table(Tcf4 = log2FoldChange>0,Chd8 = tmp[ind,'log2FoldChange']>0)))
  print(t1)
  fisher.test(t1)
})


