## separately analyze Maher, Philpot, and Sweatt RNAseq datas using DESeq2
source('../DESeq2_functions.R') #work-horse of differential expression

##############################################
# load phenotype data and RPKM expression data
load('/dcl01/lieber/ajaffe/Brady/mouseRNAseq/rawCounts_mega_dataset_nov14_n110.rda')
rownames(pd) = pd$SAMPLE_ID

########################################
# remove outliers mouse 31 and p21 mice,
ind = which(pd$Age !='p21' &rownames(pd) !='Sample_3_Adult_HT_090614_5_C76VFACXX')
pd = pd[ind,]
geneCounts = geneCounts[,ind]
pd$Age[pd$Age=='P1'] = 'p1'
pd$Age = droplevels(pd$Age)
with(pd,table(Age,Line))

#####################
# split by age groups
indList = split(seq(nrow(pd)), paste0(pd$Line,'.',pd$Age))
names(indList)

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

pdf('plots/DESeq2_MA_plots_mega_tcf4_separate.pdf')
  for(n in names(resGene)){
    plotMA(resGene[[n]],main = paste(n,'MA Plot'),ylim = c(-1,1))
  }
dev.off()

#################
# save everything
library(WriteXLS)
WriteXLS(sigGeneList, ExcelFileName = 'tables/mega_tcf4_separate_DE_table_DESeq2.xls',row.names=T)
save(sigGeneList,outGeneList,file = 'rdas/mega_tcf4_separate_DE_objects_DESeq2.rda')
save(geneDds,resGene, file = '/dcl01/lieber/ajaffe/Brady/mouseRNAseq/mega_tcf4_separate_DESeq2_svaAdj.rda')

#################
# make some plots
indP1 = grep('p1',names(outGeneList),value = T)
indAdult = grep('Adult',names(outGeneList),value = T)

library(gplots)
pdf('plots/venn_separate_mega_tcf4_degs.pdf')
geneVenn = plot(venn(lapply(outGeneList[indP1],function(g) rownames(g)[which(g$pvalue<.01 & !is.na(g$padj))]),small =1.5,show.plot = F),cex = 2)
geneVenn = plot(venn(lapply(outGeneList[indAdult[-6]],function(g) rownames(g)[which(g$pvalue<.01 & !is.na(g$padj))]),small =1.5,show.plot = F),cex = 2)
dev.off()

#############################
## replication rate by model
sapply(indP1,function(i) {
  x = with(outGeneList[[i]],rownames(outGeneList[[i]])[which(padj<.05 & !is.na(padj))])
  y = unique(unlist(lapply(outGeneList[indP1[!indP1 %in% i]],
                           function(g) rownames(g)[g$pvalue< 0.01])))
    sum(x %in% y)/length(x)*100})
#  Del.p1 Maher.p1  Nest.p1 R579W.p1 
#     13       10       17       15       

sapply(indAdult,function(i) {
  x = with(outGeneList[[i]],rownames(outGeneList[[i]])[which(padj<.05 & !is.na(padj))])
  y = unique(unlist(lapply(outGeneList[indAdult[!indAdult %in% i]],
                           function(g) rownames(g)[g$pvalue< 0.01])))
  sum(x %in% y)/length(x)*100})

#Act.Adult    Del.Adult  Maher.Adult   Nest.Adult  R579W.Adult Sweatt.Adult 
#70           86           47          57          81           29      

#######################
## replication heat map
library(Heatplus)
library(jaffelab)
library(RColorBrewer)


#########################################################
# plot replication p1 gene heat maps, in 2 of 6 models
p1Genes = unique(unlist(lapply(outGeneList[indP1],function(g) rownames(g)[which(g$padj<.05 & !is.na(g$padj))])))
p1Genes = p1Genes[sapply(p1Genes,function(g) sum(sapply(indP1,function(i) g %in% rownames(sigGeneList[[i]]))))>1]
p1Dat = do.call('cbind',lapply(outGeneList[indP1], function(g) g[p1Genes,c('log2FoldChange')]))
p1Dat[is.na(p1Dat)] = 0
colnames(p1Dat) = ss(colnames(p1Dat),'\\.')

pdf('plots/replication_heatmap_mega_tcf4_p1.pdf',height = 4,width = 4)
pheatmap(p1Dat,colorRampPalette(brewer.pal(n = 11, name ="PiYG"))(length(seq(-.6,.6,0.1))),
                    border_color = NA,ylab = 'Log2 Fold-change',legend = F,
         labels_row = outGeneList[[1]][p1Genes,'Symbol'],breaks = seq(-.6,.6,0.1))
dev.off()


#########################################################
# plot replication adult gene heat maps, in 4 of 6 models
adultGenes = unique(unlist(lapply(outGeneList[indAdult],function(g) rownames(g)[which(g$padj<.05 & !is.na(g$padj))])))
adultGenes = adultGenes[sapply(adultGenes,function(g) sum(sapply(indAdult,function(i) g %in% rownames(sigGeneList[[i]]))))>4]
adultDat = do.call('cbind',lapply(outGeneList[indAdult], function(g) g[adultGenes,c('log2FoldChange')]))
adultDat[is.na(adultDat)] = 0
colnames(adultDat) = ss(colnames(adultDat),'\\.')
maxFC = max(abs(adultDat))+.01

pdf('plots/replication_heatmap_mega_tcf4_adult.pdf',height = 6,width = 4)
pheatmap(adultDat,colorRampPalette(brewer.pal(n = 11, name ="PiYG"))(length(seq(-maxFC,maxFC,0.1))),
         legend = F,breaks = length(seq(-maxFC,maxFC,0.1)),
         border_color = NA,ylab = 'Log2 Fold-change',labels_row = outGeneList[[1]][adultGenes,'Symbol'])
dev.off()




