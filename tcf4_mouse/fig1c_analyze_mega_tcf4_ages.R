## mega analysis pooling Maher, Philpot, and Sweatt RNAseq datas using DESeq2
source('../DESeq2_functions.R') #work-horse of differential expression
library(ggplot2)
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
indList = split(seq(nrow(pd)), pd$Age)

##############################
# create and run DESeq objects
geneDds <- lapply(indList,function(i) {
  DESeq2(countData = geneCounts[,i], colData = pd[i,], 
         design = ~Genotype+Line+totalAssignedGene,sva = TRUE,parallel=TRUE)})

############################################
# get DE results, and fold-change PTHS v. WT
resGene <- lapply(geneDds,function(g) DESeq2::results(g,contrast = c('Genotype','HT','WT'), alpha=0.05))

sapply(resGene,function(g) sum(g$padj < 0.05, na.rm=TRUE))

outGeneList <- lapply(resGene,function(g) {
  outGene <- as.data.frame(g)
  outGene = outGene[order(outGene$padj,outGene$pvalue),]
  outGene = cbind(outGene,geneMap[rownames(outGene),])
})

sigGeneList = lapply(outGeneList,function(g){
  sigGene = g[which(g$pvalue<.01),]
})

postscript('plots/DESeq2_MA_plots_mega_tcf4_ages.eps',height = 4,width = 3.5)
  plotMA(resGene[[1]], main="P1 Genes MA plot", ylim=c(-.5,.5))
  plotMA(resGene[[2]], main="Adult Genes MA plot", ylim=c(-.5,.5))
dev.off()

#################
# save everything
library(WriteXLS)
WriteXLS(sigGeneList, ExcelFileName = 'tables/mega_tcf4_ages_DE_table_DESeq2.xls',row.names=T)
save(sigGeneList,outGeneList,file = 'rdas/mega_tcf4_ages_DE_objects_DESeq2.rda')
save(geneDds, file = '/dcl01/lieber/ajaffe/Brady/mouseRNAseq/mega_tcf4_ages_DESeq2_svaAdj.rda')

#################
# make some plots
cleaningY = function(y, mod, P=ncol(mod)) {
  Hat=solve(t(mod)%*%mod)%*%t(mod)
  beta=(Hat%*%t(y))
  cleany=y-t(as.matrix(mod[,-c(1:P)])%*%beta[-c(1:P),])
  return(cleany)
}
library(beeswarm)

for(a in c("p1",'Adult')){
  g = geneDds[[a]]
  modGene = model.matrix(design(g),data = colData(g))
  sigGene = sigGeneList[[a]]
  geneClean = cleaningY(counts(g,normalized = T)[rownames(sigGene),],mod = modGene,P=2)
  line = levels(colData(g)$Line)
  
  pdf(paste0('plots/mega_dataset_deg_plots_',a,'.pdf'),width = 6)
  par(mar = c(7,4,2,2))
  for (i in seq(min(100,nrow(sigGene)))){
    beeswarm(geneClean[i,]~Genotype,data = colData(g),
            main = paste(a,'-',sigGene[i,c('Symbol')]),
            ylab = 'Normalized Counts | SVs',pwpch = as.numeric(colData(g)$Line))
    mtext(paste('Fold-change =',signif(2^sigGene[i,c('log2FoldChange')],3),
                '\nFDR =',signif(sigGene[i,c('padj')],3)),1)
    legend('topright',legend = levels(colData(g)$Line),pch = seq(length(levels(colData(g)$Line))))
    }
  dev.off()
}

beeswarm(geneClean[i,]~Genotype+Line,data = colData(g),
         main = paste(a,'-',sigGene[i,c('Symbol')]),pwcol = Genotype,
         ylab = 'Normalized Counts | SVs',pwpch = as.numeric(colData(g)$Line))
boxplot(geneClean[i,]~Genotype+Line,data = colData(g),add = T)
mtext(paste('Fold-change =',signif(2^sigGene[i,c('log2FoldChange')],3),
            '\nFDR =',signif(sigGene[i,c('padj')],3)),1)
legend('topright',legend = levels(colData(g)$Line),pch = seq(length(levels(colData(g)$Line))))

library(gplots)
pdf('plots/venn_ages_mega_tcf4_degs.pdf')
geneVenn = plot(venn(lapply(outGeneList,function(g) rownames(g)[which(g$padj<.05 & !is.na(g$padj))]),small =1.5,show.plot = F),cex = 2)
dev.off()

intersection = Reduce(intersect,lapply(sigGeneList,rownames))
genes = sigGeneList[[1]][intersection,-c(1:6)]
p1 = sigGeneList[[1]][intersection,c(1:6)]
names(p1) = paste0(names(p1),'_p1')
adult = sigGeneList[[2]][intersection,c(1:6)]
names(adult) = paste0(names(adult),'_Adult')
sigGeneList[['Both']]= cbind(genes,p1,adult)

WriteXLS(sigGeneList, ExcelFileName = 'tables/mega_tcf4_ages_DE_table_DESeq2.xls',row.names=T)


################################
# overlap of DEGs at P1 and Adult
labs =  Reduce(intersect,lapply(outGeneList,rownames))
tmp = sapply(outGeneList,function(x) with(x[labs,],padj<0.05 & !is.na(padj)))
(tt = table(P1 = tmp[, 1],Adult = tmp[, 2]))
fisher.test(tt)


########################
# output data for CIBERSORT
write.table(geneCounts,sep = '\t',quote = FALSE,file = 'tables/rawGeneCounts_mega_tcf4.txt')


###################
# alignment details
load('../tcf4_mouse/rdas/mega_tcf4_pheno.rda')

median(pd$totalMapped[pd$Line=='Maher' & pd$Age !='p21']/1e6)
quantile(pd$totalMapped[pd$Line=='Maher' & pd$Age !='p21']/1e6)
median(pd$totalAssignedGene[pd$Line=='Maher' & pd$Age !='p21']*100)
quantile(pd$totalAssignedGene[pd$Line=='Maher' & pd$Age !='p21']*100)

median(pd$totalMapped[! pd$Line %in% c('Maher','Sweatt')]/1e6)
quantile(pd$totalMapped[!pd$Line %in% c('Maher','Sweatt')]/1e6)
median(pd$totalAssignedGene[!pd$Line %in% c('Maher','Sweatt')]*100)
quantile(pd$totalAssignedGene[!pd$Line %in% c('Maher','Sweatt')]*100)