# compare tcf4 DE t-statistic from mecp2 and pten DE t-statistics
getOR = function(x) x[1,1]/x[2,1]/x[1,2]*x[2,2]
getKappa = function(x) (x[1,1]+x[2,2])/sum(x)

library(jaffelab)
library(reshape2)
library(ggplot2)
options(stringsAsFactors = F)

######################################
# load TCF4, MECP2, and PTEN mice data
load('tcf4_mouse/rdas/mega_dataset_DE_objects_DESeq2.rda',envir = tcf4 <- new.env())
load('mecp2/rdas/mecp2_DE_objects_DESeq2.rda',envir = mecp2 <- new.env())
load('pten/rdas/pten_DE_objects_DESeq2.rda',envir = pten <- new.env())

# take only DEG in mega TCF4 data
geneMap = with(tcf4$outGene,tcf4$outGene[padj<0.05 & !is.na(padj),7:13])
labs = rownames(geneMap)

###############################################
# compare TCF4 vs. MECP2 log2FC anti-enrichment
outGene1 = cbind(geneMap,TCF4_log2FC= tcf4$outGene[labs,2], mecp2$outGene[labs,c(2,5)])
outGene1 = with(outGene1,outGene1[pvalue < 0.05 & Symbol!= 'Mecp2',])
outGene1 = outGene1[complete.cases(outGene1),]
dim(outGene1) #154 out of 1370 genes DE in MeCP2 

# test for negative enrichment
(t1 = with(outGene1,table(TCF4=sign(TCF4_log2FC),MECP2= sign(log2FoldChange))))
with(outGene1,cor.test(x = TCF4_log2FC,y = log2FoldChange,method = 'spearman'))
1-getKappa(t1) 


##############################################
# compare PTEN vs. TCF4 logFC (anti) enrichment
outGene2 = cbind(geneMap,TCF4_log2FC= tcf4$outGene[labs,2],pten$outGene[labs,c(2,5)])
outGene2 = with(outGene2,outGene2[pvalue < 0.05,])
outGene2 = outGene2[complete.cases(outGene2),]
dim(outGene2) #666 out of 1370 genes DE in PTEN

# test for negative enrichment
(t2 = with(outGene2,table(TCF4=sign(TCF4_log2FC),PTEN= sign(log2FoldChange))))
with(outGene2,cor.test(x = TCF4_log2FC,y = log2FoldChange,method = 'spearman'))
1-getKappa(t2) #negative enrichment

#####################################
# plot log2fc TCF4 vs. PTEN and MECP2
datList = list(MECP2=outGene1,PTEN = outGene2)
WriteXLS::WriteXLS(datList,ExcelFileName = 'tcf4_mouse/tables/stable8_mecp2_pten_comparison.xls')
dat = do.call('rbind',datList)
dat$Mutation = ss(rownames(dat),'\\.')

pdf('tcf4_mouse/plots/fig5_asd_mice_compare_logFC.pdf',height = 5,width = 9)
ggplot(data = dat,aes(x = TCF4_log2FC,y = log2FoldChange,co))+
  geom_point(alpha = .25) + facet_wrap(~Mutation)+
  ggtitle('')+ xlab('Tcf4 log2 fold change') + ylab('Mouse Mutation log2 fold change')+
  geom_vline(aes(xintercept = 0,colour= 'red',linetype = 'dashed')) + 
  geom_hline(aes(yintercept = 0,colour= 'red',linetype = 'dashed'))+
  theme(strip.background = element_rect(fill = "white", colour = "white"),
        axis.title = element_text(family="Helvetica", face="bold", size=14),
        axis.text = element_text(family="Helvetica", face="bold", size=12),
        legend.position="none")
dev.off()