# check TCF4 DEGs with human ASD from Psych ENCODE dataset
getOR = function(x) x[1,1]/x[2,1]/x[1,2]*x[2,2]
getKappa = function(x) (x[1,1]+x[2,2])/sum(x)

library(WriteXLS)
library(biomaRt)
library(data.table)
library(ggplot2)
library(jaffelab)
library(irr)

#get Ensembl mouse to human genes
ensembl=useMart("ensembl",dataset = 'mmusculus_gene_ensembl')
MMtoHG = getBM(attributes = c('ensembl_gene_id','hsapiens_homolog_ensembl_gene'),mart = ensembl)

#load human data, 'out' is the data frame
load('/dcl01/lieber/ajaffe/psychENCODE_Data/UCLA_R01MH094714/rdas/ASD_DE_gene.rda')
indFC = intersect(grep('log2FC',names(out),value = T), grep('Qual',names(out),value = T))
indPval = intersect(grep('pval',names(out),value = T), grep('Qual',names(out),value = T))

##############################
# check mouse tcf4 DEGs w/ psych encode data
load("tcf4_mouse/rdas/mega_dataset_DE_objects_DESeq2.rda") #mouse data
outGene$hsapien_homolog = MMtoHG$hsapiens_homolog_ensembl_gene[match(rownames(outGene),MMtoHG$ensembl_gene_id)]
outGene = outGene[grepl('ENSG',outGene$hsapien_homolog),] #take only genes w/ human homolog

# no enforcement mice or human
genes = outGene
genes$mouseFCsign = paste0(sign(genes$log2FoldChange),'_mm')
outFCrep = sign(out[genes$hsapien_homolog,indFC])
lapply(colnames(outFCrep),function(i) addmargins(table(outFCrep[,i],genes$mouseFCsign)))
sapply(colnames(outFCrep),function(i) fisher.test(outFCrep[,i],genes$mouseFCsign)[1:3])

# human pval < 0.05
outFCrep[out[genes$hsapien_homolog,indPval] > 0.05] = NA
lapply(colnames(outFCrep),function(i) addmargins(table(outFCrep[,i],genes$mouseFCsign)))
sapply(colnames(outFCrep),function(i) fisher.test(outFCrep[,i],genes$mouseFCsign)[1:3])

# mice padjusted < 0.05
genes = outGene[outGene$padj<0.05 & !is.na(outGene$padj),]
#genes = genes[order(genes$fdr_sva),]
genes$mouseFCsign = paste0(sign(genes$log2FoldChange),'_mm')
outFCrep = sign(out[genes$hsapien_homolog,indFC])
lapply(colnames(outFCrep),function(i) addmargins(table(outFCrep[,i],genes$mouseFCsign)))
sapply(colnames(outFCrep),function(i) fisher.test(outFCrep[,i],genes$mouseFCsign)[1:3])

# human pval < 0.05, and mice padjust < 0.05
outFCrep[out[genes$hsapien_homolog,indPval] > 0.05] = NA
lapply(colnames(outFCrep),function(i) addmargins(table(outFCrep[,i],genes$mouseFCsign)))
sapply(colnames(outFCrep),function(i) fisher.test(outFCrep[,i],genes$mouseFCsign)[1:3])

##########################################
# save table of replicated human ASD genes
outList = lapply(seq(length(indFC)),function(x){
  ret = outGene[,c('Symbol','log2FoldChange','padj','hsapien_homolog')]
  ret[,indFC[x]] = out[ret$hsapien_homolog,indFC[x]]
  ret[,indPval[x]] = out[ret$hsapien_homolog,indPval[x]]
  ret = ret[complete.cases(ret),]
  ret$sameFCsign = sign(ret[,indFC[x]]) == sign(ret$log2FoldChange)
  ret = ret[ret[,indPval[x]]<.05 & ret$padj < .05,]
  #ret = ret[ ret$padj < .01,]
  ret = ret[order(apply(log(ret[,c('padj',indPval[x])]),1,mean,na.rm = T)),]
  return(ret)
})
names(outList) = ss(indFC,'\\.')
type=c("BA 41,42,22", "BA 9","Vermis")
dat = rbindlist(outList)
names(dat) = c("Symbol", "log2FoldChange_mouse", "padj_mouse", "hsapien_homolog", "log2FoldChange_human", 
               "pvalue_human", "sameFCsign")
dat$type = rep(paste0(type,' (N=',sapply(outList,nrow),')'),times = sapply(outList,nrow))
WriteXLS(dat,ExcelFileName = 'tables/stable7_psychENCODE_human_asd_replication.xls')

#####################################
# plot human asd FC v tcf4 mouse
pdf('tcf4_mouse/plots/fig6b_human_asd_psychENCODE_fc_tcf4_mouse/plots.pdf',w = 9.5,h = 4)
ggplot(data=dat,aes(x=log2FoldChange_mouse, y=log2FoldChange_human, fill=type)) +
  geom_point(aes(x=log2FoldChange_mouse, y=log2FoldChange_human, fill=type),pch=21)+
  theme_bw(base_size = 14, base_family = "Helvetica")  + facet_wrap(~type,nrow =1,shrink=T) + 
  xlab("Mouse log2 fold-change") + ylab('Human log2 fold change') +
  ylim(c(-1,1)) + xlim(c(-1,1)) +
  geom_vline(aes(xintercept = 0,colour= 'red',linetype = 'dashed')) + 
  geom_hline(aes(yintercept = 0,colour= 'red',linetype = 'dashed'))+ 
  geom_abline(aes(slope = 1,intercept = 0, colour= 'black',linetype = 'dashed'))+
  theme(strip.background = element_rect(fill = "white", colour = "white"),
        legend.position = "none")+
  theme(axis.title.y = element_text(family="Helvetica", face="bold", size=14),
        axis.text.y = element_text(family="Helvetica", face="bold", size=12),
        axis.text.x = element_text(family="Helvetica", face="bold", size=12),
        legend.position="none")
dev.off()

#####################################
# plot human asd FC v tcf4 mouse
lapply(split(dat,dat$type),function(x) fisher.test(sign(x$log2FoldChange_mouse),sign(x$log2FoldChange_human)))
lapply(split(dat,dat$type),function(x) with(x,cor.test(log2FoldChange_mouse,log2FoldChange_human,method = 'spearman')))
lapply(split(dat,dat$type),function(x) getKappa(table(sign(x$log2FoldChange_mouse),sign(x$log2FoldChange_human))))
