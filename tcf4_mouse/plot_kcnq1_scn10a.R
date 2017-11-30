source('../DESeq2_functions.R') #work-horse of differential expression
library(ggplot2)

##############################################
# load phenotype data and RPKM expression data
load('/dcl01/lieber/ajaffe/Brady/mouseRNAseq/rpkmCounts_tcf4_mouse_OCT20_n36.rda')
load('rdas/pheno.rda',envir = dat<- new.env())
load('rdas/mouse_tcf4_DE_objects_DESeq2.rda')
pd = cbind(dat$pd[match(pd$SAMPLE_ID,dat$pd$FileID),],pd)
rownames(pd) = pd$SAMPLE_ID
with(pd,table(FileID,SAMPLE_ID))

########################################
# remove outliers mouse 31 and p21 mice,
ind = which(rownames(pd) !='Sample_3_Adult_HT_090614_5_C76VFACXX')
pd = pd[ind,]
geneRpkm = geneRpkm[,ind]
pd$Age[pd$Age=='P1'] = 'p1'
pd$Age = droplevels(pd$Age)

#############################################
#
genes = c('Scn10a','Kcnq1')
whichGenes = sapply(paste0('^',genes,'$'),grep,x = geneMap$Symbol)
names(whichGenes) = genes

pdf('plots/scn10a_kcnq1_maher_tcf4_mouse.pdf')
ggplot(data = cbind(pd,Expr = geneRpkm[whichGenes['Scn10a'],]))+ 
  geom_boxplot(aes(x = Age,y = Expr,fill = Genotype))+
  ggtitle('Scn10a Gene Expression') + ylab('FPKM') + xlab('Age')

ggplot(data = cbind(pd,Expr = geneRpkm[whichGenes['Kcnq1'],]))+ 
  geom_boxplot(aes(x = Age,y = Expr,fill = Genotype))+
  ggtitle('Kcnq1 Gene Expression') + ylab('FPKM') + xlab('Age')
dev.off()

(lmScn10a = summary(lm(Expr~Genotype+Age+Genotype:Age,
              data= cbind(pd, Expr = log2(geneRpkm[whichGenes['Scn10a'],]+1)))))
(lmKcnq1 = summary(lm(Expr~Genotype+Age+Genotype:Age,
                      data= cbind(pd, Expr = log2(geneRpkm[whichGenes['Kcnq1'],]+1)))))
