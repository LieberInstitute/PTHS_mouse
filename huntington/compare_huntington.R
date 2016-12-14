# compare disease and mouse model DE with TCF4 mutation
library(biomaRt)
library(jaffelab)

#################################
#get Ensembl mouse to human genes
ensembl=useMart("ensembl",dataset = 'mmusculus_gene_ensembl')
MMtoHG = getBM(attributes = c('ensembl_gene_id','hsapiens_homolog_ensembl_gene'),mart = ensembl)

###################
# load all the data
load('rdas/huntington_DE_objects_DESeq2.rda', envir = huntington <-new.env())
load('/dcl01/lieber/ajaffe/psychENCODE_Data/UCLA_R01MH094714/rdas/ASD_DE_gene.rda')
load('../pten/rdas/pten_DE_objects_DESeq2.rda',envir = pten <-new.env())
load('../mecp2/rdas/mecp2_DE_objects_DESeq2.rda',envir = mecp2 <-new.env())
load('../ube3a/rdas/ube3a_DE_objects_DESeq2.rda',envir = ube3a <-new.env())
load('../mef2c/rdas/mef2c_DE_objects_DESeq2.rda',envir = mef2c <-new.env())
load('../parkinson/rdas/parkinson_DE_objects_DESeq2.rda',envir = parkinson <-new.env())
load('../tcf4_mouse/rdas/mega_dataset_DE_objects_DESeq2.rda')

##########################################
# human ENSEMBL gene labels on mouse genes
outGene$hsapien_homolog = MMtoHG$hsapiens_homolog_ensembl_gene[
  match(rownames(outGene),MMtoHG$ensembl_gene_id)]
outGene = outGene[grepl('ENSG',outGene$hsapien_homolog),] #take only genes w/ human homolog
#outGene = with(outGene, outGene[padj < 0.05 & !is.na(padj),]) #take only genes DEGs

##############################
# compare human ASD DE log FC
n = c('Qual_ba41_42_22','Qual_ba9','Qual_vermis')
#ASD_auditory cortex 0.581031   0.030590  18.994   <2e-16 ***
#ASD_DLPFC 0.724040   0.039403  18.375  < 2e-16 ***
#ASD_vermis  0.51542    0.11826   4.358 3.16e-05 ***
lapply(n,function(n){
  nFC = paste0(n,'.log2FC'); nP = paste0(n,'.pval')
  out2 = cbind(outGene,out[outGene$hsapien_homolog,c(nFC,nP)])
  out2 = out2[out2[,nP]< 0.05,]
  summary(lm(as.formula(paste('log2FoldChange~',nFC)),data = out2))
})

####################################################
# compare mouse Parkinson, Pten, MeCP2, Ube3a, Mef2c

#Parkinson -0.769902   0.053771 -14.318   <2e-16 ***
summary(lm(outGene$log2FoldChange~log2FoldChange,data = parkinson$outGene[rownames(outGene),],
           subset = pvalue < 0.05))
plot(outGene$log2FoldChange~log2FoldChange,data = parkinson$outGene[rownames(outGene),],
     subset = pvalue < 0.05)

#Pten -0.46366    0.01571 -29.512  < 2e-16 ***
summary(lm(outGene$log2FoldChange~log2FoldChange,data = pten$outGene[rownames(outGene),],
           subset = Symbol != 'Pten' & pvalue < 0.05))
plot(outGene$log2FoldChange~log2FoldChange,data = pten$outGene[rownames(outGene),],
     subset = Symbol != 'Pten' & pvalue < 0.05)

#MeCP2 -0.26915    0.05405  -4.980 1.71e-06 ***
summary(lm(outGene$log2FoldChange~log2FoldChange,data = mecp2$outGene[rownames(outGene),],
           subset = Symbol != 'Mecp2' & pvalue < 0.05))
plot(outGene$log2FoldChange~log2FoldChange,data = mecp2$outGene[rownames(outGene),],
     subset = Symbol != 'Mecp2' & pvalue < 0.05)

#Ube3a -0.08676    0.16220  -0.535    0.596
summary(lm(outGene$log2FoldChange~log2FoldChange,data = ube3a$outGene[rownames(outGene),],
           subset = Symbol != 'Ube3a' & pvalue < 0.05))
plot(outGene$log2FoldChange~log2FoldChange,data = ube3a$outGene[rownames(outGene),],
     subset = Symbol != 'Ube3a' & pvalue < 0.05)

#Mef2c -0.042360   0.013562  -3.124 0.001915 ** 
summary(lm(outGene$log2FoldChange~log2FoldChange,data = mef2c$outGene[rownames(outGene),],
           subset = Symbol != 'Mef2c' & pvalue < 0.05))
plot(outGene$log2FoldChange~log2FoldChange,data = mef2c$outGene[rownames(outGene),],
     subset = Symbol != 'Mef2c' & pvalue < 0.05)

##############################
# compare Huntington DE log FC
#Huntington -1.487e-01  2.304e-02  -6.457 2.85e-10 ***
summary(lm(outGene$log2FoldChange~log2FoldChange,data = huntington$outGene[outGene$hsapien_homolog,],
           subset = pvalue < 0.05))
plot(outGene$log2FoldChange~log2FoldChange,data = huntington$outGene[outGene$hsapien_homolog,],
     subset = pvalue < 0.05)


#######################
# make pca of disorders
lfc = cbind(TCF4 = outGene$log2FoldChange,
            MECP2 = mecp2$outGene[rownames(outGene),'log2FoldChange'],
            PTEN = pten$outGene[rownames(outGene),'log2FoldChange'],
            PARK = parkinson$outGene[rownames(outGene),'log2FoldChange'],
            HUNT = huntington$outGene[outGene$hsapien_homolog,'log2FoldChange'],
            #MEF2C = mef2c$outGene[rownames(outGene),'log2FoldChange'],
            UBE3A = ube3a$outGene[rownames(outGene),'log2FoldChange'],
            Auditory = out[outGene$hsapien_homolog,'Qual_ba41_42_22.log2FC'],
            DLPFC = out[outGene$hsapien_homolog,'Qual_ba9.log2FC'],
            Vermis = out[outGene$hsapien_homolog,'Qual_vermis.log2FC'])
pval = cbind(TCF4 = outGene$pvalue,
            MECP2 = mecp2$outGene[rownames(outGene),'pvalue'],
            PTEN = pten$outGene[rownames(outGene),'pvalue'],
            PARK = parkinson$outGene[rownames(outGene),'pvalue'],
            HUNT = huntington$outGene[outGene$hsapien_homolog,'pvalue'],
            #MEF2C = mef2c$outGene[rownames(outGene),'pvalue'],
            UBE3A = ube3a$outGene[rownames(outGene),'pvalue'],
            Auditory = out[outGene$hsapien_homolog,'Qual_ba41_42_22.pval'],
            DLPFC = out[outGene$hsapien_homolog,'Qual_ba9.pval'],
            Vermis = out[outGene$hsapien_homolog,'Qual_vermis.pval'])
ind = apply(pval,1,function(p) any (p < 0.01 & !is.na(p)))
labs = colnames(lfc)
lfc = lfc[complete.cases(lfc) & ind,]
pca = prcomp(t(lfc))
varexp = (pca$sdev)^2 / sum(pca$sdev^2)

plot(pca$x, pch = '',xlab = paste('PC1',signif(varexp[1],2)*100,'% Var Explained'),
     ylab = paste('PC2',signif(varexp[2],2)*100,'% Var Explained'))
text(pca$x[,1],pca$x[,2],labs, pos = 3,cex = .75)



library(pheatmap)
library(RColorBrewer)
pheatmap(lfc,color = colorRampPalette(rev(brewer.pal(n = 11, name = "PiYG")))(100),cluster_rows = TRUE)