library(WriteXLS)
library(readxl)
library(biomaRt)
library(pheatmap)
library(RColorBrewer)
library(jaffelab)
library(reshape2)
library(parallel)
library(ggplot2)
library(beeswarm)
library(nlme)
library(lmerTest)
cleaningY = function(y, mod, P=ncol(mod)) {
  Hat=solve(t(mod)%*%mod)%*%t(mod)
  beta=(Hat%*%t(y))
  cleany=y-t(as.matrix(mod[,-c(1:P)])%*%beta[-c(1:P),])
  return(cleany)
}

reverselog_trans <- function(base = exp(1)) {
  trans <- function(x) -log(x, base)
  inv <- function(x) base^(-x)
  trans_new(paste0("reverselog-", format(base)), trans, inv, 
            log_breaks(base = base), 
            domain = c(1e-100, Inf))
}

#################################
#get Ensembl mouse to human genes
ensembl=useMart("ensembl",dataset = 'mmusculus_gene_ensembl')
MMtoHG = getBM(attributes = c('ensembl_gene_id','hsapiens_homolog_ensembl_gene'),mart = ensembl)

###########################
#load human ASD/Dup15q data
load('asd/rdas/asd_DE_objects_DESeq2_Qual.rda',envir = asd <- new.env())
asd$outGeneList = lapply(asd$outGeneList[grep('Qual',names(asd$outGeneList))],function(g){
  tmpRowNames = ss(rownames(g),'\\.')
  ind = !duplicated(tmpRowNames)
  tmp = g[ind,]; rownames(tmp) = tmpRowNames[ind]
  return(tmp)
})

##########################
# load mouse tcf4 DEG data
load("tcf4_mouse/rdas/mega_tcf4_ages_DE_objects_DESeq2.rda") #PTHS mouse data
outGeneList = lapply(outGeneList, function(g) {
  g$hsapien_homolog = MMtoHG$hsapiens_homolog_ensembl_gene[match(rownames(g),MMtoHG$ensembl_gene_id)]
  g[grepl('ENSG',g$hsapien_homolog),]})
load('mecp2/rdas/mecp2_DE_objects_DESeq2.rda',envir = mecp2 <- new.env()) #Rett mouse data
mecp2$outGene = with(mecp2$outGene, mecp2$outGene[padj<0.05 & !is.na(padj),])
load('pten/rdas/pten_ages_DE_objects_DESeq2.rda',envir = pten <- new.env()) #PTEN mouse data
pten$outGene = with(pten$outGeneList[['Adult']], pten$outGeneList[['Adult']][padj<0.05 & !is.na(padj),])

############################################
# load 34 convergent genes from mouse models
load('tcf4_mouse/rdas/overlap_tcf4_mecp2_pten.rda')
outGene1 = outGeneList[['Adult']][sigGenes,] #myelination genes
outGene1 = outGene1[complete.cases(outGene1),]

outGene2 = outGeneList[['Adult']] #all PTHS mouse DEGs
outGene2 = outGene2[outGene2$padj<0.05 & !is.na(outGene2$padj),]
outGene2 = outGene2[complete.cases(outGene1),]

# genes
outGene3 = outGeneList[['Adult']][!rownames(outGeneList[['Adult']]) %in% rownames(mecp2$outGene) &
                                    !rownames(outGeneList[['Adult']]) %in% rownames(pten$outGene) &
                                    outGeneList[['Adult']]$padj <0.05 & !is.na(outGeneList[['Adult']]$padj),]

# number of 34 genes DE in human ASD/Dup15q
sapply(asd$outGeneList,function(g) with(g,sum(pvalue<0.05 & rownames(g) %in% outGene1$hsapien_homolog)))

######################################
# load in expression RPKM of human ASD
if(FALSE){
  load('/dcl01/lieber/ajaffe/Brady/asd/rpkmCounts_asd_dec9_n204.rda')
  save(geneRpkm, pd, geneMap,file = '/dcl01/lieber/ajaffe/Brady/asd/geneRpkm_asd_dec9_n204.rda')
} else {
  load('/dcl01/lieber/ajaffe/Brady/asd/geneRpkm_asd_dec9_n204.rda')
}
load('asd/rdas/pheno.rda',envir = dat<-new.env())
pd = cbind(pd,dat$pd)
all.equal(pd$SampleID,pd$SAMPLE_ID) #samples line up
rownames(pd) = pd$SAMPLE_ID
pd$Diagnosis = factor(ifelse(pd$Detailed.Diagnosis =="Chromosome 15q Duplication Syndrome",
                             'Dup15',pd$Diagnosis), levels = c("CTL", "ASD",'Dup15'))

########################################
# residualize on confounding and qSVs
# load qSV's
load('asd/rdas/qSVAs-geschwind_asd.rdas')
mod = model.matrix(~Diagnosis + totalAssignedGene + Sequencing.Batch + 
                     Brain.Bank + RIN + Age + Sex, data =pd)
modAdj = cbind(mod,qSVs)

pd = pd[order(pd$Diagnosis),]
#yGene = cleaningY(log2(geneRpkm+1),modAdj,3)
yGene = log2(geneRpkm+1)
yGene = yGene[,rownames(pd)]
geneRpkm = geneRpkm[,rownames(pd)]

######################
# pca on the CAG genes
whichGenes = rownames(geneMap)[match(outGene1$hsapien_homolog,geneMap$ensemblID)]
whichGenes = whichGenes[!is.na(whichGenes)]

#pca1 = prcomp(t(yGene[whichGenes,]))
pca1 = prcomp(t(log2(geneRpkm+1)[whichGenes,]))
# pca1resid = prcomp(t(yGene[whichGenes,]))
## clean the PCs, not expression
pca1resid = t(cleaningY(t(pca1$x), modAdj, P=3))

postscript('asd/plots/asd_myelination_genes_pca.eps',height=5,width=4)
boxplot(pca1resid[,1]~pd$Diagnosis,main = 'Covergent Mouse Gene Set PCA',xlab = '',outline=FALSE,
        ylab= paste0('Eigengene (',signif(summary(pca1)$importance[2,1]*100,2),'% of Variance) | Model'))
col = brewer.pal(3,'Set1')
beeswarm(pca1resid[,1]~pd$Diagnosis,corral = 'wrap',add = T,pch = 20, 
	pwcol = col[as.numeric(factor(pd$Region))])
dev.off()
summary(lm(pca1$x[,1]~Diagnosis+Region+ totalAssignedGene + Sequencing.Batch + 
             Brain.Bank + RIN + Age + Sex + PC1+ PC2+ PC3+ PC4+ PC5,data = cbind(pd,qSVs)))

summary(lmer(pca1$x[,1]~Diagnosis+Region+ totalAssignedGene + Sequencing.Batch + 
             Brain.Bank + RIN + Age + Sex + PC1+ PC2+ PC3+ PC4+ PC5 + (1|Brain.ID),data = cbind(pd,qSVs)))

summary(lm(pca1$x[,1]~Diagnosis + Region+ totalAssignedGene + Sequencing.Batch + 
             Brain.Bank + RIN + Age + Sex,data = cbind(pd,qSVs)))



###############################################
# pca on the all PTHS DEGs in human ASD and 15q
whichGenes = rownames(geneMap)[match(outGene2$hsapien_homolog,geneMap$ensemblID)]
whichGenes = whichGenes[!is.na(whichGenes)]

pca2 = prcomp(t(log2(geneRpkm+1)[whichGenes,]))
# pca2resid = prcomp(t(yGene[whichGenes,]))
## clean the PCs, not expression
pca2resid = t(cleaningY(t(pca2$x), modAdj, P=3))

postscript('asd/plots/pths_mouse_genes_in_asd_pca.eps',height=5,width=4)
boxplot(pca2resid[,1]~pd$Diagnosis,main = 'All PTHS Gene Set PCA',xlab = '',outline = FALSE,
        ylab= paste0('Eigengene (',signif(summary(pca2)$importance[2,1]*100,2),'% of Variance) | Model'))
col = brewer.pal(3,'Set1')
beeswarm(pca2resid[,1]~pd$Diagnosis,corral = 'wrap',
	add = TRUE,pch = 20, pwcol= col[as.numeric(factor(pd$Region))])
legend('right',pch = 20, col = col, legend = toupper(levels(factor(pd$Region))))
dev.off()


summary(lm(pca2$x[,1]~Diagnosis + Region+ totalAssignedGene + Sequencing.Batch + 
             Brain.Bank + RIN + Age + Sex + PC1+ PC2+ PC3+ PC4+ PC5,data = cbind(pd,qSVs)))

summary(lmer(pca2$x[,1]~Diagnosis+Region+ totalAssignedGene + Sequencing.Batch + 
               Brain.Bank + RIN + Age + Sex + PC1+ PC2+ PC3+ PC4+ PC5 + (1|Brain.ID),data = cbind(pd,qSVs)))

summary(lm(pca2$x[,1]~Diagnosis + Region+ totalAssignedGene + Sequencing.Batch + 
             Brain.Bank + RIN + Age + Sex,data = cbind(pd,qSVs)))

###############################################
# pca on the uniquely PTHS DEGs in human ASD and 15q
# sfig 6
whichGenes = rownames(geneMap)[match(outGene3$hsapien_homolog,geneMap$ensemblID)]
whichGenes = whichGenes[!is.na(whichGenes)]

pca3 = prcomp(t(log2(geneRpkm+1)[whichGenes,]))
# pca3resid = prcomp(t(yGene[whichGenes,]))
## clean the PCs, not expression
pca3resid = t(cleaningY(t(pca3$x), modAdj, P=3))

postscript('asd/plots/uniquely_pths_mouse_genes_in_asd_pca.eps',height=5,width=4)
boxplot(pca3resid$x[,1]~pd$Diagnosis,main = 'Only PTHS Gene Set PCA',xlab = '',outline = FALSE,
        ylab= paste0('Eigengene (',signif(summary(pca3resid)$importance[2,1]*100,2),'% of Variance) | Model'))
beeswarm(pca3resid$x[,1]~pd$Diagnosis,corral = 'wrap',add = TRUE,pch = 20,
	pwcol= col[as.numeric(factor(pd$Region))])
dev.off()

summary(lm(pca3$x[,1]~Diagnosis + Region+ totalAssignedGene + Sequencing.Batch + 
             Brain.Bank + RIN + Age + Sex + PC1+ PC2+ PC3+ PC4+ PC5,data = cbind(pd,qSVs)))

summary(lm(pca3$x[,1]~Diagnosis + Region+ totalAssignedGene + Sequencing.Batch + 
             Brain.Bank + RIN + Age + Sex,data = cbind(pd,qSVs)))




#4/10/2018
###############################
# the 9 genes in the GO as well
findMe = c('Plp1','Cldn11','Ugt8a','Mbp','Mag','Ermn','Cnp','Aspa','Gap43')
ind = which(outGene1$Symbol %in% findMe)
whichGenes = rownames(geneMap)[match(outGene1$hsapien_homolog[ind],geneMap$ensemblID)]
whichGenes = whichGenes[!is.na(whichGenes)]

pca4 = prcomp(t(log2(geneRpkm+1)[whichGenes,]))
# pca4resid = prcomp(t(yGene[whichGenes,]))
## clean the PCs, not expression
pca4resid = t(cleaningY(t(pca4$x), modAdj, P=3))

postscript('asd/plots/asd_9GO_genes_pca.eps',height=5,width=4)
boxplot(pca4resid[,1]~pd$Diagnosis,main = '9 GO CAG Gene Set PCA',xlab = '',outline=FALSE,
        ylab= paste0('Eigengene (',signif(summary(pca4resid)$importance[2,1]*100,2),'% of Variance) | Model'))
col = brewer.pal(3,'Set1')
beeswarm(pca4resid$x[,1]~pd$Diagnosis,corral = 'wrap',add = T,pch = 20, pwbg = col[as.numeric(factor(pd$Region))])
dev.off()
summary(lm(pca4$x[,1]~Diagnosis + Region+ totalAssignedGene + Sequencing.Batch + 
             Brain.Bank + RIN + Age + Sex + PC1+ PC2+ PC3+ PC4+ PC5,data = cbind(pd,qSVs)))



#####################
# heatmap of 34 genes
whichGenes = rownames(geneMap)[match(outGene1$hsapien_homolog,geneMap$ensemblID)]
whichGenes = whichGenes[!is.na(whichGenes)]
convergeMat = t(apply(yGene[whichGenes,],1,function(x) (x-mean(x))/sd(x)))
annoColor  = list(Diagnosis = c(CTL = '#f0f0f0',ASD= '#ef8a62',`Dup15` ='#67a9cf'),
                  Region = c(`ba9` = '#e41a1c',`ba41-42-22` = '#377eb8',`vermis` = '#4daf4a'))
breaksList = seq(-6, 6, by = .1)


pdf('asd/plots/asd_heatmap_covergent_genes.pdf',height = 6,width = 8)
ind = which(pd$Region=='ba9')
pheatmap(convergeMat[,ind],ylab = 'log 2 (RPKM+1) | qSV',legend = T,breaks = breaksList,
         #cluster_cols = F,
         color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(breaksList)), 
         border_color = NA, labels_col = rep('',ncol(yGene)), labels_row = geneMap[whichGenes,'Symbol'],
         annotation_col = pd[ind,c('Diagnosis','Region')],annotation_color = annoColor)

ind = which(pd$Region=='ba41-42-22')
pheatmap(convergeMat[,ind],ylab = 'log 2 (RPKM+1) | qSV',legend = T,breaks = breaksList,
         #cluster_cols = F,
         color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(breaksList)), 
         border_color = NA, labels_col = rep('',ncol(yGene)), labels_row = geneMap[whichGenes,'Symbol'],
         annotation_col = pd[ind,c('Diagnosis','Region')],annotation_color = annoColor)

ind = which(pd$Region=='vermis')
pheatmap(convergeMat[,ind],ylab = 'log 2 (RPKM+1) | qSV',legend = T,breaks = breaksList,
         #cluster_cols = F,
         color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(breaksList)), 
         border_color = NA, labels_col = rep('',ncol(yGene)), labels_row = geneMap[whichGenes,'Symbol'],
         annotation_col = pd[ind,c('Diagnosis','Region')],annotation_color = annoColor)
dev.off()

ind = which(pd$Region=='ba9')
#ind = which(pd$Region=='ba41-42-22')
cluster     = hclust(dist(t(yGene[whichGenes,ind])), method="ward.D2")
plot(cluster, label = pd$Diagnosis[ind])



