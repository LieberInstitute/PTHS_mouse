### philpot mice differential expression analysis 
# qsub -V -l mf=10G,h_vmem=15G,h_stack=256M -cwd -b y R CMD BATCH analyze_philpot.R
cleaningY = function(y, mod, P=ncol(mod)) {
  Hat=solve(t(mod)%*%mod)%*%t(mod)
  beta=(Hat%*%t(y))
  cleany=y-t(as.matrix(mod[,-c(1:P)])%*%beta[-c(1:P),])
  return(cleany)
}

library(limma)
library(sva)
library(ggplot2)
library(jaffelab)
library(WriteXLS)

dir = "/dcl01/lieber/ajaffe/Brady/philpot"
##############################################
# load phenotype data and RPKM expression data
load('./rdas/pheno.rda',envir = dat<-new.env())
load('/dcl01/lieber/ajaffe/Brady/philpot/rpkmCounts_philpot_OCT20_n58.rda')
pd = cbind(dat$pd,pd)
all(seq(nrow(pd))==unlist(sapply(pd$SAMPLE_ID,grep,pd$'BGI file name'))) #samples line up
rownames(pd) = pd$'File rename/Unique identifier'

geneMap$coord = paste0(geneMap$Chr, ":", geneMap$Start , "-", geneMap$End)
exonMap$coord = paste0(exonMap$Chr, ":", exonMap$Start, "-", exonMap$End)

## expressed genes and exons
meanGeneExprs = rowMeans(geneRpkm)
meanExonExprs = rowMeans(exonRpkm)
meanJxnExprs = rowMeans(jRpkm)

rpkmCuts = c(0,0.1,0.2,0.5,1,2,5,10,50)
names(rpkmCuts) = paste0("RPKM > ", rpkmCuts)
numGeneExprs = sapply(rpkmCuts, function(x) sum(meanGeneExprs > x))
numGeneExprs
numExonExprs = sapply(rpkmCuts, function(x) sum(meanExonExprs > x))
numGeneByExonExprs = sapply(rpkmCuts, function(x)
  length(unique(exonMap$Geneid[meanExonExprs > x])))
numJxnExprs = sapply(rpkmCuts, function(x) sum(meanJxnExprs > x))
numJxnExprs

table(jMap$code[meanJxnExprs > 0.2])

##############################
# filter expression levels ###
# gene
gIndex = which(rowMeans(geneRpkm) > 0.1)
geneRpkm = geneRpkm[gIndex,]
geneMap = geneMap[gIndex,]

# exon
eIndex = which(rowMeans(exonRpkm) > 0.1)
exonRpkm = exonRpkm[eIndex,]
exonMap = exonMap[eIndex,]

# junction
jIndex=which(rowMeans(jRpkm) > 1 & jMap$code != "Novel")
jRpkm = jRpkm[jIndex,]
jMap = jMap[jIndex]

dim(geneRpkm) ; dim(exonRpkm) ; dim(jRpkm)

######################
# hierarchal clustering
yGene = log2(geneRpkm+1)
colnames(yGene) = rownames(pd)
hc = hclust(dist(t(yGene)))
plot(hc)

#############################
# Exploratory analysis
yGene = log2(geneRpkm+1)
pcaGene = prcomp(t(yGene))

plot(pcaGene$x, pch=21,bg=as.numeric(pd$Age),main = "Gene expression PCA", cex.lab = 1.25)
legend("bottomright", legend = levels(pd$Age),pch =21, pt.bg = 1:3,cex= 1.25)

plot(pcaGene$x, pch=21,bg=as.numeric(pd$Genotype),main = "Gene expression PCA", cex.lab = 1.25)
legend("bottomright", legend = levels(pd$Genotype),pch =21, pt.bg = 1:3,cex= 1.25)

#mitochondia map rate correlate with PC1
#repurification date clusters mito-map rate
plot(pd$mitoRate, pcaGene$x[,1],pch=21,bg=as.numeric(pd$Repurify),
     ylab = "PC1", xlab = "Mitochondria Map Rate")
legend("topright", legend = c("Not repufied on 2016-03-22","Repurified on 2016-03-22"), 
       pch =21, pt.bg = c(0:1))

plot(pd$mitoRate, pcaGene$x[,1],pch=21,bg=as.numeric(pd$Age),
     ylab = "PC1", xlab = "Mitochondria Map Rate")
legend("topright", legend = levels(pd$Age), pch =21, pt.bg = c(1:4))

###################
# statistical model
# adjust for mitochondria map rate
mod = model.matrix(~Genotype+Age+mitoRate+Repurify, data=pd)

###################
# Gene DE, with SVA
yGene = log2(geneRpkm+1)
svaGene = sva(yGene, mod) 
modGene = cbind(mod, svaGene$sv)

fitGene = lmFit(yGene, modGene)
ebGene = ebayes(fitGene)
outGene = data.frame(Symbol = geneMap$Symbol,meanRpkm = rowMeans(geneRpkm),Coordinates=geneMap$coord)
for (i in c(2:5)) {
  tmp = data.frame(log2FC = fitGene$coef[,i], pval = ebGene$p[,i], tval = ebGene$t[,i])
  tmp$fdr = p.adjust(tmp$pval, "fdr")
  print(levels(pd$Geno)[i])
  print(sum(tmp$fdr < 0.05))
  names(tmp) = paste0(names(tmp),"_",levels(pd$Geno)[i])
  outGene = cbind(outGene,tmp)
}
outGene = outGene[order(
  apply(outGene[,paste0("pval_",c(levels(pd$Geno)[2:5])),],1,
        function(x) mean(log10(x)))),]
sigGene = outGene[apply(
  outGene[,paste0("pval_",c(levels(pd$Geno)[2:5])),],1,
  function(x) any(x < 0.01)),]

#####################
# exon DE, with SVA
yExon = log2(exonRpkm+1)
svaExon = sva(yExon, mod) # 5 SV's
modExon = cbind(mod, svaExon$sv)

fitExon = lmFit(yExon, modExon)
ebExon = ebayes(fitExon)
outExon = data.frame(Symbol = exonMap$Symbol,GeneID = exonMap$Geneid,
                     Location = exonMap$coord, meanRpkm = rowMeans(exonRpkm))
for (i in c(2:5)) {
  tmp = data.frame(log2FC = fitExon$coef[,i], pval = ebExon$p[,i], tval = ebExon$t[,i])
  tmp$fdr = p.adjust(tmp$pval, "fdr")
  print(levels(pd$Geno)[i])
  print(sum(tmp$fdr < 0.05))
  names(tmp) = paste0(names(tmp),"_",levels(pd$Geno)[i])
  outExon = cbind(outExon,tmp)
}
outExon = outExon[order(
  apply(outExon[,paste0("pval_",c(levels(pd$Geno)[2:5])),],1,
        function(x) mean(log10(x)))),]
sigExon = outExon[apply(
  outExon[,paste0("pval_",c(levels(pd$Geno)[2:5])),],1,
  function(x) any(x < 0.01)),]

######################
### junction RP10Ms
txIDs = sapply(jMap$ensemblTx, paste, collapse=";")
yJxn = as.matrix(log2(jRpkm+1))
svaJxn = sva(yJxn, mod) 
modJxn = cbind(mod, svaJxn$sv)

fitJxn = lmFit(yJxn, modJxn)
ebJxn = ebayes(fitJxn)
outJxn = data.frame(Symbol = jMap$newGeneSymbol,GeneID = jMap$newGeneID,txCode = jMap$code,	
                    meanRp10M = rowMeans(jRpkm), Coordinates=ss(rownames(jRpkm),"\\("),
                    TxIDs = txIDs)
for (i in c(2:5)) {
  tmp = data.frame(log2FC = fitJxn$coef[,i], pval = ebJxn$p[,i], tval = ebJxn$t[,i])
  tmp$fdr = p.adjust(tmp$pval, "fdr")
  print(levels(pd$Geno)[i])
  print(sum(tmp$fdr < 0.05))
  names(tmp) = paste0(names(tmp),"_",levels(pd$Geno)[i])
  outJxn = cbind(outJxn,tmp)
}
outJxn = outJxn[order(
  apply(outJxn[,paste0("pval_",c(levels(pd$Geno)[2:5])),],1,
        function(x) mean(log10(x)))),]
sigJxn = outJxn[apply(
  outJxn[,paste0("pval_",c(levels(pd$Geno)[2:5])),],1,
  function(x) any(x < 0.01)),]

#######################################
# save all the differential expressions
WriteXLS(list(Gene = sigGene,Exon= sigExon, Junction = sigJxn, phenotype = pd),
         ExcelFileName = 'tables/philpot_DE_table.xls')
save(svaGene, svaExon, svaJxn, file="rdas/philpot_sva_objects.rda")
save(pd, outGene, outExon, outJxn, file="rdas/philpot_DE_objects.rda")


