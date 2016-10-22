### Pten mouse differential expression analysis 
# qsub -V -l mf=10G,h_vmem=15G,h_stack=256M -cwd -b y R CMD BATCH analyze_pten.R
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

dir = "/dcl01/lieber/ajaffe/Brady/mouseRNAseq/TilotPTEN"

##############################################
# load phenotype data and RPKM expression data
load('./rdas/pheno.rda',envir = dat<-new.env())
load('/dcl01/lieber/ajaffe/Brady/mouseRNAseq/TilotPTEN/rpkmCounts_pten_OCT20_n18.rda')
pd = cbind(dat$pd,pd) #combine gene count data with experimental conditions
all.equal(pd$FileID,pd$SAMPLE_ID) #samples line up

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
hc = hclust(dist(t(yGene)))
plot(hc)

#############################
# Exploratory analysis
yGene = log2(geneRpkm+1)
pcaGene = prcomp(t(yGene))

plot(pcaGene$x, pch=21,bg=as.numeric(pd$Age),main = "Gene expression PCA", cex.lab = 1.25)
legend("right", legend = levels(pd$Age),pch =21, pt.bg = 1:3,cex= 1.25)

plot(pcaGene$x, pch=21,bg=as.numeric(pd$Genotype),main = "Gene expression PCA", cex.lab = 1.25)
legend("right", legend = levels(pd$Genotype),pch =21, pt.bg = 1:3,cex= 1.25)

plot(pcaGene$x[,3]~pd$totalAssignedGene,bg=as.numeric(pd$Age), pch=21,
     xlab = "Gene Assignment Rate", ylab = "PC2",
     main = "PC2 v. Gene Assignment Rate", cex.lab = 1.25)
legend("bottomright", legend = levels(pd$Age),pch =21, pt.bg = 1:3,cex= 1.25)

###################
# statistical model
# adjust for linear relationship w/ PC2 and gene assignment rate
mod = model.matrix(~Genotype+Age+totalAssignedGene, data=pd)

###################
# Gene DE, with SVA
yGene = log2(geneRpkm+1)
svaGene = sva(yGene, mod) 
modGene = cbind(mod, svaGene$sv)

fitGene = lmFit(yGene, modGene)
ebGene = ebayes(fitGene)
# take only homozygous DE, 3rd column of model matrix
outGene = data.frame(Symbol = geneMap$Symbol, log2FC = fitGene$coef[,3],  
                     pval = ebGene$p[,3], tval = ebGene$t[,3],
                     meanRpkm = rowMeans(geneRpkm),Location=geneMap$coord)
outGene = outGene[order(outGene$pval),]
outGene$fdr = p.adjust(outGene$pval, "fdr")
outGene$EntrezID = geneMap[rownames(outGene),c("EntrezID")]
sigGene = outGene[outGene$pval < 0.01,]
sum(sigGene$fdr<0.05) #219 genes FDR significant

#########################
# exon DE, with SVA
yExon = log2(exonRpkm+1)
svaExon = sva(yExon, mod)
modExon = cbind(mod, svaExon$sv)

fitExon = lmFit(yExon, modExon)
ebExon = ebayes(fitExon)
outExon = data.frame(Symbol = exonMap$Symbol,GeneID = exonMap$Geneid,
                     log2FC = fitExon$coef[,3], pval = ebExon$p[,3], tval = ebExon$t[,3],
                     Location = exonMap$coord, meanRpkm = rowMeans(exonRpkm))
outExon = outExon[order(outExon$pval),]
outExon$fdr = p.adjust(outExon$pval, "fdr")
outExon$EntrezID = exonMap[rownames(outExon),c("EntrezID")]
sigExon = outExon[outExon$pval < 0.01,]
sum(sigExon$fdr<0.05) #1 exon FDR significant

######################
### junction RP10Ms
txIDs = sapply(jMap$ensemblTx, paste, collapse=";")
yJxn = as.matrix(log2(jRpkm+1))
svaJxn = sva(yJxn, mod)
modJxn = cbind(mod, svaJxn$sv)

fitJxn = lmFit(yJxn, modJxn)
ebJxn = ebayes(fitJxn)
outJxn = data.frame(Symbol = jMap$newGeneSymbol,GeneID = jMap$newGeneID,
                    log2FC = fitJxn$coef[,3], pval = ebJxn$p[,3], tval = ebJxn$t[,3],
                    txCode = jMap$code,	meanRp10M = rowMeans(jRpkm),
                    Location=ss(rownames(jRpkm),"\\("),TxIDs =txIDs)
outJxn = outJxn[order(outJxn$pval),]
outJxn$fdr = p.adjust(outJxn$pval, "fdr")
outJxn$EntrezID = geneMap[as.character(outJxn$GeneID),c("EntrezID")]
sigJxn = outJxn[outJxn$pval < 0.01,]
sum(sigJxn$fdr<.05) #0 DE junctions

#####################
# save everything! :)
WriteXLS(list(Gene = sigGene,Exon= sigExon, Junction = sigJxn),
         ExcelFileName = 'tables/pten_DE_table.xls')
save(svaGene, svaExon, svaJxn, file="rdas/pten_sva_objects.rda")
save(pd,outGene, outExon, outJxn, file="rdas/pten_DE_objects.rda")

