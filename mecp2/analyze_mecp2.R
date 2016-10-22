### Maher mouse differential expression analysis 
# qsub -V -l mf=10G,h_vmem=15G,h_stack=256M -cwd -b y R CMD BATCH analyze_mecp2.R
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

dir = "/dcl01/lieber/ajaffe/Brady/mouseRNAseq/Rett"

##############################################
# load phenotype data and RPKM expression data
load('./rdas/pheno.rda',envir = dat<-new.env())
load('/dcl01/lieber/ajaffe/Brady/mouseRNAseq/Rett/rpkmCounts_mecp2_OCT20_n6.rda')
pd = cbind(dat$pd,pd)
all.equal(pd$SampleID,pd$SAMPLE_ID) #samples line up

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

plot(pcaGene$x, pch=21,bg=as.numeric(pd$Genotype),main = "Gene expression PCA", cex.lab = 1.25)
legend("right", legend = levels(pd$Genotype),pch =21, pt.bg = 1:3,cex= 1.25)
#separate by genotype

###################
# statistical model
mod = model.matrix(~Genotype, data=pd)

#############################################
# Gene DE, no need to use surrogate variables
yGene = log2(geneRpkm+1)
fitGene = lmFit(yGene, mod)
ebGene = ebayes(fitGene)
# take only homozygous DE, 3rd column of model matrix
outGene = data.frame(Symbol = geneMap$Symbol, log2FC = fitGene$coef[,2],
                     pval = ebGene$p[,2], tval = ebGene$t[,2],
                     meanRpkm = rowMeans(geneRpkm),Location=geneMap$coord)
outGene = outGene[order(outGene$pval),]
outGene$fdr = p.adjust(outGene$pval, "fdr")
outGene$EntrezID = geneMap[rownames(outGene),c("EntrezID")]
sigGene = outGene[outGene$pval < 0.01,]
sum(sigGene$fdr<0.05) #only MeCP2 is FDR significant

#########################
# exon DE, with SVA
yExon = log2(exonRpkm+1)

fitExon = lmFit(yExon, mod)
ebExon = ebayes(fitExon)
outExon = data.frame(Symbol = exonMap$Symbol,GeneID = exonMap$Geneid,
                     log2FC = fitExon$coef[,2], pval = ebExon$p[,2], tval = ebExon$t[,2],
                     Location = exonMap$coord, meanRpkm = rowMeans(exonRpkm))
outExon = outExon[order(outExon$pval),]
outExon$fdr = p.adjust(outExon$pval, "fdr")
outExon$EntrezID = exonMap[rownames(outExon),c("EntrezID")]
sigExon = outExon[outExon$pval < 0.01,]
sum(sigExon$fdr<0.05) #14 exons FDR significant

######################
### junction RP10Ms
txIDs = sapply(jMap$ensemblTx, paste, collapse=";")
yJxn = as.matrix(log2(jRpkm+1))

fitJxn = lmFit(yJxn, mod)
ebJxn = ebayes(fitJxn)
outJxn = data.frame(Symbol = jMap$newGeneSymbol,GeneID = jMap$newGeneID,
                    log2FC = fitJxn$coef[,2], pval = ebJxn$p[,2], tval = ebJxn$t[,2],
                    txCode = jMap$code,	meanRp10M = rowMeans(jRpkm),
                    Location=ss(rownames(jRpkm),"\\("),TxIDs =txIDs)
outJxn = outJxn[order(outJxn$pval),]
outJxn$fdr = p.adjust(outJxn$pval, "fdr")
outJxn$EntrezID = geneMap[as.character(outJxn$GeneID),c("EntrezID")]
sigJxn = outJxn[outJxn$pval < 0.01,]
sum(sigJxn$fdr<.05) #3 DE junctions

#####################
# save everything! :)
WriteXLS(list(Gene = sigGene,Exon= sigExon, Junction = sigJxn),
         ExcelFileName = 'tables/mecp2_DE_table.xls')
save(pd, outGene, outExon, outJxn, file="rdas/mecp2_DE_objects.rda")

