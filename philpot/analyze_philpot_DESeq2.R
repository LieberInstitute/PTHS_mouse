### philpot mice differential expression analysis 
# qsub -V -l mf=200G,h_vmem=220G,h_stack=256M -cwd -b y R CMD BATCH analyze_philpot_DESeq2.R
source('../DESeq2_functions.R') #work-horse of differential expression

##############################################
# load phenotype data and RPKM expression data
load('./rdas/pheno.rda',envir = dat<-new.env())
load('/dcl01/lieber/ajaffe/Brady/philpot/rawCounts_philpot_OCT20_n58.rda')
pd = cbind(dat$pd,pd)
all(seq(nrow(pd))==unlist(sapply(pd$SAMPLE_ID,grep,pd$'BGI file name'))) #samples line up
rownames(pd) = pd$'File rename/Unique identifier'

##########################
# conform to DESeq2 data types
jCounts = as.matrix(as.data.frame(jCounts))
jIndex=which(jMap$code != "Novel")
jCounts = jCounts[jIndex,]
jMap = jMap[jIndex]
colnames(jCounts) = pd$SAMPLE_ID
rownames(pd) = pd$SAMPLE_ID

##############################
# create and run DESeq objects
geneDds <- DESeq2(countData = geneCounts, colData = pd, design = ~Genotype+Age,sva = TRUE,parallel=TRUE)
rm(geneCounts); gc()
exonDds <- DESeq2(countData = exonCounts, colData = pd, design = ~Genotype+Age,sva = TRUE,parallel=TRUE)
rm(exonCounts); gc()
jxnDds <- DESeq2(countData = jCounts, colData = pd, design = ~Genotype+Age,sva = TRUE,parallel=TRUE)
rm(jCounts); gc()

############################################################
# get DE results, and fold-change of each het mouse genotype
outGene = data.frame(row.names = rownames(geneDds),geneMap[rownames(geneDds),])
outExon = data.frame(row.names = rownames(exonDds),exonMap[rownames(exonDds),])
outJxn = data.frame(row.names = rownames(jxnDds),as.data.frame(jMap)[rownames(jxnDds),])

for (i in seq(2,5)) {
  het = levels(pd$Genotype)[i]
  resGene <- results(geneDds,contrast = c('Genotype',het,'WT'),alpha=0.05) 
  resExon <- results(exonDds,contrast =  c('Genotype',het,'WT'),alpha=0.05) 
  resJxn <- results(jxnDds,contrast =  c('Genotype',het,'WT'),alpha=0.05) 
  
  g <- as.data.frame(resGene); names(g) = paste0(names(g),'_',het)
  e <- as.data.frame(resExon); names(e) = paste0(names(e),'_',het)
  j <- as.data.frame(resJxn); names(j) = paste0(names(j),'_',het)
  
  outGene = cbind(outGene,g)
  outExon = cbind(outExon,e)
  outJxn = cbind(outJxn,j)
}
sapply(outGene[,paste0("padj_",c(levels(pd$Geno)[2:5])),],function(x)sum(x<.05,na.rm = T))
sapply(outExon[,paste0("padj_",c(levels(pd$Geno)[2:5])),],function(x)sum(x<.05,na.rm = T))
sapply(outJxn[,paste0("padj_",c(levels(pd$Geno)[2:5])),],function(x)sum(x<.05,na.rm = T))

outGene = outGene[order(apply(outGene[,paste0("pvalue_",c(levels(pd$Geno)[2:5])),],1,function(x) mean(log10(x)))),]
sigGene = outGene[apply(outGene[,paste0("pvalue_",c(levels(pd$Geno)[2:5])),],1,function(x) any(x < 0.01)),]

outExon = outExon[order(apply(outExon[,paste0("pvalue_",c(levels(pd$Geno)[2:5])),],1,function(x) mean(log10(x)))),]
sigExon = outExon[apply(outExon[,paste0("pvalue_",c(levels(pd$Geno)[2:5])),],1,function(x) any(x < 0.01)),]

outJxn = outJxn[order(apply(outJxn[,paste0("pvalue_",c(levels(pd$Geno)[2:5])),],1,function(x) mean(log10(x)))),]
sigJxn = outJxn[apply(outJxn[,paste0("pvalue_",c(levels(pd$Geno)[2:5])),],1,function(x) any(x < 0.01)),]

#######################################
# save all the differential expressions
library(WriteXLS)
WriteXLS(list(Gene = sigGene,Exon= sigExon, Junction = sigJxn, phenotype = pd),
         ExcelFileName = 'tables/philpot_DE_table.xls')
save(svaGene, svaExon, svaJxn, file="rdas/philpot_sva_objects.rda")
save(pd, outGene, outExon, outJxn, file="rdas/philpot_DE_objects_DESeq2.rda")
save(geneDds,exonDds,jxnDds, file = '/dcl01/lieber/ajaffe/Brady/philpot/philpot_DESeq2_svaAdj.rda')


