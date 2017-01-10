### asd mouse differential expression analysis 
# qsub -V -l mf=10G,h_vmem=15G,h_stack=256M -cwd -b y R CMD BATCH analyze_asd_DESeq2.R
source('../DESeq2_functions.R') #work-horse of differential expression

##############################################
# load phenotype data and RPKM expression data
load('./rdas/pheno.rda',envir = dat<-new.env())
load('/dcl01/lieber/ajaffe/Brady/asd/rawCounts_asd_dec9_n187.rda')
pd = cbind(pd,dat$pd)
all.equal(pd$SampleID,pd$SAMPLE_ID) #samples line up
rownames(pd) = pd$SAMPLE_ID

########################
# split by brain regions
indList = split(seq(nrow(pd)), pd$Region)

pd$Diagnosis = factor(pd$Diagnosis, levels = c("CTL", "ASD"))

#########################
# quality control metrics
library(beeswarm)
par(mar = c(2,5,3,2))
pdf('plots/asd_qc_plots.pdf')
boxplot(pd$totalAssignedGene~pd$Diagnosis+pd$Region,xaxt ='n',
        main = 'Gene Assignment Rate',ylab= 'Gene Assignment Rate')
axis(side = 1, at = 1:3*2-.5,labels = levels(factor(pd$Region)))
beeswarm(pd$totalAssignedGene~pd$Diagnosis+pd$Region,
         add = T,pch = 21,pwbg = pd$Diagnosis)
legend('topright',legend = c('CTL','ASD'),pch =21, pt.bg = 1:2)

boxplot(pd$mitoRate~pd$Diagnosis+pd$Region,xaxt ='n',
        main = 'Mito Rate', ylab = 'Mito Rate')
axis(side = 1, at = 1:3*2-.5,labels = levels(factor(pd$Region)))
beeswarm(pd$mitoRate~pd$Diagnosis+pd$Region,
         add = T,pch = 21,pwbg = pd$Diagnosis)
legend('topright',legend = c('CTL','ASD'),pch =21, pt.bg = 1:2)
dev.off()

##############################
# create and run DESeq objects
geneDdsListAdj <-lapply(indList,function(i){ 
  DESeq2(countData = geneCounts[,i], colData = pd[i,], 
         design = ~Diagnosis + totalAssignedGene +
           Sequencing.Batch + Brain.Bank + RIN + Age + Sex,
         sva = TRUE,parallel=TRUE)})

################################################################
# get DE results, and fold-change homozygous mutant v. wild-type
resGene <- lapply(geneDdsListAdj,function(g) results(g,contrast = c('Diagnosis','CTL','ASD'), alpha=0.05))

sapply(resGene,function(g) sum(g$padj < 0.05, na.rm=TRUE))

outGeneList <- lapply(resGene,function(g) {
  outGene <- as.data.frame(g)
  outGene = outGene[order(outGene$padj,outGene$pvalue),]
  outGene = cbind(outGene,geneMap[rownames(outGene),])
})

sigGeneList = lapply(outGeneList,function(g){
  sigGene = g[which(g$pvalue<.01),]
})

pdf('plots/DESeq2_MA_plots_asd.pdf')
plotMA(resGene[[1]], main="ba41-42-22 Gene MA plot", ylim=c(-2,2))
plotMA(resGene[[2]], main="ba9 Gene MA plot", ylim=c(-2,2))
plotMA(resGene[[3]], main="vermis Gene MA plot", ylim=c(-2,2))
dev.off()

#################
# save everything
library(WriteXLS)
WriteXLS(sigGeneList, ExcelFileName = 'tables/asd_DE_table_DESeq2.xls',row.names=T)
save(outGeneList, file = 'rdas/asd_DE_objects_DESeq2.rda')
save(geneDdsListAdj, file = '/dcl01/lieber/ajaffe/Brady/asd/asd_DESeq2_Adj.rda')
#cat(rownames(pd),file = '/dcl01/lieber/ajaffe/Brady/asd/samples.txt',sep = '\n')
