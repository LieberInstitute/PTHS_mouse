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

##############################
# create and run DESeq objects
geneDdsListAdj <-lapply(indList,function(i){ 
  DESeq2(countData = geneCounts[,i], colData = pd[i,], 
         design = ~Diagnosis + totalAssignedGene +
           Sequencing.Batch + Brain.Bank + RIN + Age + Sex,
         sva = FALSE,parallel=TRUE)})

################################################################
# get DE results, and fold-change homozygous mutant v. wild-type
resGene <- lapply(geneDdsListAdj,function(g) results(g,contrast = c('Diagnosis','CTL','ASD'), alpha=0.05))

sapply(resGene,function(g) sum(g$padj < 0.05, na.rm=TRUE))

outGeneList <- lapply(resGene,function(g) as.data.frame(g))
labs = Reduce(intersect,lapply(outGeneList,rownames))
outGene = do.call('cbind',lapply(outGeneList, function(g) g[labs,]))
outGene = cbind(outGene,geneMap[rownames(outGene),])
sigGene = outGene[apply(outGene[,grep('pvalue',names(outGene))],1,function(x) any(x<0.01)),]
sigGene = outGene[order(apply(sigGene[,grep('pvalue',names(sigGene))],1,mean)),]


pdf('plots/DESeq2_MA_plots_asd.pdf')
plotMA(resGene[[1]], main="ba41-42-22 Gene MA plot", ylim=c(-2,2))
plotMA(resGene[[2]], main="ba9 Gene MA plot", ylim=c(-2,2))
plotMA(resGene[[3]], main="vermis Gene MA plot", ylim=c(-2,2))
dev.off()

#################
# save everything
library(WriteXLS)
WriteXLS(list(Gene = sigGene,pd), ExcelFileName = 'tables/asd_DE_table_DESeq2.xls',row.names=T)
save(outGene, file = 'rdas/asd_DE_objects_DESeq2.rda')
save(geneDdsListAdj, file = '/dcl01/lieber/ajaffe/Brady/asd/asd_DESeq2.rda')
cat(rownames(pd),file = '/dcl01/lieber/ajaffe/Brady/asd/samples.txt',sep = '\n')
