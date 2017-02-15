### PTEN mouse differential expression analysis using DESeq2
# qsub -V -l mf=100G,h_vmem=150G,h_stack=256M -cwd -b y R CMD BATCH analyze_pten_DESeq2.R
source('../DESeq2_functions.R') #work-horse of differential expression

##############################################
# load phenotype data and RPKM expression data
load('./rdas/pheno.rda',envir = dat<-new.env())
load('/dcl01/lieber/ajaffe/Brady/mouseRNAseq/TilotPTEN/rawCounts_pten_OCT20_n18.rda')
pd = cbind(dat$pd,pd)
all.equal(pd$FileID,pd$SAMPLE_ID) #samples line up
rownames(pd) = pd$SAMPLE_ID
pd$Age = factor(pd$Age,levels = levels(pd$Age),labels = c('P14','Adult'))

#####################
# split by age groups
indList = split(seq(nrow(pd)), pd$Age)
names(indList)

##############################
# create and run DESeq objects
geneDds <- lapply(indList,function(i) {
  DESeq2(countData = geneCounts[,i], colData = pd[i,], 
         design = ~Genotype,sva = TRUE,parallel=TRUE)})

############################################
# get DE results, and fold-change PTHS v. WT
resGene <- c(lapply(geneDds,results,contrast = c('Genotype','Pten wt/m3m4','Pten wt/wt'), alpha=0.05),
             lapply(geneDds,results,contrast = c('Genotype','Pten m3m4/m3m4','Pten wt/wt'), alpha=0.05))

names(resGene) = paste0(names(resGene),'.',rep(c('Het','Mut'),each =2))
sapply(resGene,function(g) sum(g$padj < 0.05, na.rm=TRUE))
# P14.Het Adult.Het   P14.Mut Adult.Mut  
#      13       111       205      2818

outGeneList <- lapply(resGene,function(g) {
  outGene <- as.data.frame(g)
  outGene = outGene[order(outGene$padj,outGene$pvalue),]
  outGene = cbind(outGene,geneMap[rownames(outGene),])
})

sigGeneList = lapply(outGeneList,function(g){
  sigGene = g[which(g$pvalue<.01),]
})

pdf('plots/DESeq2_MA_plots_pten_ages.pdf')
plotMA(resGene[['P14.Het']], main="P14 Het MA plot", ylim=c(-1,1))
plotMA(resGene[['Adult.Het']], main="Adult Het MA plot", ylim=c(-2,2))
plotMA(resGene[['P14.Mut']], main="P14 Mut MA plot", ylim=c(-1,1))
plotMA(resGene[['Adult.Mut']], main="Adult Mut MA plot", ylim=c(-2,2))
dev.off()

#################
# save everything
library(WriteXLS)
WriteXLS(sigGeneList, ExcelFileName = 'tables/pten_DE_ages_table_DESeq2.xls',row.names=T)
save(sigGeneList,outGeneList,file = 'rdas/pten_ages_DE_objects_DESeq2.rda')
save(geneDds,resGene, file = '/dcl01/lieber/ajaffe/Brady/mouseRNAseq/TilotPTEN/pten_ages_DESeq2_svaAdj.rda')

################
# make some plots
genes = unique(unlist(sapply(outGeneList,function(g) rownames(g)[g$padj < 0.05 & !is.na(g$padj)])))

dat = do.call('cbind',lapply(outGeneList[c(1,3,2,4)],function(g) g[genes,c('stat')]))
dat = dat[complete.cases(dat),]
pairs(dat)


