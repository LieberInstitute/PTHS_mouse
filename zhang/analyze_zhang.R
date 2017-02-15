# look at TCF4 expression relative to cell types
source('../DESeq2_functions.R')
library(DESeq2)
library(jaffelab)

###############
# load the data
load('./rdas/pheno.rda',envir = dat <- new.env())
load('/dcl01/lieber/ajaffe/Brady/zhang/rawCounts_zhang_nov1_n17.rda')
pd = cbind(pd,dat$pd)
table(pd$SAMPLE_ID==pd$SampleID)

##################
# relabel cell types
pd$Cell = as.character(pd$Cell)
pd$Cell = factor(pd$Cell,levels = c('<not provided>','Astrocyte','endothelial cells', 'microglia',
                                    'neuron','oligodendrocyte precursor cells',
                                    'newly formed oligodendrocytes','myelinating oligodendrocytes'),
                 labels = c('Cortex','Astrocyte','Endothelial','Microglia','Neuron','Oligo_Pre',
                            'Oligo_New','Oligo_Mye'))

########################
# take only TRAP samples
ind = which(pd$Cell != 'Cortex')
pd = pd[ind,]; pd$Cell = droplevels(pd$Cell)
geneCounts = geneCounts[,ind]

##########################################
# split DE of each cell type v. all others
cells = levels(pd$Cell); names(cells) <- levels(pd$Cell)
geneDds = lapply(cells,function(cell){
  tmp = pd; tmp$isCell = factor(ifelse(pd$Cell==cell,'Yes','No'),levels= c('No','Yes'))
  return(DESeq2(geneCounts,tmp,~isCell, sva = TRUE,parallel=TRUE))
})

###################################################
# find genes with differentially expressed
resGene <- lapply(geneDds,results,contrast = c('isCell','Yes','No'), alpha=0.05) 
sapply(resGene,function(g) sum(g$padj < 0.05, na.rm=TRUE))

outGeneList <- lapply(resGene,function(g) {
  outGene <- as.data.frame(g)
  outGene = outGene[order(outGene$padj,outGene$pvalue),]
  outGene = cbind(outGene,geneMap[rownames(outGene),])
})

sigGeneList = lapply(outGeneList,function(g){
  sigGene = with(g, g[which(padj<.05 & !is.na(padj)),])})
sapply(sigGeneList,nrow) # upregulated in these cell types

#tmpGenes = unlist(lapply(sigGeneList,rownames))
#tmpGenes = unique(tmpGenes[duplicated(tmpGenes)]) # no overlap
#sigGeneList = lapply(sigGeneList,function(g){
#  g[which(!(rownames(g) %in% tmpGenes)),]})
#sapply(sigGeneList,nrow) #uniquely upregulate in only these cell types

pdf('plots/DESeq2_MA_plots_zhang_celltypes.pdf')
for(n in names(resGene)) plotMA(resGene[[n]], main=paste(n,"MA plot"), ylim=c(-10,10))
dev.off()

#################
# save everything
library(WriteXLS)
WriteXLS(sigGeneList, ExcelFileName = 'tables/zhang_celltypes_DE_table_DESeq2.xls',row.names=T)
save(sigGeneList,outGeneList,file = 'rdas/zhang_celltypes_DE_objects_DESeq2.rda')
save(geneDds,resGene, file = '/dcl01/lieber/ajaffe/Brady/zhang/zhang_celltypes_DESeq2_svaAdj.rda')






#################################
# check CIBERSORT signature genes
options(stringsAsFactors = F)
dat = read.delim('tables/zhang_CIBERSORT_signature.txt',row.names = 1)
dat = dat[,apply(dat,2,function(x) all(!is.na(x)))]

#############################################
# choose cell-type by column of max expression
cellSet = split(rownames(dat),names(dat)[apply(dat,1,which.max)])
names(cellSet) = c("Astrocyte","Endothelial","Microglia", "Oligo_Mye",
                   "Neuron","Oligo_New", "Oligo_Pre" )

tmp = names(cellSet);names(tmp) = names(cellSet)
lapply(tmp,function(cell) 
  table(cellSet[[cell]] %in% rownames(sigGeneList[[cell]])))
cellSet = lapply(tmp,function(cell) 
  cellSet[[cell]][cellSet[[cell]] %in% rownames(sigGeneList[[cell]])])
save(cellSet,file = 'rdas/mouse_brain_celltypes_geneset_DESeq2.rda')






save(cellSet,file = 'rdas/mouse_brain_celltypes_geneset_DESeq2.rda')

################
# check overlap
