# load in CIBERSORT signature gene table


options(stringsAsFactors = F)
dat = read.delim('tables/zhang_CIBERSORT_signature.txt',row.names = 1)
dat = dat[,apply(dat,2,function(x) all(!is.na(x)))]

#############################################
# choose cell-type by column of max expression
cellLabs = data.frame(genes = rownames(dat),
                      labs = names(dat)[apply(dat,1,which.max)],
                      var = apply(dat,1,function(x) x[which.max(x)])/apply(dat,1,sum))
cellLabs$labs = factor(cellLabs$labs, levels = c('Astrocyte','Endothelial','Microglia','Neuron',
                                                 'Pre.Oligo','New.Oligo','Mye.Oligo'),
                       labels = c('Astrocyte','Endothelial','Microglia','Neuron',
                                  'Oligo_Progenitor','Oligo_New','Oligo_Mye'))

library(beeswarm)
pdf('plots/cell_types.pdf')
boxplot(var~labs,data = cellLabs,xlab = 'Cell Type',
        ylab = "Cell type % of expression")
beeswarm(var~labs,data = cellLabs,add = T,corral = 'wrap')

cellSet = split(cellLabs$genes,cellLabs$labs)
save(cellSet,file = 'rdas/mouse_brain_celltypes_geneset.rda')


################################
#get Ensembl mouse to human genes
library(biomaRt)
ensembl=useMart("ensembl",dataset = 'mmusculus_gene_ensembl')
MMtoHG = getBM(attributes = c('ensembl_gene_id','hsapiens_homolog_ensembl_gene'),mart = ensembl)
hgMart = getBM(attributes = c('ensembl_gene_id','illumina_humanht_12_v3','entrezgene'),
               mart = useMart("ensembl",dataset = 'hsapiens_gene_ensembl'))

##########################################
# make human homolog zhang expression table
tmp = MMtoHG$hsapiens_homolog_ensembl_gene[match(rownames(dat),MMtoHG$ensembl_gene_id)]
ind = which(tmp!=''& !duplicated(tmp))
tmp2 = dat[ind,]
rownames(tmp2) = tmp[ind]
write.table(tmp2,sep = '\t',quote = FALSE,file = 'tables/signature_zhang_hg_homolog.txt')

