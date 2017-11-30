# look at TCF4 expression relative to cell types
library(DESeq2)
library(beeswarm)
library(jaffelab)
library(ggplot2)

###############
# load the data
load('./rdas/pheno.rda',envir = dat <- new.env())
load('/dcl01/lieber/ajaffe/Brady/zhang/rpkmCounts_zhang_nov1_n17.rda')
pd = cbind(pd,dat$pd)
table(pd$SAMPLE_ID==pd$SampleID)

##################
# relabel cell types
pd$Cell = as.character(pd$Cell)
pd$Cell = factor(pd$Cell,levels = c('<not provided>','endothelial cells','neuron','Astrocyte',
                                    'oligodendrocyte precursor cells',
                                    'newly formed oligodendrocytes','myelinating oligodendrocytes',
                                    'microglia'),
                 labels = c('Cortex','Endothelial','Neuron','Astrocyte','Oligo.Pre',
                            'Oligo.New','Oligo.Mye','Microglia'))

##########################
# plot TCF4 expression in these cells
ind = which(geneMap$Symbol =='Tcf4')

postscript('plots/cell-specific_tcf4_expression.eps',height = 2.5,width = 2.5)
ggplot(data =cbind(rpkm = geneRpkm[ind,],pd),aes(x = Cell, y = rpkm,fill = Cell))+
  geom_point(pch =21) +guides(FALSE) +theme(axis.text.x = element_text(angle = 30, hjust = 1))+
  scale_fill_brewer(palette = 'Set1')+ylab('FPKM')+ylim(c(0,max(geneRpkm[ind,])))+
  guides(fill = FALSE)
dev.off()


############################################
# format data for CIBERSORT signature matrix
ind = 
write.table(geneCounts[,ind],sep = '\t',quote = FALSE,file = 'tables/rawGeneCounts_zhang.txt')
mod = 2-t(model.matrix(~pd$Cell)[ind,-1])
colnames(mod) = pd$SampleID[ind]
rownames(mod) = ss(rownames(mod),'pd\\$Cell',2)
write.table(mod,sep = '\t',quote = FALSE,file = 'tables/sample_reference_zhang.txt')
save(geneDds,file = '/dcl01/lieber/ajaffe/Brady/zhang/DESeq_zhang.rda')


