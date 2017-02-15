# look at TCF4 expression relative to cell types
library(DESeq2)
library(beeswarm)
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
pd$Cell = factor(pd$Cell,levels = c('<not provided>','endothelial cells','neuron','Astrocyte',
                                    'oligodendrocyte precursor cells',
                                    'newly formed oligodendrocytes','myelinating oligodendrocytes',
                                    'microglia'),
                 labels = c('Cortex','Endothelial','Neuron','Astrocyte','Oligo_Pre',
                            'Oligo_New','Oligo_Mye','Microglia'))

##########################
# plot TCF4 expression in these cells
ind = which(geneMap$Symbol =='Tcf4')

pdf('plots/cell-specific_tcf4_expression.pdf',height = 4,width = 3)
par(mar = c(6,2,2,1),cex.lab = 1 ,font = 2)
beeswarm(geneRpkm[ind,]~pd$Cell, ylab = 'FPKM',xlab = '',pch = seq(length(levels(pd$Cell))),
         main = 'Cell-specific Tcf4 Expression',las =2)
#legend('topright',legend = levels(pd$Cell),pch = seq(length(levels(pd$Cell))),cex = .6)
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


