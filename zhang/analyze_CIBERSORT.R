library(reshape2)
library(ggplot2)
options(stringsAsFactors = F)
######################
# load the cibersort and pca cell composition data
load('../tcf4_mouse/rdas/pheno.rda')
# CIBERSORT
dat = read.csv('tables/tcf4_mouse_CIBERSORT.csv')[,1:8]
# imputed from coverage matrix
dat2 = read.csv('../tcf4_mouse/tables/composition_estimate_cells.csv')[,c(1,19:25)]
dat2 = dat2[,c('X',"endothelial_cells", "neuron","Astrocyte", "oligodendrocyte_precursor_cells", 
               "newly_formed_oligodendrocytes","myelinating_oligodendrocytes","microglia")]
colnames(dat2) = colnames(dat)

dat3 = cbind(melt(dat,value.name = 'CIBERSORT',variable.name='Cell'),Jaffe_method=melt(dat2)$value)
dat3 = cbind(dat3,pd[dat3$Input.Sample,])

ggplot(data = dat3, aes(x = CIBERSORT,y = Jaffe_method,colour = Genotype,shape = Age)) + geom_point()+
  facet_wrap(~Cell,scales = 'free')


####################################
# mean expression of signature genes
library(DESeq2)
sig = read.delim('tables/sample_reference_zhang.rawGeneCounts_zhang.bm.K999.0.txt',row.names = 1)
load( '/dcl01/lieber/ajaffe/Brady/zhang/DESeq_zhang.rda', envir = zhang<-new.env())
load('/dcl01/lieber/ajaffe/Brady/mouseRNAseq/mega_dataset_DESeq2_svaAdj.rda',
     envir = mega<- new.env())
zhang$geneDds <- estimateSizeFactors(zhang$geneDds)


pd = data.frame(SAMPLE_ID = c(colData(zhang$geneDds)$SAMPLE_ID,colData(mega$geneDds)$SAMPLE_ID),
                Data = c(rep(c('zhang'),17),as.character(colData(mega$geneDds)$Line)))
nCounts = list(counts(zhang$geneDds,normalize = T),  counts(mega$geneDds,normalize = T))
labs = intersect(x = Reduce(intersect,lapply(nCounts,rownames)),y = rownames(sig))
nCounts = do.call(cbind,lapply(nCounts,function(x) x[labs,]))

pdf('plots/cell_type_specific_mean_expression.pdf')
  plot(colMeans(nCounts),main= 'Mean Expression',pch = 21, bg = factor(pd$Data),
          xlab = 'Dataset',ylab = 'Normalized Counts',log = 'y')
  legend('topright', legend = levels( factor(pd$Data)), pt.bg = seq(length(levels( factor(pd$Data)))),
         pch = 21)
dev.off()

################################
# stats

logit = function(p) log(p/(1-p))
indList = split(seq(nrow(pd)),pd$Age)


stat = sapply(indList,function(i){
  simplify2array(apply(dat[i,],2,function(x) t.test(x~pd$Genotype[i])))['statistic',]
})
p.value = sapply(indList,function(i){
  simplify2array(apply(dat[i,],2,function(x) t.test(x~pd$Genotype[i])))['p.value',]
})

padj = p.adjust(p.value,method= 'fdr')

