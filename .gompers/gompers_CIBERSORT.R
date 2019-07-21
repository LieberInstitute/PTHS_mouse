library(DESeq2)
library(beeswarm)
library(jaffelab)
library(ggplot2)
library(reshape2)

###############
# load the data
load('rdas/pheno_gompers.rda',envir = dat <-new.env())
pd = subset(pd, Age %in% c('P0','P56'))
load('/dcl01/lieber/ajaffe/Brady/chd8_het/rawCounts_Gompers_Chd8_het_n19.rda')
pd = cbind(metrics,dat$pd)
table(pd$SAMPLE_ID==pd$SampleID)
id_vars = c('Genotype','Age')

############################################
# format data for CIBERSORT signature matrix
ind = which(!duplicated(ss(rownames(geneCounts),'\\.')))
geneCounts = geneCounts[ind, ]
rownames(geneCounts) = ss(rownames(geneCounts),'\\.')
write.table(geneCounts,sep = '\t',quote = FALSE,file = 'tables/rawGeneCounts_gompers.txt')

############################################
# load CIBERSORTed samples  
dat = read.csv('tables/chd8_gompers_CIBERSORT.csv',row.names = 1)[,1:7]

###############################
# reformat into long data table
datLong = melt(cbind(dat,pd[rownames(dat),id_vars]),
               id.vars = id_vars,variable.name = 'Celltype',
               value.name = "Fraction")
datLong$Type = paste0(datLong$Celltype,'.',datLong$Age)

pdf('plots/gompers_chd8_cibersort.pdf')
ggplot(data = datLong, aes(x = Age, y = Fraction, fill = Genotype)) +
  geom_boxplot(position = 'dodge') + facet_wrap(~Celltype, scales="free_y") + 
  ggtitle('CIBERSORTed Chd8 mouse celltypes') + 
  ylab('Percent')+ xlab('Age')
dev.off()

