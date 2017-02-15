library(reshape2)
library(ggplot2)
library(nnet)
library(ordinal)

options(stringsAsFactors = F)
######################
# load the cibersort and pca cell composition data
load('../tcf4_mouse/rdas/mega_tcf4_pheno.rda')
pd$Region = ifelse(pd$Line =='Maher','mPFC',
            ifelse(pd$Line =='Sweatt','CA1','Wholebrain'))
dat = read.csv('tables/mega_tcf4_CIBERSORT.csv',row.names = 1)[,1:7]
id_vars = c('Genotype','Age','Line','FileID','Region')

###############################
# reformat into long data table
datLong = melt(cbind(dat,pd[rownames(dat),id_vars]),
               id.vars = id_vars,variable.name = 'Celltype',
               value.name = "Fraction")
datLong = datLong[datLong$Age != 'p21',]
datLong$Age = droplevels(datLong$Age)
datLong$Type = paste0(datLong$Celltype,'.',datLong$Age)

##############################################
# fit multinomial model for shift in p1 brains
indList = split(seq(nrow(datLong)),list(datLong$Celltype,datLong$Age))
coefs = sapply(indList, function(ii){
  tmp = summary(lm(formula = Fraction~Genotype +Line,data = datLong[ii,]))
  tmp$coefficients[2,'Estimate']
})
pvals = sapply(indList, function(ii){
  tmp = summary(lm(formula = Fraction~Genotype +Line,data = datLong[ii,]))
  tmp$coefficients[2,'Pr(>|t|)']
})
pvals[pvals<0.05]
coefs[pvals<0.05]
sigtype = names(coefs)[pvals<0.05]

postscript('plots/mega_tcf4_celltype_plots.eps',width=10,height = 2.5)
ggplot(aes(x = Region,y = Fraction,fill = Genotype),
       data=datLong[datLong$Type %in% sigtype,])+
  geom_boxplot()+scale_fill_manual(values = c('gray50','red'))+
  geom_point(position = position_jitterdodge())+
  facet_wrap(~Age+Celltype,scales = 'free',nrow = 1)+
  xlab('Tissue')+ylab('Proportion')
dev.off()















####################################
# mean expression of signature genes
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











datMulti = do.call('rbind',lapply(seq(nrow(dat)),function(ii){
  tmp = dat[ii,]
  tmp2 =data.frame(Celltype = rep(names(tmp), round(tmp*100)))
  tmp2$FileID= rownames(dat)[ii]
  return(tmp2)
}))
datMulti = cbind(datMulti,pd[datMulti$FileID,id_vars])
datMulti$Celltype = factor(datMulti$Celltype)
datMulti$Celltype = relevel(datMulti$Celltype, ref = 'Neuron')

##############################################
# fit multinomial model for shift in p1 brains
tmp = subset(datMulti,Age== 'p1')
tmp$Celltype = droplevels(tmp$Celltype)
testP1 <- multinom(Celltype ~ Genotype + Line,tmp)
sumTestP1 = summary(testP1)
z <- sumTestP1$coefficients[,2]/sumTestP1$standard.errors[,2] # genotype column
p <- (1 - pnorm(abs(z), 0, 1))*2
genoCoef = exp(coef(testP1)[,2])
genoCoef[p< 0.05]

##############################################
# fit multinomial model for shift in adult brains
testAdult <- multinom(Celltype ~ Genotype + Line, data = subset(datMulti,Age=='Adult'))
sumTestAdult = summary(testAdult)
z <- sumTestAdult$coefficients[,2]/sumTestAdult$standard.errors[,2] # genotype column
p <- (1 - pnorm(abs(z), 0, 1))*2
genoCoef = exp(coef(testAdult)[,2])
genoCoef[p< 0.05]

