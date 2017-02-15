library(reshape2)
library(ggplot2)
library(nnet)

options(stringsAsFactors = F)
#######################################
# load the cibersort and phenotype data
load('rdas/pheno.rda')

dat = read.delim('tables/celltype_pten_CIBERSORT.txt',row.names = 1)[,1:7]
id_vars = c('Genotype', 'Age')

###############################
# reformat into long data table
datLong = melt(cbind(dat,pd[match(rownames(dat),pd$FileID),id_vars]),
               id.vars = id_vars,variable.name = 'Celltype',
               value.name = "Fraction")
datLong$Type = paste0(datLong$Celltype,'.',datLong$Age)

##############################################
# fit multinomial model for shift in p1 brains
indList = split(seq(nrow(datLong)),list(datLong$Celltype,datLong$Age))
coefs = sapply(indList, function(ii){
  tmp = summary(lm(formula = Fraction~Genotype,data = datLong[ii,]))
  tmp$coefficients[2,'Estimate']
})
pvals = sapply(indList, function(ii){
  tmp = summary(lm(formula = Fraction~Genotype,data = datLong[ii,]))
  tmp$coefficients[2,'Pr(>|t|)']
})
pvals[pvals<0.05]
coefs[pvals<0.05]
sigtype = names(coefs)[pvals<0.05]


coefs = sapply(indList, function(ii){
  tmp = summary(lm(formula = Fraction~Genotype,data = datLong[ii,]))
  tmp$coefficients[3,'Estimate']
})
pvals = sapply(indList, function(ii){
  tmp = summary(lm(formula = Fraction~Genotype,data = datLong[ii,]))
  tmp$coefficients[3,'Pr(>|t|)']
})
pvals[pvals<0.05]
coefs[pvals<0.05]
sigtype = unique(names(coefs)[pvals<0.05])

############
# make plots
pdf('plots/pten_celltypes_CIBERSORT.pdf',width = 10,height = 10)
ggplot(data=datLong[datLong$Type%in%sigtype ,],
       aes(x=Genotype,y = Fraction,fill = Genotype))+
  geom_boxplot(position = position_dodge(width = .5),outlier.shape=NA) +
  geom_point(pch=21, position=position_jitterdodge(dodge.width = .5)) +
  facet_wrap(~Celltype,nrow=1,scales = 'free')+xlab('Age')+
  ylab('Proportion') +scale_fill_manual(values = c("gray","pink", "red"))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(data=datLong,aes(x=Genotype,y = Fraction,fill = Genotype))+
  geom_boxplot(position = position_dodge(width = .5),outlier.shape=NA) +
  geom_point(pch=21, position=position_jitterdodge(dodge.width = .5)) +
  facet_wrap(~Celltype,scales = 'free')+xlab('Age')+
  ylab('Proportion') +scale_fill_manual(values = c("gray","pink", "red"))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()






##############################################
# fit multinomial model for shift in adult brains
datMulti = do.call('rbind',lapply(seq(nrow(dat)),function(ii){
  tmp = dat[ii,]
  tmp2 =data.frame(Celltype = rep(names(tmp), round(tmp*100)))
  tmp2$FileID= rownames(dat)[ii]
  return(tmp2)
}))
datMulti = cbind(datMulti,pd[match(datMulti$FileID,pd$FileID),id_vars])
datMulti$Celltype = factor(datMulti$Celltype)
datMulti$Celltype = relevel(datMulti$Celltype,ref = 'Endothelial')

##############################################
# fit multinomial model for shift in adult brains
indList = split(seq(nrow(datMulti)),list(datMulti$Age))
test <- lapply(indList,function(ii){
  multinom(Celltype ~ Genotype, data = datMulti[ii,])})
sumTest <- lapply(indList,function(ii){
  summary(multinom(Celltype ~ Genotype, data = datMulti[ii,]))})
z <- lapply(sumTest, function(x) x$coefficients[,2]/x$standard.errors[,2]) # genotype column
p <- lapply(z,function(x) (1 - pnorm(abs(x), 0, 1))*2)

coef = lapply(test,function(x) exp(coef(x)[,2]))
lapply(names(z), function(x) coef[[x]][p[[x]]<0.05])

pdf('plots/pten_celltypes_CIBERSORT.pdf',width = 8,height = 3)
ggplot(data=datLong[datLong$Celltype=='Mye.Oligo',],
       aes(x=Genotype,y = Fraction,fill = Genotype))+
  geom_boxplot(position = position_dodge(width = .5),outlier.shape=NA) +
  geom_point(pch=21, position=position_jitterdodge(dodge.width = .5)) +
  facet_wrap(~Age+Celltype,nrow=1,scales = 'free')+xlab('Age')+
  ylab('Proportion') +scale_fill_manual(values = c("gray","pink", "red"))
dev.off()









