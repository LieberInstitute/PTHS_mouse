library(reshape2)
library(ggplot2)
library(nnet)

options(stringsAsFactors = F)
#######################################
# load the cibersort and phenotype data
load('rdas/pheno.rda')
pd$Diagnosis2 = factor(ifelse(pd$Detailed.Diagnosis =="Chromosome 15q Duplication Syndrome",
                             'Dup15',pd$Diagnosis), levels = c("CTL", "ASD",'Dup15'))
pd$Diagnosis = factor(pd$Diagnosis,levels = c('CTL',"ASD"))
#pd$Region = factor(pd$Region,levels = c('ba9','ba41-42-22','vermis'),
 #                  labels = c('Frontal','Temporal','Cerebellum'))


dat = read.delim('tables/asd_CIBERSORT_cell_population.txt',row.names = 1)[,1:7]
id_vars = c('Diagnosis','Diagnosis2','Sequencing.Batch','Region',
            'Brain.Bank', 'RIN', 'Age', 'Sex')

###############################
# reformat into long data table
datLong = melt(cbind(dat,pd[match(rownames(dat),pd$SampleID),id_vars]),
               id.vars = id_vars,variable.name = 'Celltype',
               value.name = "Fraction")
datLong$Type = paste0(datLong$Celltype,'.',datLong$Region)

##############################################
# fit multinomial model for shift in p1 brains
indList = split(seq(nrow(datLong)),list(datLong$Celltype,datLong$Region))
coefs = sapply(indList, function(ii){
  tmp = summary(lm(formula = Fraction~Diagnosis+Sequencing.Batch+Brain.Bank + 
                     RIN+Age+Sex,data = datLong[ii,]))
  tmp$coefficients[2,'Estimate']
})
pvals = sapply(indList, function(ii){
  tmp = summary(lm(formula = Fraction~Diagnosis+Sequencing.Batch+Brain.Bank + 
                     RIN+Age+Sex,data = datLong[ii,]))
  tmp$coefficients[2,'Pr(>|t|)']
})
pvals[pvals<0.05]
coefs[pvals<0.05]
sigtype = names(coefs)[pvals<0.05]

if(FALSE){
coefs = sapply(indList, function(ii){
  tmp = summary(lm(formula = Fraction~Diagnosis+Sequencing.Batch+
                     Brain.Bank+RIN+Age+Sex,data = datLong[ii,]))
  tmp$coefficients[3,'Estimate']
})
pvals = sapply(indList, function(ii){
  tmp = summary(lm(formula = Fraction~Diagnosis+Sequencing.Batch+
                     Brain.Bank+RIN+Age+Sex,data = datLong[ii,]))
  tmp$coefficients[3,'Pr(>|t|)']
})
pvals[pvals<0.05]
coefs[pvals<0.05]
sigtype = unique(c(sigtype),names(coefs)[pvals<0.05])
}

#######
# ANOVA
indList = split(seq(nrow(datLong)),list(datLong$Celltype))
pvals = sapply(indList, function(ii){
  tmp = anova(lm(formula = Fraction~Diagnosis2+Region+Sequencing.Batch+Brain.Bank + 
                     RIN+Age+Sex,data = datLong[ii,]))
  tmp[1,'Pr(>F)']
})
pvals[pvals<0.05]


pvals = sapply(indList, function(ii){
  tmp = summary(lm(formula = Fraction~Diagnosis2+Region+Sequencing.Batch+Brain.Bank + 
                     RIN+Age+Sex,data = datLong[ii,]))
  tmp$coefficients[2,'Pr(>|t|)']
})
pvals[pvals<0.05]




sigtype = names(pvals[pvals<0.05])
############
# make plots
pdf('plots/asd_celltypes_CIBERSORT_ast_new_pre_micro.pdf',width = 9,height = 3)
ggplot(data=datLong[datLong$Celltype %in% sigtype,], 
       aes(x=Region,y = Fraction,fill = Diagnosis2))+
  geom_boxplot(position = position_dodge(width = .75),outlier.shape=NA) +
  #geom_point(pch=21, position=position_jitterdodge(dodge.width = .75)) +
  facet_wrap(~Celltype,nrow = 1, scales = 'free')+xlab('Age')+
  ylab('Proportion') + scale_fill_manual(values = c('white','gray','black'))
dev.off()







datMulti = do.call('rbind',lapply(seq(nrow(dat)),function(ii){
  tmp = dat[ii,]
  tmp2 =data.frame(Celltype = rep(names(tmp), round(tmp*100)))
  tmp2$FileID= rownames(dat)[ii]
  return(tmp2)
}))
datMulti = cbind(datMulti,pd[match(datMulti$FileID,pd$SampleID),id_vars])
datMulti$Celltype = factor(datMulti$Celltype)
datMulti$Celltype = relevel(datMulti$Celltype,ref = 'Endothelial')

##############################################
# fit multinomial model for shift in adult brains
tmp = subset(datMulti,Region =='vermis')
indList = split(seq(nrow(datMulti)),list(datMulti$Region))
test <- lapply(indList,function(ii) {
  multinom(Celltype ~ Diagnosis+Sequencing.Batch+
             Brain.Bank+RIN+Age+Sex, data = datMulti[ii,])})
sumTest <- lapply(indList,function(ii) {
  summary(multinom(Celltype ~ Diagnosis+Sequencing.Batch+
             Brain.Bank+RIN+Age+Sex, data = datMulti[ii,]))})
z <- lapply(sumTest, function(x) x$coefficients[,2]/x$standard.errors[,2]) # genotype column
p <- lapply(z,function(x) (1 - pnorm(abs(x), 0, 1))*2)
  
coef = lapply(test,function(x) exp(coef(x)[,2]))
lapply(names(z), function(x) coef[[x]][p[[x]]<0.05])

pdf('plots/asd_celltypes_CIBERSORT.pdf',width = 10,height = 3)
ggplot(data=datLong[datLong$Celltype%in%sigtype &datLong$Region=='ba9',],
       aes(x=Diagnosis,y = Fraction,fill = Diagnosis))+
  geom_boxplot(position = position_dodge(width = .5),outlier.shape=NA) +
  geom_point(pch=21, position=position_jitterdodge(dodge.width = .5)) +
  facet_wrap(~Region+Celltype,nrow=1,scales = 'free')+xlab('Age')+
  ylab('Proportion') +scale_fill_manual(values = c("lightyellow","red", "blue"))
dev.off()


datLong$Type = paste0(datLong$Celltype,'.',datLong$Region)

##############################################
# fit multinomial model for shift in p1 brains
coefs = sapply(indList, function(ii){
  tmp = summary(lm(formula = Fraction~Diagnosis+Sequencing.Batch+Brain.Bank + 
                     RIN+Age+Sex,data = datLong[ii,]))
  tmp$coefficients[2,'Estimate']
})
pvals = sapply(indList, function(ii){
  tmp = summary(lm(formula = Fraction~Diagnosis+Sequencing.Batch+Brain.Bank + 
                     RIN+Age+Sex,data = datLong[ii,]))
  tmp$coefficients[2,'Pr(>|t|)']
})
