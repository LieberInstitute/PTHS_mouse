### asd mouse differential expression analysis 
# qsub -V -l mf=10G,h_vmem=15G,h_stack=256M -cwd -b y R CMD BATCH analyze_asd_DESeq2.R
library(sva)
source('../DESeq2_functions.R') #work-horse of differential expression

##############################################
# load phenotype data and RPKM expression data
load('./rdas/pheno.rda',envir = dat<-new.env())
load('/dcl01/lieber/ajaffe/Brady/asd/geneCounts_asd_dec9_n204.rda')
load('rdas/qSVAs-geschwind_asd.rdas')
pd = cbind(pd,dat$pd)
all.equal(pd$SampleID,pd$SAMPLE_ID) #samples line up
rownames(pd) = pd$SAMPLE_ID
pd$Diagnosis = factor(ifelse(pd$Detailed.Diagnosis =="Chromosome 15q Duplication Syndrome",
                             'Dup15',pd$Diagnosis), levels = c("CTL", "ASD",'Dup15'))

#############################
# samples by region and age
with(pd, table(Diagnosis, Region))
with(pd, table(Diagnosis, Brain.Bank,Region))
with(pd, table(Diagnosis, Region,Sequencing.Batch))
with(pd, table(Diagnosis, Region,Brain.Bank))

#########################
# quality control metrics
library(beeswarm)
par(mar = c(2,5,3,2))
pdf('plots/asd_qc_plots.pdf')
boxplot(pd$totalAssignedGene~pd$Diagnosis+pd$Region,xaxt ='n',
        main = 'Gene Assignment Rate',ylab= 'Gene Assignment Rate')
axis(side = 1, at = 1:3*2-.5,labels = levels(factor(pd$Region)))
beeswarm(pd$totalAssignedGene~pd$Diagnosis+pd$Region,
         add = T,pch = 21,pwbg = pd$Diagnosis)
legend('topright',legend = c('CTL','ASD'),pch =21, pt.bg = 1:2)

boxplot(pd$mitoRate~pd$Diagnosis+pd$Region,xaxt ='n',
        main = 'Mito Rate', ylab = 'Mito Rate')
axis(side = 1, at = 1:3*2-.5,labels = levels(factor(pd$Region)))
beeswarm(pd$mitoRate~pd$Diagnosis+pd$Region,
         add = T,pch = 21,pwbg = pd$Diagnosis)
legend('topright',legend = c('CTL','ASD'),pch =21, pt.bg = 1:2)
dev.off()


########################
# split by brain regions
indList = split(seq(nrow(pd)), pd$Region)

##############################################
# create and run DESeq objects with covariates
geneDdsAdj <-lapply(indList,function(i){ 
  DESeq2(countData = geneCounts[,i], colData = pd[i,], 
         design = ~Diagnosis + totalAssignedGene +
           Sequencing.Batch + Brain.Bank + RIN + Age + Sex,
         sva = FALSE,parallel=TRUE)})

########################################
# create and run DESeq objects with qSVs
geneDdsQual <-lapply(indList,function(i){ 
  mod = model.matrix(~Diagnosis + totalAssignedGene +
  Sequencing.Batch + Brain.Bank + RIN + Age + Sex, data =pd[i,])
  qSVs = sva::qsva(degMat[,i], mod) # identify quality surrogate variables
  DESeq2(countData = geneCounts[,i], colData = cbind(pd[i,], qSVs), 
         design = as.formula(paste('~',paste(c('Diagnosis','totalAssignedGene',
      'Sequencing.Batch','Brain.Bank','RIN','Age','Sex',colnames(qSVs)), 
      collapse= "+"))), sva = FALSE,parallel=TRUE)})

################################################################
# get DE results stat adjusted, and fold-change ASD v. CTL patients
resGeneAdj <- c(lapply(geneDdsAdj,results,contrast = c('Diagnosis','ASD','CTL'), alpha=0.05),
                lapply(geneDdsAdj,results,contrast = c('Diagnosis','Dup15','CTL'), alpha=0.05))
names(resGeneAdj) = paste0(names(resGeneAdj),rep(c('_ASD','_Dup15'),each = 3))
sapply(resGeneAdj,function(g) sum(g$padj < 0.05, na.rm=TRUE))
#ba41-42-22_ASD          ba9_ASD       vermis_ASD  
#          1886             1036              322             
#ba41-42-22_Dup15     ba9_Dup15     vermis_Dup15 
#           3666            577               23 

resGeneQual <- c(lapply(geneDdsQual,results,contrast = c('Diagnosis','ASD','CTL'), alpha=0.05),
                 lapply(geneDdsQual,results,contrast = c('Diagnosis','Dup15','CTL'), alpha=0.05))
names(resGeneQual) = paste0(names(resGeneQual),rep(c('_ASD','_Dup15'),each = 3))
sapply(resGeneQual,function(g) sum(g$padj < 0.05, na.rm=TRUE))
#ba41-42-22_ASD          ba9_ASD       vermis_ASD 
#             5               99               10              
#ba41-42-22_Dup15 ba9_Dup15     vermis_Dup15 
#             654       423               87

outGeneList <- c(lapply(resGeneAdj,function(g) {
  outGene <- as.data.frame(g)
  outGene = outGene[order(outGene$padj,outGene$pvalue),]
  outGene = cbind(outGene,geneMap[rownames(outGene),])
}), lapply(resGeneQual,function(g) {
  outGene <- as.data.frame(g)
  outGene = outGene[order(outGene$padj,outGene$pvalue),]
  outGene = cbind(outGene,geneMap[rownames(outGene),])
}))
names(outGeneList) = paste0(names(outGeneList),rep(c('_Adj','_Qual'),each = 6))

sigGeneList = lapply(outGeneList,function(g){
  sigGene = g[which(g$pvalue<.01),]
})

pdf('plots/DESeq2_MA_plots_asd_Adj.pdf')
for (n in names(resGeneAdj)){
  plotMA(resGeneAdj[[n]], main=paste('Adj',n,"Gene MA plot"), ylim=c(-2,2))
}
for (n in names(resGeneQual)){
  plotMA(resGeneQual[[n]], main=paste('Qual',n,"Gene MA plot"), ylim=c(-2,2))
}
dev.off()

#################
# save everything
library(WriteXLS)
WriteXLS(sigGeneList, ExcelFileName ='tables/asd_DE_table_DESeq2_Adj.xls',row.names=T)
save(outGeneList,sigGeneList, file = 'rdas/asd_DE_objects_DESeq2_Adj.rda')
save(geneDdsAdj,geneDdsQual,resGeneAdj,resGeneQual, 
     file = '/dcl01/lieber/ajaffe/Brady/asd/asd_DESeq2_Adj.rda')



