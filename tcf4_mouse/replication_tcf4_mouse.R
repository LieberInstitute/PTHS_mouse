# find replication in sweatt dataset and philpot dataset
library(WriteXLS)

###############################################
# load discovery genes and subset to fdr < 0.05
load('rdas/mouse_tcf4_DE_objects_DESeq2.rda',envir = discovery<- new.env())
discovery = subset(discovery$outGene,padj< 0.05 & !is.na(padj))
query = discovery$Symbol
names(discovery) = paste0(names(discovery), '_Discovery')

###############################################
# load replicate dataset philpot and sweatt
load('../sweatt/rdas/sweatt_DE_objects_DESeq2.rda',envir = sweatt<- new.env())
load('../philpot/rdas/philpot_DE_objects_DESeq2.rda',envir = philpot<- new.env())

#rename sweatt data
names(sweatt$outGene) = paste0(names(sweatt$outGene),'_CA1')

#combine datasets sweatt with philpot
outGene = cbind(philpot$outGene, sweatt$outGene[rownames(philpot$outGene),1:6])

#find each gene in other data
gIndex = unlist(sapply(paste0("^",query,"$"), grep, outGene$Symbol))
sigGene = outGene[gIndex,]
sigGene = cbind(discovery,sigGene[,8:37])

##############################################################
# determine if genes that replicated in all or any other model
Line = c("Act",'Nest', 'Del', 'R579W', 'CA1')
sigGene$replicatedAny = sapply(rownames(sigGene),function(g) {
  repPval = sigGene[g,paste0('pvalue_',Line)]< 0.05
  repFC = sign(sigGene[g,paste0('log2FoldChange_',Line)])== sign(sigGene[g,'log2FoldChange_Discovery'])
  return(any(repPval & repFC))})
sigGene$replicatedAll = sapply(rownames(sigGene),function(g) {
  repPval = sigGene[g,paste0('pvalue_',Line)]< 0.05
  repFC = sign(sigGene[g,paste0('log2FoldChange_',Line)])== sign(sigGene[g,'log2FoldChange_Discovery'])
  return(all(repPval & repFC))})
sigGene$numReplicated = sapply(rownames(sigGene),function(g) {
  repPval = sigGene[g,paste0('pvalue_',Line)]< 0.05
  repFC = sign(sigGene[g,paste0('log2FoldChange_',Line)])== sign(sigGene[g,'log2FoldChange_Discovery'])
  return(sum(repPval & repFC))})
sum(sigGene$replicatedAny) #36 
sum(sigGene$replicatedAll) #7

#########################
# save for plotting sfig3
save(sigGene,file = "rdas/mouse_tcf4_replication.rda")

############################
# replication by mouse model
apply(sigGene[,grep('pvalue',names(sigGene))],2,function(x) sum(x<0.05))

############################
# genes that replicated in all 
# Mc4r, Atp2b4, Tmem44, Enpp6, Tmem88b, Gm4211, and Mog
sigGene$Symbol_Discovery[sigGene$replicatedAll] 

#############################################################
# genes that replicated in any and how many models replicated 
sigGene[sigGene$replicatedAny,c('Symbol_Discovery','numReplicated')]
sum(sigGene$numReplicated>=4)

###########################
# save replicated gene list!
WriteXLS(sigGene, ExcelFileName = "tables/stable3_DE_feature_replication.xls")
