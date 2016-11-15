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
names(sweatt$outExon) = paste0(names(sweatt$outExon),'_CA1')
names(sweatt$outJxn) = paste0(names(sweatt$outJxn),'_CA1')

#combine datasets sweatt with philpot
outGene = cbind(philpot$outGene, sweatt$outGene[rownames(philpot$outGene),1:6])
e1 = paste(philpot$outExon$Chr, philpot$outExon$Start, philpot$outExon$End)
e2 = paste(sweatt$outExon$Chr, sweatt$outExon$Start, sweatt$outExon$End)
outExon = cbind(philpot$outExon, sweatt$outExon[match(e1,e2),1:6])
outJxn = cbind(philpot$outJxn, sweatt$outJxn[rownames(philpot$outJxn),1:6])

#find each gene/exon/junction replicated in other data
gIndex = unlist(sapply(paste0("^",query,"$"), grep, outGene$Symbol))
eIndex = unlist(sapply(paste0("^",query,"$"), grep, outExon$Symbol))
jIndex = unlist(sapply(paste0("^",query,"$"), grep, outJxn$ensemblSymbol))

sigGene = outGene[gIndex,]
sigGene = cbind(discovery,sigGene[,8:37])
sigExon = outExon[eIndex,]
sigJxn = outJxn[jIndex,]

#########################
# save for plotting sfig3
save(sigGene,sigExon,sigJxn,file = "rdas/mouse_tcf4_replication.rda")

#######################################################
# filter to gene/exon/junction replicated in any genotype
sigGene = sigGene[apply(
  sigGene[,grep('pvalue',names(sigGene))],1,
  function(x) any(x < 0.05)),]
sigExon = sigExon[apply(
  sigExon[,grep('pvalue',names(sigExon))],1,
  function(x) any(x < 0.05)),]
sigJxn = sigJxn[apply(
  sigJxn[,grep('pvalue',names(sigJxn))],1,
  function(x) any(x < 0.05)),]

############################
# replication by mouse model
apply(sigGene[,grep('pvalue',names(sigGene))],2,function(x) sum(x<0.05))
apply(sigExon[,grep('pvalue',names(sigExon))],2,function(x){
  suppressWarnings(sum(tapply(x,sigExon$Symbol,min,na.rm = TRUE)>0.05))
})
apply(sigJxn[,grep('pvalue',names(sigJxn))],2,function(x) {
 suppressWarnings(sum(tapply(x,sigJxn$ensemblSymbol,min,na.rm = TRUE)>0.05))
})

###########################
# save replicated gene list!
WriteXLS(list(Gene = sigGene,Exon = sigExon, Junction = sigJxn),
         ExcelFileName = "tables/stable4_DE_feature_replication.xls")
