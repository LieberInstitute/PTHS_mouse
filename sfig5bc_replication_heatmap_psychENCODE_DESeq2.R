# plot tcf4 mouse data with human ASD and 15qDuplication RNAseq
library(WriteXLS)
library(biomaRt)
library(pheatmap)
library(RColorBrewer)
library(jaffelab)
library(reshape2)
library(parallel)
library(ggplot2)
library("scales")
reverselog_trans <- function(base = exp(1)) {
  trans <- function(x) -log(x, base)
  inv <- function(x) base^(-x)
  trans_new(paste0("reverselog-", format(base)), trans, inv, 
            log_breaks(base = base), 
            domain = c(1e-100, Inf))
}

getOR = function(x) x[1,1]/x[2,1]/x[1,2]*x[2,2]
getKappa = function(x) (x[1,1]+x[2,2])/sum(x)


enrichGene = function(genes,genesets,univ){
  fishtest = data.frame(t(sapply(genesets,function(x)
    fisher.test(table(Query = univ %in% genes, Set = univ %in% x),alternative = 'g'))))
  fishtest = fishtest[,c('p.value','estimate')]
  fishtest$estimate = unlist(fishtest$estimate)
  fishtest$padj = p.adjust(fishtest$p.value,'fdr')
  fishtest$enriched = with(fishtest, padj< 0.05)
  fishtest$N = sapply(genesets,function(x)sum(univ %in% genes & univ %in% x))
  return(fishtest)
}

#################################
#get Ensembl mouse to human genes
ensembl=useMart("ensembl",dataset = 'mmusculus_gene_ensembl')
MMtoHG = getBM(attributes = c('ensembl_gene_id','hsapiens_homolog_ensembl_gene'),mart = ensembl)

###########################
#load human ASD/Dup15q data
load('asd/rdas/asd_DE_objects_DESeq2_Qual.rda',envir = asd <- new.env())
asd$outGeneList = lapply(asd$outGeneList[grep('Qual',names(asd$outGeneList))],function(g){
  tmpRowNames = ss(rownames(g),'\\.')
  ind = !duplicated(tmpRowNames)
  tmp = g[ind,]; rownames(tmp) = tmpRowNames[ind]
  return(tmp)
})

############################################
# check mouse tcf4 DEGs w/ psych encode data
load("tcf4_mouse/rdas/mega_tcf4_ages_DE_objects_DESeq2.rda") #mouse data
outGeneList = lapply(outGeneList, function(g) {
  g$hsapien_homolog = MMtoHG$hsapiens_homolog_ensembl_gene[match(rownames(g),MMtoHG$ensembl_gene_id)]
  g[grepl('ENSG',g$hsapien_homolog),]})

##########
# plot p1 genes
p1Genes =  with(outGeneList[['p1']], hsapien_homolog[which(padj<.01 & !is.na(padj))])
indP1 = sapply(p1Genes, function(n){ 
  tmp = sapply(asd$outGeneList, function(g) g[n,c('pvalue')])
  any(tmp<0.05)}) & !is.na(p1Genes)
p1Genes = p1Genes[indP1]

p1Dat = do.call('cbind',lapply(asd$outGeneList, function(g) g[p1Genes,c('log2FoldChange')]))
p1Dat = cbind(p1Dat,Mouse = with(outGeneList[['p1']], log2FoldChange[which(padj<.01 & !is.na(padj))])[indP1])
p1Dat[is.na(p1Dat)] = 0
tmp = colnames(p1Dat)
annotation = data.frame(Case = ifelse(grepl('Dup',tmp),'15qDup',
                                  ifelse(grepl('ASD',tmp),'ASD','PTHS')),
                    Region = ifelse(grepl('ba9',tmp),'Frontal Ctx',
                                    ifelse(grepl('ba41',tmp),'Temporal Ctx',
                                    ifelse(grepl('vermis',tmp),'Cerebellum','Mouse Brain'))),
                    row.names = tmp)
annoColor  = list(Case = c(ASD= '#a6cee3',`15qDup` ='#1f78b4',PTHS = '#e31a1c'),
                        Region = c(`Frontal Ctx` = 'black',`Temporal Ctx` = 'gray50',
                                   `Cerebellum` = 'lightgray',`Mouse Brain` = 'lightyellow'))

pdf('tcf4_mouse/plots/asd_heatmap_mega_tcf4_p1.pdf',height = 4,width = 4)
pheatmap(p1Dat,colorRampPalette(brewer.pal(n = 11, name ="PiYG"))(length(seq(-.9,.9,0.1))),
         border_color = NA,ylab = 'Log2 Fold-change',legend = F,labels_col = rep('',ncol(p1Dat)),
         labels_row = with(outGeneList[['p1']], Symbol[which(padj<.01 & !is.na(padj))])[indP1],
         breaks = seq(-.9,.9,0.1),annotation_col = annotation,annotation_color = annoColor)
dev.off()


##################
# plot adult genes
adultGenes =  with(outGeneList[['Adult']], hsapien_homolog[which(padj<.01 & !is.na(padj))])
indAdult = sapply(adultGenes, function(n){ 
  tmp = sapply(asd$outGeneList, function(g) g[n,c('pvalue')])
  sum(tmp<0.05,na.rm = T)>1}) & !is.na(adultGenes)
adultGenes = adultGenes[indAdult]

adultDat = do.call('cbind',lapply(asd$outGeneList, function(g) g[adultGenes,c('log2FoldChange')]))
adultDat = cbind(adultDat,Mouse = with(outGeneList[['Adult']], 
                                 log2FoldChange[which(padj<.01 & !is.na(padj))])[indAdult])
adultDat[is.na(adultDat)] = 0
tmp = colnames(adultDat)

pdf('tcf4_mouse/plots/asd_heatmap_mega_tcf4_adult.pdf',height = 7,width = 4.5)
pheatmap(adultDat,colorRampPalette(brewer.pal(n = 11, name ="PiYG"))(length(seq(-1.1,1.1,0.1))),
         border_color = NA,ylab = 'Log2 Fold-change',legend = F,labels_col = rep('',ncol(p1Dat)),
         labels_row = rep('',nrow(adultDat)),breaks = seq(-1.1,1.1,0.1),
         annotation_col = annotation,annotation_color = annoColor)
dev.off()



##############################
# replication enrichment tests
univ = unique(unlist(lapply(asd$outGeneList,rownames)))
genesets = lapply(outGeneList,function(g) g$hsapien_homolog[g$padj<.01 & !is.na(g$padj)])
dat = lapply(lapply(asd$outGeneList[c(1,4,2,5,3,6)],function(g) rownames(g)[g$pvalue<0.05]),
       enrichGene,genesets = genesets,univ = univ)

#####
# plot replication p-values
datLong = do.call('rbind',dat)
datLong$Region = factor(ss(rownames(datLong),'_'),levels = c('ba9','ba41-42-22','vermis'),
                        labels = c('Frontal Ctx','Temporal Ctx','Cerebellum'))
datLong$Diagnosis = ss(rownames(datLong),'_',2)
datLong$Age = factor(ss(rownames(datLong),'\\.',2),levels = c('p1','Adult'),
                     labels = c('P1 DEGs','Adult DEGs'))

pdf('tcf4_mouse/plots/asd_mega_tcf4_enrichment_barplot.pdf',height = 2.5,width = 2.5)
ggplot(data =subset(datLong,Age=='Adult DEGs'),aes(x = Region,y = padj,fill = Diagnosis))+
  geom_bar(stat = "identity",position = "dodge") +# facet_wrap(~Age,scales = 'free')+
  scale_y_continuous(trans=reverselog_trans(10))+
  scale_fill_manual(values = c('gray','black'),guide=FALSE) +
  geom_hline(yintercept = 0.05,linetype = 2,color = 'red')
dev.off()

datLong = datLong[datLong$Age =='Adult DEGs',]
rownames(datLong)=ss(rownames(datLong),'\\.')

#########################
# logFC concordance rates, correlations, enrichment odds ratios
if (TRUE){
nullKappas = lapply(asd$outGeneList,function(g){
  ecdf(unlist(lapply(1:1000,function(i){
    cat(i)
    g[,1:6] = apply(g[,1:6],2,sample)
    tmp = outGeneList[['Adult']]
    tmp[,1:6] = apply(tmp[,1:6],2,sample)
    tmp = with(tmp,tmp[padj<.01,])
    tmp2 = with(tmp,cbind(M = log2FoldChange,g[hsapien_homolog,c('log2FoldChange','pvalue')]))
    with(tmp2[tmp2$pvalue<0.05,], getKappa(table(M>0,log2FoldChange>0)))
  })))
})
save(nullKappas,file = 'tcf4_mouse/rdas/psychENCODE_asd_nullKappas.rda')
} else{
  load('tcf4_mouse/rdas/psychENCODE_asd_nullKappas.rda')
}

datLong$kappas =  sapply(rownames(datLong),function(n){
  g = asd$outGeneList[[n]]
  tmp = outGeneList[['Adult']]
  tmp = with(tmp,tmp[padj<.01,])
  tmp2 = with(tmp,cbind(M = log2FoldChange,g[hsapien_homolog,c('log2FoldChange','pvalue')]))
  with(tmp2[tmp2$pvalue<0.05,], getKappa(table(M>0,log2FoldChange>0)))})

datLong$kappaPvalue =  sapply(rownames(datLong),function(n){
  tmp = nullKappas[[n]]
  1-tmp(datLong[n,'kappas'])
  })

datLong$rhos = sapply(asd$outGeneList,function(g){
  tmp = outGeneList[['Adult']]
  tmp = with(tmp,tmp[padj<.01,])
  tmp2 = with(tmp,cbind(M = log2FoldChange,g[hsapien_homolog,c('log2FoldChange','pvalue')]))
  tmp2 = tmp2[complete.cases(tmp2),]
  with(tmp2[tmp2$pvalue<0.05,], cor(x =M,y = log2FoldChange,method = 'spearman'))})

datLong$corPval = p.adjust(sapply(asd$outGeneList,function(g){
  tmp = outGeneList[['Adult']]
  tmp = with(tmp,tmp[padj<.01,])
  tmp2 = with(tmp,cbind(M = log2FoldChange,g[hsapien_homolog,c('log2FoldChange','pvalue')]))
  tmp2 = tmp2[complete.cases(tmp2),]
  with(tmp2[tmp2$pvalue<0.05,], cor.test(x =M,y = log2FoldChange,method = 'spearman')$p.value)}),'fdr')
datLong[datLong$corPval<0.05,]

pdf('tcf4_mouse/plots/asd_tcf4_concordance_plots.pdf',height = 2.5,width = 2.5)
ggplot(data = datLong,aes(x = Region,y = kappas,fill = Diagnosis))+
  geom_bar(stat = "identity",position = "dodge") + ylim(c(0,1))+
  scale_fill_manual(values = c('gray','black'),guide=FALSE) + 
  geom_hline(yintercept = 0.5,linetype = 2,color = 'red')

#ggplot(data = datLong,aes(x = Region,y = rhos,fill = Diagnosis))+
#  geom_bar(stat = "identity",position = "dodge") + ylim(c(-1,1))+
#  scale_fill_manual(values = c('gray','black')) + 
#  geom_hline(yintercept = 0,linetype = 2,color = 'red')
dev.off()
save(datLong,file = 'tcf4_mouse/rdas/psychENCODE_asd_enrichment_concordance.rda')

########################
# make venndiagram plots
library(gplots)
genesets = lapply(outGeneList,function(g) g$hsapien_homolog[g$padj<.01 & !is.na(g$padj)])
Regions = c('ASD','Dup15')
pdf('tcf4_mouse/plots/asd_tcf4_overlap_venndiagrams.pdf')
for (reg in Regions){
#  tmpList = c(genesets[1],lapply(asd$outGeneList[grep(reg,names(asd$outGeneList))],
#                                   function(g) rownames(g)[g$pvalue<0.01]))
#  geneVenn = plot(venn(tmpList,small =1.5,show.plot = F),cex = 2,main = paste(reg,'p1'))
  
  tmpList = c(genesets[2],lapply(asd$outGeneList[grepl(reg,names(asd$outGeneList)) & 
                                                   !grepl('vermis',names(asd$outGeneList))],
                                 function(g) rownames(g)[g$pvalue<0.05]))
  names(tmpList) = c("Tcf4 Mouse", "Temporal Ctx", "Frontal Ctx")
  geneVenn = plot(venn(tmpList,small =1.5,show.plot = F),cex = 2,main = paste(reg,'Adult'))
}
dev.off()



