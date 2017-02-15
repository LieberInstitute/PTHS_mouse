# plot tcf4 mouse data with human ASD and 15qDuplication RNAseq
library(WriteXLS)
library(biomaRt)
library(pheatmap)
library(RColorBrewer)
library(jaffelab)
library(reshape2)
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
  return(fishtest)
}

#################################
#get Ensembl mouse to human genes
ensembl=useMart("ensembl",dataset = 'mmusculus_gene_ensembl')
MMtoHG = getBM(attributes = c('ensembl_gene_id','hsapiens_homolog_ensembl_gene'),mart = ensembl)

###########################
#load human ASD/Dup15q data
load('asd/rdas/asd_DE_objects_DESeq2_Adj.rda',envir = asd <- new.env())
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
p1Genes =  with(outGeneList[['p1']], hsapien_homolog[which(padj<.05 & !is.na(padj))])
indP1 = sapply(p1Genes, function(n){ 
  tmp = sapply(asd$outGeneList, function(g) g[n,c('pvalue')])
  any(tmp<0.05)}) & !is.na(p1Genes)
p1Genes = p1Genes[indP1]

p1Dat = do.call('cbind',lapply(asd$outGeneList, function(g) g[p1Genes,c('log2FoldChange')]))
p1Dat = cbind(p1Dat,Mouse = with(outGeneList[['p1']], log2FoldChange[which(padj<.05 & !is.na(padj))])[indP1])
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
         labels_row = with(outGeneList[['p1']], Symbol[which(padj<.05 & !is.na(padj))])[indP1],
         breaks = seq(-.9,.9,0.1),annotation_col = annotation,annotation_color = annoColor)
dev.off()


##################
# plot adult genes
adultGenes =  with(outGeneList[['Adult']], hsapien_homolog[which(padj<.05 & !is.na(padj))])
indAdult = sapply(adultGenes, function(n){ 
  tmp = sapply(asd$outGeneList, function(g) g[n,c('pvalue')])
  sum(tmp<0.05,na.rm = T)>1}) & !is.na(adultGenes)
adultGenes = adultGenes[indAdult]

adultDat = do.call('cbind',lapply(asd$outGeneList, function(g) g[adultGenes,c('log2FoldChange')]))
adultDat = cbind(adultDat,Mouse = with(outGeneList[['Adult']], 
                                 log2FoldChange[which(padj<.05 & !is.na(padj))])[indAdult])
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
genesets = lapply(outGeneList,function(g) g$hsapien_homolog[g$padj<0.05 & !is.na(g$padj)])
dat = lapply(lapply(asd$outGeneList[c(1,4,2,5,3,6)],function(g) rownames(g)[g$pvalue<0.01]),
       enrichGene,genesets = genesets,univ = univ)

#########################
# logFC concordance rates, correlations, enrichment odds ratios
pvalThresh = 10^-(seq(from = 0.1,to = 2,by = .1))

ors = data.frame(pval = pvalThresh, sapply(asd$outGeneList,function(g){
  tmp = outGeneList[['Adult']]
  tmp2 = with(tmp,cbind(M = log2FoldChange,g[hsapien_homolog,c('log2FoldChange','pvalue')]))
  sapply(pvalThresh,function(p) getOR(table(tmp$padj<0.05,tmp2$pvalue<p)))
}))
orsLong = melt(ors,id.var = c('pval'))
orsLong$Region = ss(as.character(orsLong$variable),'_')
orsLong$Diagnosis = factor(ss(as.character(orsLong$variable),'_',2))

kappas =  data.frame(pval = pvalThresh, sapply(asd$outGeneList,function(g){
  tmp = outGeneList[['Adult']]
  tmp = with(tmp,tmp[padj<0.05,])
  tmp2 = with(tmp,cbind(M = log2FoldChange,g[hsapien_homolog,c('log2FoldChange','pvalue')]))
  sapply(pvalThresh,function(p) with(tmp2[tmp2$pvalue<p,], getKappa(table(M>0,log2FoldChange>0))))
}))
kappaLong = melt(kappas,id.var = c('pval'))
kappaLong$Region = ss(as.character(kappaLong$variable),'_')
kappaLong$Diagnosis = factor(ss(as.character(kappaLong$variable),'_',2))

rhos = data.frame(pval = pvalThresh, sapply(asd$outGeneList,function(g){
  tmp = outGeneList[['Adult']]
  tmp = with(tmp,tmp[padj<0.05,])
  tmp2 = with(tmp,cbind(M = log2FoldChange,g[hsapien_homolog,c('log2FoldChange','pvalue')]))
  tmp2 = tmp2[complete.cases(tmp2),]
  sapply(pvalThresh,function(p) with(tmp2[tmp2$pvalue<p,], 
                                           cor(x =M,y = log2FoldChange,method = 'spearman')))
}))
rhoLong = melt(rhos,id.var = c('pval'))
rhoLong$Region = ss(as.character(rhoLong$variable),'_')
rhoLong$Diagnosis = factor(ss(as.character(rhoLong$variable),'_',2))


pdf('tcf4_mouse/plots/asd_tcf4_concordance_plots.pdf',height = 4,width = 4)
ggplot(data = orsLong,aes(x = pval,y = value,colour = Region,linetype = Diagnosis))+
  geom_line(cex = 1)+geom_vline(xintercept = 0.05,colour = 'red',linetype =3)+
  xlab('P-value threshold')+ggtitle('Odds Ratios')+scale_y_log10()+
  scale_x_continuous(trans=reverselog_trans(10))+scale_linetype_manual(values = c(2,1))+
  scale_color_brewer(palette = 'Set1')+ylab('')

ggplot(data = kappaLong,aes(x = pval,y = value,colour = Region,linetype = Diagnosis))+
  geom_line(cex = 1)+geom_vline(xintercept = 0.05,colour = 'red',linetype =3)+
  xlab('P-value threshold')+ggtitle('Concordance Rate')+
  scale_x_continuous(trans=reverselog_trans(10))+scale_linetype_manual(values = c(2,1))+
  scale_color_brewer(palette = 'Set1')+ylab('')

ggplot(data = rhoLong,aes(x = pval,y = value,colour = Region,linetype = Diagnosis))+
  geom_line(cex = 1)+geom_vline(xintercept = 0.05,colour = 'red',linetype =3)+
  xlab('P-value threshold')+ggtitle('Correlation')+ylim(c(-max(rhoLong$value),max(rhoLong$value)))+
  scale_x_continuous(trans=reverselog_trans(10))+scale_linetype_manual(values = c(2,1))+
  scale_color_brewer(palette = 'Set1')+ylab('')
dev.off()

########################
# make venndiagram plots
library(gplots)
genesets = lapply(outGeneList,function(g) g$hsapien_homolog[g$padj<0.05 & !is.na(g$padj)])
Regions = c('ba9','ba41-42-22','vermis')
pdf('tcf4_mouse/plots/asd_tcf4_overlap_venndiagrams.pdf')
for (reg in Regions){
  tmpList = c(genesets[1],lapply(asd$outGeneList[grep(reg,names(asd$outGeneList))],
                                    function(g) rownames(g)[g$pvalue<0.05]))
  geneVenn = plot(venn(tmpList,small =1.5,show.plot = F),cex = 2,main = paste(reg,'p1'))
  
  tmpList = c(genesets[2],lapply(asd$outGeneList[grep(reg,names(asd$outGeneList))],
                                 function(g) rownames(g)[g$pvalue<0.05]))
  geneVenn = plot(venn(tmpList,small =1.5,show.plot = F),cex = 2,main = paste(reg,'Adult'))
}
dev.off()



#####################################
# gene set analysis in these datasets
library(clusterProfiler)
library(org.Hs.eg.db)

########################################################
# background is expressed human genes with mouse homologs
univ = as.character(out[outGene$hsapien_homolog,c('EntrezID')]) 
gList = split(out[dat$hsapien_homolog,c('EntrezID')],dat$type)
sapply(gList, length)

############################
# enrich all the GO terms :D
compareGoMf = compareCluster(gList, fun = "enrichGO", universe = univ, ont = "MF",OrgDb=org.Hs.eg.db, 
                             qvalueCutoff = 0.05, pvalueCutoff = 0.05, readable = TRUE)
compareGoBp = compareCluster(gList, fun = "enrichGO", universe = univ, ont = "BP",OrgDb=org.Hs.eg.db, 
                             qvalueCutoff = 0.05, pvalueCutoff = 0.05, readable = TRUE)
compareGoCc = compareCluster(gList, fun = "enrichGO", universe = univ, ont = "CC",OrgDb=org.Hs.eg.db, 
                             qvalueCutoff = 0.05, pvalueCutoff = 0.05, readable = TRUE)
compareGo = lapply(list(compareGoMf=compareGoMf,compareGoBp=compareGoBp,compareGoCc=compareGoCc),simplify)

###################
# save the go table
goTable = lapply(list(compareGoMf=compareGoMf,compareGoBp=compareGoBp,compareGoCc=compareGoCc),as.data.frame)
save(compareGoMf,compareGoBp,compareGoCc,file = 'tcf4_mouse/rdas/psychENCODE_human_asd_GO_terms.rda')
WriteXLS(goTable,ExcelFileName = 'tcf4_mouse/tables/stable8_psychENCODE_human_asd_GO_compareCluster.xls')
