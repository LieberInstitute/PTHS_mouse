# check TCF4 DEGs with human ASD from Psych ENCODE dataset
getOR = function(x) x[1,1]/x[2,1]/x[1,2]*x[2,2]
getKappa = function(x) (x[1,1]+x[2,2])/sum(x)

library(WriteXLS)
library(biomaRt)
library(data.table)
library(ggplot2)
library(jaffelab)
library(irr)

#################################
#get Ensembl mouse to human genes
ensembl=useMart("ensembl",dataset = 'mmusculus_gene_ensembl')
MMtoHG = getBM(attributes = c('ensembl_gene_id','hsapiens_homolog_ensembl_gene'),mart = ensembl)

#########################################
#load human data, 'out' is the data frame
load('asd/rdas/asd_DE_objects_DESeq2_Adj.rda',envir = asd <- new.env())
out = do.call('cbind',lapply(asd$outGeneList,function(g) 
  g[Reduce('intersect',lapply(asd$outGeneList,rownames)),]))
tmp = ss(Reduce('intersect',lapply(asd$outGeneList,rownames)),'\\.')
out = out[!duplicated(tmp),]
rownames(out) = tmp[!duplicated(tmp)]
indFC = grep('log2FoldChange',names(out),value = T)
indPval = grep('pvalue',names(out),value = T)

############################################
# check mouse tcf4 DEGs w/ psych encode data
load("tcf4_mouse/rdas/mega_tcf4_ages_DE_objects_DESeq2.rda") #mouse data
outGeneList = lapply(outGeneList, function(g) {
  g$hsapien_homolog = MMtoHG$hsapiens_homolog_ensembl_gene[match(rownames(g),MMtoHG$ensembl_gene_id)]
  g[grepl('ENSG',g$hsapien_homolog),]})
outGene = outGeneList[['Adult']]
#take only genes w/ human homolog

#####################################################
# get ORs and and p-values for concordance enrichment
for (n in names(outGeneList)){
  g = outGeneList[[n]]; print(n)
  lapply(indPval,function(n){ #enrichment of DEG
    print(n)
    (tt = table(Human = out[g$hsapien_homolog,n] < 0.05,Mouse = g$padj< 0.05))
    print(getOR(tt))
    fisher.test(tt)
  })
  
  # no enforcement mice or human
  g$mouseFCsign = paste0(sign(g$log2FoldChange),'_mm')
  outFCrep = sign(out[g$hsapien_homolog,indFC])
  lapply(colnames(outFCrep),function(i) addmargins(table(outFCrep[,i],g$mouseFCsign)))
  sapply(colnames(outFCrep),function(i) fisher.test(outFCrep[,i],g$mouseFCsign)[1:3])

  # human pval < 0.05
  outFCrep[out[g$hsapien_homolog,indPval] > 0.05] = NA
  lapply(colnames(outFCrep),function(i) addmargins(table(outFCrep[,i],g$mouseFCsign)))
  sapply(colnames(outFCrep),function(i) fisher.test(outFCrep[,i],g$mouseFCsign)[1:3])
  
  # mice padjusted < 0.05
  g = g[g$padj<0.05 & !is.na(g$padj),]
  g$mouseFCsign = paste0(sign(g$log2FoldChange),'_mm')
  outFCrep = sign(out[g$hsapien_homolog,indFC])
  lapply(colnames(outFCrep),function(i) addmargins(table(outFCrep[,i],g$mouseFCsign)))
  sapply(colnames(outFCrep),function(i) fisher.test(outFCrep[,i],g$mouseFCsign)[1:3])
  
  # human pval < 0.05, and mice padjust < 0.05
  outFCrep[out[g$hsapien_homolog,indPval] > 0.05] = NA
  lapply(colnames(outFCrep),function(i) addmargins(table(outFCrep[,i],g$mouseFCsign)))
  sapply(colnames(outFCrep),function(i) fisher.test(outFCrep[,i],g$mouseFCsign)[1:3])
}

##########################################
# save table of replicated human ASD genes
outList = lapply(seq(length(indFC)),function(x){
  ret = outGene[,c('Symbol','log2FoldChange','padj','hsapien_homolog')]
  ret[,indFC[x]] = out[ret$hsapien_homolog,indFC[x]]
  ret[,indPval[x]] = out[ret$hsapien_homolog,indPval[x]]
  ret = ret[complete.cases(ret),]
  ret$sameFCsign = sign(ret[,indFC[x]]) == sign(ret$log2FoldChange)
  ret = ret[ret[,indPval[x]]<.05 & ret$padj < .05,]
  #ret = ret[ ret$padj < .01,]
  ret = ret[order(apply(log(ret[,c('padj',indPval[x])]),1,mean,na.rm = T)),]
  return(ret)
})
names(outList) = names(asd$outGeneList)
type=names(asd$outGeneList)
dat = rbindlist(outList)
names(dat) = c("Symbol", "log2FoldChange_mouse", "padj_mouse", "hsapien_homolog", "log2FoldChange_human", 
               "pvalue_human", "sameFCsign")
dat$type = rep(paste0(type,' (N=',sapply(outList,nrow),')'),times = sapply(outList,nrow))
dat$region = ss(dat$type,'_')
dat$case = ss(dat$type,'_',2)
dat$adj = ss(dat$type,'_',3)
#WriteXLS(dat,ExcelFileName = 'tcf4_mouse/tables/psychENCODE_human_asd_replication.xls')

#####################################
# plot human asd FC v tcf4 mouse
pdf('tcf4_mouse/plots/human_asd_psychENCODE_fc.pdf',w = 12,h = 10)
ggplot(data=subset(dat,case =='ASD'),aes(x=log2FoldChange_mouse, y=log2FoldChange_human, fill=region)) +
  geom_point(alpha = .6,pch=21)+facet_wrap(adj~region,ncol =3,shrink=T) + 
  theme_bw(base_size = 14, base_family = "Helvetica")  + 
  xlab("Mouse log2 fold-change") + ylab('Human log2 fold change') +
  geom_vline(xintercept = 0,colour= 'red',linetype = 'dashed') + 
  geom_hline(yintercept = 0,colour= 'red',linetype = 'dashed')+ 
  theme(strip.background = element_rect(fill = "white", colour = "white"),
        axis.title = element_text(family="Helvetica", face="bold", size=14),
        axis.text = element_text(family="Helvetica", face="bold", size=12),
        legend.position="none")
ggplot(data=subset(dat,case =='Dup15'),aes(x=log2FoldChange_mouse, y=log2FoldChange_human, fill=region)) +
  geom_point(alpha = .6,pch=21)+facet_wrap(adj~region,ncol =3,shrink=T) + 
  theme_bw(base_size = 14, base_family = "Helvetica")  + 
  xlab("Mouse log2 fold-change") + ylab('Human log2 fold change') +
  geom_vline(xintercept = 0,colour= 'red',linetype = 'dashed') + 
  geom_hline(yintercept = 0,colour= 'red',linetype = 'dashed') + 
  theme(strip.background = element_rect(fill = "white", colour = "white"),
        axis.title = element_text(family="Helvetica", face="bold", size=14),
        axis.text = element_text(family="Helvetica", face="bold", size=12),
        legend.position="none")
dev.off()

##########################################################
# fisher exact tests, correlation tests, rho's and kappa's
lapply(split(dat,dat$type),function(x) fisher.test(sign(x$log2FoldChange_mouse),sign(x$log2FoldChange_human)))
lapply(split(dat,dat$type),function(x) with(x,cor.test(log2FoldChange_mouse,log2FoldChange_human,method = 'spearman')))
lapply(split(dat,dat$type),function(x) getKappa(table(sign(x$log2FoldChange_mouse),sign(x$log2FoldChange_human))))

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
