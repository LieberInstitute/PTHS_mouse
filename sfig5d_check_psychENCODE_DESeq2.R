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
ensembl2=useMart("ensembl",dataset = 'hsapiens_gene_ensembl')
MMtoHG = getBM(attributes = c('ensembl_gene_id','hsapiens_homolog_ensembl_gene'),mart = ensembl)
HGEntrez = getBM(attributes = c('ensembl_gene_id','entrezgene'),mart = ensembl2)

#########################################
#load human data, 
load('asd/rdas/asd_DE_objects_DESeq2_Qual.rda',envir = asd <- new.env())
out = do.call('cbind',lapply(asd$outGeneList,function(g) 
  g[Reduce('intersect',lapply(asd$outGeneList,rownames)),]))
tmp = ss(Reduce('intersect',lapply(asd$outGeneList,rownames)),'\\.')
out = out[!duplicated(tmp),]
rownames(out) = tmp[!duplicated(tmp)]
indFC = grep('log2FoldChange',names(out),value = T)
indPval = grep('pval',names(out),value = T)


############################################
# check mouse tcf4 DEGs w/ psych encode data
load("tcf4_mouse/rdas/mega_tcf4_ages_DE_objects_DESeq2.rda") #mouse data
outGeneList = lapply(outGeneList, function(g) {
  g$hsapien_homolog = MMtoHG$hsapiens_homolog_ensembl_gene[match(rownames(g),MMtoHG$ensembl_gene_id)]
  g$hsapien_entrez = HGEntrez$entrezgene[match(g$hsapien_homolog,HGEntrez$ensembl_gene_id)]
  g[grepl('ENSG',g$hsapien_homolog),]})
outGene = outGeneList[['Adult']]
#take only genes w/ human homolog

#####################################################
# get ORs and and p-values for concordance enrichment
for (n in names(outGeneList)){
  g = outGeneList[[n]]; print(n)
  (lapply(indPval,function(n){ #enrichment of DEG
    print(n)
    (tt = table(Human = out[g$hsapien_homolog,n] < 0.05,Mouse = g$padj< 0.05))
    print(getOR(tt))
    fisher.test(tt)
  }))
  
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
  }

##########################################
# save table of replicated human ASD genes
outList = lapply(seq(length(indFC)),function(x){
  ret = outGene[,c('Symbol','log2FoldChange','padj','hsapien_homolog','hsapien_entrez')]
  ret[,indFC[x]] = out[ret$hsapien_homolog,indFC[x]]
  ret[,indPval[x]] = out[ret$hsapien_homolog,indPval[x]]
  ret = ret[complete.cases(ret),]
  ret$sameFCsign = sign(ret[,indFC[x]]) == sign(ret$log2FoldChange)
  ret = ret[ret[,indPval[x]]<.05 & ret$padj < .01,]
  #ret = ret[ ret$padj < .01,]
  ret = ret[order(apply(log(ret[,c('padj',indPval[x])]),1,mean,na.rm = T)),]
  return(ret)
})
names(outList) = names(asd$outGeneList)
type=names(asd$outGeneList)
dat = rbindlist(outList)
names(dat) = c("Symbol", "log2FoldChange_mouse", "padj_mouse", "hsapien_homolog",'hsapien_entrez', "log2FoldChange_human", 
               "pvalue_human", "sameFCsign")
dat$type = rep(paste0(type,' (N=',sapply(outList,nrow),')'),times = sapply(outList,nrow))
dat$region = ss(dat$type,'_')
dat$case = ss(dat$type,'_',2)
dat$adj = ss(dat$type,'_',3)
WriteXLS(dat,ExcelFileName = 'tcf4_mouse/tables/psychENCODE_human_asd_replication.xls')

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

dat = dat[grep('Qual',dat$adj),]
pdf('tcf4_mouse/plots/human_asd_psychENCODE_fc_qual.pdf',w = 10,h = 4)
ggplot(data=dat,aes(x=log2FoldChange_mouse, y=log2FoldChange_human, fill=case)) +
  geom_point(alpha = .6,pch=21)+facet_wrap(~region,ncol =3,shrink=T) + 
  theme_bw(base_size = 14, base_family = "Helvetica")  + 
  ylim(low = -1, high= 1)+
  xlab("Mouse log2 fold-change") + ylab('Human log2 fold change') +
  geom_vline(xintercept = 0,colour= 'red',linetype = 'dashed') + 
  geom_hline(yintercept = 0,colour= 'red',linetype = 'dashed')+ 
  scale_fill_manual(values = c('white','black'))+
  theme(strip.background = element_rect(fill = "white", colour = "white"),
        axis.title = element_text(family="Helvetica", face="bold", size=14),
        axis.text = element_text(family="Helvetica", face="bold", size=12))
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
univ = unique(as.character(HGEntrez$entrezgene[match(rownames(out),HGEntrez$ensembl_gene_id)]))
table(dat$region, dat$case)

############################
# enrich all the GO terms :D
#compareKegg = compareCluster(fun ='enrichKEGG',hsapien_entrez ~ region+case,data = subset(dat,region != 'vermis'),
#                             universe = univ,organism = "hsa", pvalueCutoff = 0.05)
compareGoMf = compareCluster(fun ='enrichGO',hsapien_entrez ~ region+case,data =  subset(dat,region != 'vermis'), 
                             universe = univ, ont = "MF",OrgDb=org.Hs.eg.db, pvalueCutoff = 0.05, readable = TRUE)
compareGoBp = compareCluster(fun ='enrichGO',hsapien_entrez ~ region+case,data =  subset(dat,region != 'vermis'), 
                             universe = univ, ont = "BP",OrgDb=org.Hs.eg.db, pvalueCutoff = 0.05, readable = TRUE)
compareGoCc = compareCluster(fun ='enrichGO',hsapien_entrez ~ region+case,data =  subset(dat,region != 'vermis'), 
                             universe = univ, ont = "CC",OrgDb=org.Hs.eg.db,pvalueCutoff = 0.05, readable = TRUE)
compareGo = lapply(list(compareGoMf=compareGoMf,compareGoBp=compareGoBp,compareGoCc = compareGoCc),dropGO,level = 1:2)

###################
# save the go table
save(compareGo,file = 'tcf4_mouse/rdas/psychENCODE_human_asd_GO_terms.rda')
WriteXLS(lapply(c(compareGo), as.data.frame),ExcelFileName = 'tcf4_mouse/tables/stable8_psychENCODE_human_asd_GO_compareCluster.xls')


############
# make plots
postscript("tcf4_mouse/plots/gene_sets_asd_overlap1.eps",width =5,height = 3)
dotplot(compareGo[[1]],x =~case, colorBy="p.adjust", font.size =10, title = "",showCategory = 3) + 
  ggplot2::facet_grid(~factor(region,levels = c('ba9','ba41-42-22'),labels = c('Frontal','Temporal')))
dev.off()

postscript("tcf4_mouse/plots/gene_sets_asd_overlap2.eps",width =6,height = 3.5)
dotplot(compareGo[[2]],x =~case, colorBy="p.adjust", font.size =10, title = "",showCategory = 3) + 
  ggplot2::facet_grid(~factor(region,levels = c('ba9','ba41-42-22'),labels = c('Frontal','Temporal')))
dev.off()

postscript("tcf4_mouse/plots/gene_sets_asd_overlap3.eps",width =4.5,height = 3.5)
dotplot(compareGo[[3]],x =~case, colorBy="p.adjust", font.size =10, title = "",showCategory = 3) + 
  ggplot2::facet_grid(~factor(region,levels = c('ba9','ba41-42-22'),labels = c('Frontal','Temporal')))
dev.off()


