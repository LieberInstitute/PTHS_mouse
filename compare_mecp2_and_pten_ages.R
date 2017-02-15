# compare tcf4 DE t-statistic from mecp2 and pten DE t-statistics
getOR = function(x) x[1,1]/x[2,1]/x[1,2]*x[2,2]
getKappa = function(x) (x[1,1]+x[2,2])/sum(x)

library(jaffelab)
library(reshape2)
library(ggplot2)
library(gplots)
library(clusterProfiler)
library(ReactomePA)
library(org.Mm.eg.db)
library(RColorBrewer)

library("scales")
reverselog_trans <- function(base = exp(1)) {
  trans <- function(x) -log(x, base)
  inv <- function(x) base^(-x)
  trans_new(paste0("reverselog-", format(base)), trans, inv, 
            log_breaks(base = base), 
            domain = c(1e-100, Inf))
}

options(stringsAsFactors = F)

######################################
# load TCF4, MECP2, and PTEN mice data
load('tcf4_mouse/rdas/mega_tcf4_ages_DE_objects_DESeq2.rda',envir = tcf4 <- new.env())
load('mecp2/rdas/mecp2_DE_objects_DESeq2.rda',envir = mecp2 <- new.env())
load('mef2c/rdas/mef2c_DE_objects_DESeq2.rda',envir = mef2c <- new.env())
load('pten/rdas/pten_ages_DE_objects_DESeq2.rda',envir = pten <- new.env())

#########################
# get all expressed genes
geneMap = with(tcf4$outGeneList[['Adult']],tcf4$outGeneList[['Adult']][,7:13])
labs = with(geneMap, rownames(geneMap)[!Symbol %in% c('Mecp2','Tcf4','Pten')])

bigOutGeneList = list(Mecp2 = mecp2$outGene[labs,],
                      Pten = pten$outGeneList[['Adult']][labs,],
                      Tcf4 = tcf4$outGeneList[['Adult']][labs,])


##################################
# Overlap of DEGs in animal models
degList = lapply(bigOutGeneList,function(g) rownames(g)[which(g$padj<.05 & !is.na(g$padj))])
sapply(degList,length)
#Mecp2  Pten  Tcf4 
#252  2815  1830
pdf('tcf4_mouse/plots/overlap_mecp2_pten_tcf4.pdf',w=4,h=4)
geneVenn = plot(venn(degList,small =1.5,show.plot = F),cex = 2)
dev.off()

#########################
# overlap of marginal DEG
sigGenes = Reduce('intersect',lapply(bigOutGeneList,function(g) rownames(g)[which(g$padj<.05 & !is.na(g$padj))]))
entrezID = as.character(geneMap[sigGenes,'EntrezID'])
univ = as.character(geneMap[labs,'EntrezID'])
GoMF = enrichGO(entrezID,ont = "MF", universe = univ, 
                OrgDb=org.Mm.eg.db, pvalueCutoff = 0.05, readable = TRUE)
GoBP = enrichGO(entrezID,ont = "BP", universe = univ, 
                OrgDb=org.Mm.eg.db, pvalueCutoff = 0.05, readable = TRUE)
GoCC = enrichGO(entrezID,ont = "CC", universe = univ, 
                OrgDb=org.Mm.eg.db, pvalueCutoff = 0.05, readable = TRUE)
GoBP = dropGO(GoBP,level = 1:2)
GoCC = dropGO(GoCC,level = 1:2)
save(GoBP,GoCC, file = 'tcf4_mouse/rdas/overlap_mecp2_tcf4_pten_genesets.rda')

dat = do.call('rbind',lapply(list(Biological_Process = GoBP,
                                  Cellular_Component = GoCC),
                             as.data.frame))
dat$Type = ss(as.character(rownames(dat)),'\\.')
dat$GeneRatio = as.numeric(ss(dat$GeneRatio,'/'))/as.numeric(ss(dat$GeneRatio,'/',2))
dat$Description = factor(dat$Description,
           levels = dat$Description[order(dat$p.adjust,decreasing = T)])
WriteXLS::WriteXLS(list(GO = dat,Genes = geneMap[sigGenes,]),'tcf4_mouse/tables/overlap_mecp2_pten_tcf4.xls')
#dat = subset(dat,p.adjust < 0.01)
pdf('tcf4_mouse/plots/overlap_gene_sets.pdf',height = 2.5,width = 5)
ggplot(data = dat, aes(fill = GeneRatio,x = Description,
        y = p.adjust)) + geom_bar(stat="identity")+
  coord_flip()+ylab('Adjusted P-value')+xlab('Gene Ontology')+
  scale_y_continuous(trans=reverselog_trans(10))+
  geom_hline(yintercept = 0.05,color = 'red',linetype = 2)+
  scale_fill_gradient(low ="lightgray", high = "black")
dev.off()


#######################################
# compare TCF4 vs. MECP2 log2FC overlap
tmPval = do.call('cbind',lapply(bigOutGeneList[c('Tcf4','Mecp2')],'[[','padj'))
tmLog2FC = do.call('cbind',lapply(bigOutGeneList[c('Tcf4','Mecp2')],'[[','log2FoldChange'))
rownames(tmLog2FC) = rownames(bigOutGeneList[['Tcf4']])
tmLog2FC = tmLog2FC[apply(tmPval<0.05 &!is.na(tmPval),1,all),]

# test for negative enrichment = 
(t1 = with(data.frame(tmPval),table(Tcf4 = Tcf4<0.05 & !is.na(Tcf4),Mecp2 = Mecp2<0.05 & !is.na(Mecp2))))
fisher.test(t1)
getOR(t1)

(t2 = with(data.frame(tmLog2FC),table(Tcf4 = sign(Tcf4) ,Mecp2 = sign(Mecp2))))
1-getKappa(t2)
with(data.frame(tmLog2FC),cor.test(Tcf4,Mecp2,method = 'spearman'))
with(data.frame(tmLog2FC),cor.test(x = Tcf4,y = Mecp2,method = 'spearman'))
plot(tmLog2FC[,2],tmLog2FC[,1],ylab = 'Tcf4 log2 Fold-change',xlab = 'Mecp2 log2 Fold-change',
     xlim = c(-.6,.6),ylim = c(-.5,.5))

##############################################
# compare PTEN vs. TCF4 logFC (anti) enrichment
tpPval = do.call('cbind',lapply(bigOutGeneList[c('Tcf4','Pten')],'[[','padj'))
tpLog2FC = do.call('cbind',lapply(bigOutGeneList[c('Tcf4','Pten')],'[[','log2FoldChange'))
rownames(tpLog2FC) = rownames(bigOutGeneList[['Tcf4']])
tpLog2FC = tpLog2FC[apply(tpPval<0.05 &!is.na(tpPval),1,all),]
  
# test for negative enrichment = 
(t1 = with(data.frame(tpPval),table(Tcf4 = Tcf4<0.05 & !is.na(Tcf4),
                                    Pten = Pten<0.05 & !is.na(Pten))))
fisher.test(t1)
getOR(t1)
(t2 = with(data.frame(tpLog2FC),table(Tcf4 = sign(Tcf4) ,Pten = sign(Pten))))
1-getKappa(t2)
with(data.frame(tpLog2FC),cor.test(Tcf4,Pten,method = 'spearman'))

with(data.frame(tpLog2FC),cor.test(x = Tcf4,y = Pten,method = 'spearman'))
plot(tpLog2FC[,2],tpLog2FC[,1],ylab = 'Tcf4 log2 Fold-change',xlab = 'Pten log2 Fold-change',
     xlim = c(-1.2,1.2),ylim = c(-.6,.6))

#####################################
# plot log2fc TCF4 vs. PTEN and MECP2
sig = Reduce('intersect',lapply(bigOutGeneList,function(g) rownames(g)[which(g$padj<0.05 & !is.na(g$padj))]))
apply(tmLog2FC,2,range)
apply(tpLog2FC,2,range)

postscript('tcf4_mouse/plots/asd_mice_compare_logFC.eps',height = 3.5,width = 4.5)
ggplot(data = data.frame(tmLog2FC),aes(y = Tcf4,x = Mecp2,colour =  rownames(tmLog2FC) %in% sig))+
  geom_point(fill = '#00cd00',pch = 21) + ggtitle('')+ 
  ylab('Tcf4 log2 Fold-Change') + xlab('Mecp2 log2 Fold-Change')+
  geom_vline(xintercept = 0,colour= 'red',linetype = 'dashed') + 
  geom_hline(yintercept = 0,colour= 'red',linetype = 'dashed') +
  scale_color_manual(values = c('FALSE'='#00cd00','TRUE'='black'))+
  theme(strip.background = element_rect(fill = "white", colour = "white"),
        legend.position="none")
ggplot(data = data.frame(tpLog2FC),aes(y = Tcf4,x = Pten,colour =  rownames(tpLog2FC) %in% sig))+
  geom_point(fill = '#0080ff',pch =21) + ggtitle('')+ 
  ylab('Tcf4 log2 Fold-Change') + xlab('Pten log2 Fold-Change')+
  geom_vline(xintercept = 0,colour= 'red',linetype = 'dashed') + 
  geom_hline(yintercept = 0,colour= 'red',linetype = 'dashed') +
  scale_color_manual(values = c('FALSE'='#0080ff','TRUE'='black'))+
  theme(strip.background = element_rect(fill = "white", colour = "white"),
        legend.position="none")
dev.off()