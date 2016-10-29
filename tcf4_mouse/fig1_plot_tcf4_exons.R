## plot TCF4 Exon fold-change in genomic location order
## plots last exon of ENSMUST00000114982 tcf4 transcript to show timeline expression
cleaningY = function(y, mod, P=ncol(mod)) {
  Hat=solve(t(mod)%*%mod)%*%t(mod)
  beta=(Hat%*%t(y))
  cleany=y-t(as.matrix(mod[,-c(1:P)])%*%beta[-c(1:P),])
  return(cleany)
}

library(jaffelab)
library(limma)
library(biomaRt)
library(ggplot2)

load('/dcl01/lieber/ajaffe/Brady/mouseRNAseq/rpkmCounts_tcf4_mouse_OCT20_n36.rda')
load("rdas/mouse_tcf4_DE_objects.rda")

#################
# filtering exons
eIndex = which(rowMeans(exonRpkm) > 0.1)
exonRpkm = exonRpkm[eIndex,]
exonMap = exonMap[eIndex,]

#########################
# remove outlier mouse 31
pd = pd[-31,]
exonRpkm = exonRpkm[,-31]

#################
# get mouse Tcf4B
ensembl = useMart("ensembl")
ensembl = useDataset("mmusculus_gene_ensembl",mart=ensembl)
sym = getBM(attributes = c('chromosome_name','exon_chrom_start','exon_chrom_end',"ensembl_transcript_id"), 
            mart=ensembl,filters = 'ensembl_transcript_id',value ='ENSMUST00000114982')

#####################
# get only Tcf4 exons
tcf4Index=which(exonMap$Symbol == "Tcf4")
tcf4Map  = exonMap[tcf4Index,]
tcf4Rpkm = exonRpkm[tcf4Index,]

tcf4Index = which(tcf4Map$Start %in% sym$exon_chrom_start& tcf4Map$End %in% sym$exon_chrom_end)
tcf4Map  = tcf4Map[tcf4Index,]
tcf4Rpkm = tcf4Rpkm[tcf4Index,]
meanExonExprs = rowMeans(tcf4Rpkm)

############################
mod = model.matrix(~Genotype+Age, data=pd)

###################
# exon DE, with SVA
load("rdas/mouse_tcf4_sva_objects.rda")
yExon = log2(tcf4Rpkm+1)
modExon = cbind(mod, svaExon$sv)

fitExon = lmFit(yExon, modExon)
outExon = toptable(fitExon,coef = 2,number = nrow(yExon),genelist = cbind(tcf4Map,meanRpkm = meanExonExprs),confint = T)
outExon = outExon[order(outExon$Start),]
#outExon$pos = round(apply(outExon[,c('Start','End')],1,mean))
outExon$pos = factor(outExon$Start)
ind = which(outExon$adj.P.Val<.05)

################################
# Exon plots with standard error
pdf("plots/fig1_tcf4_mouse_exons.pdf",height = 3.5, width = 4)
ylim = c(min(outExon[,c('CI.L','CI.R')])-.2,max(outExon[,c('CI.L','CI.R')])+.2)
plot(outExon$pos,outExon$logFC,xlab = 'chr18',ylab = 'Log2 Fold-change',
     ylim = ylim,main = 'Tcf4 Exon Fold-change', xaxt = 'n',cex.lab = 1.25)
arrows(as.numeric(outExon$pos),outExon$CI.L,as.numeric(outExon$pos),outExon$CI.R,
       code=3,length=0.04,angle=90,col='black')
v1 = seq(1,nrow(outExon),4)
v2 = outExon$pos[v1]
axis(side = 1,at = v1,labels = v2)
abline(h = 0,col = 'red')
text(x = ind,y= outExon$CI.R[ind]+.1,labels = '*',cex = 1.5)

################################
# plot last coded coded exon chr18:69682683-69682823, e643016
tcf4Exon = 2^cleaningY(yExon[rownames(outExon),],modExon, P = 4)-1
pd$ExonExpr = tcf4Exon[nrow(outExon),]
pd$fdr = outExon$adj.P.Val[nrow(outExon)]
ggplot(data=pd) +
  geom_boxplot(aes(x=Age, y=ExonExpr, fill=Genotype),position = position_dodge(width = .5),outlier.shape=NA) +
  theme_bw(base_size = 14, base_family = "Helvetica")  + 
  geom_point(pch=21,aes(x=Age, y=ExonExpr, fill=Genotype), position=position_jitterdodge(dodge.width = .5))+
  scale_fill_manual(values = c("gray","red"),name = "",breaks = c("WT","HT")) + ylim(c(0,max(pd$ExonExpr))) +
  scale_x_discrete(labels =c("P1","P21","Adult")) + 
  ylab("Adjusted FPKM | SVs")+xlab("") + theme(axis.text.x=element_text(angle=45,hjust = 1)) +
  theme(axis.title.y = element_text(family="Helvetica", face="bold", size=16),
        axis.text.y = element_text(family="Helvetica", face="bold", size=14),
        axis.text.x = element_text(family="Helvetica", face="bold", size=14))
dev.off()
