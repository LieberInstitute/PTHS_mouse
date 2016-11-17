## plot TCF4 Exon fold-change in genomic location order
## plots last exon of ENSMUST00000114985 tcf4 transcript to show timeline expression
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
library(DESeq2)

load( '/dcl01/lieber/ajaffe/Brady/mouseRNAseq/mouse_tcf4_DESeq2_svaAdj.rda')
load("tcf4_mouse/rdas/mouse_tcf4_DE_objects_DESeq2.rda")

######################################
# summary of TCF4 expression patterns
outGene[outGene$Symbol=='Tcf4',]
tcf4Index = which(outExon$Symbol=='Tcf4' & 
                    outExon$padj< 0.05 & 
                    !is.na(outExon$padj) &
                    outExon$log2FoldChange <0)
apply(outExon[tcf4Index,1:6],2,summary)

#################
# summary of DEGs
gIndex = which(outGene$Symbol!='Tcf4' & 
                    outGene$padj< 0.05 & 
                    !is.na(outGene$padj))
sigGene = outGene[gIndex,]
dim(sigGene) #42 DEGs aside from TCF4
table(sign(sigGene$log2FoldChange)) # 21 up, 21 down

#################
# get rlog counts
yExon = rlog(exonDds)
yExon = assays(yExon)[[1]]

#################
# get mouse Tcf4B
ensembl = useMart("ensembl")
ensembl = useDataset("mmusculus_gene_ensembl",mart=ensembl)
sym = getBM(attributes = c('chromosome_name','exon_chrom_start','exon_chrom_end',"ensembl_transcript_id"), 
            mart=ensembl,filters = 'ensembl_transcript_id',value ='ENSMUST00000114985')

#####################
# get only Tcf4 exons
tcf4Index = rownames(outExon)[which(outExon$Symbol == "Tcf4")]
tcf4Map = outExon[tcf4Index,]
tcf4Counts = yExon[tcf4Index,]

tcf4Index = which(tcf4Map$Start %in% sym$exon_chrom_start & 
                    tcf4Map$End %in% sym$exon_chrom_end)
tcf4Map  = tcf4Map[tcf4Index,]
tcf4Counts = tcf4Counts[tcf4Index,]
meanExonExprs = rowMeans(tcf4Counts)

############################
pd = as.data.frame(colData(exonDds))
modExon = model.matrix(design(exonDds),data = pd)

###################
# exon DE, with SVA
outExon2 = tcf4Map
outExon2$CI.L = with(outExon2,log2FoldChange -lfcSE)
outExon2$CI.R = with(outExon2,log2FoldChange +lfcSE)
outExon2 = outExon2[order(outExon2$Start),]
#outExon$pos = round(apply(outExon[,c('Start','End')],1,mean))
outExon2$pos = as.numeric(factor(outExon2$Start))
ind = which(outExon2$padj<.05)

#############################################################
# Exon plots with standard error, mice TCF4 on forward strand
pdf("tcf4_mouse/plots/sfig1_tcf4_transcript_2_exons.pdf",height = 8, width = 10)
ylim = c(min(outExon2[,c('CI.L','CI.R')])-.1,max(outExon2[,c('CI.L','CI.R')])+.1)
outExon2$col = ifelse(outExon2$padj<.05,'red' ,'black')
par(cex = 1.5)
plot(outExon2$pos,outExon2$log2FoldChange,xlab = 'Exon #',ylab = 'Log2 Fold-change',
     ylim = ylim,main = 'Tcf4 Transcript 1 Exon Fold-change',cex.lab = 1.25,bg = outExon2$col,pch = 21)
arrows(as.numeric(outExon2$pos),outExon2$CI.L,as.numeric(outExon2$pos),outExon2$CI.R,
       code=3,length=0.04,angle=90,col='black')
v1 = seq(1,nrow(outExon2),4)
v2 = outExon2$pos[v1]
abline(h = 0,col = 'red')
text(x = ind,y= outExon2$CI.R[ind]+.1,labels = '*',cex = 1.5)
legend('bottomleft', legend = '* FDR < 1e-10',bty = "n")
dev.off()

################################
# plot last coded coded exon chr18:69683077-69684528, e672980
tcf4Exon = 2^cleaningY(yExon[rownames(outExon),],modExon, P = 4)-1
pd$ExonExpr = tcf4Exon['e673537',]
pd$fdr = outExon['e673537',c('padj')]

pdf("tcf4_mouse/plots/fig1_tcf4_last_exon.pdf",height = 3.5, width = 3.5)
ggplot(data=pd) +
  geom_boxplot(aes(x=Age, y=ExonExpr, fill=Genotype),position = position_dodge(width = .5),outlier.shape=NA) +
  theme_bw(base_size = 14, base_family = "Helvetica")  + 
  geom_point(pch=21,aes(x=Age, y=ExonExpr, fill=Genotype), position=position_jitterdodge(dodge.width = .5))+
  scale_fill_manual(values = c("gray","red"),name = "",breaks = c("WT","HT")) +
  scale_x_discrete(labels =c("P1","P21","Adult")) + 
  ylab("Adjusted Counts | SVs")+xlab("") + 
  theme(axis.title.y = element_text(family="Helvetica", face="bold", size=14),
        axis.text.y = element_text(family="Helvetica", face="bold", size=12),
        axis.text.x = element_text(family="Helvetica", face="bold", size=12),
        legend.position="none")
dev.off()
