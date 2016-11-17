## make DE volcano plot and top DE gene expression plots over development
cleaningY = function(y, mod, P=ncol(mod)) {
  Hat=solve(t(mod)%*%mod)%*%t(mod)
  beta=(Hat%*%t(y))
  cleany=y-t(as.matrix(mod[,-c(1:P)])%*%beta[-c(1:P),])
  return(cleany)
}

library(jaffelab)
library(ggplot2)
library(reshape2)
library(DESeq2)

load( '/dcl01/lieber/ajaffe/Brady/mouseRNAseq/mouse_tcf4_DESeq2_svaAdj.rda')
load("tcf4_mouse/rdas/mouse_tcf4_DE_objects_DESeq2.rda")

####################################
# gene expression summary statistics
gIndex = which(outGene$pvalue < 0.01)
sigGene = outGene[gIndex,]
sum(sigGene$pvalue < 0.01, na.rm = T) - sum(sigGene$padj < 0.05, na.rm = T) #405
table(sign(sigGene$log2FoldChange))

######################
# get rlog counts
yGene = rlog(geneDds)
yGene = assays(yGene)[[1]]
plot(hclust(dist(t(yGene))))

####################
# extract model matrix
pd = as.data.frame(colData(geneDds))
modGene = model.matrix(design(geneDds),data = pd)

#######################
# gene expression plots
sigGene = outGene[which(outGene$padj < 0.05 & !is.na(outGene$padj )),]
sigGene = sigGene[-which(sigGene$Symbol =='Tcf4'),]
sigGene$Symbol = factor(sigGene$Symbol)

geneClean = 2^cleaningY(yGene[rownames(sigGene)[1:6],], modGene, P=4)-1
dat = cbind(pd,t(geneClean))
names(dat) = c(names(pd),as.character(sigGene$Symbol[1:6]))
dat.melt = melt(data = dat, id = names(pd),value.name = "FPKM",stringsAsFactors = FALSE)
dat.melt$variable = factor(dat.melt$variable,levels = levels(sigGene$Symbol))

pdf('tcf4_mouse/plots/fig2_expression_plots.pdf',width=5, height=6)
par(mar=c(2,2,5,1), cex.lab=1.5, font.main = 2,cex.main = 2,cex.axis = 1.2,font.lab=2)
ggplot(data=dat.melt,aes(x=Age, y=FPKM, fill=Genotype)) +
  geom_boxplot(aes(x=Age, y=FPKM, fill=Genotype), position = position_dodge(width = .5),outlier.shape=NA) +
  theme_bw(base_size = 14, base_family = "Helvetica")  + facet_wrap(~variable,scales ="free_y",nrow =3) + 
  geom_point(pch=21, position=position_jitterdodge(dodge.width = .5)) +
  scale_fill_manual(values=c("gray","red"),breaks= c("WT","Het"),labels = c("","")) +
  scale_x_discrete(labels =c("P1","P21","Adult")) + ylab("Adjusted Counts | SVs") +
  xlab("") + theme(strip.background = element_rect(fill = "white", colour = "white"),
                   text = element_text(family="Helvetica", face="bold", size=16),
                   legend.position = "none")

dev.off()


############################
# plot supplemental figures
geneClean = 2^cleaningY(yGene[rownames(sigGene)[-c(1:6)],], modGene, P=4)-1
dat = cbind(pd,t(geneClean))
names(dat) = c(names(pd),as.character(sigGene$Symbol[-c(1:6)]))
dat.melt = melt(data = dat, id = names(pd),value.name = "FPKM")
dat.melt$variable = factor(as.character(dat.melt$variable))

pdf('tcf4_mouse/plots/sfig2_expression_plots.pdf',width = 8.5, height = 11)
ggplot(data=dat.melt[as.numeric(dat.melt$variable) %in% c(1:18),],aes(x=Age, y=FPKM, fill=Genotype)) +
  geom_boxplot(aes(x=Age, y=FPKM, fill=Genotype), position = position_dodge(width = .5),outlier.shape=NA) +
  theme_bw(base_size = 14, base_family = "Helvetica")  + facet_wrap(~variable,scales ="free_y",ncol =3) + 
  geom_point(pch=21, position=position_jitterdodge(dodge.width = .5)) +
  scale_fill_manual(values=c("gray","red"),breaks= c("WT","Het"),labels = c("","")) +
  scale_x_discrete(labels =c("P1","P21","Adult")) + ylab("Adjusted Counts | SVs") +
  xlab("") + theme(axis.text.x=element_text(angle=45,hjust = 1),
                   strip.background = element_rect(fill = "white", colour = "white"),
                   text = element_text(family="Helvetica", size=12),
                   legend.position = "none")

ggplot(data=dat.melt[as.numeric(dat.melt$variable) %in% c(19:36),],aes(x=Age, y=FPKM, fill=Genotype)) +
  geom_boxplot(aes(x=Age, y=FPKM, fill=Genotype), position = position_dodge(width = .5),outlier.shape=NA) +
  theme_bw(base_size = 14, base_family = "Helvetica")  + facet_wrap(~variable,scales ="free_y",ncol =3) + 
  geom_point(pch=21, position=position_jitterdodge(dodge.width = .5)) +
  scale_fill_manual(values=c("gray","red"),breaks= c("WT","Het"),labels = c("","")) +
  scale_x_discrete(labels =c("P1","P21","Adult")) + ylab("Adjusted Counts | SVs") +
  xlab("") + theme(axis.text.x=element_text(angle=45,hjust = 1),
                   strip.background = element_rect(fill = "white", colour = "white"),
                   text = element_text(family="Helvetica", size=12),
                   legend.position = "none")
dev.off()

###############
# volcano plots
outGene$col = '#000000' #black
outGene$col[outGene$pvalue<0.01] = '#ff69b480' #pink
outGene$col[outGene$padj<0.05] = '#ff0000' #red

pdf("tcf4_mouse/plots/fig2_DE_gene_MA_plot.pdf",width=5, height=6)
par(mar = c(5,5,2,2),cex.lab = 1.5,cex.axis = 1.2,font.lab=2)
with(subset(outGene,pvalue < 0.01),
     plot(baseMean,log2FoldChange,pch = 20, col = col,
          xlab = 'Normalized Mean Counts', ylab = 'Log2 Fold-change', 
          log = 'x', ylim = c(-.6, .6)))
with(subset(outGene,pvalue > 0.01),
     points(baseMean,log2FoldChange,pch = '.',col = '#00000050',cex = 2))
abline(h = 0, lty = 1, col = 'red')
legend('topright', legend = c('FDR < 0.05','pvalue < 0.01'), pch = 21, pt.bg = c('#ff0000','#ff69b480'))
dev.off()
