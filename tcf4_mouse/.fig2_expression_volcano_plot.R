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

load('/dcl01/lieber/ajaffe/Brady/mouseRNAseq/rpkmCounts_tcf4_mouse_OCT20_n36.rda')
load("rdas/mouse_tcf4_DE_objects.rda")
load("rdas/mouse_tcf4_sva_objects.rda")

######################
# throw out mouse 31
geneRpkm = geneRpkm[,-31]
yGene = log2(geneRpkm+1)

#######################
# gene expression plots
mod = model.matrix(~Genotype+Age,data = pd)
modGene = cbind(mod, svaGene$sv)
sigGene = outGene[outGene$fdr < 0.1,]

geneClean = 2^cleaningY(yGene[rownames(sigGene),], modGene, P=4)-1
dat = cbind(pd,t(geneClean))
names(dat) = c(names(pd),as.character(sigGene$Symbol[-7]))
dat.melt = melt(data = dat, id = names(pd),value.name = "FPKM")
dat.melt$variable = factor(dat.melt$variable)

par(mar=c(2,2,5,1), cex.lab=1.5, font.main = 2,cex.main = 2,cex.axis = 1.2,font.lab=2)
ggplot(data=dat.melt,aes(x=Age, y=FPKM, fill=Genotype)) +
  geom_boxplot(aes(x=Age, y=FPKM, fill=Genotype), position = position_dodge(width = .5),outlier.shape=NA) +
  theme_bw(base_size = 14, base_family = "Helvetica")  + facet_wrap(~variable,scales ="free_y",nrow =3) + 
  geom_point(pch=21, position=position_jitterdodge(dodge.width = .5)) +
  scale_fill_manual(values=c("gray","red"),breaks= c("WT","Het"),labels = c("","")) +
  scale_x_discrete(labels =c("P1","P21","Adult")) + ylab("Adjusted FPKM | SVs") +
  xlab("") + theme(axis.text.x=element_text(angle=45,hjust = 1),
                   strip.background = element_rect(fill = "white", colour = "white"),
                   text = element_text(family="Helvetica", face="bold", size=16),
                   legend.position = "none")

###############
# volcano plots


