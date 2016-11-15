# make gene set enrichment plots for GO terms

library(lattice)
library(clusterProfiler)
library(RColorBrewer)
library(jaffelab)

load("tcf4_mouse/rdas/gene_sets_tcf4_mouse.rda")
tmp = do.call('rbind',lapply(compareGo,as.data.frame))
tmp$Type = rep(c("MF", "BP",'CC'), sapply(c(compareGo),function(x) nrow(as.data.frame(x))))
tmp$PlotValue = -log10(tmp$p.adjust)
tmp$SetSize = as.numeric(ss(as.character(tmp$BgRatio), "/"))

## filter
tab = table(tmp$ID)
tab = tab[tab>1]
tmp = tmp[(tmp$ID %in% names(tab)),]
tmp$Cluster=  factor(tmp$Cluster, levels = c("Gene","Exon","Junction"))
t = tapply(tmp$PlotValue,tmp$Description,mean)
tmp$Description = factor(tmp$Description, levels = names(t)[order(t,decreasing = F)])
n = max(as.numeric(tmp$Description))
tmp = tmp[as.numeric(tmp$Description) %in% seq(n-20,n),]

## plot
pdf("tcf4_mouse/plots/fig3_GO_terms_heatmap.pdf",w=8)
#range(tmp$PlotValue)
theSeq = seq(0,17,by=0.1) 
my.col <- colorRampPalette(c("white","darkblue"))(length(theSeq))
print(levelplot(PlotValue ~ Cluster + Description, 
                data= tmp, at = theSeq,pretty=TRUE, 
                col.regions = my.col, 
                scales=list(y=list(cex=1.25), x=list(rot=90, cex=1.25)),
                ylab = "", xlab = ""))
dev.off()


