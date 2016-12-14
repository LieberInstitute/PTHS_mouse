# compare DEGs in parkinson mouse model w/ TCF4 mouse models
getOR = function(x) x[1,1]/x[2,1]/x[1,2]*x[2,2]
getKappa = function(x) (x[1,1]+x[2,2])/sum(x)

library(reshape2)
library(ggplot2)
options(stringsAsFactors = F)

######################################
# load TCF4, MECP2, and PTEN mice data
load('../tcf4_mouse/rdas/mega_dataset_DE_objects_DESeq2.rda',envir = tcf4 <- new.env())
load('../mecp2/rdas/mecp2_DE_objects_DESeq2.rda',envir = mecp2 <- new.env())
load('../pten/rdas/pten_DE_objects_DESeq2.rda',envir = pten <- new.env())
load('../parkinson/rdas/parkinson_DE_objects_DESeq2.rda',envir = parkinson <- new.env())

# gene names
geneMap = tcf4$outGene[,7:13]
labs = rownames(geneMap)

##################################
# PCA of ASDs and parkinson log2FC
col = 2
outGene = data.frame(TCF4= tcf4$outGene[labs,col],PTEN = pten$outGene[labs,col],
                     MECP2 =mecp2$outGene[labs,col],PARK = parkinson$outGene[labs,col])
#outGene = -log10(outGene[complete.cases(outGene),])
outGene = outGene[complete.cases(outGene),]
plot(hclust(dist(t(outGene))))

pca = prcomp(t(outGene))
plot(pca$x,pch =21, bg = 1:4)


ggplot(data = outGene1,aes(x = TCF4_log2FC,y = log2FoldChange))+
  geom_point(alpha = .5) +
  ggtitle('')+ xlab('Tcf4 log2 fold change') + ylab('Mouse Mutation log2 fold change')+
  geom_vline(aes(xintercept = 0,colour= 'red',linetype = 'dashed')) + 
  geom_hline(aes(yintercept = 0,colour= 'red',linetype = 'dashed'))+
  theme(strip.background = element_rect(fill = "white", colour = "white"),
        axis.title = element_text(family="Helvetica", face="bold", size=14),
        axis.text = element_text(family="Helvetica", face="bold", size=12),
        legend.position="none")

cat(outGene1$Symbol, file = '../parkinson/tables/parkinsons.txt')