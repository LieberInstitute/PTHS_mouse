enrichCellType = function(genes,univ){
  load('/users/bphan/tcf4/PTHS_mouse/zhang/rdas/mouse_brain_celltypes_geneset_DESeq2.rda',
       envir = cell <-new.env())
  fishtest = data.frame(t(sapply(cell$cellSet,function(x)
    fisher.test(table(Query = univ %in% genes, Set = univ %in% x),alternative = 'g'))))
  fishtest = fishtest[,c('p.value','estimate')]
  fishtest$estimate = unlist(fishtest$estimate)
  fishtest$padj = p.adjust(fishtest$p.value,'fdr')
  fishtest$enriched = with(fishtest, padj< 0.05)
  fishtest$percent = sapply(cell$cellSet,function(x) sum(genes %in% x)/length(x))
  return(fishtest)
}

library(jaffelab)
load('rdas/mega_tcf4_ages_DE_objects_DESeq2.rda', envir = mega <- new.env())
names(mega$outGeneList) = paste0(rep('Pooled.'),names(mega$outGeneList))

load('rdas/mega_tcf4_separate_DE_objects_DESeq2.rda')
outGeneList = c(outGeneList,mega$outGeneList)
univ = unique(unlist(lapply(outGeneList,rownames)))
datList = lapply(lapply(outGeneList,function(g) rownames(g)[g$padj<0.05&!is.na(g$padj)]),
       enrichCellType,univ = univ)

##################
# plot enrichment
library(ggplot2)
dat = do.call('rbind',datList)
dat$percent[!dat$enriched] = 0
dat$Line = factor(ss(rownames(dat),'\\.'),levels = 
                    c('R579W','Del','Nest','Act','Sweatt','Maher','Pooled'))
dat$Age = factor(ss(rownames(dat),'\\.',2),levels = c('p1','Adult'))
dat$Celltype = factor(ss(rownames(dat),'\\.',3),levels = 
                        rev(c('Astrocyte','Endothelial','Microglia','Neuron',
                        'Oligo_Pre','Oligo_New','Oligo_Mye')))
dat$fdrSignif = ifelse(dat$padj < 0.001, "***",
                ifelse(dat$padj < 0.01, "**",
                ifelse(dat$padj < 0.05, "*",'')))

pdf('plots/cellTypes_mega_tcf4_heatmap.pdf',height = 6,width = 8)
ggplot(data = dat,aes(fill=percent, y = Celltype, x = Line)) +
  geom_tile(colour = 'black') +xlab('Mouse Model') + ylab('Cell Type') + 
  geom_text(colour = 'black', aes(label = fdrSignif)) +
  facet_wrap(~Age,scales = 'free')+ scale_fill_gradient(low="white", high="black",
           guide = guide_legend(title = "Gene Set Ratio"))+
  theme(strip.background = element_rect(fill = "white", colour = "white"),
        axis.text.x = element_text(angle = 30, hjust = 1))
dev.off()

