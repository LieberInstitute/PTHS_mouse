## GO Term analysis for mega analysis pooling Maher, Philpot, and Sweatt mouse RNAseq
# qsub -V -l mf=10G,h_vmem=15G,h_stack=256M -cwd -b y R CMD BATCH gene_sets_mega_dataset.R
getSymbol = function(entrezList, geneMap){
  sapply(entrezList,function(e){
    a = unlist(strsplit(e,'/'))
    s = geneMap$Symbol[match(a,geneMap$EntrezID)]
    return(paste(s,collapse = '/'))
  })
}

library(jaffelab)
library(WriteXLS)
library(clusterProfiler)
library(ReactomePA)
library(org.Mm.eg.db)

load("rdas/mega_tcf4_ages_DE_objects_DESeq2.rda")

######################################
# make genelist of Log2FC, name is EntrezID, order decreasing
GSEList = lapply(outGeneList,function(outGene){
  ord = order(outGene$log2FoldChange,decreasing = T)
  sigStats = outGene$log2FoldChange[ord]
  names(sigStats) = as.character(outGene$EntrezID[ord])
  
  ##################################
  # find GO and KEGG term enrichment
  compareKegg = gseKEGG(sigStats, organism = "mmu", nPerm = 1000, minGSSize    = 120,pvalueCutoff = 0.05)
  compareReactome = gsePathway(sigStats, organism = "mouse", nPerm = 1000, minGSSize    = 120,pvalueCutoff = 0.05)
  #compareMKegg = gseMKEGG(sigStats, organism = "mmu", nPerm = 1000, minGSSize    = 120,pvalueCutoff = 0.05)
  compareGoMf = gseGO(sigStats, ont = "MF",OrgDb=org.Mm.eg.db,  nPerm = 1000, minGSSize    = 120,pvalueCutoff = 0.05)
  compareGoBp = gseGO(sigStats, ont = "BP",OrgDb=org.Mm.eg.db,  nPerm = 1000, minGSSize    = 120,pvalueCutoff = 0.05)
  compareGoCc = gseGO(sigStats, ont = "CC",OrgDb=org.Mm.eg.db,  nPerm = 1000, minGSSize    = 120,pvalueCutoff = 0.05)
  
  ######
  # save
  GOList =lapply(list(Molecular_Function=compareGoMf,
                      Biological_Process=compareGoBp,
                      Cellular_Compartment = compareGoCc,
                      KEGG_Term = compareKegg,
                      Reactome = compareReactome),as.data.frame)
  GO = do.call('rbind',GOList)
  GO$Category = ss(rownames(GO),'\\.')
  GO$core_enrichment = getSymbol(GO$core_enrichment,outGene)
  return(GO)
})
save(GSEList, file = "rdas/gene_sets_mega_tcf4_ages.rda")
WriteXLS(GSEList,ExcelFileName = "tables/gsea_mega_tcf4_ages.xls")


#########################
# define background genes
univ = unlist(sapply(outGeneList,'[[','EntrezID'))
univ = as.character(unique(univ[!is.na(univ)]))

######################################
# make list of genes, exons, junctions
sigStats = do.call('rbind',lapply(outGeneList,function(x) x[x$pvalue < .01,c("EntrezID",'log2FoldChange')]))
sigStats$Age = factor(ss(as.character(rownames(sigStats)),'\\.'),levels = c('p1','Adult'))
sigStats$Dir = factor(ifelse(sigStats$log2FoldChange>0,'Upregulated','Downregulated'))

##################################
# find GO and KEGG term enrichment
compareKegg = compareCluster(fun ='enrichKEGG',EntrezID ~ Age+Dir,data = sigStats, universe = univ,organism = "mmu", pvalueCutoff = 0.05)
compareGoMf = compareCluster(fun ='enrichGO',EntrezID ~ Age+Dir,data = sigStats, universe = univ, ont = "MF",OrgDb=org.Mm.eg.db, pvalueCutoff = 0.05, readable = TRUE)
compareGoBp = compareCluster(fun ='enrichGO',EntrezID ~ Age+Dir,data = sigStats, universe = univ, ont = "BP",OrgDb=org.Mm.eg.db, pvalueCutoff = 0.05, readable = TRUE)
compareGoCc = compareCluster(fun ='enrichGO',EntrezID ~ Age+Dir,data = sigStats, universe = univ, ont = "CC",OrgDb=org.Mm.eg.db,pvalueCutoff = 0.05, readable = TRUE)
compareGo = lapply(list(compareGoMf=compareGoMf,compareGoBp=compareGoBp,compareGoCc = compareGoCc),dropGO,level = 1:2)

######
# save
save(compareKegg, compareGo, file = "rdas/gene_sets_mega_dataset.rda")
WriteXLS(lapply(c(compareGo,compareKegg=compareKegg), as.data.frame),ExcelFileName = "tables/gene_sets_mega_dataset_ages.xls")

compareGo[[1]] = dropGO(compareGo[[1]],term = 'GO:0001227')


############
# make plots
postscript("plots/gene_sets_mega_dataset1.eps",width =5,height = 3)
dotplot(compareGo[[1]],x =~factor(Age,levels = c('p1','Adult')),
        colorBy="p.adjust", font.size =10, title = "",showCategory = 3) + 
  ggplot2::facet_grid(~Dir)
dev.off()

postscript("plots/gene_sets_mega_dataset2.eps",width =6,height = 3.5)
dotplot(compareGo[[2]],x =~factor(Age,levels = c('p1','Adult')), 
        colorBy="p.adjust", font.size =10, title = "",showCategory = 3) + 
  ggplot2::facet_grid(~Dir)
dev.off()

postscript("plots/gene_sets_mega_dataset3.eps",width =6,height = 3.5)
dotplot(compareGo[[3]],x =~factor(Age,levels = c('p1','Adult')),
        colorBy="p.adjust", font.size =10,title = "",showCategory = 3) + 
  ggplot2::facet_grid(~Dir)
dev.off()


