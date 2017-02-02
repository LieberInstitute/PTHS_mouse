## GO Term analysis for Maher tcf4 mouse
# qsub -V -l mf=10G,h_vmem=15G,h_stack=256M -cwd -b y R CMD BATCH gene_sets_tcf4_mouse.R
getSymbol = function(entrezList, geneMap){
  sapply(entrezList,function(e){
    a = unlist(strsplit(e,'/'))
    s = geneMap$Symbol[match(a,geneMap$EntrezID)]
    return(paste(s,collapse = '/'))
  })
}
library(jaffelab)
library(ReactomePA)
library(WriteXLS)
library(clusterProfiler)
library(org.Mm.eg.db)

#load("rdas/mouse_tcf4_DE_objects_DESeq2.rda")
load('rdas/mouse_tcf4_ages_DE_objects_DESeq2.rda')

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
save(GSEList, file = "rdas/gene_sets_tcf4_mouse_ages.rda")
WriteXLS(GSEList,ExcelFileName = "tables/gene_sets_tcf4_mouse_ages.xls")























#########################
# define background genes
#univ = outGene$EntrezID
univ = unlist(lapply(outGeneList,'[[','EntrezID'))
univ = as.character(unique(univ[!is.na(univ)]))

######################################
# make list of genes, exons, junctions
#sigStats = outGene[outGene$pvalue < 0.01,c("log2FoldChange", "EntrezID")]
#sigStats = sigStats[!is.na(gList),]
#gList = as.character(gList)
#length(gList)
sigList = lapply(outGeneList,function(outGene) outGene[outGene$pvalue < 0.01,c("log2FoldChange", "EntrezID")])
gList = endoapply(lapply(sigList,'[[','EntrezID'),function(x) as.character(x[!is.na(x)]))
sapply(gList, length)

##################################
# find GO and KEGG term enrichment
#compareKegg = enrichKEGG(gList, universe = univ,organism = "mmu", qvalueCutoff = 0.05,pvalueCutoff = 0.1)
#compareGoMf = enrichGO(gList, universe = univ, ont = "MF",OrgDb=org.Mm.eg.db, qvalueCutoff = 0.05, pvalueCutoff = 0.1, readable = TRUE)
#compareGoBp = enrichGO(gList, universe = univ, ont = "BP",OrgDb=org.Mm.eg.db, qvalueCutoff = 0.05, pvalueCutoff = 0.1, readable = TRUE)
#compareGoCc = enrichGO(gList, universe = univ, ont = "CC",OrgDb=org.Mm.eg.db, qvalueCutoff = 0.05, pvalueCutoff = 0.1, readable = TRUE)
#compareGo = lapply(list(compareGoMf=compareGoMf,compareGoBp=compareGoBp,compareGoCc = compareGoCc),simplify)

compareKegg = compareCluster(gList, fun ='enrichKEGG', universe = univ,organism = "mmu", qvalueCutoff = 0.05,pvalueCutoff = 0.1)
compareGoMf = compareCluster(gList, fun ='enrichGO', universe = univ, ont = "MF",OrgDb=org.Mm.eg.db, qvalueCutoff = 0.05, pvalueCutoff = 0.1, readable = TRUE)
compareGoBp = compareCluster(gList, fun ='enrichGO', universe = univ, ont = "BP",OrgDb=org.Mm.eg.db, qvalueCutoff = 0.05, pvalueCutoff = 0.1, readable = TRUE)
compareGoCc =compareCluster(gList, fun ='enrichGO', universe = univ, ont = "CC",OrgDb=org.Mm.eg.db, qvalueCutoff = 0.05, pvalueCutoff = 0.1, readable = TRUE)
compareGo = lapply(list(compareGoMf=compareGoMf,compareGoBp=compareGoBp,compareGoCc = compareGoCc),simplify)
############
# make plots
#pdf("plots/gene_sets_tcf4_mouse.pdf")
#dotplot(compareKegg,colorBy="qvalue", font.size =8,title = "KEGG Terms DE Group Comparisons")
#dotplot(compareGo[[1]],colorBy="qvalue", font.size =8, title = "Enriched GO-MF Terms")
#dotplot(compareGo[[2]],colorBy="qvalue", font.size =8, title = "Enriched GO-BP Terms")
#dotplot(compareGo[[3]],colorBy="qvalue", font.size =8, title = "Enriched GO-CC Terms")
#dev.off()

######
# save
save(compareKegg, compareGo, file = "rdas/gene_sets_tcf4_mouse_ages.rda")
WriteXLS(lapply(c(compareGo,compareKegg=compareKegg), as.data.frame),ExcelFileName = "tables/gene_sets_tcf4_mouse_ages.xls")

