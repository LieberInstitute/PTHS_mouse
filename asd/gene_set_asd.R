# gene ontology analysis
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
library(org.Hs.eg.db)

load('rdas/asd_DE_objects_DESeq2_Adj.rda')
outGeneList = outGeneList[grep('Qual',names(outGeneList))]

#########################
# define background genes
univ = unlist(sapply(outGeneList,'[[','EntrezID'))
univ = as.character(unique(univ[!is.na(univ)]))

######################################
# make list of genes, exons, junctions
sigStats = do.call('rbind',lapply(outGeneList,function(x) x[x$pvalue < .05,c("EntrezID",'log2FoldChange')]))
sigStats$Diagnosis = factor(ss(as.character(rownames(sigStats)),'_',2),levels = c('ASD','Dup15'))
sigStats$Region = factor(ss(as.character(rownames(sigStats)),'_'),levels = c('ba41-42-22','ba9','vermis'))

##################################
# find GO and KEGG term enrichment
# compareKegg = compareCluster(fun ='enrichKEGG',EntrezID ~ Diagnosis+Region,data = sigStats, universe = univ,organism = "mmu", pvalueCutoff = 0.05)
compareGoMf = compareCluster(fun ='enrichGO',EntrezID ~ Diagnosis+Region,data = sigStats, universe = univ, ont = "MF",OrgDb=org.Hs.eg.db, pvalueCutoff = 0.05, readable = TRUE)
compareGoBp = compareCluster(fun ='enrichGO',EntrezID ~ Diagnosis+Region,data = sigStats, universe = univ, ont = "BP",OrgDb=org.Hs.eg.db, pvalueCutoff = 0.05, readable = TRUE)
compareGoCc = compareCluster(fun ='enrichGO',EntrezID ~ Diagnosis+Region,data = sigStats, universe = univ, ont = "CC",OrgDb=org.Hs.eg.db,pvalueCutoff = 0.05, readable = TRUE)
compareGo = lapply(list(compareGoMf=compareGoMf,compareGoBp=compareGoBp,compareGoCc = compareGoCc),dropGO,level = 1:2 )

######
# save
save(compareGo, file = "rdas/gene_sets_asd.rda")
WriteXLS(lapply(compareGo, as.data.frame),ExcelFileName = "tables/gene_sets_asd.xls")
