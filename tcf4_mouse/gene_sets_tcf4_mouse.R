## GO Term analysis for Maher tcf4 mouse
# qsub -V -l mf=10G,h_vmem=15G,h_stack=256M -cwd -b y R CMD BATCH gene_sets_tcf4_mouse.R
library(jaffelab)
library(WriteXLS)
library(clusterProfiler)
library(org.Mm.eg.db)

load("rdas/mouse_tcf4_DE_objects_DESeq2.rda")

#########################
# define background genes
univ = outGene$EntrezID
univ = as.character(unique(univ[!is.na(univ)]))

######################################
# make list of genes, exons, junctions
sigStats = outGene[outGene$pvalue < 0.01,c("log2FoldChange", "EntrezID")]
sigStats = sigStats[!is.na(sigStats$EntrezID),]
sigStats$EntrezID = as.character(sigStats$EntrezID)
length(sigStats$EntrezID)

##################################
# find GO and KEGG term enrichment
compareKegg = enrichKEGG(sigStats$EntrezID, universe = univ,organism = "mmu", qvalueCutoff = 0.05,pvalueCutoff = 0.1)
compareGoMf = enrichGO(sigStats$EntrezID, universe = univ, ont = "MF",OrgDb=org.Mm.eg.db, qvalueCutoff = 0.05, pvalueCutoff = 0.1, readable = TRUE)
compareGoBp = enrichGO(sigStats$EntrezID, universe = univ, ont = "BP",OrgDb=org.Mm.eg.db, qvalueCutoff = 0.05, pvalueCutoff = 0.1, readable = TRUE)
compareGoCc = enrichGO(sigStats$EntrezID, universe = univ, ont = "CC",OrgDb=org.Mm.eg.db, qvalueCutoff = 0.05, pvalueCutoff = 0.1, readable = TRUE)
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
save(compareKegg, compareGo, file = "rdas/gene_sets_tcf4_mouse.rda")
WriteXLS(lapply(c(compareGo,compareKegg=compareKegg), as.data.frame),ExcelFileName = "tables/gene_sets_tcf4_mouse.xls")
rm(list = ls())


######################################
# make genelist of Log2FC, name is EntrezID, order decreasing
ord = order(outGene$log2FoldChange,decreasing = T)
sigStats = outGene$log2FoldChange[ord]
names(sigStats) = as.character(outGene$EntrezID[ord])

##################################
# find GO and KEGG term enrichment
compareKegg = gseKEGG(sigStats, organism = "mmu", nPerm = 1000, minGSSize    = 120,pvalueCutoff = 0.05)
compareGoMf = gseGO(sigStats, ont = "MF",OrgDb=org.Mm.eg.db,  nPerm = 1000, minGSSize    = 120,pvalueCutoff = 0.05)
compareGoBp = gseGO(sigStats, ont = "BP",OrgDb=org.Mm.eg.db,  nPerm = 1000, minGSSize    = 120,pvalueCutoff = 0.05)
compareGoCc = gseGO(sigStats, ont = "CC",OrgDb=org.Mm.eg.db,  nPerm = 1000, minGSSize    = 120,pvalueCutoff = 0.05)


######
# save
GOList =lapply(list(Molecular_Function=compareGoMf,Biological_Process=compareGoBp,Cellular_Compartment = compareGoCc,KEGG_Term = compareKegg),as.data.frame)
GO = do.call('rbind',GOList)
GO$Category = ss(rownames(GO),'\\.')
save(GO,compareKegg,compareGoMf,compareGoBp,compareGoCc, file = "rdas/gene_sets_tcf4_mouse.rda")
WriteXLS(GO,ExcelFileName = "tables/gene_sets_tcf4_mouse.xls")

