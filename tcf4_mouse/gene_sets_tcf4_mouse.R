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
pdf("plots/gene_sets_tcf4_mouse.pdf")
dotplot(compareKegg,colorBy="qvalue", font.size =8,title = "KEGG Terms DE Group Comparisons")
dotplot(compareGo[[1]],colorBy="qvalue", font.size =8, title = "Enriched GO-MF Terms")
dotplot(compareGo[[2]],colorBy="qvalue", font.size =8, title = "Enriched GO-BP Terms")
dotplot(compareGo[[3]],colorBy="qvalue", font.size =8, title = "Enriched GO-CC Terms")
dev.off()

######
# save
save(compareKegg, compareGo, file = "rdas/gene_sets_tcf4_mouse.rda")
WriteXLS(lapply(c(compareGo,compareKegg=compareKegg), as.data.frame),ExcelFileName = "tables/gene_sets_tcf4_mouse.xls")
rm(list = ls())

