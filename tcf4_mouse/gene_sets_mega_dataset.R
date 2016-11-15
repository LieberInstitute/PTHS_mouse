## GO Term analysis for mega analysis pooling Maher, Philpot, and Sweatt mouse RNAseq
# qsub -V -l mf=10G,h_vmem=15G,h_stack=256M -cwd -b y R CMD BATCH gene_sets_mega_dataset.R
library(jaffelab)
library(WriteXLS)
library(clusterProfiler)
library(org.Mm.eg.db)

load("rdas/mega_dataset_DE_objects_DESeq2.rda")

#########################
# define background genes
univ = outGene$EntrezID
univ = as.character(unique(univ[!is.na(univ)]))

######################################
# make list of genes, exons, junctions
sigStats = outGene[outGene$padj < .05,c("log2FoldChange", "EntrezID")]
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
pdf("plots/sfig4_gene_sets_mega_dataset.pdf")
dotplot(compareKegg,colorBy="qvalue", font.size =8,title = "KEGG Terms DE Group Comparisons")
dotplot(compareGo[[1]],colorBy="qvalue", font.size =8, title = "Enriched GO-MF Terms")
dotplot(compareGo[[2]],colorBy="qvalue", font.size =8, title = "Enriched GO-BP Terms")
dotplot(compareGo[[3]],colorBy="qvalue", font.size =8, title = "Enriched GO-CC Terms")
dev.off()

######
# save
save(compareKegg, compareGo, file = "rdas/gene_sets_mega_dataset.rda")
WriteXLS(lapply(c(compareGo,compareKegg=compareKegg), as.data.frame),ExcelFileName = "tables/stable6_gene_sets_mega_dataset.xls")

