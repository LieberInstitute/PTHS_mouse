## GO Term analysis for Maher tcf4 mouse
# qsub -V -l mf=10G,h_vmem=15G,h_stack=256M -cwd -b y R CMD BATCH gene_sets_tcf4_mouse.R
library(jaffelab)
library(WriteXLS)
library(clusterProfiler)
library(org.Mm.eg.db)

load("rdas/mouse_tcf4_DE_objects_DESeq2.rda")

#########################
# define background genes
outJxn$EntrezID = outGene[outJxn$newGeneID,c('EntrezID')]
univ = c(outGene$EntrezID,outExon$EntrezID, outJxn$EntrezID)
univ = as.character(unique(univ[!is.na(univ)]))

######################################
# make list of genes, exons, junctions
statList = list(Gene = outGene, Exon = outExon, Junction= outJxn)
sigStatList = lapply(statList, function(x) x[x$pvalue < 0.01,c("log2FoldChange", "EntrezID")])
sigStats = do.call("rbind", sigStatList)
sigStats$type = ss(rownames(sigStats), "\\.")
sigStats = sigStats[!is.na(sigStats$EntrezID),]
sigStats = sigStats[!duplicated(sigStats[,c("EntrezID", "type")]),]
sigStats$EntrezID = as.character(sigStats$EntrezID)
tapply(sigStats$EntrezID,sigStats$type,length) #997 exons,382 genes, 924 junctions
gList = split(sigStats$EntrezID, sigStats$type)

##################################
# find GO and KEGG term enrichment
compareKegg = compareCluster(gList,fun = 'enrichKEGG',universe = univ,organism = "mmu", qvalueCutoff = 0.05,pvalueCutoff = 0.1)
compareGoMf = compareCluster(gList, fun = "enrichGO", universe = univ, ont = "MF",OrgDb=org.Mm.eg.db, qvalueCutoff = 0.05, pvalueCutoff = 0.1, readable = TRUE)
compareGoBp = compareCluster(gList, fun = "enrichGO", universe = univ, ont = "BP",OrgDb=org.Mm.eg.db, qvalueCutoff = 0.05, pvalueCutoff = 0.1, readable = TRUE)
compareGoCc = compareCluster(gList, fun = "enrichGO", universe = univ, ont = "CC",OrgDb=org.Mm.eg.db, qvalueCutoff = 0.05, pvalueCutoff = 0.1, readable = TRUE)
compareGo = lapply(list(compareGoMf=compareGoMf,compareGoBp=compareGoBp,compareGoCc = compareGoCc),simplify)

############
# make plots
pdf("plots/gene_sets_tcf4_mouse.pdf")
plot(compareKegg,colorBy="qvalue", font.size =8,title = "KEGG Terms DE Group Comparisons")
plot(compareGo[[1]],colorBy="qvalue", font.size =8, title = "Enriched GO-MF Terms")
plot(compareGo[[2]],colorBy="qvalue", font.size =8, title = "Enriched GO-BP Terms")
plot(compareGo[[3]],colorBy="qvalue", font.size =8, title = "Enriched GO-CC Terms")
dev.off()

######
# save
save(compareKegg, compareGo, file = "rdas/gene_sets_tcf4_mouse.rda")
WriteXLS(lapply(c(compareGo,compareKegg=compareKegg), as.data.frame),ExcelFileName = "tables/gene_sets_tcf4_mouse.xls")
rm(list = ls())

