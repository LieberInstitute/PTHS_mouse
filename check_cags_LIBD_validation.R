##
library(jaffelab)
library(GEOquery)
library(biomaRt)
library(RColorBrewer)

theData = getGEO("GSE102741")[[1]]
download.file("https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE102741&format=file&file=GSE102741%5Flog2RPKMcounts%2Erda%2Egz",
	destfile = "asd/rdas/ASD_GEO_data.rda")
load("asd/rdas/ASD_GEO_data.rda")

## expression data
p = geneRpkm2
pd = pData(theData)

colnames(pd)[54:64] = ss(colnames(pd)[54:64], ":")
colnames(pd)[54:64] = gsub(" ", "_", colnames(pd)[54:64] )

pd$Dx = factor(ifelse(pd$disease_status == "Healthy control", "CONT", "ASD"),
	levels = c("CONT", "ASD"))

	
#### do CAGs
load("tcf4_mouse/rdas/overlap_tcf4_mecp2_pten.rda", verbose=TRUE)
ensembl=useMart("ensembl",dataset = 'mmusculus_gene_ensembl')
MMtoHG = getBM(attributes = c('ensembl_gene_id',
	'hsapiens_homolog_ensembl_gene', 'hsapiens_homolog_associated_gene_name'),mart = ensembl)
cagIndex = rownames(p) %in% MMtoHG$hsapiens_homolog_ensembl_gene[MMtoHG$ensembl_gene_id %in% sigGenes]

eigenGene = prcomp(t(p[cagIndex,]))$x[,1]
pcaVars = getPcaVars(prcomp(t(log2(p[cagIndex,]+1))))[1]
loads = prcomp(t(p[cagIndex,]))$rot[,1]
names(loads) = MMtoHG$hsapiens_homolog_associated_gene_name[match(names(loads),
	MMtoHG$hsapiens_homolog_ensembl_gene)] 
boxplot(eigenGene ~ pd$Dx)

f = lm(eigenGene ~ Dx + as.numeric(rin) + 
	as.numeric(percentexonicmapping) + Sex, data=pd)
summary(f)

pcs = prcomp(t(p))
fPC = lm(eigenGene ~ Dx + as.numeric(rin) + 
	as.numeric(percentexonicmapping) + Sex + pcs$x[,1:9], data=pd)
summary(fPC)

mod = model.matrix(~Dx + as.numeric(rin) + 
	as.numeric(percentexonicmapping) + Sex + pcs$x[,1:9], data=pd)
ff = summary(lm(eigenGene ~ mod - 1))

cleanEigen = cleaningY(matrix(eigenGene, nr=1), mod, P=2)
pdf("cagBoxplots_ASD_LIBD.pdf",h=6,w=7)
par(mar=c(5,6,2,2), cex.axis=1.6,cex.lab=1.6)
palette(brewer.pal(5,"Set1"))
boxplot(cleanEigen[1,] ~ pd$Dx, outline = FALSE,
	ylim = range(cleanEigen[1,]),
	ylab = paste0("CAG Eigengene: ", pcaVars, "% Var Expl | Model"))
points(cleanEigen[1,] ~ jitter(as.numeric(factor(pd$Dx)), amount=0.15),
	pch = 21, bg=factor(pd$Dx))
legend("top", paste0("p=", signif(ff$coef[2,4],3)),cex=1.5)
dev.off()

colnames(pd)[54:64] = ss(colnames(pd)[54:64], ":")

###################
## check counts ###


## do the RNA-based estiamtion
load("asd/rdas/singleCell_iPSC_quake_coefEsts_calibration_Zscale_adultOnly.rda")

# project
yExprs_Z = scale(p[ss(rownames(coefEsts),"\\."),])

cellCounts_Exprs = minfi:::projectCellType(yExprs_Z,coefEsts)
cellCounts_Exprs = as.data.frame(cellCounts_Exprs)

boxplot(cellCounts_Exprs$Oligodendrocytes ~ pd$Dx,
	ylab = "Oligo Proportion")
	
fOligo = summary(lm(cellCounts_Exprs$Oligodendrocytes ~ mod -1 ))
fOPC = summary(lm(cellCounts_Exprs$OPC ~ mod -1 ))
	
## boxplots ###
pdf("asd/plots/boxplot_of_cell_counts_asd.pdf")
palette(brewer.pal(5,"Set1"))
par(mar=c(5,6,3,2), cex.axis=2,cex.lab=2)
boxplot(cellCounts_Exprs$Oligodendrocytes ~ pd$Dx,
	ylab = "Oligo RNA Fraction", outline=FALSE,
	ylim = range(cellCounts_Exprs$Oligodendrocytes))
points(cellCounts_Exprs$Oligodendrocytes ~ 
	jitter(as.numeric(pd$Dx), amount=0.15),pch =21,bg=pd$Dx)
legend("topright", paste0("p=", signif(fOligo$coef[2,4],3)),cex=1.5)
boxplot(cellCounts_Exprs$OPC ~ pd$Dx,
	ylab = "OPC RNA Fraction", outline=FALSE,
	ylim = range(cellCounts_Exprs$OPC))
points(cellCounts_Exprs$OPC ~ 
	jitter(as.numeric(pd$Dx), amount=0.15),pch =21,bg=pd$Dx)
legend("topleft", paste0("p=", signif(fOPC$coef[2,4],3)),cex=1.5)
dev.off()
	
## boxplots, clean ###
countsAdj = as.data.frame(t(cleaningY(t(cellCounts_Exprs), mod, P=2)))
pdf("asd/plots/boxplot_of_cell_counts_asd_clean.pdf",w=5)
palette(brewer.pal(5,"Set1"))
par(mar=c(5,6,3,2), cex.axis=2,cex.lab=2)
boxplot(countsAdj$Oligodendrocytes ~ pd$Dx,
	ylab = "Oligo RNA Fraction | Model", outline=FALSE,
	ylim = range(countsAdj$Oligodendrocytes))
points(countsAdj$Oligodendrocytes ~ 
	jitter(as.numeric(pd$Dx), amount=0.15),pch =21,bg=pd$Dx)
legend("topright", paste0("p=", signif(fOligo$coef[2,4],3)),cex=1.5)
boxplot(countsAdj$OPC ~ pd$Dx,
	ylab = "OPC RNA Fraction | Model", outline=FALSE,
	ylim = range(countsAdj$OPC))
points(countsAdj$OPC ~ 
	jitter(as.numeric(pd$Dx), amount=0.15),pch =21,bg=pd$Dx)
legend("topleft", paste0("p=", signif(fOPC$coef[2,4],3)),cex=1.5)
dev.off()
