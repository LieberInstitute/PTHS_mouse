###
library(minfi)
library(SummarizedExperiment)
library(recount)
library(lmerTest)
library(jaffelab)
library(biomaRt)
library(RColorBrewer)

# load data
load("asd/rdas/rse_gene_Parikshak_ASD_human_n205.Rdata")
rse_gene = rse_gene[,-156] # drop flagged sample

## do the RNA-based estiamtion
load("asd/rdas/singleCell_iPSC_quake_coefEsts_calibration_Zscale_adultOnly.rda")
yExprs = log2(getRPKM(rse_gene, "Length")+1)

## pheno
pd = read.delim("asd/tables/ASD_sample_and_RNAseq_metadata.txt", as.is=TRUE)
pd = pd[match(colnames(rse_gene), pd$Sample.Name),]
pd$Diagnosis = factor(ifelse(pd$Detailed.Diagnosis =="Chromosome 15q Duplication Syndrome",
                             'Dup15',pd$Diagnosis), levels = c("CTL", "ASD",'Dup15'))
pd = cbind(pd, colData(rse_gene))

# load qSV's
load('asd/rdas/qSVAs-geschwind_asd.rdas')
qSVs = qSVs[ss(pd$Fastq.file.names, ".R1", 1, fixed=TRUE),]
mod = model.matrix(~Diagnosis + Region + totalAssignedGene + Sequencing.Batch + 
                     Brain.Bank + RIN + Age + Sex, data =pd)
modAdj = cbind(mod,qSVs)
			
###############
### DO CAGS ###
load("tcf4_mouse/rdas/overlap_tcf4_mecp2_pten.rda", verbose=TRUE)
ensembl=useMart("ensembl",dataset = 'mmusculus_gene_ensembl')
MMtoHG = getBM(attributes = c('ensembl_gene_id','hsapiens_homolog_ensembl_gene','hsapiens_homolog_associated_gene_name'),mart = ensembl)
cagIndexG = ss(rownames(yExprs),"\\.") %in% MMtoHG$hsapiens_homolog_ensembl_gene[MMtoHG$ensembl_gene_id %in% sigGenes]

# pcaOverall = prcomp(t(yExprs))
eigenGeneAsd = prcomp(t(yExprs[cagIndexG,]))$x[,1]
pcaVarsAsd = getPcaVars(prcomp(t(yExprs[cagIndexG,])))[1]
coefAsdCAG = summary(lmer(eigenGeneAsd ~  pd$Diagnosis + pd$Region + (1|pd$Brain.ID)))
loads =  prcomp(t(yExprs[cagIndexG,]))$rot[,1]
names(loads) = MMtoHG$hsapiens_homolog_associated_gene_name[match(ss(names(loads),"\\."),
	MMtoHG$hsapiens_homolog_ensembl_gene)] 

## make plots, regress out full model
coefAsdCAG_full = summary(lmer(eigenGeneAsd ~ mod - 1 + (1|pd$Brain.ID)))
cleanEigen_full = cleaningY(matrix(eigenGeneAsd, nr=1), mod, P=3)
cleanEigen_full = cleanEigen_full*-1 # for consistency in plots
table(sign(loads*-1)) # for consistency in plots

pdf("asd/plots/cagBoxplots_ASD_geschwind_cleaned.pdf",h=6,w=7)
par(mar=c(5,6,2,2), cex.axis=1.6,cex.lab=1.6)
palette(brewer.pal(5,"Set1"))
boxplot(cleanEigen_full[1,] ~ pd$Diagnosis , outline = FALSE,
	ylim = range(cleanEigen_full[1,]),
	ylab = paste0("CAG Eigengene: ", pcaVarsAsd, "% Var Expl | Model"))
points(cleanEigen_full[1,] ~ jitter(as.numeric(factor(pd$Diagnosis)), amount=0.15),
	pch = 20+as.numeric(factor(pd$Region)), bg=factor(pd$Diagnosis))
legend("bottomright", paste0("ASD p=", signif(coefAsdCAG_full$coef[2,5],3)),
	cex=1.5,bty="n")
dev.off()

### and then qSVA
coefAsdCAG_full_qsva = summary(lmer(eigenGeneAsd ~ modAdj - 1 + (1|pd$Brain.ID)))
cleanEigen_full_qsva = cleaningY(matrix(eigenGeneAsd, nr=1), modAdj, P=3)
cleanEigen_full_qsva = cleanEigen_full_qsva*-1 # for consistency in plots

pdf("asd/plots/cagBoxplots_ASD_geschwind_cleaned_qsva.pdf",h=6,w=7)
par(mar=c(5,6,2,2), cex.axis=1.6,cex.lab=1.6)
palette(brewer.pal(5,"Set1"))
boxplot(cleanEigen_full_qsva[1,] ~ pd$Diagnosis , outline = FALSE,
	ylim = range(cleanEigen_full_qsva[1,]),
	ylab = paste0("CAG Eigengene: ", pcaVarsAsd, "% Var Expl | Model"))
points(cleanEigen_full_qsva[1,] ~ jitter(as.numeric(factor(pd$Diagnosis)), amount=0.15),
	pch = 20+as.numeric(factor(pd$Region)), bg=factor(pd$Diagnosis))
legend("topright", paste0("ASD p=", signif(coefAsdCAG_full_qsva$coef[2,5],3)),
	cex=1.5)
dev.off()

###################
## all pths genes #
load("tcf4_mouse/rdas/mega_tcf4_ages_DE_objects_DESeq2.rda") #PTHS mouse data
outGeneList = lapply(outGeneList, function(g) {
  g$hsapien_homolog = MMtoHG$hsapiens_homolog_ensembl_gene[match(rownames(g),MMtoHG$ensembl_gene_id)]
  g[grepl('ENSG',g$hsapien_homolog),]})
outGeneAdult = outGeneList$Adult[outGeneList$Adult$padj < 0.05,]
tcf4IndexG = ss(rownames(yExprs),"\\.") %in% MMtoHG$hsapiens_homolog_ensembl_gene[MMtoHG$ensembl_gene_id %in% rownames(outGeneAdult)]

## do PCA
eigenGeneTcf4 = prcomp(t(yExprs[tcf4IndexG,]))$x[,1]
pcaVarsTcf4 = getPcaVars(prcomp(t(yExprs[tcf4IndexG,])))[1]
coefTcf4CAG = summary(lmer(eigenGeneTcf4 ~  pd$Diagnosis + pd$Region + (1|pd$Brain.ID)))
loadsTcf4 =  prcomp(t(yExprs[tcf4IndexG,]))$rot[,1]
names(loadsTcf4) = MMtoHG$hsapiens_homolog_associated_gene_name[match(ss(names(loadsTcf4),"\\."),
	MMtoHG$hsapiens_homolog_ensembl_gene)] 
	

## make plots, regress out full model
coefTcf4CAG_full = summary(lmer(eigenGeneTcf4 ~ mod - 1 + (1|pd$Brain.ID)))
coefTcf4CAG_full_qsva = summary(lmer(eigenGeneTcf4 ~ modAdj - 1 + (1|pd$Brain.ID)))

## regress out region
cleanEigenTcf4 = cleaningY(matrix(eigenGeneTcf4, nr=1), mod, P=2)

pdf("asd/plots/tcf4Boxplots_ASD_geschwind_cleaned.pdf",h=6,w=7)
par(mar=c(5,6,2,2), cex.axis=1.6,cex.lab=1.6)
palette(brewer.pal(5,"Set1"))
boxplot(cleanEigenTcf4[1,] ~ pd$Diagnosis , outline = FALSE,
	ylim = range(cleanEigenTcf4[1,]),
	ylab = paste0("CAG Eigengene: ", pcaVarsAsd, "% Var Expl | Model"))
points(cleanEigenTcf4[1,] ~ jitter(as.numeric(factor(pd$Diagnosis)), amount=0.15),
	pch = 20+as.numeric(factor(pd$Region)), bg=factor(pd$Diagnosis))
legend("bottomright", paste0("ASD p=", signif(coefTcf4CAG_full$coef[2,5],3)),
	cex=1.5,bty="n")
dev.off()

#################################################
## AND CELL TYPES, projection ##
yExprs_Z = scale(yExprs[rownames(coefEsts),])
yExprs_Z = yExprs_Z[,!is.na(colSums(yExprs_Z))]

cellCounts_Exprs = minfi:::projectCellType(yExprs_Z,coefEsts)
cellCounts_Exprs = as.data.frame(cellCounts_Exprs)
			 
# boxplot(cellCounts_Exprs$Oligodendrocytes ~ pd$Diagnosis*pd$Region,
	# ylab = "Oligo Proportion")

## linear model
lapply(cellCounts_Exprs, function(y) {
	signif(summary(lm(y ~ modAdj - 1))$coef,2)
})

####################################
## linear mixed effects modeling  ##
####################################

## just region
coefListLme_reg = lapply(cellCounts_Exprs, function(y) {
	signif(summary(lmer(y ~ pd$Diagnosis + pd$Region + (1|pd$Brain.ID)))$coef,2)
})
coefAsd_reg = t(sapply(coefListLme_reg, function(x) x[2,]))
coef22_reg = t(sapply(coefListLme_reg, function(x) x[3,]))
coefAsd_reg

## linear mixed effects model
coefListLme_qsva = lapply(cellCounts_Exprs, function(y) {
	signif(summary(lmer(y ~ modAdj - 1 + (1|pd$Brain.ID)))$coef,2)
})
coefAsd = t(sapply(coefListLme_qsva, function(x) x[2,]))
coef22 = t(sapply(coefListLme_qsva, function(x) x[3,]))
coefAsd

coefListLme_obs = lapply(cellCounts_Exprs, function(y) {
	signif(summary(lmer(y ~ mod - 1 + (1|pd$Brain.ID)))$coef,2)
})
coefAsd_obs = t(sapply(coefListLme_obs, function(x) x[2,]))
coef22_obs = t(sapply(coefListLme_obs, function(x) x[3,]))
coefAsd_obs

# ## just cortical p-values
# ctxInd = which(pd$Region != "vermis")
# coefListLmeCtx = lapply(cellCounts_Exprs, function(y) {
	# signif(summary(lmer(y[ctxInd] ~ modAdj[ctxInd,-5] - 1 + (1|pd$Brain.ID[ctxInd])))$coef,2)
# })
# t(sapply(coefListLmeCtx, function(x) x[2,]))
# t(sapply(coefListLmeCtx, function(x) x[3,]))

# ## just cerebellum p-values
# crbInd = which(pd$Region == "vermis")
# coefListLmeCrb = lapply(cellCounts_Exprs, function(y) {
	# signif(summary(lmer(y[crbInd] ~ mod[crbInd,]+ 
		# (1|factor(pd$Brain.ID[crbInd]))))$coef,2)
# })
# t(sapply(coefListLmeCrb, function(x) x[2,]))
# t(sapply(coefListLmeCrb, function(x) x[3,]))

	
## boxplots, clean ###
countsAdj = as.data.frame(t(cleaningY(t(cellCounts_Exprs), modAdj, P=3)))

pdf("asd/plots/boxplot_of_cell_counts_asd_geschwind_clean.pdf",h=6,w=7)
palette(brewer.pal(5,"Set1"))
par(mar=c(5,6,3,2), cex.axis=2,cex.lab=2)
boxplot(countsAdj$Oligodendrocytes ~ pd$Diagnosis,
	ylab = "Oligo RNA Fraction | Model", outline=FALSE,
	ylim = range(countsAdj$Oligodendrocytes))
points(countsAdj$Oligodendrocytes ~ 
	jitter(as.numeric(pd$Diagnosis), amount=0.15),pch =21,bg=pd$Diagnosis)
legend("topright", paste0("ASD p=", signif(coefAsd[4,5],3)),cex=1.5)
boxplot(countsAdj$OPC ~ pd$Diagnosis,
	ylab = "OPC RNA Fraction | Model", outline=FALSE,
	ylim = range(countsAdj$OPC))
points(countsAdj$OPC ~ 
	jitter(as.numeric(pd$Diagnosis), amount=0.15),pch =21,bg=pd$Diagnosis)
legend("topleft", paste0("ASD p=", signif(coefAsd[2,5],3)),cex=1.5)
dev.off()


#######
## big boxplot
library(tidyr)
cellCounts_Exprs$SampleID = rownames(cellCounts_Exprs)
long = gather(cellCounts_Exprs, CellType, Proportion, -SampleID )
long$Dx = pd$Diagnosis[match(long$SampleID, pd$SAMPLE_ID)]
long$Region = pd$Region[match(long$SampleID, pd$SAMPLE_ID)]
long$Label = paste0(long$CellType, ":", long$Dx)

long$Label = factor(long$Label, levels = c("Oligodendrocytes:CTL",
	"Oligodendrocytes:ASD","Oligodendrocytes:Dup15", 
	"OPC:CTL",  "OPC:ASD", "OPC:Dup15", 
	"Microglia:CTL", "Microglia:ASD","Microglia:Dup15",
	"Neurons:CTL", "Neurons:ASD", "Neurons:Dup15", 
	"Astrocytes:CTL", "Astrocytes:ASD" ,  "Astrocytes:Dup15" ,  
    "Endothelial:CTL", "Endothelial:ASD", "Endothelial:Dup15"))
    
pdf("asd/plots/cellType_changes.pdf",h=5,w=10)
par(mar=c(7,6,2,2), cex.axis=1.6,cex.lab=1.6)
palette(brewer.pal(5,"Set1"))
boxplot(Proportion ~ Label, data=long, ylim = c(0,1),
	names = rep(c("CTL", "ASD","Dup15"), times=6), 
	outline=FALSE,ylab = "RNA class proportion",las=3,
	xlim = c(0.75,18.25))
points(Proportion ~ jitter(as.numeric(Label),amount=0.15),
	data = long, pch = 21, bg = factor(Region))
text(x = seq(2,18,by=3), y = 0.95, cex=1.5,
	c("Oligo", "OPC", "Microglia", "Neuron", "Astro", "Endothelial"))
abline(v = seq(3.5,15.5,by=3), lty=2)
dev.off()

longList = split(long, long$CellType)
lapply(longList, function(dat) {
	signif(summary(lm(Proportion ~ Dx+ Region,
		data=dat))$coef,3)
})


##############################
## schizophrenia comparison ##
##############################

coefEstsSz = as.matrix(read.csv("/users/ajaffe/Lieber/Projects/RNAseq/firstRnaSeqPaper/revision/BrainSeq_Phase1_DLPFC_darmanis_hg19-based_cellComp.csv",
	row.names=1))
load("/dcs01/ajaffe/Brain/DLPFC_PolyA/szControl/data/rpkmCounts_szControlDLPFC.rda")
geneRpkmBrain = geneRpkm # just gene

geneExprsBrain = log2(geneRpkmBrain+1)
cellCountsSz_Exprs = minfi:::projectCellType(geneExprsBrain[rownames(coefEstsSz),], coefEstsSz)
cellCountsSz_Exprs = as.data.frame(cellCountsSz_Exprs)

pd = cbind(pd, cellCountsSz_Exprs)

## filter people
keepIndex = which(pd$Age > 17 & 
	pd$totalAssignedGene > 0.5 & 
	pd$mappingRate > 0.7 & pd$RIN > 6 & 
	pd$Age < 80 & pd$snpPC2 < 0.03)
pdSz = pd[keepIndex,]
geneExprsBrainSz = geneExprsBrain[,keepIndex]
## model
mod = model.matrix(~ Dx + Age + Sex + 
	snpPC1 + snpPC5 + snpPC6 + snpPC9 + snpPC10 +
	mitoRate + RIN + totalAssignedGene, data=pdSz)
coefSz = t(apply(pdSz[,57:64], 2, function(y) {
	summary(lm(y ~ mod - 1))$coef[2,]
}))


## big boxplot
cellCountsSz_Exprs = cellCountsSz_Exprs[keepIndex,]
library(tidyr)
cellCountsSz_Exprs$SampleID = rownames(cellCountsSz_Exprs)
longSz = gather(cellCountsSz_Exprs, CellType, Proportion, -SampleID )
longSz$Dx = ifelse(pd$Dx[match(longSz$SampleID, pd$RNum)] == "Control", "CONT", "SCZD")
longSz = longSz[!grepl("Fetal", longSz$CellType),]
longSz$Label = paste0(longSz$CellType, ":", longSz$Dx)


longSz$Label = factor(longSz$Label, levels = c("Oligodendrocytes:CONT",
	"Oligodendrocytes:SCZD", "OPC:CONT",  "OPC:SCZD",
	"Microglia:CONT", "Microglia:SCZD",
	"Neurons:CONT", "Neurons:SCZD", 
	"Astrocytes:CONT", "Astrocytes:SCZD" ,  
    "Endothelial:CONT", "Endothelial:SCZD"))
    
pdf("asd/plots/cellType_changes_schizophrenia.pdf",h=5,w=10)
par(mar=c(7,6,2,2), cex.axis=1.6,cex.lab=1.6)
palette(brewer.pal(5,"Set2"))
boxplot(Proportion ~ Label, data=longSz, ylim = c(0,1.05),
	names = rep(c("CONT", "SCZD"), times=6), 
	outline=FALSE,ylab = "RNA class proportion",las=3,
	xlim = c(0.75,12.25))
points(Proportion ~ jitter(as.numeric(Label),amount=0.15),
	data = longSz, pch = 21, bg = factor(Dx))
text(x = seq(1.5,12,by=2), y = 1.03, cex=1.5,
	c("Oligo", "OPC", "Microglia", "Neuron", "Astro", "Endothelial"))
abline(v = seq(2.5,12.5,by=2), lty=2)
dev.off()


#### do CAGs

cagIndex = rownames(geneExprsBrain) %in% MMtoHG$hsapiens_homolog_ensembl_gene[MMtoHG$ensembl_gene_id %in% sigGenes]
eigenGeneSz = prcomp(t(geneExprsBrainSz[cagIndex,]))$x[,1]
pcaVarsSz = getPcaVars(prcomp(t(geneExprsBrainSz[cagIndex,])))[1]
coefSzCAG = summary(lm(eigenGeneSz ~ mod - 1))

cleanEigen = cleaningY(matrix(eigenGeneSz, nr=1), mod, P=2)

pdf("asd/plots/cagBoxplots_schizophrenia.pdf",h=6,w=6)
par(mar=c(5,6,2,2), cex.axis=1.6,cex.lab=1.6)
palette(brewer.pal(5,"Set2"))
boxplot(cleanEigen[1,] ~ pdSz$Dx, outline = FALSE,
	ylim = range(cleanEigen[1,]),
	ylab = paste0("CAG Eigengene: ", pcaVarsSz, "% Var Expl (Adj)"))
points(cleanEigen[1,] ~ jitter(as.numeric(factor(pdSz$Dx)), amount=0.15),
	pch = 21, bg=factor(pdSz$Dx))
legend("bottomright", paste0("p=", signif(coefSzCAG$coef[2,4],3)),cex=1.5)
dev.off()

### Asd