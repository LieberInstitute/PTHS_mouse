###
library(GenomicRanges)
library(rafalib)
library(sva)
library(limma)

getPcaVars = function(pca)  signif(((pca$sdev)^2)/(sum((pca$sdev)^2)),3)*100

load("/dcl01/lieber/ajaffe/psychENCODE_Data/UCLA_R01MH094714/rpkmCounts_asd.rda")

plot(pd$Proportion.of.exonic.reads, pd$totalAssignedGene)

## filter
gIndex = which(rowMeans(geneRpkm) > 0.1)
geneRpkm = geneRpkm[gIndex,]
geneMap = geneMap[gIndex,]

# #### degradation files
# degFiles = paste0("Degrade/", 
	# pd$Sample.Name, "_degradeStats.txt")
# all(file.exists(degFiles))
# degCov = sapply(degFiles, function(x) {
	# cat(".")
	# read.delim(pipe(paste("cut -f10", x)), 
		# as.is=TRUE)$sum
# })
# colnames(degCov) = pd$Sample.Name
# degCov = degCov/100 # read length
# bg = matrix(rep(pd$totalMapped/80e6), nc = nrow(pd), 
	# nr = nrow(degCov),	byrow=TRUE)
# degCovAdj = degCov/bg
# save(degCovAdj, file="rdas/degradation_mat_UCLA_ASD_RiboZero_v2.rda")
load("dcl01/lieber/ajaffe/psychENCODE_Data/UCLA_R01MH094714/rdas/degradation_mat_UCLA_ASD_RiboZero_v2.rda")

## do PCA
degPca = prcomp(t(log2(degCovAdj+1)))
degPcaVars = getPcaVars(degPca)
degPcaVars[1:10]

anovaList = apply(degPca$x[,1:10], 2, function(x) {
	anova(lm(x ~ pd$Diagnosis + pd$Region))
})
anovaPval = t(sapply(anovaList, function(x) x[1:2,5]))
colnames(anovaPval) = c("Dx","Region")

boxplot(degPca$x[,4] ~ pd$Diagnosis*pd$Region)

genePca= prcomp(t(log2(geneRpkm+1)))
genePcaVars = getPcaVars(genePca)
genePcaVars[1:10]

## degradation features
mypar(3,3)
for(i in 1:9){
	plot(pd$totalAssignedGene, degPca$x[,i], 
		pch = 21, bg = factor(pd$Region),
		main = paste0("PC",i, ": ", degPcaVars[i],
			"% Var Expl"), ylab="", xlab="Assign Rate")
}

for(i in 1:9){
	plot(pd$totalAssignedGene, genePca$x[,i], 
		pch = 21, bg = factor(pd$Region),
		main = paste0("PC",i, ": ", genePcaVars[i],
			"% Var Expl"), ylab="", xlab="Assign Rate")
}

mypar()
plot(degPca$x[,1] ~ pd$RIN, pch = 21, bg = factor(pd$Region),
	xlab= "RIN", ylab="deg PC1") 
plot(pd$totalAssignedGene, pd$Median.3prime.Bias.picard)
plot(pd$totalAssignedGene, pd$Median.5prime.Bias.picard)
plot(pd$totalAssignedGene, pd$Median.5to3prime.Bias.picard)
plot(pd$totalAssignedGene, pd$Proportion.of.exonic.reads)

##################
## DE ######
rIndexes = split(1:nrow(pd), pd$Region)
pd$Diagnosis = factor(pd$Diagnosis, levels = c("CTL", "ASD"))

statListAdj = lapply(rIndexes, function(ii) {
	cat(".")
	mod = model.matrix(~Diagnosis + totalAssignedGene +
		Sequencing.Batch + Brain.Bank + RIN + Age + Sex, data =pd[ii,])
	fit = lmFit(log2(geneRpkm[,ii]+1) , mod)
	eb = ebayes(fit)
	data.frame(log2FC = fit$coef[,2], tstat =eb$t[,2],
		pval = eb$p[,2], fdr = p.adjust(eb$p[,2], "fdr"))
})
cor(sapply(statListAdj, function(x) x$log2FC))
colSums(sapply(statListAdj, function(x) x$fdr) < 0.05)

statListQual = lapply(rIndexes, function(ii) {
	cat(".")
	mod = model.matrix(~Diagnosis + totalAssignedGene +
		Sequencing.Batch + Brain.Bank + RIN + Age + Sex, data =pd[ii,])
	degPca = prcomp(t(log2(degCovAdj[,ii]+1))) # do PCA
	k = sva::num.sv(log2(degCovAdj[,ii]+1), mod) # in sva package
	qSVs = degPca$x[,1:k] # identify quality surrogate variables
	fit = lmFit(log2(geneRpkm[,ii]+1) , cbind(mod, qSVs))
	eb = ebayes(fit)
	data.frame(log2FC = fit$coef[,2], tstat =eb$t[,2],
		pval = eb$p[,2], fdr = p.adjust(eb$p[,2], "fdr"))
})
cor(sapply(statListQual, function(x) x$log2FC))
colSums(sapply(statListQual, function(x) x$fdr) < 0.05)

plot(statListAdj[[1]]$log2FC, statListQual[[1]]$log2FC )

statsAdj = do.call("cbind", statListAdj)
colnames(statsAdj) = paste0("Adj_", colnames(statsAdj))
statsQual = do.call("cbind", statListQual)
colnames(statsQual) = paste0("Qual_", colnames(statsQual))
stats = cbind(statsAdj, statsQual)
colnames(stats) = gsub("-", "_", colnames(stats))

out = cbind(stats, geneMap)
save(out, file="rdas/ASD_DE_gene.rda")

### compare to tissue on a bench ###
load("/users/ajaffe/Lieber/Projects/RnaQualPaper/Degradation/rdas/degradeStats_byLibrary_final.rda")
degradeStatsGene = degradeStats[rownames(stats),c(2,5,8,10:16,18)]

cor(stats[,seq(1,21,by=4)], degradeStatsGene$ribo_timeFC)
cor(stats[,seq(2,22,by=4)], degradeStatsGene$ribo_Tstat)

mypar(2,3)
for(i in seq(1,21,by=4)) {
	plot(stats[,i], degradeStatsGene$ribo_timeFC,
		pch = 21, bg="grey", ylab="Time log2FC",
		xlab= "ASD log2FC")
}
for(i in seq(2,22,by=4)) {
	plot(stats[,i], degradeStatsGene$ribo_Tstat,
		pch = 21, bg="grey", ylab="Time T-stat",
		xlab= "ASD T-stat")
}

or = function(mat) mat[1,1]/mat[1,2]/mat[2,1]*mat[2,2]
sapply(seq(2,22,by=4), function(i) {
	sIndex = stats[,(i+1)] < 0.01
	tt = table(sign(stats[sIndex,i]), 
		sign(degradeStatsGene$ribo_Tstat[sIndex]))
	or(tt[c("-1","1"),c("-1","1")])
})

###### which significant
pvalCuts = c(0.05, 0.01, 0.005, 0.001, 1e-4, 1e-5)
names(pvalCuts) = paste0("p<", pvalCuts)

#### adjustment models
sapply(pvalCuts, function(x) {
		sum(stats$Adj_ba41_42_22.pval < x & 
			sign(stats$Adj_ba41_42_22.log2FC) == 
				sign(stats$Adj_ba9.log2FC),
			na.rm=TRUE) / sum(stats$Adj_ba41_42_22.pval < x)
})
sapply(pvalCuts, function(x) {
	sum(stats$Qual_ba41_42_22.pval < x & 
		sign(stats$Qual_ba41_42_22.log2FC) == 
			sign(stats$Qual_ba9.log2FC),
		na.rm=TRUE) / sum(stats$Qual_ba41_42_22.pval < x)
})

## and sig
sapply(pvalCuts, function(x) {
		sum(stats$Adj_ba41_42_22.pval < x & 
			sign(stats$Adj_ba41_42_22.log2FC) == 
				sign(stats$Adj_ba9.log2FC) & 
			stats$Adj_ba9.pval < 0.05,
			na.rm=TRUE) / sum(stats$Adj_ba41_42_22.pval < x)
})

sapply(pvalCuts, function(x) {
	sum(stats$Qual_ba41_42_22.pval < x & 
		sign(stats$Qual_ba41_42_22.log2FC) == 
			sign(stats$Qual_ba9.log2FC) & 
		stats$Qual_ba9.pval < 0.05,
		na.rm=TRUE) / sum(stats$Qual_ba41_42_22.pval < x)
})
