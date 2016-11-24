# use WGCNA to find common gene groups in syndromic ASD
cleaningY = function(y, mod, P=ncol(mod)) {
  Hat=solve(t(mod)%*%mod)%*%t(mod)
  beta=(Hat%*%t(y))
  cleany=y-t(as.matrix(mod[,-c(1:P)])%*%beta[-c(1:P),])
  return(cleany)
}
library(WGCNA)
library(DESeq2)
library(reshape2)
library(jaffelab)
library(RColorBrewer)
library(beeswarm)
options(stringsAsFactors = FALSE)
enableWGCNAThreads()

################
# load the count 
load('/dcl01/lieber/ajaffe/Brady/mouseRNAseq/mega_asd_mice_DESeq2_svaAdj.rda')
geneCounts = counts(geneDds, normalize = TRUE)
pd = as.data.frame(colData(geneDds))

###########################################################
# regress out surrogate variables and gene assignment rates
modGene = model.matrix(design(geneDds),data = pd)
ind = which(grepl('SV',colnames(modGene))|grepl('Age',colnames(modGene)) | colnames(modGene)=='totalAssignedGene')
modGene = cbind(modGene[,-ind],modGene[,ind]) #ncol = 13
geneClean = cleaningY(geneCounts,modGene,P = 10)

######################################
#need to transpose expression dataframe
datExpr = as.data.frame(t(geneClean))

#########################
#check for missing values 
gsg = goodSamplesGenes(datExpr)
if (!gsg$allOK) datExpr = datExpr[gsg$goodSamples, gsg$goodGenes]

###################
#check for outliers
sampleTree = hclust(dist(datExpr), method = "centroid");
par(mar = c(0,4,2,0),cex = 0.6)
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)
abline(h = 100, col = "red")

#########################################
# Determine power from network topology analysis function
powers = c(c(1:10), seq(from = 12, to=20, by=2))
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)

#####################
# Plot power analysis
pdf('plots/wgcna_mega_asd_mice_power_analysis.pdf',height = 5, width = 9)
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()

###########################
# run WGCNA 1 step function
net = blockwiseModules(datExpr, power = 7, #determined from scale-free topology
                       TOMType = "unsigned", minModuleSize = 30,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = FALSE, verbose = 3)

####################
# plot WGCNA modules
pdf('plots/wgcna_mega_asd_mice_plot.pdf',width = 12, height = 9)
# Convert labels to colors for plotting
mergedColors = labels2colors(net$colors)
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors", dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()

###############
# save the data
moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs;
geneTree = net$dendrograms[[1]];
save(MEs, moduleLabels, moduleColors, geneTree, file = "rdas/mega_asd_mice_wgcna_ageless.RData")

###########################################
# test sample modules for general genotypes
cols = c('SAMPLE_ID',"Genotype",'Mutation','Age')
meCols = names(MEs)
meCols = meCols[order(as.numeric(ss(meCols,'ME',2)))]

pvalues = sapply(meCols,function(c){
pd2 = cbind(pd[,cols],ME = MEs[,c])
return(t.test(ME~Genotype,data = pd2)$p.value)
})

pvalues[which(pvalues< 0.05)]
padj = p.adjust(pvalues,method = 'fdr')
padj[which(padj< 0.05)]

cols = names(pvalues[which(pvalues< 0.05)])
n = length(unique(colData(geneDds)$Mutation))
col = brewer.pal(n,'Paired')
par(mar = c(7,4,2,2))
pdf('plots/wgcna_module_eigengenes.pdf')
for (i in cols){
  boxplot(MEs[,i]~Genotype*Mutation,data = pd,
          main = i,las = 2, ylab = 'Eigengene',col = rep(col,each = 2))
  beeswarm(MEs[,i]~Genotype*Mutation,data = pd,
           pch = 21,pwbg = as.numeric(pd$Genotype)-1,add = T)
}
dev.off()


