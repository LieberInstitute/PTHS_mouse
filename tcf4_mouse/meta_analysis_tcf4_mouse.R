# compute meta analysis stouffer z-values, from mega TCF4 analysis, mef2c, mecp2, and pten knockouts
library(jaffelab)
library(reshape2)
library(RColorBrewer)
library(lattice)

###################
# load the datasets
load('rdas/mega_dataset_DE_objects_DESeq2.rda',envir = tcf4<-new.env())
load('../pten/rdas/pten_DE_objects_DESeq2.rda',envir = pten<-new.env())
load('../mef2c/rdas/mef2c_DE_objects_DESeq2.rda',envir = mef2c<-new.env())
load('../mecp2/rdas/mecp2_DE_objects_DESeq2.rda',envir = mecp2<-new.env())

###############################
# extract and merge t-statistics
geneMap = tcf4$outGene[,-c(1:6)]
geneMap = geneMap[!(geneMap$Symbol %in% c('Tcf4', 'Pten', 'Mef2c','Mecp2')),]
labs = rownames(geneMap)
outStats= data.frame(TCF4 = tcf4$outGene[labs,4],
                PTEN = pten$outGene[labs,4],
                #MEF2C = mef2c$outGene[labs,4],
                MECP2 = mecp2$outGene[labs,4])
rownames(outStats) = labs

#############################################
# drop missing rows, sort by most significant
gIndex = which(complete.cases(outStats))
geneMap = geneMap[gIndex,]
outStats = outStats[gIndex,]

gIndex = order(apply(outStats,1,function(x) abs(mean(x))),decreasing = T)
geneMap = geneMap[gIndex,]
outStats = outStats[gIndex,]

######################################
# make heatmap of top concordant genes
dat = melt(cbind(ENS = rownames(outStats)[1:20],outStats[1:20,]),id.vars = 'ENS',value.name = 'stat')
dat$Symbol = geneMap[dat$ENS,c('Symbol')]
dat$variable = factor(dat$variable,levels = c('TCF4', 'PTEN', "MECP2"))
dat$Symbol = factor(dat$Symbol)

#range(dat$stat)
theSeq = seq(-13,13,by=0.01) 
my.col <- colorRampPalette(brewer.pal(11,"PiYG"))(length(theSeq))
print(levelplot(stat ~ variable + Symbol, 
                data= dat, at = theSeq,pretty=TRUE,
                col.regions = my.col, scales=list(y=list(cex=1.15), 
                                                  x=list(rot=30, cex=1.15)),
                ylab = "", xlab = ""))



############################
# order by absolute stats
gIndex = order(apply(outStats,1,function(x) mean(abs(x))),decreasing = T)
geneMap = geneMap[gIndex,]
outStats = outStats[gIndex,]

######################################
# make heatmap of top concordant genes
dat = melt(cbind(ENS = rownames(outStats)[1:20],outStats[1:20,]),id.vars = 'ENS',value.name = 'stat')
dat$Symbol = geneMap[dat$ENS,c('Symbol')]
dat$variable = factor(dat$variable,levels = c('TCF4', 'PTEN', "MECP2"))
dat$Symbol = factor(dat$Symbol)

#range(dat$stat)
theSeq = seq(-13,13,by=0.01) 
my.col <- colorRampPalette(brewer.pal(11,"PiYG"))(length(theSeq))
print(levelplot(stat ~ variable + Symbol, 
                data= dat, at = theSeq,pretty=TRUE,
                col.regions = my.col, scales=list(y=list(cex=1.15), 
                                                  x=list(rot=30, cex=1.15)),
                ylab = "", xlab = ""))


###################
# null permutations
if(FALSE){
zNull= unlist(parallel::mclapply(seq(100),function(x){
  stat = apply(outStats[,1:3],2,sample)
  z = rowMeans(abs(stat))/sqrt(ncol(stat))
  return(z)
},mc.cores = parallel::detectCores()),use.names = F)
zNull = zNull[order(zNull)]
zecdf = ecdf(zNull)
save(zecdf,file = '/dcl01/lieber/ajaffe/Brady/mouseRNAseq/meta_null_stats.rda')
} else{
  load('/dcl01/lieber/ajaffe/Brady/mouseRNAseq/meta_null_stats.rda')
}

###############################################
# Find meta p-value using absolute t-statistics
outStats$zscore = apply(abs(outStats[,1:3]),1,function(x) mean(x))/sqrt(3)
outStats$pvalue = ifelse(outStats$zscore>0,1-zecdf(outStats$zscore),zecdf(outStats$zscore))
outStats$padj = p.adjust(outStats$pvalue,method = 'fdr')
outStats = outStats[order(outStats$pvalue),]
sum(outStats$padj<0.05,na.rm = T)


