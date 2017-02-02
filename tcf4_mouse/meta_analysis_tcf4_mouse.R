# compute meta analysis stouffer z-values, from mega TCF4 analysis, mef2c, mecp2, and pten knockouts
library(jaffelab)
library(reshape2)
library(RColorBrewer)
library(lattice)
library(GenABEL)
library(qqman)

###################
# load the datasets
load('rdas/mega_dataset_DE_objects_DESeq2.rda',envir = tcf4<-new.env())
load('../pten/rdas/pten_DE_objects_DESeq2.rda',envir = pten<-new.env())
#load('../mef2c/rdas/mef2c_DE_objects_DESeq2.rda',envir = mef2c<-new.env())
load('../mecp2/rdas/mecp2_DE_objects_DESeq2.rda',envir = mecp2<-new.env())
#load("../ube3a/rdas/ube3a_DE_objects_DESeq2.rda",envir = ube3a<-new.env())

###############################
# extract and merge t-statistics
geneMap = tcf4$outGene[,-c(1:6)]
geneMap = geneMap[!(geneMap$Symbol %in% c('Tcf4','Mef2c','Pten','Mecp2','Ube3a')),]
labs = rownames(geneMap)
outStats= data.frame(TCF4 = tcf4$outGene[labs,4],
                PTEN = pten$outGene[labs,4],
                #MEF2C = mef2c$outGene[labs,4],
                MECP2 = mecp2$outGene[labs,4],
                #UBE3A = ube3a$outGene[labs,4],
                baseMean = tcf4$outGene[labs,1],geneMap)
rownames(outStats) = labs
outSE= data.frame(TCF4 = tcf4$outGene[labs,3], 
                  PTEN = pten$outGene[labs,3], 
                  #MEF2C = mef2c$outGene[labs,3], 
                  MECP2 = mecp2$outGene[labs,3])#, UBE3A = ube3a$outGene[labs,3])

#############################################
# drop missing rows, sort by most significant
gIndex = which(complete.cases(outStats))
geneMap = geneMap[gIndex,]
outStats = outStats[gIndex,]
outSE = outSE[gIndex,]

ind = which(abs(outStats$TCF4) > 3)
signif(cor(outStats[ind,1:3],use="comp"),2)
pairs(outStats[ind,1:3])

###################
# null permutations
W = 1/apply(outSE,2,mean,na.rm = T) # weighting by average standard error
pNull= unlist(parallel::mclapply(seq(300),function(x){
  stat = apply(outStats[,1:3],2,sample)
  z = rowSums(stat*matrix(W,nrow = nrow(stat),ncol = ncol(stat),byrow = T))/sqrt(sum(W^2))
  #z = rowSums(stat)/sqrt(6)
  return(2*pnorm(-abs(z)))
},mc.cores = parallel::detectCores()),use.names = F)
pNull = pNull[order(pNull)]
pecdf = ecdf(pNull)

###############################################
# Find meta p-value using absolute t-statistics
outStats$zscore = rowSums(outStats[,1:3]*matrix(W,nrow = nrow(outStats),ncol = 5,byrow = T))/sqrt(sum(W^2))
#outStats$pvalue = ifelse(outStats$zscore>0,1-zecdf(outStats$zscore),zecdf(outStats$zscore))
outStats$pvalue = pecdf(2*pnorm(-abs(outStats$zscore)))
#outStats$pvalue = 2*pnorm(-abs(outStats$zscore))
outStats$padj = p.adjust(outStats$pvalue,method = 'fdr')
outStats = outStats[order(outStats$pvalue),]
sum(outStats$padj<0.05,na.rm = T)

##############################################
#parametric pvalue qqplot and inflation lambda
qq(2*pnorm(-abs(outStats$zscore)), main = "Q-Q plot of meta parametric p-values")
estlambda(2*pnorm(-abs(outStats$zscore)), plot = TRUE, proportion = .5 )

#permutation pvalue qqplot and inflaction lambda
qq(outStats$pvalue, main = "Q-Q plot of meta permutation p-values")
estlambda(outStats$pvalue, plot = TRUE, proportion = .5 )

tmp = -log10(2*pnorm(-abs(outStats$zscore)))+log10(outStats$pvalue)
plot(outStats$baseMean, tmp,main = 'Diff(-log10P) v. Mean Exprs',col = '#00000050',
     log= 'xy',xlab = 'Mean Expression Counts',ylab = 'Difference of -log10 P',pch = 20)
#write.csv(outStats,file = 'tables/meta_analysis_genes.csv')





















if(FALSE){
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
}
