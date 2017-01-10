# find replication in sweatt dataset and philpot dataset
# plot replicated genes
library(jaffelab)
library(lattice)
library(RColorBrewer)
library(reshape2)
library(WriteXLS)

###############################################
# load discovery genes and subset to fdr < 0.05
load('tcf4_mouse/rdas/mouse_tcf4_ages_DE_objects_DESeq2.rda',envir = discovery<- new.env())
discovery = lapply(discovery$outGeneList,function(x) subset(x,padj< 0.05 & !is.na(padj)))[c(1,3)]
queryList = lapply(discovery,rownames)
discovery = lapply(discovery,function(x) {
  x[x$padj<0.05 & !is.na(x$padj),]
  names(x)= paste0('Discovery.',names(x))
  x
})

###############################################
# load replicate dataset philpot and sweatt
load('sweatt/rdas/sweatt_DE_objects_DESeq2.rda',envir = sweatt<- new.env())
load('philpot/rdas/philpot_DE_by_age_objects_DESeq2.rda',envir = philpot<- new.env())

####################
# P1 DEG replication
ind = grep('P1',names(philpot$outGeneList))
sigGeneP1 = do.call('cbind',lapply(philpot$outGeneList[ind],function(x) x[queryList[['p1']],1:6]))
names(sigGeneP1) = ss(names(sigGeneP1),'P1.',2)

##############################################################
# determine if genes that replicated in all or any other model
Line = c("Act",'Nest', 'Del', 'R579W', 'CA1')
dis = discovery[['p1']]
sigGeneP1$replicatedAny = sapply(rownames(sigGeneP1),function(g) {
  repPval = sigGeneP1[g,grep('pvalue',names(sigGeneP1))]< 0.05
  repFC = sign(sigGeneP1[g,grep('log2FoldChange',names(sigGeneP1))])== 
    sign(dis[g,'Discovery.log2FoldChange'])
  return(any(repPval & repFC))})
sigGeneP1$replicatedAll = sapply(rownames(sigGeneP1),function(g) {
  repPval = sigGeneP1[g,grep('pvalue',names(sigGeneP1))]< 0.05
  repFC = sign(sigGeneP1[g,grep('log2FoldChange',names(sigGeneP1))])== 
    sign(dis[g,'Discovery.log2FoldChange'])
  return(all(repPval & repFC))})
sigGeneP1$numReplicated = sapply(rownames(sigGeneP1),function(g) {
  repPval = sigGeneP1[g,grep('pvalue',names(sigGeneP1))]< 0.05
  repFC = sign(sigGeneP1[g,grep('log2FoldChange',names(sigGeneP1))])== 
    sign(dis[g,'Discovery.log2FoldChange'])
  return(sum(repPval & repFC))})
sigGeneP1 = cbind(dis,sigGeneP1)
sum(sigGeneP1$replicatedAny,na.rm = T) #25 
sum(sigGeneP1$replicatedAll,na.rm = T) #4


#######################
# Adult DEG replication
names(sweatt$outGene) = paste0('Adult.CA1.',names(sweatt$outGene))

#combine datasets sweatt with philpot
ind = grep('Adult',names(philpot$outGeneList))
sigGeneAdult = cbind(do.call('cbind',lapply(philpot$outGeneList[ind],function(x) x[queryList[['Adult']],1:6])), sweatt$outGene[queryList[['Adult']],1:6])
names(sigGeneAdult) = ss(names(sigGeneAdult),'Adult.',2)

dis = discovery[['Adult']]
sigGeneAdult$replicatedAny = sapply(rownames(sigGeneAdult),function(g) {
  repPval = sigGeneAdult[g,grep('pvalue',names(sigGeneAdult))]< 0.05
  repFC = sign(sigGeneAdult[g,grep('log2FoldChange',names(sigGeneAdult))])== 
    sign(dis[g,'Discovery.log2FoldChange'])
  return(any(repPval & repFC))})
sigGeneAdult$replicatedAll = sapply(rownames(sigGeneAdult),function(g) {
  repPval = sigGeneAdult[g,grep('pvalue',names(sigGeneAdult))]< 0.05
  repFC = sign(sigGeneAdult[g,grep('log2FoldChange',names(sigGeneAdult))])== 
    sign(dis[g,'Discovery.log2FoldChange'])
  return(all(repPval & repFC))})
sigGeneAdult$numReplicated = sapply(rownames(sigGeneAdult),function(g) {
  repPval = sigGeneAdult[g,grep('pvalue',names(sigGeneAdult))]< 0.05
  repFC = sign(sigGeneAdult[g,grep('log2FoldChange',names(sigGeneAdult))])== 
    sign(dis[g,'Discovery.log2FoldChange'])
  return(sum(repPval & repFC))})
sigGeneAdult = cbind(dis,sigGeneAdult)
sum(sigGeneAdult$replicatedAny,na.rm = T) #82 
sum(sigGeneAdult$replicatedAll,na.rm = T) #28


#########################
# save for plotting sfig3
save(sigGeneP1,sigGeneAdult,file = "tcf4_mouse/rdas/mouse_tcf4_replication_ages.rda")

############################
# replication by mouse model
tt = apply(sigGeneP1[,grep('pvalue',names(sigGeneP1))],2,function(x) sum(x<0.05,na.rm = T))
names(tt) = ss(names(tt),'\\.')
tt
tt =apply(sigGeneAdult[,grep('pvalue',names(sigGeneAdult))],2,function(x) sum(x<0.05,na.rm = T))
names(tt) = ss(names(tt),'\\.')
tt

############################
# genes that replicated in all 
# Mc4r, Atp2b4, Tmem44, Enpp6, Tmem88b, Gm4211, and Mog
sigGeneAdult$Discovery.Symbol[sigGeneAdult$replicatedAll] 

#############################################################
# genes that replicated in any and how many models replicated 
sigGeneAdult[sigGeneAdult$replicatedAny,c('Discovery.Symbol','numReplicated')]
sum(sigGeneAdult$numReplicated>=4,na.rm = T)

###########################
# save replicated gene list!
WriteXLS(sigGeneAdult, ExcelFileName = "tables/stable3_DE_feature_replication.xls")

##########################
# load replication dataset
load('tcf4_mouse/rdas/mouse_tcf4_replication.rda')

##########################
# find stouffer's z score
N = c(35,8,18,17,15,16)
ind = grep('stat',names(sigGeneAdult))

#not including discovery
z_unweighted = rowSums(sigGeneAdult[,ind[-1]])/sqrt(5)
rep = abs(z_unweighted)>2
sum(rep)
#############################
# make tstats into long format
dat = sigGeneP1[,-c(7:11,13,32,33,34)]
dat = melt(dat,id.var = c("Discovery.Symbol"))
dat$label = paste0(dat$Discovery.Symbol,".",ss(as.character(dat$variable),"\\."))
dat = dat[dat$Discovery.Symbol!='',]
dat$type = ss(as.character(dat$variable),"\\.",2)
dat = dcast(dat, label ~ type)
dat$Symbol = ss(dat$label,"\\.")
dat$Line = ss(dat$label,"\\.",2)
dat$Symbol = factor(dat$Symbol,levels = unique(sigGeneP1$Discovery.Symbol[order(rowMeans(sigGeneP1[,grep("stat",names(sigGeneP1))],na.rm=T))]))
labs =  c("Discovery","Act","Nest", "Del", "R579W",'CA1')
dat$Line = droplevels(factor(dat$Line,levels =labs))

###################################
# plot tstats from each mouse model
pdf(file='tcf4_mouse/plots/replication_tstats_heatmap_p1.pdf', width=6,height = 9)
# range(dat$stat,na.rm=T)
theSeq = seq(-8,8,by=0.01) 
my.col <- colorRampPalette(brewer.pal(11,"PiYG"))(length(theSeq))
print(levelplot(stat ~ Symbol+ Line , 
                data= dat, at = theSeq,pretty=T,column.values=NA,
                col.regions = my.col, scales=list(x=list(rot=60, cex=1),y = list(cex=1.15)),
                ylab = "", xlab = ""))
dev.off()

#############################
# make tstats into long format
dat = sigGeneAdult[,-c(7:11,13,44,45,46)]
dat = melt(dat,id.var = c("Discovery.Symbol"))
dat$label = paste0(dat$Discovery.Symbol,".",ss(as.character(dat$variable),"\\."))
dat = dat[dat$Discovery.Symbol!='',]
dat$type = ss(as.character(dat$variable),"\\.",2)
dat = dcast(dat, label ~ type)
dat$Symbol = ss(dat$label,"\\.")
dat$Line = ss(dat$label,"\\.",2)
dat$Symbol = factor(dat$Symbol,levels = unique(sigGeneAdult$Discovery.Symbol[order(rowMeans(sigGeneAdult[,grep("stat",names(sigGeneAdult))],na.rm=T))]))
labs =  c("Discovery","Act","Nest", "Del", "R579W",'CA1')
dat$Line = droplevels(factor(dat$Line,levels =labs))

###################################
# plot tstats from each mouse model
pdf(file='tcf4_mouse/plots/replication_tstats_heatmap_adult.pdf', width=6,height = 9)
# range(dat$stat,na.rm=T)
theSeq = seq(-13,13,by=0.01) 
my.col <- colorRampPalette(brewer.pal(11,"PiYG"))(length(theSeq))
print(levelplot(stat ~ Symbol+ Line , 
                data= dat, at = theSeq,pretty=T,column.values=NA,
                col.regions = my.col, scales=list(x=list(rot=60, cex=.7),y = list(cex=1.15)),
                ylab = "", xlab = ""))
dev.off()