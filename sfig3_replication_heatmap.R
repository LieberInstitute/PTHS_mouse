# plot replicated genes
library(jaffelab)
library(lattice)
library(RColorBrewer)
library(reshape2)

##########################
# load replication dataset
load('tcf4_mouse/rdas/mouse_tcf4_replication.rda')

##########################
# find stouffer's z score
N = c(35,8,18,17,15,16)
ind = grep('stat',names(sigGene))

#not including discovery
z_unweighted = rowSums(sigGene[,ind[-1]])/sqrt(5)
rep = abs(z_unweighted)>2
sum(rep)

#############################
# make tstats into long format
dat = sigGene[,-c(7:11,13)]
dat = melt(dat,id.var = c("Symbol_Discovery"))
dat$label = paste0(dat$Symbol,"_",ss(as.character(dat$variable),"_",2))
dat$type = ss(as.character(dat$variable),"_")
dat = dcast(dat, label ~ type)
dat$Symbol = ss(dat$label,"_")
dat$Line = ss(dat$label,"_",2)
dat$Symbol = factor(dat$Symbol,levels = c(sigGene$Symbol_Discovery[order(sigGene$stat_Discovery)][-2],'Tcf4'))
labs =  c("Discovery","Act","Nest", "Del", "R579W",'CA1')
dat$Line = factor(dat$Line,levels =labs)
dat$Line2 = factor(paste0(dat$Line,' (N=',N[as.numeric(dat$Line)],')'),levels= paste0(labs,' (N=',N,')'))

###################################
# plot tstats from each mouse model
pdf(file='tcf4_mouse/plots/sfig3_replication_tstats_heatmap.pdf', width=6,height = 9)
# range(dat$stat)
theSeq = seq(-16,16,by=0.01) 
my.col <- colorRampPalette(brewer.pal(11,"PiYG"))(length(theSeq))
print(levelplot(stat ~ Line2 + Symbol, 
                data= dat, at = theSeq,pretty=TRUE,
                col.regions = my.col, scales=list(y=list(cex=1.15), 
                                                  x=list(rot=30, cex=1.15)),
                ylab = "", xlab = ""))
dev.off()
