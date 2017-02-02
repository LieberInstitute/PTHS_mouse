library(pSI)
library(parallel)
library(WriteXLS)
library(ggplot2)
library(jaffelab)
library("scales")
reverselog_trans <- function(base = exp(1)) {
  trans <- function(x) -log(x, base)
  inv <- function(x) base^(-x)
  trans_new(paste0("reverselog-", format(base)), trans, inv, 
            log_breaks(base = base), 
            domain = c(1e-100, Inf))
}

load('rdas/csea_human_mouse_pSI.rda')
load("rdas/mega_tcf4_ages_DE_objects_DESeq2.rda")

# genes for enrichment
gList = lapply(outGeneList, function(x) toupper(x$Symbol[x$padj< 0.05 & !is.na(x$padj)]))
sapply(gList,length)

enrList = mclapply(gList, fisher.iteration, pSIs = mouse$psi.out,
                 background="data.set",mc.cores = 2)

sigEnrList = lapply(enrList, function(x) x[rowSums(x < 0.05) > 0,])
sigEnrList = lapply(sigEnrList,function(x) cbind(labs[rownames(x),],x))

WriteXLS(sigEnrList,ExcelFileName ="tables/CSEA_enrichment_mega_tcf4_all_ages.xlsx",
         row.names= TRUE)

############################
# make data into long format
dat = do.call('rbind',sigEnrList)
dat = dat[sapply(rownames(dat),stringr::str_count,pattern = '\\.')==2,]
dat$Age = factor(ss(as.character(rownames(dat)),'\\.'),levels = c('p1','Adult'))
datLong = reshape2::melt(dat,id.vars = c('Description',"Age",'Region','Cell','Region'),value.name = 'padj',variable.name = 'pSI')
datLong$pSI = as.numeric(ss(as.character(datLong$pSI),'-'))
datLong$pSI[datLong$pSI==1] = 1e-4
datLong = datLong[datLong$padj<0.05,]


############################
# plot CSEA for p1 and adult
ggplot(data = datLong,aes(y = padj,x = Age,colour = pSI,shape = Age,size = .5))+
  geom_point(position=position_jitterdodge(dodge.width = .5))+
  facet_wrap(Region~Cell,scales = 'free')+xlab('Region')+ylab('Adjusted P-value')+
  geom_hline(yintercept =0.05,colour = 'red',linetype = 2)+
  scale_y_continuous(trans = reverselog_trans(base = 10))+
  scale_fill_continuous(low="black", high="blue", limits=c(0,1))



