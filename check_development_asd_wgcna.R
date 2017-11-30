## check with Parikshak human development WGCNA data and Voineagu human ASD WGCNA data
getOR = function(x) x[1,1]/x[2,1]/x[1,2]*x[2,2]
library(WriteXLS)
library(biomaRt)
library(lattice)
library(reshape2)
library(jaffelab)
library(readxl)

################################
#get Ensembl mouse to human genes
ensembl=useMart("ensembl",dataset = 'mmusculus_gene_ensembl')
MMtoHG = getBM(attributes = c('ensembl_gene_id','hsapiens_homolog_ensembl_gene'),mart = ensembl)
hgMart = getBM(attributes = c('ensembl_gene_id','illumina_humanht_12_v3','entrezgene'),
               mart = useMart("ensembl",dataset = 'hsapiens_gene_ensembl'))

###################################
# load Parikshak supplementary data
dat = data.frame(read_excel('tcf4_mouse/tables/mmc1.xlsx',skip = 2),stringsAsFactors = F)
rownames(dat) = dat$"ENSEMBL.GENE.ID" 

############################
# load mega TCF4 DEG dataset
load('tcf4_mouse/rdas/mega_tcf4_ages_DE_objects_DESeq2.rda')
outGene = outGeneList[['Adult']]

outGene$hsapien_homolog = MMtoHG$hsapiens_homolog_ensembl_gene[match(rownames(outGene),MMtoHG$ensembl_gene_id)]
outGene = outGene[grepl('ENSG',outGene$hsapien_homolog),] #take only genes w/ human homolog
outGene$module = dat[outGene$hsapien_homolog,c('Module.Label')]
outGene = outGene[complete.cases(outGene),]

#########################################
## module enrichment of genes padj < 0.05
modules = unique(outGene$module)
modules = modules[modules != '-']
modules = modules[order(as.numeric(ss(modules,'M',2)))]
pval = sapply(modules,function(m) fisher.test(table(outGene$padj<.05,outGene$module==m))$p.value)
OR = sapply(modules,function(m) fisher.test(table(outGene$padj<.05,outGene$module==m))$estimate)
adj.pval = p.adjust(pval,'fdr')
N =  sapply(modules, function(m) sum(outGene$padj<.05& outGene$module==m))
N_max = sapply(modules, function(m) sum(outGene$module==m))
sum(adj.pval<0.05)

############################
# hold onto enriched modules
modules[adj.pval<.05 & OR >1] #"M2" "M13"
tmp1 = data.frame(logPval = -log10(adj.pval),OR = OR, Feat = 'Parikshak2013',
                  Module = paste0('dev',modules), N = N, N_max = N_max)

#################################################
# load Voineagu huma ASD WGCNA supplementary data
dat = read.csv('tcf4_mouse/tables/nature10110-s3.csv',stringsAsFactors = F)
dat = dat[order(as.numeric(ss(dat$Probe.ID,'_',2))),]
dat$X = NULL
# assign module membership by which.max and if membership greater than 0.7, from Voineagu paper
dat$Module.Label = apply(dat[,3:20],1,function(x) ifelse(any(x>.7),paste0('M',which.max(x)),'-'))

#################################################
# load microarray probes to link to ENSEMBL genes
probes = read.delim('tcf4_mouse/tables/GPL10558-50081.txt',skip = 30,stringsAsFactors = F)
probes = probes[!duplicated(probes$ID) & !is.na(ss(probes$ID,'_',2)),]
rownames(probes) = probes$ID

dat$EntrezID = probes[dat$Probe.ID,'Entrez_Gene_ID']
dat$hg_ensembl = hgMart[match(dat$Probe.ID,hgMart$illumina_humanht_12_v3),'ensembl_gene_id']
dat = dat[!is.na(dat$hg_ensembl) & !duplicated(dat$hg_ensembl),]
rownames(dat) = dat$hg_ensembl

##################################
# reload mouse data to start fresh
load('tcf4_mouse/rdas/mega_tcf4_ages_DE_objects_DESeq2.rda')
outGene = outGeneList[['Adult']]
outGene$hsapien_homolog = MMtoHG$hsapiens_homolog_ensembl_gene[match(rownames(outGene),MMtoHG$ensembl_gene_id)]
outGene = outGene[grepl('ENSG',outGene$hsapien_homolog),] #take only genes w/ human homolog
outGene$module = dat[outGene$hsapien_homolog,c('Module.Label')]
outGene = outGene[complete.cases(outGene),]

####################
## module enrichment
modules = unique(outGene$module)
modules = modules[modules != '-']
modules = modules[order(as.numeric(ss(modules,'M',2)))]
pval = sapply(modules,function(m) fisher.test(table(outGene$padj<.05,outGene$module==m))$p.value)
OR = sapply(modules,function(m) fisher.test(table(outGene$padj<.05,outGene$module==m))$estimate)
adj.pval = p.adjust(pval,'fdr')
N =  sapply(modules, function(m) sum(outGene$padj<.05& outGene$module==m))
N_max = sapply(modules, function(m) sum(outGene$module==m))

####################
# save Voineagu data
tmp2 = data.frame(logPval = -log10(adj.pval),OR = OR, Feat = 'Voineagu',
                  Module = paste0('asd',modules), N = N, N_max = N_max)
modules[adj.pval<.05 & OR >1] #"M5"  "M14"

#############################################
# load new Parikshak 2016 supplementary data
dat = data.frame(read_excel('tcf4_mouse/tables/077057-3.xlsx',skip = 2,sheet = 'TableS2a'),stringsAsFactors = F)
rownames(dat) = dat$"EnsemblID" 

##################################
# reload mouse data to start fresh
load('tcf4_mouse/rdas/mega_tcf4_ages_DE_objects_DESeq2.rda')
outGene = outGeneList[['Adult']]
outGene$hsapien_homolog = MMtoHG$hsapiens_homolog_ensembl_gene[match(rownames(outGene),MMtoHG$ensembl_gene_id)]
outGene = outGene[grepl('ENSG',outGene$hsapien_homolog),] #take only genes w/ human homolog
outGene$module = dat[outGene$hsapien_homolog,c('ModuleLabel')]
outGene = outGene[complete.cases(outGene),]
outGene$module = paste0('M',outGene$module)

####################
## module enrichment
modules = unique(outGene$module)
modules = modules[modules != '-']
modules = modules[order(as.numeric(ss(modules,'M',2)))]
modules = modules[modules != 'M11'] #non coexpressed module
pval = sapply(modules,function(m) fisher.test(table(outGene$padj<.05,outGene$module==m))$p.value)
OR = sapply(modules,function(m) fisher.test(table(outGene$padj<.05,outGene$module==m))$estimate)
adj.pval = p.adjust(pval,'fdr')
N =  sapply(modules, function(m) sum(outGene$padj<.05 & outGene$module==m))
N_max = sapply(modules, function(m) sum(outGene$module==m))

####################
# save Parikshak 2016 data
tmp3 = data.frame(logPval = -log10(adj.pval),OR = OR, Feat = 'Parikshak2016',
                  Module = paste0('ctx',modules), N = N, N_max = N_max)
modules[adj.pval<.05 & OR >1] #"M1"  "M3" "M16"

##########################################
# combine module enrichment from both papers
tmp = rbind(tmp2,tmp1,tmp3)
rownames(tmp) = tmp$Module
#module 7 doesn't exist in Parikshak, filler values log10(p) = 0, OR = 1
#tmp = rbind(tmp,data.frame( logPval =0, OR=1,Feat= 'Parikshak', Module = 'M7')) 
tmp$pval = 10^-tmp$logPval
tmp = tmp[order(tmp$pval),]
tmp[with(tmp,which(pval < 0.05 & OR >1 )),]
signif(tmp$pval[with(tmp,which(pval < 0.05 & OR >1 ))],3)


###############################
# plot the two modules together
pdf("tcf4_mouse/plots/fig6a_wgcna_module_enrichment.pdf",w=9,h =2.5)
theSeq = seq(0,20,by=0.001) 
my.col <- colorRampPalette(c("white","darkblue"))(length(theSeq))
print(levelplot(logPval ~ Module + Feat,
                data= tmp, at = theSeq,pretty=TRUE,
                col.regions = my.col, scales=list(y=list(cex=1.25), 
                                                  x=list(cex=1)),
                ylab = "", xlab = ""))
dev.off()
