### huntington human brain differential expression analysis 
# qsub -V -l mf=200,h_vmem=250G,h_stack=256M -cwd -b y R CMD BATCH analyze_PTHS_human_DESeq2.R
source('/users/bphan/tcf4/PTHS_mouse/DESeq2_functions.R') #work-horse of differential expression
source("/users/ajaffe/Lieber/lieber_functions_aj.R")
library(sva)

##############################################
# load phenotype data and RPKM expression data
load('./rdas/pheno.rda',envir = dat<-new.env())
load('/dcl01/lieber/ajaffe/Brady/huntington/rawCounts_huntington_nov28_n69.rda')
pd = cbind(dat$pd,pd)
pd$AgeGroup = ifelse(pd$Age <45, 'A',
             ifelse(pd$Age < 60, 'B',
             ifelse(pd$Age< 75,'C','D')))
table(pd$AgeGroup, pd$Case)
all(pd$"SampleID" == pd$SAMPLE_ID) #true
rownames(pd) = pd$SAMPLE_ID

########################################
# get degradation-prone region coverages 
if(FALSE){
fn = paste0("/dcl01/lieber/ajaffe/Brady/huntington/degradation/",pd$Run_s, "_degradeStats_polyA.txt")
table(file.exists(fn))
degCov = sapply(fn, function(x) {
  cat(".")
  read.delim(pipe(paste("cut -f10", x)), as.is=TRUE)$sum
})
names(degCov) = pd$Run_s
#library size and read-length normalized, 101 bp
degCovAdj = degCov/101/matrix(pd$totalMapped/80e6,nrow = nrow(degCov),
                              ncol = nrow(pd),byrow = T)
save(degCovAdj, file="rdas/degradation_list_huntington.rda")
} else {
  load('rdas/degradation_list_huntington.rda')
}

##########
# get qSVs
qSVA = qsva(degCovAdj,model.matrix(~pd$Case))
n.sv = ncol(qSVA)
pd = cbind.data.frame(pd,qSVA)

t.test(pd$totalAssignedGene~pd$Case)

##############################
# create and run DESeq objects
formula = as.formula(paste0('~Case+AgeGroup+totalAssignedGene+PC',paste0(seq(n.sv),collapse ='+PC')))
geneDds <- DESeq2(countData = geneCounts, colData = pd, sva = FALSE, parallel=TRUE,design = formula)

############################################
# get DE results, and fold-change PTHS v. WT
resGene <- results(geneDds,contrast = c('Case','CTL','HNT'),alpha=0.05) 
sum(resGene$padj < 0.05, na.rm=TRUE)

outGene <- as.data.frame(resGene[order(resGene$padj,resGene$pvalue),])
outGene = cbind(outGene,geneMap[rownames(outGene),])
sigGene = outGene[which(outGene$padj<0.05 &!is.na(outGene$padj)),]

pdf('plots/DESeq2_MA_plots_huntington.pdf')
plotMA(resGene, main="Gene MA plot", ylim=c(-2,2))
dev.off()

#################
# save everything
library(WriteXLS)
WriteXLS(list(Gene = sigGene,pheno = pd), ExcelFileName = 'tables/huntington_DE_table_DESeq2.xls',row.names=T)
save(outGene,file = 'rdas/huntington_DE_objects_DESeq2.rda')
save(geneDds, file = '/dcl01/lieber/ajaffe/Brady/huntington/huntington_DESeq2.rda')


