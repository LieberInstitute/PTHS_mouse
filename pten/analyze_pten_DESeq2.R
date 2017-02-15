### PTEN mouse differential expression analysis using DESeq2
# qsub -V -l mf=100G,h_vmem=150G,h_stack=256M -cwd -b y R CMD BATCH analyze_pten_DESeq2.R
source('../DESeq2_functions.R') #work-horse of differential expression

##############################################
# load phenotype data and RPKM expression data
load('./rdas/pheno.rda',envir = dat<-new.env())
load('/dcl01/lieber/ajaffe/Brady/mouseRNAseq/TilotPTEN/rawCounts_pten_OCT20_n18.rda')
pd = cbind(dat$pd,pd)
all.equal(pd$FileID,pd$SAMPLE_ID) #samples line up

##########################
# conform to DESeq2 data types
jCounts = as.matrix(as.data.frame(jCounts))
jIndex=which(jMap$code != "Novel")
jCounts = jCounts[jIndex,]
jMap = jMap[jIndex]
colnames(jCounts) = pd$SAMPLE_ID
rownames(pd) = pd$SAMPLE_ID

##############################
# create and run DESeq objects
geneDds <- DESeq2(countData = geneCounts, colData = pd, design = ~Genotype+Age,sva = TRUE,parallel=TRUE)
rm(geneCounts); gc()
exonDds <- DESeq2(countData = exonCounts, colData = pd, design = ~Genotype+Age,sva = TRUE,parallel=TRUE)
rm(exonCounts); gc()
jxnDds <- DESeq2(countData = jCounts, colData = pd, design = ~Genotype+Age,sva = TRUE,parallel=TRUE)
rm(jCounts); gc()

################################################################
# get DE results, and fold-change homozygous mutant v. wild-type
resGene <- results(geneDds,contrast = c('Genotype','Pten m3m4/m3m4','Pten wt/wt'),alpha=0.05) 
resExon <- results(exonDds,contrast = c('Genotype','Pten m3m4/m3m4','Pten wt/wt'),alpha=0.05) 
resJxn <- results(jxnDds,contrast = c('Genotype','Pten m3m4/m3m4','Pten wt/wt'),alpha=0.05) 

sum(resGene$padj < 0.05, na.rm=TRUE)
sum(resExon$padj < 0.05, na.rm=TRUE)
sum(resJxn$padj < 0.05, na.rm=TRUE)

outGene <- as.data.frame(resGene[order(resGene$padj,resGene$pvalue),])
outGene = cbind(outGene,geneMap[rownames(outGene),])
sigGene = outGene[which(outGene$padj<.05),]

outExon <- as.data.frame(resExon[order(resExon$padj,resExon$pvalue),])
outExon = cbind(outExon,exonMap[rownames(outExon),])
sigExon = outExon[which(outExon$padj<.05),]

outJxn <- as.data.frame(resJxn[order(resJxn$padj,resJxn$pvalue),])
outJxn = cbind(outJxn, as.data.frame(jMap)[rownames(outJxn),])
sigJxn = outJxn[which(outJxn$padj<.05),]

pdf('plots/DESeq2_MA_plots_pten.pdf')
plotMA(resGene, main="Gene MA plot", ylim=c(-.75,.75))
plotMA(resExon, main="Exon MA plot", ylim=c(-.75,.75))
plotMA(resJxn, main="Junction MA plot", ylim=c(-.75,.75))
dev.off()

#################
# save everything
library(WriteXLS)
WriteXLS(list(Gene = sigGene,Exon = sigExon,Junction = sigJxn), ExcelFileName = 'tables/pten_DE_table_DESeq2.xls',row.names=T)
save(outGene,outExon,outJxn,file = 'rdas/pten_DE_objects_DESeq2.rda')
save(geneDds,exonDds,jxnDds, file = '/dcl01/lieber/ajaffe/Brady/mouseRNAseq/TilotPTEN/pten_DESeq2_svaAdj.rda')

#################################
# format data to run on CIBERSORT
write.table(geneCounts,sep = '\t',quote = FALSE,file = 'tables/rawGeneCounts_pten.txt')

