## Check mega tcf4 dataset DEGs with human asd and animal model genes from SFARI
lookfor= function(this,inThat) { #gene name matching
  this = toupper(this); inThat = toupper(inThat);
  tmp = sapply(this,function(x) grep(paste0('^',x,'$'),inThat))
  return(sapply(tmp,function(x) ifelse(length(x)==0,NA,x[1])))}

library(jaffelab)
library(biomaRt)
library(WriteXLS)

options(stringsAsFactors = F)


##################################
# get Ensembl mouse to human genes
ensembl = useMart("ensembl",dataset = 'mmusculus_gene_ensembl')
MMtoHG = getBM(attributes = c('ensembl_gene_id','hsapiens_homolog_ensembl_gene'),
               mart = ensembl)

##############################
# load the mega_tcf4_mice data
load('tcf4_mouse/rdas/mega_tcf4_ages_DE_objects_DESeq2.rda')

outGene = outGeneList[['Adult']]
# subset mouse genes to those w/ human homolog
outGene$hsapien_homolog = MMtoHG$hsapiens_homolog_ensembl_gene[match(rownames(outGene),MMtoHG$ensembl_gene_id)]
# keep only expressed genes (DESeq2 will not adjust p-value if not expressed enough)
#outGene = outGene[!is.na(outGene$padj),] #take only genes w/ human homolog
outGene = outGene[grepl('ENSG',outGene$hsapien_homolog),] #take only genes w/ human homolog
dim(outGene) #16303 background genes
sigGene = with(outGene,outGene[which(padj < 0.05 &!is.na(padj)),])
dim(sigGene) #1725 DEGs with human homolog

########################
# load SFARI human genes
humanSFARI = read.csv('tcf4_mouse/tables/SFARI-Gene_genes_export30-11-2017.csv')
humanSFARI = with(humanSFARI, humanSFARI[order(gene.score),])
humanSFARI = cbind(humanSFARI,outGene[lookfor(humanSFARI$gene.symbol,outGene$Symbol),])
humanSFARI= humanSFARI[!is.na(humanSFARI$Symbol),]
humanSFARI = humanSFARI[!duplicated(humanSFARI),]
nrow(humanSFARI) # 886 genes used

#########################
# load SFARI mouse models
mouseSFARI = read.csv('tcf4_mouse/tables/SFARI-Gene_animal-genes_export30-11-2017.csv')
mouseSFARI = with(mouseSFARI,mouseSFARI[model.species=='Mus musculus',])
mouseSFARI = cbind(mouseSFARI,outGene[lookfor(mouseSFARI$gene.symbol,outGene$Symbol),])
mouseSFARI= mouseSFARI[!is.na(mouseSFARI$Symbol),]
mouseSFARI = mouseSFARI[order(match(mouseSFARI$Symbol,humanSFARI$Symbol)),]
mouseSFARI = mouseSFARI[!duplicated(mouseSFARI),]
nrow(mouseSFARI) # 234 expressed in mouse tcf4 dataset

###########################
# list of DEG in mouse SFARI
# 67 genes differentially expressed and in SFARI 
(t1 = with(outGene,table(inMouseSFARI = Symbol %in% mouseSFARI$Symbol,inDEG = padj < 0.05)))
fisher.test(t1) # OR = 3.141425, p-value = 1.231e-12
ind1 = with(sigGene, which(mouseSFARI$Symbol %in% Symbol))

#######################################
# list of DEG in scored human SFARI list
# 170 human SFARI ASD-linked genes in our list of DEGs
(t2 = with(outGene,table(inHumanSFARI = Symbol %in% humanSFARI$Symbol,inDEG = padj < 0.05)))
fisher.test(t2) # OR = 1.899831, p-value = 1.261e-11
ind2 = with(sigGene, which(humanSFARI$Symbol %in% Symbol))

####################
# save the two lists
sfariList = list(Scored_Human_SFARI_Genes = humanSFARI[ind2,],
                 Mouse_Model_Genes = mouseSFARI[ind1,])
WriteXLS(sfariList,ExcelFileName = 'tcf4_mouse/tables/mega_tcf4_adult_SFARI_asd_gene_list_20173011.xls')

### venn diagram
library(limma)
vennMat = matrix(FALSE, nrow = nrow(outGene), nc = 3)
colnames(vennMat) = c("Tcf4", "SFARI_Hs", "SFARI_Mm")
rownames(vennMat) = rownames(outGene)
vennMat[rownames(vennMat) %in% rownames(sigGene),1] = TRUE 
vennMat[outGene$hsapien_homolog %in% humanSFARI$hsapien_homolog,2] = TRUE 
vennMat[outGene$Symbol %in% mouseSFARI$Symbol,3] = TRUE 
pdf("./tcf4_mouse/plots/vennDiagram_tcf4_20171130.pdf")
vennDiagram(vennCounts(vennMat),cex=1.5)
dev.off()