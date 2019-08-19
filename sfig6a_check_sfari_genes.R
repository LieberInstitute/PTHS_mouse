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
dim(outGene) #16123 background genes
sigGene = with(outGene,outGene[which(padj < 0.05 &!is.na(padj)),])
dim(sigGene) #1711 DEGs with human homolog

########################
# load SFARI human genes
humanSFARI = read.csv('tcf4_mouse/tables/SFARI-Gene_genes_08-06-2019release_08-08-2019export.csv')
humanSFARI = with(humanSFARI, humanSFARI[order(gene.score),])
humanSFARI = cbind(humanSFARI,outGene[lookfor(humanSFARI$gene.symbol,outGene$Symbol),])
humanSFARI= humanSFARI[!is.na(humanSFARI$Symbol),]
humanSFARI = humanSFARI[!duplicated(humanSFARI),]
nrow(humanSFARI) # 992 genes used

#########################
# load SFARI mouse models
mouseSFARI = read.csv('tcf4_mouse/tables/SFARI-Gene_animal-genes_08-06-2019release_08-08-2019export.csv')
mouseSFARI = with(mouseSFARI,mouseSFARI[model.species=='Mus musculus',])
mouseSFARI = cbind(mouseSFARI,outGene[lookfor(mouseSFARI$gene.symbol,outGene$Symbol),])
mouseSFARI= mouseSFARI[!is.na(mouseSFARI$Symbol),]
mouseSFARI = mouseSFARI[order(match(mouseSFARI$Symbol,humanSFARI$Symbol)),]
mouseSFARI = mouseSFARI[!duplicated(mouseSFARI),]
nrow(mouseSFARI) # 255 expressed in mouse tcf4 dataset

###########################
# list of DEG in mouse SFARI
# 72 genes differentially expressed and in SFARI 
(t1 = with(outGene,table(inMouseSFARI = Symbol %in% mouseSFARI$Symbol,inDEG = padj < 0.05 & !is.na(padj))))
fisher.test(t1) # OR = 3.42, p-value = 2.874e-15
ind1 = with(sigGene, which(mouseSFARI$Symbol %in% Symbol))

#######################################
# list of DEG in scored human SFARI list
# 190 human SFARI ASD-linked genes in our list of DEGs
(t2 = with(outGene,table(inHumanSFARI = Symbol %in% humanSFARI$Symbol,inDEG = padj < 0.05 & !is.na(padj))))
fisher.test(t2) # OR = 2.116, p-value = 2.2e-16
ind2 = with(sigGene, which(humanSFARI$Symbol %in% Symbol))

##########################
## sensitivity for gene score
humanSFARI_strict = humanSFARI[humanSFARI$gene.score %in% 1:3 | 
								humanSFARI$syndromic == 1,] 
(t2_strict = with(outGene,table(inHumanSFARI = Symbol %in% humanSFARI_strict$Symbol,inDEG = padj < 0.05 & !is.na(padj))))
fisher.test(t2_strict) # OR = 2.95988, p-value = 2.2e-16
				
####################
# save the two lists
sfariList = list(Scored_Human_SFARI_Genes = humanSFARI[ind2,],
                 Mouse_Model_Genes = mouseSFARI[ind1,])
WriteXLS(sfariList,ExcelFileName = 'tcf4_mouse/tables/mega_tcf4_adult_SFARI_asd_gene_list_20190814.xls')




####################
# Sander et al, Neuron 2015 enrichments
sander = readxl::read_xlsx('tcf4_mouse/tables/Sanders_et_al-Neuron_2015-Table-S6_20190814.xlsx', sheet = 'TADA_Scores')
sander = sander[sander$tadaFdrAscSscExomeSscAgpSmallDel <0.1, ] #select the 65 small del genes FDR <0.1
sander$RefSeqGeneName = snakecase::to_upper_camel_case(sander$RefSeqGeneName)
sander = cbind(sander,outGene[lookfor(sander$RefSeqGeneName,outGene$Symbol),])
sander= sander[!is.na(sander$Symbol),]
sander = sander[!duplicated(sander),]

nrow(sander) # 61 genes used
sum(sander$Symbol %in% outGene$Symbol)

(t3 = with(outGene,table(inHumanSFARI = Symbol %in% sander$Symbol,inDEG = padj < 0.05 & !is.na(padj))))
fisher.test(t3) # OR = 5.92, p-value = 7.332e-10
ind3 = with(sigGene, which(sander$Symbol %in% Symbol))

####################
# Satterstrom et al, bioRxiv 2018 enrichments
satterstrom = readxl::read_xlsx('tcf4_mouse/tables/Satterstrom_et_al-bioRxiv_Table-S4_20190814.xlsx', sheet = '102 ASD FDR≤0.1 genes')
satterstrom$hugoGene = snakecase::to_upper_camel_case(satterstrom$hugoGene)
satterstrom = cbind(satterstrom,outGene[lookfor(satterstrom$hugoGene,outGene$Symbol),])
satterstrom= satterstrom[!is.na(satterstrom$Symbol),]
satterstrom = satterstrom[!duplicated(satterstrom),]

nrow(satterstrom) # 97 genes used
sum(satterstrom$Symbol %in% outGene$Symbol)

(t4 = with(outGene,table(inHumanSFARI = Symbol %in% satterstrom$Symbol,inDEG = padj < 0.05 & !is.na(padj))))
fisher.test(t4) # OR = 4.83, p-value = 2.325e-11
ind4 = with(sigGene, which(satterstrom$hugoGene %in% Symbol))




### venn diagram
library(limma)
vennMat = matrix(FALSE, nrow = nrow(outGene), nc = 3)
colnames(vennMat) = c("Tcf4", "SFARI_Hs", "SFARI_Mm")
rownames(vennMat) = rownames(outGene)
vennMat[rownames(vennMat) %in% rownames(sigGene),1] = TRUE 
vennMat[outGene$hsapien_homolog %in% humanSFARI$hsapien_homolog,2] = TRUE 
vennMat[outGene$Symbol %in% mouseSFARI$Symbol,3] = TRUE 
pdf("./tcf4_mouse/plots/vennDiagram_tcf4_20190814.pdf")
vennDiagram(vennCounts(vennMat),cex=1.5)
dev.off()




