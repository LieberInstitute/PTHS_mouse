##########################
# CSEA for Maher tcf4 mice
library(pSI)
library(WriteXLS)
if(FALSE){
  source("/users/ajaffe/Lieber/Projects/KeriM/csea_functions.R") 
  labs = read.delim('tables/CSEA_mouse_cell_types.csv',stringsAsFactors = F)
  rownames(labs) = toupper(labs$Label); labs[,1] = NULL
  mouse = lapply(mouse,function(x) {
    names(x) = toupper(names(x)); 
    rownames(x) = toupper(rownames(x)); 
    x})
  human = lapply(mouse,function(x) {
    names(x) = toupper(names(x)); 
    rownames(x) = toupper(rownames(x)); 
    x})
  save(labs,mouse,human,file = 'rdas/csea_human_mouse_pSI.rda')
} else{
  load('rdas/csea_human_mouse_pSI.rda')
}
load('rdas/mouse_tcf4_ages_DE_objects_DESeq2.rda')




# genes for enrichment
gList = lapply(outGeneList, function(x) toupper(x$Symbol[x$padj< 0.05 & !is.na(x$padj)]))
sapply(gList,length)

enrList = lapply(gList, fisher.iteration, pSIs = mouse$psi.out,
                 background="data.set")

sigEnrList = lapply(enrList, function(x) x[rowSums(x < 0.05) > 0,])
sigEnrList = lapply(sigEnrList,function(x) cbind(Description = labs[rownames(x),1],x))

WriteXLS(sigEnrList,ExcelFileName ="tables/CSEA_enrichment_tcf4_mouse_all_ages.xlsx",
         row.names= TRUE)
