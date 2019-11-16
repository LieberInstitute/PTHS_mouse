## pool Maher, Philpot, and Sweatt RNAseq datasets together
library(jaffelab)

#############################
# load all the phenotype data
load('./rdas/pheno.rda',envir = maherDat<-new.env())
load('../philpot/rdas/pheno.rda',envir = philpotDat<-new.env())

##########################
# clean up phenotype data
maherDat$pd$Line = 'Maher'
philpotDat$pd$FileID = ss(philpotDat$pd$'BGI file name','.fq.gz')
philpotDat$pd$Genotype = ifelse(grepl('Mutant',philpotDat$pd$'General genotype'),'HT', 'WT')

#########
# combine
cols = c('FileID', 'Line', 'Genotype','Age')
pd = rbind(maherDat$pd[,cols],philpotDat$pd[,cols])

##################
# load gene counts
load('/dcl01/lieber/ajaffe/Brady/mouseRNAseq/rawCounts_tcf4_mouse_OCT20_n36.rda', envir = maherCounts<-new.env())
load('/dcl01/lieber/ajaffe/Brady/philpot/rawCounts_philpot_OCT20_n58.rda', envir = philpotCounts<-new.env())

####################
# add processed data
cols = Reduce(intersect,lapply(list(maherCounts$pd,philpotCounts$pd),names))
pd = cbind(pd,rbind(maherCounts$pd[,cols],philpotCounts$pd[,cols]))

####################
# add the FASTQ paths
searchDirs = c('/dcl01/lieber/ajaffe/Brady/mouseRNAseq/FASTQ2','/dcl01/lieber/ajaffe/Brady/philpot/FASTQ')

fqFiles = unlist(sapply(paste0(searchDirs ), list.files, pattern = '.gz', full.names = T, recursive = T))
names(fqFiles) = NULL

R1 = fqFiles[grep('R1|_1.fq', fqFiles)]
R2 = fqFiles[grep('R2', fqFiles)]

pd$R1 = sapply(pd$FileID,  grep, R1,value = T)
tmp = unlist(sapply(pd$FileID,  grep, R2,value = T))
pd$R2 = NULL
pd[names(tmp), 'R2'] = tmp
  
write.csv(pd, 'tables/maher_philpot_PTHS_mouse_rnaseq_pheno.csv', quote = F, row.names = F)
