## pool Maher, Philpot, and Sweatt RNAseq datasets together
library(jaffelab)

#############################
# load all the phenotype data
load('./rdas/pheno.rda',envir = maherDat<-new.env())
load('../philpot/rdas/pheno.rda',envir = philpotDat<-new.env())
load('../sweatt/rdas/pheno.rda',envir = sweattDat<-new.env())

##########################
# clean up phenotype data
maherDat$pd$Line = 'Maher'
sweattDat$pd$Line = 'Sweatt'
sweattDat$pd$Age = 'Adult'
philpotDat$pd$FileID = ss(philpotDat$pd$'BGI file name','.fq.gz')
philpotDat$pd$Genotype = ifelse(grepl('Mutant',philpotDat$pd$'General genotype'),'HT', 'WT')

#########
# combine
cols = c('FileID', 'Line', 'Genotype','Age')
pd = rbind(maherDat$pd[,cols],philpotDat$pd[,cols],sweattDat$pd[,cols])

##################
# load gene counts
load('/dcl01/lieber/ajaffe/Brady/mouseRNAseq/rawCounts_tcf4_mouse_OCT20_n36.rda', envir = maherCounts<-new.env())
load('/dcl01/lieber/ajaffe/Brady/philpot/rawCounts_philpot_OCT20_n58.rda', envir = philpotCounts<-new.env())
load('/dcl01/lieber/ajaffe/Brady/sweatt/rawCounts_sweatt_OCT20_n16.rda', envir = sweattCounts<-new.env())

####################
# add processed data
cols = Reduce(intersect,lapply(list(maherCounts$pd,philpotCounts$pd,sweattCounts$pd),names))
pd = cbind(pd,rbind(maherCounts$pd[,cols],philpotCounts$pd[,cols],sweattCounts$pd[,cols]))
#all.equal(pd$FileID,pd$SAMPLE_ID) #TRUE

###############################################
# intersect expressed gene/exon/junction counts
mGene = Reduce(intersect,list(rownames(maherCounts$geneMap),
                              rownames(philpotCounts$geneMap),
                              rownames(sweattCounts$geneMap)))
mExon = Reduce(intersect,list(with(maherCounts$exonMap,paste(Chr,Start,End))
                              ,with(philpotCounts$exonMap,paste(Chr,Start,End)),
                              with(sweattCounts$exonMap,paste(Chr,Start,End))))
mJxn = Reduce(intersect,list(names(maherCounts$jMap),
                             names(philpotCounts$jMap),
                             names(sweattCounts$jMap)))

###################
# merge annotations
geneMap = maherCounts$geneMap[mGene,]
indExon = match(mExon,with(maherCounts$exonMap,paste(Chr,Start,End)))
exonMap = maherCounts$exonMap[indExon,]
jMap = maherCounts$jMap[mJxn,]

#############
# merge counts
indExon1 = match(mExon,with(maherCounts$exonMap,paste(Chr,Start,End)))
indExon2 = match(mExon,with(philpotCounts$exonMap,paste(Chr,Start,End)))
indExon3 = match(mExon,with(sweattCounts$exonMap,paste(Chr,Start,End)))

geneCounts = cbind(maherCounts$geneCounts[mGene,],
                   philpotCounts$geneCounts[mGene,],
                   sweattCounts$geneCounts[mGene,])

exonCounts = cbind(maherCounts$exonCounts[indExon1,],
                philpotCounts$exonCounts[indExon2,],
                sweattCounts$exonCounts[indExon3,])

jCounts = cbind(maherCounts$jCounts[mJxn,],
                   philpotCounts$jCounts[mJxn,],
                   sweattCounts$jCounts[mJxn,])

save(pd, geneCounts,geneMap,exonCounts,exonMap,jCounts,jMap,file = '/dcl01/lieber/ajaffe/Brady/mouseRNAseq/rawCounts_mega_dataset_nov14_n110.rda')