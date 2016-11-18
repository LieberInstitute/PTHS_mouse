## pool Tcf4, Mef2c, MeCP2, Pten, Ube3a datasets together
library(jaffelab)

#############################
# load all the phenotype data
load('./rdas/pheno.rda',envir = maherDat<-new.env())
load('../philpot/rdas/pheno.rda',envir = philpotDat<-new.env())
load('../sweatt/rdas/pheno.rda',envir = sweattDat<-new.env())
load('../mef2c/rdas/pheno.rda',envir = mef2cDat<-new.env())
load('../mecp2/rdas/pheno.rda',envir = mecp2Dat<-new.env())
load('../pten/rdas/pheno.rda',envir = ptenDat<-new.env())
load('../ube3a/rdas/pheno.rda',envir = ube3aDat<-new.env())

ptenIndex = which(ptenDat$pd$Genotype != 'Pten wt/m3m4') # exclude pten het mice

##########################
# clean up phenotype data
maherDat$pd$Line = 'Maher'
sweattDat$pd$Line = 'Sweatt'
mef2cDat$pd$Line = 'Mef2c'
mecp2Dat$pd$Line = 'Mecp2'
ptenDat$pd$Line = 'Pten'
ube3aDat$pd$Line = 'Ube3a'

sweattDat$pd$Age = 'Adult'
mef2cDat$pd$Age = 'Adult'
mecp2Dat$pd$Age = 'Adult'
ptenDat$pd$Age = ifelse(grepl('6-weeks old at sacrifice',ptenDat$pd$Age),'p14','Adult')

mef2cDat$pd$FileID = mef2cDat$pd$SampleID
mecp2Dat$pd$FileID = mecp2Dat$pd$SampleID
philpotDat$pd$FileID = ss(philpotDat$pd$'BGI file name','.fq.gz')
ube3aDat$pd$FileID = ss(ube3aDat$pd$'BGI file name','.fq.gz')

sweattDat$pd$Genotype = ifelse(grepl('HT',sweattDat$pd$Genotype),'MT', 'WT')
maherDat$pd$Genotype = ifelse(grepl('HT',maherDat$pd$Genotype),'MT', 'WT')
philpotDat$pd$Genotype = ifelse(grepl('Mutant',philpotDat$pd$'General genotype'),'MT', 'WT')
ube3aDat$pd$Genotype = ifelse(grepl('Mutant',ube3aDat$pd$'General genotype'),'MT', 'WT')
mef2cDat$pd$Genotype = ifelse(as.numeric(mef2cDat$pd$Genotype)==1,'WT', 'MT')
mecp2Dat$pd$Genotype = ifelse(grepl('MeCP2 KO',mecp2Dat$pd$Genotype),'MT', 'WT')
ptenDat$pd$Genotype = ifelse(grepl('Pten m3m4/m3m4',ptenDat$pd$Genotype),'MT', 'WT')

#########
# combine
cols = c('FileID', 'Line', 'Genotype','Age')
pd = rbind(maherDat$pd[-31,cols],philpotDat$pd[,cols],sweattDat$pd[,cols],
          mecp2Dat$pd[,cols], mef2cDat$pd[,cols], ptenDat$pd[ptenIndex,cols], ube3aDat$pd[,cols])
pd$Age[pd$Age=='P1']= 'p1'

################################
# load TCF4 mutation gene counts
load('/dcl01/lieber/ajaffe/Brady/mouseRNAseq/rawCounts_tcf4_mouse_OCT20_n36.rda', envir = maherCounts<-new.env())
load('/dcl01/lieber/ajaffe/Brady/philpot/rawCounts_philpot_OCT20_n58.rda', envir = philpotCounts<-new.env())
load('/dcl01/lieber/ajaffe/Brady/sweatt/rawCounts_sweatt_OCT20_n16.rda', envir = sweattCounts<-new.env())

###############################################
# load other syndromic ASD mutation gene counts
load('/dcl01/lieber/ajaffe/Brady/mouseRNAseq/Rett/rawCounts_mecp2_OCT20_n6.rda', envir = mecp2Counts<-new.env())
load('/dcl01/lieber/ajaffe/Brady/mef2c/rawCounts_MEF2C_OCT27_n6.rda', envir = mef2cCounts<-new.env())
load('/dcl01/lieber/ajaffe/Brady/mouseRNAseq/TilotPTEN/rawCounts_pten_OCT20_n18.rda', envir = ptenCounts<-new.env())
load('/dcl01/lieber/ajaffe/Brady/ube3a/rawCounts_ube3a_nov17_n8.rda', envir = ube3aCounts<-new.env())


####################
# add processed data
cols = Reduce(intersect,lapply(list(maherCounts$pd,philpotCounts$pd,sweattCounts$pd,
         mecp2Counts$pd, mef2cCounts$pd,ptenCounts$pd ,ube3aCounts$pd),names))
pd = cbind(pd,rbind(maherCounts$pd[-31,cols],philpotCounts$pd[,cols],sweattCounts$pd[,cols],
           mecp2Counts$pd[,cols], mef2cCounts$pd[,cols], ptenCounts$pd[ptenIndex,cols], ube3aCounts$pd[,cols]))
#all.equal(pd$FileID,pd$SAMPLE_ID) #TRUE

###############################################
# intersect expressed gene/exon/junction counts
mGene = Reduce(intersect,list(rownames(maherCounts$geneMap),rownames(philpotCounts$geneMap),
                              rownames(sweattCounts$geneMap),rownames(mecp2Counts$geneMap),
                              rownames(mef2cCounts$geneMap),rownames(ptenCounts$geneMap),
                              rownames(ube3aCounts$geneMap)))

##############################
# merge annotations and counts
geneMap = maherCounts$geneMap[mGene,]

geneCounts = cbind(maherCounts$geneCounts[mGene,-31],
                   philpotCounts$geneCounts[mGene,],
                   sweattCounts$geneCounts[mGene,],
                   mecp2Counts$geneCounts[mGene,],
                   mef2cCounts$geneCounts[mGene,],
                   ptenCounts$geneCounts[mGene,ptenIndex],
                   ube3aCounts$geneCounts[mGene,])
#all.equal(colnames(geneCounts),pd$SAMPLE_ID) #TRUE
save(pd, geneCounts,geneMap, file = '/dcl01/lieber/ajaffe/Brady/mouseRNAseq/rawCounts_mega_dataset_asd_mice_nov17_n141.rda')
