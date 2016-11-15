# load coverage data (bigWig files) into derfinder and find expressed regions
# qsub -V -l mf=400G,h_vmem=450G,h_stack=256M -m e -M badoi@jhu.edu -cwd -b y R CMD BATCH load_coverage.R
library(derfinder)
library(GenomicRanges)
library(jaffelab)

######################
# load phenotype data and add mapped information
load('./rdas/pheno.rda')
pd = cbind(pd,read.csv('/dcl01/lieber/ajaffe/Brady/zhang/annotated_pd.csv'))

######################
# load BigWig files into full coverage
files = paste0('/dcl01/lieber/ajaffe/Brady/zhang/Coverage/',pd$SampleID,'.bw')
names(files) = pd$SampleID
# all(file.exists(files)) #TRUE
fullCov = fullCoverage(files = files, chrs = paste0("chr", c(1:19,"X","Y")),totalMapped = rep(1, length(files)),
                       verbose = FALSE,targetSize = 10e6,mc.cores = 6) #normalize to 10M reads per sample

########################
# make region matrices
pd$sumMapped = pd$totalMapped + pd$mitoMapped
regionMat = regionMatrix(fullCov, cutoff = 5, L = 101, verbose = FALSE,returnBP=FALSE, 
                         totalMapped = pd$sumMapped,mc.cores=6)
save(regionMat, file="/dcl01/lieber/ajaffe/Brady/zhang/overallRegionMat.rda")

##############################
# extract regions and coverage
regions = unlist(GRangesList(lapply(regionMat, '[[', 'regions')))
names(regions) = NULL
regionMat = do.call("rbind", lapply(regionMat, '[[', 'coverageMatrix'))
#colnames(regionMat) = pd$SampleID

#############################
# region length filter > 50bp
lenIndex = which(width(regions) > 50)
regions = regions[lenIndex,]
regionMat = regionMat[lenIndex,]
save(regionMat,regions, file="/dcl01/lieber/ajaffe/Brady/zhang/filteredRegionMat.rda")



