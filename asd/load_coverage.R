# load coverage data (bigWig files) into derfinder and find expressed regions
# qsub -V -l mf=400G,h_vmem=450G,h_stack=256M -m e -M badoi@jhu.edu -cwd -b y R CMD BATCH load_coverage.R
library(derfinder)
library(GenomicRanges)
library(jaffelab)

######################
# load phenotype data and add mapped information
load('./rdas/pheno.rda')
pd = cbind(pd,read.csv('/dcl01/lieber/ajaffe/Brady/asd/annotated_pd.csv'))

######################
# load BigWig files into full coverage
files = pd$bamFile
names(files) = pd$SampleID
# all(file.exists(files)) #TRUE
#normalize to 40M reads per sample
fullCov = fullCoverage(files = files, chrs = paste0("chr", c(1:19,"X","Y")),
                       totalMapped = rep(1, length(files)),
                       verbose = TRUE,targetSize = 40e6,
                       mc.cores = parallel::detectCores()) 

########################
# make region matrices
pd$sumMapped = pd$totalMapped + pd$mitoMapped
regionMat = regionMatrix(fullCov, cutoff = 5, L = 50, verbose = TRUE,returnBP=FALSE, 
                         totalMapped = pd$sumMapped,mc.cores=parallel::detectCores())
save(regionMat, file="/dcl01/lieber/ajaffe/Brady/asd/overallRegionMat.rda")

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
save(regionMat,regions, file="/dcl01/lieber/ajaffe/Brady/asd/filteredRegionMat.rda")



