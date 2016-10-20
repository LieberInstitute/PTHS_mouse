# make phenotype data for pten mouse 
library(jaffelab)
library(parallel)

# find reads and appropriately rename
Dir = "/dcl01/lieber/ajaffe/Brady/mouseRNAseq/TilotPTEN/FASTQ"
#reads = list.files(Dir, pattern="fastq.gz$", recur=TRUE, full.names=TRUE)
#thecall = paste('cp',reads,gsub('_','_R',reads))
#mclapply(thecall,system,mc.cores = 18)
reads = list.files(Dir, pattern="R[1-2].fastq.gz$", full.names=TRUE)

## make phenotype data
sampleID = ss(ss(reads, "/", 9),"_",1)
bySample = split(reads,sampleID)
bySample = t(sapply(bySample, function(x) {
  left = paste(x[grep("_1.",x)],collapse=",")
  right = paste(x[grep("_2.",x)],collapse=",")
  c(left,right)}))

t = read.delim('./tables/SraRunTable_tilot.txt')
t = t[unlist(sapply(rownames(bySample),grep,t$Run_s)),]

pd = data.frame(SampleID = paste0("Mouse", seq(along=split(reads,sampleID))), 
                FileID = rownames(bySample), stringsAsFactors=FALSE)
pd$Genotype = factor(t$genotype_s,levels = c("Pten wt/wt",'Pten wt/m3m4','Pten m3m4/m3m4'))
pd$Age = factor(t$age_s,levels = c('2-weeks old at sacrifice','6-weeks old at sacrifice'))
save(pd,file = './tables/pheno.rda')

write.table(pd$FileID,file = "/dcl01/lieber/ajaffe/Brady/mouseRNAseq/TilotPTEN/SAMPLE_IDs.txt",
            row.names = F,col.names = F,quote = F)