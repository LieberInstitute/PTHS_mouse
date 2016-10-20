### Maher mouse Genotype phenotype data
library(jaffelab)
library(parallel)

Dir = "/dcl01/lieber/ajaffe/Brady/mouseRNAseq/FASTQ"
reads = list.files(Dir, pattern="001.fastq.gz$", recur=TRUE, full.names=TRUE)

sampleID = gsub("-","_",ss(reads, "/", 8))
bySample = split(reads,sampleID)
bySample = t(sapply(bySample, function(x) {
  left = paste(x[grep("_R1_",x)],collapse=" ")
  right = paste(x[grep("_R2_",x)],collapse=" ")
  c(left,right)}))

#####################################
# merge files and put in FASTQ folder
write.table(names(split(reads,sampleID)),file = paste0('/dcl01/lieber/ajaffe/Brady/mouseRNAseq/SAMPLE_IDs.txt'),
            col.names = F,row.names = F,quote = F)
R1 = paste0(Dir,'/',names(split(reads,sampleID)),'_R1.fastq.gz')
R2 = paste0(Dir,'/',names(split(reads,sampleID)),'_R2.fastq.gz')
thecall1 = paste('cat',bySample[,1],' > ',R1)
thecall2 = paste('cat',bySample[,2],' > ',R2)
#mclapply(c(thecall1,thecall2),system,mc.cores = 18)

###########################
# make phenotype data
pd = data.frame(SampleID = paste0("Mouse", 1:nrow(bySample)),
                FileID = rownames(bySample), 
                Genotype = ss(rownames(bySample), "_", 5),
                Age = ss(rownames(bySample), "_", 4),
                FlowCell = ss(rownames(bySample), "_", 6),
                Date = ss(rownames(bySample), "_", 2),
                stringsAsFactors=FALSE)

### clean up flowCells
ind = grep("C76VFACXX",pd$FileID)
pd$Genotype[ind] = ss(pd$FileID[ind], "_", 4)
pd$Age[ind] = ss(pd$FileID[ind], "_", 3)
pd$Date[ind] = ss(pd$FileID[ind], "_", 5)
pd$FlowCell[ind] = "C76VFACXX"
pd$Genotype[pd$Genotype == "H"] = "HT"

ind = which(ss(pd$FileID, "_", 5) == "H5TJTBBXX") 
indH = grep("Het",ss(pd$FileID, "_", 4))
pd$Genotype[indH] = "HT"
indW = grep("WT",ss(pd$FileID, "_", 4))
pd$Genotype[indW] = "WT"
indH = grep("Adult",ss(pd$FileID, "_", 4))
pd$Age[indH] = "Adult"
indW = grep("p1",ss(pd$FileID, "_", 4))
pd$Age[indW] = "p1"
pd$FlowCell[ind] = "H5TJTBBXX"
pd$Genotype[pd$Genotype == "H"] = "HT"

# factorize phenotype data
pd$Genotype = factor(pd$Genotype, levels=c("WT", "HT"))
pd$Age = factor(pd$Age,levels = c("p1","p21","Adult"))
save(pd, file = './rdas/pheno.rda')

