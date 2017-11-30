# mog average intensity of dentate gyrus in adult mice
library(readxl)
library(lme4)
library(lmerTest)
library(beeswarm)
options(stringsAsFactors = F)

#############################
# read in data and key
dat =data.frame(read_excel('MOG_ICC_data.xlsx'))
key =data.frame(read_excel("MOG_Blind_Key.xlsx"))

#########################
#names(dat)
#[1] "Animal" "Slice"  "Batch"  "Mean"   "Min"    "Max" 

###########
#names(key)
#[1] "ACTUAL ANIMAL ID "     "LABELED AS FOR BLIND " "GENOTYPE" 
key$Batch = rep(1:2,each = 8)

dat$Genotype = factor(key[match(toupper(dat$Animal),toupper(key[,2])),3],
                      levels = c('W','H'),labels = c("WT",'Het'))
dat$Batch = key[match(toupper(dat$Animal),toupper(key[,2])),4]

summary(glmer(data = dat,Median~Genotype+Batch+Area+(1|Animal),family = poisson))
summary(lmer(data = dat,Mean~Genotype+Batch+Area+(1|Animal)))

summary(glm(data = dat,Median~Genotype,family = poisson))
summary(lm(data = dat,Mean~Genotype))

boxplot(Mean/255~Genotype+Batch, data = dat,main = 'Adult mice DG Mog ICC',
        ylab= 'Normalized Intensity')
beeswarm(Mean/255~Genotype+Batch,data = dat,add = T,pch = 21,pwbg = factor(Animal))


boxplot(Median/255~Genotype+Batch, data = dat,main = 'Adult mice DG Mog ICC',
        ylab= 'Normalized Intensity',log = 'y')
beeswarm(Median/255~Genotype+Batch,data = dat,add = T,pch = 21,pwbg = factor(Animal))