# mog average intensity of dentate gyrus in adult mice
library(readxl)
library(lme4)
library(lmerTest)
library(beeswarm)
options(stringsAsFactors = F)

#############################
# read in data and key
dat =read_excel('MOG_ICC_data.xlsx')
key =read_excel("MOG_Blind_Key.xlsx")

#########################
#names(dat)
#[1] "Animal" "Slice"  "Mean"   "Min"    "Max" 

###########
#names(key)
#[1] "ACTUAL ANIMAL ID "     "LABELED AS FOR BLIND " "GENOTYPE" 

dat$Genotype = factor(key[match(toupper(dat$Animal),toupper(key[,2])),3],
                      levels = c('W','H'),labels = c("WT",'Het'))

summary(lmer(data = dat,Mean~Genotype+(1|Animal)))
boxplot(Mean/255~Genotype, data = dat,main = 'Adult mice DG Mog ICC',
        ylab= 'Normalized Intensity',log = 'y')
beeswarm(Mean/255~Genotype,data = dat,add = T,pch = 21,pwbg = factor(Animal))