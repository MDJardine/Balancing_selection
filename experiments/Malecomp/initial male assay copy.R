setwd("/Users/michaeljardine/Desktop/DTP/Datasets/Malecomp")

### analsis for block 1 and been altered and partailly destroed while fiing to make it wok properly for block 2
### it will actaully be ver simailr to block 2 apart from the notes 
### paper print out from this block is sill availble and can be used as reference



##### Analysis and prep for Block 2 (August 2018) #####


comp <-read.csv("malecompassayblock2FULL.csv")

## PACKAGES
library(Rmisc)
library(lme4)
library(aod)
library(lsmeans)
library(ggplot2)

## read in data and check

## change Line to a factor
comp$Line <- factor(comp$Line)
## change Round to a factor
comp$Round <- factor(comp$Round)
## change Block to a factor
comp$Block <- factor(comp$Block)

str(comp)

## remove NA

comp <- na.omit(comp)

str(comp)
summary(comp)



## change pupa to a 1/0 score to see who the winner is
comp$winner <- ifelse(comp$Pupa %in% c("T"), 0, 1)

## winner: competitor = 0 (Tub), focal = 1 (Norm)



## some simple plots
par(mfrow=c(2,2))
plot(winner ~ Allele, data=comp)
plot(winner ~ Round, data=comp)
plot(winner ~ Genotype, data=comp)
plot(winner ~ Line, data=comp)
## these don't seem greatly useful, probs they use the median and we only have 1/0 data - lsmeans plots may be better later

## some summary statistics for means and errors 
summarySE(comp, measurevar=c("winner"), groupvars=c("Allele"))
summarySE(comp, measurevar=c("winner"), groupvars=c("Line"))
summarySE(comp, measurevar=c("winner"), groupvars=c("Genotype"))
summarySE(comp, measurevar=c("winner"), groupvars=c("Linegenotype"))
summarySE(comp, measurevar=c("winner"), groupvars=c("Allelegenotype"))

## combine the line and the genotype to create 5B and 5D 
comp$Linegenotype <- interaction(comp$Line, comp$Genotype)
str(comp)
summarySE(comp, measurevar=c("winner"), groupvars=c("Linegenotype"))
## comnine for genotype with introgressed part
comp$Allelegenotype <- interaction(comp$Allele, comp$Genotype)
str(comp)
## combine these so that we have a knowledge of the introgessed background, the line and the genotype of each vial
comp$Allelelinegenotype <- interaction(comp$Allele, comp$Linegenotype)
##is this final combination necessary? - probs not dosen't explain anything more than the previous ones and won't be as useful in the model


lines6=data.frame(line=c("05", "08", "16", "19", "24", "95"),
                  mean=c(0.291, 0.538, 0.365, 0.387, 0.705, 0.431),
                  lower=c(0.248, 0.486, 0.312, 0.336, 0.653, 0.369), 
                  upper=c(0.333, 0.59, 0.417, 0.438, 0.757, 0.493))
lines6
str(lines6)


ggplot() + 
  geom_pointrange(data=lines6, mapping=aes(x=line, y=mean, ymin=lower, ymax=upper)) + #width=0.2, size=1, color="blue", fill="white", shape=22) + 
  theme(axis.title.y = element_text(size=20), 
        axis.title.x = element_text(size=20), 
        axis.text.x=element_text(size=16, colour="black"), 
        axis.text.y=element_text(size=16, colour="black")) +
  scale_y_continuous(name="Proportion of matings by focal males") +
  coord_cartesian(ylim=c(0.0, 1.0))




#### plot out the means of the12 lines to eyeball for patterns
#lines12=data.frame(line=c("05B","05D", "08B", "08D", "16B", "16D", "19B", "19D", "24B", "24D", "95B", "95D"),
 #                  mean=c(0.521, 0.412, 0.625, 0.611, 0.409, 0.353, 0.42, 0.667, 0.7, 0.714, 1, 0.333),
  #                 lower=c(0.415, 0.289, 0.524, 0.493, 0.302, 0.234, 0.309, 0.553, 0.547, 0.613, 1, 0.228),
   #                upper=c(0.627, 0.535, 0.726, 0.729, 0.516, 0.427, 0.531, 0.781, 0.853, 0.815, 1, 0.438))
lines12
## plot mean +/- SE with ggplot
ggplot() + 
  geom_pointrange(data=lines12, mapping=aes(x=line, y=mean, ymin=lower, ymax=upper)) + #width=0.2, size=1, color="blue", fill="white", shape=22) + 
  theme(axis.title.y = element_text(size=20), 
        axis.title.x = element_text(size=20), 
        axis.text.x=element_text(size=16, colour="black"), 
        axis.text.y=element_text(size=16, colour="black")) +
  scale_y_continuous(name="Proportion of matings by focal males") +
  coord_cartesian(ylim=c(0.0, 1.0))

## background part: theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), 
## panel.background=element_blank(), axis.line=element_line(colour="black")) + 


### What we wnat to know if how the verison of the fruitless allele affestc the copeitive performace of the males
### The way we have meausred compeitve performance is the phenotype of the pupa, whcih directly reveals the phenotype of the father
### We have six lines, 3 male alleles and 3 with females alleles - these are basically replicates
### we also have two diffeent genotypes, B and D, for each line. This provides information on the dominace of the allele
### the assay has been carried out over sereval days use to logistics. There may be variation in the data set depending on day (Round) - random factor


### out 1st model will simply look at the i there are differences ebtween those flies with the male or female allele, taking into account the variation acros lines and day

GLMM1 <- glmer(winner ~ Allele + (1 | Round) + (1 | Line), data=comp, family=binomial(link=logit))

res <- resid(GLMM1)
ran <- ranef(GLMM1)
fit <- fitted(GLMM1)
par(mfrow=c(1,4))
hist(ran$Round[,1], main='')
plot(res ~ fit)
qqnorm(res, main='')
qqline(res)

summary(GLMM1)
anova(GLMM1)
wald.test(b=fixef(GLMM1), Sigma=vcov(GLMM1), Terms = 2, df=1)

### appears to be no difference
### even if there was one, from the summar stats above it is likely to be in the opposite direction to the one expected
### we can try other arabgements of he random factors to get better diagnostics ut this is unlikely to alter the main outcome

### the flies differ in their genotype (B or D) so we will inclide this in the model also to account for this variation 

GLMM2 <- glmer(winner ~ Allele*Genotype + (1 | Round) +(1 | Line), data=comp, family=binomial(link=logit))

res <- resid(GLMM2)
ran <- ranef(GLMM2)
fit <- fitted(GLMM2)
par(mfrow=c(1,3))
hist(ran$Round[,1], main='')
plot(res ~ fit)
qqnorm(res, main='')
qqline(res)

summary(GLMM2)
anova(GLMM2)

wald.test(b=fixef(GLMM2), Sigma=vcov(GLMM2), Terms = 2, df=1)
wald.test(b=fixef(GLMM2), Sigma=vcov(GLMM2), Terms = 3, df=1)
wald.test(b=fixef(GLMM2), Sigma=vcov(GLMM2), Terms = 4, df=1)

###somewthat mixed results 
### from first look appears still no affect of the fruitless allele
### High F-value for genotype indictaing a srtong affect, males with D genoype are less successful
### also moderatly high F for interction, males with M allele and D genotype are particularly bad.

### becomes a bit confusing when looking at the different metrics.
### althought genotype has a higher F and the interqaction, it does not come out as significant in the summary and wald tests, althog it would on a traditional p-value caluclator
### this is rpobably something to dicusss with Max and kevin next week

### agai we can change up the random factos to find the best arrangement but maybe leave for now.




##### Combine blocks 1 and 2 #####

Tcomp <- read.csv("blocks1and2.csv")

str(Tcomp)

## change Line to a factor
Tcomp$Line <- factor(Tcomp$Line)
## change Round to a factor
Tcomp$Round <- factor(Tcomp$Round)
## change Block to a factor
Tcomp$Block <- factor(Tcomp$Block)

str(Tcomp)

## remove NA
Tcomp <- na.omit(Tcomp)

str(Tcomp)
summary(Tcomp)

## change pupa to a 1/0 score to see who the winner is
Tcomp$winner <- ifelse(Tcomp$Pupa %in% c("T"), 0, 1)

## winner: competitor = 0 (Tub), focal = 1 (Norm)

str(Tcomp)
summary(Tcomp)

## some simple plots
par(mfrow=c(2,2))
plot(winner ~ Allele, data=comp)
plot(winner ~ Round, data=comp)
plot(winner ~ Genotype, data=comp)
plot(winner ~ Line, data=comp)
## these don't seem greatly useful, probs they use the median and we only have 1/0 data - lsmeans plots may be better later

## combine the line and the genotype to create 5B and 5D 
Tcomp$Linegenotype <- interaction(Tcomp$Line, Tcomp$Genotype)
str(Tcomp)
## comnine for genotype with introgressed part
Tcomp$Allelegenotype <- interaction(Tcomp$Allele, Tcomp$Genotype)
str(Tcomp)
## combine these so that we have a knowledge of the introgessed background, the line and the genotype of each vial
Tcomp$Allelelinegenotype <- interaction(Tcomp$Allele, Tcomp$Linegenotype)
str(Tcomp)
##is this final combination necessary? - probs not dosen't explain anything more than the previous ones and won't be as useful in the model

## some summary statistics for means and errors 
summarySE(Tcomp, measurevar=c("winner"), groupvars=c("Allele"))
summarySE(Tcomp, measurevar=c("winner"), groupvars=c("Line"))
summarySE(Tcomp, measurevar=c("winner"), groupvars=c("Genotype"))
summarySE(Tcomp, measurevar=c("winner"), groupvars=c("Linegenotype"))
summarySE(Tcomp, measurevar=c("winner"), groupvars=c("Allelegenotype"))

## these summary statistics show simialr paterns to those from block 2 - maybe because B2 out numers B1 b more than 2:1


### What we wnat to know if how the verison of the fruitless allele affects the competitive performace of the males
### The way we have meausred competitive performance is the phenotype of the pupa, which directly reveals the phenotype of the father
### We have six lines, 3 male alleles and 3 with females alleles - these are basically replicates
### we also have two different genotypes, B and D, for each line. This provides information on the dominace of the allele
### the assay has been carried out over sereval days use to logistics. There may be variation in the data set depending on day (Round) - random factor
### the same goes for the two blocks

### out 1st model will simply look at the i there are differences ebtween those flies with the male or female allele, taking into account the variation acros lines and day

GLMM3 <- glmer(winner ~ Allele + (1 | Round) + (1 | Line) + (1 | Block), data=Tcomp, family=binomial(link=logit))

res <- resid(GLMM3)
ran <- ranef(GLMM3)
fit <- fitted(GLMM3)
par(mfrow=c(1,4))
hist(ran$Round[,1], main='')
plot(res ~ fit)
qqnorm(res, main='')
qqline(res)

summary(GLMM3)
anova(GLMM3)
wald.test(b=fixef(GLMM3), Sigma=vcov(GLMM1), Terms = 2, df=1)

### appears to be no difference
### even if there was one, from the summary stats above it is likely to be in the opposite direction to the one expected
### we can try other arrangements of the random factors to get better diagnostics but this is unlikely to alter the main outcome

### the flies differ in their genotype (B or D) so we will inclide this in the model also to account for this variation 

GLMM4 <- glmer(winner ~ Allele*Genotype + (1 | Round) + (1 | Line) + (1 | Block), data=Tcomp, family=binomial(link=logit))

res <- resid(GLMM4)
ran <- ranef(GLMM4)
fit <- fitted(GLMM4)
par(mfrow=c(1,3))
hist(ran$Round[,1], main='')
plot(res ~ fit)
qqnorm(res, main='')
qqline(res)

summary(GLMM4)
anova(GLMM4)

wald.test(b=fixef(GLMM2), Sigma=vcov(GLMM4), Terms = 2, df=1)
wald.test(b=fixef(GLMM2), Sigma=vcov(GLMM4), Terms = 3, df=1)
wald.test(b=fixef(GLMM2), Sigma=vcov(GLMM4), Terms = 4, df=1)

###somewthat mixed results 
### from first look appears still no affect of the fruitless allele
### High F-value for genotype indictaing a srtong affect, males with D genoype are less successful
### also moderatly high F for interction, males with M allele and D genotype are particularly bad.


### suggested to alter the random effects structure of the model from max

GLMM5 <- glmer(winner ~ Allele*Genotype + (1 | Round) + (1 + Genotype | Line) + (1 | Block), data=Tcomp, family=binomial(link=logit))

res <- resid(GLMM5)
ran <- ranef(GLMM5)
fit <- fitted(GLMM5)
par(mfrow=c(1,3))
hist(ran$Round[,1], main='')
plot(res ~ fit)
qqnorm(res, main='')
qqline(res)

summary(GLMM5)
anova(GLMM5)

wald.test(b=fixef(GLMM5), Sigma=vcov(GLMM5), Terms = 2, df=1)
wald.test(b=fixef(GLMM5), Sigma=vcov(GLMM5), Terms = 3, df=1)
wald.test(b=fixef(GLMM5), Sigma=vcov(GLMM5), Terms = 4, df=1)

### becomes a bit confusing when looking at the different metrics.
### althought genotype has a higher F and the interqaction, it does not come out as significant in the summary and wald tests, althog it would on a traditional p-value caluclator

### max suggested using the pbkrtest package to help with this

library(pbkrtest)







######## remove the poorly performing 95 line. does this help? #####

fivelines <- read.csv("malecompassayno95.csv")


## change Line to a factor
fivelines$Line <- factor(fivelines$Line)
## change Round to a factor
fivelines$Round <- factor(fivelines$Round)
## change Block to a factor
fivelines$Block <- factor(fivelines$Block)

## remove NA
fivelines <- na.omit(fivelines)

## change pupa to a 1/0 score to see who the winner is
fivelines$winner <- ifelse(fivelines$Pupa %in% c("Tub"), 0, 1)

## winner: competitor = 0 (Tub), focal = 1 (Norm)
str(fivelines)

lines5 <- glmer(winner ~ Line*Genotype + (1 | Round), data=fivelines, family=binomial(link=logit))

res <- resid(lines5)
ran <- ranef(lines5)
fit <- fitted(lines5)
par(mfrow=c(1,3))
hist(ran$Round[,1], main='')
plot(res ~ fit)
qqnorm(res, main='')
qqline(res)

summary(lines5)
anova(lines5)
wald.test(b=fixef(lines5), Sigma=vcov(lines5), Terms = 2, df=4)
wald.test(b=fixef(lines5), Sigma=vcov(lines5), Terms = 3, df=1)
wald.test(b=fixef(lines5), Sigma=vcov(lines5), Terms = 4, df=4)

