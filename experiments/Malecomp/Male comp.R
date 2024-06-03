###################################################################################
###################################################################################
###################################################################################
##### All blocks, 1, 2 and 3 ##########################


### INDEX
# data organisation and modification: 13 - 72
# GLMMs testing for the effect of fru allele and genotype on mating success; 73 - 152
# parametric bootsrapping: 153 - 180
# plots: 181 - 226 (end)

### set working directory
setwd("/Users/michaeljardine/Desktop/DTP/Data and analysis/fruitless/Malecomp")

## PACKAGES
library(Rmisc)
library(lme4)
library(aod)
library(lsmeans)
library(ggplot2)
library(pbkrtest)
library(car)


Tcomp <- read.csv("ALL_3_BLOCKS_ALT.csv")
summary(Tcomp)
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

### change alleles fom M and F to S (short) and L (long)
levels(Tcomp$Allele) <- c("S", "L")

## change pupa to a 1/0 score to see who the winner is
Tcomp$winner <- ifelse(Tcomp$Pupa %in% c("T"), 0, 1)
## winner: competitor = 0 (Tub), focal = 1 (Norm)
str(Tcomp)
summary(Tcomp)

## combine the line and the genotype to create 5B/5D, 8B/8D etc - this is for the plots later on 
Tcomp$Linegenotype <- interaction(Tcomp$Line, Tcomp$Genotype)
str(Tcomp)

## combine for genotype with introgressed part
Tcomp$Allelegenotype <- interaction(Tcomp$Allele, Tcomp$Genotype)
Tcomp$Lineblock <- interaction(Tcomp$Line, Tcomp$Block)

str(Tcomp)
##is this final combination necessary? - probs not dosen't explain anything more than the previous ones and won't be as useful in the model


## some summary statistics for means and errors 
summarySE(Tcomp, measurevar=c("winner"), groupvars=c("Allele"))
summarySE(Tcomp, measurevar=c("winner"), groupvars=c("Line"))
summarySE(Tcomp, measurevar=c("winner"), groupvars=c("Genotype"))
summarySE(Tcomp, measurevar=c("winner"), groupvars=c("Linegenotype"))
summarySE(Tcomp, measurevar=c("winner"), groupvars=c("Allelegenotype"))
summarySE(Tcomp, measurevar=c("winner"), groupvars=c("Lineblock"))

Tcomp$linegenotype <- interaction(Tcomp$Line, Tcomp$Genotype)
Tcomp$BLG <- interaction(Tcomp$linegenotype, Tcomp$Block)
summarySE(Tcomp, measurevar=c("winner"), groupvars=c("BLG"))



#### MODELS ####

### What we wnat to know if how the verison of the fruitless allele affects the competitive performace of the males
### The way we have meausred competitive performance is the phenotype of the pupa, which directly reveals the fru allele possed by the father
### We have six lines, 3 with short fru alleles and 3 with long fru alleles - these are basically replicates
### we also have two different genotypes, B and D, for each line. This provides information on the dominace of the allele
### the assay has been carried out over sereval days due to logistics. There may be variation in the data set depending on day (Round) - random factor
### the same goes for the three blocks


### out 1st model will simply look at the if there are differences between those flies with the male or female allele, taking into account the variation across lines and day
### this has block, line and round as random effects - not taking into account any experimental desgin structure
### since the variable of interest is binary we will use a binomial error distribution and logit link function


TGLMM1.1 <- glmer(winner ~ Allele + (1 | Line) + (1 | Block) + (1 | Round), data=Tcomp, family=binomial(link=logit))

res <- resid(TGLMM1.1)
ran <- ranef(TGLMM1.1)
fit <- fitted(TGLMM1.1)
par(mfrow=c(1,5))
hist(ran$Line[,1], main='')
hist(ran$Round[,1], main='')
hist(ran$Block[,1], main='')
plot(res ~ fit)
qqnorm(res, main='')
qqline(res)

summary(TGLMM1.1)
anova(TGLMM1.1)
Anova(TGLMM1.1)
### with just the flies fru allele included then there is no effect on mating success


### the flies differ in their genotype (B or D) so we will inclide this in the model also to account for this variation 
### keep the simple random effect for now
GLMM2.1 <- glmer(winner ~ Allele*Genotype +  (1 | Line) + (1 | Block) + (1 | Round), data=Tcomp, family=binomial(link=logit))

res <- resid(GLMM2.1)
ran <- ranef(GLMM2.1)
fit <- fitted(GLMM2.1)
par(mfrow=c(1,5))
hist(ran$Line[,1], main='')
hist(ran$Round[,1], main='')
hist(ran$Block[,1], main='')
plot(res ~ fit)
qqnorm(res, main='')
qqline(res)
summary(GLMM2.1)
anova(GLMM2.1)
Anova(GLMM2.1)
### no real difference, but cements the effect of the genotype


### it was mentioned before that we may want to account fo the fact that the design is imbalanced, e.g. each line does not have both fru alleles

GLMM2.3 <- glmer(winner ~ Allele*Genotype +  (1 + Allele | Line) + (1 | Block) + (1 | Round), data=Tcomp, family=binomial(link=logit))
### warning messages but still works
res <- resid(GLMM2.3)
ran <- ranef(GLMM2.3)
fit <- fitted(GLMM2.3)
par(mfrow=c(1,5))
hist(ran$Line[,1], main='')
hist(ran$Round[,1], main='')
hist(ran$Block[,1], main='')
plot(res ~ fit)
qqnorm(res, main='')
qqline(res)
summary(GLMM2.3)
anova(GLMM2.3)
Anova(GLMM2.3)
## very simialr results but with error messages on the model
anova(GLMM2.1, GLMM2.3)
## doesn't appear to be justified leave with more basic random effects structure

#### parametric bootstrapping ####
# based on model GLMM2.1

comp1 <- glmer(winner ~ 1 +  (1 | Line) + (1 | Block) + (1 | Round), data=Tcomp, family=binomial(link=logit))
comp2 <- glmer(winner ~ Allele +  (1 | Line) + (1 | Block) + (1 | Round), data=Tcomp, family=binomial(link=logit))
comp3 <- glmer(winner ~ Genotype +  (1 | Line) + (1 | Block) + (1 | Round), data=Tcomp, family=binomial(link=logit))
comp4 <- glmer(winner ~ Allele + Genotype +  (1 | Line) + (1 | Block) + (1 | Round), data=Tcomp, family=binomial(link=logit))
comp5 <- glmer(winner ~ Allele*Genotype +  (1 | Line) + (1 | Block) + (1 | Round), data=Tcomp, family=binomial(link=logit))

MC.allele <- PBmodcomp(comp4, comp3, nsim=1000)
MC.allele
MC.chro <- PBmodcomp(comp4, comp2, nsim=1000)
MC.chro
MC.inter <- PBmodcomp(comp5, comp4, nsim=1000)
MC.inter

MC.inter2 <- PBmodcomp(comp5, comp4, nsim=1000)
MC.inter2


#### plots ####

### same as with eggs we would like ones with allele, line and allele by genotype


d1=data.frame(Allele=c("Long","Short"), mean=c(0.436, 0.3881), lower=c(-0.054, -0.0921), upper=c(0.926, 0.868))
d1$allele <- c("L", "S")


ggplot() + 
  geom_pointrange(data=d1, mapping=aes(x=Allele, y=mean, ymax=upper, ymin=lower, col=allele), size=1.5) + 
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold")) +
  scale_y_continuous(name="proportion of matings by focal males") 



### winner by line
d2=data.frame(Line=c("05", "08", "16", "19", "24", "95"), mean=c(0.32, 0.5, 0.34, 0.33, 0.58, 0.39), lower=c(-.015, 0, -0.13, -0.14, 0.09, -0.1), upper=c(0.79, 1, 0.81, 0.8, 1.07, 0.88)) 
d2$allele <- c("S", "S", "S", "L", "L", "L")

ggplot() + 
  geom_pointrange(data=d2, mapping=aes(x=Line, y=mean, ymax=upper, ymin=lower, col=allele), size=1.5) + 
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold")) +
  scale_y_continuous(name="proportion of matings by focal males") 


### winner by allele and genotype
d3=data.frame(Allelegenotype=c("Short, B", "Short, D", "Long, B", "Long, D"), mean=c(0.47, 0.3, 0.47, 0.4), lower=c(-0.03, -0.16, -0.03, -0.1), upper=c(0.97, 0.76, 0.97, 0.9))
d3$allele <- c("S", "S", "L", "L")

ggplot() + 
  geom_pointrange(data=d3, mapping=aes(x=Allelegenotype, y=mean, ymax=upper, ymin=lower, col=allele), size=1.5) + 
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold")) +
  scale_y_continuous(name="proportion of matings by focal males") 


### winner by genotype
d4=data.frame(Genotype=c("Balancer", "Deletion"), mean=c(0.47, 0.35), lower=c(-0.03, -0.13), upper=c(0.97, 0.83))
d4$geno <- c("B", "D")

ggplot() + 
  geom_pointrange(data=d4, mapping=aes(x=Genotype, y=mean, ymax=upper, ymin=lower, col=geno), size=1.5) + 
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold")) +
  scale_y_continuous(name="proportion of matings by focal males") 

### Filip's code






### making plots from new data frames cause I'm not very good

Abars <- read.csv("allelebars.csv")
str(Abars)

bullshit <- ggplot(aes(y=winsfrac, x=Line, fill=Allele), data=Abars) +
  geom_bar(stat="identity", colour="black", position="dodge") +
  geom_errorbar(aes(ymin=winsfrac-se, ymax=winsfrac+se),
                  width=.2,
                  position=position_dodge(.9)) +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
          panel.background=element_blank(), axis.line=element_line(colour="black")) +
  theme(axis.title.y = element_text(size=16),
          axis.title.x = element_text(size=18),
          axis.text.x=element_text(size=15, colour="black"),
          axis.text.y=element_text(size=15)) +
  scale_y_continuous(name="Proportion of matings by focal males") +
  scale_x_discrete(labels=c("M05", "M08", "M16", "M19", "M24", "M95"))
bullshit

# WHY HAVE THE AXIS GONE?? THIS ALL WORKED  

Ibars <- read.csv("INTER BARS.csv")
str(Ibars)
Ibars$linegenotype <- interaction(Ibars$line, Ibars$geno)
str(Ibars)


bullshit2 <- ggplot(aes(y=winfrac, x=line, fill=allele), data=Ibars) +
  geom_bar(stat="identity", colour="black", position="dodge") +
  geom_errorbar(aes(ymin=winfrac-se, ymax=winfrac+se),
                width=.2,
                position=position_dodge(.9)) +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        panel.background=element_blank(), axis.line=element_line(colour="black")) +
  theme(axis.title.y = element_text(size=16),
        axis.title.x = element_text(size=18),
        axis.text.x=element_text(size=15, colour="black"),
        axis.text.y=element_text(size=15)) +
  scale_y_continuous(name="Proportion of matings by focal males") +
  scale_x_discrete(labels=c("M05", "M08", "M16", "M19", "M24", "M95")) +
  facet_wrap(~geno,scales="free_x")
bullshit2


