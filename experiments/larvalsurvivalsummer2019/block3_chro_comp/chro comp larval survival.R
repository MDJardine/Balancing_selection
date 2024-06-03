setwd("/Users/michaeljardine/Desktop/DTP/Data and analysis/fruitless/larvalsurvivalsummer2019/block3_chro_comp")
Larval.surv.chro <- read.csv("Larval survival chro compliment complete.csv")
str(Larval.surv.chro)
summary(Larval.surv.chro)


## packages
library('lme4')
library('ggplot2')
library('Rmisc')
library('car')
library('chemCal')
library("RColorBrewer")
library("ggsci")
library("plyr")


#### Look at data and genrate summary stats ####


## histograms of males, females and totals
par(mfrow=c(1,3))
hist(Larval.surv.chro$total)
hist(Larval.surv.chro$males)
hist(Larval.surv.chro$females)
## not completly normal, sex specific looks zero skewed - try logging?
hist(log(Larval.surv.chro$males))
hist(log(Larval.surv.chro$females))
## not much better

shapiro.test(Larval.surv.chro$total)
shapiro.test(log(Larval.surv.chro$total))


## also look at the proportions of the two sexes
par(mfrow=c(1,3))
hist(sqrt(Larval.surv.chro$propsurv))
hist(Larval.surv.chro$propmale)
hist(Larval.surv.chro$propfem)
## all look fairly normal


## calculating means to help interperate model results later
# firest create an interaction term between allele and genotype
Larval.surv.chro$allelechro <- interaction(Larval.surv.chro$allele, Larval.surv.chro$chro)

# means for total surviving offspring
summarySE(Larval.surv.chro, measurevar=c("total"), groupvars=c("allele"))
# L = 13.42, S = 8.875

summarySE(Larval.surv.chro, measurevar=c("total"), groupvars=c("chro"))
# B = 10.24444, d = 12.55556

summarySE(Larval.surv.chro, measurevar=c("total"), groupvars=c("allelechro"))
# L.B = 12.22, S.B = 7.775, L.D = 14.62, S.D = 9.75

summarySE(Larval.surv.chro, measurevar=c("total"), groupvars=c("Line"))
# L1 = 11.38, L2  = 12.815, L3 = 15.5278, S1 = 8.454, S2 = 10, S3 = 8.625

# means for proportion surviving
summarySE(Larval.surv.chro, measurevar=c("propsurv"), groupvars=c("allele"))
# L = 0.5368, S = 0.3550

summarySE(Larval.surv.chro, measurevar=c("propsurv"), groupvars=c("chro"))
# B = 0.40978, D = 0.5022

summarySE(Larval.surv.chro, measurevar=c("propsurv"), groupvars=c("allelechro"))
# L.B = 0.4888, S.B = 0.311, L.D = 0.5848, S.D = 0.399

summarySE(Larval.surv.chro, measurevar=c("propsurv"), groupvars=c("Line"))
# L1 = 0.4554, L2  = 0.513, L3 = 0.621, S1 = 0.338, S2 = 0.4, S3 = 0.345

# means for proportion male - sexratio
summarySE(Larval.surv.chro, measurevar=c("propmale"), groupvars=c("allele"))
# L = 0.476, S = 0.4948

summarySE(Larval.surv.chro, measurevar=c("propmale"), groupvars=c("chro"))
# B = 0.4717, D = 0.4973

summarySE(Larval.surv.chro, measurevar=c("propmale"), groupvars=c("allelechro"))
# L.B = 0.4756, S.B = 0.4663, L.D = 0.4765, S.D = 0.5233

summarySE(Larval.surv.chro, measurevar=c("propmale"), groupvars=c("Line"))
# L1 = 0.455, l2  = 0.4902, L3 = 0.4769, S1 = 0.4963, S2 = 0.5022, S3 = 0.4813

## Before startign the models create another interaction term
# 2 measures of each vial are present in the data set since offspring split at the pupa stage
Larval.surv.chro$vial_line <- interaction(Larval.surv.chro$Line, Larval.surv.chro$vial)
# include this as a random variable in the models


#################################################################################################
#### testing for differnences in number of surviving offsrping ####

## does the number of flies that reach adulthood vary depending in the fru allele it has????

## construct models to look at the total number of surviving offspring
## 2 models, 1 with normal error and theother with poisson - standard for count data
## look at full model with interaction term

M.total.chro.1 <- lmer(total ~ allele*chro + (1 | Line) + (1 | vial_line), data=Larval.surv.chro)

res <- resid(M.total.chro.1)
ran <- ranef(M.total.chro.1)
fit <- fitted(M.total.chro.1)
par(mfrow=c(1,4))
hist(ran$Line[,1], main='')
hist(ran$vial_line[,1], main='')
plot(res ~ fit)
qqnorm(res, main='')
qqline(res)
summary(M.total.chro.1)
anova(M.total.chro.1)
Anova(M.total.chro.1)
# an affect of allele and chro

## try poisson for count data 
M.total.chro.2 <- glmer(total ~ allele*chro + (1 | Line) + (1 | vial_line), data=Larval.surv.chro, family = poisson(link = log))

res <- resid(M.total.chro.2)
ran <- ranef(M.total.chro.2)
fit <- fitted(M.total.chro.2)
par(mfrow=c(1,4))
hist(ran$Line[,1], main='')
hist(ran$vial_line[,1], main='')
plot(res ~ fit)
qqnorm(res, main='')
qqline(res)
summary(M.total.chro.2)
anova(M.total.chro.2)
Anova(M.total.chro.2)

anova(M.total.chro.1, M.total.chro.2)
##not actually any better - AIC higher, diagnostic plots messier

## try logging data innstead?

M.total.chro.3 <- lmer(log(total) ~ allele*chro + (1|Line) + (1 | vial_line), data=Larval.surv.chro)

res <- resid(M.total.chro.3)
ran <- ranef(M.total.chro.3)
fit <- fitted(M.total.chro.3)
par(mfrow=c(1,4))
hist(ran$Line[,1], main='')
hist(ran$vial_line[,1], main='')
plot(res ~ fit)
qqnorm(res, main='')
qqline(res)
summary(M.total.chro.3)
anova(M.total.chro.3)
Anova(M.total.chro.3)
## plots look worse

anova(M.total.chro.1, M.total.chro.2, M.total.chro.3)
## but AIC is much less and loglikelihood is much higher - why?
## logging affecting values but not their fit?

## maybe better to keep as logged?

## futher model accounts for nested nature of the design

M.total.chro.4 <- lmer(log(total) ~ allele*chro + (1 | allele:Line) + (1 | vial_line), data=Larval.surv.chro)

res <- resid(M.total.chro.4)
ran <- ranef(M.total.chro.4)
fit <- fitted(M.total.chro.4)
par(mfrow=c(1,4))
hist(ran$"allele:Line"[,1], main='')
hist(ran$vial_line[,1], main='')
plot(res ~ fit)
qqnorm(res, main='')
qqline(res)
summary(M.total.chro.4)
anova(M.total.chro.4)
Anova(M.total.chro.4)

anova(M.total.chro.3, M.total.chro.4)
## pretty much eaxtly the same - probably btter due to the structur of the experiment.

## now drop the interaction term

M.total.chro.5 <- lmer(log(total) ~ allele + chro + (1 | allele:Line) + (1 | vial_line), data=Larval.surv.chro)

res <- resid(M.total.chro.5)
ran <- ranef(M.total.chro.5)
fit <- fitted(M.total.chro.5)
par(mfrow=c(1,4))
hist(ran$'allele:Line'[,1], main='')
hist(ran$vial_line[,1], main='')
plot(res ~ fit)
qqnorm(res, main='')
qqline(res)
summary(M.total.chro.5)
anova(M.total.chro.5)
Anova(M.total.chro.5)

## S allele flies have fewer surviving offspring
## D chromosome flies have more surviving offspring
## max prefers doing this parametric bootstrapping  so lets just do it
## based on M.total.chro.4 and 5 

tot1 <- lmer(log(total) ~ 1 + (1 | allele:Line) + (1 | vial_line), data=Larval.surv.chro)
tot2 <- lmer(log(total) ~ allele + (1 | allele:Line) + (1 | vial_line), data=Larval.surv.chro)
tot3 <- lmer(log(total) ~ chro + (1 | allele:Line) + (1 | vial_line), data=Larval.surv.chro)
tot4 <- lmer(log(total) ~ allele + chro + (1 | allele:Line) + (1 | vial_line), data=Larval.surv.chro)
tot5 <- lmer(log(total) ~ allele*chro + (1 | allele:Line) + (1 | vial_line), data=Larval.surv.chro)

total.off.allele <- PBmodcomp(tot4, tot3, nsim=1000)
total.off.allele
total.off.chro <- PBmodcomp(tot4, tot2, nsim=1000)
total.off.chro
total.off.inter <- PBmodcomp(tot5, tot4)
total.off.inter


########################################################################################################
#### Testing for differenecs in the proportion of eggs surviving to adult hood #### 

## tring with the proportion surving from both sexes - basically jsut the same measurement as before as same starting numbers
## but not sure if this response will be more intuitive and comparable when writing up later
## therefore just doing final 2 models version from the tests above which I'll still call models 4 and 5 incase we add 1-3 later

M.propsurv.chro.4 <- lmer(propsurv ~ allele*chro + (1 | allele:Line) + + (1 | vial_line), data=Larval.surv.chro)
res <- resid(M.propsurv.chro.4)
ran <- ranef(M.propsurv.chro.4)
fit <- fitted(M.propsurv.chro.4)
par(mfrow=c(1,4))
hist(ran$`allele:Line`[,1], main='')
hist(ran$vial_line[,1], main='')
plot(res ~ fit)
qqnorm(res, main='')
qqline(res)
summary(M.propsurv.chro.4)
anova(M.propsurv.chro.4)
Anova(M.propsurv.chro.4)
## looks ok

## non important interacytion term - drop
M.propsurv.chro.5 <- lmer(propsurv ~ allele + chro + (1 | allele:Line) + (1 | vial_line), data=Larval.surv.chro)
res <- resid(M.propsurv.chro.5)
ran <- ranef(M.propsurv.chro.5)
fit <- fitted(M.propsurv.chro.5)
par(mfrow=c(1,4))
hist(ran$`allele:Line`[,1], main='')
hist(ran$vial_line[,1], main='')
plot(res ~ fit)
qqnorm(res, main='')
qqline(res)
summary(M.propsurv.chro.5)
anova(M.propsurv.chro.5)
Anova(M.propsurv.chro.5)

anova(M.propsurv.chro.4, M.propsurv.chro.5) # models pretty similar

## same as before - affect of both allele and chromosome
## smaller proportion of S allele flies survive to adult hood
## more D flies survive to adult hood

pro1 <- lmer(propsurv ~ 1 + (1 | allele:Line) + (1 | vial_line), data=Larval.surv.chro)
pro2 <- lmer(propsurv ~ allele + (1 | allele:Line) + (1 | vial_line), data=Larval.surv.chro)
pro3 <- lmer(propsurv ~ chro + (1 | allele:Line) + (1 | vial_line), data=Larval.surv.chro)
pro4 <- lmer(propsurv ~ allele + chro + (1 | allele:Line) + (1 | vial_line), data=Larval.surv.chro)
pro5 <- lmer(propsurv ~ allele*chro + (1 | allele:Line) + + (1 | vial_line), data=Larval.surv.chro)
             
pro.surv.allele <- PBmodcomp(pro4, pro3, nsim=1000)
pro.surv.allele
pro.surv.chro <- PBmodcomp(pro4, pro2, nsim=1000)
pro.surv.chro
pro.surv.inter <- PBmodcomp(pro5, pro4)
pro.surv.inter

### extra tests splitting by balancer genotype 
## analyse within each of the these balancer genotypes
balancersplitsurv <- split(Larval.surv.chro, Larval.surv.chro$chro)
survB <- balancersplitsurv[[1]]
survD <- balancersplitsurv[[2]]
str(survB)
str(survD)

## B only models and tests
survB1 <- lmer(propsurv ~ 1 + (1 | allele:Line), data=survB)
survB2 <- lmer(propsurv ~ allele + (1 | allele:Line), data=survB)

MC.allele.Bonly <- PBmodcomp(survB2, survB1, nsim=1000)
MC.allele.Bonly

## D only models and tests
survD1 <- lmer(propsurv ~ 1 + (1 | allele:Line), data=survD)
survD2 <- lmer(propsurv ~ allele + (1 | allele:Line), data=survD)

MC.allele.Donly <- PBmodcomp(survD2, survD1, nsim=1000)
MC.allele.Donly




###########################################################################
                  #### sexspefic survival rates ####

                        #### sex ratio ####

# one of the assumptions made in this experiment is that the starting sex ratio is the same
# but the adult sex ratio may alter if there are diffrrences in the relative survival of the 2 sexes
# sexratio is recorded as the proportion of the survivng adults that are male (Larval.surv.chro$propmale)

par(mfrow=c(1,3))
hist(Larval.surv.chro$propmale)
hist(sqrt(Larval.surv.chro$propmale))
hist(log(Larval.surv.chro$propmale)) 
# normal is best
shapiro.test(Larval.surv.chro$propmale)
shapiro.test(sqrt(Larval.surv.chro$propmale))
shapiro.test(log(Larval.surv.chro$propmale))

## seesm fine to use propmale alone?

# since its a proportion porbably best to use glmer and binomial logit link
M.sexratio.1 <- lmer(propmale ~ allele*chro + (1 | allele:Line) + (1 | vial_line), data=Larval.surv.chro)

res <- resid(M.sexratio.1)
ran <- ranef(M.sexratio.1)
fit <- fitted(M.sexratio.1)
par(mfrow=c(1,4))
hist(ran$`allele:Line`[,1], main='')
hist(ran$vial_line[,1], main='')
plot(res ~ fit)
qqnorm(res, main='')
qqline(res)
summary(M.sexratio.1)
anova(M.sexratio.1)
Anova(M.sexratio.1)
## no affect of anything

## try sqrt??
M.sexratio.2 <- lmer(sqrt(propmale) ~ allele*chro + (1 | allele:Line) + (1 | vial_line), data=Larval.surv.chro)

res <- resid(M.sexratio.2)
ran <- ranef(M.sexratio.2)
fit <- fitted(M.sexratio.2)
par(mfrow=c(1,4))
hist(ran$`allele:Line`[,1], main='')
hist(ran$vial_line[,1], main='')
plot(res ~ fit)
qqnorm(res, main='')
qqline(res)
summary(M.sexratio.2)
anova(M.sexratio.2)
Anova(M.sexratio.2)

### again not any difference
anova(M.sexratio.1, M.sexratio.2)
## actually seems much better from this



## drop interaction term
M.sexratio.3 <- lmer(sqrt(propmale) ~ allele+chro + (1 | allele:Line) + (1 | vial_line), data=Larval.surv.chro)

res <- resid(M.sexratio.3)
ran <- ranef(M.sexratio.3)
fit <- fitted(M.sexratio.3)
par(mfrow=c(1,4))
hist(ran$`allele:Line`[,1], main='')
hist(ran$vial_line[,1], main='')
plot(res ~ fit)
qqnorm(res, main='')
qqline(res)
summary(M.sexratio.3)
anova(M.sexratio.3)
Anova(M.sexratio.3)

# no sex ratio difference due to fru allele or chromosomal compliment

SR1 <- lmer(sqrt(propmale) ~ 1 + (1 | allele:Line) + (1 | vial_line), data=Larval.surv.chro)
SR2 <- lmer(sqrt(propmale) ~ allele+ + (1 | allele:Line) + (1 | vial_line), data=Larval.surv.chro)
SR3 <- lmer(sqrt(propmale) ~ chro + (1 | allele:Line) + (1 | vial_line), data=Larval.surv.chro)
SR4 <- lmer(sqrt(propmale) ~ allele+chro + (1 | allele:Line) + (1 | vial_line), data=Larval.surv.chro)
SR5 <- lmer(sqrt(propmale) ~ allele*chro + (1 | allele:Line) + (1 | vial_line), data=Larval.surv.chro)

SR.allele <- PBmodcomp(SR4, SR3, nsim=1000)
SR.allele
SR.chro <- PBmodcomp(SR4, SR2, nsim=1000)
SR.chro
SR.inter <- PBmodcomp(SR5, SR4)
SR.inter


### extra tests splitting by balancer genotype 
## analyse within each of the these balancer genotypes
balancersplitSR <- split(Larval.surv.chro, Larval.surv.chro$chro)
SRB <- balancersplitSR[[1]]
SRD <- balancersplitSR[[2]]
str(SRB)
str(SRD)

## B only models and tests
SRB1 <- lmer(sqrt(propmale) ~ 1 + (1 | allele:Line), data=SRB)
SRB2 <- lmer(sqrt(propmale) ~ allele + (1 | allele:Line), data=SRB)

MC.allele.Bonly <- PBmodcomp(SRB2, SRB1, nsim=1000)
MC.allele.Bonly

## D only models and tests
SRD1 <- lmer(sqrt(propmale) ~ 1 + (1 | allele:Line), data=SRD)
SRD2 <- lmer(sqrt(propmale) ~ allele + (1 | allele:Line), data=SRD)

MC.allele.Donly <- PBmodcomp(SRD2, SRD1, nsim=1000)
MC.allele.Donly





            #### relative sex-specifc survival ####


## however what we really want to know is if the version of the fru allele affects the egg to adult survival of the two sexes differently 
## for this then we need to take into account the relative survival of both sexes in each line
## this is done by dividing the proportio of males who survived by the number of females that also survived.
Larval.surv.chro$sexrelsurv <- Larval.surv.chro$propexpmale/Larval.surv.chro$propexpfem
## this metric means that if sexrelsurv > 1 then more males survive relative to females, and if < 1 than more females survive realtive to males.

par(mfrow=c(1,3))
hist(Larval.surv.chro$sexrelsurv)
hist(log(Larval.surv.chro$sexrelsurv))
hist(sqrt(Larval.surv.chro$sexrelsurv))
# logging looks much better


# four zeros (rows 72, 81, 85, 90)  -remove as logging looks better and will probs be necessary for decent model fit
Larval.surv.chro.Edit <- Larval.surv.chro[-c(72, 81, 85, 90), ]
# no infinties in this data set
str(Larval.surv.chro.Edit) # seesm to have worked, 176 rows


## replot the data
par(mfrow=c(1,3))
hist(Larval.surv.chro.Edit$sexrelsurv)
hist(log(Larval.surv.chro.Edit$sexrelsurv))
hist(sqrt(Larval.surv.chro.Edit$sexrelsurv))
# logging looks much better
shapiro.test(log(Larval.surv.chro.Edit$sexrelsurv)) 

summarySE(Larval.surv.chro.Edit, measurevar=c("sexrelsurv"), groupvars=c("allele"))
# L = 1.129, S = 1.396

summarySE(Larval.surv.chro.Edit, measurevar=c("sexrelsurv"), groupvars=c("chro"))
# B = 1.222, D = 1.266

summarySE(Larval.surv.chro, measurevar=c("sexrelsurv"), groupvars=c("allelechro"))

summarySE(Larval.surv.chro.Edit, measurevar=c("sexrelsurv"), groupvars=c("Line"))
# L1 = 0.9443, L2 = 1.2965, L3 = 1.0858 S1 = 1.3725, S2 = 1.46885, S3 = 1.3714  
# only 14 S3!!!!! -issues


## models to test for difference in the relative survival for the two sexes - using edited data set
# first a simple model - no logging
M.relsexsurv.chro.1 <- lmer(sexrelsurv ~ allele*chro + (1 | allele:Line) + (1 | vial_line), data=Larval.surv.chro.Edit)
res <- resid(M.relsexsurv.chro.1)
ran <- ranef(M.relsexsurv.chro.1)
fit <- fitted(M.relsexsurv.chro.1)
par(mfrow=c(1,4))
hist(ran$`allele:Line`[,1], main='')
hist(ran$vial_line[,1], main='')
plot(res ~ fit)
qqnorm(res, main='')
qqline(res)
summary(M.relsexsurv.chro.1)
anova(M.relsexsurv.chro.1)
Anova(M.relsexsurv.chro.1)
## diagnostic plots look awful


M.relsexsurv.chro.2 <- lmer(log(sexrelsurv) ~ allele*chro + (1 | allele:Line) + (1 | vial_line), data=Larval.surv.chro.Edit)
res <- resid(M.relsexsurv.chro.2)
ran <- ranef(M.relsexsurv.chro.2)
fit <- fitted(M.relsexsurv.chro.2)
par(mfrow=c(1,4))
hist(ran$`allele:Line`[,1], main='')
hist(ran$vial_line[,1], main='')
plot(res ~ fit)
qqnorm(res, main='')
qqline(res)
summary(M.relsexsurv.chro.2)
anova(M.relsexsurv.chro.2)
Anova(M.relsexsurv.chro.2)
# diagnostic plots look much better 
# allele is now a borderline significant affect

## compare the two models 
anova(M.relsexsurv.chro.1, M.relsexsurv.chro.2)
## model two seems much better


# remove intercation term - see if allele effect p< 0.05 
M.relsexsurv.chro.3 <- lmer(log(sexrelsurv) ~ allele + chro + (1 | allele:Line) + (1 | vial_line), data=Larval.surv.chro.Edit)
res <- resid(M.relsexsurv.chro.3)
ran <- ranef(M.relsexsurv.chro.3)
fit <- fitted(M.relsexsurv.chro.3)
par(mfrow=c(1,4))
hist(ran$`allele:Line`[,1], main='')
hist(ran$vial_line[,1], main='')
plot(res ~ fit)
qqnorm(res, main='')
qqline(res)
summary(M.relsexsurv.chro.3)
anova(M.relsexsurv.chro.3)
Anova(M.relsexsurv.chro.3)
# again a borderline significant term

###
##### seems bad to throw away data
## add 0.001 tio every sexrelsurv entry - will now be able to log 

## create new column in original dataset
Larval.surv.chro$Add_sexrelsurv <- Larval.surv.chro$sexrelsurv + 0.001

## summary stats for plots
summarySE(Larval.surv.chro, measurevar=c("sexrelsurv"), groupvars=c("allele"))

summarySE(Larval.surv.chro, measurevar=c("sexrelsurv"), groupvars=c("chro"))

summarySE(Larval.surv.chro, measurevar=c("sexrelsurv"), groupvars=c("Line"))

summarySE(Larval.surv.chro, measurevar=c("sexrelsurv"), groupvars=c("allelechro"))


## can now do the models the same as before but using this altered full data set

M.relsexsurv.chro.4 <- lmer(Add_sexrelsurv ~ allele*chro + (1 | allele:Line) + (1 | vial_line), data=Larval.surv.chro, control = lmerControl(optimizer = "bobyqa"))
res <- resid(M.relsexsurv.chro.4)
ran <- ranef(M.relsexsurv.chro.4)
fit <- fitted(M.relsexsurv.chro.4)
par(mfrow=c(1,4))
hist(ran$`allele:Line`[,1], main='')
hist(ran$vial_line[,1], main='')
plot(res ~ fit)
qqnorm(res, main='')
qqline(res)
summary(M.relsexsurv.chro.4)
anova(M.relsexsurv.chro.4)
Anova(M.relsexsurv.chro.4)
## diagnostic plots look awful


M.relsexsurv.chro.5 <- lmer(log(Add_sexrelsurv) ~ allele*chro + (1 | allele:Line) + (1 | vial_line), data=Larval.surv.chro, control = lmerControl(optimizer = "bobyqa"))
res <- resid(M.relsexsurv.chro.5)
ran <- ranef(M.relsexsurv.chro.5)
fit <- fitted(M.relsexsurv.chro.5)
par(mfrow=c(1,4))
hist(ran$`allele:Line`[,1], main='')
hist(ran$vial_line[,1], main='')
plot(res ~ fit)
qqnorm(res, main='')
qqline(res)
summary(M.relsexsurv.chro.5)
anova(M.relsexsurv.chro.5)
Anova(M.relsexsurv.chro.5)
# diagnostic plots look much better 
# interaction and chro are now mariginal affects

anova(M.relsexsurv.chro.4, M.relsexsurv.chro.5) 
## when this doen the first time there is a warning message involving 'bobyqa -- a trust region step failed to reduce q'
## solution to this to add 'control = lmerControl(optimizer = "bobyqa")' to each of the models
## this didn't help but keep them in for now

## interaction is almost significant - check direction
summarySE(Larval.surv.chro, measurevar=c("Add_sexrelsurv"), groupvars=c("allelechro"))
# L.B = 1.0675, S.B = 1.294, L.D = 1.1923, S.D = 1.3606


## parametrci bootstrapping for p-vaues
RSS1 <- lmer(log(Add_sexrelsurv) ~ 1 + (1 | allele:Line) + (1 | vial_line), data=Larval.surv.chro, control = lmerControl(optimizer = "bobyqa"))
RSS2 <- lmer(log(Add_sexrelsurv) ~ allele + (1 | allele:Line) + (1 | vial_line), data=Larval.surv.chro, control = lmerControl(optimizer = "bobyqa"))
RSS3 <- lmer(log(Add_sexrelsurv) ~ chro + (1 | allele:Line) + (1 | vial_line), data=Larval.surv.chro, control = lmerControl(optimizer = "bobyqa"))
RSS4 <- lmer(log(Add_sexrelsurv) ~ allele+chro + (1 | allele:Line) + (1 | vial_line), data=Larval.surv.chro, control = lmerControl(optimizer = "bobyqa"))
RSS5 <- lmer(log(Add_sexrelsurv) ~ allele*chro + (1 | allele:Line) + (1 | vial_line), data=Larval.surv.chro, control = lmerControl(optimizer = "bobyqa"))

RSS.allele <- PBmodcomp(RSS4, RSS3)
RSS.allele
RSS.chro <- PBmodcomp(RSS4, RSS2)
RSS.chro
RSS.inter <- PBmodcomp(RSS5, RSS4)
RSS.inter




#############################################
##### PLOTS for larval survival #####


## plot of lines split by chro for total surviving offspring
dat_vlinesS <- data.frame(chro=c("B", "D"), xval=c(7.775, 9.975))
dat_vlinesL <- data.frame(chro=c("B", "D"), xval=c(12.22, 14.62))

dfS <- data.frame(chro=c("B", "D"), x1=c(3.6, 3.6), x2 = c(6.5, 6.5), y1=c(7.775, 9.975), y2=c(7.775, 9.975))
dfL <- data.frame(chro=c("B", "D"), x1=c(0.5, 0.5), x2 = c(3.4, 3.4), y1=c(12.22, 14.62), y2=c(12.22, 14.62))

total.surving.offspring_points <- ggplot(aes(y=total, x=Line,fill=allele), data=Larval.surv.chro)+
  geom_boxplot(outlier.shape = NA)+
  geom_point(aes(fill = allele), size = 2, shape = 21, position = position_jitterdodge(0.8))+
  geom_segment(aes(x=x1, y=y1, xend=x2, yend=y2), data=dfS, inherit.aes = FALSE, linetype=2, size=1.5)+
  geom_segment(aes(x=x1, y=y1, xend=x2, yend=y2), data=dfL, inherit.aes = FALSE, linetype=2, size=1.5)+
  ylab("No. of offspring surving to adulthood")+
  xlab("Line")+
  theme_bw()+
  theme(axis.text = element_text(size=30),
        axis.title.x = element_text(size=45), axis.title.y = element_text(size=45), legend.position = "none",strip.text=element_text(size=30))+
  theme(axis.title.x = element_text(vjust=-0.4), axis.title.y = element_text(vjus =2.5, hjust = 0.5))+
  theme(plot.margin = unit(c(0.5,0.1,0.7,0.7), "cm"))+
  theme(legend.position="none")+
  theme(legend.position="none")+
  scale_fill_manual(values=c("#CA3542", "#AECBC9"))+
  facet_wrap(~chro,scales="free_x")
total.surving.offspring_points
ggsave("Plots/total_offspring_surving_points.png")


## alternative verison using comments from reviewers
total.surving.offspring_points_alt <- ggplot(aes(y=total, x=Line, col=allele), data=Larval.surv.chro)+
  geom_boxplot(outlier.shape = NA)+
  geom_point(aes(col = allele), size = 2, shape = 21, position = position_jitterdodge(0.8))+
  geom_segment(aes(x=x1, y=y1, xend=x2, yend=y2), data=dfS, inherit.aes = FALSE, linetype=2, size=1.5)+
  geom_segment(aes(x=x1, y=y1, xend=x2, yend=y2), data=dfL, inherit.aes = FALSE, linetype=2, size=1.5)+
  ylab("No. of offspring surving to adulthood")+
  xlab("Line")+
  theme_bw()+
  theme(axis.text = element_text(size=30),
        axis.title.x = element_text(size=45), axis.title.y = element_text(size=45), legend.position = "none",strip.text=element_text(size=30))+
  theme(axis.title.x = element_text(vjust=-0.4), axis.title.y = element_text(vjus =2.5, hjust = 0.5))+
  theme(plot.margin = unit(c(0.5,0.1,0.7,0.7), "cm"))+
  theme(legend.position="none")+
  theme(legend.position="none")+
  scale_color_manual(values=c("#CA3542", "#7A8ED9"))+
  facet_wrap(~chro,scales="free_x")
total.surving.offspring_points_alt
ggsave("Plots/total.surving.offspring_points_alt.png")

### fix x axis from comments 

##/S
levels(Larval.surv.chro$Line)
levels(Larval.surv.chro$Line) <- c("L1/S","L2/S","L3/S", "S1/S", "S2/S", "S3/S")
levels(Larval.surv.chro$Line)
# NOW PLOT
total.surving.offspring_points_XS <- ggplot(aes(y=total, x=Line, col=allele), data=Larval.surv.chro)+
  geom_boxplot(outlier.shape = NA)+
  geom_point(aes(col = allele), size = 2, shape = 21, position = position_jitterdodge(0.8))+
  geom_segment(aes(x=x1, y=y1, xend=x2, yend=y2), data=dfS, inherit.aes = FALSE, linetype=2, size=1.5)+
  geom_segment(aes(x=x1, y=y1, xend=x2, yend=y2), data=dfL, inherit.aes = FALSE, linetype=2, size=1.5)+
  ylab("No. of offspring surving to adulthood")+
  xlab("Line")+
  theme_bw()+
  theme(axis.text = element_text(size=30),
        axis.title.x = element_text(size=45), axis.title.y = element_text(size=45), legend.position = "none",strip.text=element_text(size=30))+
  theme(axis.title.x = element_text(vjust=-0.4), axis.title.y = element_text(vjus =2.5, hjust = 0.5))+
  theme(plot.margin = unit(c(0.5,0.1,0.7,0.7), "cm"))+
  theme(legend.position="none")+
  theme(legend.position="none")+
  scale_color_manual(values=c("#CA3542", "#7A8ED9"))+
  facet_wrap(~chro,scales="free_x")
total.surving.offspring_points_XS
ggsave("Plots/total.surving.offspring_points_XS.png")

## /-
levels(Larval.surv.chro$Line)
levels(Larval.surv.chro$Line) <- c("L1/-","L2/-","L3/-", "S1/-", "S2/-", "S3/-")
levels(Larval.surv.chro$Line)
## AGAIN PLOT
total.surving.offspring_points_Xdash <- ggplot(aes(y=total, x=Line, col=allele), data=Larval.surv.chro)+
  geom_boxplot(outlier.shape = NA)+
  geom_point(aes(col = allele), size = 2, shape = 21, position = position_jitterdodge(0.8))+
  geom_segment(aes(x=x1, y=y1, xend=x2, yend=y2), data=dfS, inherit.aes = FALSE, linetype=2, size=1.5)+
  geom_segment(aes(x=x1, y=y1, xend=x2, yend=y2), data=dfL, inherit.aes = FALSE, linetype=2, size=1.5)+
  ylab("No. of offspring surving to adulthood")+
  xlab("Line")+
  theme_bw()+
  theme(axis.text = element_text(size=30),
        axis.title.x = element_text(size=45), axis.title.y = element_text(size=45), legend.position = "none",strip.text=element_text(size=30))+
  theme(axis.title.x = element_text(vjust=-0.4), axis.title.y = element_text(vjus =2.5, hjust = 0.5))+
  theme(plot.margin = unit(c(0.5,0.1,0.7,0.7), "cm"))+
  theme(legend.position="none")+
  theme(legend.position="none")+
  scale_color_manual(values=c("#CA3542", "#7A8ED9"))+
  facet_wrap(~chro,scales="free_x")
total.surving.offspring_points_Xdash
ggsave("Plots/total.surving.offspring_points_Xdash.png")


## plot for proportion surviving
dat_vlinesS <- data.frame(chro=c("B", "D"), xval=c(0.311, 0.399))
dat_vlinesL <- data.frame(chro=c("B", "D"), xval=c(0.4888, 0.5848))


proportion.surving.offspring <- ggplot(aes(y=propsurv, x=Line,fill=allele), data=Larval.surv.chro)+
  geom_boxplot(outlier.shape = NA)+
  geom_hline(aes(yintercept=xval), data=dat_vlinesS, size=1.5, linetype=2)+
  geom_hline(aes(yintercept=xval), data=dat_vlinesL, size=1.5, linetype=3)+
  ylab("Proportion of eggs reaching adulthood")+
  xlab("Line")+
  theme_bw()+
  theme(legend.title=element_blank(),axis.text = element_text(size=30),
        axis.title = element_text(size=40),legend.position = "none",strip.text=element_text(size=30))+
  theme(legend.position="none")+
  scale_fill_manual(values=c("#CA3542", "#AECBC9"))+
  facet_wrap(~chro,scales="free_x")
proportion.surving.offspring
ggsave("Plots/proportion_offspring_surviving.png")

## plot for sex ratio
# dat_vlinesS <- data.frame(chro=c("B", "D"), xval=c(0.466, 0.523))
# dat_vlinesL <- data.frame(chro=c("B", "D"), xval=c(0.476, 0.476))

dfS <- data.frame(chro=c("B", "D"), x1=c(3.6, 3.6), x2 = c(6.5, 6.5), y1=c(0.466, 0.523), y2=c(0.466, 0.523))
dfL <- data.frame(chro=c("B", "D"), x1=c(0.5, 0.5), x2 = c(3.4, 3.4), y1=c(0.476, 0.476), y2=c(0.476, 0.476))

sex.ratio_points <- ggplot(aes(y=propmale, x=Line,fill=allele), data=Larval.surv.chro)+
  geom_boxplot(outlier.shape = NA)+
  geom_point(aes(fill = allele), size = 2, shape = 21, position = position_jitterdodge(0.8))+
  geom_segment(aes(x=x1, y=y1, xend=x2, yend=y2), data=dfS, inherit.aes = FALSE, linetype=2, size=1.5)+
  geom_segment(aes(x=x1, y=y1, xend=x2, yend=y2), data=dfL, inherit.aes = FALSE, linetype=2, size=1.5)+
  ylab("Sex ratio")+
  xlab("Line")+
  theme_bw()+
  theme(axis.text = element_text(size=30),
        axis.title.x = element_text(size=45), axis.title.y = element_text(size=45), legend.position = "none",strip.text=element_text(size=30))+
  theme(axis.title.x = element_text(vjust=-0.4), axis.title.y = element_text(vjus =2.5, hjust = 0.5))+
  theme(plot.margin = unit(c(0.5,0.1,0.7,0.7), "cm"))+
  theme(legend.position="none")+
  scale_fill_manual(values=c("#CA3542", "#AECBC9"))+
  facet_wrap(~chro,scales="free_x")
sex.ratio_points
ggsave("Plots/sex_ratio_surviving_offspring_points.png")

## alternatiove version based on comments from the reviewers
sex.ratio_points_alt <- ggplot(aes(y=propmale, x=Line,col=allele), data=Larval.surv.chro)+
  geom_boxplot(outlier.shape = NA)+
  geom_point(aes(col = allele), size = 2, shape = 21, position = position_jitterdodge(0.8))+
  geom_segment(aes(x=x1, y=y1, xend=x2, yend=y2), data=dfS, inherit.aes = FALSE, linetype=2, size=1.5)+
  geom_segment(aes(x=x1, y=y1, xend=x2, yend=y2), data=dfL, inherit.aes = FALSE, linetype=2, size=1.5)+
  ylab("Sex ratio")+
  xlab("Line")+
  theme_bw()+
  theme(axis.text = element_text(size=30),
        axis.title.x = element_text(size=45), axis.title.y = element_text(size=45), legend.position = "none",strip.text=element_text(size=30))+
  theme(axis.title.x = element_text(vjust=-0.4), axis.title.y = element_text(vjus =2.5, hjust = 0.5))+
  theme(plot.margin = unit(c(0.5,0.1,0.7,0.7), "cm"))+
  theme(legend.position="none")+
  scale_color_manual(values=c("#CA3542", "#7A8ED9"))+
  facet_wrap(~chro,scales="free_x")
sex.ratio_points_alt
ggsave("Plots/sex_ratio_surviving_offspring_points_alt.png")


### alternative x axis

##/S
levels(Larval.surv.chro$Line)
levels(Larval.surv.chro$Line) <- c("L1/S","L2/S","L3/S", "S1/S", "S2/S", "S3/S")
levels(Larval.surv.chro$Line)
## plot
sex.ratio_points_XS <- ggplot(aes(y=propmale, x=Line,col=allele), data=Larval.surv.chro)+
  geom_boxplot(outlier.shape = NA)+
  geom_point(aes(col = allele), size = 2, shape = 21, position = position_jitterdodge(0.8))+
  geom_segment(aes(x=x1, y=y1, xend=x2, yend=y2), data=dfS, inherit.aes = FALSE, linetype=2, size=1.5)+
  geom_segment(aes(x=x1, y=y1, xend=x2, yend=y2), data=dfL, inherit.aes = FALSE, linetype=2, size=1.5)+
  ylab("Sex ratio")+
  xlab("Line")+
  theme_bw()+
  theme(axis.text = element_text(size=30),
        axis.title.x = element_text(size=45), axis.title.y = element_text(size=45), legend.position = "none",strip.text=element_text(size=30))+
  theme(axis.title.x = element_text(vjust=-0.4), axis.title.y = element_text(vjus =2.5, hjust = 0.5))+
  theme(plot.margin = unit(c(0.5,0.1,0.7,0.7), "cm"))+
  theme(legend.position="none")+
  scale_color_manual(values=c("#CA3542", "#7A8ED9"))+
  facet_wrap(~chro,scales="free_x")
sex.ratio_points_XS

##/-
levels(Larval.surv.chro$Line)
levels(Larval.surv.chro$Line) <- c("L1/-","L2/-","L3/-", "S1/-", "S2/-", "S3/-")
levels(Larval.surv.chro$Line)
#plot
sex.ratio_points_dash <- ggplot(aes(y=propmale, x=Line,col=allele), data=Larval.surv.chro)+
  geom_boxplot(outlier.shape = NA)+
  geom_point(aes(col = allele), size = 2, shape = 21, position = position_jitterdodge(0.8))+
  geom_segment(aes(x=x1, y=y1, xend=x2, yend=y2), data=dfS, inherit.aes = FALSE, linetype=2, size=1.5)+
  geom_segment(aes(x=x1, y=y1, xend=x2, yend=y2), data=dfL, inherit.aes = FALSE, linetype=2, size=1.5)+
  ylab("Sex ratio")+
  xlab("Line")+
  theme_bw()+
  theme(axis.text = element_text(size=30),
        axis.title.x = element_text(size=45), axis.title.y = element_text(size=45), legend.position = "none",strip.text=element_text(size=30))+
  theme(axis.title.x = element_text(vjust=-0.4), axis.title.y = element_text(vjus =2.5, hjust = 0.5))+
  theme(plot.margin = unit(c(0.5,0.1,0.7,0.7), "cm"))+
  theme(legend.position="none")+
  scale_color_manual(values=c("#CA3542", "#7A8ED9"))+
  facet_wrap(~chro,scales="free_x")
sex.ratio_points_dash



## relative survival of the sexes
dat_vlinesS <- data.frame(chro=c("B", "D"), xval=c(1.294, 1.359))
dat_vlinesL <- data.frame(chro=c("B", "D"), xval=c(1.067, 1.192))

relative.sex.surv <- ggplot(aes(y=sexrelsurv, x=Line,fill=allele), data=Larval.surv.chro)+
  geom_boxplot(outlier.shape = NA)+
  geom_hline(aes(yintercept=xval), data=dat_vlinesS, size=1.5, linetype=2)+
  geom_hline(aes(yintercept=xval), data=dat_vlinesL, size=1.5, linetype=3)+
  ylab("Relative survival of the sexes")+
  xlab("Line")+
  theme_bw()+
  theme(legend.title=element_blank(),axis.text = element_text(size=25),
        axis.title = element_text(size=25),legend.position = "none",strip.text=element_text(size=25))+
  theme(legend.position="none")+
  coord_cartesian(ylim =c(0, 3))+
  scale_fill_manual(values=c("#CA3542", "#AECBC9"))+
  facet_wrap(~chro,scales="free_x")
relative.sex.surv
ggsave("Plots/relative_survival_of_each_sex.png")

######################################################################################################
                  #### DEVELOPMENT TIME ####
# read it in
setwd("/Users/michaeljardine/Desktop/DTP/Data and analysis/fruitless/larvalsurvivalsummer2019/block3_chro_comp")
devtime <- read.csv("expanded_development_time_data.csv")
head(devtime)
summary(devtime)
str(devtime)

## first thing to make vial a factor, unique to each vial inthe experiment
devtime$line_chro <- interaction(devtime$lines, devtime$chro)
devtime$vial <- interaction(devtime$line_chro, devtime$vial)

str(devtime)

levels(devtime$vial)

### distribution of data
par(mfrow=c(1,2))
hist(devtime$eclosion_date)
hist(log(devtime$eclosion_date))
# logging maybe better but not huge variation anyway

hist(devtime$total_from_vial)
hist(log(devtime$total_from_vial))
## looks good by its self

## caluctaing avergaes and errors 


summarySE(devtime, measurevar=c("eclosion_date"), groupvars=c("sex"))
# F = 10.26, M = 10.48

summarySE(devtime, measurevar=c("eclosion_date"), groupvars=c("lines"))
# L1 = 10.31, L2 = 10.08, L3 = 10.54, S1 = 10.6, S2 = 10.275, S3 = 10.34

summarySE(devtime, measurevar=c("eclosion_date"), groupvars=c("chro"))
# B = 10.22, D = 10.49

summarySE(devtime, measurevar=c("eclosion_date"), groupvars=c("allele"))
# L = 10.32, S = 10.46

## create intraction between allele and sex
devtime$allele_sex <- interaction(devtime$allele, devtime$sex)
summarySE(devtime, measurevar=c("eclosion_date"), groupvars=c("allele_sex"))
# L.F = 10.23, S.F = 10.33, L.M = 10.43, S.M = 10.59

# interaction for line and chro
devtime$allele_chro <- interaction(devtime$allele, devtime$chro)
summarySE(devtime, measurevar=c("eclosion_date"), groupvars=c("allele_chro"))
# L.B = 10.19, S.B = 10.3, L.D = 10.44, S.D = 10.58
# B quicker than D?

# interaction for sex and chro
devtime$sex_chro <- interaction(devtime$sex, devtime$chro)
summarySE(devtime, measurevar=c("eclosion_date"), groupvars=c("sex_chro"))

#### sex difference in evelopment time ####

## females usually develop quicker than males  -does that occur here
## control for the different lines and for variation in vials 
## and for the chromosome compliment 

dev.sex.M1 <- lmer(eclosion_date ~ sex + (1 | vial) + (1 | lines) + (1| chro), data = devtime)

res <- resid(dev.sex.M1)
ran <- ranef(dev.sex.M1)
fit <- fitted(dev.sex.M1)
par(mfrow=c(1,5))
hist(ran$vial[,1], main='')
hist(ran$lines[,1], main='')
hist(ran$chro[,1], main='')
plot(res ~ fit)
qqnorm(res, main='')
qqline(res)
summary(dev.sex.M1)
anova(dev.sex.M1)
Anova(dev.sex.M1)
## despite the small difference - there is support for a differenec between males and females 

## but was better before when logged
dev.sex.M2 <- lmer(log(eclosion_date) ~ sex + (1 | vial) + (1 | lines) + (1| chro), data = devtime)

res <- resid(dev.sex.M2)
ran <- ranef(dev.sex.M2)
fit <- fitted(dev.sex.M2)
par(mfrow=c(1,5))
hist(ran$vial[,1], main='')
hist(ran$lines[,1], main='')
hist(ran$chro[,1], main='')
plot(res ~ fit)
qqnorm(res, main='')
qqline(res)
summary(dev.sex.M2)
anova(dev.sex.M2)
Anova(dev.sex.M2)
## result is the same and difficult to see from the diagnostic plots

anova(dev.sex.M1, dev.sex.M2)

### max said that a poisson distribution might actually be better for the model
dev.sex.M3 <- glmer(eclosion_date ~ sex + (1 | vial) + (1 | lines) + (1| chro), data = devtime, family = poisson(link = log))

res <- resid(dev.sex.M3)
ran <- ranef(dev.sex.M3)
fit <- fitted(dev.sex.M3)
par(mfrow=c(1,5))
hist(ran$vial[,1], main='')
hist(ran$lines[,1], main='')
hist(ran$chro[,1], main='')
plot(res ~ fit)
qqnorm(res, main='')
qqline(res)
summary(dev.sex.M3)
anova(dev.sex.M3)
Anova(dev.sex.M3)
# no longer a significant effect for sex
# check the two models

anova(dev.sex.M2, dev.sex.M3)
## logging better than poisson

#### an effect of fru on sex-specific development time? ####

dev.fru.sex.M1 <- lmer(log(eclosion_date) ~ allele*sex*chro + (1 | vial) + (1 | lines), data = devtime)

res <- resid(dev.fru.sex.M1)
ran <- ranef(dev.fru.sex.M1)
fit <- fitted(dev.fru.sex.M1)
par(mfrow=c(1,4))
hist(ran$vial[,1], main='')
hist(ran$lines[,1], main='')
plot(res ~ fit)
qqnorm(res, main='')
qqline(res)
summary(dev.fru.sex.M1)
anova(dev.fru.sex.M1)
Anova(dev.fru.sex.M1)
## drop 3 way interaction

### first lets tests this with the poisson distribution
dev.fru.sex.M1.1 <- glmer(eclosion_date ~ allele*sex*chro + (1 | vial) + (1 | lines), data = devtime, family = poisson(link = log))

res <- resid(dev.fru.sex.M1.1)
ran <- ranef(dev.fru.sex.M1.1)
fit <- fitted(dev.fru.sex.M1.1)
par(mfrow=c(1,4))
hist(ran$vial[,1], main='')
hist(ran$lines[,1], main='')
plot(res ~ fit)
qqnorm(res, main='')
qqline(res)
summary(dev.fru.sex.M1.1)
anova(dev.fru.sex.M1.1)
Anova(dev.fru.sex.M1.1)

# compare those two full models
anova(dev.fru.sex.M1, dev.fru.sex.M1.1)
## the first model using lmer and a logged response variable is better - keep this form of the model. 

dev.fru.sex.M2 <- lmer(log(eclosion_date) ~ allele*sex + allele*chro + sex*chro + (1 | vial) + (1 | lines), data = devtime)

res <- resid(dev.fru.sex.M2)
ran <- ranef(dev.fru.sex.M2)
fit <- fitted(dev.fru.sex.M2)
par(mfrow=c(1,4))
hist(ran$vial[,1], main='')
hist(ran$lines[,1], main='')
plot(res ~ fit)
qqnorm(res, main='')
qqline(res)
summary(dev.fru.sex.M2)
anova(dev.fru.sex.M2)
Anova(dev.fru.sex.M2)
## drop allele by chro

dev.fru.sex.M3 <- lmer(log(eclosion_date) ~ allele*sex + sex*chro + (1 | vial) + (1 | lines), data = devtime)

res <- resid(dev.fru.sex.M3)
ran <- ranef(dev.fru.sex.M3)
fit <- fitted(dev.fru.sex.M3)
par(mfrow=c(1,4))
hist(ran$vial[,1], main='')
hist(ran$lines[,1], main='')
plot(res ~ fit)
qqnorm(res, main='')
qqline(res)
summary(dev.fru.sex.M3)
anova(dev.fru.sex.M3)
Anova(dev.fru.sex.M3)
## drop allele by sex

dev.fru.sex.M4 <- lmer(log(eclosion_date) ~ sex*chro + allele + (1 | vial) + (1 | lines), data = devtime)

res <- resid(dev.fru.sex.M4)
ran <- ranef(dev.fru.sex.M4)
fit <- fitted(dev.fru.sex.M4)
par(mfrow=c(1,4))
hist(ran$vial[,1], main='')
hist(ran$lines[,1], main='')
plot(res ~ fit)
qqnorm(res, main='')
qqline(res)
summary(dev.fru.sex.M4)
anova(dev.fru.sex.M4)
Anova(dev.fru.sex.M4)
## doesn't change the model much and M3 had the interaction of interest. 

anova(dev.fru.sex.M1, dev.fru.sex.M2, dev.fru.sex.M3, dev.fru.sex.M4)


### bootsrapping for p-values
## model including factor and mostcomplicated model without factor

## allele
DT90 <- lmer(log(eclosion_date) ~ sex*chro + (1 | vial) + (1 | lines), data = devtime)
DT91 <- lmer(log(eclosion_date) ~ allele + sex*chro + (1 | vial) + (1 | lines), data = devtime)
# sex
DT92 <- lmer(log(eclosion_date) ~ allele*chro + (1 | vial) + (1 | lines), data = devtime)
DT93 <- lmer(log(eclosion_date) ~ allele*chro + sex + (1 | vial) + (1 | lines), data = devtime)
# chro
DT94 <- lmer(log(eclosion_date) ~ allele*sex + (1 | vial) + (1 | lines), data = devtime)
DT95 <- lmer(log(eclosion_date) ~ allele*sex + chro + (1 | vial) + (1 | lines), data = devtime)
# inetractions
DT96 <- lmer(log(eclosion_date) ~ allele*chro + sex*chro + (1 | vial) + (1 | lines), data = devtime)
DT97 <- lmer(log(eclosion_date) ~ allele*sex +  sex*chro + (1 | vial) + (1 | lines), data = devtime)
DT98 <- lmer(log(eclosion_date) ~ allele*sex + allele*chro + (1 | vial) + (1 | lines), data = devtime)
DT99 <- lmer(log(eclosion_date) ~ allele*sex + allele*chro + sex*chro + (1 | vial) + (1 | lines), data = devtime)
DT100 <- lmer(log(eclosion_date) ~ allele*sex*chro + (1 | vial) + (1 | lines), data = devtime)


DevT.A <- PBmodcomp(DT91, DT90, nsim=1000)
DevT.A
DevT.S <- PBmodcomp(DT93, DT92, nsim=1000)
DevT.S
DevT.C <- PBmodcomp(DT95, DT94, nsim=1000)
DevT.C
DevT.AS <- PBmodcomp(DT99, DT96, nsim=1000)
DevT.AS
DevT.AC <- PBmodcomp(DT99, DT97, nsim=1000)
DevT.AC
DevT.SC <- PBmodcomp(DT99, DT98, nsim=1000)
DevT.SC
DevT.3way <- PBmodcomp(DT100, DT99, nsim=1000)
DevT.3way

### extra tests splitting by balancer genotype 
## analyse within each of the these balancer genotypes
balancersplitDT <- split(devtime, devtime$chro)
DTB <- balancersplitDT[[1]]
DTD <- balancersplitDT[[2]]
str(DTB)
str(DTD)

## B only models and tests
DTB1 <- lmer(log(eclosion_date) ~ allele + (1 | vial) + (1 | lines), data = DTB)
DTB2 <- lmer(log(eclosion_date) ~ sex + (1 | vial) + (1 | lines), data = DTB)
DTB3 <- lmer(log(eclosion_date) ~ allele + sex + (1 | vial) + (1 | lines), data = DTB)
DTB4 <- lmer(log(eclosion_date) ~ allele*sex + (1 | vial) + (1 | lines), data = DTB)

MC.sex.Bonly <- PBmodcomp(DTB3, DTB1, nsim=1000)
MC.sex.Bonly
MC.allele.Bonly <- PBmodcomp(DTB3, DTB2, nsim=1000)
MC.allele.Bonly
MC.inter.Bonly <- PBmodcomp(DTB4, DTB3, nsim=1000)
MC.inter.Bonly


## D only models and tests
DTD1 <- lmer(log(eclosion_date) ~ allele + (1 | vial) + (1 | lines), data = DTD)
DTD2 <- lmer(log(eclosion_date) ~ sex + (1 | vial) + (1 | lines), data = DTD)
DTD3 <- lmer(log(eclosion_date) ~ allele + sex + (1 | vial) + (1 | lines), data = DTD)
DTD4 <- lmer(log(eclosion_date) ~ allele*sex + (1 | vial) + (1 | lines), data = DTD)

MC.sex.Donly <- PBmodcomp(DTD3, DTD1, nsim=1000)
MC.sex.Donly
MC.allele.Donly <- PBmodcomp(DTD3, DTD2, nsim=1000)
MC.allele.Donly
MC.inter.Donly <- PBmodcomp(DTD4, DTD3, nsim=1000)
MC.inter.Donly



#####################
### so there is no affect of fru allele on the development time of these flies

#### relationship between development time and number of fleis from vial ####

# is there a correlation?
cor(devtime$eclosion_date, devtime$total_from_vial)
## very weak - worth looking at - maybe 

# model
dev.total.M1 <- lmer(log(eclosion_date) ~ total_from_vial + (1 | vial) + (1 | lines) + (1| chro), data = devtime)

res <- resid(dev.total.M1)
ran <- ranef(dev.total.M1)
fit <- fitted(dev.total.M1)
par(mfrow=c(1,5))
hist(ran$vial[,1], main='')
hist(ran$lines[,1], main='')
hist(ran$chro[,1], main='')
plot(res ~ fit)
qqnorm(res, main='')
qqline(res)
summary(dev.total.M1)
anova(dev.total.M1)
Anova(dev.total.M1)

# include interaction with sex since this is so important

dev.total.M2 <- lmer(log(eclosion_date) ~ total_from_vial*sex + (1 | vial) + (1 | lines) + (1| chro), data = devtime)

res <- resid(dev.total.M2)
ran <- ranef(dev.total.M2)
fit <- fitted(dev.total.M2)
par(mfrow=c(1,5))
hist(ran$vial[,1], main='')
hist(ran$lines[,1], main='')
hist(ran$chro[,1], main='')
plot(res ~ fit)
qqnorm(res, main='')
qqline(res)
summary(dev.total.M2)
anova(dev.total.M2)
Anova(dev.total.M2)

## nope- sex is significant as it was before - nothing else is 

## bootsrapping for p-values
vial.total.1 <- lmer(log(eclosion_date) ~ 1 + (1 | vial) + (1 | lines) + (1| chro), data = devtime)
vial.total.2 <- lmer(log(eclosion_date) ~ total_from_vial + (1 | vial) + (1 | lines) + (1| chro), data = devtime)
vial.total.3 <- lmer(log(eclosion_date) ~ sex + (1 | vial) + (1 | lines) + (1| chro), data = devtime)
vial.total.4 <- lmer(log(eclosion_date) ~ total_from_vial+sex + (1 | vial) + (1 | lines) + (1| chro), data = devtime)
vial.total.5 <- lmer(log(eclosion_date) ~ total_from_vial*sex + (1 | vial) + (1 | lines) + (1| chro), data = devtime)

VT.total.from.vial <- PBmodcomp(vial.total.4, vial.total.3)
VT.total.from.vial
VT.sex <- PBmodcomp(vial.total.4, vial.total.2)
VT.sex
VT.inter <- PBmodcomp(vial.total.5, vial.total.4)
VT.inter

##################################################
##### plots #####

# point and line  plot would be better
line.plots <- summarySE(devtime, measurevar=c("eclosion_date"), groupvars=c("line_chro"))
line.plots$chro <- c(rep('B', 6), rep('D', 6))
line.plots$allele <- c(rep('L', 3), rep('S', 3), rep('L', 3), rep('S', 3))
line.plots$line <- c('L1', 'L2', 'L3', 'S1', 'S2', 'S3', 'L1', 'L2', 'L3', 'S1', 'S2', 'S3')

dfS <- data.frame(chro=c("B", "D"), x1=c(3.6, 3.6), x2 = c(6.5, 6.5), y1=c(10.304, 10.584), y2=c(10.304, 10.584))
dfL <- data.frame(chro=c("B", "D"), x1=c(0.5, 0.5), x2 = c(3.4, 3.4), y1=c(10.19, 10.47), y2=c(10.19, 10.47))

dev.time.line <- ggplot(aes(x=line, y=eclosion_date, color=allele), data=line.plots)+
  geom_point(size=8)+
  geom_errorbar(aes(ymin=(eclosion_date-se), ymax=(eclosion_date+se)), size=3.5)+
  geom_segment(aes(x=x1, y=y1, xend=x2, yend=y2), data=dfS, inherit.aes = FALSE, linetype=2, size=1.5)+
  geom_segment(aes(x=x1, y=y1, xend=x2, yend=y2), data=dfL, inherit.aes = FALSE, linetype=2, size=1.5)+
  ylab("Development time (days)")+
  xlab("Line")+
  theme_bw()+
  theme(axis.text = element_text(size=30),
        axis.title.x = element_text(size=45), axis.title.y = element_text(size=45), legend.position = "none",strip.text=element_text(size=30))+
  theme(axis.title.x = element_text(vjust=-0.4), axis.title.y = element_text(vjus =2.5, hjust = 0.5))+
  theme(plot.margin = unit(c(0.5,0.1,0.7,0.7), "cm"))+
  scale_color_manual(values=c("#CA3542", "#AECBC9"))+
  facet_wrap(~chro,scales="free_x")
dev.time.line
ggsave("Plots/development_time_line_point_plot.png")

## suggestion to produce 2 plots split by sex since we know that sex is a major factor in development time
split.by.sex <- split(devtime, devtime$sex)
devtime.Male <- split.by.sex[[2]]
devtime.Female <- split.by.sex[[1]]

### male
line.plots.M <- summarySE(devtime.Male, measurevar=c("eclosion_date"), groupvars=c("line_chro"))
line.plots.M$chro <- c(rep('B', 6), rep('D', 6))
line.plots.M$allele <- c(rep('L', 3), rep('S', 3), rep('L', 3), rep('S', 3))
line.plots.M$line <- c('L1', 'L2', 'L3', 'S1', 'S2', 'S3', 'L1', 'L2', 'L3', 'S1', 'S2', 'S3')

summarySE(devtime.Male, measurevar=c("eclosion_date"), groupvars=c("allele_chro"))
dfS <- data.frame(chro=c("B", "D"), x1=c(3.6, 3.6), x2 = c(6.5, 6.5), y1=c(10.39, 10.72), y2=c(10.39, 10.72))
dfL <- data.frame(chro=c("B", "D"), x1=c(0.5, 0.5), x2 = c(3.4, 3.4), y1=c(10.28, 10.55), y2=c(10.28, 10.55))

dev.time.line.Male <- ggplot(aes(x=line, y=eclosion_date, color=allele), data=line.plots.M)+
  geom_point(size=8)+
  geom_errorbar(aes(ymin=(eclosion_date-se), ymax=(eclosion_date+se)), size=3.5)+
  geom_segment(aes(x=x1, y=y1, xend=x2, yend=y2), data=dfS, inherit.aes = FALSE, linetype=2, size=1.5)+
  geom_segment(aes(x=x1, y=y1, xend=x2, yend=y2), data=dfL, inherit.aes = FALSE, linetype=2, size=1.5)+
  ylab("Development time (days)")+
  xlab("Line")+
  theme_bw()+
  theme(axis.text = element_text(size=30),
        axis.title.x = element_text(size=45), axis.title.y = element_text(size=45), legend.position = "none",strip.text=element_text(size=30))+
  theme(axis.title.x = element_text(vjust=-0.4), axis.title.y = element_text(vjus =2.5, hjust = 0.5))+
  theme(plot.margin = unit(c(0.5,0.1,0.7,0.7), "cm"))+
  scale_color_manual(values=c("#CA3542", "#7A8ED9"))+
  facet_wrap(~chro,scales="free_x")
dev.time.line.Male
ggsave("Plots/development_time_line_point_plot_male.png")



##/S
line.plots.M$line <- c('L1/S', 'L2/S', 'L3/S', 'S1/S', 'S2/S', 'S3/S', 'L1/S', 'L2/S', 'L3/S', 'S1/S', 'S2/S', 'S3/S')
## plot
dev.time.line.Male.XS <- ggplot(aes(x=line, y=eclosion_date, color=allele), data=line.plots.M)+
  geom_point(size=8)+
  geom_errorbar(aes(ymin=(eclosion_date-se), ymax=(eclosion_date+se)), size=3.5)+
  geom_segment(aes(x=x1, y=y1, xend=x2, yend=y2), data=dfS, inherit.aes = FALSE, linetype=2, size=1.5)+
  geom_segment(aes(x=x1, y=y1, xend=x2, yend=y2), data=dfL, inherit.aes = FALSE, linetype=2, size=1.5)+
  ylab("Development time (days)")+
  xlab("Line")+
  theme_bw()+
  theme(axis.text = element_text(size=30),
        axis.title.x = element_text(size=45), axis.title.y = element_text(size=45), legend.position = "none",strip.text=element_text(size=30))+
  theme(axis.title.x = element_text(vjust=-0.4), axis.title.y = element_text(vjus =2.5, hjust = 0.5))+
  theme(plot.margin = unit(c(0.5,0.1,0.7,0.7), "cm"))+
  scale_color_manual(values=c("#CA3542", "#7A8ED9"))+
  facet_wrap(~chro,scales="free_x")
dev.time.line.Male.XS

##/-
line.plots.M$line <- c('L1/-', 'L2/-', 'L3/-', 'S1/-', 'S2/-', 'S3/-', 'L1/-', 'L2/-', 'L3/-', 'S1/-', 'S2/-', 'S3/-')
## plot
dev.time.line.Male.Xdash <- ggplot(aes(x=line, y=eclosion_date, color=allele), data=line.plots.M)+
  geom_point(size=8)+
  geom_errorbar(aes(ymin=(eclosion_date-se), ymax=(eclosion_date+se)), size=3.5)+
  geom_segment(aes(x=x1, y=y1, xend=x2, yend=y2), data=dfS, inherit.aes = FALSE, linetype=2, size=1.5)+
  geom_segment(aes(x=x1, y=y1, xend=x2, yend=y2), data=dfL, inherit.aes = FALSE, linetype=2, size=1.5)+
  ylab("Development time (days)")+
  xlab("Line")+
  theme_bw()+
  theme(axis.text = element_text(size=30),
        axis.title.x = element_text(size=45), axis.title.y = element_text(size=45), legend.position = "none",strip.text=element_text(size=30))+
  theme(axis.title.x = element_text(vjust=-0.4), axis.title.y = element_text(vjus =2.5, hjust = 0.5))+
  theme(plot.margin = unit(c(0.5,0.1,0.7,0.7), "cm"))+
  scale_color_manual(values=c("#CA3542", "#7A8ED9"))+
  facet_wrap(~chro,scales="free_x")
dev.time.line.Male.Xdash

### female
line.plots.F <- summarySE(devtime.Female, measurevar=c("eclosion_date"), groupvars=c("line_chro"))
line.plots.F$chro <- c(rep('B', 6), rep('D', 6))
line.plots.F$allele <- c(rep('L', 3), rep('S', 3), rep('L', 3), rep('S', 3))
line.plots.F$line <- c('L1', 'L2', 'L3', 'S1', 'S2', 'S3', 'L1', 'L2', 'L3', 'S1', 'S2', 'S3')

summarySE(devtime.Female, measurevar=c("eclosion_date"), groupvars=c("allele_chro"))
dfS <- data.frame(chro=c("B", "D"), x1=c(3.6, 3.6), x2 = c(6.5, 6.5), y1=c(10.215, 10.43), y2=c(10.215, 10.43))
dfL <- data.frame(chro=c("B", "D"), x1=c(0.5, 0.5), x2 = c(3.4, 3.4), y1=c(10.106, 10.33), y2=c(10.106, 10.33))

dev.time.line.Female <- ggplot(aes(x=line, y=eclosion_date, color=allele), data=line.plots.F)+
  geom_point(size=8)+
  geom_errorbar(aes(ymin=(eclosion_date-se), ymax=(eclosion_date+se)), size=3.5)+
  geom_segment(aes(x=x1, y=y1, xend=x2, yend=y2), data=dfS, inherit.aes = FALSE, linetype=2, size=1.5)+
  geom_segment(aes(x=x1, y=y1, xend=x2, yend=y2), data=dfL, inherit.aes = FALSE, linetype=2, size=1.5)+
  ylab("Development time (days)")+
  xlab("Line")+
  theme_bw()+
  theme(axis.text = element_text(size=30),
        axis.title.x = element_text(size=45), axis.title.y = element_text(size=45), legend.position = "none",strip.text=element_text(size=30))+
  theme(axis.title.x = element_text(vjust=-0.4), axis.title.y = element_text(vjus =2.5, hjust = 0.5))+
  theme(plot.margin = unit(c(0.5,0.1,0.7,0.7), "cm"))+
  scale_color_manual(values=c("#CA3542", "#7A8ED9"))+
  facet_wrap(~chro,scales="free_x")
dev.time.line.Female
ggsave("Plots/development_time_line_point_plot_female.png")

##/S
line.plots.F$line <- c('L1/S', 'L2/S', 'L3/S', 'S1/S', 'S2/S', 'S3/S', 'L1/S', 'L2/S', 'L3/S', 'S1/S', 'S2/S', 'S3/S')
## plot
dev.time.line.Female.XS <- ggplot(aes(x=line, y=eclosion_date, color=allele), data=line.plots.F)+
  geom_point(size=8)+
  geom_errorbar(aes(ymin=(eclosion_date-se), ymax=(eclosion_date+se)), size=3.5)+
  geom_segment(aes(x=x1, y=y1, xend=x2, yend=y2), data=dfS, inherit.aes = FALSE, linetype=2, size=1.5)+
  geom_segment(aes(x=x1, y=y1, xend=x2, yend=y2), data=dfL, inherit.aes = FALSE, linetype=2, size=1.5)+
  ylab("Development time (days)")+
  xlab("Line")+
  theme_bw()+
  theme(axis.text = element_text(size=30),
        axis.title.x = element_text(size=45), axis.title.y = element_text(size=45), legend.position = "none",strip.text=element_text(size=30))+
  theme(axis.title.x = element_text(vjust=-0.4), axis.title.y = element_text(vjus =2.5, hjust = 0.5))+
  theme(plot.margin = unit(c(0.5,0.1,0.7,0.7), "cm"))+
  scale_color_manual(values=c("#CA3542", "#7A8ED9"))+
  facet_wrap(~chro,scales="free_x")
dev.time.line.Female.XS

##/-
line.plots.F$line <- c('L1/-', 'L2/-', 'L3/-', 'S1/-', 'S2/-', 'S3/-', 'L1/-', 'L2/-', 'L3/-', 'S1/-', 'S2/-', 'S3/-')
## plot
dev.time.line.Female.dash <- ggplot(aes(x=line, y=eclosion_date, color=allele), data=line.plots.F)+
  geom_point(size=8)+
  geom_errorbar(aes(ymin=(eclosion_date-se), ymax=(eclosion_date+se)), size=3.5)+
  geom_segment(aes(x=x1, y=y1, xend=x2, yend=y2), data=dfS, inherit.aes = FALSE, linetype=2, size=1.5)+
  geom_segment(aes(x=x1, y=y1, xend=x2, yend=y2), data=dfL, inherit.aes = FALSE, linetype=2, size=1.5)+
  ylab("Development time (days)")+
  xlab("Line")+
  theme_bw()+
  theme(axis.text = element_text(size=30),
        axis.title.x = element_text(size=45), axis.title.y = element_text(size=45), legend.position = "none",strip.text=element_text(size=30))+
  theme(axis.title.x = element_text(vjust=-0.4), axis.title.y = element_text(vjus =2.5, hjust = 0.5))+
  theme(plot.margin = unit(c(0.5,0.1,0.7,0.7), "cm"))+
  scale_color_manual(values=c("#CA3542", "#7A8ED9"))+
  facet_wrap(~chro,scales="free_x")
dev.time.line.Female.dash




### sex difference in devlopemnt time
dat_vlinesF <- data.frame(chro=c("B", "D"), xval=c(10.14, 10.363))
dat_vlinesM <- data.frame(chro=c("B", "D"), xval=c(10.323, 10.616))

sex_plots <- summarySE(devtime, measurevar=c("eclosion_date"), groupvars=c("sex_chro"))
sex_plots$chro <- c(rep('B', 2), rep('D', 2))
sex_plots$sex <- c('F', 'M', 'F', 'M')
sex_plots

dev.time.sex <- ggplot(aes(x=sex, y=eclosion_date, color=sex), data=sex_plots)+
  geom_point(size=8)+
  geom_errorbar(aes(ymin=(eclosion_date-se), ymax=(eclosion_date+se), size=3.5))+
  ylab("Development time (days)")+
  xlab("Sex")+
  theme_bw()+
  theme(axis.text = element_text(size=30),
        axis.title.x = element_text(size=45), axis.title.y = element_text(size=45), legend.position = "none",strip.text=element_text(size=30))+
  theme(axis.title.x = element_text(vjust=-0.4), axis.title.y = element_text(vjus =2.5, hjust = 0.5))+
  theme(plot.margin = unit(c(0.5,0.1,0.7,0.7), "cm"))+
  scale_colour_manual(values=c("#E69F00", "#009E73"))+
  facet_wrap(~chro,scales="free_x")
dev.time.sex
ggsave("Plots/development_time_sex_point_plot.png")



## create box plot version of the abive plot based on suggestions from reviewer 1
dev.time.sex.box <- ggplot(aes(y=eclosion_date, x=sex, fill=sex), data=devtime)+
  geom_boxplot(outlier.shape = NA)+
  geom_point(aes(fill = sex), size = 2, shape = 21, position = position_jitterdodge(0.8))+
  ylab("Development time (days)")+
  xlab("Sex")+
  theme_bw()+
  theme(axis.text = element_text(size=30),
        axis.title.x = element_text(size=45), axis.title.y = element_text(size=45), legend.position = "none",strip.text=element_text(size=30))+
  theme(axis.title.x = element_text(vjust=-0.4), axis.title.y = element_text(vjus =2.5, hjust = 0.5))+
  theme(plot.margin = unit(c(0.5,0.1,0.7,0.7), "cm"))+
  scale_fill_manual(values=c("#E69F00", "#009E73"))+
  facet_wrap(~chro,scales="free_x")
dev.time.sex.box
ggsave("Plots/development_time_sex_box_plot.png")
## doesn't look good
## discrete and limited spread of data conributes to uninformative plot hwne takin gthe same approach as other plots in the paper




## scatter plot development time with number of flies from each vial


time_and_total<- ggplot(aes(x=total_from_vial, y=eclosion_date), data=devtime) +
  geom_smooth(method=lm, formula = y~x, se=T, level=0.99, fullrange=T, colour="black") +
  geom_jitter(width = 0.2, height = 0.1, shape=1)+
  ylab("Development time(days)")+
  xlab("Total number of flies eclosing per vial")+
  theme_bw()+
  theme(legend.title=element_blank(),axis.text = element_text(size=25),
        axis.title = element_text(size=30), strip.text=element_text(size=25))
time_and_total
ggsave("Plots/dev_time_versus_number_per_vial.png")

mean(devtime$total_from_vial)

## same plot but with 1 point per vial
time_and_total_plots <- summarySE(devtime, measurevar=c("total_from_vial"), groupvars=c("vial"))
times <- summarySE(devtime, measurevar=c("eclosion_date"), groupvars=c("vial"))
time_and_total_plots$eclosion_date <- times$eclosion_date

time_and_total_1_per_vial<- ggplot(aes(x=total_from_vial, y=eclosion_date), data=time_and_total_plots) +
  geom_point(shape=1, col=1)+
  geom_smooth(method=lm, formula = y~x, se=T, level=0.99, fullrange=T, colour="black") +
  geom_jitter(width = 0.15, shape=1)+
  ylab("Development time(days)")+
  xlab("Total number of flies eclosing per vial")+
  theme_bw()+
  theme(legend.title=element_blank(),axis.text = element_text(size=25),
        axis.title = element_text(size=30), strip.text=element_text(size=25))
time_and_total_1_per_vial
ggsave("Plots/dev_time_versus_number_per_vial.1_point_per_vial.png")


plot.dat <- ggplot_build(time_and_total_1_per_vial)$data[[2]]
x1 = plot.dat[2,1]
x2 = plot.dat[12,1]
y1 = plot.dat[2,2]
y2 = plot.dat[12,2]
m = (y2-y1)/(x2-x1)
m



geom_violin(trim=FALSE)
#### violin plot
line.plots.M <- summarySE(devtime.Male, measurevar=c("eclosion_date"), groupvars=c("line_chro"))
line.plots.M$chro <- c(rep('B', 6), rep('D', 6))
line.plots.M$allele <- c(rep('L', 3), rep('S', 3), rep('L', 3), rep('S', 3))
line.plots.M$line <- c('L1', 'L2', 'L3', 'S1', 'S2', 'S3', 'L1', 'L2', 'L3', 'S1', 'S2', 'S3')

summarySE(devtime.Male, measurevar=c("eclosion_date"), groupvars=c("allele_chro"))
dfS <- data.frame(chro=c("B", "D"), x1=c(3.6, 3.6), x2 = c(6.5, 6.5), y1=c(10.39, 10.72), y2=c(10.39, 10.72))
dfL <- data.frame(chro=c("B", "D"), x1=c(0.5, 0.5), x2 = c(3.4, 3.4), y1=c(10.28, 10.55), y2=c(10.28, 10.55))

data_summary <- function(x) {
  m <- mean(x)
  ymin <- m-sd(x)
  ymax <- m+sd(x)
  return(c(y=m,ymin=ymin,ymax=ymax))
}

### change plots conpletely
##/S
levels(devtime.Male$lines)
devtime.Male$lines <- revalue(devtime.Male$lines, c("L1"="L1/S", "L2"="L2/S", "L3"="L3/S", "S1"="S1/S", "S2"="S2/S", "S3"="S3/S"))
levels(devtime.Male$lines)
## plot
stat_sum_scatter_male.S <- ggplot(aes(x=lines, y=eclosion_date, color=allele), data=devtime.Male)+
  geom_point(size=1)+
  geom_jitter(width=0.2, height=0.05)+
  stat_summary(fun.data=data_summary, lwd=1.7, size=2)+
  geom_segment(aes(x=x1, y=y1, xend=x2, yend=y2), data=dfS, inherit.aes = FALSE, linetype=2, size=1.5)+
  geom_segment(aes(x=x1, y=y1, xend=x2, yend=y2), data=dfL, inherit.aes = FALSE, linetype=2, size=1.5)+
  ylab("Development time (days)")+
  xlab("Line")+
  theme_bw()+
  theme(axis.text = element_text(size=30),
        axis.title.x = element_text(size=45), axis.title.y = element_text(size=45), legend.position = "none",strip.text=element_text(size=30))+
  theme(axis.title.x = element_text(vjust=-0.4), axis.title.y = element_text(vjus =2.5, hjust = 0.5))+
  theme(plot.margin = unit(c(0.5,0.1,0.7,0.7), "cm"))+
  scale_colour_manual(values=c("#CA3542", "#7A8ED9"))+
  facet_wrap(~chro,scales="free_x")
stat_sum_scatter_male.S

## reset
split.by.sex <- split(devtime, devtime$sex)
devtime.Male <- split.by.sex[[2]]
devtime.Female <- split.by.sex[[1]]

##/-
levels(devtime.Male$lines)
devtime.Male$lines <- revalue(devtime.Male$lines, c("L1"="L1/-", "L2"="L2/-", "L3"="L3/-", "S1"="S1/-", "S2"="S2/-", "S3"="S3/-"))
levels(devtime.Male$lines)
# Plot
stat_sum_scatter_male.dash <- ggplot(aes(x=lines, y=eclosion_date, color=allele), data=devtime.Male)+
  geom_point(size=1)+
  geom_jitter(width=0.2, height=0.05)+
  stat_summary(fun.data=data_summary, lwd=1.7, size=2)+
  geom_segment(aes(x=x1, y=y1, xend=x2, yend=y2), data=dfS, inherit.aes = FALSE, linetype=2, size=1.5)+
  geom_segment(aes(x=x1, y=y1, xend=x2, yend=y2), data=dfL, inherit.aes = FALSE, linetype=2, size=1.5)+
  ylab("Development time (days)")+
  xlab("Line")+
  theme_bw()+
  theme(axis.text = element_text(size=30),
        axis.title.x = element_text(size=45), axis.title.y = element_text(size=45), legend.position = "none",strip.text=element_text(size=30))+
  theme(axis.title.x = element_text(vjust=-0.4), axis.title.y = element_text(vjus =2.5, hjust = 0.5))+
  theme(plot.margin = unit(c(0.5,0.1,0.7,0.7), "cm"))+
  scale_colour_manual(values=c("#CA3542", "#7A8ED9"))+
  facet_wrap(~chro,scales="free_x")
stat_sum_scatter_male.dash

### female plot
summarySE(devtime.Female, measurevar=c("eclosion_date"), groupvars=c("allele_chro"))
dfS <- data.frame(chro=c("B", "D"), x1=c(3.6, 3.6), x2 = c(6.5, 6.5), y1=c(10.215, 10.43), y2=c(10.215, 10.43))
dfL <- data.frame(chro=c("B", "D"), x1=c(0.5, 0.5), x2 = c(3.4, 3.4), y1=c(10.106, 10.33), y2=c(10.106, 10.33))
## /S
levels(devtime.Female$lines)
devtime.Female$lines <- revalue(devtime.Female$lines, c("L1"="L1/S", "L2"="L2/S", "L3"="L3/S", "S1"="S1/S", "S2"="S2/S", "S3"="S3/S"))
levels(devtime.Female$lines)
#plot
stat_sum_scatter_fem.S <- ggplot(aes(x=lines, y=eclosion_date, color=allele), data=devtime.Female)+
  geom_point(size=1)+
  geom_jitter(width=0.2, height=0.05)+
  stat_summary(fun.data=data_summary, lwd=1.7, size=2)+
  geom_segment(aes(x=x1, y=y1, xend=x2, yend=y2), data=dfS, inherit.aes = FALSE, linetype=2, size=1.5)+
  geom_segment(aes(x=x1, y=y1, xend=x2, yend=y2), data=dfL, inherit.aes = FALSE, linetype=2, size=1.5)+
  ylab("Development time (days)")+
  xlab("Line")+
  theme_bw()+
  theme(axis.text = element_text(size=30),
        axis.title.x = element_text(size=45), axis.title.y = element_text(size=45), legend.position = "none",strip.text=element_text(size=30))+
  theme(axis.title.x = element_text(vjust=-0.4), axis.title.y = element_text(vjus =2.5, hjust = 0.5))+
  theme(plot.margin = unit(c(0.5,0.1,0.7,0.7), "cm"))+
  scale_colour_manual(values=c("#CA3542", "#7A8ED9"))+
  facet_wrap(~chro,scales="free_x")
stat_sum_scatter_fem.S

# reset
split.by.sex <- split(devtime, devtime$sex)
devtime.Male <- split.by.sex[[2]]
devtime.Female <- split.by.sex[[1]]

## /-
levels(devtime.Female$lines)
devtime.Female$lines <- revalue(devtime.Female$lines, c("L1"="L1/-", "L2"="L2/-", "L3"="L3/-", "S1"="S1/-", "S2"="S2/-", "S3"="S3/-"))
levels(devtime.Female$lines)
#plot
stat_sum_scatter_fem.dash <- ggplot(aes(x=lines, y=eclosion_date, color=allele), data=devtime.Female)+
  geom_point(size=1)+
  geom_jitter(width=0.2, height=0.05)+
  stat_summary(fun.data=data_summary, lwd=1.7, size=2)+
  geom_segment(aes(x=x1, y=y1, xend=x2, yend=y2), data=dfS, inherit.aes = FALSE, linetype=2, size=1.5)+
  geom_segment(aes(x=x1, y=y1, xend=x2, yend=y2), data=dfL, inherit.aes = FALSE, linetype=2, size=1.5)+
  ylab("Development time (days)")+
  xlab("Line")+
  theme_bw()+
  theme(axis.text = element_text(size=30),
        axis.title.x = element_text(size=45), axis.title.y = element_text(size=45), legend.position = "none",strip.text=element_text(size=30))+
  theme(axis.title.x = element_text(vjust=-0.4), axis.title.y = element_text(vjus =2.5, hjust = 0.5))+
  theme(plot.margin = unit(c(0.5,0.1,0.7,0.7), "cm"))+
  scale_colour_manual(values=c("#CA3542", "#7A8ED9"))+
  facet_wrap(~chro,scales="free_x")
stat_sum_scatter_fem.dash
