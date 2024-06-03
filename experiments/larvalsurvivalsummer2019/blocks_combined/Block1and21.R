### Using all data from blocks 1 and 2 of laravl survival experiments peformed by Charlotte in Summer 2019
### Testing to see if there is any difference between lines with the L or S allele and the sex ratio of surviving offspring.
### additional tests for number of survivors and relative survival rates

library('lme4')
library('ggplot2')

setwd("/Users/michaeljardine/Desktop/DTP/Data and analysis/larvalsurvivalsummer2019")
Larval.surv <- read.csv("Block1 and 2 combined.csv")
str(Larval.surv)
summary(Larval.surv)


## histograms of males, females and totals
par(mfrow=c(1,3))
hist(Larval.surv$total)
hist(Larval.surv$males)
hist(Larval.surv$females)
## males normal, fmelaes less so, totals bunched to thr right, ways to sort?

## also look at the proportions of the two sexes
par(mfrow=c(1,3))
hist(Larval.surv$propsurv)
hist(Larval.surv$propmale)
hist(Larval.surv$propfem)
## resonably normal looking apart from the totla proportion surving - fix
par(mfrow=c(1,2))
hist(log(Larval.surv$propsurv))  ## not good
hist(sqrt(Larval.surv$propsurv)) ## also not amazing - better off left alone?


#### testing for differnences in number of surviving offsrping ####

## construct models to look at the totla number of surviving offspring
## 2 models, 1 with normal error and theother with poisson - standard for count data
M.total.1 <- lmer(total ~ allele + (1|Line) + (1 | Block), data=Larval.surv)
res <- resid(M.total.1)
ran <- ranef(M.total.1)
fit <- fitted(M.total.1)
par(mfrow=c(1,3))
hist(ran$line[,1], main='')
plot(res ~ fit)
qqnorm(res, main='')
qqline(res)
summary(M.total.1)
anova(M.total.1)
Anova(M.total.1)

M.total.2 <- glmer(total ~ allele + (1|Line) + (1 | Block), data=Larval.surv, family=poisson(link=log))
res <- resid(M.total.2)
ran <- ranef(M.total.2)
fit <- fitted(M.total.2)
par(mfrow=c(1,3))
hist(ran$line[,1], main='')
plot(res ~ fit)
qqnorm(res, main='')
qqline(res)
summary(M.total.2)
anova(M.total.2)
Anova(M.total.2)

## no difference between the 2 alleles in the number of offspring surviving to adulthood.
## output of models not amazing and 1st does not converge - review later but move on for now

## tring witht he proportion surving from both sexes - basically jsut the same measurement as before as smae starting numbers
M.propsurv.1 <- lmer(propsurv ~ allele + (1|Line) + (1 | Block), data=Larval.surv)
res <- resid(M.propsurv.1)
ran <- ranef(M.propsurv.1)
fit <- fitted(M.propsurv.1)
par(mfrow=c(1,4))
hist(ran$Line[,1], main='')
hist(ran$Block[,1], main='')
plot(res ~ fit)
qqnorm(res, main='')
qqline(res)
summary(M.propsurv.1)
anova(M.propsurv.1)
Anova(M.propsurv.1)

M.propsurv.2 <- glmer(propsurv ~ allele + (1|Line) + (1 | Block), data=Larval.surv, family=binomial(link=logit))
res <- resid(M.propsurv.2)
ran <- ranef(M.propsurv.2)
fit <- fitted(M.propsurv.2)
par(mfrow=c(1,4))
hist(ran$Line[,1], main='')
hist(ran$Block[,1], main='')
plot(res ~ fit)
qqnorm(res, main='')
qqline(res)
summary(M.propsurv.2)
anova(M.propsurv.2)
Anova(M.propsurv.2)
## as expected no differnece in the proportion surving for offsrping with either allele  - as expected 
## binomial model probably no the best to use - not binary response variable -get rid of?

### relative survival rates of all offsrping ####
X.propsurv <- mean(Larval.surv$propsurv)

Larval.surv$relsurv <- Larval.surv$propsurv/X.propsurv ## mean of resurv = 1

## repeat the same models as before 
M.relsurv.1 <- lmer(relsurv ~ allele + (1 | Line) + (1 | Block), data=Larval.surv)
res <- resid(M.relsurv.1)
ran <- ranef(M.relsurv.1)
fit <- fitted(M.relsurv.1)
par(mfrow=c(1,4))
hist(ran$Line[,1], main='')
hist(ran$Block[,1], main='')
plot(res ~ fit)
qqnorm(res, main='')
qqline(res)
summary(M.relsurv.1)
anova(M.relsurv.1)
Anova(M.relsurv.1)
# no difference but fit is not great

# looking at ditribution 
par(mfrow=c(1,3))
hist(Larval.surv$relsurv)
hist(log(Larval.surv$relsurv))  
hist(sqrt(Larval.surv$relsurv))


Larval.surv$fixrelsurv <- (Larval.surv$relsurv/(1-Larval.surv$relsurv))
hist(log(Larval.surv$fixrelsurv)) 

#### sexspefic survival rates ####

## sex ratio ##
# one of the assumptions made in this experiment is that the starting sex ratio is the same
# but the adult sex ratio may alter if there are diffrrences int he realtive survival of the 2 sexes
# sexratio is recorded as the proportion of the survivng adults that are male (Larval.surv$propmale)
par(mfrow=c(1,3))
hist(Larval.surv$propmale)
hist(sqrt(Larval.surv$propmale))
hist(log(Larval.surv$propmale)) 

# since its a proportion porbably best to use glmer and binomial logit link
M.sexratio.1 <- lmer(propmale ~ allele + (1 | Line) + (1 | Block), data=Larval.surv)
res <- resid(M.sexratio.1)
ran <- ranef(M.sexratio.1)
fit <- fitted(M.sexratio.1)
par(mfrow=c(1,4))
hist(ran$Line[,1], main='')
hist(ran$Block[,1], main='')
plot(res ~ fit)
qqnorm(res, main='')
qqline(res)
summary(M.sexratio.1)
anova(M.sexratio.1)
Anova(M.sexratio.1)

M.sexratio.2 <- lmer(sqrt(propmale) ~ allele + (1 | Line) + (1 | Block), data=Larval.surv)
res <- resid(M.sexratio.2)
ran <- ranef(M.sexratio.2)
fit <- fitted(M.sexratio.2)
par(mfrow=c(1,4))
hist(ran$Line[,1], main='')
hist(ran$Block[,1], main='')
plot(res ~ fit)
qqnorm(res, main='')
qqline(res)
summary(M.sexratio.2)
anova(M.sexratio.2)
Anova(M.sexratio.2)
# no sex ratio difference due to fru allele - even poorer when sqrt

## however what we really want to know is if the version of the fru allele affects the egg to adult survival of the two sexes differently 
## for this then we need to take into account the relative survival of both sexes in each line
## this is done by dividing the proportio of males who survived by the number of females that also survived.
Larval.surv$sexrelsurv <- Larval.surv$propexpsurvmale/Larval.surv$propexpsurvfem
## this metric means that if sexrelsurv > 1 then more males survive relative to females, and if < 1 than more females survive realtive to males.

# one issue from this is row 106 which was very poor and only produced 3 males and no females. The above calucation thereofre gives an infintity (3/0)
# deleting this row will solve the problem for now

Larval.survEdit <- Larval.surv[-c(106), ]
# saved this as a new data frame LsurvEdit incase we want to use the full data set for other purposes

par(mfrow=c(1,3))
hist(Larval.survEdit$sexrelsurv)
hist(log(Larval.survEdit$sexrelsurv))
hist(sqrt(Larval.survEdit$sexrelsurv))
# logging looks much better

## now repeat models using this combined value 
# first a simple model
M.relsexsurv.1 <- lmer(sexrelsurv ~ allele + (1 | Line) + (1 | Block), data=Larval.survEdit)
res <- resid(M.relsexsurv.1)
ran <- ranef(M.relsexsurv.1)
fit <- fitted(M.relsexsurv.1)
par(mfrow=c(1,4))
hist(ran$Line[,1], main='')
hist(ran$Block[,1], main='')
plot(res ~ fit)
qqnorm(res, main='')
qqline(res)
summary(M.relsexsurv.1)
anova(M.relsexsurv.1)
Anova(M.relsexsurv.1)

# seems fine but logging looked better when looking at full data
Larval.survEdit$logsexrelsurv <- log(Larval.survEdit$sexrelsurv)


M.relsexsurv.2 <- lmer(logsexrelsurv ~ allele + (1 | Line) + (1 | Block), data=Larval.survEdit)
res <- resid(M.relsexsurv.2)
ran <- ranef(M.relsexsurv.2)
fit <- fitted(M.relsexsurv.2)
par(mfrow=c(1,4))
hist(ran$Line[,1], main='')
hist(ran$Block[,1], main='')
plot(res ~ fit)
qqnorm(res, main='')
qqline(res)
summary(M.relsexsurv.2)
anova(M.relsexsurv.2)
Anova(M.relsexsurv.2)
# dosen't change the outcome of the test of the allele but a much reduced residual standard error.
AIC(M.relsexsurv.1, M.relsexsurv.2)
# AIC is much lower in M3.2 -> logging the response variable seems juistified.

