## The fru alllee that a fly has does not affect its survival from egg to adult hood when homozygous
## but it does when it is homozygous, or with a balancer chromosome 
## so I take the relac=vent data from the chro compiment experiment and include it with the previous 2 blocks as a 3rd block

setwd("/Users/michaeljardine/Desktop/DTP/Data and analysis/larvalsurvivalsummer2019/3_blocks_comb")
Larval.surv.B3 <- read.csv("Block1and2 plus chro comp data.csv")
str(Larval.surv.B3)
summary(Larval.surv.B3)

## packages
library('lme4')
library('ggplot2')
library('Rmisc')
library('car')


#### plots and summary stats ####

## histograms of males, females and totals
par(mfrow=c(1,3))
hist(Larval.surv.B3$total)
hist(Larval.surv.B3$males)
hist(Larval.surv.B3$females)
## all bunched to the left, less so for the totals (higher means)

# try logging?
hist(log(Larval.surv.B3$total))
hist(log(Larval.surv.B3$males))
hist(log(Larval.surv.B3$females))
# maybe better for the two sexes but not great for totals

# or sqrt?
hist(sqrt(Larval.surv.B3$total))
hist(sqrt(Larval.surv.B3$males))
hist(sqrt(Larval.surv.B3$females))
## they actually look much better - irratting that ther is no consistancy between all the combination sof the blocks in the normailsation of the data
## when it comes to testing the number of offspring will include models of all 3 and comapre

## also look at the proportions of the two sexes
par(mfrow=c(1,3))
hist(Larval.surv.B3$propsurv)
hist(Larval.surv.B3$propmale)
hist(Larval.surv.B3$propfem)
## resonably normal looking apart from the total proportion surving - fix
par(mfrow=c(1,2))
hist(log(Larval.surv.B3$propsurv))  ## not good
hist(sqrt(Larval.surv.B3$propsurv)) ## also not amazing but ok - test later

## calculating means to help interperate model results later

# means for total surviving offspring
summarySE(Larval.surv.B3, measurevar=c("total"), groupvars=c("allele"))
# L = 22.032, S = 16.66

summarySE(Larval.surv.B3, measurevar=c("total"), groupvars=c("Line"))
# L1  = 20.56, L2 = 21.013, L3  = 24.35, S1 = 17.62, S2 = 22.32, S3 = 9.11111

summarySE(Larval.surv.B3, measurevar=c("total"), groupvars=c("Block"))
# B1 = 23.133, B2 = 28.88, B3 = 11.4
# why was B3 so poor? - this was the line that was outcrossed? - something wrong with the food during this block?

# means for proportion surviving
summarySE(Larval.surv.B3, measurevar=c("propsurv"), groupvars=c("allele"))
# L = 0.44065, S = 0.3332

summarySE(Larval.surv.B3, measurevar=c("propsurv"), groupvars=c("Line"))
# L1  = 0.41125, L2 = 0.4203, L3  = 0.487, S1 = 0.352, S2 = 0.446, S3 = 0.1822

summarySE(Larval.surv.B3, measurevar=c("propsurv"), groupvars=c("Block"))
# B1 = 0.463, B2 = 0.578, B3 = 0.228

# means for proportion male - sexratio
summarySE(Larval.surv.B3, measurevar=c("propmale"), groupvars=c("allele"))
# L = 0.4855, S = 0.4966

summarySE(Larval.surv.B3, measurevar=c("propmale"), groupvars=c("Line"))
# L1  = 0.473, L2 = 0.4905, L3  = 0.4916, S1 = 0.4899, S2 = 0.4908, S3 = 0.5131

summarySE(Larval.surv.B3, measurevar=c("propmale"), groupvars=c("Block"))
# B1 = 0.5148, B2 = 0.4745, B3 = 0.48


#################################################################################################
#### testing for differences in number of surviving offsrping ####

## construct models to look at the totla number of surviving offspring
## 2 models, 1 with normal error and theother with poisson - standard for count data
M.totalB3.1 <- lmer(total ~ allele + (1|allele:Line) + (1 | Block), data=Larval.surv.B3)
res <- resid(M.totalB3.1)
ran <- ranef(M.totalB3.1)
fit <- fitted(M.totalB3.1)
par(mfrow=c(1,4))
hist(ran$`allele:Line`[,1], main='')
hist(ran$Block[,1], main="")
plot(res ~ fit)
qqnorm(res, main='')
qqline(res)
summary(M.totalB3.1)
anova(M.totalB3.1)
Anova(M.totalB3.1)
## no affect of allele - plots look quite good

M.totalB3.2 <- lmer(log(total) ~ allele + (1|allele:Line) + (1 | Block), data=Larval.surv.B3)
res <- resid(M.totalB3.2)
ran <- ranef(M.totalB3.2)
fit <- fitted(M.totalB3.2)
par(mfrow=c(1,4))
hist(ran$`allele:Line`[,1], main='')
hist(ran$Block[,1], main="")
plot(res ~ fit)
qqnorm(res, main='')
qqline(res)
summary(M.totalB3.2)
anova(M.totalB3.2)
Anova(M.totalB3.2)
# allele closer to 0.05 but plots actaully look worse

M.totalB3.3 <- lmer(sqrt(total) ~ allele + (1|allele:Line) + (1 | Block), data=Larval.surv.B3)
res <- resid(M.totalB3.3)
ran <- ranef(M.totalB3.3)
fit <- fitted(M.totalB3.3)
par(mfrow=c(1,4))
hist(ran$`allele:Line`[,1], main='')
hist(ran$Block[,1], main="")
plot(res ~ fit)
qqnorm(res, main='')
qqline(res)
summary(M.totalB3.3)
anova(M.totalB3.3)
Anova(M.totalB3.3)
# allele still not involved - plots lokk better than logging

## compare models
anova(M.totalB3.1, M.totalB3.2, M.totalB3.3)

## the second model has much lower AIc and much higehr logliklihood, despite looking the worst from the plots
## best to go with this model as the numbers are so much better.

#### repeat the above but with proportion surviving ####




