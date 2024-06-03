setwd("/Users/michaeljardine/Desktop/DTP/Data and analysis/larvalsurvivalsummer2019")
Lsurv <- read.csv("Larvalsurvivalblock1.csv")
str(Lsurv)
summary(Lsurv)


### histograms of data

#numbers of surving mlaes, females and total
par(mfrow=c(1,3))
hist(Lsurv$males)
hist(Lsurv$females)
hist(Lsurv$total)

#proportions of males and females
par(mfrow=c(1,2))
hist(Lsurv$propmale)
hist(Lsurv$propfem)

# proportions of survival
par(mfrow=c(1,3))
hist(Lsurv$propsurv)
hist(Lsurv$propexpsurvmale)
hist(Lsurv$propexpsurvfem)

## all historgrams look reasonably normal

## some simple models 


#### effect of fru on proportion of surviving larave ####
# does the proportion of surviving offspring vary between lines with differnet fru alleles  
M1.1 <- glmer(propsurv ~ allele +  (1 | Line), data=Lsurv, family=binomial(link=logit))
summary(M1.1)
anova(M1.1)
Anova(M1.1)
# warning messages appear
# no realtionship between fru allele and the proportion of eggs reaching adulthood

# really we want to compare the realtive survival of the lines - do this by dividing the proportion survining by the average survival rate across all the vials 
# mean propsurv = 0.4627
Lsurv$relsurv <- Lsurv$propsurv/0.4627

## what happens if we run a simialr model again?
M2.1 <- lmer(Lsurv$relsurv ~ allele +  (1 | Line), data=Lsurv)
summary(M2.1)
anova(M2.1)
Anova(M2.1)
 # there appears to be no difference in the relative survival due to the fru allele 



#### sex specific larval survival ####
## however what we really want to know is if the version of the fru allele affects the egg to adult survival of the two sexes differently 
## for this then we need to take into account the relative survival of both sexes in each line
## this is done by dividing the proportio of males who survived by the number of females that also survived.
Lsurv$sexrelsurv <- Lsurv$propexpsurvmale/Lsurv$propexpsurvfem
## this metric means that if sexrelsurv > 1 then more males survive relative to females, and if < 1 than more females survive realtive to males.

# one issue from this is row 106 which was very poor and only produced 3 males and no females. The above calucation thereofre gives an infintity (3/0)
# deleting this row will solve the problem for now

LsurvEdit <- Lsurv[-c(106), ]
# saved this as a new data frame LsurvEdit incase we want to use the full data set for other purposes

par(mfrow=c(1,2))
hist(LsurvEdit$sexrelsurv)
hist(log(LsurvEdit$sexrelsurv))
# logging looks much better


## now repeat models using this combined value 
# first a simple model
M3.1 <- lm(sexrelsurv ~ allele, data=LsurvEdit)
par(mfrow= c(2,2))
plot(M3.1)

summary(M3.1)
anova(M3.1)

# seems fine but logging looked better when looking at full data
M3.2 <- lm(log(sexrelsurv) ~ allele, data=LsurvEdit)
par(mfrow= c(2,2))
plot(M3.2)

summary(M3.2)
anova(M3.2)
# dosen't change the outcome of the test of the allele but a much reduced residual standard error.
AIC(M3.1, M3.2)
# AIC is much lower in M3.2 -> logging the response variable seems juistified.


## but we have multiple lines which must be takedn into account and included as a random factor
## therefore we build a liner mixed model 
M3.3 <- lmer(log(sexrelsurv) ~ allele + (1 | Line), data=LsurvEdit)

res <- resid(M3.3)
ran <- ranef(M3.3)
fit <- fitted(M3.3)
par(mfrow=c(1,3))
hist(ran$line[,1], main='')
plot(res ~ fit)
qqnorm(res, main='')
qqline(res)

summary(M3.3)
Anova(M3.3)
# again no affect of the fru allele of the relative survival of males and females.
# is this the best method of analysing this? 
# maybe need to and the simualtions as Max did before? Or wald tests? 

## we can now look at other possble ideas for differences between the lines: e.g. sex ratio. Also we would like some overall stats to compare e.g overall numbers of survival.

#### adult sex ratio ####

# one of the assumptions made in this experiment is that the starting sex ratio is the same
# while we can do nothing about that and from this it appears the reltive survival of the two sexes is the smae regardless of genotype probs best to check sex ratio its self

# sexratio is recorded as the proportion of the survivng adults that are male (Lsurv$propmale)

hist(Lsurv$propmale)
hist(sqrt(Lsurv$propmale))

# since its a proportion porbably best to use glmer and binomial logit link
M4.1 <- glmer(propmale ~ allele + (1 | Line), data=LsurvEdit, family=binomial(link=logit))

res <- resid(M4.1)
ran <- ranef(M4.1)
fit <- fitted(M4.1)
par(mfrow=c(1,3))
hist(ran$line[,1], main='')
plot(res ~ fit)
qqnorm(res, main='')
qqline(res)

summary(M4.1)
Anova(M4.1)
# dosen't seem to any affect of the fru allele on sex ratio

# check models as problems of convergenece
M4.2 <- lmer(propmale ~ allele + (1 | Line), data=LsurvEdit)

res <- resid(M4.2)
ran <- ranef(M4.2)
fit <- fitted(M4.2)
par(mfrow=c(1,3))
hist(ran$line[,1], main='')
plot(res ~ fit)
qqnorm(res, main='')
qqline(res)

summary(M4.2)
Anova(M4.2)

AIC(M4.1, M4.2)
anova(M4.1, M4.2)
# so the second model seems better? no affect of allele on sex ratio anyway so not a huge issue?

#### number of surviving offspring ####

# we can do this in three ways: no. of males, no. of females, and total number.

# do a simple model of each of these three.
# count data so check if log values may be better? - do bot and compare

# total

M5.1 <- lmer(total ~ allele + (1 | Line), data=Lsurv)

res <- resid(M5.1)
ran <- ranef(M5.1)
fit <- fitted(M5.1)
par(mfrow=c(1,3))
hist(ran$line[,1], main='')
plot(res ~ fit)
qqnorm(res, main='')
qqline(res)

summary(M5.1)
anova(M5.1)
Anova(M5.1)

# try the logging of the data?
M5.2 <- lmer(log(total) ~ allele + (1 | Line), data=LsurvEdit)

res <- resid(M5.2)
ran <- ranef(M5.2)
fit <- fitted(M5.2)
par(mfrow=c(1,3))
hist(ran$line[,1], main='')
plot(res ~ fit)
qqnorm(res, main='')
qqline(res)

summary(M5.2)
anova(M5.2)
Anova(M5.2)

AIC(M5.1, M5.2)
# maybe logging looks better?
# no affect of fru allele

### males 
M6.1 <- lmer(males ~ allele + (1 | Line), data=Lsurv)

res <- resid(M6.1)
ran <- ranef(M6.1)
fit <- fitted(M6.1)
par(mfrow=c(1,3))
hist(ran$line[,1], main='')
plot(res ~ fit)
qqnorm(res, main='')
qqline(res)

summary(M6.1)
anova(M6.1)
Anova(M6.1)

# try the logging of the data?
M6.2 <- lmer(log(males) ~ allele + (1 | Line), data=Lsurv)

res <- resid(M6.2)
ran <- ranef(M6.2)
fit <- fitted(M6.2)
par(mfrow=c(1,3))
hist(ran$line[,1], main='')
plot(res ~ fit)
qqnorm(res, main='')
qqline(res)

summary(M6.2)
anova(M6.2)
Anova(M6.2)

AIC(M6.1, M6.2)
## logging improves the fitting of the model - logginf seems to wok for smaller counts data

### females
# remove row 106 due to zero femlaes recorded and log(0) = Inf
M7.1 <- lmer(females ~ allele + (1 | Line), data=LsurvEdit)

res <- resid(M7.1)
ran <- ranef(M7.1)
fit <- fitted(M7.1)
par(mfrow=c(1,3))
hist(ran$line[,1], main='')
plot(res ~ fit)
qqnorm(res, main='')
qqline(res)

summary(M7.1)
anova(M7.1)
Anova(M7.1)

# try the logging of the data?
M7.2 <- lmer(log(females) ~ allele + (1 | Line), data=LsurvEdit)

res <- resid(M7.2)
ran <- ranef(M7.2)
fit <- fitted(M7.2)
par(mfrow=c(1,3))
hist(ran$line[,1], main='')
plot(res ~ fit)
qqnorm(res, main='')
qqline(res)

summary(M7.2)
anova(M7.2)
Anova(M7.2)

AIC(M7.1, M7.2)
# again logging seems justified.
# no differences between larval survival due to the fru allele 
