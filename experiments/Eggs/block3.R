setwd("/Users/michaeljardine/Desktop/DTP/Datasets/Eggs")

### Packages
library(lme4)
library(Rmisc)
library(aod)
library(pbkrtest)


#### block 3 part 1 ####
eggs3 <- read.csv("blockthree101018.csv", header = TRUE)
str(eggs3)


### convert line to a factor with 6 levels
eggs3$line <- factor(eggs3$line)
str(eggs3)
### all good

### we have two numbers produced for the number eggs, predicted and corrected
x = eggs3$Predicted.count.
y = eggs3$Corrected.count.
par(mfrow=c(1,1))
plot(x, y)
### since these correlate exactly, it dosen't matter which one we use
### I will use the predicted count for consistancy

### histograms to check for normality
par(mfrow=c(2,2))
hist(eggs3$Predicted.count.)
### bunched mostly towards zero try logging?
hist(log(eggs3$Predicted.count.))
### dosen't really help, now bunched towards the right -square rooting??
hist(sqrt(eggs3$Predicted.count.))
### that actuall looks really good


par(mfrow=c(3, 1))
hist(eggs3$Corrected.count.)
hist(log(eggs3$Corrected.count.))
hist(sqrt(eggs3$Predicted.count.))
### almost exactly the same patterns with the corrected count



### summary statistics for reference, relating to the plots created in JMP
# some summary statistics for means and errors 
summarySE(eggs3, measurevar=c("Predicted.count."), groupvars=c("allele"))
summarySE(eggs3, measurevar=c("Predicted.count."), groupvars=c("line"))
summarySE(eggs3, measurevar=c("Predicted.count."), groupvars=c("gentype"))

eggs3$allelegenotype <- interaction(eggs3$allele, eggs3$gentype)
summarySE(eggs3, measurevar=c("Predicted.count."), groupvars=c("allelegenotype"))



### Models
### we want to look at the differenecs in the number of eggs laid between the lines of flies with different fruitless alleles
### flies also vary in which line they belong to (repliate of allele) and their genotype
### 2 models: 1) egg number explaiend purely by allele; 2) include allele and genotype (B/D) and thier iteraction

### these models will be linear mixed models

LMM3.1 <- lmer(Predicted.count. ~ allele + (1 | line), data=eggs3)

res <- resid(LMM3.1)
ran <- ranef(LMM3.1)
fit <- fitted(LMM3.1)
par(mfrow=c(1,3))
hist(ran$line[,1], main='')
plot(res ~ fit)
qqnorm(res, main='')
qqline(res)

summary(LMM3.1)
anova(LMM3.1)

wald.test(b=fixef(LMM3.1), Sigma=vcov(LMM3.1), Terms = 2, df=1)

### 
### 

### 

SQLMM3.1 <- lmer(sqrt(Predicted.count.) ~ allele + (1 | line), data=eggs3)

res <- resid(SQLMM3.1)
ran <- ranef(SQLMM3.1)
fit <- fitted(SQLMM3.1)
par(mfrow=c(1,3))
hist(ran$line[,1], main='')
plot(res ~ fit)
qqnorm(res, main='')
qqline(res)

summary(SQLMM3.1)
anova(SQLMM3.1)

wald.test(b=fixef(SQLMM3.1), Sigma=vcov(SQLMM3.1), Terms = 2, df=1)


### looks perhaps a litle better compare
AIC(LMM3.1, SQLMM3.1)
anova(LMM3.1, SQLMM3.1)
## sqrt justified


### we want to account for the variation caused by the fly's genotpe as well

LMM3.2 <- lmer(Predicted.count. ~ allele*gentype + (1 | line), data=eggs3)

res <- resid(LMM3.2)
ran <- ranef(LMM3.2)
fit <- fitted(LMM3.2)
par(mfrow=c(1,3))
hist(ran$line[,1], main='')
plot(res ~ fit)
qqnorm(res, main='')
qqline(res)

summary(LMM3.2)
anova(LMM3.2)

wald.test(b=fixef(LMM3.2), Sigma=vcov(LMM3.2), Terms = 2, df=1)
wald.test(b=fixef(LMM3.2), Sigma=vcov(LMM3.2), Terms = 3, df=1)
wald.test(b=fixef(LMM3.2), Sigma=vcov(LMM3.2), Terms = 4, df=1)

### again do decernable affect of allele
### high and modertae F-values for genotype and an interaction
### problems maintained with the different metrics, what is correct??

### lets see what sqrt does this time

SQLMM3.2 <- lmer(sqrt(Predicted.count.) ~ allele*gentype + (1 | line), data=eggs3)

res <- resid(SQLMM3.2)
ran <- ranef(SQLMM3.2)
fit <- fitted(SQLMM3.2)
par(mfrow=c(1,3))
hist(ran$line[,1], main='')
plot(res ~ fit)
qqnorm(res, main='')
qqline(res)

summary(SQLMM3.2)
anova(SQLMM3.2)

wald.test(b=fixef(SQLMM3.2), Sigma=vcov(SQLMM3.2), Terms = 2, df=1)
wald.test(b=fixef(SQLMM3.2), Sigma=vcov(SQLMM3.2), Terms = 3, df=1)
wald.test(b=fixef(SQLMM3.2), Sigma=vcov(SQLMM3.2), Terms = 4, df=1)

### COMPARE THE TWO MODELS
AIC(LMM3.2, SQLMM3.2)
anova(LMM3.2, SQLMM3.2)

### again the sqrt appears justified

### can we compare the 4 models??
AIC(LMM3.1, SQLMM3.1, LMM3.2, SQLMM3.2)
anova(LMM3.1, SQLMM3.1, LMM3.2, SQLMM3.2)

### not sure what this tells us but mabe that includeing even more terms is justified


#### block 3 part 2 ####


### although this is still block three we'll all in 4 just now until we work out how to combine it with the rest of block 3
### first we'll do the analysis as though it is anotht block and if the prgram has worked well

eggs4 <- read.csv("blockthree121018.csv", header = TRUE)
str(eggs4)
summary(eggs4)

### convert line to a factor with 6 levels
eggs4$line <- factor(eggs4$line)
str(eggs4)
### all good

### we have two numbers produced for the number eggs, predicted and corrected
x = eggs4$Predicted.count.
y = eggs4$Corrected.count.
par(mfrow=c(1,1))
plot(x, y)
### since these correlate exactly, it dosen't matter which one we use
### I will use the predicted count for consistancy

### histograms to check for normality
par(mfrow=c(2,2))
hist(eggs4$Predicted.count.)
### bunched mostly towards zero try logging?
hist(log(eggs4$Predicted.count.))
### dosen't really help, now bunched towards the right -square rooting??
hist(sqrt(eggs4$Predicted.count.))
### that actuall looks really good


par(mfrow=c(3, 1))
hist(eggs4$Corrected.count.)
hist(log(eggs4$Corrected.count.))
hist(sqrt(eggs4$Predicted.count.))
### almost exactly the same patterns with the corrected count



### summary statistics for reference, relating to the plots created in JMP
# some summary statistics for means and errors 
summarySE(eggs4, measurevar=c("Predicted.count."), groupvars=c("allele"))
summarySE(eggs4, measurevar=c("Predicted.count."), groupvars=c("line"))
summarySE(eggs4, measurevar=c("Predicted.count."), groupvars=c("genotype"))

eggs4$allelegenotype <- interaction(eggs4$allele, eggs4$genotype)
summarySE(eggs4, measurevar=c("Predicted.count."), groupvars=c("allelegenotype"))



### Models
### we want to look at the differenecs in the number of eggs laid between the lines of flies with different fruitless alleles
### flies also vary in which line they belong to (repliate of allele) and their genotype
### 2 models: 1) egg number explaiend purely by allele; 2) include allele and genotype (B/D) and thier iteraction

### these models will be linear mixed models

LMM4.1 <- lmer(Predicted.count. ~ allele + (1 | line), data=eggs4)

res <- resid(LMM4.1)
ran <- ranef(LMM4.1)
fit <- fitted(LMM4.1)
par(mfrow=c(1,3))
hist(ran$line[,1], main='')
plot(res ~ fit)
qqnorm(res, main='')
qqline(res)

summary(LMM4.1)
anova(LMM4.1)

wald.test(b=fixef(LMM4.1), Sigma=vcov(LMM4.1), Terms = 2, df=1)

### 
### 

### 

SQLMM4.1 <- lmer(sqrt(Predicted.count.) ~ allele + (1 | line), data=eggs4)

res <- resid(SQLMM4.1)
ran <- ranef(SQLMM4.1)
fit <- fitted(SQLMM4.1)
par(mfrow=c(1,3))
hist(ran$line[,1], main='')
plot(res ~ fit)
qqnorm(res, main='')
qqline(res)

summary(SQLMM4.1)
anova(SQLMM4.1)

wald.test(b=fixef(SQLMM4.1), Sigma=vcov(SQLMM4.1), Terms = 2, df=1)


### looks perhaps a litle better compare
AIC(LMM4.1, SQLMM4.1)
anova(LMM4.1, SQLMM4.1)
## sqrt justified


### we want to account for the variation caused by the fly's genotpe as well

LMM4.2 <- lmer(Predicted.count. ~ allele*genotype + (1 | line), data=eggs4)

res <- resid(LMM4.2)
ran <- ranef(LMM4.2)
fit <- fitted(LMM4.2)
par(mfrow=c(1,3))
hist(ran$line[,1], main='')
plot(res ~ fit)
qqnorm(res, main='')
qqline(res)

summary(LMM4.2)
anova(LMM4.2)

wald.test(b=fixef(LMM4.2), Sigma=vcov(LMM4.2), Terms = 2, df=1)
wald.test(b=fixef(LMM4.2), Sigma=vcov(LMM4.2), Terms = 3, df=1)
wald.test(b=fixef(LMM4.2), Sigma=vcov(LMM4.2), Terms = 4, df=1)

### 
### 
### 

### lets see what sqrt does this time

SQLMM4.2 <- lmer(sqrt(Predicted.count.) ~ allele*genotype + (1 | line), data=eggs4)

res <- resid(SQLMM4.2)
ran <- ranef(SQLMM4.2)
fit <- fitted(SQLMM4.2)
par(mfrow=c(1,3))
hist(ran$line[,1], main='')
plot(res ~ fit)
qqnorm(res, main='')
qqline(res)

summary(SQLMM4.2)
anova(SQLMM4.2)

wald.test(b=fixef(SQLMM4.2), Sigma=vcov(SQLMM4.2), Terms = 2, df=1)
wald.test(b=fixef(SQLMM4.2), Sigma=vcov(SQLMM4.2), Terms = 3, df=1)
wald.test(b=fixef(SQLMM4.2), Sigma=vcov(SQLMM4.2), Terms = 4, df=1)

### COMPARE THE TWO MODELS
AIC(LMM4.2, SQLMM4.2)
anova(LMM4.2, SQLMM4.2)

### again the sqrt appears justified

### can we compare the 4 models??
AIC(LMM4.1, SQLMM4.1, LMM4.2, SQLMM4.2)
anova(LMM4.1, SQLMM4.1, LMM4.2, SQLMM4.2)

### not sure what this tells us but maybe that including even more terms is justified
