setwd("/Users/michaeljardine/Desktop/DTP/Datasets/Eggs")

### Packages
library(lme4)
library(Rmisc)
library(aod)

eggs <- read.csv("blockone300518.csv")
str(eggs)


### convert line to a factor with 6 levels
eggs$line <- factor(eggs$line)
str(eggs)
### all good

### we have two numbers produced for the number eggs, predicted and corrected
x = eggs$Predicted.count.
y = eggs$Corrected.count.
cor(x, y)
### since these correlate exactly, it dosen't matter which one we use
### I will use the predicted count for consistancy

### histograms to check for normality
hist(eggs$Predicted.count.)
### looks fairly normal, a slight righthand tail, try logging?
hist(log(eggs$Predicted.count.))
### definately dosen't help, continue with the plain data

hist(eggs$Corrected.count.)
hist(log(eggs$Corrected.count.))

### summary statistics for reference, relating to the plots created in JMP
# some summary statistics for means and errors 
summarySE(eggs, measurevar=c("Predicted.count."), groupvars=c("allele"))
summarySE(eggs, measurevar=c("Predicted.count."), groupvars=c("line"))
summarySE(eggs, measurevar=c("Predicted.count."), groupvars=c("genotype"))
summarySE(eggs, measurevar=c("Predicted.count."), groupvars=c("age"))


### Models
### we want to look att eh differenecs in the number of eggs liad between the different lines, genotypes and alleles, while also controlling for variatrion due to age.

#### Full model ####

# we'll first do a model including all factors and an interaction between line and genotype
# this should work within a simple linear model

full <- lm(Predicted.count. ~ allele + line*genotype + age, data=eggs)

par(mfrow=c(2,2))
plot(full)
### looks reasonably good, slight devaition within the Q-Q plot at the right-hand end
summary(full)
anova(full)

### so it appears that there is a variation in egg number due to line.
### however in the breakdown of the lines this may only be due to the very low numbers from M16
### no differences due to other factors, although allele is close


#### remove genotype from the model ####

### since genotype expplained the least in the full model I'll run the modle again from the next one
### this also removes the interaction between line and genotype but this explined the seconf least amount of variation.

minusgenotype <- lm(Predicted.count. ~ allele + line + age, data=eggs)

par(mfrow=c(2,2))
plot(minusgenotype)
### looks reasonably good, slight devaition within the Q-Q plot at the right-hand end
summary(minusgenotype)
anova(minusgenotype)


### again only a difference for M16 due to it laying many fewer eggs
### in this and the full model no values have been calculated for M95. Why? Lack of adaquate replication??
### the amount of variatce explained by the version of the introgressed allele is simialr in both models

### neither model is great overall, describing only ~9% of the data

### Compare models, justfy dropping 
AIC(full, minusgenotype)

### the minus genotype model has a lower AIC score (not by much) largly due to dropping 6 extra degrees of freedom

#### But really line is not of specific interest as the two triplets of lines are replicates for the fruitless allele intrgressed ####
### nowa new model will be made that will use line as a random factor to look at variation due to the two alleles
### this will also include genotype, both as a factors and also an inclusion in as econd model as part of the random variaance component


### 1st mixed model - random intercept for line

ranline <- lmer(Predicted.count. ~ allele*genotype + (1 | line), data=eggs)

res <- resid(ranline)
ran <- ranef(ranline)
fit <- fitted(ranline)
par(mfrow=c(1,3))
hist(ran$line[,1], main='')
plot(res ~ fit)
qqnorm(res, main='')
qqline(res)

summary(ranline)
anova(ranline)
wald.test(b=fixef(ranline), Sigma=vcov(ranline), Terms = 2, df=1)
wald.test(b=fixef(ranline), Sigma=vcov(ranline), Terms = 3, df=1)
wald.test(b=fixef(ranline), Sigma=vcov(ranline), Terms = 4, df=1)

### and now model with random genotype compoent as well

ranlinegeno <- lmer(Predicted.count. ~ allele*genotype + (1 + genotype | line), data=eggs)

res <- resid(ranlinegeno)
ran <- ranef(ranlinegeno)
fit <- fitted(ranlinegeno)
par(mfrow=c(1,3))
hist(ran$line[,1], main='')
plot(res ~ fit)
qqnorm(res, main='')
qqline(res)

summary(ranlinegeno)
anova(ranlinegeno)
wald.test(b=fixef(ranlinegeno), Sigma=vcov(ranlinegeno), Terms = 2, df=1)
wald.test(b=fixef(ranlinegeno), Sigma=vcov(ranlinegeno), Terms = 3, df=1)
wald.test(b=fixef(ranlinegeno), Sigma=vcov(ranlinegeno), Terms = 4, df=1)

### again no effect of the allele on fecundity



#### Just allele ####

justallele <- lm(Predicted.count. ~ allele, data=eggs)

plot(justallele)

summary(justallele)
anova(justallele)

#### age only ####

ageonly <- lm(Predicted.count. ~ age, data=eggs)

plot(ageonly)

summary(ageonly)
anova(ageonly)

#### how about just the two together?? ####

ageallele <- lm(Predicted.count. ~ allele + age, data=eggs)

plot(ageallele)

summary(ageallele)
anova(ageallele)

### no combination gives there being an effect of allele on the number of eggs laid or and effect of age.
### however we will keep age in the main models as it is good to control for this
### also in the main two models allele is almost sigificant (p~0.09) so maybe with a larger sample size (especially for M95) we will see an effect of allele too




#### block 2 ####

setwd("/Users/michaeljardine/Desktop/DTP/Datasets/Eggs")

### Packages
library(lme4)
library(Rmisc)
library(aod)
library(pbkrtest)

eggs2 <- read.csv("blocktwo060818.csv")
str(eggs2)


### convert line to a factor with 6 levels
eggs2$line <- factor(eggs2$line)
str(eggs2)
### all good

### we have two numbers produced for the number eggs, predicted and corrected
x = eggs2$Predicted.count.
y = eggs2$Corrected.count.
cor(x, y)
### since these correlate exactly, it dosen't matter which one we use
### I will use the predicted count for consistancy

### histograms to check for normality
par(mfrow=c(2,2))
hist(eggs2$Predicted.count.)
### bunched mostly towards zero try logging?
hist(log(eggs2$Predicted.count.))
### dosen't really help, now bunched towards the right -square rooting??
hist(sqrt(eggs2$Predicted.count.))
### that actuall looks really good


par(mfrow=c(3, 1))
hist(eggs$Corrected.count.)
hist(log(eggs$Corrected.count.))
hist(sqrt(eggs2$Predicted.count.))
### almost exactly the same patterns with the corrected count



### summary statistics for reference, relating to the plots created in JMP
# some summary statistics for means and errors 
summarySE(eggs2, measurevar=c("Predicted.count."), groupvars=c("allele"))
summarySE(eggs2, measurevar=c("Predicted.count."), groupvars=c("line"))
summarySE(eggs2, measurevar=c("Predicted.count."), groupvars=c("genotype"))

### Models
### we want to look at the differenecs in the number of eggs laid between the lines of flies with different fruitless alleles
### flies also vary in which line they belong to (repliate of allele) and their genotype
### 2 models: 1) egg number explaiend purely by allele; 2) include allele and genotype (B/D) and thier iteraction

### these models will be linear mixed models

LMM1 <- lmer(Predicted.count. ~ allele + (1 | line), data=eggs2)

res <- resid(LMM1)
ran <- ranef(LMM1)
fit <- fitted(LMM1)
par(mfrow=c(1,3))
hist(ran$line[,1], main='')
plot(res ~ fit)
qqnorm(res, main='')
qqline(res)

summary(LMM1)
anova(LMM1)

wald.test(b=fixef(LMM1), Sigma=vcov(LMM1), Terms = 2, df=1)

### no affect of allele on the number of eggs laid
### would have been opposite direction to that expected anyway

### sqrt looked good on hist before and diagnostic plots could be better - worth trying?

SQLMM1 <- lmer(sqrt(Predicted.count.) ~ allele + (1 | line), data=eggs2)

res <- resid(SQLMM1)
ran <- ranef(SQLMM1)
fit <- fitted(SQLMM1)
par(mfrow=c(1,3))
hist(ran$line[,1], main='')
plot(res ~ fit)
qqnorm(res, main='')
qqline(res)

summary(SQLMM1)
anova(SQLMM1)

wald.test(b=fixef(SQLMM1), Sigma=vcov(SQLMM1), Terms = 2, df=1)


### looks perhaps a litle better compare
AIC(LMM1, SQLMM1)
anova(LMM1, SQLMM1)
## sqrt justified


### we now want to account for the variation caused by the fly's genotpe as well

LMM2 <- lmer(Predicted.count. ~ allele*genotype + (1 | line), data=eggs2)

res <- resid(LMM2)
ran <- ranef(LMM2)
fit <- fitted(LMM2)
par(mfrow=c(1,3))
hist(ran$line[,1], main='')
plot(res ~ fit)
qqnorm(res, main='')
qqline(res)

summary(LMM2)
anova(LMM2)

wald.test(b=fixef(LMM2), Sigma=vcov(LMM2), Terms = 2, df=1)
wald.test(b=fixef(LMM2), Sigma=vcov(LMM2), Terms = 3, df=1)
wald.test(b=fixef(LMM2), Sigma=vcov(LMM2), Terms = 4, df=1)

### again do decernable affect of allele
### high and modertae F-values for genotype and an interaction
### problems maintained with the different metrics, what is correct??

### lets see what sqrt does this time

SQLMM2 <- lmer(sqrt(Predicted.count.) ~ allele*genotype + (1 | line), data=eggs2)

res <- resid(SQLMM2)
ran <- ranef(SQLMM2)
fit <- fitted(SQLMM2)
par(mfrow=c(1,3))
hist(ran$line[,1], main='')
plot(res ~ fit)
qqnorm(res, main='')
qqline(res)

summary(SQLMM2)
anova(SQLMM2)

wald.test(b=fixef(SQLMM2), Sigma=vcov(SQLMM2), Terms = 2, df=1)
wald.test(b=fixef(SQLMM2), Sigma=vcov(SQLMM2), Terms = 3, df=1)
wald.test(b=fixef(SQLMM2), Sigma=vcov(SQLMM2), Terms = 4, df=1)

### COMPARE THE TWO MODELS
AIC(LMM2, SQLMM2)
anova(LMM2, SQLMM2)

### again the sqrt appears justified

### can we compare the 4 models??
AIC(LMM1, SQLMM1, LMM2, SQLMM2)
anova(LMM1, SQLMM1, LMM2, SQLMM2)

### not sure what this tells us but mabe that includeing even more terms is justified


#### block 3 part 1 ####

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
summary(eggs4)
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
