setwd("/Users/michaeljardine/Desktop/DTP/Datasets/Eggs")

### Packages
library(lme4)
library(Rmisc)
library(aod)

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


#######################################
##### Blocks 1 and 2 - all data ######
  
### this will follow the same basic code as that for the second block
### the data set will be termed 'Teggs' for total eggs
### we first need to correctt he parts of the dta set for factros and plot the data for patterns as before

Teggs <- read.csv("alleggblocks.csv")
str(Teggs)


### convert line to a factor with 6 levels
Teggs$line <- factor(Teggs$line)
str(Teggs)
Teggs$block <- factor(Teggs$block)
str(Teggs)
### all good

### we have two numbers produced for the number eggs, predicted and corrected
x = eggs2$Predicted.count.
y = eggs2$Corrected.count.
cor(x, y)
### since these correlate exactly, it dosen't matter which one we use
### I will use the predicted count for consistancy

### histograms to check for normality
par(mfrow=c(2,2))
hist(Teggs$Predicted.count.)
### bunched mostly towards zero try logging?
hist(log(Teggs$Predicted.count.))
### dosen't really help, now bunched towards the right -square rooting??
hist(sqrt(Teggs$Predicted.count.))
### that actuall looks really good


par(mfrow=c(3, 1))
hist(Teggs$Corrected.count.)
hist(log(Teggs$Corrected.count.))
hist(sqrt(Teggs$Predicted.count.))
### almost exactly the same patterns with the corrected count


# some summary statistics for means and errors 
summarySE(Teggs, measurevar=c("Predicted.count."), groupvars=c("allele"))
summarySE(Teggs, measurevar=c("Predicted.count."), groupvars=c("line"))
summarySE(Teggs, measurevar=c("Predicted.count."), groupvars=c("genotype"))
summarySE(Teggs, measurevar=c("Predicted.count."), groupvars=c("block"))


Teggs$allelegenotype <- interaction(Teggs$allele, Teggs$genotype)
summarySE(Teggs, measurevar=c("Predicted.count."), groupvars=c("allelegenotype"))

### Models
### we want to look at the differenecs in the number of eggs laid between the lines of flies with different fruitless alleles
### flies also vary in which line they belong to (repliate of allele) and their genotype
### we also have 2 block which may var independantly so include as a random factor in all models
### 2 models: 1) egg number explaiend purely by allele; 2) include allele and genotype (B/D) and thier iteraction

### these models will be linear mixed models

TLMM1 <- lmer(Predicted.count. ~ allele + (1 | line) + (1 | block), data=Teggs)

res <- resid(TLMM1)
ran <- ranef(TLMM1)
fit <- fitted(TLMM1)
par(mfrow=c(1,3))
hist(ran$line[,1], main='')
plot(res ~ fit)
qqnorm(res, main='')
qqline(res)

summary(TLMM1)
anova(TLMM1)

wald.test(b=fixef(TLMM1), Sigma=vcov(TLMM1), Terms = 2, df=1)

### no affect of allele on the number of eggs laid
### would have been opposite direction to that expected anyway

### sqrt looked good on hist before and diagnostic plots could be better - worth trying?

TSQLMM1 <- lmer(sqrt(Predicted.count.) ~ allele + (1 | line) + (1 | block), data=Teggs)

res <- resid(TSQLMM1)
ran <- ranef(TSQLMM1)
fit <- fitted(TSQLMM1)
par(mfrow=c(1,3))
hist(ran$line[,1], main='')
plot(res ~ fit)
qqnorm(res, main='')
qqline(res)

summary(TSQLMM1)
anova(TSQLMM1)

wald.test(b=fixef(TSQLMM1), Sigma=vcov(TSQLMM1), Terms = 2, df=1)


### looks perhaps a litle better compare
AIC(TLMM1, TSQLMM1)
anova(TLMM1, TSQLMM1)
## sqrt justified


### we now want to account for the variation caused by the fly's genotpe as well

TLMM2 <- lmer(Predicted.count. ~ allele*genotype + (1 | line) + (1 | block), data=Teggs)

res <- resid(TLMM2)
ran <- ranef(TLMM2)
fit <- fitted(TLMM2)
par(mfrow=c(1,3))
hist(ran$line[,1], main='')
plot(res ~ fit)
qqnorm(res, main='')
qqline(res)

summary(TLMM2)
anova(TLMM2)

wald.test(b=fixef(TLMM2), Sigma=vcov(TLMM2), Terms = 2, df=1)
wald.test(b=fixef(TLMM2), Sigma=vcov(TLMM2), Terms = 3, df=1)
wald.test(b=fixef(TLMM2), Sigma=vcov(TLMM2), Terms = 4, df=1)

### again do decernable affect of allele
### high and modertae F-values for genotype and an interaction
### problems maintained with the different metrics, what is correct??

### lets see what sqrt does this time

TSQLMM2 <- lmer(sqrt(Predicted.count.) ~ allele*genotype + (1 | line) + (1 | block), data=Teggs)

res <- resid(TSQLMM2)
ran <- ranef(TSQLMM2)
fit <- fitted(TSQLMM2)
par(mfrow=c(1,3))
hist(ran$line[,1], main='')
plot(res ~ fit)
qqnorm(res, main='')
qqline(res)

summary(TSQLMM2)
anova(TSQLMM2)

wald.test(b=fixef(TSQLMM2), Sigma=vcov(TSQLMM2), Terms = 2, df=1)
wald.test(b=fixef(TSQLMM2), Sigma=vcov(TSQLMM2), Terms = 3, df=1)
wald.test(b=fixef(TSQLMM2), Sigma=vcov(TSQLMM2), Terms = 4, df=1)

### COMPARE THE TWO MODELS
AIC(TLMM2, TSQLMM2)
anova(TLMM2, TSQLMM2)

### again the sqrt appears justified

### can we compare the 4 models??
AIC(LMM1, SQLMM1, LMM2, SQLMM2)
anova(LMM1, SQLMM1, LMM2, SQLMM2)


### alternate random efffects structure as suggested b Max the line also dependant on the gneotpe as lines only have one genotype

TLMM3 <- lmer(Predicted.count. ~ allele*genotype + (1 + genotype | line) + (1 | block), data=Teggs)

res <- resid(TLMM3)
ran <- ranef(TLMM3)
fit <- fitted(TLMM3)
par(mfrow=c(1,3))
hist(ran$line[,1], main='')
plot(res ~ fit)
qqnorm(res, main='')
qqline(res)

summary(TLMM3)
anova(TLMM3)

wald.test(b=fixef(TLMM3), Sigma=vcov(TLMM3), Terms = 2, df=1)
wald.test(b=fixef(TLMM3), Sigma=vcov(TLMM3), Terms = 3, df=1)
wald.test(b=fixef(TLMM3), Sigma=vcov(TLMM3), Terms = 4, df=1)

### and again with the square rooting

TSQLMM3 <- lmer(sqrt(Predicted.count.) ~ allele*genotype + (1 + genotype| line) + (1 | block), data=Teggs)

res <- resid(TSQLMM3)
ran <- ranef(TSQLMM3)
fit <- fitted(TSQLMM3)
par(mfrow=c(1,3))
hist(ran$line[,1], main='')
plot(res ~ fit)
qqnorm(res, main='')
qqline(res)

summary(TSQLMM3)
anova(TSQLMM3)

wald.test(b=fixef(TSQLMM3), Sigma=vcov(TSQLMM3), Terms = 2, df=1)
wald.test(b=fixef(TSQLMM3), Sigma=vcov(TSQLMM3), Terms = 3, df=1)
wald.test(b=fixef(TSQLMM3), Sigma=vcov(TSQLMM3), Terms = 4, df=1)

### COMPARE THE TWO MODELS
AIC(TLMM3, TSQLMM3)
anova(TLMM3, TSQLMM3)

### justifys the squre rooting again 


### can we compare the 6 models??
AIC(LMM1, SQLMM1, LMM2, SQLMM2, TLMM3, TSQLMM3)
anova(LMM1, SQLMM1, LMM2, SQLMM2, TLMM3, TSQLMM3)

### sort off....



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
