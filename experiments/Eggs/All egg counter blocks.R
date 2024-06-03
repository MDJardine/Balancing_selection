#######################################
##### Blocks 1, 2 and 3 - ALL DATA ######

### the data set will be termed 'Teggs' for total eggs
### we first need to correct the parts of the data set for factors and plot the data for patterns 
### here the secod round of block 3 will be treated as 'block 4' for ease at the moment

setwd("/Users/michaeljardine/Desktop/DTP/Datasets/Eggs")

Teggs <- read.csv("All_egg_counts_blocks_1-3.csv")
str(Teggs)


### convert line to a factor with 6 levels
Teggs$line <- factor(Teggs$line)
str(Teggs)
### convert block to factor
Teggs$block <- factor(Teggs$block)
str(Teggs)
summary(Teggs)

### we have two numbers produced for the number eggs, predicted and corrected
x = eggs2$Predicted.count.
y = eggs2$Corrected.count.
cor(x, y)
### since these correlate exactly, it dosen't matter which one we use
### I will use the predicted count for consistancy

### histograms to check for normality

par(mfrow=c(3,1))
hist(Teggs$Predicted.count.)
### bunched mostly towards zero try logging?
hist(log(Teggs$Predicted.count.))
### dosen't really help, now bunched towards the right - square rooting?
hist(sqrt(Teggs$Predicted.count.))
### that actuall looks really good


# some summary statistics for means and errors 
summarySE(Teggs, measurevar=c("Predicted.count."), groupvars=c("allele"))
summarySE(Teggs, measurevar=c("Predicted.count."), groupvars=c("line"))
summarySE(Teggs, measurevar=c("Predicted.count."), groupvars=c("genotype"))
summarySE(Teggs, measurevar=c("Predicted.count."), groupvars=c("block"))

### later models will have an intercation between allele and genotype so we will get summary stats for this too
Teggs$allelegenotype <- interaction(Teggs$allele, Teggs$genotype)
summarySE(Teggs, measurevar=c("Predicted.count."), groupvars=c("allelegenotype"))

### Models
### we want to look at the differenecs in the number of eggs laid between the lines of flies with different fruitless alleles
### flies also vary in which line they belong to (repliate of allele) and their genotype
### we also have 2 block which may var independantly so include as a random factor in all models
### 2 models: 1) egg number explaiend purely by allele; 2) include allele and genotype (B/D) and thier interaction

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

### Appears to be an effect of the allele on the number of eggs laid

### sqrt looked good on hist before and diagnostic plots could be better
### helped improve all the block when looked at individually

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

### now allele does not effect the number of eggs laid


### looks the smae really
AIC(TLMM1, TSQLMM1)
anova(TLMM1, TSQLMM1)
## sqrt justified if looking at these



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
TLMM2
### reasonable p-values which ma mean somethign depending on the degrees of freedom

wald.test(b=fixef(TLMM2), Sigma=vcov(TLMM2), Terms = 2, df=1)
wald.test(b=fixef(TLMM2), Sigma=vcov(TLMM2), Terms = 3, df=1)
wald.test(b=fixef(TLMM2), Sigma=vcov(TLMM2), Terms = 4, df=1)

### wald tests seem no longer to work when there is more than 1 term - they did 2 years ago?
### may need to do pbktest to get significance tests for the models 

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

### Decrease in the F value of both the allele and allele x genotype terms

wald.test(b=fixef(TSQLMM2), Sigma=vcov(TSQLMM2), Terms = 2, df=1)
wald.test(b=fixef(TSQLMM2), Sigma=vcov(TSQLMM2), Terms = 3, df=1)
wald.test(b=fixef(TSQLMM2), Sigma=vcov(TSQLMM2), Terms = 4, df=1)
### no idea what these mean now .... (remove?)



### COMPARE THE TWO MODELS
AIC(TLMM2, TSQLMM2)
anova(TLMM2, TSQLMM2)

### again the sqrt appears justified


### alternate random efffects structure as suggested by Max the line also dependant on the gneotpe as lines only have one genotype

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

### seems to result in a much stronger afect of the allele on the number of eggs laid
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

### again a decrease in F values with square rooting with the square rooting

### COMPARE THE TWO MODELS
AIC(TLMM3, TSQLMM3)
anova(TLMM3, TSQLMM3)



### parametirc bootstraping ......


ParB1 <- PBmodcomp(TLMM3, TSQLMM3, nsim=500)

summary(ParB1)
