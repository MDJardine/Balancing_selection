#######################################
##### Blocks 1, 2 and 3 - ALL DATA ######



setwd("/Users/michaeljardine/Desktop/DTP/Data and analysis/Eggs")

TM6B <- read.csv("TM6B_comparison.csv")
str(TM6B)


### convert line to a factor with 6 levels
Teggs$line <- factor(Teggs$line)
str(Teggs)
### convert block to factor
Teggs$block <- factor(Teggs$block)
str(Teggs)
summary(Teggs)

### convert M and F alleles to long and short ( L and S) respectively
levels(Teggs$allele) <- c("S", "L")


### we have two numbers produced for the number eggs, predicted and corrected
x = TM6B$Predicted.count.
y = TM6B$Corrected.count.
plot(cor(x, y))
### since these correlate exactly, it dosen't matter which one we use
### I will use the predicted count for consistancy

### histograms to check for normality

par(mfrow=c(3,1))
hist(TM6B$Predicted.count.)
### bunched mostly towards zero try logging?
hist(log(TM6B$Predicted.count.))
### dosen't really help, now bunched towards the right - square rooting?
hist(sqrt(TM6B$Predicted.count.))
### that actuall looks really good


# some summary statistics for means and errors 
summarySE(Teggs, measurevar=c("Predicted.count."), groupvars=c("allele"))
summarySE(Teggs, measurevar=c("Predicted.count."), groupvars=c("line"))
summarySE(Teggs, measurevar=c("Predicted.count."), groupvars=c("genotype"))
summarySE(Teggs, measurevar=c("Predicted.count."), groupvars=c("block"))

### later models will have an intercation between allele and genotype so we will get summary stats for this too
Teggs$allelegenotype <- interaction(Teggs$allele, Teggs$genotype)
summarySE(Teggs, measurevar=c("Predicted.count."), groupvars=c("allelegenotype"))

Teggs$linegenotype <- interaction(Teggs$line, Teggs$genotype)
Teggs$BLG <- interaction(Teggs$linegenotype, Teggs$block)
summarySE(Teggs, measurevar=c("Predicted.count."), groupvars=c("BLG"))


summarySE(Teggs, measurevar=c("Predicted.count."), groupvars=c("linegenotype"))
Teggs$blockline <- interaction(Teggs$line, Teggs$block)
summarySE(Teggs, measurevar=c("Predicted.count."), groupvars=c("blockline"))


#### Models ####
### we want to look at the differenecs in the number of eggs laid between the lines of flies with different fruitless alleles
### flies also vary in which line they belong to (repliate of allele) and their genotype
### we also have 2 block which may var independantly so include as a random factor in all models
### 2 models: 1) egg number explaiend purely by allele; 2) include allele and genotype (B/D) and thier interaction

### these models will be linear mixed models

TLMM2 <- lm(Predicted.count. ~ allele*genotype, data=Teggs)

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
Anova(TLMM2)

### Appears to be an effect of the allele on the number of eggs laid

### sqrt looked good on hist before and diagnostic plots could be better
### helped improve all the block when looked at individually

TSQLMM2 <- lm(sqrt(Predicted.count.) ~ allele*genotype, data=Teggs)

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
Anova(TSQLMM2)

wald.test(b=fixef(TSQLMM1), Sigma=vcov(TSQLMM1), Terms = 2, df=1)

### now allele does not effect the number of eggs laid


### looks the smae really
AIC(TLMM1, TSQLMM1)
anova(TLMM1, TSQLMM1)
## sqrt justified if looking at these



### we now want to account for the variation caused by the fly's genotpe as well

TLMM2 <- lm(Predicted.count. ~ allele*genotype + (1 | line) + (1 | block), data=Teggs)

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
Anova(TLMM2)

### all of allele, genotype and thier intercation seem to affect the number of eggs laid 


### lets see what sqrt does this time

TSQLMM2 <- lm(sqrt(Predicted.count.) ~ allele*genotype + (1 | line) + (1 | block), data=Teggs)

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
Anova(TSQLMM2)

### when using sqrt allele no longer affects the number of eggs



### COMPARE THE TWO MODELS
AIC(TLMM2, TSQLMM2)
anova(TLMM2, TSQLMM2)

### again the sqrt appears justified even id not giving the answer we want


### alternate random efffects structure as suggested by Max the line also dependant on the allele, as the structure of the experiment is not balanced
### is this the corrc way to do this, I am very confused with how I should (or if I should) acount for this

TLMM3 <- lmer(Predicted.count. ~ allele*genotype + (1 + allele | line) + (1 | block), data=Teggs)
### warning message

res <- resid(TLMM3)
ran <- ranef(TLMM3)
fit <- fitted(TLMM3)
par(mfrow=c(1,3))
hist(ran$line[,1], main='')
plot(res ~ fit)
qqnorm(res, main='')
qqline(res)

summary(TLMM3)
anova(TLMM3, test="F")
Anova(TLMM3)



### seems to result in a much stronger afect of the allele on the number of eggs laid
### and again with the square rooting

TSQLMM3 <- lmer(sqrt(Predicted.count.) ~ allele*genotype + (1 + allele | line) + (1 | block), data=Teggs)
# warning message

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
Anova(TSQLMM3)


### as before the affect of the flies fru allele drops out when using sqrt

### COMPARE THE TWO MODELS
AIC(TLMM3, TSQLMM3)
anova(TLMM3, TSQLMM3)



### parametirc bootstraping ...... ####

### multiple models

EG.PB1 <- lmer(Predicted.count. ~ allele + (1 | line) + (1 | block), data=Teggs )
EG.PB2 <- update(EG.PB1, .~. -allele)
ParB1 <- PBmodcomp(EG.PB1, EG.PB2, nsim=500)
summary(ParB1)
### sems to have worked, not sure what it means



EG.PB3 <- lmer(Predicted.count. ~ genotype + (1 | line) + (1 | block), data=Teggs) 
EG.PB4 <- update(EG.PB4, .~. -genotype)
ParB2 <- PBmodcomp(EG.PB3, EG.PB4, nsim=500)
summary(ParB2)
### again seems to work, not sure what it means


### more complicated model

EG.PB3 <- lmer(Predicted.count. ~ allele*genotype + (1 | line) + (1 | block), data=Teggs) 
### drop genotype
ParB3 <- PBmodcomp(EG.PB5, EG.PB4, nsim=500)
summary(ParB3)
### ERROR!

### drop allele
ParB4 <- PBmodcomp(EG.PB5, EG.PB2, nsim=500)
summary(ParB4)
### ERROR!

### so when trying to make models more complicated there are errors!

### also all these models have different random effects structures to those models I used above to look at the effects they have on egg production
### this is because when we use these models witht he parametirc bootstrap it gets really messy
### e.g.

TLMM3.drop.allele <- update(TLMM3, .~. -allele)
### warning but ok

ParB5 <- PBmodcomp(TLMM3, TLMM3.drop.allele, nsim=500)


### waning mesages go on for ever!


### so the big problems seems to be: 
### 1) what random error structure are we using? what is appropriate for the set up of our experiment? and are we happy to use this treatment even if we loose the signifcance of the terms we thought we had
### 2) how can the parametric bootsrap be used to produce significance test for more complicated models?
### 3) what does the output tell us?


#### plots ####


### eggs b allele
qplot( x=allele , y=Predicted.count. , data=Teggs , geom=c("boxplot"),
       xlab = "Fruitless Allele", ylab = "Total number of eggs") +
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold")) + 
  stat_summary(fun.y=mean, geom="point", shape=20, size=10, color="blue", fill="blue")


qplot( x=allele , y=Predicted.count. , data=Teggs , geom=c("boxplot", "jitter"), fill=allele,
       xlab = "Fruitless Allele", ylab = "Total number of eggs")  +
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold")) + 
  stat_summary(fun.y=mean, geom="point", shape=20, size=10, color="blue", fill="blue")


### eggs by line
qplot( x=line , y=Predicted.count. , data=Teggs , geom=c("boxplot"), fill=allele,
       xlab = "Line", ylab = "Total number of eggs")  +
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold")) + 
  stat_summary(fun.y=mean, geom="point", shape=20, size=10, color="blue", fill="blue")

qplot( x=line , y=Predicted.count. , data=Teggs , geom=c("boxplot", "jitter"), fill=allele,
       xlab = "Line", ylab = "Total number of eggs") +
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold")) + 
  stat_summary(fun.y=mean, geom="point", shape=20, size=10, color="blue", fill="blue")


### allele and genotype
qplot( x=allelegenotype , y=Predicted.count. , data=Teggs , geom=c("boxplot"), fill=allele,
       xlab = "Allele / Genotype combination", ylab = "Total number of eggs") +
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold")) + 
  stat_summary(fun.y=mean, geom="point", shape=20, size=10, color="blue", fill="blue")

qplot( x=allelegenotype , y=Predicted.count. , data=Teggs , geom=c("boxplot", "jitter"), fill=allele,
       xlab = "Allele / Genotype combination", ylab = "Total number of eggs") +
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold")) + 
  stat_summary(fun.y=mean, geom="point", shape=20, size=10, color="blue", fill="blue")


ggplot(aes(y=Predicted.count., x=allele), data=TM6B) +
  + geom_boxplot()
+ stat_summary(fun.y=mean, geom="point", shape=1, size=4) +
  + theme(axis.title.y = element_text(size=18),
          + axis.title.x = element_text(size=18),
          + axis.text.x=element_text(size=15, colour="black"),
          + axis.text.y=element_text(size=15)) +
  + scale_y_continuous(name="Total number of eggs laid") +
  + scale_x_discrete(name="Fruitless Allele", limits=c("F", "M"), labels=c("F", "M"))

### lines split by genotype
qplot( x=linegenotype , y=Predicted.count. , data=Teggs , geom=c("boxplot"), fill=allele,
       xlab = "Line", ylab = "Total number of eggs")  +
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold")) + 
  stat_summary(fun.y=mean, geom="point", shape=20, size=10, color="blue", fill="blue")


bp7 <- ggplot(aes(y=Predicted.count., x=line, fill=genotype), data=Teggs) +
  geom_boxplot()
bp7
bp7 + theme(axis.title.y = element_text(size=18),
            axis.title.x = element_text(size=18),
            axis.text.x=element_text(size=15, colour="black"),
            axis.text.y=element_text(size=15)) +
  scale_y_continuous(name="Total number of eggs laid") +
  scale_x_discrete(name="Line") +
  scale_fill_manual(values=c("grey80", "goldenrod2"))



### box plots in filip's style
ggplot(aes(y=Predicted.count., x=line, col=allele), data=Teggs)+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(width=0.3)+
  ylab("Number of eggs laid")+
  xlab("Line")+
  theme_bw()+
  theme(legend.title=element_blank(),axis.text = element_text(size=25),axis.title = element_text(size=25),legend.position = "none",strip.text=element_text(size=25))+
  facet_wrap(~allele,scales="free_x")



ggplot(aes(y=Predicted.count., x=linegenotype, col=allele), data=Teggs)+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(width=0.3)+
  ylab("Number of eggs laid")+
  xlab("Line")+
  theme_bw()+
  theme(legend.title=element_blank(),axis.text = element_text(size=25),axis.title = element_text(size=25),legend.position = "none",strip.text=element_text(size=25))+
  facet_wrap(~genotype,scales="free_x") +
  scale_x_discrete(labels=c("M05", "M08", "M16", "M19", "M24", "M95"))




##### finding means for line/block comparisons ####

TeggsA <- read.csv("All_egg_counts_blocks_1-3_combine blocks3-4.csv")
str(Teggs)


### convert line to a factor with 6 levels
TeggsA$line <- factor(Teggs$line)
str(TeggsA)
### convert block to factor
TeggsA$block <- factor(TeggsA$block)
str(TeggsA)
summary(TeggsA)

### convert M and F alleles to long and short ( L and S) respectively
levels(TeggsA$allele) <- c("S", "L")

TeggsA$blockline <- interaction(TeggsA$line, TeggsA$block)
summarySE(TeggsA, measurevar=c("Predicted.count."), groupvars=c("blockline"))


