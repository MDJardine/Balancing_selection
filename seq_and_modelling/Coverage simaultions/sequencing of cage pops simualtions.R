### binomial models of popution smapling and pool sequencing ####
setwd("/Users/michaeljardine/Desktop/DTP/Data and analysis/coverage_simulations")


# rbinom(n, size, prob)

# n = number of onservations - number of times that the process is repeated 
# size = number of trials - number of samples - number of chromosomes taken from the population
# prob = probability of success - the true allele frequency within our poputlations

## Max's example: S = 0.4, sample 48 flies (96 chro), repeat 1000 times
## for all fdreq of S 0.1 - 0.5 
NS.5 <- rbinom(n=1000, size=96, prob=0.5)
NS.4 <- rbinom(n=1000, size=96, prob=0.4)
NS.3 <- rbinom(n=1000, size=96, prob=0.3)
NS.2 <- rbinom(n=1000, size=96, prob=0.2)
NS.1 <- rbinom(n=1000, size=96, prob=0.1)

par(mfrow=c(1,5))
hist(NS.5/96)
hist(NS.4/96)
hist(NS.3/96)
hist(NS.2/96)
hist(NS.1/96)

## works

## coverage of x50 
NSS.5.50 <- rbinom(n=1000, size=50, prob=NS.5/96)
NSS.4.50 <- rbinom(n=1000, size=50, prob=NS.4/96)
NSS.3.50 <- rbinom(n=1000, size=50, prob=NS.3/96)
NSS.2.50 <- rbinom(n=1000, size=50, prob=NS.2/96)
NSS.1.50 <- rbinom(n=1000, size=50, prob=NS.1/96)

FSS.5.50 <- NSS.5.50/50
FSS.4.50 <- NSS.4.50/50
FSS.3.50 <- NSS.3.50/50
FSS.2.50 <- NSS.2.50/50
FSS.1.50 <- NSS.1.50/50

par(mfrow=c(1,5))
hist(FSS.5.50)
hist(FSS.4.50)
hist(FSS.3.50)
hist(FSS.2.50)
hist(FSS.1.50)

## repeat for 60
NSS.5.60 <- rbinom(n=1000, size=60, prob=NS.5/96)
NSS.4.60 <- rbinom(n=1000, size=60, prob=NS.4/96)
NSS.3.60 <- rbinom(n=1000, size=60, prob=NS.3/96)
NSS.2.60 <- rbinom(n=1000, size=60, prob=NS.2/96)
NSS.1.60 <- rbinom(n=1000, size=60, prob=NS.1/96)

FSS.5.60 <- NSS.5.60/60
FSS.4.60 <- NSS.4.60/60
FSS.3.60 <- NSS.3.60/60
FSS.2.60 <- NSS.2.60/60
FSS.1.60 <- NSS.1.60/60

par(mfrow=c(1,5))
hist(FSS.5.60)
hist(FSS.4.60)
hist(FSS.3.60)
hist(FSS.2.60)
hist(FSS.1.60)

## repeat for 70
NSS.5.70 <- rbinom(n=1000, size=70, prob=NS.5/96)
NSS.4.70 <- rbinom(n=1000, size=70, prob=NS.4/96)
NSS.3.70 <- rbinom(n=1000, size=70, prob=NS.3/96)
NSS.2.70 <- rbinom(n=1000, size=70, prob=NS.2/96)
NSS.1.70 <- rbinom(n=1000, size=70, prob=NS.1/96)

FSS.5.70 <- NSS.5.70/70
FSS.4.70 <- NSS.4.70/70
FSS.3.70 <- NSS.3.70/70
FSS.2.70 <- NSS.2.70/70
FSS.1.70 <- NSS.1.70/70

par(mfrow=c(1,5))
hist(FSS.5.70)
hist(FSS.4.70)
hist(FSS.3.70)
hist(FSS.2.70)
hist(FSS.1.70)

## repeat for 80
NSS.5.80 <- rbinom(n=1000, size=80, prob=NS.5/96)
NSS.4.80 <- rbinom(n=1000, size=80, prob=NS.4/96)
NSS.3.80 <- rbinom(n=1000, size=80, prob=NS.3/96)
NSS.2.80 <- rbinom(n=1000, size=80, prob=NS.2/96)
NSS.1.80 <- rbinom(n=1000, size=80, prob=NS.1/96)

FSS.5.80 <- NSS.5.80/80
FSS.4.80 <- NSS.4.80/80
FSS.3.80 <- NSS.3.80/80
FSS.2.80 <- NSS.2.80/80
FSS.1.80 <- NSS.1.80/80

par(mfrow=c(1,5))
hist(FSS.5.80)
hist(FSS.4.80)
hist(FSS.3.80)
hist(FSS.2.80)
hist(FSS.1.80)

## repeat for 90
NSS.5.90 <- rbinom(n=1000, size=90, prob=NS.5/96)
NSS.4.90 <- rbinom(n=1000, size=90, prob=NS.4/96)
NSS.3.90 <- rbinom(n=1000, size=90, prob=NS.3/96)
NSS.2.90 <- rbinom(n=1000, size=90, prob=NS.2/96)
NSS.1.90 <- rbinom(n=1000, size=90, prob=NS.1/96)

FSS.5.90 <- NSS.5.90/90
FSS.4.90 <- NSS.4.90/90
FSS.3.90 <- NSS.3.90/90
FSS.2.90 <- NSS.2.90/90
FSS.1.90 <- NSS.1.90/90

par(mfrow=c(1,5))
hist(FSS.5.90)
hist(FSS.4.90)
hist(FSS.3.90)
hist(FSS.2.90)
hist(FSS.1.90)

## repeat at 100 for 0.1 - 0.5
NSS.5.100 <- rbinom(n=1000, size=100, prob=NS.5/96)
NSS.4.100 <- rbinom(n=1000, size=100, prob=NS.4/96)
NSS.3.100 <- rbinom(n=1000, size=100, prob=NS.3/96)
NSS.2.100 <- rbinom(n=1000, size=100, prob=NS.2/96)
NSS.1.100 <- rbinom(n=1000, size=100, prob=NS.1/96)

FSS.5.100 <- NSS.5.100/100
FSS.4.100 <- NSS.4.100/100
FSS.3.100 <- NSS.3.100/100
FSS.2.100 <- NSS.2.100/100
FSS.1.100 <- NSS.1.100/100

par(mfrow=c(1,5))
hist(FSS.5.100)
hist(FSS.4.100)
hist(FSS.3.100)
hist(FSS.2.100)
hist(FSS.1.100)

### now take the results from some of these models and put into a data frame to make the plots

# start small with a data frame of x50 coverage
coveragex50 <- data.frame(FSS.5.50, FSS.4.50, FSS.3.50, FSS.2.50, FSS.1.50)

# not sure how to manipulate this data frame atm
# export to excel and sort there
setwd("/Users/michaeljardine/Desktop/DTP/Data and analysis/coverage_simulations")
write.csv(coveragex50, file = "coveragex50.csv")

combx50 <- read.csv("coveragex50_comb.csv")
str(combx50)
summary(combx50)

par(mfrow=c(1,1))
plot(SimfreqS ~ PopS, data=combx50)

# another one for x100
coveragex100 <- data.frame(FSS.5.100, FSS.4.100, FSS.3.100, FSS.2.100, FSS.1.100)

# not sure how to manipulate this data frame atm
# export to excel and sort there
setwd("/Users/michaeljardine/Desktop/DTP/Data and analysis/coverage_simulations")
write.csv(coveragex100, file = "coveragex100.csv")

combx100 <- read.csv("coveragex100_comb.csv")
str(combx100)
summary(combx100)

par(mfrow=c(1,1))
plot(SimfreqS ~ PopS, data=combx100)

# add coverage variable
x50 <- rep(50, 5000)
combx50$coverage <- x50

x100 <- rep(100, 5000)
combx100$coverage <- x100

# can we combine the two data sets?
comb_x50_x100 <- rbind(combx50, combx100)
str(comb_x50_x100)
comb_x50_x100$PopS <- factor(comb_x50_x100$PopS)
comb_x50_x100$coverage <- factor(comb_x50_x100$coverage)
str(comb_x50_x100)


## can we plot this?
library(ggplot2)
ggplot(aes(x=PopS, y=SimfreqS, color=coverage), data=comb_x50_x100) +
  geom_point(shape=1, position=position_jitter()) +
  scale_colour_hue(l=50) +
  geom_smooth(method = lm, se=TRUE, fullrange=TRUE)

library(plyr)
cdata <- ddply(comb_x50_x100, c("PopS", "coverage"), summarise,
               N    = length(SimfreqS),
               mean = mean(SimfreqS),
               sd   = sd(SimfreqS),
               se   = sd / sqrt(N),
               ci  = 1.96*(sd/(sqrt(N)))
)
cdata

## using se
pd <- position_dodge(0.2)
ggplot(aes(x=PopS, y=mean, color=coverage), data=cdata) + 
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.1, position=pd) +
  geom_line(position=pd) +
  geom_point(position=pd) +
  geom_smooth(method=lm, formula=y~x, se=TRUE, fullrange=TRUE, size=2)

x50vsx100 <- ggplot(aes(x=PopS, y=mean, color=coverage), data=cdata)+ 
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.3, position=pd)+
  geom_point(aes(size=1), position=pd)+
  ylab("Frequency of S determined by sequencing")+
  xlab("True population frequency of S")+
  theme_bw()+
  theme(legend.title=element_blank(),axis.text = element_text(size=25),
        axis.title = element_text(size=30),legend.position = "none",strip.text=element_text(size=25))+
  scale_y_continuous(limits=c(0, 0.7))

x50vsx100

## using sd
x50vsx100.sd <- ggplot(aes(x=PopS, y=mean, color=coverage), data=cdata)+ 
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.3, position=pd)+
  geom_point(aes(size=1), position=pd)+
  ylab("Frequency of S determined by sequencing")+
  xlab("True population frequency of S")+
  theme_bw()+
  theme(legend.title=element_blank(),axis.text = element_text(size=25),
        axis.title = element_text(size=30),legend.position = "none",strip.text=element_text(size=25))+
  scale_y_continuous(limits=c(0, 0.7))+
  theme(legend.position="right")

x50vsx100.sd


## using confidence intervals
x50vsx100.ci <- ggplot(aes(x=PopS, y=mean, color=coverage), data=cdata)+ 
  geom_errorbar(aes(ymin=mean-ci, ymax=mean+ci), width=.3, position=pd)+
  geom_point(aes(size=1), position=pd)+
  ylab("Frequency of S determined by sequencing")+
  xlab("True population frequency of S")+
  theme_bw()+
  theme(legend.title=element_blank(),axis.text = element_text(size=25),
        axis.title = element_text(size=30),legend.position = "none",strip.text=element_text(size=25))+
  scale_y_continuous(limits=c(0, 0.7))

x50vsx100.ci
