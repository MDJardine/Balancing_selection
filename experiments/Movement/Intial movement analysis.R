setwd("/Users/michaeljardine/Desktop/DTP/Data and analysis/Movement")

move <- read.csv("initial movement assay.csv")
str(move)


### convert line to a factor with 6 levels
move$Vial <- factor(move$Vial)
str(move)

summary(move)

par(mfrow=c(3,1))
### histograms of mean time
hist(move$Mean.time)
# logging
hist(log(move$Mean.time))
#square rooting
hist(sqrt(move$Mean.time))



## compare the two lines of flies for their times to reach the top of the vial as a measure of movement

m1 <- t.test(move$Mean.time ~ move$Line)
summary(m1)
m1

### linear model 

m2 <- lm(Mean.time ~ Line, data=move)

par(mfrow=c(2,2))
plot(m2)


par(mfrow=c(1,1))
boxplot(Mean.time ~ Line, data=move, notch=T)

library(ggplot2)

ggplot(aes(y=Mean.time, x=Line, fill=Line), data=move)+
  geom_boxplot()+
  ylab("Time to reach top (s)")+
  xlab("Line")+
  theme_bw()+
  theme(legend.title=element_blank(),axis.text = element_text(size=25),axis.title = element_text(size=25),legend.position = "none",strip.text=element_text(size=25))
