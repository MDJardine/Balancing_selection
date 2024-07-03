### plot code from filip
ggplot(aes(y=Predicted.count.,x=line,col=genotype), data=Teggs)+
  geom_boxplot(outlier.shape = NA)+
  ylab("Number of eggs laid")+
  xlab("Line")+
  theme_bw()+
  theme(legend.title=element_blank(),axis.text = element_text(size=25),axis.title = element_text(size=25),legend.position = "none",strip.text=element_text(size=25))+
  scale_x_discrete(labels=c("L1", "L2", "L3", "S1", "S2", "S3"))+
  theme(legend.position="right")



### final plots

E.end <- ggplot(aes(y=Predicted.count.,x=line,fill=allele), data=Teggs)+
  geom_boxplot(outlier.shape = NA)+
  ylab("Number of eggs laid")+
  xlab("Line")+
  theme_bw()+
  theme(legend.title=element_blank(),axis.text = element_text(size=25),
        axis.title = element_text(size=30),legend.position = "none",strip.text=element_text(size=25))+
  theme(legend.position="right")+
  facet_wrap(~genotype,scales="free_x")
E.end

### add lines

egglines <- read.csv("linesdf.csv")
egglines$linegenotype <- interaction(egglines$Algen, egglines$geno)


dat_vlinesS <- data.frame(genotype=c("B", "D"), xval=c(26.03, 29.67))
dat_vlinesL <- data.frame(genotype=c("B", "D"), xval=c(23.52, 23.57))


E.end2 <- ggplot(aes(y=Predicted.count.,x=line,fill=allele), data=Teggs)+
  geom_boxplot(outlier.shape = NA)+
  geom_hline(aes(yintercept=xval), data=dat_vlinesS, col="#F8766D", size=1.5, linetype=2)+
  geom_hline(aes(yintercept=xval), data=dat_vlinesL, col="#00BFC4", size=1.5, linetype=3)+
  ylab("Number of eggs laid")+
  xlab("Line")+
  theme_bw()+
  theme(legend.title=element_blank(),axis.text = element_text(size=25),
        axis.title = element_text(size=30),legend.position = "none",strip.text=element_text(size=25))+
  theme(legend.position="right")+
  facet_wrap(~genotype,scales="free_x")
E.end2




??linetype




### plot codes from masters

bp7 <- ggplot(aes(y=total, x=population, fill=allele), data=Teggs) +
  geom_boxplot()



+
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        panel.background=element_blank(), axis.line=element_line(colour="black"))
bp7 + theme(axis.title.y = element_text(size=18),
            axis.title.x = element_text(size=18),
            axis.text.x=element_text(size=15, colour="black"),
            axis.text.y=element_text(size=15)) +
  scale_y_continuous(name="Total number of eggs laid") +
  scale_x_discrete(name="Population") +
  scale_fill_manual(values=c("grey80", "goldenrod2"))
