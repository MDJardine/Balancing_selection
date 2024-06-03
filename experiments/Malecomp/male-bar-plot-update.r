setwd("/Users/michaeljardine/Desktop/poster stuff to transfer")

interbars <- read.csv("INTER BARS.csv")
interbars$Line

thing2 <- ggplot(aes(y=winfrac, x=Line, fill=Genotype), data=interbars) +
  geom_bar(stat="identity", colour="black", position="dodge") +
  geom_errorbar(aes(ymin=winfrac-se, ymax=winfrac+se),
                width=.2,
                position=position_dodge(.9)) +
  theme_bw()+
  theme(legend.title=element_blank(),axis.text = element_text(size=30),axis.title = element_text(size=30),legend.position = "none",strip.text=element_text(size=30))+
  scale_y_continuous(name="Proportion of matings by focals")+
  theme(legend.position="right")
thing2






C.end <- ggplot(aes(y=winfrac, x=Line, fill=allele), data=interbars)+
  geom_bar(stat="identity", colour="black", position="dodge")+
  geom_errorbar(aes(ymin=winfrac-se, ymax=winfrac+se),
                width=.2,
                position=position_dodge(.9))+
  ylab("Number of eggs laid")+
  xlab("Line")+
  theme_bw()+
  theme(legend.title=element_blank(),axis.text = element_text(size=25),
        axis.title = element_text(size=30),legend.position = "none",strip.text=element_text(size=25))+
  theme(legend.position="right")+
  scale_fill_manual(values=c("#00BFC4", "#F8766D"))+
  facet_wrap(~Genotype,scales="free_x")
C.end

### add lines
C.dat_vlinesS <- data.frame(Genotype=c("B", "D"), xval=c(0.4687, 0.4))
C.dat_vlinesL <- data.frame(Genotype=c("B", "D"), xval=c(0.4695, 0.3))

C.end2 <- ggplot(aes(y=winfrac, x=Line, fill=allele), data=interbars)+
  geom_bar(stat="identity", colour="black", position="dodge")+
  geom_errorbar(aes(ymin=winfrac-se, ymax=winfrac+se),
                width=.2,
                position=position_dodge(.9))+
  geom_hline(aes(yintercept=xval), data=C.dat_vlinesS, col="#F8766D", size=1.5, linetype=2)+
  geom_hline(aes(yintercept=xval), data=C.dat_vlinesL, col="#00BFC4", size=1.5, linetype=3)+
  ylab("Proportion of mating  focal males")+
  xlab("Line")+
  theme_bw()+
  theme(legend.title=element_blank(),axis.text = element_text(size=25),
        axis.title = element_text(size=30),legend.position = "none",strip.text=element_text(size=25))+
  theme(legend.position="right")+
  scale_fill_manual(values=c("#00BFC4", "#F8766D"))+
  facet_wrap(~Genotype,scales="free_x")
C.end2






#####ggplot(aes(x=line, y=winfrac, fill=geno), data=interbars)+
geom_errorbar(aes(ymin=winfrac-se, ymax=winfrac+se),
              + width=.2,
              + position=position_dodge(.9))+
  theme_bw()
theme(legend.title=element_blank(),axis.text = element_text(size=25),axis.title = element_text(size=25),legend.position = "none",strip.text=element_text(size=25))+
  ylab("Proportion of matings by focal males")+
  xlab("Line")

+
  theme_bw()

+
  + geom_bar(position=position_dodge(), stat="identity", colour="black")+
  
  
  +
  + theme(axis.title.y = element_text(size=16),
          + axis.title.x = element_text(size=18),
          + axis.text.x=element_text(size=15, colour="black"),
          + axis.text.y=element_text(size=15)) +
  + theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
          + panel.background=element_blank(), axis.line=element_line(colour="black")) +
  + scale_y_continuous(name="")

+
  + scale_x_discrete(labels=c("1", "2", "3", "4", "5", " ", " ", " ", "6", "7", "8", "9", "10")) +
  + scale_fill_manual(values=c("red1", "cornflowerblue"))
