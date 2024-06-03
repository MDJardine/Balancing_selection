### final plots
setwd("")

Teggs <- read.csv("All_egg_counts_blocks_1-3.csv")


E.end <- ggplot(aes(y=Predicted.count.,x=line,fill=allele), data=Teggs)+
  geom_boxplot(outlier.shape = NA)+
  ylab("Number of eggs laid")+
  xlab("Line")+
  theme_bw()+
  theme(legend.title=element_blank(),axis.text = element_text(size=25),
        axis.title = element_text(size=30),legend.position = "none",strip.text=element_text(size=25))+
  theme(legend.position="right")+
  scale_fill_manual(values=c("gold", "firebrick3"))+
  facet_wrap(~genotype,scales="free_x")
E.end

### add lines
egglines <- read.csv("linesdf.csv")
egglines$linegenotype <- interaction(egglines$Algen, egglines$geno)


dat_vlinesS <- data.frame(genotype=c("B", "D"), xval=c(26.03, 29.67))
dat_vlinesL <- data.frame(genotype=c("B", "D"), xval=c(23.52, 23.57))


E.end2 <- ggplot(aes(y=Predicted.count.,x=line,fill=allele), data=Teggs)+
  geom_boxplot(outlier.shape = NA)+
  geom_hline(aes(yintercept=xval), data=dat_vlinesS, size=1.5, linetype=1)+
  geom_hline(aes(yintercept=xval), data=dat_vlinesL, size=1.5, linetype=3)+
  ylab("Number of eggs laid")+
  xlab("Line")+
  theme_bw()+
  theme(legend.title=element_blank(),axis.text = element_text(size=25),
        axis.title = element_text(size=25),legend.position = "none",strip.text=element_text(size=25))+
  theme(legend.position="none")+
  scale_fill_manual(values=c("#CA3542", "#AECBC9"))+
  facet_wrap(~genotype,scales="free_x")
E.end2








