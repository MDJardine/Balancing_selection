setwd("/Users/michaeljardine/Desktop/DTP/Data and analysis/fruitless/Malecomp")


library(ggplot2)

interbars <- read.csv("INTER BARS.csv")

### basic plot 
thing2 <- ggplot(aes(y=winfrac, x=line, fill=geno), data=interbars) +
  geom_bar(stat="identity", colour="black", position="dodge") +
  geom_errorbar(aes(ymin=winfrac-se, ymax=winfrac+se),
                width=.2,
                position=position_dodge(.9)) +
  theme_bw()+
  theme(legend.title=element_blank(),axis.text = element_text(size=30),axis.title = element_text(size=30),legend.position = "none",strip.text=element_text(size=30))+
  scale_y_continuous(name="Proportion of matings by focals")+
  theme(legend.position="right")
thing2


### use faceting to split the B and D flies
C.end <- ggplot(aes(y=winfrac, x=line, fill=allele), data=interbars)+
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
  facet_wrap(~geno,scales="free_x")
C.end

### add lines plus change colour an legend
C.dat_vlinesS <- data.frame(geno=c("B", "D"), xval=c(0.4687, 0.4))
C.dat_vlinesL <- data.frame(geno=c("B", "D"), xval=c(0.4695, 0.3))

C.end2 <- ggplot(aes(y=winfrac, x=line, fill=allele), data=interbars)+
  geom_bar(stat="identity", colour="black", position="dodge")+
  geom_errorbar(aes(ymin=winfrac-se, ymax=winfrac+se),
                width=.2,
                position=position_dodge(.9))+
  geom_hline(aes(yintercept=xval), data=C.dat_vlinesS, size=1.5, linetype=2)+
  geom_hline(aes(yintercept=xval), data=C.dat_vlinesL, size=1.5, linetype=3)+
  ylab("Proportion of mating  focal males")+
  xlab("Line")+
  theme_bw()+
  theme(legend.title=element_blank(),axis.text = element_text(size=25),
        axis.title = element_text(size=30),legend.position = "none",strip.text=element_text(size=25))+
  theme(legend.position="none")+
  scale_fill_manual(values=c("firebrick3", "gold"))+
  facet_wrap(~geno,scales="free_x")
C.end2
ggsave("maleplot_red_yellow.png")

maleplot <- ggplot(aes(y=winfrac, x=line, fill=allele), data=interbars)+
  geom_bar(stat="identity", colour="black", position="dodge")+
  geom_errorbar(aes(ymin=winfrac-se, ymax=winfrac+se),
                width=.2,
                position=position_dodge(.9))+
  geom_hline(aes(yintercept=xval), data=C.dat_vlinesS, size=1.5, linetype=2)+
  geom_hline(aes(yintercept=xval), data=C.dat_vlinesL, size=1.5, linetype=3)+
  ylab("Proportion of mating  focal males")+
  xlab("Line")+
  theme_bw()+
  theme(axis.text = element_text(size=30),
        axis.title = element_text(size=40),legend.position = "none",strip.text=element_text(size=30))+
  scale_fill_manual(values=c("#CA3542", "#AECBC9"))+
  facet_wrap(~geno,scales="free_x")
maleplot
ggsave("male_plot_red_blue.png")


maleplot <- ggplot(aes(y=winfrac, x=line, fill=allele), data=interbars)+
  geom_bar(stat="identity", colour="black", position="dodge")+
  geom_errorbar(aes(ymin=winfrac-se, ymax=winfrac+se),
                width=.2,
                position=position_dodge(.9))+
  geom_hline(aes(yintercept=xval), data=C.dat_vlinesS, size=1.5, linetype=2)+
  geom_hline(aes(yintercept=xval), data=C.dat_vlinesL, size=1.5, linetype=3)+
  ylab("Proportion of matings obtained by focal males")+
  xlab("Line")+
  theme_bw()+
  theme(axis.text = element_text(size=30),
        axis.title = element_text(size=40),legend.position = "none",strip.text=element_text(size=30))+
  scale_fill_manual(values=c("#CA3542", "#AECBC9"))+
  facet_wrap(~geno,scales="free_x")
maleplot
ggsave("male_plot_red_blue.png")


maleplot <- ggplot(aes(y=winfrac, x=line, fill=allele), data=interbars)+
  geom_bar(stat="identity", colour="black", position="dodge")+
  geom_errorbar(aes(ymin=winfrac-se, ymax=winfrac+se),
                width=.2,
                position=position_dodge(.9))+
  geom_hline(aes(yintercept=xval), data=C.dat_vlinesS, size=1.5, linetype=2)+
  geom_hline(aes(yintercept=xval), data=C.dat_vlinesL, size=1.5, linetype=3)+
  ylab("Proportion of matings obtained by focal males")+
  xlab("Line")+
  theme_bw()+
  theme(axis.text = element_text(size=30),
        axis.title = element_text(size=40),legend.position = "none",strip.text=element_text(size=30))+
  scale_fill_manual(values=c("#CA3542", "#AECBC9"))+
  facet_wrap(~geno,scales="free_x")
maleplot


C.dat_vlinesS <- data.frame(geno=c("B", "D"), xval=c(0.4687, 0.4))
C.dat_vlinesL <- data.frame(geno=c("B", "D"), xval=c(0.4695, 0.3))

dfS <- data.frame(geno=c("B", "D"), x1=c(3.6, 3.6), x2 = c(6.5, 6.5), y1=c(0.4687, 0.4), y2=c(0.4687, 0.4))
dfL <- data.frame(geno=c("B", "D"), x1=c(0.5, 0.5), x2 = c(3.4, 3.4), y1=c(0.4695, 0.3), y2=c(0.4695, 0.3))

maleplot_segment <- ggplot(aes(y=winfrac, x=line, fill=allele), data=interbars)+
  geom_bar(stat="identity", colour="black", position="dodge")+
  geom_errorbar(aes(ymin=winfrac-se, ymax=winfrac+se),
                width=.2,
                position=position_dodge(.9))+
  geom_segment(aes(x=x1, y=y1, xend=x2, yend=y2), data=dfS, inherit.aes = FALSE, linetype=2, size=1.5)+
  geom_segment(aes(x=x1, y=y1, xend=x2, yend=y2), data=dfL, inherit.aes = FALSE, linetype=2, size=1.5)+
  ylab("Proportion of matings by focal males")+
  xlab("Line")+
  theme_bw()+
  theme(axis.text = element_text(size=30),
        axis.title.x = element_text(size=45), axis.title.y = element_text(size=45), legend.position = "none",strip.text=element_text(size=30))+
  theme(axis.title.x = element_text(vjust=-0.4), axis.title.y = element_text(vjus =2.5, hjust = 0.5))+
  theme(plot.margin = unit(c(0.5,0.1,0.7,0.7), "cm"))+
  scale_fill_manual(values=c("#CA3542", "#AECBC9"))+
  facet_wrap(~geno,scales="free_x")
maleplot_segment 
ggsave("Plots/maleplot_red_blue_segment_lines.png")


## reviewer suggestion to box plot like the other plots


dfS <- data.frame(geontype=c("B", "D"), x1=c(3.6, 3.6), x2 = c(6.5, 6.5), y1=c(0.6324, 0.5466), y2=c(0.6324, 0.5466))
dfL <- data.frame(genotype=c("B", "D"), x1=c(0.5, 0.5), x2 = c(3.4, 3.4), y1=c(0.5801, 0.4668), y2=c(0.5801, 0.4668))

maleplot_segment <- ggplot(aes(y=winner, x=Line, fill=Allele), data=Tcomp)+
  geom_boxplot(outlier.shape = NA)+
  geom_segment(aes(x=x1, y=y1, xend=x2, yend=y2), data=dfS, inherit.aes = FALSE, linetype=2, size=1.5)+
  geom_segment(aes(x=x1, y=y1, xend=x2, yend=y2), data=dfL, inherit.aes = FALSE, linetype=2, size=1.5)+
  ylab("Proportion of matings by focal males")+
  xlab("Line")+
  theme_bw()+
  theme(axis.text = element_text(size=30),
        axis.title.x = element_text(size=45), axis.title.y = element_text(size=45), legend.position = "none",strip.text=element_text(size=30))+
  theme(axis.title.x = element_text(vjust=-0.4), axis.title.y = element_text(vjus =2.5, hjust = 0.5))+
  theme(plot.margin = unit(c(0.5,0.1,0.7,0.7), "cm"))+
  scale_fill_manual(values=c("#CA3542", "#AECBC9"))+
  facet_wrap(~genotype,scales="free_x")
maleplot_segment 
ggsave("Plots/maleplot_red_blue_segment_lines.png")
## well that looks awful


## Other reviewer comment to use outlines of plots ratehr than fill
dfS <- data.frame(geno=c("B", "D"), x1=c(3.6, 3.6), x2 = c(6.5, 6.5), y1=c(0.4687, 0.4), y2=c(0.4687, 0.4))
dfL <- data.frame(geno=c("B", "D"), x1=c(0.5, 0.5), x2 = c(3.4, 3.4), y1=c(0.4695, 0.3), y2=c(0.4695, 0.3))

maleplot_segment_2 <- ggplot(aes(y=winfrac, x=line, color=allele), data=interbars)+
  geom_bar(stat="identity", fill="white", position="dodge")+
  geom_errorbar(aes(ymin=winfrac-se, ymax=winfrac+se),
                width=.2,
                position=position_dodge(.9), color="black")+
  geom_segment(aes(x=x1, y=y1, xend=x2, yend=y2), data=dfS, inherit.aes = FALSE, linetype=2, size=1.5)+
  geom_segment(aes(x=x1, y=y1, xend=x2, yend=y2), data=dfL, inherit.aes = FALSE, linetype=2, size=1.5)+
  ylab("Proportion of matings by focal males")+
  xlab("Line")+
  theme_bw()+
  theme(axis.text = element_text(size=30),
        axis.title.x = element_text(size=45), axis.title.y = element_text(size=45), legend.position = "none",strip.text=element_text(size=30))+
  theme(axis.title.x = element_text(vjust=-0.4), axis.title.y = element_text(vjus =2.5, hjust = 0.5))+
  theme(plot.margin = unit(c(0.5,0.1,0.7,0.7), "cm"))+
  scale_color_manual(values=c("#CA3542", "#AECBC9"))+
  facet_wrap(~geno,scales="free_x")
maleplot_segment_2
ggsave("Plots/maleplot_red_blue_segment_lines_2.png")

##again colours are quite weak
maleplot_segment_2.alt <- ggplot(aes(y=winfrac, x=line, color=allele), data=interbars)+
  geom_bar(stat="identity", fill="white", position="dodge", size=1.2)+
  geom_errorbar(aes(ymin=winfrac-se, ymax=winfrac+se),
                width=.2,
                position=position_dodge(.9), color="black")+
  geom_segment(aes(x=x1, y=y1, xend=x2, yend=y2), data=dfS, inherit.aes = FALSE, linetype=2, size=1.5)+
  geom_segment(aes(x=x1, y=y1, xend=x2, yend=y2), data=dfL, inherit.aes = FALSE, linetype=2, size=1.5)+
  ylab("Proportion of matings by focal males")+
  xlab("Line")+
  theme_bw()+
  theme(axis.text = element_text(size=30),
        axis.title.x = element_text(size=45), axis.title.y = element_text(size=45), legend.position = "none",strip.text=element_text(size=30))+
  theme(axis.title.x = element_text(vjust=-0.4), axis.title.y = element_text(vjus =2.5, hjust = 0.5))+
  theme(plot.margin = unit(c(0.5,0.1,0.7,0.7), "cm"))+
  scale_color_manual(values=c("#CA3542", "#7a8ed9"))+
  facet_wrap(~geno,scales="free_x")
maleplot_segment_2.alt
ggsave("Plots/maleplot_red_blue_segment_lines_2.alt.png")
