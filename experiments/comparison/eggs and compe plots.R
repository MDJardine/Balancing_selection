#### eggs/block comapred to comp/block ####

library(lme4)
library(ggplot2)
library(RColorBrewer)
library(car)

setwd("/Users/michaeljardine/Desktop/DTP/Data and analysis/fruitless/comparison")

## colours for these plots come from colour brewer - which we find the hexcodes from colur palates 
display.brewer.all()
display.brewer.pal(n = 9, name = 'YlOrRd')
brewer.pal(n = 9, name = "YlOrRd")
display.brewer.pal(n = 9, name = 'Blues')
brewer.pal(n = 9, name = "Blues")



pallete1 <- c("#BD0026", "#FC4E2A", "#FEB24C", "#08306B", "#2171B5", "#6BAED6")

##### single point for each line / chro combo ####
# this collapses sibling information from 3 blocks into 1 single number
## plotting across blocks for genral pattern
## but does not show tradeoff in male and female fitness
# due to just 1 example from each line will not be possible to do statistical analysis
combeggcomp <- read.csv("eggs_comp_comparison_combine_blocks.csv")
str(eggcomp)

# make plot the same way as before
combine_P1 <- ggplot(aes(x=comp, y=eggs, color=line), data=combeggcomp)+
  geom_point(aes(size=3))+
  geom_errorbar(aes(ymin=eggsmin, ymax=eggsmax))+
  geom_errorbarh(aes(xmin=compmin, xmax=compmax))+
  ylab("Number of eggs laid")+
  xlab("Proportion of matings by focal males")+
  theme_bw()+
  theme(legend.title=element_blank(),axis.text = element_text(size=25),
        axis.title = element_text(size=30),strip.text=element_text(size=25))+
  scale_colour_manual(values=pallete1)+
  facet_wrap(~chro,scales="free_x")
combine_P1
ggsave("Plots/single_point_per_line.png")
# differnt patterns appear when combining plots
## not much pattern in the B plot but in the D one all S lines high mating and high fecindity
## this corrleates with the findings from the two assays which showwed that the S lines were fitter


## when doing these experiments we used sibling for the male and female assays 
## we there fore have 3 measures for each line covering both assays
## we know tgere is no antagonistic relationship due to the fru alleles but
## there could be an intra generational tradeoff for male and female fitness

eggcomp <- read.csv("eggs_comp_comparison.csv")
str(eggcomp)
## 36 rows of data - 3 for each line / chromsome combination

## change block to a factor
eggcomp$block <- factor(eggcomp$block)
str(eggcomp)

#### plotting with 36 rows ##### 

# we want to have a compariosn of male competion and egg putput for the lines 
# the intial data set has 36 rows, representing the 6 M lines, the 2 chromosome compliments and the 3 blocks of the experiment

## first plot will add a point for each row

P1 <- ggplot(aes(x=comp, y=eggs), data=eggcomp)+
  geom_point(aes())+
  geom_errorbar(aes(ymin=eggsmin, ymax=eggsmax))+
  geom_errorbarh(aes(xmin=compmin, xmax=compmax))
P1
# that works


# now color by line and facet by the chromosomal compliment. also increase size of points 
P2 <- ggplot(aes(x=comp, y=eggs, color=line), data=eggcomp)+
geom_point(aes(size=3))+
  geom_errorbar(aes(ymin=eggsmin, ymax=eggsmax))+
  geom_errorbarh(aes(xmin=compmin, xmax=compmax))+
  facet_wrap(~chro,scales="free_x")
P2
# also looks good

# now add asthetics to make it look nice
P3 <- ggplot(aes(x=comp, y=eggs, color=line), data=eggcomp)+
  geom_point(aes(size=3))+
  geom_errorbar(aes(ymin=eggsmin, ymax=eggsmax))+
  geom_errorbarh(aes(xmin=compmin, xmax=compmax))+
  ylab("Number of eggs laid")+
  xlab("Proporion of matings by focal males")+
  theme_bw()+
  theme(legend.title=element_blank(),axis.text = element_text(size=25),
        axis.title = element_text(size=30),strip.text=element_text(size=25))+
  scale_colour_manual(values=pallete1)+
  facet_wrap(~chro,scales="free_x")
P3
ggsave("Plots/split_chro_no_lines.png")
# looks good
# this current plot has three points for each Mline/chromosome combination for which we have both male and female fitness data = three blocks


## it would be nice to introduce some best fit lines to see trends
#### line plotting ####

### simple comparison of eggs by comp all together
one_line_for_all <- ggplot(aes(y=eggs, x=comp), data=eggcomp)+
  geom_point(aes(size=3))+
  geom_smooth(method=lm, formula = y~x, se=T, fullrange=TRUE, size=2)+
  ylab("Number of eggs laid")+
  xlab("Proportion of matings by focal males")+
  theme_bw()+
  theme(legend.title=element_blank(),axis.text = element_text(size=25),
        axis.title = element_text(size=30),legend.position = "none",strip.text=element_text(size=25))+
  scale_y_continuous(limits=c(0, 60))+
  scale_x_continuous(limits=c(0, 1))+
  theme(legend.position="right")
one_line_for_all
ggsave("Plots/one_line.png")
# should be 36 points
## find gradient
plot.dat <- ggplot_build(one_line_for_all)$data[[2]]
x1 = plot.dat[2,1]
x2 = plot.dat[12,1]
y1 = plot.dat[2,2]
y2 = plot.dat[12,2]
m = (y2-y1)/(x2-x1)
m

### simple first
one_line_split_chro <- ggplot(aes(y=eggs, x=comp), data=eggcomp)+
  geom_point(aes(size=3))+
  geom_smooth(method=lm, formula = y~x, se=T, fullrange=TRUE, size=2)+
  ylab("Number of eggs laid")+
  xlab("Proportion of matings by focal males")+
  theme_bw()+
  theme(legend.title=element_blank(),axis.text = element_text(size=22),
        axis.title = element_text(size=30),legend.position = "none",strip.text=element_text(size=25))+
  scale_y_continuous(limits=c(0, 60))+
  scale_x_continuous(limits=c(0, 1))+
  theme(legend.position="right")+
  facet_wrap(~chro,scales="free_x")
one_line_split_chro
ggsave("Plots/one_line_split_chro.png")
### still negative realtionships for both but they vary, B is more negative 
### insert value for these corellations:



#### best fit per allele ####
### combine lines with each allele togther - 1 line per allele
pallete2 <- c("#FC4E2A", "#2171B5")

one.line.per.allele <- ggplot(aes(y=eggs, x=comp, color=allele), data=eggcomp)+
  geom_point(aes(size=3))+
  geom_smooth(method=lm, formula = y~x, se=F, fullrange=TRUE, size=2)+
  ylab("Number of eggs laid")+
  xlab("Proportion of matings by focal males")+
  theme_bw()+
  theme(legend.title=element_blank(),axis.text = element_text(size=25),
        axis.title = element_text(size=30),legend.position = "none",strip.text=element_text(size=25))+
  scale_y_continuous(limits=c(0, 50))+
  scale_x_continuous(limits=c(0, 1))+
  scale_colour_manual(values=pallete2)+
  theme(legend.position="right")
one.line.per.allele
ggsave("Plots/one_line_per_allele.png")
## 2 negative lines, 18 points per line

plot.dat <- ggplot_build(one.line.per.allele)$data[[2]]
# L gradient first
x1 = plot.dat[2,2]
x2 = plot.dat[12,2]
y1 = plot.dat[2,3]
y2 = plot.dat[12,3]
m = (y2-y1)/(x2-x1)
m
# S gradient
x1 = plot.dat[82,2]
x2 = plot.dat[92,2]
y1 = plot.dat[82,3]
y2 = plot.dat[92,3]
m = (y2-y1)/(x2-x1)
m

## also show this difference by splitting by chro as beofre
one.line.per.allele.split.chro <- ggplot(aes(y=eggs, x=comp, color=allele), data=eggcomp)+
  geom_point(aes(size=3))+
  geom_smooth(method=lm, formula = y~x, se=FALSE, fullrange=TRUE, size=2)+
  ylab("Number of eggs laid")+
  xlab("Proportion of matings by focal males")+
  theme_bw()+
  theme(legend.title=element_blank(),axis.text = element_text(size=22),
        axis.title = element_text(size=30),legend.position = "none",strip.text=element_text(size=25))+
  scale_y_continuous(limits=c(0, 50))+
  scale_x_continuous(limits=c(0, 1))+
  scale_colour_manual(values=pallete2)+
  theme(legend.position="right")+
  facet_wrap(~chro,scales="free_x")
one.line.per.allele.split.chro
ggsave("Plots/one_lineper_allele_split_chro.png")
## caluclate the slopes of these lines


#### best fit per line ####
### eggs by comp for each line only - no splitting by chro
comparison_per_line <- ggplot(aes(y=eggs, x=comp, color=line), data=eggcomp)+
  geom_point(aes(size=3))+
  stat_smooth(method=lm, formula = y~x, se=F, fullrange=T, size=2)+
  ylab("Number of eggs laid")+
  xlab("Proportion of matings by focal males")+
  theme_bw()+
  theme(legend.title=element_blank(),axis.text = element_text(size=25),
        axis.title = element_text(size=30),legend.position = "none",strip.text=element_text(size=25))+
  scale_y_continuous(limits=c(0, 60))+
  scale_x_continuous(limits=c(0, 1))+
  scale_colour_manual(values=pallete1)+
  theme(legend.position="right")
comparison_per_line
ggsave("Plots/comparison_per_line.png")
# 6 negative lines, 6 points per line

### caluclate the slopes of these lines
plot.dat <- ggplot_build(comparison_per_line)$data[[2]]
# L1
m = (plot.dat[12,3]-plot.dat[2,3])/(plot.dat[12,2]-plot.dat[2,2])
m
# L2
m = (plot.dat[92,3]-plot.dat[82,3])/(plot.dat[92,2]-plot.dat[82,2])
m
# L3
m = (plot.dat[172,3]-plot.dat[162,3])/(plot.dat[172,2]-plot.dat[162,2])
m
# S1
m = (plot.dat[252,3]-plot.dat[242,3])/(plot.dat[252,2]-plot.dat[242,2])
m
# S2
m = (plot.dat[332,3]-plot.dat[322,3])/(plot.dat[332,2]-plot.dat[322,2])
m
# S3
m = (plot.dat[412,3]-plot.dat[402,3])/(plot.dat[412,2]-plot.dat[402,2])
m


### option of shaping each point by block 
one.line.per.line.blockshape <- ggplot(aes(y=eggs, x=comp, color=line), data=eggcomp)+
  geom_point(aes(size=3, shape=block))+
  geom_smooth(method=lm, formula = y~x, se=FALSE, fullrange=TRUE, size=2)+
  ylab("Number of eggs laid")+
  xlab("Proportion of matings by focal males")+
  theme_bw()+
  theme(legend.title=element_blank(),axis.text = element_text(size=25),
        axis.title = element_text(size=30),legend.position = "none",strip.text=element_text(size=25))+
  scale_y_continuous(limits=c(0, 60))+
  scale_x_continuous(limits=c(0, 1))+
  scale_colour_manual(values=pallete1)+
  theme(legend.position="right")
one.line.per.line.blockshape
ggsave("Plots/comparison_per_line_blockshape.png")
# with lines, works as only 2 points per line, no so useful


### best fit per line and split by chro

### colour plots by line as in the first plot
one.line.per.line.split.chro <- ggplot(aes(y=eggs, x=comp, color=line), data=eggcomp)+
  geom_point(aes(size=3))+
  geom_smooth(method=lm, formula = y~x, se=FALSE, fullrange=T, size=1)+
  ylab("Number of eggs laid")+
  xlab("Proporion of matings by focal males")+
  theme_bw()+
  theme(legend.title=element_blank(),axis.text = element_text(size=25),
        axis.title = element_text(size=30),legend.position = "none",strip.text=element_text(size=25))+
  scale_y_continuous(limits=c(0, 60))+
  scale_x_continuous(limits=c(0, 1))+
  theme(legend.position="right")+
  scale_colour_manual(values=pallete1)+
  facet_wrap(~chro,scales="free_x")
one.line.per.line.split.chro
ggsave("Plots/per_line_split_chro.png")
## all but S3 D is still negative - probelms with lines here because of only 3 points per line (12 lines in total across the plots)
## bizzarely it seems tat the best fit lines for S2 and S3 are approximatly identical and over lap so much that S2 is invisible 




#### Models to tests for real differences ####

## start simple
hist(eggcomp$eggs)
hist(log(eggcomp$eggs))
## logging doesn't look much better

m1 <- lm(eggs ~ comp, data=eggcomp)
par(mfrow=c(2,2))
plot(m1)
summary(m1)
anova(m1)

# try logging as we did before
m1.1 <- lm(log(eggs) ~ comp, data=eggcomp)
par(mfrow=c(2,2))
plot(m1.1)
summary(m1.1)
# not much better
# keep with un-logged


m2 <- lm(eggs ~ comp + line + chro + allele, data = eggcomp)
summary(m2)
anova(m2)
## over fitted - remove line as replicate of allele 

m3 <- lm(eggs ~ comp + chro + allele, data = eggcomp)
par(mfrow=c(2,2))
plot(m3)
summary(m3)
anova(m3)
Anova(m3)


m4 <- lmer(eggs ~ comp + chro + allele + (1 | block), data = eggcomp)
summary(m4)
anova(m4)
Anova(m4)
## well nothing is significant now

## what about line instead?
m5 <- lmer(eggs ~ comp + chro + allele + (1 | line), data = eggcomp)
summary(m5)
anova(m5)
Anova(m5)

#### model of interest
## links to one line per allele.png
m6 <- lm(eggs ~ comp*allele, data=eggcomp)
summary(m6)
anova(m6)

### link to "comparison_per_line.png"
m7 <- lmer(eggs ~ comp*allele + (1 |line), data=eggcomp)
summary(m7)
anova(m7)
Anova(m7)

## para bootsrapping for p-values
comparison1 <- lmer(eggs ~ 1 + (1 | line), data=eggcomp)
comparison2 <- lmer(eggs ~ allele + (1 | line), data=eggcomp)
comparison3 <- lmer(eggs ~ comp + (1 | line), data=eggcomp)
comparison4 <- lmer(eggs ~ comp + allele + (1 | line), data=eggcomp)
comparison5 <- lmer(eggs ~ comp*allele + (1 | line), data=eggcomp)

comparison.allele <- PBmodcomp(comparison4, comparison3)
comparison.allele
comparison.comp <- PBmodcomp(comparison4, comparison2)
comparison.comp
comparison.inter <- PBmodcomp(comparison5, comparison4)
comparison.inter
