### problems with the IHA ageing dataset extraction
## work done for initially was severely flawed 
## fixed dataset manually
## most of the same code will work

#### create data set ####
lines <- c(rep('L1',65), rep('L2',95), rep('L3',90), rep('S1',110), rep('S2',50), rep('S3',40), rep('L1',65), rep('L2',95), rep('L3',90), rep('S1',110), rep('S2',50), rep('S3',40))
balancertype <- c(rep('B',450), rep('D', 450))
replicares <- c(rep(c(1:13), each=5), rep(c(1:19), each=5), rep(c(1:18), each=5), rep(c(1:22), each=5), rep(c(1:10), each=5), rep(c(1:8), each =5), rep(c(1:13), each=5), rep(c(1:19), each=5), rep(c(1:18), each=5), rep(c(1:22), each=5), rep(c(1:10), each=5), rep(c(1:8), each =5))

# create data frame
data <- data.frame(lines = lines, balancer = balancertype, replicate = replicares)
str(data)


line <- c(rep('L1',88), rep('L2',86), rep('L3',83), rep('L1',82), rep('L2',62), rep('L3',83), rep('S1',94), rep('S2',98), rep('S3',88), 
          rep('S1',83), rep('S2',75), rep('S3',65), rep('L1',68), rep('L2',89), rep('L3',90), rep('L1',86), rep('L2',57), rep('L3',73), 
          rep('S1',81), rep('S2',88), rep('S3',76), rep('S1',70), rep('S2',71), rep('S3',80))
sex <- c(rep('F',987), rep('M',929))
chro <- c(rep('B',257), rep('D',227), rep('B',280), rep('D',223), rep('B',247), rep('D',216), rep('B',245), rep('D',221))
censor <- c(rep(0,83), rep(1,5), rep(0,78), rep(1,8), rep(0,77), rep(1,6), rep(0,76), rep(1,6), rep(0,56), rep(1,6), 
            rep(0,76), rep(1,7), rep(0,82), rep(1,12), rep(0,93), rep(1,5), rep(0,82), rep(1,6), rep(0,74), rep(1,9),
            rep(0,68), rep(1,7), rep(0,54), rep(1,11), rep(0,54), rep(1,14), rep(0,72), rep(1,17), rep(0,70), rep(1,20), 
            rep(0,76), rep(1,10), rep(0,43), rep(1,14), rep(0,61), rep(1,12), rep(0,68), rep(1,13), rep(0,74), rep(1,14), 
            rep(0,66), rep(1,10), rep(0,62), rep(1,8), rep(0,54), rep(1,17), rep(0,60), rep(1,20))
age.prep <- data.frame(line = line, sex = sex, chro = chro, censor = censor)
str(age.prep)

setwd("/Users/michaeljardine/Desktop/DTP/Data and analysis/fruitless/ageing")
write.csv(age.prep, "ageing.manual.csv")

#### read manually prepeared dta and check ####
setwd("/Users/michaeljardine/Desktop/DTP/Data and analysis/fruitless/ageing")
# packages needed
library(survival)
library(plyr)
library(survminer)
library(RColorBrewer)
library(ggsci)

agefru <- read.csv("ageing.manual.csv")
summary(agefru)
str(agefru)
agefru$allele.chro <- interaction(agefru$allele, agefru$chro)
agefru$allele.sex <- interaction(agefru$allele, agefru$sex)
agefru$chro.sex <- interaction(agefru$chro, agefru$sex)
str(agefru)


#### big model with no splitting ####
surv.total.1 <- coxph(Surv(day, censor) ~ chro*allele*sex, data=agefru)
surv.total.1
cox.zph(surv.total.1)
# each variable ok but global < 0.045
ggcoxzph(cox.zph(surv.total.1)) 
summary(surv.total.1)
anova(surv.total.1)
## lots going on as before

## remove 3-way ineraction
surv.total.2 <- coxph(Surv(day, censor) ~ allele*chro + allele*sex + chro*sex, data=agefru)
surv.total.2
cox.zph(surv.total.2)
# < 0.05 for chro by allle and global
ggcoxzph(cox.zph(surv.total.2)) 
summary(surv.total.2)
anova(surv.total.2)

## as before the effect fo the chro comp seesm to confuse and outeay the results
## split by chro and see 
chro_split <- split(agefru, agefru$chro)
agefruB <- chro_split[[1]]
agefruD <- chro_split[[2]]

#### D flies ####
summary(agefruD)
str(agefruD)

D.surv.1 <- coxph(Surv(day, censor) ~ allele*sex, data=agefruD)
D.surv.1
cox.zph(D.surv.1)
# each variable ok but global not
ggcoxzph(cox.zph(D.surv.1))
summary(D.surv.1)
anova(D.surv.1)
 
#### B flies ####
B.surv.1 <- coxph(Surv(day, censor) ~ allele*sex, data=agefruB)
B.surv.1
cox.zph(B.surv.1)
# sex and global not proportional
ggcoxzph(cox.zph(B.surv.1))
summary(B.surv.1)
anova(B.surv.1)
## everything seesm impotant now


#### also split by sex
sex_split <- split(agefru, agefru$sex)
agefruF <- sex_split[[1]]
agefruM <- sex_split[[2]]

### Male flies ####
M.surv.1 <-(coxph(Surv(day, censor) ~ chro*allele, data=agefruM))
M.surv.1
cox.zph(M.surv.1)
## interaction borderline - global fine
ggcoxzph(cox.zph(M.surv.1))
summary(M.surv.1)
anova(M.surv.1)


#### Female flies ####
F.surv.1 <- coxph(Surv(day, censor) ~ chro*allele , data=agefruF)
F.surv.1
cox.zph(F.surv.1)
## all seems fine
ggcoxzph(cox.zph(F.surv.1))
summary (F.surv.1)
anova(F.surv.1)


#####################################
#### PLOTS #####
#####################################
#####################################

####  all flies plots #####
## PLOTS
setwd("/Users/michaeljardine/Desktop/DTP/Data and analysis/fruitless/ageing/")

## different types of these: 
# 1) Kaplan-meier / survival curves
# 2) hazard ratio plots

# 1) Kaplan-Meier survival curves

## work by creating a survival object for each factor and then plotting it
## Sex
fit.total.sex <- survfit(Surv(day, censor) ~ sex, data=agefru)
S.Curve.total.Sex<- ggsurvplot(fit.total.sex, data=agefru, pval=F, palette = c('#e69f00', '#009e73'),
                               xlab="Time in days", censor.size=6, size=1.25, font.x=c(40), font.y=c(40),
                               font.tickslab=c(30), legend.labs=c("Female", "Male"), legend.title="")
S.Curve.total.Sex$plot <-S.Curve.total.Sex$plot + theme(legend.text = element_text(size=25))
S.Curve.total.Sex
ggsave("Plots/Survival_curve_total_sex.png")

##allele
fit.total.allele <- survfit(Surv(day, censor) ~ allele, data=agefru)
S.curv.total.allele <- ggsurvplot(fit.total.allele, data=agefru, pval=F, palette = c('#ca3542', '#7a8ed9'),
                                  xlab="Time in days", censor.size=6, size=1.25, font.x=c(40), font.y=c(40),
                                  font.tickslab=c(30), legend.labs=c("L", "S"), legend.title="")
S.curv.total.allele$plot <-S.curv.total.allele$plot + theme(legend.text = element_text(size=25))
S.curv.total.allele
ggsave("Plots/Survival_curve_total_allele.png")

##chro
fit.total.chro <- survfit(Surv(day, censor) ~ chro, data=agefru)
S.curv.total.chro <- ggsurvplot(fit.total.chro, data=agefru, pval=F, palette = c('#ffd24d', '#b06bc7'), 
                                xlab="Time in days", censor.size=6, size=1.25, font.x=c(40), font.y=c(40),
                                font.tickslab=c(30), legend.labs=c("Balancer", "Deletion"), legend.title="")
S.curv.total.chro$plot <-S.curv.total.chro$plot + theme(legend.text = element_text(size=25))
S.curv.total.chro
ggsave("Plots/Survival_curve_total_chro.png")

## allele x chro
fit.total.A.C <- survfit(Surv(day, censor) ~ allele.chro, data=agefru)
S.curv.total.AxC <- ggsurvplot(fit.total.A.C, data=agefru, pval=F, palette = c('jama'),
                               xlab="Time in days", censor.size=6, size=1.25, font.x=c(40), font.y=c(40),
                               font.tickslab=c(30), legend.labs=c("L / balancer", "L / deletion", "S / balancer", 'S / deletion'), legend.title="")  
S.curv.total.AxC$plot <-S.curv.total.AxC$plot + theme(legend.text = element_text(size=25))
S.curv.total.AxC
ggsave("Plots/Survival_curve_total_allelexchro.png")

## allele x sex
fit.total.A.S<- survfit(Surv(day, censor) ~ allele.sex, data=agefru)
S.curv.total.AxS <- ggsurvplot(fit.total.A.S, data=agefru, pval=F, palette = c('jama'),
                               xlab="Time in days", censor.size=6, size=1.25, font.x=c(40), font.y=c(40),
                               font.tickslab=c(30), legend.labs=c("L / female", "L / male", "S / female", "S / male"), legend.title="")
S.curv.total.AxS$plot <- S.curv.total.AxS$plot +  theme(legend.text = element_text(size=25))
S.curv.total.AxS
ggsave("Plots/Survival_curve_total_allelexsex.png")

## chro x sex
fit.total.C.S <- survfit(Surv(day, censor) ~ chro.sex, data=agefru)
S.curv.total.CxS <- ggsurvplot(fit.total.C.S, data=agefru, pval=F, palette = c('jama'),
                               xlab="Time in days", censor.size=6, size=1.25, font.x=c(40), font.y=c(40),
                               font.tickslab=c(30), legend.labs=c("Female balancer", "Male balancer", "Female deletion", "Male deletion"), legend.title="")
S.curv.total.CxS$plot <- S.curv.total.CxS$plot +  theme(legend.text = element_text(size=25))
S.curv.total.CxS
ggsave("Plots/Survival_curve_total_chroxsex.png")

#####################################
####  D flies plots #####
# sex
fit.D.sex <- survfit(Surv(day, censor) ~ sex, data=agefruD)
S.Curve.D.Sex<- ggsurvplot(fit.D.sex, data=agefruD, pval=F, palette = c('#e69f00', '#009e73'),
                           xlab="Time in days", censor.size=6, size=1.25, font.x=c(40), font.y=c(40),
                           font.tickslab=c(30), legend.labs=c("Female", "Male"), legend.title="")
S.Curve.D.Sex$plot <- S.Curve.D.Sex$plot + theme(legend.text = element_text(size=25))
S.Curve.D.Sex
ggsave("Plots/Survival_curve_D_sex.png")

# allele
fit.D.allele <- survfit(Surv(day, censor) ~ allele, data=agefruD)
S.curv.D.allele <- ggsurvplot(fit.D.allele, data=agefruD, pval=F, palette = c('#ca3542', '#7a8ed9'),
                              xlab="Time in days", censor.size=6, size=1.25, font.x=c(40), font.y=c(40),
                              font.tickslab=c(30), legend.labs=c("L", "S"), legend.title="")
S.curv.D.allele$plot <- S.curv.D.allele$plot + theme(legend.text = element_text(size=25))
S.curv.D.allele
ggsave("Plots/Survival_curve_D_allele.png")

# allele x sex
fit.D.A.S<- survfit(Surv(day, censor) ~ allele.sex, data=agefruD)
S.curv.D.AxS <- ggsurvplot(fit.D.A.S, data=agefruD, pval=F, palette = 'jama',
                           xlab="Time in days", censor.size=6, size=1.25, font.x=c(40), font.y=c(40),
                           font.tickslab=c(30), legend.labs=c("L / female", "S / female", "L / male", "S / male"), legend.title="")
S.curv.D.AxS$plot <- S.curv.D.AxS$plot + theme(legend.text = element_text(size=25))
S.curv.D.AxS
ggsave("Plots/Survival_curve_D_allelexsex.png")

## redo this plot for paper with consistant colours from other plots
S.curv.D.AxS.paper <- ggsurvplot(fit.D.A.S, data=agefruD, pval=F, palette = c('#CA3542', '#aecbc9', '#CA3542', '#aecbc9'),
                            linetype = c(1, 1, 4, 4), xlab="Time in days", censor.size=6, size=1.25, font.x=c(40), font.y=c(40),
                           font.tickslab=c(30), legend.labs=c("L / female", "S / female", "L / male", "S / male"), legend.title="")
S.curv.D.AxS.paper$plot <- S.curv.D.AxS.paper$plot + theme(legend.text = element_text(size=40))
S.curv.D.AxS.paper
ggsave("Plots/Survival_curve_D_allelexsex_paper.png")

## original colours for these plots are difficult to see so here is an alternative
### USE THIS ONE!!!
S.curv.D.AxS.paper2 <- ggsurvplot(fit.D.A.S, data=agefruD, pval=F, palette = c('#CA3542', '#7a8ed9', '#CA3542', '#7a8ed9'),
                                 linetype = c(1, 1, 6, 6), xlab="Time in days", censor.size=6, size=1.8, font.x=c(40), font.y=c(40),
                                 font.tickslab=c(30), legend.title="")
S.curv.D.AxS.paper2$plot <- S.curv.D.AxS.paper2$plot + 
  theme(axis.title.x = element_text(size=50), axis.title.y = element_text(size=45))+
  theme(axis.title.y = element_text(vjus =1.4))
S.curv.D.AxS.paper2
ggsave("Plots/Survival_curve_D_allelexsex_paper2.png")

### dotted line
S.curv.D.AxS.paper.dots <- ggsurvplot(fit.D.A.S, data=agefruD, pval=F, palette = c('#CA3542', '#7a8ed9', '#CA3542', '#7a8ed9'),
                                  linetype = c(1, 1, 3, 3), xlab="Time in days", censor.size=6, size=1.8, font.x=c(40), font.y=c(40),
                                  font.tickslab=c(30), legend.labs=c("L/female", "S/female", "L/male", "S/male"), legend.title="")
S.curv.D.AxS.paper.dots$plot <- S.curv.D.AxS.paper.dots$plot + 
  theme(legend.text = element_text(size=55), axis.title.x = element_text(size=50), axis.title.y = element_text(size=45))+
  theme(axis.title.y = element_text(vjus =1.4))
S.curv.D.AxS.paper.dots
ggsave("Plots/Survival_curve_D_allelexsex_paper_dots.png")


### 4 colour e
S.curv.D.AxS.paper.4cols <- ggsurvplot(fit.D.A.S, data=agefruD, pval=F, palette = c('#d55d68', '#a2b0e4', '#a22a35', '#526cce'),
                                      linetype = c(1, 1, 1, 1), xlab="Time in days", censor.size=6, size=1.8, font.x=c(40), font.y=c(40),
                                      font.tickslab=c(30), legend.labs=c("L/female", "S/female", "L/male", "S/male"), legend.title="")
S.curv.D.AxS.paper.4cols$plot <- S.curv.D.AxS.paper.4cols$plot + 
  theme(legend.text = element_text(size=55), axis.title.x = element_text(size=50), axis.title.y = element_text(size=45))+
  theme(axis.title.y = element_text(vjus =1.4))
S.curv.D.AxS.paper.4cols
ggsave("Plots/Survival_curve_D_allelexsex_paper_4cols.png")

#####################################
####  B flies plots #####
# sex
fit.B.sex <- survfit(Surv(day, censor) ~ sex, data=agefruB)
S.Curve.B.Sex<- ggsurvplot(fit.B.sex, data=agefruB, pval=F, palette = c('#e69f00', '#009e73'),
                           xlab="Time in days", censor.size=6, size=1.25, font.x=c(40), font.y=c(40),
                           font.tickslab=c(30), xlim=c(0,100), break.x.by=25, legend.labs=c("Female", "Male"), legend.title="")
S.Curve.B.Sex$plot <- S.Curve.B.Sex$plot + theme(legend.text = element_text(size=25))
S.Curve.B.Sex
ggsave("Plots/Survival_curve_B_sex.png")

# allele
fit.B.allele <- survfit(Surv(day, censor) ~ allele, data=agefruB)
S.curv.B.allele <- ggsurvplot(fit.B.allele, data=agefruB, pval=F, palette = c('#ca3542', '#7a8ed9'),
                              xlab="Time in days", censor.size=6, size=1.25, font.x=c(40), font.y=c(40),
                              font.tickslab=c(30), xlim=c(0,100), break.x.by=25, legend.labs=c("L", "S"), legend.title="")
S.curv.B.allele$plot <- S.curv.B.allele$plot + theme(legend.text = element_text(size=25))
S.curv.B.allele
ggsave("Plots/Survival_curve_B_allele.png")

# allele x sex
fit.B.A.S<- survfit(Surv(day, censor) ~ allele.sex, data=agefruB)
S.curv.B.AxS <- ggsurvplot(fit.B.A.S, data=agefruB, pval=F, palette = 'jama',
                           xlab="Time in days", censor.size=6, size=1.25, font.x=c(40), font.y=c(40),
                           font.tickslab=c(30), xlim=c(0,100), break.x.by=25, legend.labs=c("L / female", "L / male", "S / female", "S / male"), legend.title="")
S.curv.B.AxS$plot <- S.curv.B.AxS$plot + theme(legend.text = element_text(size=25))
S.curv.B.AxS
ggsave("Plots/Survival_curve_B_allelexsex.png")

## redo this plot for paper with consistant colours from other plots
S.curv.B.AxS.paper <- ggsurvplot(fit.B.A.S, data=agefruD, pval=F, palette = c('#CA3542', '#aecbc9', '#CA3542', '#aecbc9'),
                                 linetype = c(1, 1, 6, 6), xlab="Time in days", censor.size=6, size=1.25, font.x=c(40), font.y=c(40),
                                 font.tickslab=c(30), xlim=c(0,100), break.x.by=25, legend.labs=c("L / female", "S / female", "L / male", "S / male"), legend.title="")
S.curv.B.AxS.paper$plot <- S.curv.B.AxS.paper$plot + theme(legend.text = element_text(size=40))
S.curv.B.AxS.paper
ggsave("Plots/Survival_curve_B_allelexsex_paper.png")

## original colours for these plots are difficult to see so here is an alternative
S.curv.B.AxS.paper2 <- ggsurvplot(fit.B.A.S, data=agefruD, pval=F, palette = c('#CA3542', '#7a8ed9', '#CA3542', '#7a8ed9'),
                                 linetype = c(1, 1, 6, 6), xlab="Time in days", censor.size=8, size=1.8, font.x=c(40), font.y=c(40),
                                 font.tickslab=c(30), xlim=c(0,100),  break.x.by=25, legend.title="")
S.curv.B.AxS.paper2$plot <- S.curv.B.AxS.paper2$plot + 
  theme(axis.title.x = element_text(size=50), axis.title.y = element_text(size=45))+
  theme(axis.title.y = element_text(vjus =1.4))
S.curv.B.AxS.paper2
ggsave("Plots/Survival_curve_B_allelexsex_paper2.png")

## dotted
S.curv.B.AxS.paper.dots <- ggsurvplot(fit.B.A.S, data=agefruD, pval=F, palette = c('#CA3542', '#7a8ed9', '#CA3542', '#7a8ed9'),
                                  linetype = c(1, 1, 3, 3), xlab="Time in days", censor.size=6, size=1.8, font.x=c(40), font.y=c(40),
                                  font.tickslab=c(30), xlim=c(0,100),  break.x.by=25, legend.labs=c("L/female", "S/female", "L/male", "S/male"), legend.title="")
S.curv.B.AxS.paper.dots$plot <- S.curv.B.AxS.paper.dots$plot + 
  theme(legend.text = element_text(size=55), axis.title.x = element_text(size=50), axis.title.y = element_text(size=45))+
  theme(axis.title.y = element_text(vjus =1.4))
S.curv.B.AxS.paper.dots
ggsave("Plots/Survival_curve_B_allelexsex_paper_dots.png")

## 4 colours
S.curv.B.AxS.paper.4cols <- ggsurvplot(fit.B.A.S, data=agefruD, pval=F, palette = c('#d55d68', '#a2b0e4', '#a22a35', '#526cce'),
                                      linetype = c(1, 1, 1, 1), xlab="Time in days", censor.size=6, size=1.8, font.x=c(40), font.y=c(40),
                                      font.tickslab=c(30), xlim=c(0,100),  break.x.by=25, legend.labs=c("L/female", "S/female", "L/male", "S/male"), legend.title="")
S.curv.B.AxS.paper.4cols$plot <- S.curv.B.AxS.paper.4cols$plot + 
  theme(legend.text = element_text(size=55), axis.title.x = element_text(size=50), axis.title.y = element_text(size=45))+
  theme(axis.title.y = element_text(vjus =1.4))
S.curv.B.AxS.paper.4cols
ggsave("Plots/Survival_curve_B_allelexsex_paper_4cols.png")

#####################################
####  Male flies plots #####

##allele
fit.M.allele <- survfit(Surv(day, censor) ~ allele, data=agefruM)
S.curv.M.allele <- ggsurvplot(fit.M.allele, data=agefruM, pval=F, palette = c('#ca3542', '#7a8ed9'),
                              xlab="Time in days", censor.size=6, size=1.25, font.x=c(40), font.y=c(40),
                             font.tickslab=c(30), xlim=c(0,100), break.x.by=25, legend.labs=c("L", "S"), legend.title="")
S.curv.M.allele$plot <- S.curv.M.allele$plot + theme(legend.text = element_text(size=25))
S.curv.M.allele
ggsave("Plots/Survival_curve_M_allele_paper2.png")

## chro
fit.M.chro <- survfit(Surv(day, censor) ~ chro, data=agefruM)
S.curv.M.chro <- ggsurvplot(fit.M.chro, data=agefruM, pval=F, palette = c('#ffd24d', '#b06bc7'), 
                            xlab="Time in days", censor.size=6, size=1.25, font.x=c(40), font.y=c(40),
                            font.tickslab=c(30), xlim=c(0,100), break.x.by=25, legend.labs=c("Balancer", "Deletion"), legend.title="")
S.curv.M.chro$plot <- S.curv.M.chro$plot + theme(legend.text = element_text(size=25))
S.curv.M.chro
ggsave("Plots/Survival_curve_M_chro.png")

## allele chro
fit.M.A.C <- survfit(Surv(day, censor) ~ allele.chro, data=agefruM)
S.curv.M.AxC <- ggsurvplot(fit.M.A.C, data=agefruM, pval=F, palette = 'jama',
                           xlab="Time in days", censor.size=6, size=1.25, font.x=c(40), font.y=c(40),
                           font.tickslab=c(30), xlim=c(0,100), break.x.by=25, legend.labs=c("L / balancer", "L / deletion", "S / balancer", 'S / deletion'), legend.title="")
S.curv.M.AxC$plot <- S.curv.M.AxC$plot + theme(legend.text = element_text(size=25))
S.curv.M.AxC
ggsave("Plots/Survival_curve_M_allelexchro.png")


#####################################
####  Female flies plots #####

##allele
fit.F.allele <- survfit(Surv(day, censor) ~ allele, data=agefruF)
S.curv.F.allele <- ggsurvplot(fit.F.allele, data=agefruF, pval=F, palette = c('#ca3542', '#7a8ed9'),
                              xlab="Time in days", censor.size=6, size=1.25, font.x=c(40), font.y=c(40),
                              font.tickslab=c(30), xlim=c(0,100), break.x.by=25, legend.labs=c("L", "S"), legend.title="")
S.curv.F.allele$plot <- S.curv.F.allele$plot + theme(legend.text = element_text(size=25))
S.curv.F.allele
ggsave("Plots/Survival_curve_F_allele.png")

## chro
fit.F.chro <- survfit(Surv(day, censor) ~ chro, data=agefruF)
S.curv.F.chro <- ggsurvplot(fit.F.chro, data=agefruF, pval=F, palette = c('#ffd24d', '#b06bc7'), 
                            xlab="Time in days", censor.size=6, size=1.25, font.x=c(40), font.y=c(40),
                            font.tickslab=c(30), xlim=c(0,100), break.x.by=25, legend.labs=c("Balancer", "Deletion"), legend.title="")
S.curv.F.chro$plot <- S.curv.F.chro$plot + theme(legend.text = element_text(size=25))
S.curv.F.chro
ggsave("Plots/Survival_curve_F_chro.png")

## allele chro
fit.F.A.C <- survfit(Surv(day, censor) ~ allele.chro, data=agefruF)
S.curv.F.AxC <- ggsurvplot(fit.F.A.C, data=agefruF, pval=F, palette = 'jama',
                           xlab="Time in days", censor.size=6, size=1.25, font.x=c(40), font.y=c(40),
                           font.tickslab=c(30), xlim=c(0,100), break.x.by=25, legend.labs=c("L / balancer", "L / deletion", "S / balancer", 'S / deletion'), legend.title="")
S.curv.F.AxC$plot <- S.curv.F.AxC$plot + theme(legend.text = element_text(size=25))
S.curv.F.AxC
ggsave("Plots/Survival_curve_F_allelexchro.png")

#####################################
####  HR plots ####

### 1 plot per model set

## all flies
term <- c("alleleS", "chroD", "SexM", 'alleleS:chroD', "alleleS:SexM", "chroD:SexM")
h.ratio <- c(1.3183, 0.5156, 0.5305, 0.6925, 0.8209, 2.6241)
min.D <- c(1.1259, 0.4395, 0.4488, 0.5701, 0.6762, 2.1535)
max.D <- c(1.5437, 0.6119, 0.627, 0.8412, 0.9966, 3.1976)
T.HR <- data.frame(Term = term, HR = h.ratio, HR.min = min.D, HR.max = max.D)
T.HR
levels(T.HR$Term)
T.HR$Term <- factor(T.HR$Term, levels = c("alleleS", "chroD", "SexM", 'alleleS:chroD', "alleleS:SexM", "chroD:SexM"))

HR.plot.T<- ggplot(aes(x=Term, y=HR, color=Term), data=T.HR)+
  geom_point(aes(size=10))+
  geom_errorbar(aes(ymin=(min.D), ymax=(max.D)), size=2)+
  geom_hline(yintercept = 1)+
  ylab("Hazard-ratio")+
  xlab("Term in Cox model")+
  theme_bw()+
  theme(legend.title=element_blank(),axis.text = element_text(size=25),
        axis.title = element_text(size=30),legend.position = "none",strip.text=element_text(size=20))+
  scale_color_jama()
HR.plot.T
ggsave("Plots/HZ/T_HR_plot.png")

# D flies
term <- c("alleleS", "sexM", 'alleleS:sexM')
h.ratio <- c(0.8702, 1.3204, 0.927)
min.D <- c(0.7153, 1.0805, 0.6963)
max.D <- c(1.059, 1.614, 1.234)
D.HR <- data.frame(Term = term, HR = h.ratio, HR.min = min.D, HR.max = max.D)
D.HR
levels(D.HR$Term)
D.HR$Term <- factor(D.HR$Term, levels = c("alleleS", "sexM", "alleleS:sexM"))

HR.plot.D<- ggplot(aes(x=Term, y=HR, color=Term), data=D.HR)+
  geom_point(aes(size=10))+
  geom_errorbar(aes(ymin=(min.D), ymax=(max.D)), size=2)+
  geom_hline(yintercept = 1)+
  ylab("Hazard-ratio")+
  xlab("Term in Cox model")+
  theme_bw()+
  theme(legend.title=element_blank(),axis.text = element_text(size=25),
        axis.title = element_text(size=30),legend.position = "none",strip.text=element_text(size=25))+
  scale_color_jama()
HR.plot.D
ggsave("Plots/HZ/D_HR_plot.png")

# B flies
term <- c("alleleS", "sexM", 'alleleS:sexM')
h.ratio <- c(1.3857, 0.5719, 0.7308)
min.D <- c(1.16, 0.4724, 0.5605)
max.D <- c(1.6554, 0.6924, 0.9527)
B.HR <- data.frame(Term = term, HR = h.ratio, HR.min = min.D, HR.max = max.D)
B.HR
levels(B.HR$Term)
B.HR$Term <- factor(B.HR$Term, levels = c("alleleS", "sexM", "alleleS:sexM"))

HR.plot.B <- ggplot(aes(x=Term, y=HR, color=Term), data=B.HR)+
  geom_point(aes(size=10))+
  geom_errorbar(aes(ymin=(min.D), ymax=(max.D)), size=2)+
  geom_hline(yintercept = 1)+
  ylab("Hazard-ratio")+
  xlab("Term in Cox model")+
  theme_bw()+
  theme(legend.title=element_blank(),axis.text = element_text(size=25),
        axis.title = element_text(size=30),legend.position = "none",strip.text=element_text(size=25))+
  scale_color_jama()
HR.plot.B
ggsave("Plots/HZ/B_HR_plot.png")

# male flies
term <- c("alleleS", "chromosomeD", 'alleleS:chromosomeD')
h.ratio <- c(1.0385, 1.3008, 0.7724)
min.D <- c(0.8539, 1.0611, 0.5798)
max.D <- c(1.263, 1.595, 1.029)
M.HR <- data.frame(Term = term, HR = h.ratio, HR.min = min.D, HR.max = max.D)
M.HR
levels(M.HR$Term)
M.HR$Term <- factor(F.HR$Term, levels = c("alleleS", "chromosomeD", "alleleS:chromosomeD"))

HR.plot.M<- ggplot(aes(x=Term, y=HR, color=Term), data=M.HR)+
  geom_point(aes(size=10))+
  geom_errorbar(aes(ymin=(min.D), ymax=(max.D)), size=2)+
  geom_hline(yintercept = 1)+
  ylab("Hazard-ratio")+
  xlab("Term in Cox model")+
  theme_bw()+
  theme(legend.title=element_blank(),axis.text = element_text(size=25),
        axis.title = element_text(size=30),legend.position = "none",strip.text=element_text(size=25))+
  scale_color_jama()
HR.plot.M
ggsave("Plots/HZ/M_HR_plot.png")

# female flies
term <- c("alleleS", "chromosomeD", 'alleleS:chromosomeD')
h.ratio <- c(1.3814, 0.5423, 0.6114)
min.D <- c(1.1565, 0.4491, 0.4685)
max.D <- c(1.65, 0.6548, 0.798)
F.HR <- data.frame(Term = term, HR = h.ratio, HR.min = min.D, HR.max = max.D)
F.HR
levels(F.HR$Term)
F.HR$Term <- factor(F.HR$Term, levels = c("alleleS", "chromosomeD", "alleleS:chromosomeD"))

HR.plot.F<- ggplot(aes(x=Term, y=HR, color=Term), data=F.HR)+
  geom_point(aes(size=10))+
  geom_errorbar(aes(ymin=(min.D), ymax=(max.D)), size=2)+
  geom_hline(yintercept = 1)+
  ylab("Hazard-ratio")+
  xlab("Term in Cox model")+
  theme_bw()+
  theme(legend.title=element_blank(),axis.text = element_text(size=25),
        axis.title = element_text(size=30),legend.position = "none",strip.text=element_text(size=25))+
  scale_color_jama()
HR.plot.F
ggsave("Plots/HZ/F_HR_plot.png")

### combined plots comparing models

# with data from combined models above
model <- c("D flies", "B flies", "Males", "Females")
Haz_R <- c(0.8702, 1.3857, 1.0385, 1.3414)
minHR <- c(0.7153, 1.16, 0.8539, 1.1565)
maxHR <- c(1.059, 1.6554, 1.2663, 1.65)
C.HR <- data.frame(Model = model, HR = Haz_R, HR.min = minHR, HR.max = maxHR)
levels(C.HR$Model)

HR.comp.1 <- ggplot(aes(x=Model, y=HR, color=Model), data=C.HR)+
  geom_point()+
  geom_errorbar(aes(ymin=(minHR), ymax=(maxHR)), size=2)+
  geom_hline(yintercept = 1)+
  ylab("Hazard-ratio of the S allele")+
  xlab("Subset of data")+
  theme_bw()+
  theme(legend.title=element_blank(),axis.text = element_text(size=25),
        axis.title = element_text(size=30),legend.position = "none",strip.text=element_text(size=25))+
  scale_color_jama()
HR.comp.1
ggsave("Plots/HZ/HR_Sallele_from_models.png")

### procure seprate HR for the S allele for each of th different subsets of data
### can produce plots simialr to waht we had before

chrosex_split <- split(agefruD, agefruD$sex)
agefruDF <- chrosex_split[[1]]
agefruDM <- chrosex_split[[2]]

chroBsex_split <- split(agefruB, agefruB$sex)
agefruBF <- chroBsex_split[[1]]
agefruBM <- chroBsex_split[[2]]

## models
DF.surv <- coxph(Surv(day, censor) ~ allele, data=agefruDF)
DM.surv <- coxph(Surv(day, censor) ~ allele, data=agefruDM)
BF.surv <- coxph(Surv(day, censor) ~ allele, data=agefruBF)
BM.surv <- coxph(Surv(day, censor) ~ allele, data=agefruBM)

# MODEL OUPUT
# D females
summary(DF.surv)
# D males
summary(DM.surv)
# B females
summary(BF.surv)
# B males
summary(BM.surv)

## plot this
term <- c("B Males", "B Females", 'D Males', "D Females")
chro <- c("B", "B", "D", "D")
sex <- c ("M", "F", "M", "F")
h.ratio <- c(1.037, 1.341, 0.8322, 0.8464)
min.D <- c(0.8517, 1.121, 0.6746, 0.6945)
max.D <- c(1.264, 1.603, 1.026, 1.031)
S.HR <- data.frame(Term = term, HR = h.ratio, HR.min = min.D, HR.max = max.D, chromosome = chro, sex = sex)
S.HR
levels(S.HR$Term)

S_HR.plots.single.term.models<- ggplot(aes(x=Term, y=HR, color=sex), data=S.HR)+
  geom_point(size=5)+
  geom_errorbar(aes(ymin=(min.D), ymax=(max.D)), size=2)+
  geom_hline(yintercept = 1)+
  ylab("H-R of S allele")+
  xlab("Sex and chromosome combination")+
  theme_bw()+
  theme(legend.title=element_blank(),axis.text = element_text(size=25),
        axis.title = element_text(size=30),legend.position = "none",strip.text=element_text(size=25))+
  scale_colour_manual(values=c("#E69F00", "#009E73"))+
  facet_wrap(~chromosome,scales="free_x")
S_HR.plots.single.term.models
ggsave("Plots/HZ/simple_model_Sallele_plot.png")




#### max's extra stuff
## only use allele

## all flies
surv.total.allele.only <- coxph(Surv(day, censor) ~ allele, data=agefru)
surv.total.allele.only
summary(surv.total.allele.only)
anova(surv.total.allele.only)

# split them up
chrosex_split <- split(agefruD, agefruD$sex)
agefruDF <- chrosex_split[[1]]
agefruDM <- chrosex_split[[2]]

chroBsex_split <- split(agefruB, agefruB$sex)
agefruBF <- chroBsex_split[[1]]
agefruBM <- chroBsex_split[[2]]

DF.surv <- coxph(Surv(day, censor) ~ allele, data=agefruDF)
DM.surv <- coxph(Surv(day, censor) ~ allele, data=agefruDM)
BF.surv <- coxph(Surv(day, censor) ~ allele, data=agefruBF)
BM.surv <- coxph(Surv(day, censor) ~ allele, data=agefruBM)

# MODEL OUPUT
# D females
summary(DF.surv)
# D males
summary(DM.surv)
# B females
summary(BF.surv)
# B males
summary(BM.surv)




## B females flies
BF.surv
summary(BF.surv)
anova(BF.surv)

# survival~sex in D data only
D.surv.sex <- coxph(Surv(day, censor) ~ sex, data=agefruD)
D.surv.sex
summary(D.surv.sex)
anova(D.surv.sex)

# survival~sex in B data only
B.surv.sex <- coxph(Surv(day, censor) ~ sex, data=agefruB)
B.surv.sex
summary(B.surv.sex)
anova(B.surv.sex)

# survival~allele in D data only
D.surv.allele <- coxph(Surv(day, censor) ~ allele, data=agefruD)
D.surv.allele
summary(D.surv.allele)
anova(D.surv.allele)

# survival~allele in B data only
B.surv.allele <- coxph(Surv(day, censor) ~ allele, data=agefruB)
B.surv.allele
summary(B.surv.allele)
anova(B.surv.allele)

# survival~allele in female data only
Fem.surv.allele <- coxph(Surv(day, censor) ~ allele, data=agefruF)
Fem.surv.allele
summary(Fem.surv.allele)
anova(Fem.surv.allele)

# survival~allele in male data only
Male.surv.allele <- coxph(Surv(day, censor) ~ allele, data=agefruM)
Male.surv.allele
summary(Male.surv.allele)
anova(Male.surv.allele)

## all by sex
surv.total.sex.only <- coxph(Surv(day, censor) ~ sex, data=agefru)
surv.total.sex.only
summary(surv.total.sex.only)
anova(surv.total.sex.only)

## all by complement
surv.total.chro.only <- coxph(Surv(day, censor) ~ chro, data=agefru)
surv.total.chro.only
summary(surv.total.chro.only)
anova(surv.total.chro.only)
