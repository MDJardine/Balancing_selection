library(ggplot2)

setwd("/Users/michaeljardine/Desktop/DTP/Data and analysis/Cage analysis/fruitless/")
fru <- read.delim("fru_rc")

str(fru)
summary(fru)

#### Time point 1 only ####

## use select_SNP_freq as the starting point
# select_SNP_frequencies <- read.csv("SNP_frequencies_near_fru.csv")
# for some reason there are later probelsm if doing this.
# run the whole prep_for_popoolation_data.R first


## get only samples from the first time point
TP1 <- select_SNP_freq[,grepl("X1.", colnames(select_SNP_freq))]
## add in the other information
TP1$MA <- select_SNP_freq$MA
str(TP1)

which(TP==1)


## creare vecotr of the major allele states for time point 1 only
maa <- length(select_SNP_freq$major_alleles)
 for(i in 1:nrow(select_SNP_freq)){
  temp <- unlist(strsplit(as.character(select_SNP_freq$major_alleles)[i], split=""))
  temp2 <- temp[which(TP==1)]
  maa[i] <- paste(temp2, collapse = "")
}
maa


## add this information to frequc data for timepoint 1
TP1$major_alleles <- maa
## add all extra information regarding position count etc
TP1 <- cbind(TP1, full_SNP_info)
### now find the polymorphic rows at time point 1

## Polymorphism loop - which rows are polymorphic?
polymorph_test <- vector(length=nrow(TP1))
for(i in 1:nrow(TP1)){
  tmp_allele_1 <- unlist(strsplit(as.character(TP1$MA)[i], split="/"))[1]
  tmp_major_alleles <- unlist(strsplit(as.character(TP1$major_alleles)[i], split=""))
  polymorph_test[i] <- sum(1*(tmp_major_alleles==tmp_allele_1))>4 & sum(1*(tmp_major_alleles==tmp_allele_1))<7 # can change the > & < numbers to alter the selection criteria
  
 # !(sum(1*(tmp_major_alleles==tmp_allele_1))==10 | sum(1*(tmp_major_alleles==tmp_allele_1))==0)
}
polymorph_test

## trim data set based on the amount of polymorphism
TP1_poly_select <- TP1[polymorph_test,]
str(TP1_poly_select)
summary(TP1_poly_select)
## produces 20 polymorphic rows  

## for some reason two rows have an NAs - poor coverage? - least variable loci anyway
TP1_poly_select <- na.omit(TP1_poly_select)
str(TP1_poly_select)
# now only 9 rows

## doesn't look good
write.csv(TP1_poly_select, "timepoint1.csv")






#### compare the two sets of 5 cages to see which time points show interesting patterns ####

test <- TP1_poly_select[3,]
TP1_poly_select[3, c(1,3,5,7,9)]
str(TP1_poly_select)

TP1_poly_select[3, c(1,3,5,7,9)]


colMeans(TP1_poly_select[, c(1,3,5,7,9)])

high.MB <- rowMeans(TP1_poly_select[, c(1,3,5,7,9)])
high.FB <- rowMeans(TP1_poly_select[, c(2,4,6,8,10)])

mean.freq.TP1 <- data.frame(high.MB, high.FB)

## take information from the former data set
info <- TP1_poly_select[,11:19]
info2 <- rbind(info, info)
mean.freq.TP1 <- cbind(mean.freq.TP1, info)


major.alleles.freqs.tp1 <- c(high.MB, high.FB)

starting.states <- c(rep("MB", 9), rep("FB", 9))

mean.freq.TP1 <- data.frame(major.alleles.freqs.tp1, starting.states)
mean.freq.TP1 <- cbind(mean.freq.TP1, info2)
str(mean.freq.TP1)
mean.freq.TP1$pos <- as.factor(mean.freq.TP1$pos)
str(mean.freq.TP1)
write.csv(mean.freq.TP1, "time_point1_mean_frequncy comaprison.csv")


TP1.comparison <- ggplot(aes(x=pos, y=major.alleles.freqs.tp1, fill=starting.states), data=mean.freq.TP1)+
  geom_bar(stat="identity", position=position_dodge())
TP1.comparison



#### histograms of the poltmorphic sites from tp1 across all time points ####
## do as a total or 
## as 5 and 5 as before? 

## take a vecotr of the positions of the 9 sites

positions.of.tp1 <- TP1_poly_select$pos
# filter larger date set by this
poly.tp1.full <- select_SNP_freq[select_SNP_freq$pos %in% positions.of.tp1, ]


SNP.1181 <- poly.tp1.full[3,]
SNP.1208 <- poly.tp1.full[4,]
SNP.1521 <- poly.tp1.full[5,]

write.csv(SNP.1181, "SNP.1181.csv")
write.csv(SNP.1208, "SNP.1208.csv")
write.csv(SNP.1521, "SNP.1521.csv")

## plan to produce 1 plot for each position
## see the difference in frequencies between the two set of 5 over time

#### plottign tp1 loci of interest ####


## first remove extra information
SNP.1181.info <- SNP.1181[,1:9]
# and save the rest as the frequencies
SNP.1181.freqs <- SNP.1181[,10:99]
# flip oreintation of this dataset
SNP.1181.freqs.t <- data.frame(t(SNP.1181.freqs))
str(SNP.1181.freqs.t)
## create cage identity vector 
Cage <- c(1:9, 1:8, 10, 1:7, 9, 10, 1:6, 8:10, 1:5, 7:10, 1:4, 6:10, 1:3, 5:10, 1, 2, 4:10, 1, 3:10, 2:10)
## add to freqs
SNP.1181.freqs.t$cage <- Cage
## timepoint sample vecotor
collection.timepoint <- c("4", "7", "5", "8", "1", "2", "6", "3", "9", "7", "5", "8", "1", "2", "6", "3", "9", "4", "5", "8", "1", "2", "6", "3", "9", "4", "7", "8", "1", "2", 
                                  "6", "3", "9", "4", "7", "5", "1", "2", "6", "3", "9", "4", "7", "5", "8", "2", "6", "3", "9", "4", "7", "5", "8", "1", "6", "3", "9", "4", "7", "5", 
                                  "8", "1", "2", "3", "9", "4", "7", "5", "8", "1", "2", "6", "9", "4", "7", "5", "8", "1", "2", "6", "3", "4", "7", "5", "8", "1", "2", "6", "3", "9")
## add the freqs
SNP.1181.freqs.t$timepoint <- collection.timepoint
## vectr of moths from start
month <- c("6", "18", "10", "22", "1", "2", "14", "4", "28", "18", "10", "22", "1", "2", "14", "4", "28", "6", "10", "22", "1", "2", "14", "4", "28", "6", "18", "22", "1", "2", 
            "14", "4", "28", "6", "18", "10", "1", "2", "14", "4", "28", "6", "18", "10", "22", "2", "14", "4", "28", "6", "18", "10", "22", "1", "14", "4", "28", "6", "18", "10", 
            "22", "1", "2", "4", "28", "6", "18", "10", "22", "1", "2", "14", "28", "6", "18", "10", "22", "1", "2", "14", "4", "6", "18", "10", "22", "1", "2", "14", "4", "28")
## add to freqs
SNP.1181.freqs.t$month <- month
str(SNP.1181.freqs.t)
SNP.1181.freqs.t$cage.name <- as.factor(SNP.1181.freqs.t$cage)
SNP.1181.freqs.t$month <- as.numeric(SNP.1181.freqs.t$month)
str(SNP.1181.freqs.t)
## name columns correctly
colnames(SNP.1181.freqs.t) <- c("frequency", "cage", "timepoint", "month", "cage_name")

## repeat for the other 2 loci of interest
### .1208
## first remove extra information
SNP.1208.info <- SNP.1208[,1:9]
# and save the rest as the frequencies
SNP.1208.freqs <- SNP.1208[,10:99]
# flip oreintation of this dataset
SNP.1208.freqs.t <- data.frame(t(SNP.1208.freqs))
## add cage info
SNP.1208.freqs.t$cage <- Cage
## add timepoint info
SNP.1208.freqs.t$timepoint <- collection.timepoint
## add month info
SNP.1208.freqs.t$month <- month
str(SNP.1208.freqs.t)
## corrections
SNP.1208.freqs.t$cage <- as.factor(SNP.1208.freqs.t$cage)
SNP.1208.freqs.t$month <- as.numeric(SNP.1208.freqs.t$month)
colnames(SNP.1208.freqs.t) <- c("frequency", "cage", "timepoint", "month")

### .1521
## first remove extra information
SNP.1521.info <- SNP.1521[,1:9]
# and save the rest as the frequencies
SNP.1521.freqs <- SNP.1521[,10:99]
# flip oreintation of this dataset
SNP.1521.freqs.t <- data.frame(t(SNP.1521.freqs))
## add cage info
SNP.1521.freqs.t$cage <- Cage
## add timepoint info
SNP.1521.freqs.t$timepoint <- collection.timepoint
## add month info
SNP.1521.freqs.t$month <- month
str(SNP.1521.freqs.t)
## corrections
SNP.1521.freqs.t$cage <- as.factor(SNP.1521.freqs.t$cage)
SNP.1521.freqs.t$month <- as.numeric(SNP.1521.freqs.t$month)
colnames(SNP.1521.freqs.t) <- c("frequency", "cage", "timepoint", "month")


#### plot the trajectories  ####
plot.SNP1181 <- ggplot(aes(y=frequency, x=month, group=cage_name), data=SNP.1181.freqs.t)+
  geom_line(aes(color=cage_name), size=1.6)+
  geom_point(size=4)+
  ggtitle("18,521,181")+
  ylab("Frequency of the L allele")+
  xlab("Months since the start of the experiment")+
  theme_bw()+
  theme(plot.title = element_text(size=30, face="bold"))+
  theme(axis.text = element_text(size=30),
        axis.title.x = element_text(size=45), 
        axis.title.y = element_text(size=45),strip.text=element_text(size=30))
plot.SNP1181

plot.SNP1208 <- ggplot(aes(y=frequency, x=month, group=cage_name), data=SNP.1208.freqs.t)+
  geom_line(aes(color=cage_name), size=1.6)+
  geom_point(size=4)+
  ggtitle("18,521,208")+
  ylab("Frequency of the L allele")+
  xlab("Months since the start of the experiment")+
  theme_bw()+
  theme(plot.title = element_text(size=30, face="bold"))+
  theme(axis.text = element_text(size=30),
        axis.title.x = element_text(size=45), 
        axis.title.y = element_text(size=45),strip.text=element_text(size=30))
plot.SNP1208

plot.SNP1521 <- ggplot(aes(y=frequency, x=month, group=cage_name), data=SNP.1521.freqs.t)+
  geom_line(aes(color=cage_name), size=1.6)+
  geom_point(size=4)+
  ggtitle("18,521,521")+
  ylab("Frequency of the L allele")+
  xlab("Months since the start of the experiment")+
  theme_bw()+
  theme(plot.title = element_text(size=30, face="bold"))+
  theme(axis.text = element_text(size=30),
        axis.title.x = element_text(size=45), 
        axis.title.y = element_text(size=45),strip.text=element_text(size=30))
plot.SNP1521


### different sort of plot

plot.SNP1181.groupmean <- ggplot(aes(x=month, y=frequency, colour=start.freq), data=SNP.1181.freqs.t)+
  geom_line()+
  geom_point()+
  geom_point(aes(x=month, y=frequency, colour=start.freq), data=SNP.1181.freqs.t)
plot.SNP1181.groupmean


#### add startign frequency information ####
## starting frequency information
start.freq <- 1:90
for(i in 1:nrow(SNP.1181.freqs.t)){
  if(SNP.1181.freqs.t$cage[i] %% 2){
    start.freq[i] <- "high L start"
  } else {
    start.freq[i] <- "low L start"
  }
  }
start.freq
SNP.1181.freqs.t$start.freq <- start.freq
SNP.1208.freqs.t$start.freq <- start.freq
SNP.1521.freqs.t$start.freq <- start.freq

### write all as csv files
write.csv(SNP.1181.freqs.t, "proxy_SNP_data/1181_metrics.csv")
write.csv(SNP.1208.freqs.t, "proxy_SNP_data/1208_metrics.csv")
write.csv(SNP.1521.freqs.t, "proxy_SNP_data/1521_metrics.csv")


#### recording metrics of fru data ####
### these are alos in a R scripts called 'functions for trajectory analysis'

## create summary table for metrics
Cage.ID <- c("Cage1", "Cage2", "Cage3", "Cage4", "Cage5", "Cage6", "Cage7", "Cage8", "Cage9", "Cage10")
Initial_state <- c("high_L", "low_L", "high_L", "low_L", "high_L", "low_L", "high_L", "low_L", "high_L", "low_L")
cage_metrics <- data.frame(Cage.ID = Cage.ID, Initial_state = Initial_state)
str(cage_metrics)

### 1. is either allele fixed? - last two time points frequcy equals 1 or 0
## aso in anther script called "functions for tajectory analysis"

fixation.fn <- function(data){
  fixation <- 2
  if(data[data$timepoint==8,]$frequency == 1 && data[data$timepoint==9,]$frequency == 1){
    fixation <- 1
  } else {
    fixation <- 0
  } |
    if(data[data$timepoint==8,]$frequency == 0 && data[data$timepoint==9,]$frequency == 0){
      fix <- 1
    } else {
      fixation <- 0
    }
  return(fixation)
}


### 2. persistance time
## how long in months from the start of the experiment until either allele fixes
persistance.time.fn <- function(data){
  persistance.time <- 22
  for(i in 2:nrow(data)){
if(persistance.time!=22) break    
if(data[data$timepoint==i-1,]$frequency == 1 && data[data$timepoint==i,]$frequency == 1){
    persistance.time <- data[data$timepoint==i-1,]$month
  } else {} |
      if(data[data$timepoint==i-1,]$frequency == 0 && data[data$timepoint==i,]$frequency == 0){
                persistance.time <- data[data$timepoint==i-1,]$month
              } else {}
    }
    return(persistance.time)
}


### 3. chnage from starting frequency
#raw diff between timepoint 9 frequency and 0.1 or 0.9 (0.9 if high L start == odd, probably, could be other way)

diff.start.fn <- function(data){
  #set difference to 0
  diff <- 0
  #if it's an even cage (remainder of cage number/2 is 0)
  if(data$cage[1] %% 2==0){
    #find the difference between the frequency at timepoint 9 and 0.1 (starting frequency)
    diff <- data[data$timepoint==9,]$frequency-0.1
  } else{
    #otherwise find the difference between the frequency at timepoint 9 and 0.9 (starting frequency)
    diff <- data[data$timepoint==9,]$frequency-0.9}
  return(diff)
}

### 4. slope from timepoint 1 to timepoint 9 - just those two, don't fit the other points
#CHANGED TO MONTH
slope.t1t9.fn <- function(data){
  #initially set the slope to 0
  slope <- 0
  #make a new data frame (sub) with only timepoints 1 and 9
  sub <- subset(data, timepoint%in%c(1,9))
  #find the slope of a linear regression through timepoints 1 and 9
  slope <- lm(frequency~month, data=sub)$coefficients[2]
  return(slope)
}

### 5. best fit line (regression slope) with all the points
#CHANGED TO MONTH
slope.all.fn <- function(data){
  #initially set the slope to 0
  slope <- 0
  #find the slope of a linear regression through all timepoints
  slope <- lm(frequency~month, data=data)$coefficients[2]
  return(slope)
}

### sum of squares distance from slope (did sum of squares error, i think)
##DEFINITELY CHECK THIS ONE
sse.fn <- function(data){
  #initially set the sum of squares error to 0
  sse <- 0
  #find the linear regression through all timepoints
  md <- lm(frequency~month, data=data)
  #sum of squares error is the sum of the squares of the residuals
  sse <- sum((md$residuals)^2)
  
  return(sse)
}

### variance of distance of points to line
#should be the variance of the residuals?
##DEFINITELY CHECK THIS ONE
var.error.fn <- function(data){
  #initially set the variance in error to 0
  var.e <- 0
  #find the linear regression through all timepoints
  md <- lm(frequency~month, data=data)
  #variance of the error is the variance of the residuals
  var.e <- var(md$residuals)
  
  return(var.e)
}

## 9. find all rows where timepoint is 9
tp9 <- subset(SNP.1181.freqs.t, timepoint==9)
## 9.1 difference in variance from tp.1 and tp.9
tp1 <- subset(SNP.1181.freqs.t, timepoint==1)
var.diff <- (var(tp1$frequency)) - (var(tp9$frequency))
var.diff

## use the functions to calculate the values for our cage data 
## one SNP at a time
for (i in 1:nrow(cage_metrics)){
  cage_metrics$fixation[i] <- fixation.fn(subset(SNP.1181.freqs.t,cage==i))
  cage_metrics$persistance_time[i] <- persistance.time.fn(subset(SNP.1181.freqs.t,cage==i))
  cage_metrics$diff.start[i] <- diff.start.fn(subset(SNP.1181.freqs.t,cage==i))
  cage_metrics$slope.t1t9[i] <- slope.t1t9.fn(subset(SNP.1181.freqs.t,cage==i))
  cage_metrics$slope.all[i] <- slope.all.fn(subset(SNP.1181.freqs.t,cage==i))
  cage_metrics$sse[i] <- sse.fn(subset(SNP.1181.freqs.t,cage==i))
  cage_metrics$var.err[i] <- var.error.fn(subset(SNP.1181.freqs.t,cage==i))
  cage_metrics$zigzag[i] <- zigzag.fn(subset(SNP.1181.freqs.t,cage==i))
  cage_metrics$endvar[i] <- var(tp9$frequency)
  cage_metrics$startvar[i] <- var(tp1$frequency)
  cage_metrics$vardiff[i] <- (var(tp1$frequency)) - (var(tp9$frequency))
}
## will create a function to do this for the simaultd data, but its in a slightly different format so keep the loop here

