## can ignore this data loading stuff
snp.data <- read.csv("1181_metrics.csv", head=TRUE)

data <- snp.data[snp.data$cage==3,]


#### recording metrics of fru data ####


## create summary table for metrics
Cage.ID <- c("Cage1", "Cage2", "Cage3", "Cage4", "Cage5", "Cage6", "Cage7", "Cage8", "Cage9", "Cage10")
Initial_state <- c("high_L", "low_L", "high_L", "low_L", "high_L", "low_L", "high_L", "low_L", "high_L", "low_L")
cage_metrics <- data.frame(Cage.ID = Cage.ID, Initial_state = Initial_state)
str(cage_metrics)

### 1. is either allele fixed? - last two time points frequcy equals 1 or 0
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
#using 0 as the default because month 22 exists. could choose something else  
  persistance.time <- 0
#start at the end to work our way backwards (i starts at maximum number - number of rows)
  i <- nrow(data)
  #put all the logical statements together in one. 
  #if the end frequency is either 0 or 1 AND the previous time point has the same frequency
    if ((data[data$timepoint==i,]$frequency==1|data[data$timepoint==i,]$frequency==0)&&data[data$timepoint==i-1,]$frequency==data[data$timepoint==i,]$frequency){
      #as long as the previous time point has the same frequency as the focal time point
      while(data[data$timepoint==i-1,]$frequency==data[data$timepoint==i,]$frequency){
        #make the persistance time equal to the month of the previous time point
        persistance.time <- data[data$timepoint==i-1,]$month
        #decrease the focal time point by 1 to check the "while" condition again
        i <- i-1
      }
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

### 7. zigzag metric - how many times does it go up vs how many times does it go down (+1 up, -1 down) - frequency has to change by 0.05
zigzag.fn <- function(data){
  #initially set the 'zigzag' measure to 0
  zz <- 0
  #for all i starting at 2 and ending at the maximum number
  for (i in 2:nrow(data)){
    #temporarily store the difference between the frequency at timepoint i and the previous timepoint
    temp <- data[data$timepoint==i,]$frequency - data[data$timepoint==i-1,]$frequency
    #if the absolute value of the difference is greater than 0.05
    if(abs(temp)>0.05){
      #change the zigzag score according to the sign of the difference. +1 if positive, -1 if negative
      zz <- zz+sign(temp)
    }
  }
  return(zz)
}




## 9. find all rows where timepoint is 9
### variance at timepoint 9 of all of them (include all odds and evens)
#put in table as same value repeated
tp9 <- subset(SNP.1181.freqs.t, timepoint==9)
## 9.1 difference in variance from tp.1 and tp.9
tp1 <- subset(SNP.1181.freqs.t, timepoint==1)
var.diff <- (var(tp1$frequency)) - (var(tp9$frequency))
var.diff



## apply all functions to a data set
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

var(TP1)

#combine all metrics into a single function that can be applied to data (real and siulated)
apply_metrics <- function(df){
  Cage.ID <- c("Cage1", "Cage2", "Cage3", "Cage4", "Cage5", "Cage6", "Cage7", "Cage8", "Cage9", "Cage10")
  metrics <- data.frame(Cage.ID = Cage.ID)
for (i in 1:nrow(cage_metrics)){
  metrics$fixation[i] <- fixation.fn(subset(df,cage==i))
  metrics$persistance_time[i] <- persistance.time.fn(subset(df,cage==i))
  metrics$diff.start[i] <- diff.start.fn(subset(df,cage==i))
  metrics$slope.t1t9[i] <- slope.t1t9.fn(subset(df,cage==i))
  metrics$slope.all[i] <- slope.all.fn(subset(df,cage==i))
  metrics$sse[i] <- sse.fn(subset(df,cage==i))
  metrics$var.err[i] <- var.error.fn(subset(df,cage==i))
  metrics$zigzag[i] <- zigzag.fn(subset(df,cage==i))
  metrics$endvar[i] <- var(df[,9])
  metrics$startvar[i] <- var(df[,1])
  metrics$vardiff[i] <- var(df[,1]) - var(df[,9])
}
}


last.tp <- sel_gen[,9]
var(last.tp)
var(sel_gen[,9])
