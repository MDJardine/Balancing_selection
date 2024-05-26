library(learnPopGen)
??learnPopGen

?drift.selection()

## neutral
m1 <- drift.selection(p0=0.9, Ne=500, w=c(1,1,1), ngen=100, nrep=5, colors=NULL)
m1

m2 <- drift.selection(p0=0.1, Ne=500, w=c(1,1,1), ngen=100, nrep=5, colors=NULL)
m2

## directionalto AA
m3 <- drift.selection(p0=0.9, Ne=500, w=c(1,1,0.9), ngen=100, nrep=5, colors=NULL)
m3
m3[1]

m4 <- drift.selection(p0=0.1, Ne=500, w=c(1,1,0.9), ngen=100, nrep=5, colors=NULL)
m4
m4[1]

## directional to aa
m5 <- drift.selection(p0=0.9, Ne=500, w=c(0.9,1,1), ngen=100, nrep=5, colors=NULL)
m5
m5[1]

m6 <- drift.selection(p0=0.1, Ne=500, w=c(0.9,1,1), ngen=100, nrep=5, colors=NULL)
m6
m6[1]

##
m7 <- drift.selection(p0=0.9, Ne=500, w=c(0.9,1,0.9), ngen=100, nrep=5, colors=NULL)
m7
m7[1]

m8 <- drift.selection(p0=0.1, Ne=500, w=c(0.9,1,0.9), ngen=100, nrep=5, colors=NULL)
m8
m8[1]

par(mfrow=c(2,2))


dirft_sim.10K <- 0
for (i in 1:10000){
  dirft_sim.10K[i] <- drift.selection(p0=0.1, Ne=3000, w=c(0.9,1,0.9), ngen=100, nrep=5, colors=NULL)
}

#############################
### neutral ###
m1 <- drift.selection(p0=0.9, Ne=500, w=c(1,1,1), ngen=56, nrep=5000, colors=NULL)
m1

list_length <- 10 
neutral_list <- vector(mode = "list", length = list_length)

for (i in 1:10){
}


  LH <- drift.selection(p0=0.9, Ne=500, w=c(1,1,1), ngen=56, nrep=5)
  LL <- drift.selection(p0=0.1, Ne=500, w=c(1,1,1), ngen=56, nrep=5)
  df5H <- data.frame(matrix(unlist(LH), nrow=length(L), byrow=TRUE))
  df5L <- data.frame(matrix(unlist(LL), nrow=length(L), byrow=TRUE))
  df10 <- rbind(df5H, df5L)
  sel_gen <- select(df10, X3, X5, X9, X13, X21, X29, X37, X45, X57)


  
pool_sim <- function(df, m.cov){
chrosample <- matrix(0, nrow = nrow(df), ncol =ncol(df))
for (i in 1:nrow(df)){
  for (j in 1:ncol(df)){
    chrosample[i,j] <- (rbinom(n=1, size=96, prob = df[i,j])/96)
  }
}
readsample <- matrix(0, nrow = nrow(df), ncol =ncol(df))
for(i in 1:nrow(chrosmaple)){
  for (j in 1:ncol(chrosmaple)){
    readsample[i,j] <- (rbinom(n=1, size=m.cov, prob = chrosmaple[i,j])/35)
  }
}
return(readsample)
}


neut_sim_1 <- pool_sim(sel_gen, 35)


### we now need to take this table of frequencie and tun it into a table similar to that
### we have for the real cage data
### this is 90 rows long and also has extra info needed for the metric calculations

## this should work
neut_sim_freqs <- as.vector(t(neut_sim_1))

## then need the other information
cage <- c(rep("1", 9), rep("2", 9), rep("3", 9), rep("4", 9), rep("5", 9), rep("6", 9),
          rep("7", 9), rep("8", 9), rep("9", 9), rep("10", 9))
time_point <- c(rep(c("1","2","3","4","5", "6","7","8","9"), 10))
month <- c(rep(c("1","2","4","6","10", "14","18","22","28"), 10))
gen <- c(rep(c("2","4","8","12","20", "28","36","44","56"), 10))
freq_start <- c(rep("High_L", 45), rep("Low_L", 45))
exp_info <- cbind(cage, time_point, month, gen, freq_start)

is.matrix(exp_info)
str(exp_info)

### take this part and put it beofore all the loops and simualtions. 
### will be too much to create each time 
## we can just call it instead
## once in this form we should be able to use the metric function the same as before

## will need another line included here that has somethign like
## Nsim <- cbind(neut_sim_freqs, exp_info)

## and then feed that into the metrics bid loop function
## 

## calculate metrics for neutral sims
## need to apadt the loop for all datasets 

