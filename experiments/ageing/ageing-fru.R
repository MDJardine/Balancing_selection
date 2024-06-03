data<-read.csv("pop_data.csv", header=T, sep=",")
species1<-data$species_1
species2<-data$species_2
year<-data$year
temp<-data$temp[1:19]

#add empty columns into the table to fill in with calculated data later
data<-cbind(data, growth_rate1="0")
data<-cbind(data, growth_rate2="0")
data<-cbind(data, log_growth_rate1="0")
data<-cbind(data, log_growth_rate2="0")

#function for calculating growth rate per year
growth_rate_per_year<-function(vec){
  growth_rate<-vector("numeric", 20)
  for(i in 1:(length(vec)-1)){
    calc<-vec[i+1]/vec[i]
    growth_rate[i]<-calc
  }
  growth_rate
}

#fill in empty columns with calculated data
data$growth_rate1[2:21]<-growth_rate_per_year(data$species_1)
data$growth_rate2[2:21]<-growth_rate_per_year(data$species_2)
data$log_growth_rate1[2:21]<-log10(data$growth_rate1[2:21])
data$log_growth_rate2[2:21]<-log10(data$growth_rate2[2:21])

#function for calculating mean of a vector
mean_func<-function(vec){
  sum_vec<-0
  for(i in 1:length(vec)){
    sum_vec<-sum(sum_vec+vec[i])
  }
  mean_vec<-sum_vec/length(vec)
  return(mean_vec)
}

#intrinsic rate of growth for both species
intrinsic_growth_rate1<-mean_func(data$log_growth_rate1)
intrinsic_growth_rate2<-mean_func(data$log_growth_rate2)

#sum of cross products function
sum_cross_prod_func<-function(vec1, vec2){
  sum_cp<-0
  for (i in 1:length(vec1)){
    cross_prod<-(vec1[i]-mean_func(vec1))*(vec2[i]-mean_func(vec2))
    sum_cp<-sum_cp+cross_prod
  }
  return(sum_cp) 
}

#sum of squares function
sum_sq_func <- function(vec){
  sum_sq <- 0
  for(i in 1:length(vec)){
    sum_sq <- sum_sq + (vec[i] - mean_func(vec) )^2
  }
  return(sum_sq)
}

#PMCC function
cor_coef_func <- function(vec1, vec2){
  cor_coef <- sum_cross_prod_func(vec1, vec2) / 
    sqrt( sum_sq_func(vec1)*sum_sq_func(vec2)) 
  return(cor_coef) 
}

#functions for line of best fit
b_func<-function(vec1, vec2){
  b<-sum_cross_prod_func(vec1, vec2)/sum_sq_func(vec1)
  return(b)
}

a_func<-function(vec1, vec2){
  a<-mean_func(vec2)-(b_func(vec1, vec2)*mean_func(vec1))
  return(a)
}

#plots for temp vs species number
plot(data$temp[1:19], data$log_growth_rate1[1:19], main="temp vs species1", xlab="Temperature difference (˚C)", ylab="Population size", pch=20)
b<-b_func(data$temp[1:19], data$log_growth_rate1[1:19])
a<-a_func(data$temp[1:19], data$log_growth_rate1[1:19])
abline(a, b)

plot(data$temp[1:19], data$log_growth_rate2[1:19], main="temp vs species2", xlab="Temperature difference (˚C)", ylab="Population size", pch=20)
b<-b_func(data$temp[1:19], data$log_growth_rate2[1:19])
a<-a_func(data$temp[1:19], data$log_growth_rate2[1:19])
abline(a, b)

#calculate correlation coefficients
PMCC1<-cor_coef_func(data$temp[1:19], data$log_growth_rate1[1:19])
PMCC2<-cor_coef_func(data$temp[1:19], data$log_growth_rate2[1:19])

#function to calculate t-stat
calc_t_func <- function(vec1, vec2) {
  r <- cor_coef_func(vec1, vec2) 
  t <- r/sqrt((1 - r^2)/(length(vec1) - 2)) 
}

#calculate t_stat for correlation coefficient
t_stat1<-calc_t_func(data$temp[1:19], data$log_growth_rate1[1:19])
t_stat2<-calc_t_func(data$temp[1:19], data$log_growth_rate2[1:19])

#calculate P-values
Pval1<-pt(t_stat1, (length(temp)-2))
Pval2<-pt(t_stat2, (length(temp)-2))

