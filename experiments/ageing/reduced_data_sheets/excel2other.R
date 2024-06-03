rm(list=ls())

library("rJava")
library("xlsxjars")
library("xlsx")

#setwd("~/Dropbox/students/IbrahimRita/analysis+plots_Rita")

args <- commandArgs(trailingOnly=T)
#args <- c("InR_TF_1.xls", "InR_TF1")

dd <- read.xlsx(args[1], 1, stringsAsFactors=F)[5:83,5:15]
ddNames <- as.character(dd[1,])
dat2 <- data.frame(dd[3:nrow(dd),])	
colnames(dat2) <- ddNames
dim(dat2)
dat2 <- dat2[complete.cases(dat2),]
dim(dat2)
#head(dat2)

deaths <- data.frame(day=dat2$day, event=as.numeric(as.character(c(dat2[,2],dat2[,3],dat2[,4],dat2[,5],dat2[,6],dat2[,7],dat2[,8],dat2[,9],dat2[,10],dat2[,11]))), condition=rep(colnames(dat2)[2:ncol(dat2)], each=nrow(dat2)))
deaths$censor <- 1
dim(deaths)
#str(deaths)
deaths <- subset(deaths, event>0)
dim(deaths)

d <- read.xlsx(args[1], 1, stringsAsFactors=F)[91:166,5:15]
colnames(d) <- ddNames
dim(d)
#head(d)
d <- d[complete.cases(d),]
dim(d)
censors <- data.frame(day=d$day, event=as.numeric(as.character(c(d[,2],d[,3],d[,4],d[,5],d[,6],d[,7],d[,8],d[,9],d[,10],d[,11]))), condition=rep(colnames(d)[2:ncol(d)], each=nrow(d)))
dim(censors)
censors$censor <- 0
censors <- subset(censors, event>0)

colnames(deaths)
colnames(censors)
allOfIt <- data.frame(rbind(deaths, censors))
write.table(allOfIt, sep="\t", file=paste(args[2], "_JMPformat.txt", sep=""), row.names=FALSE, quote=F)

jmp2r <- function(x){
	for(i in 1:nrow(x)){
#		print(paste("Processing row", i))
		xx <- x[i,]
		if(xx$event[1] > 0){
		xxx <- matrix(rep(as.matrix(xx), xx[1,"event"]), ncol=ncol(xx), byrow=T)
		if(exists("a")){a <- rbind(a, xxx)}else{(a <- xxx)}
	}
	}
		a <- as.data.frame(a)
	for(i in 1:ncol(a)){
		colnames(a)[i] <- colnames(x)[i]
		if(class(a[,i]) != class(x[,i])) {a[,i] <- as.character(a[,i]); class(a[,i]) <- class(x[,i])}
	}	

a$censor <- ifelse(a$censor==0, 1, 0)
a <- a[,!colnames(a) %in% c("count", "deathIsZero_censorIsOne")]
a
}

newDat <- jmp2r(allOfIt)
newDat$censor <- factor(newDat$censor)
levels(newDat$censor) <- rev(levels(newDat$censor))
dim(newDat)

write.table(newDat, sep="\t", file=paste(args[2], "_Rformat.txt", sep=""), row.names=FALSE, quote=F)