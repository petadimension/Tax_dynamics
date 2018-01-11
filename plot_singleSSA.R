rm( list=ls(all=TRUE) )

library(ggplot2)

is_skip <- "yes"
skipnum <- 500
tscale <- 0.00065

data <- read.table("Taxsim.txt",sep="\t",header=TRUE)
data <- data[seq(1,nrow(data),skipnum),]
timev <- data[,1]*tscale
numv <- data[,2]
labelv <- rep("Tax",length=nrow(data))
df <- data.frame(time=timev,number=numv,variable=labelv)

sfn <- "bad.eps"
postscript(sfn,horizontal=FALSE,onefile=FALSE,paper="special",height=9,width=9)
#sfn <- "Tax_SSA_singlerun.png"
#png(sfn,height=600,width=600)
theme_set(theme_gray(base_size=18))
plt <- ggplot() + 
geom_line(data=df,aes(x=time,y=number,color=variable),lwd=1) +
ylim(c(0.01,max(numv))) + 
xlab("time (hours)") + ylab("Number of molecules") + labs(title="Tax dynamics")
print( plt )
dev.off()

data <- read.table("Taxsim.txt",sep="\t",header=TRUE)
zeroidx <- which( data[,2]==100 )
print( data[zeroidx,1]*tscale )
dummy <- c(0,0)
if( is_skip=="yes" ) {
	tmp1 <- data[zeroidx[1]:zeroidx[4],]
	tmp2 <- data[zeroidx[4]:zeroidx[7],]
	tmp3 <- data[zeroidx[10]:zeroidx[11],]
	notmp1 <- data[1:zeroidx[1],]
	notmp2 <- data[zeroidx[7]:zeroidx[10],]
	notmp3 <- data[zeroidx[11]:nrow(data),]

	dummy <- rbind(dummy,notmp1)
	tmp1 <- tmp1[seq(1,nrow(data),skipnum),]
	dummy <- rbind(dummy,tmp1)
	tmp2 <- tmp2[seq(1,nrow(data),skipnum),]
	dummy <- rbind(dummy,tmp2)
	dummy <- rbind(dummy,notmp2)
	tmp3 <- tmp3[seq(1,nrow(data),skipnum),]
	dummy <- rbind(dummy,tmp3)
	dummy <- rbind(dummy,notmp3)
#	for( i in 1:(length(zeroidx)-1) ) {
#		tmp <- data[zeroidx[i]:zeroidx[i+1],]
#		if( (data[zeroidx[i+1],1]-data[zeroidx[i],1])*tscale>10.0 ) {
#			tmp <- tmp[seq(1,nrow(data),skipnum),]
#			dummy <- rbind(dummy,tmp)
#		} else {
#			dummy <- rbind(dummy,tmp)
#		}
#	}
}
data <- dummy[-1,]

timev <- data[,1]*tscale
numv <- data[,2]
labelv <- rep("Tax",length=nrow(data))
df <- data.frame(time=timev,number=numv,variable=labelv)

sfn <- "Tax_SSA_singlerun.eps"
postscript(sfn,horizontal=FALSE,onefile=FALSE,paper="special",height=9,width=9)
#sfn <- "Tax_SSA_singlerun.png"
#png(sfn,height=600,width=600)
theme_set(theme_gray(base_size=18))
plt <- ggplot() + 
geom_line(data=df,aes(x=time,y=number,color=variable),lwd=1) +
ylim(c(0.01,max(numv))) + 
xlab("time (hours)") + ylab("Number of molecules") + labs(title="Tax dynamics")
print( plt )
dev.off()
