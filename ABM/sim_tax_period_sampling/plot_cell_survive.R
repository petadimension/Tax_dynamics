rm( list=ls(all=TRUE) )

library(ggplot2)
library(scales)

seedn <- 123
NinitNC <- 104 # total number of cells per well
NinitKD <- 63 # total number of cells per well

Tmin <- 3.0 # start of experiment
Tmax <- 20.0 # 4.0 # (hours) end of experiment
PTmax <- 19.0
ptime0 <- seq(Tmin,Tmin+1,by=0.1) # time point for plot
ptime1 <- seq(Tmin+1,Tmax,by=0.1) # time point for plot
ptime <- c(ptime0,ptime1)

## experimental data ##
# shNC cell count #
dpNC <- read.table(file="../data/MOCKcell.csv",sep=",",header=F)
dpNC <- cbind(dpNC,rep("shNC",length=nrow(dpNC)))
colnames(dpNC) <- c("days","values","label")
# Tax KD cell count #
dpKD <- read.table(file="../data/TaxKDcell.csv",sep=",",header=F)
dpKD <- cbind(dpKD,rep("TaxKD",length=nrow(dpKD)))
colnames(dpKD) <- c("days","values","label")
# bind #
dp <- rbind(dpNC,dpKD)
ov20 <- which( dp[,1]>20.0 )
dp <- dp[-ov20,]
ud3 <- which( dp[,1]<3.0 )
dp <- dp[-ud3,]
# shNC cell count #
gfpNC <- read.table(file="../data/MOCKGFP.csv",sep=",",header=F)
gfpNC <- cbind(gfpNC,rep("shNC-GFP+",length=nrow(gfpNC)))
colnames(gfpNC) <- c("days","values","label")
# Tax KD cell count #
gfpKD <- read.table(file="../data/TaxKDGFP.csv",sep=",",header=F)
gfpKD <- cbind(gfpKD,rep("TaxKD-GFP+",length=nrow(gfpKD)))
colnames(gfpKD) <- c("days","values","label")
# bind #
dpGFP <- rbind(gfpNC,gfpKD)
ov20 <- which( dpGFP[,1]>20.0 )
dpGFP <- dpGFP[-ov20,]
ud3 <- which( dpGFP[,1]<3.0 )
dpGFP <- dpGFP[-ud3,]

## computed MOCK cell count ##
rfn <- paste("TaxOnOffApo_num",NinitNC,"seed",seedn,".txt",sep="")
tmp <- read.table(rfn,sep="\t",header=TRUE)
stime <- tmp[,1]
Total_NC <- tmp[,2]+tmp[,3]+tmp[,4]
# choose the nearest time point #
tn <- length(ptime)
ids <- rep(0,length=tn)
for( i in 1:tn ) {
    ids[i] <- which( abs(stime-ptime[i])==min(abs(stime-ptime[i])) )
}
#nonskipids <- 1:(ids[1]-1)
#ids <- c(nonskipids,ids)
#tn <- length(ids)
#ptime <- c(tmp[nonskipids,1],ptime)
## cell count data frame ##
dfNC <- data.frame(days=ptime,values=Total_NC[ids],label=rep("shNC",length=tn))
ovPTmax <- which( dfNC$days>PTmax )
dfNC <- dfNC[-ovPTmax,]

tmp2 <- tmp[ids,2]
zeroidx <- which( tmp2<1.0 )
tmp2[zeroidx] <- 1
dfNCOn <- data.frame(days=ptime,values=tmp2,label=rep("NC-On",length=tn))
dfNCOff <- data.frame(days=ptime,values=tmp[ids,3],label=rep("NC-Off",length=tn))
dfNCApo <- data.frame(days=ptime,values=tmp[ids,4],label=rep("NC-Apo",length=tn))
dfNCAll <- rbind(dfNCOn,dfNCOff,dfNCApo)
ovPTmax <- which( dfNCAll$days>PTmax )
dfNCAll <- dfNCAll[-ovPTmax,]
## Tax-on ratio data frame ##
OnRatio <- data.frame(days=ptime,values=tmp[ids,2]/Total_NC[ids],label=rep("TaxOnRatio",length=tn))
ovPTmax <- which( OnRatio$days>PTmax )
OnRatio <- OnRatio[-ovPTmax,]

sfn <- paste("ABM_shNCAll.eps",sep="")
postscript(sfn,horizontal=FALSE,onefile=FALSE,paper="special",height=9,width=9)
theme_set(theme_gray(base_size=18))
plt <- ggplot() + 
#geom_point(data=dp,aes(x=days,y=log(values,base=10),shape=label),size=12) +
geom_path(data=dfNC,aes(x=days,y=log(values,base=10)),lwd=2) +
geom_path(data=dfNCAll,aes(x=days,y=log(values,base=10),linetype=label,color=label),lwd=2) +
scale_y_continuous(breaks=0:7,labels=trans_format('+',math_format(10^.x))) +
#coord_trans(y="log10")+ scale_y_continuous(trans=log10_trans(),breaks=trans_breaks("log10",function(x) 10^x),labels=trans_format("log10",math_format(10^.x))) +
xlab("time (days)") + ylab("Number of shNC TaxOn / Off / Apoptotic cells") + labs(title="Cell level simulation") #+ 
#ylim(c(0.0,6.5)) #+ xlim(c(Tmin,Tmax))
print( plt )
dev.off()

## computed TaxKD cell count ##
rfn <- paste("TaxOnOffApo_num",NinitKD,"seed",seedn,"_TaxKD.txt",sep="")
tmp <- read.table(rfn,sep="\t",header=TRUE)
stime <- tmp[,1]
Total_KD <- tmp[,2]+tmp[,3]+tmp[,4]
# choose the nearest time point #
tn <- length(ptime)
ids <- rep(0,length=tn)
for( i in 1:tn ) {
    ids[i] <- which( abs(stime-ptime[i])==min(abs(stime-ptime[i])) )
}
#nonskipids <- 1:(ids[1]-1)
#ids <- c(nonskipids,ids)
#tn <- length(ids)
#ptime <- c(tmp[nonskipids,1],ptime0)

tmp2 <- tmp[ids,2]
zeroidx <- which( tmp2<1.0 )
tmp2[zeroidx] <- 1
dfKD <- data.frame(days=ptime,values=Total_KD[ids],label=rep("TaxKD",length=tn))
dfKDOn <- data.frame(days=ptime,values=tmp2,label=rep("KD-On",length=tn))
dfKDOff <- data.frame(days=ptime,values=tmp[ids,3],label=rep("KD-Off",length=tn))
dfKDApo <- data.frame(days=ptime,values=tmp[ids,4],label=rep("KD-Apo",length=tn))
dfKDAll <- rbind(dfKDOn,dfKDOff,dfKDApo)
## merge ##
dfmix <- rbind(dfNC,dfKD)

sfn <- paste("ABM_TaxKDAll.eps",sep="")
postscript(sfn,horizontal=FALSE,onefile=FALSE,paper="special",height=9,width=9)
theme_set(theme_gray(base_size=18))
plt <- ggplot() + 
#geom_point(data=dp,aes(x=days,y=log(values,base=10),shape=label),size=12) +
geom_path(data=dfKD,aes(x=days,y=log(values,base=10)),lwd=2) +
geom_path(data=dfKDAll,aes(x=days,y=log(values,base=10),linetype=label,color=label),lwd=2) +
scale_y_continuous(breaks=0:5,labels=trans_format('+',math_format(10^.x))) +
xlab("time (days)") + ylab("Number of TaxKD TaxOn / Off / Apoptotic cells") + labs(title="Cell level simulation") #+
#ylim(c(0.0,5.0))
print( plt )
dev.off()

## GFP+ percentage data frame ##
modratio <- 0.455
rmidx <- which( ptime>PTmax )
dfGFPNC <- data.frame(days=ptime[-rmidx],values=modratio*100*dfNC$values/(dfNC$values),label=rep("shNC-GFP+",length=(tn-length(rmidx))))
dfGFPKD <- data.frame(days=ptime[-rmidx],values=100*dfKD$values[-rmidx]/(dfNC$values+dfKD$values[-rmidx]),label=rep("TaxKD-GFP+",length=(tn-length(rmidx))))
## merge ##
dfGFP <- rbind(dfGFPNC,dfGFPKD)

## plot cell count (log-scale) ##
sfn <- paste("ABM_cellcount.eps",sep="")
postscript(sfn,horizontal=FALSE,onefile=FALSE,paper="special",height=9,width=9)
#png(sfn,width=800,height=800)
theme_set(theme_gray(base_size=18))
plt <- ggplot() +
geom_point(data=dp,aes(x=days,y=log(values,base=10),shape=label),size=12) +
#geom_point(data=dfmix,aes(x=days,y=log(values,base=10),shape=label,color=label),size=2) +
geom_path(data=dfmix,aes(x=days,y=log(values,base=10),linetype=label,color=label),lwd=2) +
scale_y_continuous(breaks=0:7,labels=trans_format('+',math_format(10^.x))) +
#ylim(c(0.0,6.5)) +
xlab("time (days)") + ylab("Number of shNC and TaxKD cells") + labs(title="Cell level simulation") #+
print( plt )
dev.off()

## plot cell count (log-scale) ##
sfn <- paste("ABM_GFP.eps",sep="")
postscript(sfn,horizontal=FALSE,onefile=FALSE,paper="special",height=9,width=9)
#png(sfn,width=800,height=800)
theme_set(theme_gray(base_size=18))
plt <- ggplot() +
geom_point(data=dpGFP,aes(x=days,y=values,shape=label),size=12) +
geom_path(data=dfGFP,aes(x=days,y=values,linetype=label,color=label),lwd=2) +
#xlim(c(Tmin,Tmax)) + 
xlab("time (days)") + ylab("Percentage of GFP+ cells") + labs(title="Cell level simulation") #+
print( plt )
dev.off()

## plot Tax-on ratio ##
sfn <- paste("TaxOnRatio.eps",sep="")
postscript(sfn,horizontal=FALSE,onefile=FALSE,paper="special",height=9,width=9)
#png(sfn,width=800,height=800)
theme_set(theme_gray(base_size=18))
plt <- ggplot() +
geom_path(data=OnRatio,aes(x=days,y=values),lwd=2,color="black") +
#xlim(c(Tmin,Tmax)) + 
ylim(c(0.0,0.15)) +
xlab("time (days)") + ylab("Ratio of Tax-on cells") + labs(title="Cell level simulation") #+
print( plt )
dev.off()
