rm( list=ls(all=TRUE) )

get_alivecell_IDs <- function( lpop ) {
	ID_vec <- lpop[[1]]
	alive_IDs <- which( is.na(ID_vec)==FALSE )
	return( alive_IDs  )
}

get_deadcell_IDs <- function( lpop ) {
	ID_vec <- lpop[[1]]
	dead_IDs <- which( is.na(ID_vec)==TRUE )
	return( dead_IDs )
}

get_TaxOn_IDs <- function( lpop ) {
	state_vec <- lpop[[2]]
	OnIDs <- which( state_vec=="TaxOn" )
	aidx <- get_alivecell_IDs(lpop)
	OnIDs <- intersect(OnIDs,aidx)
	return( OnIDs )
}

get_TaxOff_IDs <- function( lpop ) {
	state_vec <- lpop[[2]]
	OffIDs <- which( state_vec=="TaxOff" )
	aidx <- get_alivecell_IDs(lpop)
	OffIDs <- intersect(OffIDs,aidx)
	return( OffIDs )
}

get_Apoptotic_IDs <- function( lpop ) {
	state_vec <- lpop[[2]]
	ApoIDs <- which( state_vec=="Apop" )
	aidx <- get_alivecell_IDs(lpop)
	ApoIDs <- intersect(ApoIDs,aidx)
	return( ApoIDs )
}

get_TaxOnhist_IDs <- function( lpop ) {
	hist_vec <- lpop[[4]]
	histIDs <- which( hist_vec=="yes" )
	aidx <- get_alivecell_IDs(lpop)
	histIDs <- intersect(histIDs,aidx)
	return( histIDs )
}

death_process <- function( lpop, tID ) {
	ID_vec <- lpop[[1]]
	ID_vec[tID] <- NA
	lpop[[1]] <- ID_vec
	return( lpop )
}

division_process <- function( lpop, mID, dID ) {
	ID_vec <- lpop[[1]]
	state_vec <- lpop[[2]]
	clock_vec <- lpop[[3]]
	hist_vec <- lpop[[4]]
	ID_vec[dID] <- dID
	state_vec[dID] <- state_vec[mID]
	hist_vec[dID] <- hist_vec[mID]
	if( state_vec[mID]=="TaxOn" ) {
		clock_vec[dID] <- sample(t_periods,1,replace=TRUE)
	} else if( state_vec[mID]=="TaxOff" ) {
		clock_vec[dID] <- sample(t_intervals,1,replace=TRUE)
	} else {
		clock_vec[dID] <- Tmax # nothing to do with apoptotic cells
	}
	lpop[[1]] <- ID_vec
	lpop[[2]] <- state_vec
	lpop[[3]] <- clock_vec
	lpop[[4]] <- hist_vec
	return( lpop )
}

newborn_registry <- function( lpop, newID ) {
	ID_vec <- lpop[[1]]
	state_vec <- lpop[[2]]
	clock_vec <- lpop[[3]]
	hist_vec <- lpop[[4]]
	ID_vec <- c(ID_vec,newID)
	state_vec <- c(state_vec,NA)
	clock_vec <- c(clock_vec,NA)
	hist_vec <- c(hist_vec,NA)
	lpop[[1]] <- ID_vec
	lpop[[2]] <- state_vec
	lpop[[3]] <- clock_vec
	lpop[[4]] <- hist_vec
	return( lpop )
}

update_clock <- function( lpop, time_step ) {
	clock_vec <- lpop[[3]]
	aidx <- get_alivecell_IDs( lpop )
	clock_vec[aidx] <- clock_vec[aidx]-time_step
	lpop[[3]] <- clock_vec
	return( lpop )
}

change_TaxOnOff <- function( lpop ) {
	state_vec <- lpop[[2]]
	clock_vec <- lpop[[3]]
	hist_vec <- lpop[[4]]
	tidx <- which( clock_vec<0 )
	if( length(tidx)>0 ) {
		for( tid in tidx ) {
			if( state_vec[tid]=="TaxOn" ) {
				state_vec[tid] <- "TaxOff"
				clock_vec[tid] <- Tmax # no transition to TaxOn #sample(t_intervals,1,replace=TRUE)
			} #else if( state_vec[tid]=="TaxOff" ) {
            #    state_vec[tid] <- "TaxOn"
            #    clock_vec[tid] <- sample(t_periods,1,replace=TRUE)
            #    hist_vec[tid] <- "yes"
			#} else if( state_vec[tid]=="Apop" ) {
            #    state_vec[tid] <- "TaxOn"
            #    clock_vec[tid] <- sample(t_periods,1,replace=TRUE)
            #    hist_vec[tid] <- "yes"
			#}
		}
	}
	lpop[[2]] <- state_vec
	lpop[[3]] <- clock_vec
	lpop[[4]] <- hist_vec
	return( lpop )
}

change_TaxOff2Apo <- function( lpop, tid ) {
	state_vec <- lpop[[2]]
	clock_vec <- lpop[[3]]
	state_vec[tid] <- "Apop"
	#clock_vec[tid] <- sample(t_periods,1,replace=TRUE)
	clock_vec[tid] <- Tmax # Apoptotic cells never come back to TaxOff
	lpop[[2]] <- state_vec
	lpop[[3]] <- clock_vec
	return( lpop )
}

seedn <- 123 #123456
set.seed(seedn) # set a seed of random number generator
Tmin <- 3.0 # start of experiment
Tmax <- 21.0 # (hours) end of experiment
min_step <- 1.0
skipnum <- 10000
Nthr <- 10000
#kbind=0.0025_kunbind=0.018_ka=5.0

# 2.995819	0	17	46	0
Tstart <- 2.995819
NonInit <- 6
NoffInit <- 1
NapoInit <- 56
Ninit <- NonInit+NoffInit+NapoInit # total number of cells per well
gon <- 0.54 # 0.52 # cell division rate of TaxOn cells (1/day)
g <- 0.68 # 0.65 # cell division rate (1/day)
delta <- 0.25 # 0.6 # additional death via apoptosis (1/day)
p <- 0.40 # 0.40 # transition from TaxOff to Apoptotic cell

is_expdata <- "no"
if( is_expdata=="yes" ) {
	t_periods <- as.numeric( read.table(file=rfn <- "../data/t_period.txt",header=T)$t_period ) ## exp
} else {
	t_periods <- read.table("optim_t_period.txt",header=F)$V1 ## sim
}

t_periods <- t_periods/24.0 # convert to day unit
t_intervals <- Tmax # no transition to TaxOn #read.table("optim_t_interval.txt",header=F)$V1
#t_intervals <- t_intervals/24.0 # convert to day unit
Onpercent <- 0.02
Apopercent <- 0.2

if( is_expdata!="yes" ) {
    over60idx <- which( t_periods*24.0>60.0 )
    if( length(over60idx)>0 ) {
    	t_periods <- t_periods[-over60idx]
    }
    less3idx <- which( t_periods*24.0<3.0 )
    if( length(less3idx)>0 ) {
    	t_periods <- t_periods[-less3idx]
    }
}


#library(MASS)
#tint <- read.table("optim_t_interval.txt",header=F)$V1/24
#res <- fitdistr(tint,"exponential",rate=1.0)
#print( res )
#png("t_interval_expfit.png",width=600,height=600)
#hist(tint,freq=FALSE)
#x <- seq(0.0,350,by=1.0)
#lines(x,dexp(x,rate=0.0157),lwd=2,col="blue")
#dev.off()
#
#res <- fitdistr(t_periods,"exponential",rate=1.0)
#print( res )
#png("t_period_expfit.png",width=600,height=600)
#hist(t_periods,freq=FALSE)
#x <- seq(0.0,2.5,by=0.1)
#lines(x,dexp(x,rate=0.93),lwd=2,col="blue")
#print( summary(t_periods) )
#print( sqrt(1.0/var(t_periods)) )
#dev.off()
#
#quit()

## Define a population as a list of dataframes ##
IDs <- 1:Ninit
#InitialOnIDs <- sample(IDs,Onpercent*Ninit,replace=FALSE)
InitialOnIDs <- 1:NonInit
tmp <- setdiff(IDs,InitialOnIDs)
InitialOffIDs <- (NonInit+1):(NonInit+NoffInit+1)
InitialApoIDs <- (NonInit+NoffInit+1+1):Ninit
states <- rep("",length=Ninit)
states[InitialOnIDs] <- "TaxOn"
states[InitialOffIDs] <- "TaxOff"
states[InitialApoIDs] <- "Apop"
clocks <- rep(0.0,length=Ninit)
clocks[InitialOnIDs] <- sample(t_periods,length(InitialOnIDs),replace=TRUE)
clocks[InitialOffIDs] <- Tmax # no transition to TaxOn #sample(t_intervals,length(InitialOffIDs),replace=TRUE)
clocks[InitialApoIDs] <- Tmax # no transition to TaxOn #sample(t_intervals,length(InitialOffIDs),replace=TRUE)
histories <- rep("",length=Ninit)
histories[InitialOnIDs] <- "yes"
histories[InitialOffIDs] <- "no"
histories[InitialApoIDs] <- "no"
dfpop <- data.frame(ID=IDs,state=states,clock=clocks,TaxOnhistory=histories,stringsAsFactors=FALSE)
lpop <- as.list(dfpop)

Nalive <- Ninit

t <- Tstart
rates <- c(0.0,0.0,0.0,0.0,0.0)
# 1: cell division of TaxOn
# 2: cell division of TaxOff
# 3: cell division of Apoptotic
# 4: cell death of Apoptotic
# 5: Transition from TaxOff to Apoptotic
rates[1] <- gon*length(InitialOnIDs) #g*Non
rates[2] <- g*length(InitialOffIDs)
rates[3] <- g*length(InitialApoIDs)
rates[4] <- delta*length(InitialApoIDs)
rates[5] <- p*length(InitialOffIDs)
trates <- rates[1]+rates[2]+rates[3]+rates[4]+rates[5]
pp <- rates/trates

## output ##
sfn <- paste("TaxOnOffApo_num",Ninit,"seed",seedn,"_TaxKD.txt",sep="")
pout <- "time\tNOn\tNOff\tApo\thist\n"
cat(pout,file=sfn,append=FALSE)

pout <- c("time (days)","# of TaxOff","# of TaxOn","# of Apoptotic","# of TaxOn history")
print( pout )
out <- c(0,0,0,0,0)
count <- 1
while(  t<Tmax ) {
	if( Nalive>0 ) {
		aliveIDs <- get_alivecell_IDs(lpop)
		deadIDs <- get_deadcell_IDs(lpop)
		Offs <- get_TaxOff_IDs(lpop)
		Ons <- get_TaxOn_IDs(lpop)
		Apos <- get_Apoptotic_IDs(lpop)
		hists <- get_TaxOnhist_IDs(lpop)
		Nalive <- length(aliveIDs)
		Noff <- length(Offs)
		Non <- length(Ons)
		Napo <- length(Apos)
		Nhist <- length(hists)
		if( Nalive>Nthr ) {
			if( count>skipnum ) {
				rates[1] <- gon*Non #g*Non
				rates[2] <- g*Noff
				rates[3] <- g*Napo
				rates[4] <- delta*Napo
				rates[5] <- p*Noff
				trates <- rates[1]+rates[2]+rates[3]+rates[4]+rates[5]
				pp <- rates/trates
				count <- 0
			} else {
				count <- count+1
			}
		} else {
			rates[1] <- gon*Non #g*Non
			rates[2] <- g*Noff
			rates[3] <- g*Napo
			rates[4] <- delta*Napo
			rates[5] <- p*Noff
			trates <- rates[1]+rates[2]+rates[3]+rates[4]+rates[5]
			pp <- rates/trates
		}
		event <- sample(c("gon","goff","gap","dapo","tOff2Apo"),1,prob=pp)
		if( event=="gon" ) {
			## division process of TaxOn cells ##
			motherID <- sample(Ons,1,replace=F) # choice of a mother cell
			if( length(deadIDs)==0 ) {
				daughterID <- Nalive+1
				lpop <- newborn_registry(lpop,daughterID)
				aliveIDs <- c(aliveIDs,daughterID) # temporary but necessary?
			} else if( length(deadIDs)==1 ) {
				daughterID <- deadIDs[1]
				deadIDs <- c() # temporary but necessary
			} else {
				daughterID <- sample(deadIDs,1,replace=F)
				deadIDs <- deadIDs[-daughterID] # becomes empty for future dying cells
			}
			lpop <- division_process(lpop,motherID,daughterID)
		} else if( event=="goff" ) {
			## division process of TaxOff cells ##
			motherID <- sample(Offs,1,replace=F) # choice of a mother cell
			if( length(deadIDs)==0 ) {
				daughterID <- Nalive+1
				lpop <- newborn_registry(lpop,daughterID)
				aliveIDs <- c(aliveIDs,daughterID) # temporary but necessary?
			} else if( length(deadIDs)==1 ) {
				daughterID <- deadIDs[1]
				deadIDs <- c() # temporary but necessary?
			} else {
				daughterID <- sample(deadIDs,1,replace=F)
				deadIDs <- deadIDs[-daughterID] # becomes empty for future dying cells
			}
			lpop <- division_process(lpop,motherID,daughterID)
		} else if( event=="gap" ) {
			## division process of Apoptotitc cells ##
			motherID <- sample(Apos,1,replace=F) # choice of a mother cell
			if( length(deadIDs)==0 ) {
				daughterID <- Nalive+1
				lpop <- newborn_registry(lpop,daughterID)
				aliveIDs <- c(aliveIDs,daughterID) # temporary but necessary?
			} else if( length(deadIDs)==1 ) {
				daughterID <- deadIDs[1]
				deadIDs <- c() # temporary but necessary?
			} else {
				daughterID <- sample(deadIDs,1,replace=F)
				deadIDs <- deadIDs[-daughterID] # becomes empty for future dying cells
			}
			lpop <- division_process(lpop,motherID,daughterID)
		} else if( event=="dapo") {
			## death process of apoptotic cells ##
			go2deathnote <- sample(Apos,1,replace=F)
			deadIDs <- c(deadIDs,go2deathnote) # registry of the dead ID
			lpop <- death_process(lpop,go2deathnote)
		} else if( event=="tOff2Apo" ) {
			## Transition from TaxOff to Apoptotic cells ##
			tID <- sample(Offs,1,replace=F) # choice of a target cell
			lpop <- change_TaxOff2Apo(lpop,tID)
		}
		## calculate time step due to death ##
		dt <- rexp(1,rate=trates)
		time_step <- dt
		## update status ##
		lpop <- change_TaxOnOff(lpop)
		lpop <- update_clock(lpop,time_step)
	} else {
		print( "No alive cells. Finish computation." )
		t <- Tmax
		next
	}
	t <- t+time_step
	out <- c(t,Non,Noff,Napo,Nhist)
	#print( out )
	pout <- sprintf("%f\t%d\t%d\t%d\t%d\n",out[1],out[2],out[3],out[4],out[5])
	cat(pout,file=sfn,append=TRUE)
}
