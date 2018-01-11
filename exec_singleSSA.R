rm( list=ls(all=TRUE) )

#rnd_seed <- 1234
rnd_seed <- 376
parms <- c(kon=3.0e-6,koff=0.01,km=0.1,kp=10.0,kbind=0.0025,kunbind=0.018,ka=5.0,dm=1.0,dp=0.125)
cmd <- paste("./singlesim.out",rnd_seed,paste(parms,collapse=" "))
system(cmd)
