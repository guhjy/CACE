#### ed sim ####
# trying to fix truncated distribution issue w Ed's code
# at this point (8/31), the plug in outperforms the IF-based estimators
# at this point (9/7), I realized that it's fine if the correctly specified plug in is more efficient
# need to figure out what comparisons to do, have as options:
#   - true plug in
#   - estimated plug in (how to estimate?)
#   - true IF
#   - estimated IF


# this is a version that doesn't actually run CACE (for speed)
getEstimatorBias <- function(data, delta, true.single, true.double, psi=1){
  alpha = matrix(c(1,1,-1,-1),nrow = 4); beta = matrix(c(1,-1,-1,1))
  x = t(as.matrix(data[,5:8]))
  y = data$y
  a = data$a
  z = data$z
  y0 = data$y0
  t = data$t

  #### true parameters ####
  bX = (t(beta)%*%x)[1,]; muZ = (1.5*sign(t(alpha)%*%x))[1,]; muT = bX + y0

  aplus = as.numeric(z + delta >= t); amin = as.numeric(z - delta >= t)
  yplus = y0 + aplus*psi*t; ymin = y0 + amin*psi*t

  muA = pnorm((z - bX))
  muAplus = pnorm(((z + delta) - bX))
  muAmin = pnorm(((z - delta) - bX))

  muY = bX*c(pnorm((z - bX))) - c(dnorm((z - bX)))#*sqrt(2)
  muYplus = bX*pnorm(((z + delta) - bX)) - dnorm(((z + delta) - bX))#*sqrt(2)
  muYmin = bX*pnorm(((z - delta) - bX)) - dnorm(((z - delta) - bX))#*sqrt(2)

  pi = dnorm(z,mean = muZ, sd = 2)
  pi_min = dnorm(z-delta,mean = muZ, sd = 2)
  pi_plus = dnorm(z+delta,mean = muZ, sd = 2)

  yplus = y0 + aplus*psi*t; ymin = y0 + amin*psi*t

  #### plug-in ####
  plugIn = mean(muYplus - muY)/mean(muAplus - muA)
  plugInDub = mean(muYplus - muYmin)/mean(muAplus - muAmin)

  # #### empirical ####
  # Empirical = mean(t[which(aplus>a)])
  # EmpiricalDub = mean(yplus[which(aplus>amin)] - ymin[which(aplus>amin)])

  #### phi single estimates ####
  phi_y = (y - muY)*(pi_min/pi) - (y - muYplus)
  phi_a = (a - muA)*(pi_min/pi) - (a - muAplus)
  IF = mean(phi_y)/mean(phi_a)

  #### phi double estimates ####
  phi_y2 = ( (y - muY)*((pi_min - pi_plus)/pi) ) + muYplus - muYmin
  phi_a2 = ( (a - muA)*((pi_min - pi_plus)/pi) ) + muAplus - muAmin
  IFdub = mean(phi_y2)/mean(phi_a2)

  # #### CACE estimates ####
  # CACE = CACE(y=y,a=a,z=z,cov=data[,5:8],delta = delta,ranger = T,type = 'single')
  # CACEdub = CACE(y=y,a=a,z=z,cov=data[,5:8],delta = delta,ranger = T,type = 'double')
  #
  # #### Parametric plug in estimates ####
  # plugInPar = CACE(y=y,a=a,z=z,cov=data[,5:8],delta = delta,ranger = T,type = 'plugin')
  # plugInParDub = CACE(y=y,a=a,z=z,cov=data[,5:8],delta = delta,ranger = T,type = 'pluginDouble')

  return(data.frame(if.bias = IF - true.single,
                    pi.bias = plugIn - true.single,
                    #CACE.bias = CACE$phi - Empirical,
                    #piPar.bias = plugInPar$phi - Empirical,
                    if2.bias = IFdub - true.double,
                    pi2.bias = plugInDub - true.double,
                    #CACE2.bias = CACEdub$phi - EmpiricalDub,
                    #piPar2.bias = plugInParDub$phi - EmpiricalDub,
                    delta = delta
                    )
  )

}

library(parallel)
# Calculate the number of cores
no_cores <- detectCores() - 1
# Initiate cluster
cl <- makeCluster(no_cores)
clusterExport(cl, "getEstimatorBias")

#make k datasets (however many simulations you want at each delta)
k = 100
N.list <- lapply((1:k), function(x) 10000)
clusterExport(cl, "N.list");clusterExport(cl, "makeSim")
data.sets <- parLapply(cl,N.list, makeSim)

# output on datasets for delta = 1
delta = 1
clusterExport(cl, "delta")
out.1 <- parLapply(cl, data.sets, function(x) getEstimatorBias(x, delta = delta))
df.1 <- rbind.fill(out.1)
mns.1 <- apply(df.1,2,mean)
sds.1 <- apply(df.1,2,sd)

# to get this file, comment out everything to make CACE, CACE2, piPar, piPar2
# run with delta = 1 over 100 simulated datasets
name = paste("Datasims/simulationsDelta",delta,"NoCACE.csv",sep = "")
write.csv(df.1, file = name)

run.all <- function(delta){
  #make k datasets (however many simulations you want at each delta)
  k = 100
  N.list <- lapply((1:k), function(x) 10000)
  clusterExport(cl, "N.list");clusterExport(cl, "makeSim")
  data.sets <- parLapply(cl,N.list, makeSim)

  # output on datasets for delta
  clusterExport(cl, "delta")
  out <- parLapply(cl, data.sets, function(x) getEstimatorBias(x, delta = delta))
  df <- rbind.fill(out)

  return(df)
}

# get for a range of deltas
delta.range <- seq(0.5,4,by = .5)
out.range <- lapply(delta.range, run.all)

# save the output
for(i in 1:length(out.range)){
  name = paste("Datasims/simulationsDelta",i,"NoCACE.csv",sep = "")
  write.csv(out.range[[i]], file = name)
}

all.out <- rbind.fill(out.range)

mns1 <- ddply(all.out, ~delta, summarize, mean.if = mean(if.bias), mean.pi = mean(pi.bias))
sds1 <- ddply(all.out, ~delta, summarize, sd.if = sd(if.bias), sd.pi = sd(pi.bias))
summarized1 <- merge(mns1,sds1)

# plot
p = ggplot(summarized1, aes(delta)) +
  geom_pointrange(data = summarized1, aes(x = delta, y = mean.pi, ymin = mean.pi - 1.96*sd.pi, ymax= mean.pi + 1.96*sd.pi, linetype = 'Plug in estimate', shape = "Plug in estimate")) +
  geom_pointrange(data = summarized1, aes(x = delta, y = mean.if, ymin = mean.if - 1.96*sd.if, ymax= mean.if + 1.96*sd.if, linetype = 'IF estimate', shape = "IF estimate")) +
  geom_hline( yintercept = 0 )+
  ylab("Bias") +
  xlab("Shift")

ggsave("Datasims/sims100NoCACE.png",plot = p, width = 7, height = 4)

#### try it with the actual estimators ####
#make k datasets (however many simulations you want at each delta)
k = 50
N.list <- lapply((1:k), function(x) 5000)
data.sets <- lapply(N.list, makeSim)

library(parallel)
# Calculate the number of cores
no_cores <- detectCores() - 1
# Initiate cluster
cl <- makeCluster(no_cores)

# output on datasets for delta = 1
library(plyr)
delta = 1
true.ests <- EstimateTrue(delta = 1); t.single <- true.ests[1]; t.double <- true.ests[2]
clusterExport(cl, "delta")
clusterExport(cl, "CACE")
clusterExport(cl, "getEstimatorBias")
clusterExport(cl, "t.single")
clusterExport(cl, "t.double")
start.time <- proc.time()[1]
out.1 <- data.frame(if.bias = rep(NA,length(data.sets)),
                    pi.bias = rep(NA,length(data.sets)),
                    CACE.bias = rep(NA,length(data.sets)),
                    piPar.bias = rep(NA,length(data.sets)),
                    delta = rep(NA,length(data.sets))
                    )
for(i in 20:length(data.sets)){
  out.1[i,] <- getEstimatorBias(data.sets[[i]], delta = 1)
}

# to get this file, leave in all single estimator strategies
# run with delta = 1 over 50 simulated datasets
name = paste("Datasims/simulationsDelta",delta,"CACE.csv",sep = "")
write.csv(out.1, file = name)

# keeps giving a connection error
out.1 <- parLapply(cl, data.sets, function(x) getEstimatorBias(x, delta = delta))
proc.time()[1] - start.time
df.1 <- rbind.fill(out.1)

# to get this file, leave in all estimator strategies
# run with delta = 1 over 2 simulated datasets
name = paste("Datasims/simulationsDeltaParallel",delta,"CACE.csv",sep = "")
write.csv(df.1, file = name)


mns.1 <- apply(df.1,2,mean)
sds.1 <- apply(df.1,2,sd)


# get for a range of deltas
delta.range <- seq(0.5,4,by = .5)
start <- proc.time()
out.range <- lapply(delta.range, run.all)
proc.time- start

# save the output
for(i in 1:length(out.range)){
  name = paste("Datasims/simulationsDeltaBig",i,"CACE.csv",sep = "")
  write.csv(out.range[[i]], file = name)
}

all.out <- rbind.fill(out.range)

mns1 <- ddply(all.out, ~delta, summarize, mean.if = mean(if.bias), mean.pi = mean(pi.bias))
sds1 <- ddply(all.out, ~delta, summarize, sd.if = sd(if.bias), sd.pi = sd(pi.bias))
summarized1 <- merge(mns1,sds1)

# plot
p = ggplot(summarized1, aes(delta)) +
  geom_pointrange(data = summarized1, aes(x = delta, y = mean.pi, ymin = mean.pi - 1.96*sd.pi, ymax= mean.pi + 1.96*sd.pi, linetype = 'Plug in estimate', shape = "Plug in estimate")) +
  geom_pointrange(data = summarized1, aes(x = delta, y = mean.if, ymin = mean.if - 1.96*sd.if, ymax= mean.if + 1.96*sd.if, linetype = 'IF estimate', shape = "IF estimate")) +
  geom_hline( yintercept = 0 )+
  ylab("Bias") +
  xlab("Shift")

ggsave("Datasims/simsBigCACE.png",plot = p, width = 7, height = 4)

#### with glm as propscore ####
dat.sim <- makeSim()

#test with just one
attach(dat.sim)
start <- proc.time()
test.out <- CACE(y=y,a=a,z=z,cov=dat.sim[,5:8],delta = 1,ranger = T,type = 'simple', split = F)
proc.time() - start

# get for a range of deltas
require(plyr)
delta.range <- seq(0,10,by = 2)[-1]
start <- proc.time()
out.range <- lapply(delta.range, run.all)
proc.time() - start

out.df <- rbind.fill(out.range)
write.csv(out.df, 'CACEsingleSims.csv')
summarized <- ddply(out.df, ~delta, summarize, IF=mean(if.bias), IFsd = sd(if.bias),
                    PI = mean(pi.bias), PIsd = sd(pi.bias),
                    CACE = mean(CACE.bias), CACEsd = sd(CACE.bias),
                    piPar = mean(piPar.bias), piParsd = sd(piPar.bias)
                    )
write.csv(summarized, 'CACEsingleSimsSummary.csv')

require(ggplot2)
p = ggplot(summarized, aes(delta)) +
  geom_pointrange(data = summarized, aes(x = delta, y = CACE, ymin = CACE - 1.96*CACEsd, ymax= CACE + 1.96*CACEsd, linetype = 'IF estimate', shape = "IF estimate")) +
  geom_pointrange(data = summarized, aes(x = delta, y = piPar, ymin = piPar - 1.96*piParsd, ymax= piPar + 1.96*piParsd, linetype = 'Plug in estimate', shape = "Plug in estimate")) +
  geom_hline( yintercept = 0 )+
  ylab("Bias") +
  xlab("Shift") 

ggsave("C:/Users/jackie/Desktop/Figures/simsSingleBigCACE.png",plot = p, width = 7, height = 4)