#### ed sim ####
# trying to fix truncated distribution issue w Ed's code
# at this point (8/31), the plug in outperforms the IF-based estimators
# at this point (9/7), I realized that it's fine if the correctly specified plug in is more efficient
# need to figure out what comparisons to do, have as options:
#   - true plug in
#   - estimated plug in (how to estimate?)
#   - true IF
#   - estimated IF

makeSim <- function(N = 10000, psi = 1){
  y0 = rnorm(N,0,1); meanx = c(0,0,0,0)
  alpha = matrix(c(1,1,-1,-1),nrow = 4); beta = matrix(c(1,-1,-1,1))

  x = t(matrix(unlist(lapply(meanx, function(x) rnorm(N,x,1))), nrow = N, byrow =T))
  z = rnorm(N, 1.5*sign(t(alpha)%*%x), sd = 2)
  t = rnorm(N, t(beta)%*%x+y0, sd = 1)
  a = as.numeric(z >= t)
  y = y0 + a*psi*t
  dat = as.data.frame(cbind(y,a,z,t,t(x),y0))

  return(dat)
}

getEstimatorOutputs <- function(data, delta, psi=1){
  alpha = matrix(c(1,1,-1,-1),nrow = 4); beta = matrix(c(1,-1,-1,1))
  attach(data)

  #### true parameters ####
  bX = (t(beta)%*%x)[1,]; muZ = (1.5*sign(t(alpha)%*%x))[1,]; muT = bX + y0

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
  detach(data)

  #### phi single estimates ####
  phi_y = (y - muY)*(pi_min/pi) - (y - muYplus)
  phi_a = (a - muA)*(pi_min/pi) - (a - muAplus)
  psihat.single = mean(phi_y)/mean(phi_a)

  #### phi double estimates ####
  phi_y2 = ( (y - muY)*((pi_min - pi_plus)/pi) ) + muYplus - muYmin
  phi_a2 = ( (a - muA)*((pi_min - pi_plus)/pi) ) + muAplus - muAmin
  psihat.double = mean(phi_y2)/mean(phi_a2)

  # #### standard deviation single ####
  # top = phi_y - psihat.single*phi_a
  # bottom = mean(phi_a2)
  # v.single = mean( ( top/bottom )^2  )/ length(y)
  #
  # #### standard deviation double ####
  # top = phi_y2 - psihat.double*phi_a2
  # bottom = mean(phi_a2)
  # v.double = mean( ( top/bottom )^2  )/ length(y)

  return(data.frame(plugIn = mean(muYplus - muY)/mean(muAplus - muA),
                    IF = psihat.single,
                    #IF.sd = v.single,
                    Empirical = mean(t[which(aplus>a)]),
                    plugInDub = mean(muYplus - muYmin)/mean(muAplus - muAmin),
                    IFdub = psihat.double,
                    #IFdub.sd = v.double,
                    EmpiricalDub = mean(yplus[which(aplus>amin)] - ymin[which(aplus>amin)])
                    )
         )
}

getEstimatorBias <- function(data, delta, psi=1){
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

  #### empirical ####
  Empirical = mean(t[which(aplus>a)])
  EmpiricalDub = mean(yplus[which(aplus>amin)] - ymin[which(aplus>amin)])

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

  return(data.frame(if.bias = IF - Empirical,
                    pi.bias = plugIn - Empirical,
                    #CACE.bias = CACE$phi - Empirical,
                    #piPar.bias = plugInPar$phi - Empirical,
                    if2.bias = IFdub - EmpiricalDub,
                    pi2.bias = plugInDub - EmpiricalDub,
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
getEstimatorBias <- function(data, delta, psi=1){
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
  #plugInDub = mean(muYplus - muYmin)/mean(muAplus - muAmin)
  
  #### empirical ####
  Empirical = mean(t[which(aplus>a)])
  #EmpiricalDub = mean(yplus[which(aplus>amin)] - ymin[which(aplus>amin)])
  
  #### phi single estimates ####
  phi_y = (y - muY)*(pi_min/pi) - (y - muYplus)
  phi_a = (a - muA)*(pi_min/pi) - (a - muAplus)
  IF = mean(phi_y)/mean(phi_a)
  
  #### phi double estimates ####
  # phi_y2 = ( (y - muY)*((pi_min - pi_plus)/pi) ) + muYplus - muYmin
  # phi_a2 = ( (a - muA)*((pi_min - pi_plus)/pi) ) + muAplus - muAmin
  # IFdub = mean(phi_y2)/mean(phi_a2)
  
  #### CACE estimates ####
  CACE = tryCatch({
    CACE(y=y,a=a,z=z,cov=data[,5:8],delta = delta,ranger = T,type = 'single', split = F) 
  } error = function(e){
    print(e)
    return(NA)
  }
  )
   
  #CACEdub = CACE(y=y,a=a,z=z,cov=data[,5:8],delta = delta,ranger = T,type = 'double', split = F)

  #### Parametric plug in estimates ####
  plugInPar = tryCatch({
    CACE(y=y,a=a,z=z,cov=data[,5:8],delta = delta,ranger = T,type = 'plugin', split = F)
  } error = function(e){
    print(e)
    return(NA)
  }
  )
  #plugInParDub = CACE(y=y,a=a,z=z,cov=data[,5:8],delta = delta,ranger = T,type = 'pluginDouble', split = F)
  
  df = data.frame(if.bias = IF - Empirical,
                    pi.bias = plugIn - Empirical,
                    CACE.bias = CACE$phi - Empirical,
                    piPar.bias = plugInPar$phi - Empirical,
                    #if2.bias = IFdub - EmpiricalDub,
                    #pi2.bias = plugInDub - EmpiricalDub,
                    #CACE2.bias = CACEdub$phi - EmpiricalDub,
                    #piPar2.bias = plugInParDub$phi - EmpiricalDub,
                    delta = delta)
  return(df)
  
}

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
clusterExport(cl, "delta")
clusterExport(cl, "CACE")
clusterExport(cl, "getEstimatorBias")
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

# have not run:
run.all <- function(delta){
  #make k datasets (however many simulations you want at each delta)
  k = 50
  N.list <- lapply((1:k), function(x) 5000)
  data.sets <- lapply(N.list, makeSim)
  
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

