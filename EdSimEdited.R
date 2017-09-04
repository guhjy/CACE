#### ed sim ####
# trying to fix truncated distribution issue w Ed's code
# at this point (8/31), the plug in outperforms the IF-based estimators

makeSim <- function(){
  rm(list = ls())
  require(truncnorm)
  today = format(Sys.Date(), "%Y%m%d")
  sig2y0 = 1;N = 10000
  y0 = rnorm(N,0,sig2y0); meanx = c(0,0,0,0)
  alpha = matrix(c(1,1,-1,-1),nrow = 4); beta = matrix(c(1,-1,-1,1)); psi =1
  delta = 2

  x = t(matrix(unlist(lapply(meanx,function(x) rnorm(N,x,1))),nrow = N,byrow =T))
  z = rnorm(N,1.5*sign(t(alpha)%*%x),sd = 2)
  t = rnorm(N,t(beta)%*%x+y0,1)
  a = as.numeric(z>=t); aplus = as.numeric(z+delta>=t); amin = as.numeric(z-delta>=t)
  y = y0 + a*psi*t; yplus = y0 + aplus*psi*t; ymin = y0 + amin*psi*t

  dat = as.data.frame(cbind(y,a,z,t,t(x)))
  train = sample(1:dim(dat)[1],8000)
  datTrain = dat[train,]; datTest = dat[-train,]

  #### true parameters ####
  bX = (t(beta)%*%x)[1,]; muZ = (1.5*sign(t(alpha)%*%x))[1,];muT = bX + y0

  muA = pnorm((z-bX));muAplus = pnorm(((z+delta)-bX));muAmin = pnorm(((z-delta)-bX))

  muY = bX*c(pnorm((z-bX))) - c(dnorm((z-bX)))#*sqrt(2)
  muYplus = bX*pnorm(((z+delta)-bX)) - dnorm(((z+delta)-bX))#*sqrt(2)
  muYmin = bX*pnorm(((z-delta)-bX)) - dnorm(((z-delta)-bX))#*sqrt(2)

  pi = dnorm(z,mean = muZ, sd = 2); pi_min = dnorm(z-delta,mean = muZ, sd = 2); pi_plus = dnorm(z+delta,mean = muZ, sd = 2)
  keep = as.numeric( (pi_plus != 0) & (pi_min != 0) )

  yplus = y0 + aplus*psi*t; ymin = y0 + amin*psi*t

  #### phi single estimates ####
  phi_y = (y - muY)*(pi_min/pi) - (y - muYplus)
  phi_a = (a - muA)*(pi_min/pi) - (a - muAplus)

  #### phi double estimates ####
  phi_y2 = ( (y - muY)*((pi_min - pi_plus)/pi) ) + muYplus - muYmin
  phi_a2 = ( (a - muA)*((pi_min - pi_plus)/pi) ) + muAplus - muAmin

  psihat = mean(phi_y2)/mean(phi_a2)

  #### standard deviation ####
  top = phi_y2 - psihat*phi_a2; bottom = mean(phi_a2)
  v = mean( ( top/bottom )^2  )/ length(y)

  return(data.frame(plugIn = mean(muYplus - muY)/mean(muAplus - muA),
                    plugInDub = mean(muYplus - muYmin)/mean(muAplus - muAmin),
                    IF = mean(phi_y)/mean(phi_a),
                    IFdub = psihat,
                    Empirical = mean(t[which(aplus>a)]),
                    EmpiricalDub = mean(yplus[which(aplus>amin)] - ymin[which(aplus>amin)])
                    )
         )
}

k = 1000
out = data.frame(plugIn = rep(NA,k),plugInDub = rep(NA,k), IF = rep(NA,k), IFdub = rep(NA,k), Empirical = rep(NA,k), EmpiricalDub = rep(NA,k))
for(i in 1:k){
  cat(i)
  out[i,] = makeSim()
}

diff1 = abs(out$plugIn - out$Empirical)
diff2 = abs(out$IF - out$Empirical)
diff3 = abs(out$plugInDub - out$EmpiricalDub)
diff4 = abs(out$IFdub - out$EmpiricalDub)
mean(diff1); mean(diff2)
mean(diff3); mean(diff4)


# # calculate the true mean E(T|Z+delta>=T>Z) from truncated normal
# Alpha = (z - muT); Beta = z+delta - muT; Z = pnorm(Beta) - pnorm(Alpha)
# truncMean = muT + (dnorm(Alpha)-dnorm(Beta))/Z; truncMean[is.infinite(truncMean)]<-NA
# Calc = mean(truncMean,na.rm = T)

