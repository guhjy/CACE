#### This is a very simplified simulation based on Ed's
#### Ed's simulation led to positivity violations for the shift estimator
#### This was caused by the truncation of the z values
#### This simulation still not working: plugin outperforming IF

makeSimSimple<-function(){
  rm(list = ls())
  today = format(Sys.Date(), "%Y%m%d")
  sig2y0 = 1;N = 10000
  y0 = rnorm(N,0,sig2y0)
  delta = 2
  psi = 1

  #sign <- sample(c(-1,1), N, replace = T)

  z <- rnorm(N,1.5,sd = 2)
  t <- rnorm(N,y0-.2,1)
  a <- as.numeric(z>=t); aplus <- as.numeric(z+delta>=t); amin <- as.numeric(z-delta>=t)
  y = y0 + a*psi*t; yplus = y0 + aplus*psi*t; ymin = y0 + amin*psi*t

  #### true parameters ####
  muA = pnorm(z+.2);muAplus = pnorm(((z+delta)+.2));muAmin = pnorm(((z-delta)+.2))
  muY = -.2*pnorm(z+.2) - dnorm(z + .2)#*sqrt(2)
  muYplus = -.2*pnorm(z+delta+.2) - dnorm(z + delta + .2)#*sqrt(2)
  muYmin = -.2*pnorm(z-delta+.2) - dnorm(z - delta + .2)#*sqrt(2)
  pi = dnorm(z,mean = 1.5, sd = 2); pi_min = dnorm(z-delta,mean = 1.5, sd = 2); pi_plus = dnorm(z+delta,mean = 1.5, sd = 2)
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

  #### algorithm estimator ####
  phiSingle = CACE(y=y,a=a,z=z,cov = rnorm(length(y)),delta = delta,ranger = T,type = 'single',split = F)
  phiDub = CACE(y=y,a=a,z=z,cov=rnorm(length(y)),delta = delta,ranger = T,type = 'double',split = F)

  return(data.frame(plugIn = mean(muYplus - muY)/mean(muAplus - muA),
                    plugInDub = mean(muYplus - muYmin)/mean(muAplus - muAmin),
                    IF = mean(phi_y)/mean(phi_a),
                    IFdub = psihat,
                    IFest = phiSingle$phi,
                    IFestdub = phiDub$phi,
                    Empirical = mean(t[which(aplus>a)]),
                    EmpiricalDub = mean(yplus[which(aplus>amin)] - ymin[which(aplus>amin)])
  )
  )
}

k = 100
out = NA
for(i in 1:k){
  print(i)
  temp = makeSimSimple()
  print(temp)
  out<-rbind(out,temp)
}

diff1 = abs(out$plugIn - out$Empirical)/out$Empirical
diff2 = abs(out$IF - out$Empirical)/out$Empirical
diff3 = abs(out$plugInDub - out$EmpiricalDub)/out$EmpiricalDub
diff4 = abs(out$IFdub - out$EmpiricalDub)/out$EmpiricalDub
mean(diff1); mean(diff2)
mean(diff3); mean(diff4)

par(mfrow = c(2,2))
hist(diff1); hist(diff2);hist(diff3);hist(diff4)
par(mfrow = c(1,1))
