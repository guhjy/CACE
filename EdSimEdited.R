#### ed sim ####
# trying to fix truncated distribution issue w Ed's code
# at this point (8/31), the plug in outperforms the IF-based estimators
# at this point (9/7), I realized that it's fine if the correctly specified plug in is more efficient
# need to figure out what comparisons to do, have as options:
#   - true plug in
#   - estimated plug in (how to estimate?)
#   - true IF
#   - estimated IF

makeSim <- function(N = 10000, delta = 2, psi = 1){
  today = format(Sys.Date(), "%Y%m%d")
  y0 = rnorm(N,0,1); meanx = c(0,0,0,0)
  alpha = matrix(c(1,1,-1,-1),nrow = 4); beta = matrix(c(1,-1,-1,1))

  # observed data
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

  #### CACE estimates ####
  CACE = CACE(y=y,a=a,z=z,cov=data[,5:8],delta = delta,ranger = T,type = 'single')
  CACEdub = CACE(y=y,a=a,z=z,cov=data[,5:8],delta = delta,ranger = T,type = 'double')

  #### Parametric plug in estimates ####
  plugInPar = CACE(y=y,a=a,z=z,cov=data[,5:8],delta = delta,ranger = T,type = 'plugin')
  plugInParDub = CACE(y=y,a=a,z=z,cov=data[,5:8],delta = delta,ranger = T,type = 'pluginDouble')

  return(data.frame(if.bias = IF - Empirical,
                    pi.bias = plugIn - Empirical,
                    CACE.bias = CACE$phi - Empirical,
                    piPar.bias = plugInPar$phi - Empirical,
                    if2.bias = IFdub - EmpiricalDub,
                    pi2.bias = plugInDub - EmpiricalDub,
                    CACE2.bias = CACEdub$phi - EmpiricalDub,
                    piPar2.bias = plugInParDub$phi - EmpiricalDub,
                    delta = delta
                    )
  )

}

delta.range <- seq(0 , 1, by = 1)
temp.function <-function(x){delta = x; data = makeSim(delta = delta, N = 10000); getEstimatorBias(data, delta = delta)}
ests <- lapply(delta.range, function(x) temp.function(x))


#### think i want to delete all of this ####
get.N.ests <- function(N, delta){
  k = N - 1
  out = makeSim(delta = delta)
  for(i in 1:k){
    out = rbind(out,makeSim())
  }
  out$PI.bias = out$plugIn - out$Empirical
  out$IF.bias = out$IF - out$Empirical
  out$PI2.bias = out$plugInDub - out$EmpiricalDub
  out$IF2.bias = out$IFdub - out$EmpiricalDub

  means = apply(out, 2, mean)
  sds = apply(out, 2, sd)
  return(cbind(means, sds))
}

big.out <- lapply(delta.range, function(x) get.N.ests(N = 500, delta = x))
pi.bias <- as.data.frame(cbind(delta.range,matrix(unlist(lapply(big.out, function(x) x['PI.bias',])),ncol = 2, byrow = T)))
if.bias <- as.data.frame(cbind(delta.range,matrix(unlist(lapply(big.out, function(x) x['IF.bias',])),ncol = 2, byrow = T)))
names(pi.bias) <- c('delta','PlugInBias', 'PlugInSD')
names(if.bias) <- c('delta','IFBias', 'IFSD')
single.set <- as.data.frame(cbind(pi.bias,if.bias))[,-4]

# using ggplot (not done, legend being dumb)
theme_set(theme_bw(base_size = 10))
p <- ggplot() +
  geom_point(data = if.bias, aes(x = delta, y = IFBias, ymin = IFBias - 1.96*IFSD, ymax= IFBias + 1.96*IFSD, colour = 'red'))+
  geom_pointrange(data = if.bias, aes(x = delta, y = IFBias, ymin = IFBias - 1.96*IFSD, ymax= IFBias + 1.96*IFSD, colour = 'red')) +
  geom_point(data = pi.bias, aes(x = delta, y = PlugInBias), colour = 'blue') +
  geom_pointrange(data = pi.bias, aes(x = delta, y = PlugInBias, ymin = PlugInBias - 1.96*PlugInSD, ymax= PlugInBias + 1.96*PlugInSD), colour = 'blue') +
  geom_hline( yintercept = 0 )+
  ylab("Bias") +
  xlab("Shift")

min.val = min(single.set$PlugInBias-1.96*single.set$PlugInSD,single.set$IFBias-1.96*single.set$IFSD)
max.val = max(single.set$PlugInBias+1.96*single.set$PlugInSD,single.set$IFBias+1.96*single.set$IFSD)
plot(PlugInBias~delta, data = single.set, ylab = "Bias", xlab = "Shift", ylim = c(min.val, max.val), pch = 19)
points(IFBias ~ delta, data = single.set, col = 'red', pch = 19)
abline(h = 0, lty = 2)
legend('topright',legend = c('Plug in', 'IF'), col = c('black', 'red'), pch = c(19,19), cex = .75)


