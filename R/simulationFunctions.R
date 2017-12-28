#' Make simulations
#'
#' @description Simulate various datasets for testing
#'
#' @param N dataset size
#' @param psi causal effect for some
#' @param t threshold for some
#' @param delta some simulations allow you to specify delta to get true effect and counterfactuals
#'
#' @return dataset

#### original ed sim ####
makeEdsim <- function(N){
  require(truncnorm)
  sig2y0 = 1
  y0 = rnorm(N,0,sig2y0); meanx = c(0,0,0,0)
  alpha = matrix(c(1,1,-1,-1),nrow = 4); beta = matrix(c(1,-1,-1,1)); psi =1
  delta = 5

  x = t(matrix(unlist(lapply(meanx,function(x) rnorm(N,x,1))),nrow = N,byrow =T))
  z = rtruncnorm(N,a=-2,b=2,mean=1.5*sign(t(alpha)%*%x),sd = 2)
  t = rnorm(N,t(beta)%*%x+y0,1)

  a = as.numeric(z>=t)
  aplus = as.numeric(z+delta>=t)
  amin = as.numeric(z-delta>=t)

  y = y0 + a*psi*t
  yplus = y0 + aplus*psi*t
  ymin = y0 + amin*psi*t

  true.single = mean((yplus - y)[which(aplus>a)])
  true.double = mean((yplus - ymin)[which(aplus>amin)])

  return(list(obs.dat = as.data.frame(cbind(y,a,z,t,t(x))),
              unobs.dat = as.data.frame(cbind(y,yplus,ymin,a,aplus,amin,t,y0)),
              true.eff = c(true.single, true.double)
              )
         )
}


#### based on Ed's paper - no truncation but tlatent t ####
makeSim <- function(N = 10000, psi = 1, delta = 1){
  y0 = rnorm(N,0,1); meanx = c(0,0,0,0)
  alpha = matrix(c(1,1,-1,-1),nrow = 4); beta = matrix(c(1,-1,-1,1))

  x = t(matrix(unlist(lapply(meanx, function(x) rnorm(N,x,1))), nrow = N, byrow =T))
  z = rnorm(N, 1.5*sign(t(alpha)%*%x), sd = 2)
  t = rnorm(N, t(beta)%*%x+y0, sd = 1)
  a = as.numeric(z >= t)
  y = y0 + a*psi*t
  aplus = as.numeric(z+delta >= t); amin = as.numeric(z-delta >= t)
  yplus = y0 + aplus*psi*t; ymin = y0 + amin*psi*t
  dat = as.data.frame(cbind(y,a,z,t,t(x),y0,yplus,ymin,aplus,amin))

  return(dat)
}

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
  #plugInDub = mean(muYplus - muYmin)/mean(muAplus - muAmin)

  #### empirical ####
  #Empirical = mean(t[which(aplus>a)])
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
    CACE(y=y,a=a,z=z,cov=data[,5:8],delta = delta,ranger = T,type = 'simple', split = F)
  }, error = function(e){
    print(e)
    return(NA)
  }
  )

  #CACEdub = CACE(y=y,a=a,z=z,cov=data[,5:8],delta = delta,ranger = T,type = 'double', split = F)

  #### Parametric plug in estimates ####
  plugInPar = tryCatch({
    CACE(y=y,a=a,z=z,cov=data[,5:8],delta = delta,ranger = T,type = 'simplePlugin', split = F)
  }, error = function(e){
    print(e)
    return(NA)
  }
  )
  #plugInParDub = CACE(y=y,a=a,z=z,cov=data[,5:8],delta = delta,ranger = T,type = 'pluginDouble', split = F)

  df = data.frame(if.bias = IF - true.single,
                  pi.bias = plugIn - true.single,
                  CACE.bias = CACE$phi - true.single,
                  piPar.bias = plugInPar$phi - true.single,
                  #if2.bias = IFdub - true.double,
                  #pi2.bias = plugInDub - true.double,
                  #CACE2.bias = CACEdub$phi - true.double,
                  #piPar2.bias = plugInParDub$phi - true.double,
                  delta = delta)
  return(df)

}

EstimateTrue <- function(delta, psi = 1, N = 1e6){
  dat <- makeSim(N)
  aplus = as.numeric( dat$a + delta >= dat$t )
  amin = as.numeric( dat$a - delta >= dat$t )
  yplus = dat$y0 + aplus*psi*dat$t
  ymin = dat$y0 + amin*psi*dat$t
  effect.single = mean( yplus[aplus > dat$a] - dat$y[aplus > dat$a])
  effect.double = mean( yplus[aplus > amin] - ymin[aplus > amin])
  return(c(effect.single,effect.double))
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

  aplus = as.numeric( a + delta >=t ); amin = as.numeric( a - delta >=t )
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

#### very simple ####
newSim <- function(N = 10000, psi = 1, t = 1){
  meanx = c(0,0,0,0)
  x = t(matrix(unlist(lapply(meanx, function(x) rnorm(N,x,1))), nrow = N, byrow =T))
  alpha = matrix(c(1,1,-1,-1),nrow = 4)
  z = rnorm(N, 1.5*sign(t(alpha)%*%x), sd = 2)
  a = as.numeric(z >= t)
  y = a*psi + rnorm(N)
  dat = as.data.frame(cbind(y,a,z,rep(t,N),t(x)))
  names(dat) <- c('y','a','z','t','x1','x2','x3','x4')
  return(dat)
}

EstimateNewSim <- function(delta, psi = 1, N = 1e6){
  dat <- newSim(N)
  aplus = as.numeric( dat$a + delta >= dat$t )
  amin = as.numeric( dat$a - delta >= dat$t )
  yplus = aplus*psi
  ymin = amin*psi
  effect.single = mean( yplus[aplus > dat$a] - dat$y[aplus > dat$a])
  effect.double = mean( yplus[aplus > amin] - ymin[aplus > amin])
  return(c(effect.single,effect.double))
}

#### simplified version of ed's (no truncation, no latent t) ####
simFunc <- function(N=5000,delta = 1, psi = 2){
  y0 = rnorm(N,0,1); meanx = c(0,0,0,0)
  alpha = matrix(c(1,1,-1,-1),nrow = 4)
  x = matrix(unlist(lapply(meanx, function(x) rnorm(N,x,1))), nrow = N, byrow =T)
  z = rnorm(N, t(alpha)%*%t(x), sd = 2)
  a = as.numeric(z >= y0)
  y = y0 + a*psi

  true.eff = psi
  true.ymean = psi*pnorm(z)
  true.amean = pnorm(z)
  true.ymean.plus = psi*pnorm(z+delta)
  true.amean.plus = pnorm(z+delta)
  true.ymean.min = psi*pnorm(z-delta)
  true.amean.min = pnorm(z-delta)
  true.z <- dnorm(z, mean = t(alpha)%*%t(x), sd = 2)
  true.z.min <- dnorm(z-delta, mean = t(alpha)%*%t(x), sd = 2)
  true.z.plus <- dnorm(z+delta, mean = t(alpha)%*%t(x), sd = 2)

  return(data.frame(y,a,z,x,true.ymean,true.amean,true.ymean.plus,true.ymean.min,
                    true.amean.plus,true.amean.min,
                    true.z, true.z.min, true.z.plus))
}

#### like simplified ed's but with ks style misspecification of covariates ####
gammaSim <- function(psi = 1, delta = 10, N = 5000, dep.x = FALSE){
  x = matrix(unlist(lapply(c(0,0,0,0), function(x) rnorm(N,x,1))), nrow = N, byrow =T)
  z = rgamma(N, rate = min(15,exp(x[,4])+abs(x[,1])+expit(x[,2])), shape = 10+abs(x[,3]*x[,4]))
  beta.x = abs(x[,4])+log(abs(1-x[,1]/x[,2])) + exp(x[,3]-x[,4])
  y0 = rnorm(N,0,1)
  a = as.numeric(z >= 75*y0 + dep.x*(beta.x))
  aplus = as.numeric(z + delta >= 75*y0 + dep.x*(beta.x))
  amin = as.numeric(z - delta >= 75*y0 + dep.x*(beta.x))
  y = y0 + a*psi
  yplus = y0 + aplus*psi
  ymin = y0 + amin*psi

  true.single = mean((yplus - y)[which(aplus > a)])
  true.double = mean((yplus - ymin)[which(aplus > amin)])

  return(list(data = data.frame(y=y,a=a,z=z,x=x,y0=y0,yplus=yplus,ymin=ymin,aplus = aplus, amin = amin),
              true = data.frame(true.single, true.double)))
}

#### similar to gammasim but effect depends on delta ####
# need the compliers to be different
gammaSimD <- function(psi = 1, delta = 10, N = 5000, dep.x = TRUE){
  x = matrix(unlist(lapply(c(0,0,0,0), function(x) rnorm(N,x,1))), nrow = N, byrow =T)
  z = rgamma(N, rate = min(15,exp(x[,4])+abs(x[,1])+expit(x[,2])), shape = 10+abs(x[,3]*x[,4]))
  beta.x = abs(x[,4])+log(abs(1-x[,1]/x[,2])) + exp(x[,3]-x[,4])
  y0 = rnorm(N,0,1)
  a = as.numeric(z/50 >= y0 + dep.x*(beta.x))
  aplus = as.numeric(z/50 + delta >= y0 + dep.x*(beta.x))
  amin = as.numeric(z/50 - delta >= y0 + dep.x*(beta.x))
  y = y0 + a*psi + dep.x*(beta.x)
  yplus = y0 + aplus*psi + dep.x*(beta.x)
  ymin = y0 + amin*psi + dep.x*(beta.x)

  true.single = mean((yplus - y)[which(aplus > a)])
  true.double = mean((yplus - ymin)[which(aplus > amin)])

  return(list(data = data.frame(y=y,a=a,z=z,x=x,y0=y0,yplus=yplus,ymin=ymin,aplus = aplus, amin = amin),
              true = data.frame(true.single, true.double)))
}
