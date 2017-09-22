#### file created 9/7/2017
#### should output simulations for the JASA paper

rm(list = ls())
require(truncnorm)
today = format(Sys.Date(), "%Y%m%d")
sig2y0 = 1;N = 10000
y0 = rnorm(N,0,sig2y0); meanx = c(0,0,0,0)
alpha = matrix(c(1,1,-1,-1),nrow = 4); beta = matrix(c(1,-1,-1,1)); psi =1
delta = .1

x = t(matrix(unlist(lapply(meanx,function(x) rnorm(N,x,1))),nrow = N,byrow =T))
muZ <- (1.5*sign(t(alpha)%*%x))[1,]; muT <- (t(beta)%*%x+y0)[1,]
z = rtruncnorm(N,a=-2,b=2,mean=muZ,sd = 2)
t = rnorm(N,muT,1)
a = as.numeric(z>=t)
y = y0 + a*psi*t

# counterfactuals
aplus = as.numeric(z+delta>=t); amin = as.numeric(z-delta>=t)
yplus = y0 + aplus*psi*t; ymin = y0 + amin*psi*t
