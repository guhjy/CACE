makeEdsim <- function(N){
  require(truncnorm)
  sig2y0 = 1
  y0 = rnorm(N,0,sig2y0); meanx = c(0,0,0,0)
  alpha = matrix(c(1,1,-1,-1),nrow = 4); beta = matrix(c(1,-1,-1,1)); psi =1
  delta = 5

  x = t(matrix(unlist(lapply(meanx,function(x) rnorm(N,x,1))),nrow = N,byrow =T))
  z = rtruncnorm(N,a=-2,b=2,mean=1.5*sign(t(alpha)%*%x),sd = 2)
  t = rnorm(N,t(beta)%*%x+y0,1)
  a = as.numeric(z>=t); aplus = as.numeric(z+delta>=t)
  y = y0 + a*psi*t; yplus = y0 + aplus*psi*t

  return(as.data.frame(cbind(y,a,z,t,t(x))))
}
