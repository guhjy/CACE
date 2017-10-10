library(doParallel)

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
  true.z <- dnorm(z, mean = t(alpha)%*%t(x), sd = 2)
  true.z.min <- dnorm(z-delta, mean = t(alpha)%*%t(x), sd = 2)

  return(data.frame(y,a,z,x,true.ymean,true.amean,true.ymean.plus,true.amean.plus,
                    true.z, true.z.min))
}


k = 2
m = 2
deltas = seq(1,5, length = m)
output = NULL
algs = list(y.est = 'random forest', a.est = 'random forest', z.est = 'glm')

foreach (j=1:length(deltas), .options.multicore=list(preschedule=TRUE)) %dopar% {
  for(i in 1:k){
    data = simFunc()
    test.pi <- single.shift.pi(y = data$y,a = data$a,z=data$z,x = data[,4:7], delta = deltas[j], algo = algs)
    test.if <- single.shift(y = data$y,a = data$a,z=data$z,x = data[,4:7], delta = deltas[j], algo = algs)
    output = rbind(output,c(true.eff, deltas[j], test.if$psi, test.pi$psi, test.if$sd, test.pi$sd))
    print(paste(deltas[j],":",i))
  }
}
