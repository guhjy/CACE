# basic test
require(ggplot2)
require(plyr)

#### fixed T, no endogeneity ####
k = 50
m = 5
deltas = seq(1,5, length = m)
output = NULL
algs = list(y.est = 'ranger', a.est = 'ranger', z.est = 'flexcode')

for(d in deltas){
  for(i in 1:k){
    N = 5000
    y0 = rnorm(N,0,1); meanx = c(0,0,0,0)
    alpha = matrix(c(1,1,-1,-1),nrow = 4); beta = matrix(c(1,-1,-1,1))
    x = matrix(unlist(lapply(meanx, function(x) rnorm(N,x,1))), nrow = N, byrow =T)
    z = rnorm(N, 1.5*sign(t(alpha)%*%t(x)), sd = 2)
    t = 1; psi = 1
    a = as.numeric(z >= t)
    y = y0 + a*psi*t
    true.eff = t*psi
    test.pi <- single.shift.pi(y = y,a = a,z=z,x = x, delta = d, algo = algs)
    test.if <- single.shift(y = y,a = a,z=z,x = x, delta = d, algo = algs)
    output = rbind(output,c(true.eff, d, test.if$psi, test.pi$psi,
                            mean(as.numeric(z+d >=t) - as.numeric(z >= t)),
                            mean(as.numeric(z+d >=t)*psi*t - as.numeric(z >= t)*psi*t)))
    print(paste(d,":",i))
  }
}

output <- as.data.frame(output)
names(output) <- c('true.eff', 'delta', 'IF', 'PI', 'denom', 'numerator')
mean(output$IF[which(output$delta==1)])
mean(output$IF[which(output$delta==5)])
mean(output$PI[which(output$delta==1)])
mean(output$PI[which(output$delta==5)])


results = data.frame(results = c(output[,3], output[,4]),
                     delta = rep(output[,2],2),
                     type = c(rep('IF', dim(output)[1]), rep('PI', dim(output)[1])))
summed <- ddply(results, .(type, delta), summarize, mean = mean(results), sd = sd(results))
summed$lower <- summed$mean - summed$sd
summed$upper <- summed$mean + summed$sd

dodge <- position_dodge(width=0.2)
ggplot(summed) + geom_point(aes(x = delta, y = mean, shape = type), position = dodge) +
  geom_errorbar(aes(x = delta, ymin = lower, ymax = upper, group = type), width = .1, position = dodge) +
  geom_hline(yintercept = 1)

require(xtable)
x <- xtable(summed[,1:4])
print.xtable(x, type='latex', file = 'NPEstTable.tex')

#### actual endogeneity ####
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

N = 5000; delta = 1; psi = 2
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

sim <- simFunc(); attach(sim)
dat <- as.data.frame(cbind(z,x))
dat.plus <- as.data.frame(cbind(z+delta,x))
names(dat.plus) <- names(dat)

estFunc <- function(data, delta = 1){
  dat = data[,3:7]
  dat.plus = cbind(data[,3]+delta, data[,4:7])
  names(dat.plus) <- names(dat)

  y.est = predict(y.mean.est(data$y, dat, 'glm'), dat, type = 'response')
  y.plus.est = predict(y.mean.est(data$y, dat, 'glm'), dat.plus, type = 'response')

  a.est = predict(a.mean.est(data$a, dat, 'glm'), dat, type = 'response')
  a.plus.est = predict(a.mean.est(data$a, dat, 'glm'), dat.plus, type = 'response')

  z.mean <- predict(glm(z~x, family = 'gaussian'))
  z.var <- mean( (z - z.mean)^2  )
  phat.gk <- sapply(data$z, function(y) (1/N)*sum(gK(sqrt( ((y - z.mean))^2/z.var ) )))
  phat.gk.min <- sapply((data$z-delta), function(y) (1/N)*sum(gK(sqrt( ((y - z.mean))^2/z.var ) )))

  return(data.frame(data$y,data$a,data$z,y.est, y.plus.est, a.est, a.plus.est, phat.gk, phat.gk.min))
}

# ymean: not so good
y.est = predict(y.mean.est(y, dat, 'glm'), dat, type = 'response')
y.plus.est = predict(y.mean.est(y, dat, 'glm'), dat.plus, type = 'response')
plot(y.est ~ true.ymean)
plot(y.plus.est ~ true.ymean.plus)

# amean: very good
a.est = predict(a.mean.est(a, dat, 'glm'), dat, type = 'response')
a.plus.est = predict(a.mean.est(a, dat, 'glm'), dat.plus, type = 'response')
plot(a.est ~ true.amean)
plot(a.plus.est ~ true.amean.plus)

# z prob: glm for means, gaussian kernel for density
# does poorly
z.mean <- predict(glm(z~x, family = 'gaussian'))
z.var <- mean( (z - z.mean)^2  )

gK <- function(x){(1/sqrt(2*pi))*exp(-(x^2)/2)}

phat.gk <- sapply(z, function(y) (1/N)*sum((1/sqrt(z.var))*gK( sqrt( ((y - z.mean))^2/z.var ) ) ) )
phat.gk.min <- sapply((z-delta), function(y) (1/N)*sum(gK(sqrt( ((y - z.mean))^2/z.var ) )))
phat.ks = density( (z - z.mean)/sqrt(z.var)  )$y

plot(phat.ks ~ true.z); abline(0,1,col = 'red')
plot((phat.gk)/sqrt(z.var) ~ true.z); abline(0,1,col = 'red')

plot(phat.gk.min ~ true.z.min)

#boxcar (ignoring glm): terrible
K<- function(m,l,h){
  sapply(apply(sqrt((m-l)^2),1,sum),function(x) as.numeric(x<=h/2))
  #sapply(sum(sqrt((m-l)^2)),function(x) as.numeric(x<=h/2))
}
phat.boxcar <- sapply(z, function(y) (1/N)*sum(K(x,y,3)/3))
plot(phat.boxcar~true.z)

#### Properly specified ests ####
# PI does perfect; IF does well
phi.PI = mean(true.ymean.plus - true.ymean)/mean(true.amean.plus - true.amean)
phi.IF = mean((true.z.min/true.z)*(y - true.ymean) - (y - true.ymean.plus))/
  mean((true.z.min/true.z)*(a - true.amean) - (a - true.amean.plus))

output = matrix(rep(NA,200), ncol = 2)
for(i in 1:100){
  dat = simFunc()
  phi.PI = mean(dat$true.ymean.plus - dat$true.ymean)/
    mean(dat$true.amean.plus - dat$true.amean)
  phi.IF = mean((dat$true.z.min/dat$true.z)*(dat$y - dat$true.ymean) - (dat$y - dat$true.ymean.plus))/
    mean((dat$true.z.min/dat$true.z)*(dat$a - dat$true.amean) - (dat$a - dat$true.amean.plus))
  output[i,] <- c(phi.PI, phi.IF)
}
apply(output, 2, mean); apply(output, 2, sd)

k = 100
m = 5
deltas = seq(1,5, length = m)
output = NULL
for(d in deltas){
  for(i in 1:k){
    true.eff = 2
    dat = simFunc()
    test.pi <- mean(dat$true.ymean.plus - dat$true.ymean)/
      mean(dat$true.amean.plus - dat$true.amean)
    test.if <- mean((dat$true.z.min/dat$true.z)*(dat$y - dat$true.ymean) - (dat$y - dat$true.ymean.plus))/
      mean((dat$true.z.min/dat$true.z)*(dat$a - dat$true.amean) - (dat$a - dat$true.amean.plus))
    output = rbind(output,c(true.eff, d, test.if, test.pi))
    print(paste(d,":",i))
  }
}

output <- as.data.frame(output)
names(output) <- c('true.eff', 'delta', 'IF', 'PI')
write.csv(output, 'correctlySpecified.csv')

results = data.frame(results = c(output[,3], output[,4]),
                     delta = rep(output[,2],2),
                     type = c(rep('IF', dim(output)[1]), rep('PI', dim(output)[1])))
summed <- ddply(results, .(type, delta), summarize, mean = mean(results), sd = sd(results))
summed$lower <- summed$mean - 2*summed$sd
summed$upper <- summed$mean + 2*summed$sd

dodge <- position_dodge(width=0.2)
ggplot(summed) + geom_point(aes(x = delta, y = mean, shape = type), position = dodge) +
  geom_errorbar(aes(x = delta, ymin = lower, ymax = upper, group = type), width = .1, position = dodge) +
  geom_hline(yintercept = true.eff, col = 'red')+
  ggtitle('Correctly specified IF vs PI')+
  xlab('Delta')+ylab('Estimate')
ggsave('correctlySpecified.png', height = 4, width = 6)

#### With parametrically estimated means/prop scores ####
# both do well, IF maybe a little better
phi.PI.est = mean(y.plus.est - y.est)/mean(a.plus.est - a.est)
phi.IF.est = mean((phat.gk.min/phat.gk)*(y - y.est) - (y - y.plus.est))/
  mean((phat.gk.min/phat.gk)*(a - a.est) - (a - a.plus.est))

output = matrix(rep(NA,200), ncol = 2)
for(i in 1:100){
  dat = estFunc(simFunc())
  phi.PI = mean(dat$y.plus.est - dat$y.est)/mean(dat$a.plus.est - dat$a.est)
  phi.IF = mean((dat$phat.gk.min/dat$phat.gk)*(dat$data.y - dat$y.est) - (dat$data.y - dat$y.plus.est))/
    mean((dat$phat.gk.min/dat$phat.gk)*(dat$data.a - dat$a.est) - (dat$data.a - dat$a.plus.est))
  output[i,] <- c(phi.PI, phi.IF)
}
apply(output, 2, mean); apply(output, 2, sd)

#### do it a bunch of times ####
k = 10
m = 5
deltas = seq(1,5, length = m)
output = NULL
algs = list(y.est = 'glm', a.est = 'glm', z.est = 'glm')

for(d in deltas){
  for(i in 1:k){
    data = simFunc()
    test.pi <- single.shift.pi(y = data$y,a = data$a,z=data$z,x = data[,4:7]
                               , delta = d, algo = algs)
    test.if <- single.shift(y = data$y,a = data$a,z=data$z,x = data[,4:7]
                            , delta = d, algo = algs)
    output = rbind(output,c(true.eff, d, test.if$psi, test.pi$psi, test.if$sd, test.pi$sd))
    print(paste(d,":",i))
  }
}

output <- as.data.frame(output)
names(output) <- c('true.eff', 'delta', 'IF', 'PI', 'IF.sd', 'PI.sd')
write.csv(output, 'ParametricEsts.csv')

results = data.frame(results = c(output[,3], output[,4]),
                     delta = rep(output[,2],2),
                     type = c(rep('IF', dim(output)[1]), rep('PI', dim(output)[1])))
summed <- ddply(results, .(type, delta), summarize, mean = mean(results), sd = sd(results))
summed$lower <- summed$mean - 2*summed$sd
summed$upper <- summed$mean + 2*summed$sd


dodge <- position_dodge(width=0.2)
ggplot(summed) + geom_point(aes(x = delta, y = mean, shape = type), position = dodge) +
  geom_errorbar(aes(x = delta, ymin = lower, ymax = upper, group = type), width = .1, position = dodge) +
  geom_hline(yintercept = true.eff, col = 'red')+
  ggtitle('Incorrectly specified parametric estimates, IF vs PI')+
  xlab('Delta')+ylab('Estimate')
ggsave('ParametricEsts.png', height = 4, width = 6)

# use ranger
k = 100
m = 5
N = 5000
true.eff = 2
deltas = seq(1,5, length = m)
output = NULL
algs = list(y.est = 'ranger', a.est = 'ranger', z.est = 'glm')

for(d in deltas){
  for(i in 1:k){
    data = simFunc()
    test.pi <- single.shift.pi(y = data$y,a = data$a,z=data$z,x = data[,4:7]
                               , delta = d, algo = algs)
    test.if <- single.shift(y = data$y,a = data$a,z=data$z,x = data[,4:7]
                            , delta = d, algo = algs)
    output = rbind(output,c(true.eff, d, test.if$psi, test.pi$psi, test.if$sd, test.pi$sd))
    print(paste(d,":",i))
  }
}

output <- as.data.frame(output)
names(output) <- c('true.eff', 'delta', 'IF', 'PI', 'IF.sd', 'PI.sd')
write.csv(output, 'RangerRangerGLMEsts.csv')

results = data.frame(results = c(output[,3], output[,4]),
                     delta = rep(output[,2],2),
                     type = c(rep('IF', dim(output)[1]), rep('PI', dim(output)[1])))
summed <- ddply(results, .(type, delta), summarize, mean = mean(results), sd = sd(results))
summed$lower <- summed$mean - 2*summed$sd
summed$upper <- summed$mean + 2*summed$sd


dodge <- position_dodge(width=0.2)
ggplot(summed) + geom_point(aes(x = delta, y = mean, shape = type), position = dodge) +
  geom_errorbar(aes(x = delta, ymin = lower, ymax = upper, group = type), width = .1, position = dodge) +
  geom_hline(yintercept = true.eff, col = 'red')+
  ggtitle('Incorrectly specified Ranger estimates, IF vs PI')+
  xlab('Delta')+ylab('Estimate')
ggsave('RangerRangerGLMEsts.png', height = 4, width = 6)


#### hide variables a la k+s ####
simFunc.alt <- function(N=5000,delta = 1, psi = 2){
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

  x1 <- cbind( exp(x[,1]/2) , 10+x[,2]/(1+exp(x[,1])) , (x[,1]*x[,3]/25 + 0.6)^3 , (x[,2]+x[,4]+20)^2)

  return(data.frame(y,a,z,x,x1,true.ymean,true.amean,true.ymean.plus,true.amean.plus,
                    true.z, true.z.min))
}



#### double shift ####
dat = simFunc()
dub.if <- double.shift(y = dat$y, a = dat$a, z = dat$z, x = dat[,4:7], delta = 1)
