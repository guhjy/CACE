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
N = 5000
y0 = rnorm(N,0,1); meanx = c(0,0,0,0)
alpha = matrix(c(1,1,-1,-1),nrow = 4)
x = matrix(unlist(lapply(meanx, function(x) rnorm(N,x,1))), nrow = N, byrow =T)
z = rnorm(N, t(alpha)%*%t(x), sd = 2)
psi = 1
a = as.numeric(z >= y0)
y = y0 + a*psi
true.eff = psi
true.ymean = psi*pnorm(z)
true.amean = pnorm(z)

dat <- as.data.frame(cbind(z,x))

# not so good
y.est = predict(y.mean.est(y, dat, 'glm'), dat, type = 'response')
plot(y.est ~ true.ymean)

# very good
a.est = predict(a.mean.est(a, dat, 'glm'), dat, type = 'response')
plot(a.est ~ true.amean)

true.z <- dnorm(z, mean = t(alpha)%*%t(x), sd = 2)
z.reg <- glm(z~x, family = 'gaussian')
z.mean <- predict(glm(z~x, family = 'gaussian'))
z.var <- mean( (z - z.mean)^2  )
p.z <- ksmooth( x = z, y = (z - z.mean)/(sqrt(z.var))  )
plot(p.z$y ~ z)

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
    z = rnorm(N, t(alpha)%*%t(x), sd = 2)
    psi = 1
    a = as.numeric(z >= y0)
    y = y0 + a*psi
    true.eff = psi
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
