#### make a harder simulation ####
rm( list = ls() )
psi = 1
delta = 1
N = 5000
x = matrix(unlist(lapply(c(0,0,0,0), function(x) rnorm(N,x,1))), nrow = N, byrow =T)
z = rgamma(N, rate = 10*exp(x[,4])*abs(x[,1]/x[,2]), shape = 2e-2*abs(x[,3]*x[,4]))
y0 = rnorm(N,0,1)
a = as.numeric(z >= 2*y0)
aplus = as.numeric(z + delta >= 2*y0)
amin = as.numeric(z - delta >= 2*y0)
y = y0 + a*psi
yplus = y0 + aplus*psi
ymin = y0 + amin*psi

true.single = mean((yplus - y)[which(aplus > a)])
true.double = mean((yplus - ymin)[which(aplus > amin)])

algs = list(y.est = 'ranger', a.est = 'ranger', z.est = 'flexcode')

test.pi <- single.shift.pi(y = y, a = a, z = z, x = x, delta = delta, algo = algs)
test.if <- single.shift(y = y, a = a, z = z, x = x, delta = delta, algo = algs)

# use ranger and flexcode
k = 10
m = 3
N = 5000
true.eff = 2
deltas = seq(1,5, length = m)
output = NULL
algs = list(y.est = 'ranger', a.est = 'ranger', z.est = 'flexcode')

for(d in deltas){
  for(i in 1:k){
    x = matrix(unlist(lapply(c(0,0,0,0), function(x) rnorm(N,x,1))), nrow = N, byrow =T)
    z = rgamma(N, rate = 10*exp(x[,4])*abs(x[,1]/x[,2]), shape = 2e-2*abs(x[,3]*x[,4]))
    y0 = rnorm(N,0,1)
    a = as.numeric(z >= 2*y0)
    aplus = as.numeric(z + delta >= 2*y0)
    amin = as.numeric(z - delta >= 2*y0)
    y = y0 + a*psi
    yplus = y0 + aplus*psi
    ymin = y0 + amin*psi

    true.eff = mean((yplus - y)[which(aplus > a)])

    test.pi <- single.shift.pi(y = y, a = a, z = z, x = x, delta = delta, algo = algs, nfolds = 1)
    test.if <- single.shift(y = y, a = a, z = z, x = x, delta = delta, algo = algs, nfolds = 1)

    output = rbind(output,c(true.eff, d, test.if$psi, test.pi$psi, test.if$sd, test.pi$sd))
    print(paste(d,":",i))
  }
}

output <- as.data.frame(output)
write.csv(output, 'simsHardFCRanger.csv', row.names = F)
names(output) <- c('true.eff', 'delta', 'IF', 'PI', 'IF.sd', 'PI.sd')

library(plyr); library(ggplot2)
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
  labs(x = 'delta',y = 'Estimates')+
  ggtitle(paste('Estimates at',k, 'simulations (Ranger and FlexCoDE)'))

ggsave('FCRangerEsts.png')
require(xtable)
x <- xtable(summed[,1:4])
print.xtable(x, type='latex', file = 'FCRangerEsts.tex')
