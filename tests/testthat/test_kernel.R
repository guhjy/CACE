### test the kernel ###
gK <- function(x){(1/sqrt(2*pi))*exp(-(x^2)/2)}

test_that('kernel can work in simple case',{
  x = rnorm(1000)
  N = length(x)
  h = sd(x)
  test.seq = seq(-2,2,by=.2)
  test.pi = sapply(test.seq, function(y) (1/length(x))*sum(gK((1/h)*sqrt((y - x)^2))) )

  expect_lt(-cor(test.pi, dnorm(test.seq)),-.8)

  plot(test.pi~dnorm(test.seq),
       xlim = c(min(test.pi, dnorm(test.seq)),max(test.pi, dnorm(test.seq))),
       ylim = c(min(test.pi, dnorm(test.seq)),max(test.pi, dnorm(test.seq))))
  abline(0,1)
})

x = rnorm(1000)
z = rnorm(length(x),x,1)
true.density = dnorm(z,mean = x, sd = 1)
mean.est = predict(glm(z~x, family = 'gaussian'))
sd.est = sqrt( mean( (mean.est - z)^2 ) )
pi.est = sapply(z, function(y) (1/length(x))*sum(gK((1/sd.est)*sqrt((y - x)^2))) )/sd.est
pi.est = sapply(z, function(y) (1/length(x))*sum( gK( y - mean.est ) ))/sd.est
pi.est = sapply(z, function(y) (1/length(x))*sum( gK( (y - mean.est)/sd.est ) ))/sd.est
plot(pi.est ~ true.density)

library(CACE)
library(FlexCoDE)
fc = z.condldens.est(z,as.matrix(x),algo = 'flexcode')
pred = predict(fc, as.matrix(x))
pihat = get_probs(z, pred$z, pred$CDE)
plot(pihat~true.density, xlim = c(0,.4), ylim = c(0,.4))


library(hdrcde)
pi.est = cde(x,z)
plot(pi.est ~ true.density)
