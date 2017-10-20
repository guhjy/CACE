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

x1 = rnorm(1000); x2 = rnorm(1000)
z = rnorm(length(x1),x1+x2,1)
true.density = dnorm(z,mean = x1+x2, sd = sd(z))
mean.est = predict(glm(z~x1+x2, family = 'gaussian'))
sd.est = sqrt( mean( (mean.est - z)^2 ) )
pi.est = sapply(z, function(y) (1/length(x))*sum(gK((1/sd.est)*sqrt((y - x)^2))/sd.est) )
pi.est = sapply(z, function(y) (1/length(x))*sum( gK( y - mean.est )/sd.est ))
pi.est = sapply(z, function(y) (1/length(x))*sum( gK( (y - mean.est)/sd.est )/sd.est ))
plot(pi.est ~ true.density)

library(CACE)
library(FlexCoDE)
fc = z.condldens.est(z,as.matrix(x),algo = 'flexcode')
pred = predict(fc, as.matrix(x))
pihat = get_probs(z, pred$z, pred$CDE)
plot(pihat~true.density, xlim = c(0,.4), ylim = c(0,.4))


# make xgrid finer
# make it adaptable to number of xs
# speed up
# average over nearby grid to get pihat
library(hdrcde)
Zst = as.numeric((z - mean.est)/sd.est)
xmat = cbind(x1,x2)
pi.reg = cde(xmat,Zst,nxmargin = 100)


pi.est = rep(NA,length(Zst))
for(i in 1:length(Zst)){
  x1rank = min(rank(append(x1[i],pi.reg$x$x1),ties.method = 'last')[1],dim(pi.reg$z)[1])
  x2rank = min(rank(append(x2[i],pi.reg$x$x2),ties.method = 'last')[1],dim(pi.reg$z)[1])
  yrank = min(rank(append(Zst[i],pi.reg$y),ties.method = 'last')[1],length(pi.reg$y))

  rg1 = c(max(x1rank-1, 1):x1rank)
  rg2 = c(max(x2rank-1, 1):x2rank)
  rg3 = c(max(yrank-1, 1):yrank)
  pi.est[i] <- mean(pi.reg$z[rg1,rg2,rg3])
}

plot(pi.est~true.density); abline(0,1,col = 'red')


get.xrank <- function(x,reg,nx){
  out = rep(NA,length(x))
  for(i in 1:length(x)){
    out[i] = min(rank(append(x[i],reg$x[[i]]),ties.method = 'last')[1],nx)
  }
  return(out)
}
format.range <- function(x){
  if(length(x)>1){
    out = paste(x[1],x[2], sep = ":")
  }
  else{out = x}
  return(out)
}

pi.est = rep(NA,length(Zst))
for(i in 1:length(Zst)){
  xrank = get.xrank(xmat[i,], pi.reg, dim(pi.reg$z)[1])
  yrank = min(rank(append(Zst[i],pi.reg$y),ties.method = 'last')[1],length(pi.reg$y))

  xrange = lapply(xrank, function(x) c(max(x-1,1):x))
  yrange = c(max(yrank-1, 1):yrank)
  temprange = xrange
  temprange[[length(temprange)+1]] = yrange
  range = Reduce(c,lapply(temprange, format.range))
  pi.est[i] <- mean(eval(parse(text =paste('pi.reg$z[',paste(range, collapse = ','),']', sep = ""))))
}

plot(pi.est~true.density); abline(0,1,col = 'red')



train = sample(1:length(z), length(z)/2)
test = c(1:length(z))[-train]
pi.reg = cde(cbind(x1,x2)[train,], Zst[train], x.margin = cbind(x1,x2)[test,])
