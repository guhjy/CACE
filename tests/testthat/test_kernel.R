### test the kernel ###
N = 1000
x1 <- rnorm(N)
x2 <- rnorm(N)
y <- 1 + 2*x1 - 3*x2 + rnorm(N)
true.dens <- dnorm(y, mean = 1 + 2*x1 - 3*x2, sd = sqrt(14))

s <- sample(c(TRUE,FALSE), length(x1), replace = TRUE)
rDens <- density(y[s])

for(i in 1:sum(!s)){
  dens.est[i] = (rDens$y[rank(append(y[!s][i],rDens$x))[1]]+rDens$y[rank(append(y[!s][i],rDens$x))[1]+1])/2
}
plot(dens.est~true.dens[!s])

# slightly better
rCondlDens <- cde(cbind(x1,x2)[s,], y[s],nxmargin = 100)
X = cbind(x1,x2)
dens.est <- rep(NA,sum(!s))
for(i in 1:sum(!s)){
  x.rank = rep(NA,length(rCondlDens$x))
  for(j in 1:length(rCondlDens$x)){
    x.rank[j] = min(rank(append(X[!s,j][i], rCondlDens$x[[j]]))[1], length(rCondlDens$x[[j]])-1)
  }
  y.rank = min(rank(append(y[!s][i], rCondlDens$y))[1], length(rCondlDens$y)-1)
  ranks = paste(c(x.rank,y.rank),collapse = ",")
  ranks1 = paste(c(x.rank,y.rank)+1,collapse = ",")
  low.dens = eval(parse(text = paste("rCondlDens$z[",ranks,"]", sep = "")))
  high.dens = eval(parse(text = paste("rCondlDens$z[",ranks1,"]", sep = "")))

  dens.est[i] = (low.dens+high.dens)/2
}
plot(dens.est~true.dens[!s])




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

# this is all half-baked: not really set up to look at 2 x's
# try out using R's density function
x = rnorm(1000)
z = rnorm(length(x), 2*x, 1)
true.density = dnorm(z, 2*x, sqrt(2))
mean.est = predict(glm(z~x, family = 'gaussian'))
sd.est = sqrt( mean(mean.est - z)^2 / length(z) )
Zst = as.numeric( (z-mean.est)/sd.est )

x1 = rnorm(1000); x2 = rnorm(1000)
z = x1 + x2 + rnorm(1000)
true.density = dnorm(z,mean = x1+x2, sd = sqrt(3))
mean.est = predict(glm(z~x1+x2, family = 'gaussian'))
sd.est = sqrt( mean( (mean.est - z)^2 ) )
Zst = as.numeric( (z - mean.est) /sd.est)

pi.est = sapply(z, function(y) (1/length(x))*sum(gK((1/sd.est)*sqrt((y - x)^2))/sd.est) )
pi.est = sapply(z, function(y) (1/length(x))*sum( gK( y - mean.est )/sd.est ))
pi.est = sapply(z, function(y) (1/length(x))*sum( gK( (y - mean.est)/sd.est )/sd.est ))

piR = density(Zst)
for(i in 1:length(Zst)){
  pi.est[i] = mean(piR$y[rank(append(Zst[i],piR$x))[1]],piR$y[rank(append(Zst[i],piR$x))[1]-1])/sd(z)
}

lims = c(min(pi.est,true.density), max(pi.est,true.density))
plot(pi.est ~ true.density, xlim = lims, ylim = lims)
abline(0,1)


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
