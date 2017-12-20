#' Double Shift Estimator Across Delta Range
#'
#' @description Estimates the effects of shifting up and down
#' across a range of values. Includes the multiplier bootstrap
#' for calculations of uniform confidence bands
#'
#' @param y your outcome variable (a vector)
#' @param a your treatment variable (a vector)
#' @param z your instrument (a vector)
#' @param x a dataframe of covariates
#' @param delta a vector of shift levels
#' @param Y.est an algorithm for estimating the means of your outcome.
#' Should be 'glm', 'superlearner' or 'ranger'. If you choose superlearner,
#' you can also specify the libraries you would like to use. Default libraries
#' are c("SL.glm","SL.randomForest","SL.polymars","SL.mean").
#' @param A.est an algorithm for estimating the means of your treatment.
#' Should be 'glm', 'superlearner' or 'ranger'. If you choose superlearner,
#' you can also specify the libraries you would like to use. Default libraries
#' are c("SL.glm","SL.randomForest","SL.polymars","SL.mean").
#' @param Z.est an algorithm for estimating the means of your treatment.
#' Should be 'glm' in which case a glm is used to estimte the mean and variance
#' and a kernel is used to estimate the density, or 'flexcode'. If you choose
#' 'flexcode', you can specify the regression function. The default is regressionFunction.NW
#' @param nfolds number of folds. Defaults to 2
#' @param zmax the upper bound on Z, default is Inf
#' @param zmin the lower bound on Z, default is -Inf
#' @param alpha the alpha level for your multiplier bootstrap confidence bands (default 0.05)
#' @param nbs the number of boostrap samples to take. Default is 10,000.
#'
#' @return a list including an estimate of the effect for each delta value and its standard deviation,
#' as well as upper and lower confidence intervals for the pointwise and uniform case.


double.shift.range <- function(y,a,z,x,delta,Y.est,A.est,Z.est,
                               nfolds = 2, zmax = Inf, zmin = -Inf, alpha = 0.05, nbs = 10000, ...){

  #move these to their own function
  N = length(y); j = length(delta)
  if(any(length(a),length(z),dim(x)[1])!=N){stop('y,a,z and x must be same length')}
  if(zmax <= zmin){stop('zmax must be bigger than zmin')}
  keep <- complete.cases(cbind(y,a,z,x))
  y <- y[keep]; a <- a[keep]; z <- z[keep]; x <- x[keep,]

  n = length(y)
  s = sample(rep(1:nfolds,ceiling(n/nfolds))[1:n])

  dat = as.data.frame(cbind(z, x))
  psihat <- sd <- matrix(rep(NA,nfolds*j), ncol = j, nrow = nfolds)

  for(i in 1:nfolds){
    train = (s!=i); test = (s==i)
    if(nfolds == 1){train <- test}

    ymean = y.mean.est(y[train],dat[train,],Y.est)
    amean = a.mean.est(a[train],dat[train,],A.est)
    zmean = z.condldens.est(z[train],x[train,],Z.est)

    mu.hat = my.pred(model = ymean, newdata = dat[test,], algo = Y.est)
    lambda.hat = my.pred(model = amean, newdata = dat[test,], algo = A.est)
    pi.hat = my.zpred(model = zmean, z = z[test], newdata = dat[test,], algo = Z.est)

    psihat[i,] = sapply(delta, function(x) psi.est(x))

  }
  # average within delta values
  est.eff = apply(psihat, 2, mean)
  sigma = sqrt(apply(psihat, 2, var))

  # pointwise confidence interval
  ll1 <- est.eff-qnorm(1-alpha/2)*sigma/sqrt(n); ul1 <- est.eff+qnorm(1-alpha/2)*sigma/sqrt(n)

  # multiplier bootstrap
  eff.mat <- matrix(rep(est.eff,n), nrow = n, byrow = T)
  sig.mat <- matrix(rep(sd,n), nrow = n, byrow = T)
  ifvals.std <- (psihat - eff.mat)/sig.mat
  mult <- matrix(2*rbinom(n*nbs, 1, .5)-1, nrow = n, ncol = nbs)
  maxvals <- sapply(1:nbs, function(col){
    max(abs(apply(mult[,col]*ifvals.std,2,sum)/sqrt(n)))
  })
  calpha <- quantile(maxvals, 1-alpha)
  eff.ll <- est.eff - calpha*sigma/sqrt(n)
  eff.ul <- est.eff + calpha*sigma/sqrt(n)

  out = list(Estimate = est.eff, SD = sigma, CI.lower.point = ll1, CI.upper.point = ul1,
             MBquantile = calpha, CI.lower.unif = ll2, CI.upper.unif = ul2)
}

# these need to be tested, esp for superlearner
my.pred <- function(model, newdata, algo){
  if(algo == 'ranger'){predict(model, newdata)$pred}
  else if(algo == 'superlearner'){predict(model, newdata, onlySL = T)$pred}
  else{predict(model, newdata, type = 'response')}
}

my.zpred <- function(model, z, newdata, algo){
  if(algo == 'glm'){
    zhat = predict(model, newdata, type = 'response')
    z.var = mean( (z - zhat)^2  )
    N = length(zhat)

    gK <- function(x){(1/sqrt(2*pi))*exp(-(x^2)/2)}
    out = sapply(z, function(y) (1/N)*sum(gK(sqrt( ((y - zhat))^2/z.var ) )))
  }
  else if(algo == 'flexcode'){
    pred = predict(model, newdata)
    out = get_probs(z, pred$z, pred$CDE)
  }
  if(length(out==0)>0){warning('0 values found for pi, setting to NA')}
  out[out==0] <- NA
  return(out)
}


psi.est <- function(delta){
  dat.plus = as.data.frame(cbind(z+delta, x))
  dat.min = as.data.frame(cbind(z-delta, x))

  mu.hat.plus = my.pred(model = ymean, newdata = dat.plus[test,], algo = Y.est)
  lambda.hat.plus = my.pred(model = amean, newdata = dat.plus[test,], algo = A.est)
  pi.hat.plus = my.zpred(model = zmean, z = z[test], newdata = dat.plus[test,], algo = Z.est)

  mu.hat.min = my.pred(model = ymean, newdata = dat.min[test,], algo = Y.est)
  lambda.hat.min = my.pred(model = amean, newdata = dat.min[test,], algo = A.est)
  pi.hat.min = my.zpred(model = zmean, z = z[test], newdata = dat.min[test,], algo = Z.est)

  xi.y.delta = ((y[test] - mu.hat)*(pi.hat.min/pi.hat) - (y[test] - mu.hat.plus))*((z[test]+delta) < zmax)
  xi.y.Mdelta = ((y[test] - mu.hat)*(pi.hat.plus/pi.hat) - (y[test] - mu.hat.min))*((z[test]-delta) > zmin)

  xi.a.delta = ((a[test] - lambda.hat)*(pi.hat.min/pi.hat) - (a[test] - lambda.hat.plus))*((z[test]+delta) < zmax)
  xi.a.Mdelta = ((a[test] - lambda.hat)*(pi.hat.plus/pi.hat) - (a[test] - lambda.hat.min))*((z[test]-delta) > zmin)

  mean(xi.y.delta - xi.y.Mdelta)/mean(xi.a.delta - xi.a.Mdelta)
}
