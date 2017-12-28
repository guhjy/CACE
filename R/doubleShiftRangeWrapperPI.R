#' Double Shift Plug in Estimator over range of deltas
#'
#' @description Wraps up the whole process of estimating phi
#' using sample splitting.
#'
#' @param y your outcome variable (a vector)
#' @param a your treatment variable (a vector)
#' @param z your instrument (a vector)
#' @param x a dataframe of covariates
#' @param delta is the amount you want to shift by (a vector)
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
#' @param nfolds defaults to 2
#' @param zmax the upper bound on Z, default is Inf
#' @param zmin the lower bound on Z, default is -Inf
#'
#' @return a list including an estimate of the effect and of its standard deviation.
#' Standard deviation is not valid in general.

double.shift.rangePI <- function(y,a,z,x,delta,Y.est,A.est,Z.est,
                               nfolds = 2, zmax = Inf, zmin = -Inf, ...){

  #move these to their own function
  N = length(y); j = length(delta)
  if(any(length(a)!=N,length(z)!=N,dim(x)[1]!=N)){stop('y,a,z and x must be same length')}
  if(zmax <= zmin){stop('zmax must be bigger than zmin')}
  keep <- complete.cases(cbind(y,a,z,x))
  y <- y[keep]; a <- a[keep]; z <- z[keep]; x <- x[keep,]

  n = length(y)
  s = sample(rep(1:nfolds,ceiling(n/nfolds))[1:n])

  dat = as.data.frame(cbind(z, x))
  psihat <- sd <- matrix(rep(NA,nfolds*j), ncol = j, nrow = nfolds)
  ifvals <- matrix(rep(NA,n*j), ncol = j, nrow = n)

  for(i in 1:nfolds){
    train = (s!=i); test = (s==i)
    if(nfolds == 1){train <- test}

    ymean = y.mean.est(y[train],dat[train,],Y.est)
    amean = a.mean.est(a[train],dat[train,],A.est)
    zmean = z.condldens.est(z[train],x[train,],Z.est)

    mu.hat = my.pred(model = ymean, newdata = dat[test,], algo = Y.est)
    lambda.hat = my.pred(model = amean, newdata = dat[test,], algo = A.est)

    dat.plus = lapply(delta, function(d) as.data.frame(cbind(z+d,x)))
    dat.min = lapply(delta, function(d) as.data.frame(cbind(z-d,x)))

    for(ii in 1:length(dat.plus)){names(dat.plus[[ii]]) <- names(dat.min[[ii]]) <- names(dat)}

    mu.hat.plus <- lapply(dat.plus, function(df) my.pred(model = ymean, newdata = df[test,], algo = Y.est))
    lambda.hat.plus <- lapply(dat.plus, function(df) my.pred(model = amean, newdata = df[test,], algo = A.est))

    mu.hat.min <- lapply(dat.min, function(df) my.pred(model = ymean, newdata = df[test,], algo = Y.est))
    lambda.hat.min <- lapply(dat.min, function(df) my.pred(model = amean, newdata = df[test,], algo = A.est))

    for(jj in 1:j){
      d = delta[jj]
      mp = mu.hat.plus[[jj]];  mm = mu.hat.min[[jj]]
      lp = lambda.hat.plus[[jj]];  lm = lambda.hat.min[[jj]]

      psihat[i,jj] = mean(mp*(z[test]+d<zmax) - mm*((z[test]-d) > zmin), na.rm = T)/
        mean(lp*(z[test]+d<zmax) - lm*((z[test]-d) > zmin), na.rm = T)
    }
  }
  # average within delta values
  est.eff = apply(psihat, 2, mean, na.rm = T)

  out = list(Estimate = est.eff, delta = delta)

  return(out)
}
