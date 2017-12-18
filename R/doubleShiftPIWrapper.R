#' Double Shift Plug in Estimator
#'
#' @description Wraps up the whole process of estimating phi
#' using sample splitting.
#'
#' @param y your outcome variable (a vector)
#' @param a your treatment variable (a vector)
#' @param z your instrument (a vector)
#' @param x a dataframe of covariates
#' @param delta is the amount you want to shift by
#' @param algo a list of three algorithms you want to use: y.est, a.est and z.est
#' @param nfolds defaults to 2
#' @param zmax the upper bound on Z, default is Inf
#' @param zmin the lower bound on Z, default is -Inf
#'
#' @return a list including an estimate of the effect and of its standard deviation.
#' Standard deviation is not valid in general.

double.shift.pi <- function(y,a,z,delta,x,data = NULL,
                         algo = list(y.est = 'glm',a.est = 'glm',z.est = 'glm'),
                         nfolds = 2,zmax = Inf, zmin = -Inf,...){
  # want to specify data frame and draw from that w/o attaching
  # should be able to do more than 2 folds

  full.dat <- cbind(y,a,z,x)
  keep <- complete.cases(full.dat)
  y <- y[keep]; a <- a[keep]; z <- z[keep]; x <- x[keep,]

  # set up data ----
  n = length(y)
  s = sample(rep(1:nfolds,ceiling(n/nfolds))[1:n])

  dat = as.data.frame(cbind(z, x))
  dat.plus = as.data.frame(cbind(z+delta,x))
  dat.min = as.data.frame(cbind(z-delta,x))
  names(dat.plus) <- names(dat.min) <- names(dat)

  psihat <- sd <- rep(NA,nfolds)

  # estimate nuisance parameters and phi ----
  for(i in 1:nfolds){
    train = (s!=i); test = (s==i)

    ymean = y.mean.est(y[train],dat[train,],algo$y.est)
    amean = a.mean.est(a[train],dat[train,],algo$a.est)

    # get predictions

    # predict y
    if(algo$y.est == 'glm' | algo$y.est == 'random forest'){
      yhat = predict(ymean, newdata = dat[test,], type = 'response')
      yhat.plus = predict(ymean, newdata = dat.plus[test,], type = 'response')
      yhat.min = predict(ymean, newdata = dat.min[test,], type = 'response')
    }
    else{
      yhat.plus = predict(ymean, dat.plus[test,])$pred
      yhat.min = predict(ymean, dat.min[test,])$pred
    }

    # predict a
    if(algo$a.est == 'glm' | algo$a.est == 'random forest'){
      ahat.plus = predict(amean, newdata = dat.plus[test,], type = 'response')
      ahat.min = predict(amean, newdata = dat.min[test,], type = 'response')
    }
    else{
      ahat.plus = predict(amean, dat.plus[test,])$pred
      ahat.min = predict(amean, dat.min[test,])$pred
    }

    # get phi
    psihat[i] = mean(yhat.plus*(z[test]+delta<zmax) - yhat.min*(z[test]-delta > zmin))/mean(ahat.plus*(z[test]+delta<zmax) - ahat.min*(z[test]-delta > zmin))
    n = length(yhat.plus)
    v = var(  ((yhat.plus*(z[test]+delta<zmax) - yhat.min*(z[test]-delta > zmin))/mean(ahat.plus*(z[test]+delta<zmax) - ahat.min*(z[test]-delta > zmin))) )
    sd[i] = sqrt(v)
  }

  # average across folds
  psihat = mean(psihat)
  sd = mean(sd)

  return(list(psi = psihat, sd = sd))

}
