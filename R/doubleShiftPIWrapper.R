#' Double Shift Estimator
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
#' @param nfolds defaults to 2, can only be 1 or 2 at this point
#'
#' @return a list including an estimate of the effect and of its standard deviation

double.shift.pi <- function(y,a,z,delta,x,data = NULL,
                         algo = list(y.est = 'glm',a.est = 'glm',z.est = 'glm'),
                         nfolds = 2,...){
  # want to specify data frame and draw from that w/o attaching
  # would like to be able to do more than 2 folds

  # set up data ----
  n = length(y)
  s1 = sample(1:n, n/2)
  s2 = c(1:n)[-s1]
  s = cbind(s1,s2)

  dat = as.data.frame(cbind(z, x))
  dat.plus = as.data.frame(cbind(z+delta,x))
  dat.min = as.data.frame(cbind(z-delta,x))
  names(dat.plus) <- names(dat.min) <- names(dat)

  psihat <- sd <- rep(NA,nfolds)

  # estimate nuisance parameters and phi ----
  for(i in 1:nfolds){
    train = s[,i]; test = s[,-i]

    ymean = y.mean.est(y[train],dat[train,],algo$y.est)
    amean = a.mean.est(a[train],dat[train,],algo$a.est)

    # get predictions

    # predict y
    if(algo$y.est == 'glm' | algo$y.est == 'random forest'){
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
    psihat[i] = mean(yhat.plus - yhat.min)/mean(ahat.plus - ahat.min)

    # get sd
    n = length(yhat.plus)
    v = mean(  ((yhat.plus - yhat.min)/mean(ahat.plus - ahat.min))^2 )/n
    sd[i] = sqrt(v)

  }

  # average across folds
  psihat = mean(psihat)
  sd = mean(sd)

  return(list(psi = psihat, sd = sd))

}
