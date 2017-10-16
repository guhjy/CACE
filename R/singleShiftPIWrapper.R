#' Single Shift Estimator
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

single.shift.pi <- function(y,a,z,delta,x,data = NULL,
                         algo = list(y.est = 'glm',a.est = 'glm',z.est = 'glm'),
                         nfolds = 2){
  # want to specify data frame and draw from that w/o attaching

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

  for(i in 1:nfolds){
    train = s[,i]; test = s[,-i]

    ymean = y.mean.est(y[train],dat[train,],algo$y.est)
    amean = a.mean.est(a[train],dat[train,],algo$a.est)

    # get predictions ----
    # predict y
    if(algo$y.est == 'glm'| algo$y.est == 'random forest'){
      yhat = predict(ymean, dat[test,], type = 'response')
      yhat.plus = predict(ymean, dat.plus[test,], type = 'response')
    }
    else{
      yhat = predict(ymean, dat[test,])$pred
      yhat.plus = predict(ymean, dat.plus[test,])$pred
    }

    # predict a
    if(algo$a.est == 'glm'| algo$a.est == 'random forest'){
      ahat = predict(amean, dat[test,], type = 'response')
      ahat.plus = predict(amean, dat.plus[test,], type = 'response')
    }
    else{
      ahat = predict(amean, dat[test,])$pred
      ahat.plus = predict(amean, dat.plus[test,])$pred
    }

    # get phi ----
    psihat[i] = mean(yhat.plus - yhat)/mean(ahat.plus - ahat)
    n = length(yhat)
    v = mean(  ((yhat.plus - yhat)/mean(ahat.plus - ahat))^2 )/n
    sd[i] = sqrt(v)
  }

  psihat = mean(psihat)
  sd = mean(sd)

  return(list(psi = psihat, sd = sd))

}
