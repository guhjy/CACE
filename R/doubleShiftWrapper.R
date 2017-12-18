#' Double Shift IF Estimator
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



double.shift <- function(y,a,z,delta,x,data = NULL,
                         algo = list(y.est = 'glm',a.est = 'glm',z.est = 'glm'),
                         nfolds = 2,
                         zmax = Inf, zmin = -Inf, ...){
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
    if(nfolds == 1){train <- test}

    ymean = y.mean.est(y[train],dat[train,],algo$y.est)
    amean = a.mean.est(a[train],dat[train,],algo$a.est)
    zmean = z.condldens.est(z[train],x[train,],algo$z.est)

    # get predictions

    # predict y
    if(algo$y.est == 'glm' | algo$y.est == 'random forest'){
      yhat = predict(ymean, newdata = dat[test,], type = 'response')
      yhat.plus = predict(ymean, newdata = dat.plus[test,], type = 'response')
      yhat.min = predict(ymean, newdata = dat.min[test,], type = 'response')
    }
    else{
      yhat = predict(ymean, dat[test,])$pred
      yhat.plus = predict(ymean, dat.plus[test,])$pred
      yhat.min = predict(ymean, dat.min[test,])$pred
    }

    # predict a
    if(algo$a.est == 'glm' | algo$a.est == 'random forest'){
      ahat = predict(amean, newdata = dat[test,], type = 'response')
      ahat.plus = predict(amean, newdata = dat.plus[test,], type = 'response')
      ahat.min = predict(amean, newdata = dat.min[test,], type = 'response')
    }
    else{
      ahat = predict(amean, dat[test,])$pred
      ahat.plus = predict(amean, dat.plus[test,])$pred
      ahat.min = predict(amean, dat.min[test,])$pred
    }

    # predict z
    if(algo$z.est == 'glm'){
      zhat <- predict(zmean, dat[test,], type = 'response')
      z.var <- mean( (z - zhat)^2  )
      N = length(zhat)

      gK <- function(x){(1/sqrt(2*pi))*exp(-(x^2)/2)}
      pihat <- sapply(z, function(y) (1/N)*sum(gK(sqrt( ((y - zhat))^2/z.var ) )))
      pihat.min <- sapply((z-delta), function(y) (1/N)*sum(gK(sqrt( ((y - zhat))^2/z.var ) )/sqrt(z.var)))
      pihat.plus <- sapply((z+delta), function(y) (1/N)*sum(gK(sqrt( ((y - zhat))^2/z.var ) )/sqrt(z.var)))
    }
    else{
      pred = predict(zmean, dat[test,])
      pihat = get_probs(z[test], pred$z, pred$CDE)
      pihat.min = get_probs((z-delta)[test], pred$z, pred$CDE)
      pihat.plus = get_probs((z+delta)[test], pred$z, pred$CDE)
    }

    # get phi
    phi_y1 = ((y[test] - yhat)*(pihat.min/pihat) - (y[test] - yhat.plus))*((z[test]+delta) < zmax)
    phi_y2 = ((y[test] - yhat)*(pihat.plus/pihat) - (y[test] - yhat.min))*((z[test]-delta) > zmin)
    phi_a1 = ((a[test] - ahat)*(pihat.min/pihat) - (a[test] - ahat.plus))*((z[test]+delta) < zmax)
    phi_a2 = ((a[test] - ahat)*(pihat.plus/pihat) - (a[test] - ahat.min))*((z[test]-delta) > zmin)

    if(length(which(pi==0))>0){warning(paste("Number of zero probability values (positivity violation):",length(which(pi==0))))}
    pos = which(pihat!=0)
    psihat[i] = mean((phi_y1-phi_y2)[pos])/mean((phi_a1 - phi_a2)[pos])

    # get sd
    n = length(phi_y1[pos])
    top = (phi_y1-phi_y2)[pos] - psihat[i]*(phi_a1 - phi_a2)[pos]
    bottom = mean((phi_a1 - phi_a2)[pos])
    v = mean( ( top/bottom )^2  )/ n
    sd[i] = sqrt(v)

  }

  # average across folds
  psihat = mean(psihat)
  sd = mean(sd)

  return(list(psi = psihat, sd = sd))

}
