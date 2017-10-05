single.shift.pi <- function(y,a,z,delta,x,data = NULL,
                         algo = list(y.est = 'glm',a.est = 'glm',z.est = 'glm'),
                         nfolds = 2){
  # want to specify data frame and draw from that w/o attaching

  # set up data ----
  n = length(y)
  train = sample(1:n, n/2)
  test = c(1:n)[-train]

  dat = as.data.frame(cbind(z, x))
  dat.plus = as.data.frame(cbind(z+delta,x))
  dat.min = as.data.frame(cbind(z-delta,x))
  names(dat.plus) <- names(dat.min) <- names(dat)

  ymean = y.mean.est(y[train],dat[train,],algo$y.est)
  amean = a.mean.est(a[train],dat[train,],algo$a.est)

  # get predictions ----
  # predict y
  if(algo$y.est == 'glm'){
    yhat = predict(ymean, dat[test,], type = 'response')
    yhat.plus = predict(ymean, dat.plus[test,], type = 'response')
  }
  else{
    yhat = predict(ymean, dat[test,])$pred
    yhat.plus = predict(ymean, dat.plus[test,])$pred
  }

  # predict a
  if(algo$a.est == 'glm'){
    ahat = predict(amean, dat[test,], type = 'response')
    ahat.plus = predict(amean, dat.plus[test,], type = 'response')
  }
  else{
    ahat = predict(amean, dat[test,])$pred
    ahat.plus = predict(amean, dat.plus[test,])$pred
  }

  # get phi ----
  psihat = mean(yhat.plus - yhat)/mean(ahat.plus - ahat)
  n = length(yhat)
  v = mean(  ((yhat.plus - yhat)/(ahat.plus - ahat))^2 )
  sd = sqrt(v)

  return(list(psi = psihat, sd = sd))

}
