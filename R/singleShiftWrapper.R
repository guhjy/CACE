single.shift <- function(y,a,z,delta,x,data = NULL,
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
  zmean = z.condldens.est(z[train],x[train,],algo$z.est)

  # get predictions ----
  # predict y
  if(algo$y.est == 'glm'){
    yhat = predict(ymean, newdata = dat[test,], type = 'response')
    yhat.plus = predict(ymean, newdata = dat.plus[test,], type = 'response')
    }
  else{
    yhat = predict(ymean, dat[test,])$pred
    yhat.plus = predict(ymean, dat.plus[test,])$pred
  }

  # predict a
  if(algo$a.est == 'glm'){
    ahat = predict(amean, newdata = dat[test,], type = 'response')
    ahat.plus = predict(amean, newdata = dat.plus[test,], type = 'response')
  }
  else{
    ahat = predict(amean, dat[test,])$pred
    ahat.plus = predict(amean, dat.plus[test,])$pred
  }

  # predict z
  if(algo$z.est == 'glm'){
    pihat = predict(zmean, dat[test,])
    pihat.min = predict(zmean, dat.min[test,])
  }
  else{
    pred = predict(zmean, dat[test,])
    pihat = get_probs(z[test], pred$z, pred$CDE)
    pihat.min = get_probs((z-delta)[test], pred$z, pred$CDE)
  }

  # get phi ----
  phi_y = (y - yhat)*(pihat.min/pihat) - (y - yhat.plus)
  phi_a = (a - ahat)*(pihat.min/pihat) - (a - ahat.plus)

  if(length(which(pi==0))>0){warning(paste("Number of zero probability values (positivity violation):",length(which(pi==0))))}
  keep = which(pihat!=0)
  psihat = mean(phi_y[keep])/mean(phi_a[keep])

  n = length(phi_y[keep])
  top = phi_y[keep] - psihat*phi_a[keep]; bottom = mean(phi_a[keep])
  v = mean( ( top/bottom )^2  )/ n
  sd = sqrt(v)

  return(list(psi = psihat, sd = sd))

}
