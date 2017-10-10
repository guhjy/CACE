single.shift <- function(y,a,z,delta,x,data = NULL,
                         algo = list(y.est = 'glm',a.est = 'glm',z.est = 'glm'),
                         nfolds = 2,...){
  # want to specify data frame and draw from that w/o attaching
  # need to make this do repeat for each half and then average

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
    zmean = z.condldens.est(z[train],x[train,],algo$z.est)

    # get predictions ----
    # predict y
    if(algo$y.est == 'glm' | algo$y.est == 'random forest'){
      yhat = predict(ymean, newdata = dat[test,], type = 'response')
      yhat.plus = predict(ymean, newdata = dat.plus[test,], type = 'response')
    }
    else{
      yhat = predict(ymean, dat[test,])$pred
      yhat.plus = predict(ymean, dat.plus[test,])$pred
    }

    # predict a
    if(algo$a.est == 'glm' | algo$a.est == 'random forest'){
      ahat = predict(amean, newdata = dat[test,], type = 'response')
      ahat.plus = predict(amean, newdata = dat.plus[test,], type = 'response')
    }
    else{
      ahat = predict(amean, dat[test,])$pred
      ahat.plus = predict(amean, dat.plus[test,])$pred
    }

    # predict z
    if(algo$z.est == 'glm'){
      zhat <- predict(zmean, dat[test,], type = 'response')
      z.var <- mean( (z - zhat)^2  )

      gK <- function(x){(1/sqrt(2*pi))*exp(-(x^2)/2)}
      pihat <- sapply(z, function(y) (1/N)*sum(gK(sqrt( ((y - zhat))^2/z.var ) )))
      pihat.min <- sapply((z-delta), function(y) (1/N)*sum(gK(sqrt( ((y - zhat))^2/z.var ) )))
    }
    else{
      pred = predict(zmean, dat[test,])
      pihat = get_probs(z[test], pred$z, pred$CDE)
      pihat.min = get_probs((z-delta)[test], pred$z, pred$CDE)
    }

    # get phi ----
    phi_y = (y[test] - yhat)*(pihat.min/pihat) - (y[test] - yhat.plus)
    phi_a = (a[test] - ahat)*(pihat.min/pihat) - (a[test] - ahat.plus)

    if(length(which(pi==0))>0){warning(paste("Number of zero probability values (positivity violation):",length(which(pi==0))))}
    keep = which(pihat!=0)
    psihat[i] = mean(phi_y[keep])/mean(phi_a[keep])

    n = length(phi_y[keep])
    top = phi_y[keep] - psihat[i]*phi_a[keep]; bottom = mean(phi_a[keep])
    v = mean( ( top/bottom )^2  )/ n
    sd[i] = sqrt(v)

  }

  psihat = mean(psihat)
  sd = mean(sd)

  return(list(psi = psihat, sd = sd))

}
