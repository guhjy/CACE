#' Wrapper
#'
#' @description Wraps up the whole process of estimating phi
#' using sample splitting. Note that a lot of problems you
#' might face at this point will be due to the data you are
#' inputting. It has trouble with factors and with intercept
#' terms.
#'
#' @param y your outcome variable (a vector)
#' @param a your treatment variable (a vector)
#' @param z your instrument (a vector)
#' @param cov a dataframe of covariates
#' @param delta is the amount you want to shift by
#' @param type is either 'double' meaning you shift both up and
#' down by delta, or 'single' meaning you shift up.
#'
#' @return an estimate of the causal effect

CACE<-function(y,a,z,cov,delta=2,ranger = F,type = 'double'){
  ptm <- proc.time()
  # split data
  data    <- as.data.frame(cbind(y,a,z,cov))
  mysplit(data)

  if(ranger == F){
    # estimate means using what you've specified
    means1  = mean_est(y=ds1[,1],a = ds1[,2],z = ds1[,3],cov = ds1[,c(4:dim(ds1)[2])])
    means2  = mean_est(y=ds2[,1],a = ds2[,2],z = ds2[,3],cov = ds2[,c(4:dim(ds2)[2])])
  }

  else{
    # estimates neans using ranger
    means1  = mean_est_ranger(y=ds1[,1],a = ds1[,2],z = ds1[,3],cov = ds1[,c(4:dim(ds1)[2])])
    means2  = mean_est_ranger(y=ds2[,1],a = ds2[,2],z = ds2[,3],cov = ds2[,c(4:dim(ds2)[2])])

  }

  # what kind of estimator do you want?
  if(type=='double'){phi_est = phi_est_double}
  if(type=='single'){phi_est = phi_est_single}
  if(type=='min'){phi_est = phi_est_min}

  # first split
  prop1   = propscore_est(y=ds1[,3],x=ds1[,c(4:dim(ds2)[2])])
  out1    = phi_est(y=ds2[,1],a = ds2[,2],z = ds2[,3],cov = ds2[,c(4:dim(ds2)[2])],ymean=ymean,amean=amean,p=prop1,delta=delta)
  phi1    = out1$phi

  # second split
  prop2   = propscore_est(y=ds2[,3],x=ds2[,c(4:dim(ds2)[2])])
  out2    = phi_est(y=ds1[,1],a = ds1[,2],z = ds1[,3],cov = ds1[,c(4:dim(ds1)[2])],ymean=ymean,amean=amean,p=prop2,delta=delta)
  phi2    = out2$phi

  # average
  phi     = .5*(phi1 + phi2)
  vr      = .5*(out1$var + out2$var)

  print(paste("Total estimation runtime:",(proc.time()-ptm)[1]))
  return(list(phi=phi,var = vr))
}
