#' Local Average Treatment Effect Estimator
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

CACE<-function(y,a,z,cov,delta=2,ranger = F,type = 'double',quiet = T){
  ptm <- proc.time()

  # what kind of estimator do you want?
  if(type=='double'){phi_est = phi_est_double}
  if(type=='single'){phi_est = phi_est_single}
  if(type=='min'){phi_est = phi_est_min}
  if(type == 'plugin'){phi_est = plugin}

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

  # first split
  cat('round 1')
  prop1   = propscore_est(y=ds1[,3],x=ds1[,c(4:dim(ds1)[2])])
  out1    = phi_est(y=ds2[,1],a = ds2[,2],z = ds2[,3],cov = ds2[,c(4:dim(ds2)[2])],ymean=ymean,amean=amean,p=prop1,delta=delta)
  phi1    = out1$phi
  num1    = out1$numerator
  den1    = out1$denominator

  # second split
  cat('round 2')
  prop2   = propscore_est(y=ds2[,3],x=ds2[,c(4:dim(ds2)[2])])
  out2    = phi_est(y=ds1[,1],a = ds1[,2],z = ds1[,3],cov = ds1[,c(4:dim(ds1)[2])],ymean=ymean,amean=amean,p=prop2,delta=delta)
  phi2    = out2$phi
  num2    = out2$numerator
  den2    = out2$denominator

  # average
  phi     = .5*(phi1 + phi2)
  sd      = .5*(out1$sd + out2$sd)
  num     = .5*(num1 + num2)
  den     = .5*(den1 + den2)

  print(paste("Total estimation runtime:",(proc.time()-ptm)[1]))
  return(list(phi=phi,sd = sd,numerator = num, denominator = den))
}
