#' Wrapper
#'
#' @description Wraps up the whole process of estimating phi
#' using sample splitting. For now the covariates are hard
#' coded, as are a bunch of things. Will fix later. I'm taking
#' out the loading step so that other people can use this, but
#' it still has to be used on the visits dataset
#'
#' @return an estimate of the causal effect

CACE<-function(y,a,z,cov,delta=2){

  # split data
  data    <- as.data.frame(cbind(y,a,z,cov))
  mysplit(data)

  # estimates on first split
  means1  = mean_est(y=ds1[,1],a = ds1[,2],z = ds1[,3],cov = ds1[,c(4:dim(ds1)[2])])
  prop1   = propscore_est(y=ds1[,1],x=ds1[,c(4:dim(cov)[2])])
  phi1    = phi_est(y=ds2[,1],a = ds2[,2],z = ds2[,3],cov = ds2[,c(4:dim(ds2)[2])]
                    ,ymean=ymean,amean=amean,p=prop1,delta=delta)

  # estimates on second split
  means2  = mean_est(y=ds2[,1],a = ds2[,2],z = ds2[,3],cov = ds2[,c(4:dim(ds2)[2])])
  prop2   = propscore_est(y=ds2[,1],x=ds2[,c(4:dim(cov)[2])])
  phi2    = phi_est(y=ds1[,1],a = ds1[,2],z = ds1[,3],cov = ds1[,c(4:dim(ds1)[2])]
                    ,ymean=ymean,amean=amean,p=prop2,delta=delta)

  # average
  phi     = mean(.5*(phi1 + phi2),na.rm = T)

  return(phi)
}
