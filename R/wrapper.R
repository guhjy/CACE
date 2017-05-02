#' Wrapper
#'
#' @description Wraps up the whole process of estimating phi
#' using sample splitting
#'
#' @return an estimate of the causal effect

wrapper<-function(){

  #load and split data
  loadem()
  mysplit(vis)

  # estimates on first split
  cov = ds1[,c(5:18,20)]
  means1  = mean_est(y=ds1$NCRecid3
                     ,a = ds1$visitslastlocyn1
                     ,z = ds1$total_time
                     ,cov = cov)
  prop1   = propscore_est(y=ds1$NCRecid3,x=cov)
  phi1    = phi_est(ds2,ymean=ymean,amean=amean,p=prop1)

  # estimates on second split
  cov = ds2[,c(5:18,20)]
  means2  = mean_est(y=ds2$NCRecid3
                     ,a = ds2$visitslastlocyn1
                     ,z = ds2$total_time
                     ,cov = cov)
  prop2   = propscore_est(y=ds2$NCRecid3,x=cov)
  phi2    = phi_est(ds1,ymean=ymean,amean=amean,p=prop2)

  # average
  phi     = mean(.5*(phi1 + phi2),na.rm = T)

  return(phi)
}
