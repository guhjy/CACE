#' Wrapper
#'
#' @description Wraps up the whole process of estimating phi
#' using sample splitting
#'
#' @return an estimate of the causal effect

wrapper<-function(){
  loadem()
  mysplit(vis)

  means1  = mean_est(ds1)
  prop1   = propscore_est(ds1)
  phi1    = phi_est(ds2,means1,prop1)

  means2  = mean_est(ds2)
  prop2   = propscore_est(ds2)
  phi2    = phi_est(ds1,means2,prop22)

  phi     = .5*(phi1 + phi2)

  return(phi)
}
