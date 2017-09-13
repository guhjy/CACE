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

CACE<-function(y,a,z,cov,delta=2,ranger = F,type = 'double',quiet = T,split = T){
  ptm <- proc.time()

  if(split == T){
    # split data
    data    <- as.data.frame(cbind(y,a,z,cov))
    mysplit(data)
    res1 = bundle(ds1,ranger = ranger,type = type, delta = delta)
    res2 = bundle(ds2,ranger = ranger,type = type, delta = delta)

    # average
    phi     = .5*(res1$phi + res2$phi)
    sd      = .5*(res1$sd + res2$sd)
    num     = .5*(res1$numerator + res2$numerator)
    den     = .5*(res1$denominator + res2$denominator)

    print(paste("Total estimation runtime:",(proc.time()-ptm)[1]))
    return(list(phi=phi,sd = sd,numerator = num, denominator = den))
  }
  else{
    data <- as.data.frame(cbind(y,a,z,cov))
    data <- data[complete.cases(data),]
    res = bundle(data,ranger = ranger,type = type, delta = delta)
    print(paste("Total estimation runtime:",(proc.time()-ptm)[1]))
    return(list(phi=res$phi,sd = res$sd,numerator = res$numerator, denominator = res$denominator))
  }
}
