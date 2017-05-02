#' Estimate Phi on half of data
#'
#' @description Take in the estimates done of the nuisance function on
#' the training data. Get predictions using these on the test data.
#' Input those into the formula for your real parameter. Output the
#' estimate of your parameter for this train/test combo.
#'
#' @param dat is a dataset that is being used as your training
#' dataset overall
#' @param m a list of your mean estimates from the mean_est function
#' @param p the conditional density estimation object from the
#' propscore_est funciton
#' @param delta the level you want to shift by
#'
#' @return an estimate of the causal effect

phi_est <- function(dat,ymean, amean, p, delta = 2){
  attach(dat)

  z <- total_time
  cov = dat[,c(5:18,20)]

  xnew = as.data.frame(cbind(z,cov))
  xnewplus = as.data.frame(cbind(z=z + delta,cov))
  xnewmin = as.data.frame(cbind(z=z - delta,cov))

  #### Y means ####
  # predicted means of y|x,z
  mu_y_xz <- predict(ymean,newdata = xnew)$pred
  mu_y_xzplus <- predict(ymean,newdata = xnewplus)$pred
  mu_y_xzmin <- predict(ymean,newdata = xnewmin)$pred


  #### A means ####
  # predicted means of a|x,z
  mu_a_xz <- predict(amean,newdata = xnew)$pred
  mu_a_xzplus <- predict(amean,newdata = xnewplus)$pred
  mu_a_xzmin <- predict(amean,newdata = xnewmin)$pred

  #### prop scores ####
  pi <- predict(p, xNew = xnew)$CDE[,1] #change this, just for testing
  pi_plus <- predict(p, xNew = xnewplus)$CDE[,1] #change this, just for testing
  pi_min <- predict(p, xNew = xnewmin)$CDE[,1] #change this, just for testing

  #need to then get: pi = prob(Z = z); pi_min = prob(Z = z-delta); pi_plus = prob(Z = z+delta)

  #### estimator ####
  phi_top = (NCRecid3 - mu_y_xz)*(pi_min - pi_plus)/pi + mu_y_xzplus - mu_y_xzmin
  phi_bot = (visitslastlocyn1 - mu_a_xz)*(pi_min - pi_plus)/pi + mu_a_xzplus - mu_a_xzmin
  phi = phi_top/phi_bot

  detach(dat)
  return(phi)
}
