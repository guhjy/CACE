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

phi_est <- function(y,a,z,cov,ymean, amean, p, delta = 2){
  print("Estimating Parameter")
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
  # need to figure out what predict actually does in this case
  pred <- predict(p, cov)
  pi <- get_probs(z,pred$z,pred$CDE)
  pi_plus <- get_probs((z+delta),pred$z,pred$CDE)
  pi_min <- get_probs((z-delta),pred$z,pred$CDE)

  #### estimator ####
  phi_top = (y - mu_y_xz)*(pi_min - pi_plus)/pi + mu_y_xzplus - mu_y_xzmin
  phi_bot = (a - mu_a_xz)*(pi_min - pi_plus)/pi + mu_a_xzplus - mu_a_xzmin
  phi = phi_top/phi_bot

  return(phi)
}
