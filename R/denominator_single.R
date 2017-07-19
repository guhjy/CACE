#' Compliers
#'
#' @description Take in the estimates done of the nuisance function on
#' the training data. Get predictions using these on the test data.
#' Input those into the formula for the denominator of the parameter.
#' Output the estimate for this train/test combo. Right now this
#' works for the ranger mean estimation. Otherwise you need to change
#' the predict function. I'm going to make this automatic at some point.
#' This function is for the single-shift estimator, meaning the contrast
#' between moving up by delta vs the current value.
#'
#' @param y is the outcome vector
#' @param a is the treatment vector
#' @param z is the instrument vector
#' @param cov is a dataframe of covariates
#' @param ymean is the mean estimates for the outcome
#' variable from the mean_est function
#' @param amean is the mean estimates for the treatment
#' variable from the mean_est function
#' @param p the conditional density estimation object from the
#' propscore_est funciton
#' @param delta the level you want to shift by
#'
#' @return an estimate of the estimated complier proportion

denominator_single <- function(y,a,z,cov,ymean, amean, p, delta = 20){
  print("Estimating Denominator")
  xnew = as.data.frame(cbind(z,cov))
  xnewplus = as.data.frame(cbind(z=z + delta,cov))

  #### A means ####
  # predicted means of a|x,z
  mu_a_xz <- predict(amean,data = xnew)$pred
  mu_a_xzplus <- predict(amean,data = xnewplus)$pred

  #### prop scores ####
  pred <- predict(p, cov)
  pi <- get_probs(z,pred$z,pred$CDE)
  pi_min <- get_probs((z-delta),pred$z,pred$CDE)

  #### estimator ####
  phi_a = (a - mu_a_xz)*pi_min/pi - (a - mu_a_xzplus)

  return(mean(phi_a))
}
