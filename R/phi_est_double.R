#' Estimate Phi on half of data
#'
#' @description Take in the estimates done of the nuisance function on
#' the training data. Get predictions using these on the test data.
#' Input those into the formula for your real parameter. Output the
#' estimate of your parameter for this train/test combo. Right now this
#' works for the ranger mean estimation. Otherwise you need to change
#' the predict function. I'm going to make this automatic at some point.
#' This function is for the double-shift estimator, meaning the contrast
#' between moving up by delta vs moving down by delta.
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
#' @return an estimate of the causal effect

phi_est_double_ranger <- function(y,a,z,cov,ymean, amean, p, delta = 20){
  print("Estimating Parameter")
  xnew = as.data.frame(cbind(z,cov))
  xnewplus = as.data.frame(cbind(z=z + delta,cov))
  xnewmin = as.data.frame(cbind(z=z - delta,cov))

  #### Y means ####
  # predicted means of y|x,z
  mu_y_xz <- predict(ymean,data = xnew)$pred
  mu_y_xzplus <- predict(ymean,data = xnewplus)$pred
  mu_y_xzmin <- predict(ymean,data = xnewmin)$pred

  #### A means ####
  # predicted means of a|x,z
  mu_a_xz <- predict(amean,data = xnew)$pred
  mu_a_xzplus <- predict(amean,data = xnewplus)$pred
  mu_a_xzmin <- predict(amean,data = xnewmin)$pred

  #### prop scores ####
  pred <- predict(p, cov)
  pi <- get_probs(z,pred$z,pred$CDE)
  pi_plus <- get_probs((z+delta),pred$z,pred$CDE)
  pi_min <- get_probs((z-delta),pred$z,pred$CDE)

  #### estimator ####
  phi_y = ( (y - mu_y_xz)*(pi_min - pi_plus)/pi ) + mu_y_xzplus - mu_y_xzmin
  phi_a = ( (a - mu_a_xz)*(pi_min - pi_plus)/pi ) + mu_a_xzplus - mu_a_xzmin
  print(paste(length(which(pi==0)),"zero probability values"));keep = which(pi!=0)
  psihat = mean(phi_y[keep])/mean(phi_a[keep])
  n = length(phi_y[keep])

  v = mean( ((phi_y[keep] - psihat*phi_a[keep])/mean(phi_a[keep]))^2  )/n
  sd = sqrt(v)

  return(list(phi = psihat, sd = sd))
}
