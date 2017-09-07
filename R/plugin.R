#' Estimate plug-in on half of data
#'
#' @description Take in the estimates done of the nuisance function on
#' the training data. Get predictions using these on the test data.
#' Input those into the plug in formula for your real parameter. Output the
#' estimate of your parameter for this train/test combo.
#' The first function is for the single-shift estimator, meaning the contrast
#' between moving up by delta vs the current value. The second is for the
#' double shift.
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

plugin_ranger <- function(y,a,z,cov,ymean, amean, p, delta = 20){
  print("Estimating Parameter")
  xnew = as.data.frame(cbind(z,cov))
  xnewplus = as.data.frame(cbind(z=z + delta,cov))

  #### Y means ####
  # predicted means of y|x,z
  mu_y_xz <- predict(ymean,data = xnew)$pred
  mu_y_xzplus <- predict(ymean,data = xnewplus)$pred

  #### A means ####
  # predicted means of a|x,z
  mu_a_xz <- predict(amean,data = xnew)$pred
  mu_a_xzplus <- predict(amean,data = xnewplus)$pred

  #### prop scores ####
  pred <- predict(p, cov)
  pi <- get_probs(z,pred$z,pred$CDE)

  #### estimator ####
  tx = mu_a_xzplus - mu_a_xz; ty = mu_y_xzplus - mu_y_xz
  psihat = mean(ty)/mean(tx)

  v = (1/mean(tx)) * (var(ty) + psihat^2*var(tx) - 2*psihat*cov(tx,ty))
  print('variance:'); print(v)

  return(list(phi = psihat, var = v))
}

plugin2_ranger <- function(y,a,z,cov,ymean, amean, p, delta = 20){
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

  # #### prop scores ####
  # pred <- predict(p, cov)
  # pi <- get_probs(z,pred$z,pred$CDE)

  #### estimator ####
  tx = mu_a_xzplus - mu_a_xzmin; ty = mu_y_xzplus - mu_y_xzmin
  psihat = mean(ty)/mean(tx)

  v = (1/mean(tx)) * (var(ty) + psihat^2*var(tx) - 2*psihat*cov(tx,ty))
  print('variance:'); print(v)

  return(list(phi = psihat, var = v))
}
