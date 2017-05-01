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
#'
#' @return an estimate of the causal effect

phi_est <- function(dat,m,p){
  attach(dat)

  z <- total_time
  delta = 2
  zplus <- z + delta
  zmin <- z - delta

  ymean = m[[1]]
  amean = m[[2]]

  #### Y means ####
  # predicted means of y|x,z
  xnew <-  as.data.frame(dat[,c(2,5:19,21)])
  mu_y_xz <- predict(ymean,newdata = xnew)$pred

  # predicted means of y|x,z+delta
  xnew <- as.data.frame(cbind(zplus,dat[,c(5:19,21)]))
  mu_y_xzplus <- predict(ymean,newdata = xnew)$pred

  # predicted means of y|x,z-delta
  xnew <- as.data.frame(cbind(zmin,dat[,c(5:19,21)]))
  mu_y_xzplus <- predict(ymean,newdata = xnew)$pred


  #### A means ####
  # predicted means of a|x,z
  xnew <-  as.data.frame(dat[,c(2,5:19,21)])
  mu_a_xz <- predict(amean,newdata = xnew)$pred
  # other.preds <- predict(ymean,newdata = xnew)$library.predict

  # predicted means of a|x,z+delta
  xnew <- as.data.frame(cbind(zplus,dat[,c(5:19,21)]))
  mu_a_xzplus <- predict(amean,newdata = xnew)$pred

  # predicted means of a|x,z-delta
  xnew <- as.data.frame(cbind(zmin,dat[,c(5:19,21)]))
  mu_a_xzplus <- predict(amean,newdata = xnew)$pred

  #### prop scores ####
  xnew <-  as.data.frame(dat[,c(5:19,21)])
  pred_pi <- predict(p, xNew = xnew)$CDE

  #need to then get: pi = prob(Z = z); pi_min = prob(Z = z-delta); pi_plus = prob(Z = z+delta)

  #### estimator ####
  phi_top = (NCRecid3 - mu_y_xz)(pi_min - pi_plus)/pi + mu_y_xzplus - mu_y_xzmin
  phi_bot = (visiteveryn - mu_a_xz)(pi_min - pi_plus)/pi + mu_a_xzplus - mu_a_xzmin
  phi = phi_top/phi_bot
}
