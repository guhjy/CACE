# preds on new dataset
# detach(ds1)
# attach(ds2)

phi_est <- function(dat){
  attach(dat)

  z <- total_time
  delta = 2
  zplus <- z + delta
  zmin <- z - delta

  #### Y means ####
  # predicted means of y|x,z
  xnew <-  as.data.frame(dat[,c(2,5:19,21)])
  mu_y_xz <- predict(ymean,newdata = xnew)$pred
  # other.preds <- predict(ymean,newdata = xnew)$library.predict

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
  pred_pi <- predict(fit, xNew = xnew)$CDE

  #### estimator ####
  phi_top = (NCRecid3 - mu_y_xz)(pi_min - pi_plus)/pi + mu_y_xzplus - mu_y_xzmin
  phi_bot = (visiteveryn - mu_a_xz)(pi_min - pi_plus)/pi + mu_a_xzplus - mu_a_xzmin
  phi = phi_top/phi_bot
}
