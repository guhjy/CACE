SL.ranger <- function (Y, X, newX, family, ...) {
  require("ranger")
  fit.rf <- ranger::ranger(Y ~ ., data=X,write.forest = T)
  pred <- predict(fit.rf,data=newX)$predictions
  fit <- list(object = fit.rf)
  out <- list(pred = pred, fit = fit)
  class(out$fit) <- c("SL.ranger")
  return(out)
}
