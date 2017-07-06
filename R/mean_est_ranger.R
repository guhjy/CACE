#' Regression Function Estimates using Ranger
#'
#' @description Use ranger to estimate the mean
#' functions for the outcome and the treatment given covariates
#' and the instrument
#'
#' @param y the outcome variable in your training set
#' @param a the treatment variable in your training set
#' @param z the instrument in your training set
#' @param cov the covariates in your training set. Should be a
#' data frame and it should not include any factors for now.
#'
#' @return gives you ymean and amean in the global envr. Note that I
#' tried to return these as a list and it gave me something weird, but
#' this works.

mean_est_ranger <- function(y,a,z,cov){
  require("ranger")
  print("Estimating Means using Ranger")
  ptm <- proc.time()

  ##### Estimate E( recid | distance, cov ) #####
  # set up data frame:
  dat1 <- as.data.frame(cbind(y,z,cov))
  dat2 <- as.data.frame(cbind(a,z,cov))

  # run the super learner for y
  (rang.res1 <- ranger::ranger(y~.,data = dat1,write.forest = T))

  # run the super learner for a
  (rang.res2 <- ranger::ranger(a~.,data = dat2,write.forest = T))

  #return(list(ymean = sl.res1,amean = sl.res2))
  assign("ymean",rang.res1, envir = .GlobalEnv)
  assign("amean",rang.res2, envir = .GlobalEnv)
  print(paste("Mean estimation runtime:",(proc.time()-ptm)[1]))
}
