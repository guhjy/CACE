#' Regression Function Estimates
#'
#' @description Use a SuperLearner to estimate the mean
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

mean_est <- function(y,a,z,cov){
  print("Estimating Means")
  ptm <- proc.time()
  library(SuperLearner)
  set.seed(87932)

  ##### Estimate E( recid | distance, cov ) #####
  # set up data frame:
  x <- as.data.frame(cbind(z,cov))

  # pick which algorithms you want
  sl.lib <- c("SL.glm"
              #,
              #"SL.ranger"
              ,"SL.randomForest"
              ,"SL.polymars"
              #,"SL.gam2"
              #,"SL.gam3"
              #,"SL.gam"
              #,"SL.glmnet"
              ,"SL.mean"
              #,"SL.knn"
              #,"SL.loess"
              #,"SL.leekasso"
  )

  # run the super learner for y
  (sl.res1 <- SuperLearner(Y=y, X=x, SL.library=sl.lib, family=gaussian()))

  # run the super learner for a
  (sl.res2 <- SuperLearner(Y=a, X=x, SL.library=sl.lib, family=binomial()))

  #return(list(ymean = sl.res1,amean = sl.res2))
  assign("ymean",sl.res1, envir = .GlobalEnv)
  assign("amean",sl.res2, envir = .GlobalEnv)
  print(paste("Mean estimation runtime:",(proc.time()-ptm)[1]))
}
