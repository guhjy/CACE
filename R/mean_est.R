#' Regression Function Estimates
#'
#' @description Use a SuperLearner to estimate the mean
#' functions for the outcome and the treatment given covariates
#' and the instrument
#'
#' @param dat is a dataset that is being used as your training
#' dataset overall (not within this algorithm)
#'
#' @return list of the means for outcome ('ymean') and treatment ('amean')

mean_est <- function(dat){
  library(SuperLearner)
  set.seed(87932)
  attach(dat)

  ##### Estimate E( recid | distance, cov ) #####
  # set up data frame:
  # outcome: recidivate w/in three years (0/1)
  # instrument: total drive time (real number)
  y <- NCRecid3
  x <- as.data.frame(dat[,c(2,5:19,21)])

  # don't include any factors
  # x$facility = as.factor(x$facility)
  # x$commit_cnty = as.factor(x$commit_cnty)

  # pick which algorithms you want
  sl.lib <- c("SL.gam"
              ,"SL.polymars"
              #,"SL.gam2"
              #,"SL.gam3"
              ,"SL.glm"
              ,"SL.glmnet"
              ,"SL.mean"
              ,"SL.knn"
              #,"SL.loess"
              #,"SL.leekasso"
  )

  # run the super learner
  (sl.res1 <- SuperLearner(Y=y, X=x, SL.library=sl.lib, family=binomial()))


  ##### Estimate E( visit | distance, cov ) #####
  # outcome: ever visited (0/1)
  # instrument: total drive time (real number)
  y <- visitseveryn
  x <- as.data.frame(dat[,c(2,5:19,21)])

  (sl.res2 <- SuperLearner(Y=y, X=x, SL.library=sl.lib, family=binomial()))
  # sl.preds2 <- predict(sl.res2)$pred
  # other.preds2 <- predict(sl.res2)$library.predict
  detach(dat)
  return(list(ymean = sl.res1,amean = sl.res2))
}

#
#
# library(SuperLearner)
# # listWrappers() #<- if you want to see what library options you have
# set.seed(87932)
# # loadem()
# # mysplit(vis)
# attach(ds1)
#
# ##### Estimate E( recid | distance, cov ) #####
# # set up data frame:
# # outcome: recidivate w/in three years (0/1)
# # instrument: total drive time (real number)
# y <- NCRecid3
# x <- as.data.frame(ds1[,c(2,5:19,21)])
#
# # don't include any factors
# # x$facility = as.factor(x$facility)
# # x$commit_cnty = as.factor(x$commit_cnty)
#
# # pick which algorithms you want
# sl.lib <- c("SL.gam"
#             #,"SL.gam2"
#             #,"SL.gam3"
#             ,"SL.glm"
#             ,"SL.glmnet"
#             ,"SL.mean"
#             ,"SL.knn"
#             #,"SL.loess"
#             #,"SL.leekasso"
#             )
#
# # run the super learner
# (sl.res1 <- SuperLearner(Y=y, X=x, SL.library=sl.lib, family=binomial()))
#
#
# ##### Estimate E( visit | distance, cov ) #####
# # outcome: ever visited (0/1)
# # instrument: total drive time (real number)
# y <- visitseveryn
# x <- as.data.frame(cbind(as.factor(facility),
#                          loslastloc,white,male,urban,priorarrests,
#                          married,violent,lsirscore,ageyrs,custody_level,
#                          numofpriorinc,mh,highschoolgrad,
#                          numofpriormisconducts,numoftotalmisconducts,
#                          as.factor(commit_cnty),CountyClass
#                          ,total_time))
#
# ok <- complete.cases(x,y)
# x <- x[ok,]
# y <- y[ok]
#
# # pick which algorithms you want
# sl.lib <- c("SL.gam"
#             #,"SL.gam2"
#             #,"SL.gam3"
#             ,"SL.glm"
#             ,"SL.glmnet"
#             ,"SL.mean"
#             ,"SL.knn"
#             #,"SL.loess"
#             #,"SL.leekasso"
# )
#
# (sl.res2 <- SuperLearner(Y=y, X=x, SL.library=sl.lib, family=binomial()))
# # sl.preds2 <- predict(sl.res2)$pred
# # other.preds2 <- predict(sl.res2)$library.predict
