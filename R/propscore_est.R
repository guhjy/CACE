#' Conditional Density Estimates
#'
#' @description Use FlexCoDE to estimate the conditional densities
#' that give the propensity scores
#'
#' @param dat is a dataset that is being used as your training
#' dataset overall (not within this algorithm)
#'
#' @return conditional density object


propscore_est <- function(dat){
  #load ann lee's conditional density estimation package
  library(digest)
  #devtools::install_github(repo = "rizbicki/FlexCoDE")
  library(FlexCoDE)

  attach(dat)
  #### Estimate pi(time | covariates)
  y <- total_time
  x <- as.data.frame(dat[,c(5:19,21)])
  n = dim(x)[1]

  # subsetting, not sure if this is right given i'm already only on one half of the data
  nTrain=round(0.7*n)
  nValidation=round(0.25*n)
  nTest=n-nTrain-nValidation

  # split data
  randomIndex=sample(1:n)
  xTrain=x[randomIndex[1:nTrain],]
  xValidation=x[randomIndex[(nTrain+1):(nTrain+nValidation)],]
  xTest=x[randomIndex[(nTrain+nValidation+1):n],]
  yTrain=y[randomIndex[1:nTrain]]
  yValidation=y[randomIndex[(nTrain+1):(nTrain+nValidation)]]
  yTest=y[randomIndex[(nTrain+nValidation+1):n]]

  # Fit nearest neighbors FlexCoDE
  fit=fitFlexCoDE(xTrain,yTrain,xValidation,yValidation,xTest,yTest,
                  nIMax = 30,regressionFunction = regressionFunction.NN)

  if(quite == F){
    fit$estimatedRisk
    print(fit)
    plot(fit,xTest,yTest)
  }
  detach(dat)
  return(fit)
}
