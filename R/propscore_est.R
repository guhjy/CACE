#' Conditional Density Estimates
#'
#' @description Use FlexCoDE to estimate the conditional densities
#' that give the propensity scores
#'
#' @param y is your outcome variable (here it should be your instrument)
#' @param x is a dataframe of covariates
#'
#' @return conditional density object


propscore_est <- function(y, x, quiet = T,regFunc = regressionFunction.NW){
  print("Estimating Propensity Scores")
  ptm <- proc.time()

  #load ann lee's conditional density estimation package
  library(digest)
  library(FlexCoDE)

  #### Estimate pi(time | covariates)

  # subsetting, not sure if this is right given i'm already only on one half of the data
  x = as.matrix(x)
  n = dim(x)[1]
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
                  nIMax = 40,regressionFunction = regFunc,
                  regressionFunction.extra=list(nCores=3),verbose = T)

  if(quiet == F){
    fit$estimatedRisk
    print(fit)
    plot(fit,xTest,yTest)
  }
  print(paste("Propensity score estimation runtime:",(proc.time()-ptm)[1]))
  return(fit)
}
