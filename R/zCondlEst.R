z.condldens.est <- function(z,x,algorithm,regFunc = regressionFunction.NW){
  if(tolower(algorithm) == 'glm'){
    df = as.data.frame(cbind(z,x))
    out <- glm(z~., data = df, family = 'gaussian')
    }

  else if(tolower(algorithm) == 'flexcode'){
    require(FlexCoDE)
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
    zTrain=z[randomIndex[1:nTrain]]
    zValidation=z[randomIndex[(nTrain+1):(nTrain+nValidation)]]
    zTest=z[randomIndex[(nTrain+nValidation+1):n]]

    # Fit nearest neighbors FlexCoDE
    out=fitFlexCoDE(xTrain,zTrain,xValidation,zValidation,xTest,zTest,
                    nIMax = 40,regressionFunction = regFunc,
                    regressionFunction.extra=list(nCores=3),verbose = T)
  }

  else{stop('Use flexcode or glm as algorithm')}

  return(out)
}
