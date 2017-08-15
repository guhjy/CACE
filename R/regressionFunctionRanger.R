# trying to make a ranger function
regressionFunction.Ranger<-function(x,responses,extra=NULL){
  #x a matrix of covariates used to train the model
  #responses matrix where each row corresponds to a row in x and each column corresponds to a different response (in practice, it will be $\phi_i(z)$, but you do not need to worry about this)

  #responses come in in a 7*n x 40 matrix--need to figure out what to do with that
  require("ranger")
  require(doParallel)
  #dat = data.frame(y = responses[,1],x = x)
  #fit.rf <- ranger::ranger(y~ ., data = dat, write.forest = T)
  #class(fit.rf) <- c('rangerReg')
  #return(fit.rf)

  # Both x and responses are matrices
  n=dim(x)[1]
  random=sample(1:n)
  nTrain=round(0.7*n)
  xTrain=x[random[1:nTrain],,drop=FALSE]
  responsesTrain=responses[random[1:nTrain],,drop=FALSE]
  xValidation=x[random[-c(1:nTrain)],,drop=FALSE]
  responsesValidation=responses[random[-c(1:nTrain)],,drop=FALSE]

  p0Vec=extra$p0Vec
  if(is.null(p0Vec))
    p0Vec=round(seq(1,ncol(xTrain),length.out = 5))

  ntree=extra$ntree
  if(is.null(ntree))
    ntree=500

  maxnodes=extra$maxnodes

  nCores=extra$nCores
  if(is.null(nCores))
    nCores=1

  cl <- parallel::makeCluster(nCores)
  doParallel::registerDoParallel(cl)


  fittedReg <- foreach(ii=1:ncol(responsesTrain), .verbose = T) %dopar% {
  #fittedReg<-for(ii in 2:ncol(responsesTrain)){
    error=rep(NA,length(p0Vec))
    dat = data.frame(y = responsesTrain[,ii,drop=FALSE],x = xTrain)
    newdat = data.frame(y= responsesValidation[,ii,drop = FALSE],x = xValidation)
    for(s in 1:length(p0Vec))
    {
      ajuste = ranger::ranger(y~.,data = dat,write.forest = T, mtry=p0Vec[s])
      #ajuste = randomForest::randomForest(x=xTrain,y=responsesTrain[,ii,drop=FALSE],mtry=p0Vec[s],importance = FALSE)
      predito = predict(ajuste, data = newdat)$predictions
      error[s]=mean((predito-responsesValidation[,ii,drop=FALSE])^2)
    }
    bestP0=p0Vec[which.min(error)]
    #ajuste = randomForest(x=xTrain,y=responsesTrain[,ii,drop=FALSE],mtry=bestP0,importance = TRUE)
    ajuste = ranger::ranger(y~.,data = dat,mtry=bestP0)
    #return(ajuste)
    object=NULL
    object$fit=ajuste
    object$importance=ajuste$importance
    object$bestP0=bestP0
    object$errors=ajuste$mse
    gc(verbose=FALSE)
    return(object)
  }

  parallel::stopCluster(cl)

  regressionObject=NULL
  regressionObject$fittedReg=fittedReg
  class(regressionObject)="FC.ranger"
  return(regressionObject)
}

predict.Ranger<-function(object,xNew,maxTerms=NULL){
    if(class(object)!="FC.ranger")
      stop("Object has wrong class, should be FC.ranger")

  if(!is.null(maxTerms))
  {
    maxTerms=min(maxTerms,length(object$fittedReg))
  } else {
    maxTerms=length(object$fittedReg)
  }

  predictedValidation=apply(as.matrix(1:maxTerms),1,function(xx)
  {
    colnames(xNew)=NULL
    predicted = predict(object$fittedReg[[xx]]$fit,data=xNew)$predictions
    return(predicted)
  })

  return(predictedValidation)
}

print.Ranger<-function(regressionObject,bestI,nameCovariates){
  cat(paste('this is nothing'))
}
