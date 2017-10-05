a.mean.est <- function(a,dat,algorithm){
  if(tolower(algorithm) == 'glm'){
    df = as.data.frame(cbind(a,dat))
    out = glm(a~., data = df, family = 'binomial')
    }

  else if(tolower(algorithm) == 'superlearner'){
    require(SuperLearner)
    sl.lib <- c("SL.glm","SL.randomForest","SL.polymars","SL.mean")
    out = SuperLearner(Y=a, X=dat, SL.library=sl.lib, family=binomial())}

  else if(tolower(algorithm) == 'ranger'){
    require(ranger)
    out = ranger::ranger(a~.,data = as.data.frame(cbind(a,dat)), write.forest = T)
  }

  else{stop('Use ranger, superlearner or glm as algorithm')}

  return(out)
}
