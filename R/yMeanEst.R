y.mean.est <- function(y,dat,algorithm,sl.lib = c("SL.glm","SL.randomForest","SL.polymars","SL.mean")){
  data = as.data.frame(cbind(y,dat))

  if(tolower(algorithm) == 'glm'){
    out = glm(y~., data = data, family = 'gaussian')
    }

  else if(tolower(algorithm) == 'superlearner'){
    require(SuperLearner)
    out = SuperLearner(Y=y, X=dat, SL.library=sl.lib, family=gaussian())
    }

  else if(tolower(algorithm) == 'ranger'){
    require(ranger)
    out = ranger::ranger(y~.,data = data, write.forest = T)
  }

  else if(tolower(algorithm) == 'random forest'){
    require(randomForest)
    out = randomForest(y~., data = data)
  }

  else{stop('Use ranger, superlearner or glm as algorithm')}

  return(out)
}
