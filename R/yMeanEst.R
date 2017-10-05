y.mean.est <- function(y,dat,algorithm){
  if(tolower(algorithm) == 'glm'){
    df = as.data.frame(cbind(y,dat))
    out = glm(y~., data = df, family = 'gaussian')
    }
  else if(tolower(algorithm) == 'superlearner'){
    require(SuperLearner)
    sl.lib <- c("SL.glm","SL.randomForest","SL.polymars","SL.mean")
    out = SuperLearner(Y=y, X=dat, SL.library=sl.lib, family=gaussian())}
  else if(tolower(algorithm) == 'ranger'){
    require(ranger)
    out = ranger::ranger(y~.,data = as.data.frame(cbind(y,dat)), write.forest = T)
  }

  else{stop('Use ranger, superlearner or glm as algorithm')}

  return(out)
}
