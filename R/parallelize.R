# this is not doing what it needs to
# it can't seem to find packages that have been loaded elsewhere
# it can only download from master, obvi :,( i'm never making a branch again

library(parallel)
no_cores <- detectCores() - 1
cl <- makeCluster(no_cores)

dummyTest<-function(dat){
  require(CACE); require('ranger')
  CACE:::mean_est_ranger(y=dat[,1],a = dat[,2],z = dat[,3],cov = dat[,c(4:dim(dat)[2])])
  #CACE:::propscore_est(y=dat[,3],x=dat[,c(4:dim(dat)[2])])
  }

getMeanProp<-function(dat,ranger = T){
  require(CACE); require('ranger')
  if(ranger == F){
    # estimate means using what you've specified
    means  = CACE:::mean_est(y=dat[,1],a = dat[,2],z = dat[,3],cov = dat[,c(4:dim(dat)[2])])
  }
  else{
    # estimates neans using ranger
    means  = CACE:::mean_est_ranger(y=dat[,1],a = dat[,2],z = dat[,3],cov = dat[,c(4:dim(dat)[2])])
  }

  # get propensity score
  prop   = CACE:::propscore_est(y=dat[,3],x=dat[,c(4:dim(dat)[2])])

  return(list(means,prop))
}

start = proc.time()
test = parLapply(cl,list(ds1,ds2), dummyTest)
end = proc.time()


stopCluster(cl)
