# this is not doing what it needs to
# it can't seem to find packages that have been loaded elsewhere

# library(parallel)
# no_cores <- detectCores() - 1
# cl <- makeCluster(no_cores)
#
# dummyTest<-function(dat){
#   require(CACE)
#   CACE:::mean_est(y=dat[,1],a = dat[,2],z = dat[,3],cov = dat[,c(4:dim(dat)[2])])
# }
#
# getMeanProp<-function(dat,ranger = ranger){
#   if(ranger == F){
#     # estimate means using what you've specified
#     means  = mean_est(y=dat[,1],a = dat[,2],z = dat[,3],cov = dat[,c(4:dim(dat)[2])])
#   }
#   else{
#     # estimates neans using ranger
#     means  = mean_est_ranger(y=dat[,1],a = dat[,2],z = dat[,3],cov = dat[,c(4:dim(dat)[2])])
#   }
#
#   # get propensity score
#   prop   = propscore_est(y=dat[,3],x=dat[,c(4:dim(dat)[2])])
#
#   return(list(means,prop))
# }
#
# test = parLapply(cl,list(ds1,ds2), dummyTest)
#
#
# stopCluster(cl)
