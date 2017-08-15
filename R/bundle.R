bundle<-function(dat,ranger,type,delta){
  # what kind of estimator do you want?

  if(ranger == T){
    if(type=='double'){phi_est = phi_est_double_ranger}
    if(type=='single'){phi_est = phi_est_single_ranger}
    if(type=='min'){phi_est = phi_est_min_ranger}
    if(type == 'plugin'){phi_est = plugin_ranger}

    means = mean_est_ranger(y=dat[,1],a = dat[,2],z = dat[,3],cov = dat[,c(4:dim(dat)[2])])
    prop = propscore_est(y=dat[,3],x=dat[,c(4:dim(dat)[2])])
    out = phi_est(y=dat[,1],a = dat[,2],z = dat[,3],cov = dat[,c(4:dim(dat)[2])],ymean=ymean,amean=amean,p=prop,delta=delta)
  }
  else{
    #prediction not done for these yet
    if(type=='double'){phi_est = phi_est_double}
    if(type=='single'){phi_est = phi_est_single}
    if(type=='min'){phi_est = phi_est_min}
    if(type == 'plugin'){phi_est = plugin}

    means = mean_est(y=dat[,1],a = dat[,2],z = dat[,3],cov = dat[,c(4:dim(dat)[2])])
    prop = propscore_est(y=dat[,3],x=dat[,c(4:dim(dat)[2])])
    out = phi_est(y=dat[,1],a = dat[,2],z = dat[,3],cov = dat[,c(4:dim(dat)[2])],ymean=ymean,amean=amean,p=prop,delta=delta)
  }

  return(out)

}
