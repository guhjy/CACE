
bundle<-function(dat,ranger,type,delta){
  # what kind of estimator do you want?

  if(ranger == T){
    if(type=='double'){phi_est = phi_est_double_ranger}
    if(type=='single'){phi_est = phi_est_single_ranger}
    if(type=='min'){phi_est = phi_est_min_ranger}
    if(type == 'plugin'){phi_est = plugin_ranger}
    if(type == 'pluginDouble'){phi_est = plugin2_ranger}
    if(type == 'simple'){phi_est = phi_est_single_simple; propscore_est = propscore_est_simple}
    if(type == 'simplePlugin'){phi_est = plugin_simple; propscore_est = propscore_est_simple}

    means = mean_est_ranger(y=dat[,1],a = dat[,2],z = dat[,3],cov = dat[,c(4:dim(dat)[2])])
    prop = propscore_est(y=dat[,3],x=dat[,c(4:dim(dat)[2])])
    out = phi_est(y=dat[,1],a = dat[,2],z = dat[,3],cov = dat[,c(4:dim(dat)[2])],ymean=ymean,amean=amean,p=prop,delta=delta)
  }
  else{
    if(type=='double'){phi_est = phi_est_double}
    if(type=='single'){phi_est = phi_est_single}
    if(type=='min'){phi_est = phi_est_min}
    if(type == 'plugin'){phi_est = plugin}
    if(type == 'pluginDouble'){print('Estimation not created yet')}
    if(type == 'simplePlugin'){print('Estimation not created yet')}
    if(type == 'simple'){phi_est = phi_est_single_simple; propscore_est = propscore_est_simple}

    means = mean_est(y=dat[,1],a = dat[,2],z = dat[,3],cov = dat[,c(4:dim(dat)[2])])

    if(type == 'plugIn' | type == 'plugInDouble'){prop = NULL}
    else{prop = propscore_est(y=dat[,3],x=dat[,c(4:dim(dat)[2])])}
    out = phi_est(y=dat[,1],a = dat[,2],z = dat[,3],cov = dat[,c(4:dim(dat)[2])],ymean=ymean,amean=amean,p=prop,delta=delta)
  }

  return(out)

}
