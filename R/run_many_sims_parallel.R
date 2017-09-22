run.all <- function(delta){
  #make k datasets (however many simulations you want at each delta)
  require('plyr')
  library(parallel)
  k = 10
  N.list <- lapply((1:k), function(x) 5000)
  data.sets <- lapply(N.list, makeSim)

  # get estimates of truth
  true.est <- EstimateTrue(delta)
  t.single <- true.est[1]; t.double <- true.est[2]
  
  out <- lapply(data.sets, function(x) getEstimatorBias(x, delta = delta, true.single = t.single, true.double = t.double))

  # # output on datasets for delta
  # clusterExport(cl, varlist = c("delta","t.single","t.double","getEstimatorBias"), envir=environment())
  # out <- parLapply(cl, data.sets, function(x) getEstimatorBias(x, delta = delta, true.single = t.single, true.double = t.double))
  df <- rbind.fill(out)

  return(df)
}
