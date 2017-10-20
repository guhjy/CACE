# helper functions for kde predictions

library(hdrcde)

get.xrank <- function(x,reg,nx){
  out = rep(NA,length(x))
  for(i in 1:length(x)){
    out[i] = min(rank(append(x[i],reg$x[[i]]),ties.method = 'last')[1],nx)
  }
  return(out)
}
format.range <- function(x){
  if(length(x)>1){
    out = paste(x[1],x[2], sep = ":")
  }
  else{out = x}
  return(out)
}

get.pihat <- function(z,xmat,mean.est,sd.est,nxmargin = 100,maxk){
  Zst = as.numeric((z - mean.est)/sd.est)
  pi.reg = cde(xmat,Zst,nxmargin = nxmargin, maxk = maxk)

  pi.est = rep(NA,length(Zst))
  for(i in 1:length(Zst)){
    xrank = get.xrank(xmat[i,], pi.reg, dim(pi.reg$z)[1])
    yrank = min(rank(append(Zst[i],pi.reg$y),ties.method = 'last')[1],length(pi.reg$y))

    xrange = lapply(xrank, function(x) c(max(x-1,1):x))
    yrange = c(max(yrank-1, 1):yrank)
    temprange = xrange
    temprange[[length(temprange)+1]] = yrange
    range = Reduce(c,lapply(temprange, format.range))
    pi.est[i] <- mean(eval(parse(text =paste('pi.reg$z[',paste(range, collapse = ','),']', sep = ""))))
  }

  return(pi.est)
}

