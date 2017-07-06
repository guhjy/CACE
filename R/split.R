#' Split data into 2
#'
#' @description Splits your dataset into 2 randomly
#'
#' @param dat your dataset
#'
#' @return creates two datasets in the global environment called ds1 and ds2
#'
#' @example

mysplit<-function(dat){
  mf = dat[complete.cases(dat),]
  n1 = round(.5*dim(mf)[1])
  indx = sample(1:dim(mf)[1])
  ds1 = mf[indx[1:n1],]
  ds2 = mf[indx[(n1+1):length(indx)],]
  assign("ds1",ds1, envir = .GlobalEnv)
  assign("ds2",ds2, envir = .GlobalEnv)
}


