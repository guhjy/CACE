#' Load data
#'
#' @description This probably won't do anyone else any good
#'
#' @return list of my two datasets

loadem <- function(){
  library(readxl)
  setwd("C:/Users/jackie/Dropbox/prisoners")
  vis <- read_excel("Visit Study - All Supervision.xlsx")
  move <- read_excel("Movements Prior to 2008 Release.xlsx")
  setwd("C:/Users/jackie/CACE")
  return(list(move,vis))
}

