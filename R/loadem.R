#' Load data
#'
#' @description This probably won't do anyone else any good
#'
#' @return list of my two datasets

loadem <- function(){
  library(readxl)
  setwd("C:/Users/jackie/Dropbox/prisoners")
  vis_all <- read_excel("Visit Study - All Supervision.xlsx")
  vis <- vis_all[,c('NCRecid3','total_time','visitseveryn',
                'facility','loslastloc','white','male','urban',
                'priorarrests','married','violent','lsirscore',
                'ageyrs','custody_level','numofpriorinc','mh',
                'highschoolgrad','numofpriormisconducts','numoftotalmisconducts',
                'commit_cnty','CountyClass')]
  mov <- read_excel("Movements Prior to 2008 Release.xlsx")
  setwd("C:/Users/jackie/CACE")
  assign("vis",vis, envir = .GlobalEnv)
  assign("move",mov, envir = .GlobalEnv)
  #return(list(move,vis))
}

