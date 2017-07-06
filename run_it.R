library(readxl)
setwd("C:/Users/jackie/Dropbox/prisoners")
vis<- read_excel("Visit Study - All Supervision.xlsx")
setwd("C:/Users/jackie/CACE")

#get covariate dataframe
t = apply(vis[,31:33],1,sum) #these release codes sum to one, leave out one
vis$commit_cnty <- as.factor(vis$commit_cnty)
vis$CountyClass <- as.factor(vis$CountyClass)
vis$servedLen = as.numeric(as.Date(as.character(vis$mov_move_date2008),"%Y%m%d") - as.Date(as.character(vis$admit),"%Y%m%d"))
mat<-model.matrix(NCRecid3~-1+visitslastlocyn1+total_time+loslastloc+white+urban+priorarrests+married+violent
                  +lsirscore+ageyrs+custody_level+numofpriorinc+prerelease+parole
                  +releasetocenter+mh+highschoolgrad+facility+numofpriormisconducts
                  +servedLen+commit_cnty+CountyClass,
                  data = vis)

# can't figure out predictions on ranger
y = vis$NCRecid3[as.numeric(rownames(mat))]
phi = CACE(y=y,a=mat[,1],z=mat[,2],cov=mat[,3:dim(mat)[2]],delta = 20)

vals = c(20,60,120,240)
out = sapply(vals,function(x) CACE(y=y,a=mat[,1],z=mat[,2],cov=mat[,3:dim(mat)[2]],delta = x))


