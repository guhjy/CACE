# kang and schafer simulated data

expit <- function(x){exp(x)/(1+exp(x))}
logit <- function(p){log(p/(1-p))}

ks_data<-function(n){
  z <- data.frame(matrix(rnorm(4*n),nrow=n,ncol=4)); colnames(z) <- paste("z",1:4,sep="")
  x <- cbind( exp(z[,1]/2) , 10+z[,2]/(1+exp(z[,1])) , (z[,1]*z[,3]/25 + 0.6)^3 , (z[,2]+z[,4]+20)^2)
  x <- as.data.frame(x); colnames(x) <- paste("x",1:4,sep="")
  pi <- 1-expit(-z[,1]+0.5*z[,2]-0.25*z[,3]-0.1*z[,4]); r <- rbinom(n,1,pi)
  mu <- 210+27.4*z[,1]+13.7*z[,2]+13.7*z[,3]+13.7*z[,4]; y <- mu + rnorm(n); y[r==0] <- 0
  dat <- data.frame(x,z,r,y); colnames(dat)[1:8] <- c(paste("x",1:4,sep=""),paste("z",1:4,sep=""))
  dat
}


