# ed sim
require(truncnorm)

N = 10000
sig2y0 = 1;
y0 = rnorm(N,0,sig2y0); meanx = c(0,0,0,0)
alpha = matrix(c(1,1,-1,-1),nrow = 4); beta = matrix(c(1,-1,-1,1)); psi =1
delta = 5

x = t(matrix(unlist(lapply(meanx,function(x) rnorm(N,x,1))),nrow = N,byrow =T))
z = rtruncnorm(N,a=-2,b=2,mean=1.5*sign(t(alpha)%*%x),sd = 2)
t = rnorm(N,t(beta)%*%x+y0,1)
a = as.numeric(z>=t); aplus = as.numeric(z+delta>=t)
y = y0 + a*psi*t; yplus = y0 + aplus*psi*t

num_compliers = length(which((z+delta)>=t))-length(which(z>=t)) #monotone in delta
Z = pnorm(-2)-pnorm(2)
meanZ = 1.5*sign(t(alpha)%*%x) + 2*(dnorm(-2)-dnorm(2))/Z

muA = pnorm((z-t(beta)%*%x)/sqrt(2));muAplus = pnorm(((z+delta)-t(beta)%*%x)/sqrt(2)); muhatA = mean(a)
muY = c(t(beta)%*%x)*c(pnorm((z-t(beta)%*%x)/sqrt(2))) - sqrt(2)*c(dnorm((z-t(beta)%*%x)));muYplus = c(t(beta)%*%x)*c(pnorm(((z+delta)-t(beta)%*%x)/sqrt(2))) - sqrt(2)*c(dnorm(((z+delta)-t(beta)%*%x))); muhatY = mean(y)
pihat = mean(z>=t); piminhat = mean((z+delta)>=t)
pi = dnorm((z-1.5*sign(t(alpha)%*%x))/2)/(2*Z); pi_min =  dnorm(((z+delta)-1.5*sign(t(alpha)%*%x))/2)/(2*Z)
phi_y = (y - muY)*(pi_min/pi) - (y - muYplus)
phi_a = (a - muA)*(pi_min/pi) - (a - muAplus)

LATEest1 = mean(muYplus - muY)/mean(muAplus - muA)
LATEest2 = mean(phi_y)/mean(phi_a)
LATEest3 = mean(t[which((z+delta)>=t & t>z)])
est = CACE(y=y,a=a,z=z,cov = as.data.frame(t(x)),delta = delta,ranger = T,type = 'single')

temp = cbind(c('delta','bottom','a-mua','a-muaplus','y-muy','y-muyplus','pimin','pi'),temp)
goodtemp = cbind(c('delta','bottom','a-mua','a-muaplus','y-muy','y-muyplus','pimin','pi'),
                 c(delta,mean(phi_a),mean(a - muA),mean(a - muAplus),mean(y-muY),mean(y-muYplus),mean(pi_min),mean(pi)))

f<- function(d){CACE(y=y,a=a,z=z,cov = as.data.frame(t(x)),delta = d,ranger = T,type = 'single')}
f2 <- function(d){mean(t[which((z+d)>=t & t>z)])}
trues = unlist(lapply(c(0:10),f2))
outRange = lapply(c(0:10),f)
outMat = as.data.frame(matrix(unlist(outRange),ncol = 4,byrow = T))
outMat = cbind(c(0:10),outMat,trues)
names(outMat)<- c('delta','estimate','sd','numerator','denominator','truth')
write.csv(outMat, file = "EdSimForest.csv")

plot(outMat$delta,outMat$truth, xlab = "delta", ylab = "estimate", ylim  = c(min(c(outMat$truth,outMat$estimate),na.rm = T),max(c(outMat$truth,outMat$estimate),na.rm = T)),main = "Regression.NW estimates on Ed simulation")
points(outMat$delta,outMat$estimate,pch = 19, col = 'red')
error.bar(outMat$delta,outMat$estimate,2*outMat$sd)