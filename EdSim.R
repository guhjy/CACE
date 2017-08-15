#### ed sim ####
rm(list = ls())
require(truncnorm)
sig2y0 = 1;N = 10000
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

dat = as.data.frame(cbind(y,a,z,t,t(x)))
train = sample(1:dim(dat)[1],8000)
datTrain = dat[train,]; datTest = dat[-train,]

#### true parameters ####
muA = pnorm((z-t(beta)%*%x)/sqrt(2));muAplus = pnorm(((z+delta)-t(beta)%*%x)/sqrt(2)); muhatA = mean(a)
muY = c(t(beta)%*%x)*c(pnorm((z-t(beta)%*%x)/sqrt(2))) - sqrt(2)*c(dnorm((z-t(beta)%*%x)));muYplus = c(t(beta)%*%x)*c(pnorm(((z+delta)-t(beta)%*%x)/sqrt(2))) - sqrt(2)*c(dnorm(((z+delta)-t(beta)%*%x))); muhatY = mean(y)
pihat = mean(z>=t); piminhat = mean((z+delta)>=t)
pi = dnorm((z-1.5*sign(t(alpha)%*%x))/2)/(2*Z); pi_min =  dnorm(((z+delta)-1.5*sign(t(alpha)%*%x))/2)/(2*Z)
phi_y = (y - muY)*(pi_min/pi) - (y - muYplus)
phi_a = (a - muA)*(pi_min/pi) - (a - muAplus)

LATEest1 = mean(muYplus - muY)/mean(muAplus - muA) #plug in
LATEest2 = mean(phi_y)/mean(phi_a) #infl-fn based
LATEest3 = mean(t[which((z+delta)>=t & t>z)])

####check mean functions####
delta = 2
#using Ranger
mnsTrain = mean_est_ranger(y=datTrain$y,a=datTrain$a,z=datTrain$z,cov=datTrain[,c(5:8)])
ymean1 <- ymean; amean1 <- amean
predDat = data.frame(z = datTest$z,V5 = datTest$V5,V6 = datTest$V6,V7 = datTest$V7,V8 = datTest$V8)
predDatplus = data.frame(z = (datTest$z+delta),V5 = datTest$V5,V6 = datTest$V6,V7 = datTest$V7,V8 = datTest$V8)
yPred = predict(ymean1,predDat);aPred = predict(amean1,predDat)
yPredplus = predict(ymean1,predDatplus); aPredplus = predict(amean1,predDatplus)

par(mfrow = c(2,2),mar = c(4,4,2.5,1))
plot(yPred[[1]],datTest$y,xlim=c(min(yPred[[1]],datTest$y),max(yPred[[1]],datTest$y)),xlab = 'predicted Y',ylab = 'true Y')
plot(aPred[[1]],datTest$a,xlim = c(min(aPred[[1]],datTest$a),max(aPred[[1]],datTest$a)),xlab = 'predicted A',ylab = 'true A')
plot(yPredplus[[1]],yplus[-train],xlim=c(min(yPredplus[[1]],yplus[-train]),max(yPredplus[[1]],yplus[-train])),xlab = 'predicted Y plus',ylab = 'true Y plus')
plot(aPredplus[[1]],aplus[-train],xlim = c(min(aPredplus[[1]],aplus[-train]),max(aPredplus[[1]],aplus[-train])),xlab = 'predicted A plus',ylab = 'true A plus')
mtext("Ranger Estimates", side = 1, line = -21, outer = TRUE)

#using SuperLearner
mnsTrain2 = mean_est(y=datTrain$y,a=datTrain$a,z=datTrain$z,cov=datTrain[,c(5:8)])
ymean2 <- ymean; amean2 <- amean
predDat = data.frame(z = datTest$z,V5 = datTest$V5,V6 = datTest$V6,V7 = datTest$V7,V8 = datTest$V8)
predDatplus = data.frame(z = (datTest$z+delta),V5 = datTest$V5,V6 = datTest$V6,V7 = datTest$V7,V8 = datTest$V8)
yPred2 = predict(ymean2,predDat); aPred2 = predict(amean2,predDat)
yPredplus2 = predict(ymean2,predDatplus); aPredplus2 = predict(amean2,predDatplus)

plot(yPred2[[1]],datTest$y,xlim=c(min(yPred2[[1]],datTest$y),max(yPred2[[1]],datTest$y)),xlab = 'predicted Y',ylab = 'true Y')
plot(aPred2[[1]],datTest$a,xlim = c(min(aPred2[[1]],datTest$a),max(aPred2[[1]],datTest$a)),xlab = 'predicted A',ylab = 'true A')
plot(yPredplus2[[1]],yplus[-train],xlim=c(min(yPredplus2[[1]],yplus[-train]),max(yPredplus2[[1]],yplus[-train])),xlab = 'predicted Y plus',ylab = 'true Y plus')
plot(aPredplus2[[1]],aplus[-train],xlim =c(min(aPredplus2[[1]],aplus[-train]),max(aPredplus2[[1]],aplus[-train])),xlab = 'predicted A plus',ylab = 'true A plus')
mtext("SuperLearner Estimates", side = 1, line = -21, outer = TRUE)
par(mfrow = c(1,1))

#### check propscore functions--not good ####
truePi = dtruncnorm(datTest$z,a=-2,b=2,mean=1.5*sign(t(alpha)%*%x[,-train]),sd = 2)
propscoreTrain = propscore_est(y = datTrain$z,x = datTrain[,c(5:8)])
predcovs = data.frame(V5 = datTest$V5,V6 = datTest$V6,V7 = datTest$V7,V8 = datTest$V8)
predZ = predict(propscoreTrain,predcovs)
predZactual = get_probs(datTest$z,predZ$z,predZ$CDE)
plot(truePi,predZactual,xlab = 'real data',ylab ='predictions',main = 'FlexCode prediction (NW)'
     ,xlim = c(min(truePi,predZactual),max(truePi,predZactual))
     ,ylim = c(min(truePi,predZactual),max(truePi,predZactual))
     )

propscoreTrain = propscore_est(y = datTrain$z,x = datTrain[,c(5:8)],regFunc = regressionFunction.SpAM)
predcovs = as.matrix(datTest[,5:8])
predZ = predict(propscoreTrain,predcovs)
predZactual = get_probs(datTest$z,predZ$z,predZ$CDE)
plot(truePi,predZactual,xlab = 'real data',ylab ='predictions',main = 'FlexCode prediction (SpAM)'
     ,xlim = c(min(truePi,predZactual),max(truePi,predZactual))
     ,ylim = c(min(truePi,predZactual),max(truePi,predZactual))
)
propscoreTrain = propscore_est(y = datTrain$z,x = datTrain[,c(5:8)],regFunc = regressionFunction.Forest)
predcovs = as.matrix(datTest[,5:8])
predZ = predict(propscoreTrain,predcovs)
predZactual = get_probs(datTest$z,predZ$z,predZ$CDE)
plot(truePi,predZactual,xlab = 'real data',ylab ='predictions',main = 'FlexCode prediction (Forest)'
     ,xlim = c(min(truePi,predZactual),max(truePi,predZactual))
     ,ylim = c(min(truePi,predZactual),max(truePi,predZactual))
)


#### check phi ####
est = CACE(y=y,a=a,z=z,cov = as.data.frame(t(x)),delta = delta,ranger = T,type = 'single',split = F)

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

lm = lm(y~a+z+V5+V6+V7+V8,data = datTrain)
predLM = predict(lm,datTest)
plot(predLM~datTest$y)
