logit <- function(x) {exp(x)/(1+exp(x))}
N = 6000
x1 = rnorm(N,3,1)
x2 = runif(N)
x3 = rbeta(N,5,100)
cov = data.frame(x1 = x1, x2 = x2, x3 = x3)
z = rbeta(N,.5*abs(x1) + .2*x2 + .1*x3,10)

delta = mean(z)/2 #shift
thr = .1 #threshold (higher ==> more compliers)

# need to set up a fn so that aplus>amin (no defiers)
f = rbinom(N,1,logit(.3*x3/(1+x2) - .2*x1))
a = ifelse(f==0, f + as.numeric(z>thr), f)
aplus = ifelse(f==0, f + as.numeric((z+delta)>thr), f)
amin  = ifelse(f==0, f + as.numeric((z-delta)>thr), f)
#mean(aplus>amin); mean(aplus<amin)

f2 = a*(.3*x1 + .4*x2 - x1) + (1-a)*(x1-x2+2*x1)
f2p = aplus*(.3*x1 + .4*x2 - x1) + (1-aplus)*(x1-x2+2*x1)
f2m = amin*(.3*x1 + .4*x2 - x1) + (1-amin)*(x1-x2+2*x1)
y = rbinom(N,1,logit(f2))
yplus = rbinom(N,1,logit(f2p))
ymin = rbinom(N,1,logit(f2m))

dat = data.frame(y=y,a=a,z=z,x1=x1,x2=x2,x3=x3)

true = mean(yplus - ymin)/mean(aplus - amin)
#phi = CACE(y=y,a=a,z=z,cov=cov,delta = delta,ranger = T)

trueSingle = mean(yplus - y)/mean(aplus - a)
phiSingle = CACE(y=y,a=a,z=z,cov=cov,delta = delta,ranger = T,type = 'single')

# compare to tsls -- weird to do since the parameter is not the usual one
first.stage = glm(a~z+x1+x2+x3,family = binomial(link=probit)) #E(A|Z,X)
inst = first.stage$fitted.values #Ahat
second.stage = glm(y~inst+x1+x2+x3,family = binomial(link=probit)) #E(Y|Ahat,X)
tsls.coef = second.stage$coefficients['inst']


# testing ranger
# r = ranger::ranger(y~z+x1+x2+x3,data = ds1, write.forest = T)
# newdat = ds2[,3:6]
# pred <- predict(r,data=newdat)$predictions

# testing propensity scores
data    <- as.data.frame(cbind(y,a,z,x1,x2,x3))
mysplit(data)
test = propscore_est(y=ds1[,3],x=ds1[,c(4:dim(ds1)[2])])
pred = predict(test, ds2[,c(4:dim(ds2)[2])])
probs = get_probs(ds2$z,pred$z,pred$CDE)

for(i in which(probs==0)){
  plot(pred$z,pred$CDE[i,],pch =19,main = i)
  abline(v=ds2$z[i])
}


## run a few times and output a curve

f <- function(x){
  tryCatch({CACE(y=y,a=a,z=z,cov=cov,delta = x)},
  error = function(e){print(paste("error is: ",e)); print(paste("val: ",x)); return(NA)})
}

vals = seq(min(z)*2,max(z)*.5, length = 5)
out = sapply(vals,function(x) f(x))

png(file = "sim_effect_curve5.png")
plot(vals,out,xlab = "delta", ylab = "phi")
dev.off()

vals = seq(min(z)*2,max(z)*.5, length = 10)
out = sapply(vals,function(x) f(x)) #started 1:17pm

png(file = "sim_effect_curve10.png")
plot(vals,unlist(as.matrix(out)[1,]),xlab = "delta", ylab = "phi")
dev.off()

###### simpler data set ########
N = 10000
x = 3 + rnorm(N)
z = exp(x/6) + rnorm(N)
delta = mean(z)/2
pi = expit(.8*z); piplus = expit(.8*(z+delta))
a = rbinom(N,1,pi); aplus = rbinom(N,1,piplus)
length(which(a>aplus))
y = 2*a + .4*x + rnorm(N); yplus = 2*(aplus) + .4*x + rnorm(N)
keep = a<=aplus

data <- data.frame(y=y[keep],a=a[keep],yplus=yplus[keep],aplus=aplus[keep],z=z[keep],x=x[keep])

true = mean(data$yplus-data$y)/mean(data$aplus - data$a)
phiSingle = CACE(y=data$y,a=data$a,z=data$z,cov=data$x,delta = delta,ranger = T,type = 'single')


mysplit(data)
test = propscore_est(y=ds1[,3],x=as.matrix(ds1$x))
pred = predict(test, ds2$x)
plot(pred$z,ds2$z)
probs = get_probs(ds2$z,pred$z,pred$CDE)

for(i in which(probs==0)){
  plot(pred$z,pred$CDE[i,],pch =19,main = i)
  abline(v=ds2$z[i])
}


dat = data.frame(y=y,a=a,z=z,x=x)

phi = CACE(y=y,a=a,z=z,cov=x,delta = delta,ranger = T)
true = mean(yplus-ymin)

# compare to tsls
first.stage = lm(a~z+x)
inst = first.stage$fitted.values
second.stage = lm(y~inst+x)
tsls.coef = second.stage$coefficients['inst']

vals = seq(min(z)*2,max(z)*.5, length = 10)
f <- function(x){
  tryCatch({CACE(y=y,a=a,z=z,cov=cov,delta = x)},
           error = function(e){print(paste("error is: ",e)); print(paste("val: ",x)); return(NA)})
}
out2 = sapply(vals,function(x) f(x))

png(file = "sim_effect_curve2.png")
plot(vals,out2$phi,xlab = "delta", ylab = "phi")
dev.off()


#### using kang and schafer data ####
dat <- ks_data(1000)
delta = .05
phi = CACE(y=dat$y,a=dat$r,z=dat$z1,cov=dat[,c(1:4)],delta = delta,ranger = T)

#### plug in ####
N = 10000
x = 3 + rnorm(N)
z = exp(x/6) + rnorm(N)
delta = mean(z)/2
pi = expit(.8*z); piplus = expit(.8*(z+delta))
plot(pi~z, pch = 19)
points(z, piplus, col = 'red', pch = 19)

a = rbinom(N,1,pi); aplus = rbinom(N,1,piplus)
length(which(a>aplus))
y = 2*a + .4*x + rnorm(N); yplus = 2*(aplus) + .4*x + rnorm(N)
keep = a<=aplus

dat <- data.frame(y=y[keep],a=a[keep],yplus=yplus[keep],aplus=aplus[keep],z=z[keep],x=x[keep])

true = mean(yplus - y)/mean(aplus - a)
phi_plug = CACE(y=dat$y,a=dat$a,z=dat$z,cov=dat$x,delta = delta,ranger = T, type = 'plugin')
phi_sing = CACE(y=dat$y,a=dat$a,z=dat$z,cov=dat$x,delta = delta,ranger = T, type = 'single')
phi_dub = CACE(y=dat$y,a=dat$a,z=dat$z,cov=dat$x,delta = delta,ranger = T, type = 'double')
