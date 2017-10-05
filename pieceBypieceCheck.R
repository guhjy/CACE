##### goal: test how well each part of the shift estimator does ####
# we use the edited ed simulation
test.sim <- makeSim()


#### propensity scores ####
# get estimates of propensity score using flexcode
props <- propscore_est(y = test.sim$z, x = test.sim[,5:8])
pred <- predict(props, test.sim[,5:8])
pi <- get_probs(test.sim$z,pred$z,pred$CDE)

#get estimates of propensity score using glm
prop.glm <- glm(z~V5+V6+V7+V8, data = test.sim)
pi.glm <- predict(prop.glm, test.sim[,5:8])

# get true propensity score
mns <- 1.5*sign(test.sim[,5] + test.sim[,6] - test.sim[,7] - test.sim[,8])
true.pi <- dnorm(test.sim$z, mns, sd = 2)

#compare with flexcode
cor(pi,true.pi)
png(filename = 'flexCodeEstvsTrue.png', height = 480, width = 480)
plot(pi~true.pi, xlim = c(0,.2), ylim = c(0,.2), xlab = 'True', ylab = 'FlexCode estimate', main = "Propensity Score Estimation")
abline(a=0,b=1, col = 'red')
dev.off()

plot(pi[mns == -1.5] ~ true.pi[mns == -1.5])
plot(pi[mns == 1.5] ~ true.pi[mns == 1.5])

#compare with glm
cor(pi.glm, true.pi); plot(pi.glm~true.pi)

#### means ####
# true means
mnT <- test.sim$V5 - test.sim$V6 - test.sim$V7 + test.sim$V8 + test.sim$y0
mnA <- pnorm(test.sim$z, mnT, 1)
mnY <- mnT*mnA

require(ranger)
# set up data frame:
dat1 <- test.sim[,-c(2,4,9)]
dat2 <- test.sim[,-c(1,4,9)]

# run ranger & get predictions
rang.res1 <- ranger::ranger(y~.,data = dat1,write.forest = T)
rang.res2 <- ranger::ranger(a~.,data = dat2,write.forest = T)
xnew <- test.sim[,c(3,5:8)]
muY <- predict(rang.res1,data = xnew)$pred
muA <- predict(rang.res2,data = xnew)$pred

# try superlearner
require(SuperLearner)
sl.lib <- c("SL.glm", "SL.randomForest","SL.polymars","SL.mean")
muA.SL <- SuperLearner(Y=test.sim$a, X=test.sim[,c(3,5:8)], SL.library=sl.lib, family=binomial())
sl.pred.a <- predict(muA.SL, test.sim[,c(3,5:8)])$pred

#run glm
glm.a <- glm(a~.,data = dat2, family = binomial('probit'))
glm.y <- glm(y~.,data = dat1)
glm.pred.a <- predict(glm.a, dat2)
glm.pred.y <- predict(glm.y, dat1)

#compare
cor(muA,mnA)
png(filename = 'RangerEstMeanAvsTrue.png', width = 480, height = 480)
plot(muA~mnA, xlab = 'True', ylab = 'Ranger Estimate', xlim = c(0,1), ylim = c(0,1))
abline(0,1,col='red')
dev.off()

cor(glm.pred.a,mnA)
png(filename = 'GLMEstMeanAvsTrue.png', width = 480, height = 480)
plot(glm.pred.a ~ mnA, xlab = 'True', ylab = 'GLM Estimate')
abline(0,1,col='red')
dev.off()

cor(muY,mnY)

png(filename = 'RangerEstMeanYvsTrue.png', width = 480, height = 480)
plot(muY~mnY, xlab = 'True', ylab = 'Ranger Estimate', xlim = c(ll,ul), ylim = c(ll,ul),
     main = "Estimates of mu(Y) vs. True mu(Y)")
abline(0,1,col='red')
dev.off()

plot(glm.pred.y~mnY)
cor(glm.pred.y,mnY)

ll = min(sl.pred.a,mnA); ul = max(sl.pred.a,mnA)
png(filename = 'SLEstMeanAvsTrue.png', width = 480, height = 480)
plot(sl.pred.a~mnA, xlim = c(ll,ul), ylim = c(ll,ul), xlab = "True", ylab = "Estimated using SuperLearner")
abline(0,1,col = 'red')
dev.off()

#### estimate single shift ###
delta = 1
mnA.plus <- pnorm(test.sim$z + delta, mnT, 1)
true.eff <- mean(mnT*(mnA.plus - mnA)) / mean(mnA.plus - mnA)

props.min <- propscore_est(y = test.sim$z - delta, x = test.sim[,5:8])
pred.min <- predict(props.min, test.sim[,5:8])
pi.min <- get_probs(test.sim$z - delta,pred.min$z,pred.min$CDE)

covs = cbind(test.sim$z+delta,test.sim[,5:8])
muA.plus <- SuperLearner(Y=test.sim$a, X=covs, SL.library=sl.lib, family=binomial())
predA.plus <- predict(muA.plus, covs)$pred

rang.res.plus <- ranger::ranger(y~.,data = test.sim[,c(1,3,5:8)],write.forest = T)
xnew <- test.sim[,c(3,5:8)]
muY.plus <- predict(rang.res.plus,data = xnew)$pred

phi_y = (test.sim$y - muY)*(pi.min/pi) - (test.sim$y - muY.plus)
phi_a = (test.sim$a - muY)*(pi.min/pi) - (test.sim$a - muY.plus)

print(paste(length(which(pi==0)),"zero probability values"));keep = which(pi!=0)
psihat = mean(phi_y[keep])/mean(phi_a[keep])
