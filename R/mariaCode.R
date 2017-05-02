mariacode <- function(){
  # true parameters
  samplesize = 1000 # test
  index = 1:samplesize
  beta = 0.5
  Y1 = rbinom(n = samplesize, size = 1, prob = beta)

  X1 = rnorm(n = samplesize, mean = 0, sd = 1) # correct model
  X2 = rnorm(n = samplesize, mean = 0, sd = 1)
  X3 = rnorm(n = samplesize, mean = 0, sd = 1)
  X4 = rnorm(n = samplesize, mean = 0, sd = 1)

  X1star = exp(X1/2) # misspecified model
  X2star = X2/(1 + exp(X1)) + 10
  X3star = (X1*X3/25 + 0.6)^3
  X4star = (X2 + X4 + 20)^2

  pi = expit(-X1+0.5*X2-0.25*X3-0.1*X4)

  A = rbinom(n = samplesize, size = 1, prob = expit(-X1+0.5*X2-0.25*X3-0.1*X4))

  # Defining my parameter
  beta = 0.5
  mu0 = beta/(1+exp(-X1+0.5*X2-0.25*X3-0.1*X4))
  mu1 = rep(beta, samplesize)
  #gamma = 1 - mu0/mu1
  gamma = expit(-X1+0.5*X2-0.25*X3-0.1*X4) # Don't we need to make it so that gamma is equal to 1-mu0/mu1??

  # generate data frame
  df = as.data.frame(cbind(index, X1, X2, X3, X4, X1star, X2star, X3star, X4star, pi, gamma, A, Y1))

  # generate Y0 conditional on combinations of values of A and Y1
  dfy11 = df[which(df$Y1==1),]
  dfy11$Y0 = rbinom(n = nrow(dfy11) , size = 1, prob = (1-gamma)) # or is it location = expit(t(PSI)*X), scale = 0?

  dfy10 = df[which(df$Y1==0),]
  dfy10$Y0 = 0

  # add Y0 to dataframe
  df_wy0 = as.data.frame(rbind(dfy11, dfy10))

  # apply consistency to get Y
  df_wy0$Y = ifelse(df_wy0$A==1, df_wy0$Y1, df_wy0$Y0)

  # ordering data so it's as it was at the beginning
  dff = df_wy0[order(df_wy0$index),]

  #And here's the same superlearner code:
  sl.lib2 <- c("SL.glm", "SL.randomForest", "SL.gam", "SL.polymars", "SL.mean")


  # Fitting model
  data0 = dff[which(dff$A==0),] # need to put data into superlearner as df
  data1 = dff[which(dff$A==1),]

  mis_N_PI_mu0SL = SuperLearner(Y=data0$Y, X=as.data.frame(cbind(data0$X1star,data0$X2star,data0$X3star,data0$X4star)),
                                SL.library = sl.lib2, family=binomial())
  mis_N_PI_mu1SL = SuperLearner(Y=data1$Y, X=as.data.frame(cbind(data1$X1star,data1$X2star,data1$X3star,data1$X4star)),
                                SL.library = sl.lib2, family=binomial())

  # Getting fitted values (after inverse link function?)
  mis_N_PI_mu0_hatSL = predict(mis_N_PI_mu0SL, newdata=as.data.frame(cbind(dff$X1star,dff$X2star,dff$X3star,dff$X4star)))$pred
  mis_N_PI_mu1_hatSL = predict(mis_N_PI_mu1SL, newdata=as.data.frame(cbind(dff$X1star,dff$X2star,dff$X3star,dff$X4star)))$pred

  # Gamma hat
  mis_N_PI_gammahatSL = (mis_N_PI_mu1_hatSL - mis_N_PI_mu0_hatSL)/mis_N_PI_mu1_hatSL

  # RMSE
  RMSE_mis_N_PISL = sqrt(  mean( (mis_N_PI_gammahatSL - gamma)^2 )  )
  RMSE_mis_N_PI = RMSE_mis_N_PISL

}
