test_that('glm y est doing ok',{
  z <- rnorm(1000)
  x <- rnorm(1000)
  y <- 3*x - .5*z + rnorm(1000)
  test.glm <- y.mean.est(y,cbind(z,x),'glm')
  predict <- predict(test.glm, data.frame(z=z, x=x))
  expect_equal(class(test.glm), c('glm', 'lm'))
  expect_lt(-cor(predict,y),-.5)
})

test_that('ranger y est doing ok',{
  x <- rnorm(1000)
  z <- rnorm(1000)
  y <- 3*x - .5*z + rnorm(1000)
  test.r <- y.mean.est(y,cbind(z,x),'Ranger')
  predict <- predict(test.r, cbind(z=z, x=x))$pred
  expect_equal(class(test.r), 'ranger')
  expect_lt(-cor(predict,y),-.5)
})

test_that('sl y est doing ok',{
  x <- rnorm(1000)
  z <- rnorm(1000)
  y <- 3*x - .5*z + rnorm(1000)
  test.sl <- y.mean.est(y,as.data.frame(cbind(z,x)),'superlearner')
  predict <- predict(test.sl, as.data.frame(cbind(z=z, x=x)))$pred
  expect_equal(class(test.sl), 'SuperLearner')
  expect_lt(-cor(predict,y),-.5)
})

test_that('algorithm error message',{
  x <- rnorm(1000)
  z <- rnorm(1000)
  y <- 3*x - .5*z + rnorm(1000)
  expect_error(y.mean.est(y,as.data.frame(cbind(z,x)),'blah'))
})


test_that('glm a est doing ok',{
  x <- rnorm(1000); z <- rnorm(1000)
  a <- as.numeric(3*x - .5*z + rnorm(1000) >=1)
  test.glm <- a.mean.est(a,cbind(z,x),'glm')
  predict <- predict(test.glm, data.frame(z=z, x=x), type = 'response')
  true.prob <- 1 - pnorm(1-3*x+.5*z)
  expect_equal(class(test.glm), c('glm', 'lm'))
  expect_lt(-cor(predict,a),-.5)
})

test_that('ranger a est doing ok',{
  x <- rnorm(1000)
  z <- rnorm(1000)
  a <- as.numeric(3*x - .5*z + rnorm(1000) >=1)
  test.r <- a.mean.est(a,cbind(z,x),'Ranger')
  predict <- predict(test.r, cbind(z=z, x=x))$pred
  expect_equal(class(test.r), 'ranger')
  expect_lt(-cor(predict,a),-.5)
})

test_that('sl a est doing ok',{
  x <- rnorm(1000)
  z <- rnorm(1000)
  a <- as.numeric(3*x - .5*z + rnorm(1000) >=1)
  test.sl <- a.mean.est(a,as.data.frame(cbind(z,x)),'superlearner')
  predict <- predict(test.sl, as.data.frame(cbind(z=z, x=x)))$pred
  expect_equal(class(test.sl), 'SuperLearner')
  expect_lt(-cor(predict,a),-.5)
})


test_that('glm z est doing ok',{
  x <- rnorm(1000)
  z <- exp(x) - x^2 + rnorm(1000)
  test.glm <- z.condldens.est(z,x,'glm')
  expect_equal(class(test.glm), c('glm', 'lm'))
  predict <- predict(test.glm, data.frame(z=z, x=x))
  expect_lt(-cor(predict,z),-.5)
})

test_that('flexcode z est doing ok',{
  x1 <- rnorm(1000,1,1)
  x2 <- rnorm(1000,1,1)
  z <- 2*x1 - x2 + rnorm(1000)
  true.prob <- dnorm(z, mean = 1, sd = sqrt(6))
  test.fc <- z.condldens.est(z,cbind(x1,x2),'flexcode')
  expect_equal(class(test.fc), 'FlexCoDE')
  predict <- predict(test.fc, cbind(x1,x2))
  probs <- get_probs(z, predict$z, predict$CDE)
  expect_lt(-cor(probs,true.prob),-.5)
})

test_that('single shift wrapper doing ok',{
  x1 <- rnorm(1000,1,1)
  x2 <- rnorm(1000,1,1)
  z <- rnorm(1000,1,1)
  a <- as.numeric(exp(x1) - x2^2 + z >= 0)
  y = a*(x1 + x2) + (2*x1 - x2) + rnorm(1000)
  test.out <- single.shift(y = y,a = a,z = z,x = cbind(x1,x2), delta = 1)
  test.out <- single.shift(y = y,a = a,z = z,x = cbind(x1,x2), delta = 1,
                           algo = list(y.est = 'glm', a.est = 'glm', z.est = 'glm'))
  test.out <- single.shift(y = y,a = a,z = z,x = cbind(x1,x2), delta = 1,
                           algo = list(y.est = 'glm', a.est = 'glm', z.est = 'flexcode'))
})

#consistently biased downwards
test_that('getting close on easy sim',{
  N = 5000
  y0 = rnorm(N,0,1); meanx = c(0,0,0,0)
  alpha = matrix(c(1,1,-1,-1),nrow = 4); beta = matrix(c(1,-1,-1,1))
  x = matrix(unlist(lapply(meanx, function(x) rnorm(N,x,1))), nrow = N, byrow =T)
  z = rnorm(N, 1.5*sign(t(alpha)%*%t(x)), sd = 2)
  t = 2; psi = 1
  a = as.numeric(z >= t)
  y = y0 + a*psi*t
  true.eff = t*psi
  test.out <- single.shift(y = y,a = a,z = z,x = x, delta = 2)
  print(test.out$psi)
  expect_true( (test.out$psi >= true.eff*.8) & (test.out$psi <= true.eff*1.2) )
})

#consistently biased downwards, slightly worse
test_that('plug in version',{
  N = 5000
  y0 = rnorm(N,0,1); meanx = c(0,0,0,0)
  alpha = matrix(c(1,1,-1,-1),nrow = 4); beta = matrix(c(1,-1,-1,1))
  x = matrix(unlist(lapply(meanx, function(x) rnorm(N,x,1))), nrow = N, byrow =T)
  z = rnorm(N, 1.5*sign(t(alpha)%*%t(x)), sd = 2)
  t = 2; psi = 1
  a = as.numeric(z >= t)
  y = y0 + a*psi*t
  true.eff = t*psi
  test.out <- single.shift.pi(y = y,a = a,z=z,x = x, delta = 2)
  print(test.out$psi)
  expect_true( (test.out$psi >= true.eff*.8) & (test.out$psi <= true.eff*1.2) )
})
