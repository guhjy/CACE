test_that('glm z est doing ok',{
  x <- rnorm(1000)
  z <- exp(x) - x^2 + rnorm(1000)
  test.glm <- z.condldens.est(z,x,'glm')
  expect_equal(class(test.glm), c('glm', 'lm'))
  # predict <- predict(test.glm, data.frame(z=z, x=x))
  # expect_lt(-cor(predict,y),-.5)
})
