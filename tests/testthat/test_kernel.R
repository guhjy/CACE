### test the kernel ###
gK <- function(x){(1/sqrt(2*pi))*exp(-(x^2)/2)}

test_that('kernel can work in simple case',{
  x = rnorm(1000)
  N = length(x)
  h = sd(x)
  test.seq = seq(-2,2,by=.2)
  test.pi = sapply(test.seq, function(y) (1/length(x))*sum(gK((1/h)*sqrt((y - x)^2))) )

  expect_lt(-cor(test.pi, dnorm(test.seq)),-.8)

  plot(test.pi~dnorm(test.seq),
       xlim = c(min(test.pi, dnorm(test.seq)),max(test.pi, dnorm(test.seq))),
       ylim = c(min(test.pi, dnorm(test.seq)),max(test.pi, dnorm(test.seq))))
  abline(0,1)
})

# dunno about this one..
test_that('kernel can work in conditional case',{
  x = rnorm(1000)
  N = length(x)
  h = sd(x)
  test.seq = x+1
  test.pi = sapply(test.seq, function(y) (1/length(x))*sum(gK((1/h)*sqrt((y - x)^2))) )

  expect_lt(-cor(test.pi, dnorm((test.seq-1))),-.8)

  plot(test.pi~dnorm((test.seq),mean = 1),
       xlim = c(min(test.pi, dnorm(test.seq)),max(test.pi, dnorm(test.seq))),
       ylim = c(min(test.pi, dnorm(test.seq)),max(test.pi, dnorm(test.seq))))
  abline(0,1)
})





