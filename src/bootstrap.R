# bootstraping methods implement for standard errors
library(mvtnorm) # for generating random number from multivariate normal


# residuals bootstraping --------------------------------------------------

# data generating
n  <- 1000
X  <- rmvnorm(n, rep(0, 4), sigma = diag(1, 4, 4))
b0 <- matrix(c(.5, 1.6, .7, 2.8), nrow = 4)

y <- X %*% b0 + rnorm(n)

# analytic form 
model <- lm(y ~ X)

# bootstraping

# step 0. get the beta hat from the data

b_hat <- solve(t(X) %*% X) %*% (t(X) %*% y)

# step 1. create the estimated residuals u_hat
y_hat <- X %*% b_hat
u_hat <- y - y_hat

# step 2. start bootstraping

b_times <- 1000
bs <- vector('list', b_times)
for(i in 1:b_times){
  
  # draw from u_hat with replacement
  b_idx <- sample(1:n, size = n, replace = T)
  u_bootstrap <- u_hat[b_idx]
  # rebuild the y_hat_b = Xb + u_b
  y_bootstrap <- X %*% b_hat + u_bootstrap
  
  # reestimate the b_hat_b
  b_hat_boostrap <- solve(t(X) %*% X) %*% (t(X) %*% y_bootstrap)
  
  # save the results
  bs[[i]] <- b_hat_boostrap
  
}

library(purrr)

# results
x1_b <- map_dbl(bs, 1)
x2_b <- map_dbl(bs, 2)
x3_b <- map_dbl(bs, 3)
x4_b <- map_dbl(bs, 4)

show_boostrap_results <- function(){
  
  message("\tEstimation\t\tStd. Error\t\tt-value")
  message('x1\t', mean(x1_b), '\t', sd(x1_b), '\t', (mean(x1_b) - 0)/ sd(x1_b))
  message('x2\t', mean(x2_b), '\t', sd(x2_b), '\t', (mean(x2_b) - 0)/ sd(x2_b))
  message('x3\t', mean(x3_b), '\t', sd(x3_b), '\t', (mean(x3_b) - 0)/ sd(x3_b))
  message('x4\t', mean(x4_b), '\t', sd(x4_b), '\t', (mean(x4_b) - 0)/ sd(x4_b))
}

# compare to the analytic results
show_boostrap_results()
summary(model)



# plot
plot(density(x4_b))
abline(v=quantile(x4_b)[2:4], lty='dashed', col =  'red')


# nonparametric bootstraping ----------------------------------------------

# data generating
n  <- 1000
X  <- rmvnorm(n, rep(0, 4), sigma = diag(1, 4, 4))
b0 <- matrix(c(.5, 1.6, .7, 2.8), nrow = 4)

y <- X %*% b0 + rnorm(n)

# bootstrap
bs2 <- vector('list', b_times)
b_times <- 1000
for (i in 1:b_times) {
  idx <- sample(1:n, 1000, replace = T)
  Xb <- X[idx, ]
  yb <- y[idx, ]
  bs2[[i]] <- lm(yb ~ Xb)$coefficients[2:5]

}

# results
x1_b2 <- map_dbl(bs2, 1)
x2_b2 <- map_dbl(bs2, 2)
x3_b2 <- map_dbl(bs2, 3)
x4_b2 <- map_dbl(bs2, 4)

show_boostrap_results2 <- function(){
  
  message("\tEstimation\t\tStd. Error\t\tt-value")
  message('x1\t', mean(x1_b2), '\t', sd(x1_b2), '\t', (mean(x1_b2) - 0)/ sd(x1_b2))
  message('x2\t', mean(x2_b2), '\t', sd(x2_b2), '\t', (mean(x2_b2) - 0)/ sd(x2_b2))
  message('x3\t', mean(x3_b2), '\t', sd(x3_b2), '\t', (mean(x3_b2) - 0)/ sd(x3_b2))
  message('x4\t', mean(x4_b2), '\t', sd(x4_b2), '\t', (mean(x4_b2) - 0)/ sd(x4_b2))
}

show_boostrap_results2()


model <- lm(y ~ X)
summary(model)

sd1 <- sqrt(mean((x1_b2 - mean(x1_b2))**2))


