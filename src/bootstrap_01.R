# Non parametric bootstrap ------------------------------------------------

x <- rnorm(50, 2, 1)
y <- x*-2.3 + rnorm(50)

R <- 1000
betas <- vector("numeric", R)

B_TIMES <- 1000

for (r in 1:R) {
  
  bootstraps <- sample(1:50, B_TIMES, replace = T)
  
  x_b <- x[bootstraps]
  y_b <- y[bootstraps]
  
  
  betas[r] <- lm(y_b ~ x_b)$coef[2]
  
  
}

plot(density(betas))

quantile(betas, c(.05, .95))
sd(betas)




# residual bootstrap ------------------------------------------------------

beta_hat <- lm(y~x)$coef[2]
errors <- y - x*beta_hat

betas <- vector('numeric', B_TIMES)
for (b in 1:B_TIMES) {
  
  u_b <- errors[sample(1:50, 50, replace = T)]
  
  y_b <- x*beta_hat + u_b
  
  betas[b] <- lm(y_b ~ x)$coef[2]
  
}

plot(density(betas))
quantile(betas, c(.05, .25, .5, .75, .95))

#      5%       25%       50%       75%       95% 
# -2.571124 -2.444133 -2.351911 -2.266474 -2.128981 




