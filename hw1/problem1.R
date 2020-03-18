library(extRemes)
library(dplyr)
# Problem 1-------------------------------------------------------------------------

# Setups ------------------------------------------------------------------
set.seed(123)
n <- 400
beta1 <- 1.0
beta2 <- -0.5

x1 <- rnorm(n, 0, 1)
x2 <- rchisq(n, df = 1)

u1 <- revd(n)
u2 <- revd(n)

# (1) data generating ---------------------------------------------------------------------

y = ifelse(x1*beta1 + u1 > x2*beta2 + u2, 1, 0)

# (2) log-likelihood function ---------------------------------------------

# first, build a distribution function called `cdf`.
# using CDf of logit distribution, since subtraction of 2 type-I distribution is logit-distributed.
cdf <- function(x1, x2, beta1, beta2){
  
  u <- x1*beta1 - x2*beta2
  res <- exp(u) / (1+exp(u))
  
  return(res)
  
}


# define a function called `log_likelihood` using {beta1, beta2} as arguments.
log_likelihood <- function(beta1, beta2){
  
  res <- (y*log(cdf(x1, x2, beta1, beta2)) + (1-y)*log(1 - cdf(x1, x2, beta1, beta2)))
  res_sum <- sum(res)

  return(res_sum)
}

# test
log_likelihood(1, -0.5)




# (3) grid search -------------------------------------------------------------

get_mle_betas <- function(from = -5, to = 5, step = 0.1){
  # expand grid
  grids <- expand.grid(seq(from, to, step), seq(from, to, step))
  
  # calculate llk for each set of beta
  grids$llk <- mapply(log_likelihood, grids$Var1, grids$Var2)
  
  # get index of the maximum value of log likelihood ratio
  index <- which.max(grids$llk)  # (0.8, -0.5)
  betas <- grids[index, ]
  # returning betas
  return(betas)
}

get_mle_betas()
get_mle_betas(-1, 1, 0.01)



# (4) BHHH ----------------------------------------------------------------
library(maxLik)

# rewrite the log likelihood function form to meet maxLik requirements
cdf <- function(x1, x2, beta1, beta2){
  
  u <- x1*beta1 - x2*beta2
  res <- exp(u) / (1+exp(u))
  
  return(res)
  
}
log_likelihood <- function(beta){
  beta1 <- beta[1]
  beta2 <- beta[2]
  
  res <- (y*log(cdf(x1, x2, beta1, beta2)) + (1-y)*log(1 - cdf(x1, x2, beta1, beta2)))

  return(res)
}

BHHH <- maxLik(log_likelihood, start = c(0, 0), method = 'BHHH')
BHHH$estimate

# R = 100, N = 400
t <- proc.time()
betas <- list()
for (i in 1:100) {
  # setups
  n <- 400
  beta1 <- 1.0
  beta2 <- -0.5
  x1 <- rnorm(n, 0, 1)
  x2 <- rchisq(n, df = 1)
  u1 <- revd(n)
  u2 <- revd(n)
  y <- ifelse(x1*beta1 + u1 > x2*beta2 + u2, 1, 0)
  
  # BHHH estimation
  BHHH <- maxLik(log_likelihood, start = c(0, 0), method = 'BHHH')
  
  betas[[i]] <- BHHH$estimate
  message(i, '-th simulation...')
}
proc.time() - t
# the computation time: 1.087s

library(purrr)

show_results <- function(betas) {
  beta1s <- map_dbl(.x = betas, .f = 1)
  beta2s <- map_dbl(betas, 2)
  
  message("The mean of beta 1 is: ", mean(beta1s), " , and the sd is: ", sd(beta1s))
  message("The mean of beta 2 is: ", mean(beta2s), " , and the sd is: ", sd(beta2s))
}

show_results(betas = betas)

# The mean of beta 1 is: 0.99400702590802 , and the sd is: 0.128792152916141
# The mean of beta 2 is: -0.50031288370689 , and the sd is: 0.091142855808259


# (4-2) Nelder Mead Method ------------------------------------------------

betas_nm <- list()
for (i in 1:100) {
  # setups
  n <- 400
  beta1 <- 1.0
  beta2 <- -0.5
  x1 <- rnorm(n, 0, 1)
  x2 <- rchisq(n, df = 1)
  u1 <- revd(n)
  u2 <- revd(n)
  y <- ifelse(x1*beta1 + u1 > x2*beta2 + u2, 1, 0)
  
  # BHHH estimation
  NM <- maxLik(log_likelihood, start = c(0, 0), method = 'NM')
  
  betas_nm[[i]] <- NM$estimate
  # message(i, '-th simulation...')
}

show_results(betas_nm)

# The mean of beta 1 is: 1.00602680797551 , and the sd is: 0.13462494456679
# The mean of beta 2 is: -0.509983180077846 , and the sd is: 0.0952679195463455


