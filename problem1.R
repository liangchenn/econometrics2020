library(extRemes)
library(dplyr)
# Problem 1-------------------------------------------------------------------------


# Setups ------------------------------------------------------------------

n <- 400
beta1 <- 1.0
beta2 <- -0.5

x1 <- rnorm(n, 0, 1)
x2 <- rchisq(n, df = 1)

u1 <- revd(n)
u2 <- revd(n)



# (1) data Generating ---------------------------------------------------------------------

data <- data.frame(x1, x2, u1, u2) %>%
  mutate(y = ifelse(x1*beta1 + u1 > x2*beta2 + u2, 1, 0))

y <- data$y

# (2) log-likelihood function ---------------------------------------------

# first, build a distribution function called `cdf`
cdf <- function(x1, x2, beta1, beta2){
  
  u <- x1*beta1 - x2*beta2
  res <- exp(u) / (1+exp(u))
  
  return(res)
  
}


# define a function called `log_likelihood` using {beta1, beta2} as arguments.
# since R supports operations be taken elementwise, function could be written as follow.
log_likelihood <- function(beta1, beta2){
  
  res <- (y*log(cdf(x1, x2, beta1, beta2)) + (1-y)*log(1 - cdf(x1, x2, beta1, beta2)))
  res_sum <- sum(res)

  return(res_sum)
}

# test
log_likelihood(1, -0.5)



# checking result with dplyr ----------------------------------------------

res <- data %>%
  mutate(prob = cdf(x1, x2, -5, -5)) %>%
  mutate(llk = y*log(prob) + (1-y)*log((1-prob)))
sum(res$llk)
# -226.8987

# grid search -------------------------------------------------------------


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


get_mle_betas(-1, 1, 0.01)
get_mle_betas()
