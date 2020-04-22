# setwd('./hw2/')
library(rstan)
library(data.table)
library(magrittr)
library(maxLik)

# Problem 1 ---------------------------------------------------------------

# `preprocessing01_get_feature_name.R` get the features names and features to use.
source('./preprocessing01_get_features_name.R')

# extract x's columns
cols <- match(feature.to.use, feature.names)

# read all 76 groups/ networks data into 2 lists
files_nx <- list.files('./problem_set_2_sample/network/', full.names = T)
files_cov <- list.files('./problem_set_2_sample/group/', full.names = T)

nx <- lapply(files_nx, fread) 
group <- lapply(files_cov, fread)

# extract gpa y and covariates x for all 76 data, and convert into matrix
nx <- lapply(nx, as.matrix)
ys <- lapply(group, function(x) x[, .SD, .SDcols=which(feature.names=='gpa')] %>% as.matrix)
xs <- lapply(group, function(x) x[, ..cols]  %>% as.matrix)


# log-likelihood function
# theta = (lambda, sigma, alpha_g, beta1', beta2')', which is a 29x1 matrix

loglik_fn <- function(theta){
  
  lambda <- theta[1]
  sigma <- theta[2]
  alpha_g  <- theta[3]
  beta_  <- theta[-3:-1]
  beta1  <- beta_[1:(length(beta_)/2)]
  beta2  <- beta_[((length(beta_)/2) + 1):length(beta_)]
  
  llk <- 0
  for (i in 1:76) {
    
    W <- nx[[i]]
    ng <- nrow(W)
    y <- ys[[i]]
    x <- xs[[i]]
    
    ep <- (diag(ng) - lambda * W) %*% y - x %*% beta1 - W %*% x %*% beta2 - alpha_g
    res <- -(ng/2)*log(2*pi) - (ng)*log(sigma) + log(abs(det(diag(ng) - lambda*W))) - (1/2/sigma**2) * (t(ep) %*% ep)
    llk <- llk + res
  }
  
  return(llk)
}

# optimization using `maxLik`-package-provided BFGS method
# set the initial value to (0, 1, 0,...,0) (29x1)
theta_0 <- c(0, 1, rep(0, 27))

# estimate
t <- proc.time()
estimate <- maxLik(loglik_fn, start = theta_0, method = 'BFGS')
proc.time() - t

# coefs and std err
summary(estimate)
# output:
#          Estimate Std. error   t value    Pr(> t)    
# lambda   0.057609   0.008394   6.863     6.72e-12 ***
# sigma    0.697826   0.010994   63.476    < 2e-16 ***
# alpha_g  3.398482   0.151347   22.455    < 2e-16 ***











