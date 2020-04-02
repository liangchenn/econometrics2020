# ridge regression in numerical method
library(dplyr)
beta <- c(1, -.4, 7)
# DGP
x1 <- rnorm(1000)
x2 <- rnorm(1000)
x3 <- rchisq(1000, 1)
e  <- rnorm(1000)

# response
y <- beta %*% matrix(c(x1, x2, x3), nrow = 3, ncol = 1000) + e
y <- y[1, ]



obj_f <- function(beta){
  
  hat <- beta %*% matrix(c(x1, x2, x3), nrow = 3, ncol = 1000)
  
  res <- sum(-(y - hat)**2)
  
  return(res)
  
}


obj_f2 <- function(beta){
  
  
  sum(-(y - beta %*% matrix(c(x1, x2, x3), nrow = 3, ncol = 1000))**2) - 0.1*beta %*% beta
  
  
}

obj_f(c(1, -.4, 7))

library(maxLik)


res <- maxLik(obj_f2, start = c(0, 0, 0))





