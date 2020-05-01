library(mvtnorm)
library(invgamma)
library(magrittr)

ones <- function(n){
  
  res <- matrix(rep(1, n))
  return(res)
  
}
zeros <- function(n){
  res <- matrix(rep(0, n))
  return(res)
}
inv  <- function(X){
  
  return(solve(X))
}


# load data information W, Y, X both are list with length of 76


load('data.RData')


# MCMC for R times (Monte Carlo)

R <- 1000
monte <- 0

while(monte <= R){
  
  monte <- monte + 1
  message(" ===================== This is ", monte, "-th times monte carlo simulation. ======================")
  # ============================= prior belief initial values =============================
  
  k <- 13     # total num of independent variables
  G <- 76     # total num of groups
  N <- 100   # total num of MCMC
  
  # beta ~ N(beta.0, B.0)  (2k by 1)
  beta.0 <- zeros(2*k)
  B.0    <- diag(2*k) * 3
  
  # sigma_e2 ~ invgamma((rho.0 + N)/2, (eta.0 + e'e)/2) ~ (eta.0 + e'e) * (1/ chisq(rho.0 + N))  (1 by 1)
  rho.0 <- 2.2
  eta.0 <- 0.1
  
  # alpha_g ~ N(Alpha.0, A.0)  (1 by 1)
  Alpha.0 <- 1
  A.0     <- 1
  
  
  # =============================       estimation store     =============================
  # T = 1:N
  
  beta_t   <- matrix(0, nrow = 2*k, ncol = N)
  alpha_t  <- matrix(0, nrow = G, ncol = N)
  sige_t   <- matrix(1, nrow = 1, ncol = N)
  lambda_t <- matrix(0, nrow = 1, ncol = N)
  
  
  
  for (i in 2:N){
    
    # =============================       posterior of lambda     =============================
    
    # propose a lambda_
    
    if(i < 3){
      
      lambda_ <- rmvnorm(1, mean = lambda_t[, (i-1)], sigma = diag(1)*0.1^2)[1]
      
    }else{
      
      a <- rmvnorm(1, mean = lambda_t[, (i-1)], sigma = matrix(lambda_t[1:(i-1)]) %>% cov()*(2.83^2))[1]
      b <- rmvnorm(1, mean = lambda_t[, (i-1)], sigma = diag(1)*0.1^2 )[1]
      
      lambda_ <- 0.95 * a + 0.05 * b
      
    }
    
    
    
    ppl_ln <- 0
    
    # V = sige_t[, (i-1)]
    
    for(g in 1:76){
      ng <- nrow(W[[g]])
      
      ZZ <- ZZZ[[g]]
      
      S1 <- ( diag(ng) - lambda_t[, (i-1)] * W[[g]] )
      S2 <- ( diag(ng) - lambda_ * W[[g]] )
      
      e1 <- S1 %*% Y[[g]] - ZZ %*% beta_t[, i] - alpha_t[g, i]
      e2 <- S2 %*% Y[[g]] - ZZ %*% beta_t[, i] - alpha_t[g, i]
      
      likelihood_1 <- det(S1) * exp(-.5 * (t(e1) %*% e1) * sige_t[, (i-1)]^(-1))
      likelihood_2 <- det(S2) * exp(-.5 * (t(e2) %*% e2) * sige_t[, (i-1)]^(-1))
      
      ppl_ln = ppl_ln + log(1+likelihood_1)- log(1+likelihood_2)
      message(g, '\t', ppl_ln)
    }
    
    ppl_ <- min(exp(ppl_ln), 1)
    
    if(runif(1) <= ppl_ | is.na(ppl_)){
      # accept propose
      lambda_t[, i] <- lambda_
      
    }else{
      # reject
      lambda_t[, i] <- lambda_t[, (i-1)]
    }
    
    
    
    # =============================       posterior of beta     =============================
    
    ZVX <- 0
    ZVY <- 0
    V   <- sige_t[, i-1]
    
    for (g in 1:G) {
      
      ng <- nrow(W[[g]])
      YY <- (diag(ng) - lambda_t[, i]*W[[g]]) %*% Y[[g]] - alpha_t[g, (i-1)]
      ZZ <- ZZZ[[g]]
      
      ZVX <- ZVX + (t(ZZ) %*% ZZ) * V^(-1)
      ZVY <- ZVY + (t(ZZ) %*% YY) * V^(-1)
      
    }
    
    B_ <- inv(inv(B.0) + ZVX)
    beta_ <- B_ %*% (solve(B.0) %*% beta.0 + ZVY)
    
    # update beta_t
    beta_t[, i] <- rmvnorm(1, mean = beta_, sigma = B_)
    
    
    # =============================       posterior of sigmae     =============================
    
    ee <- 0
    total_n <- 0
    
    for(g in 1:76){
      
      ng <- nrow(W[[g]])
      ZZ <-  ZZZ[[g]]
      e <- (diag(ng) - lambda_t[, i] * W[[g]]) %*% Y[[g]] - ZZ %*% beta_t[, i] - alpha_t[g, (i-1)]*ones(ng)
      
      ee <- ee + t(e) %*% e
      total_n <- total_n + nrow(W[[g]])
    }
    
    # update sigma_e2_t
    sige_t[, i] <- (eta.0 + ee) / rchisq(1, df = rho.0 + total_n)
    
    
    
    
    # =============================       posterior of alpha     =============================
    
    for (g in 1:76){
      
      ng <- nrow(W[[g]])
      
      # Rg <- inv(solve(A.0) + sige_t[, i]^(-1) * t(ones(ng)) %*% ones(ng))
      dd <- ( Alpha.0^(-1) + (sige_t[, i]^(-1)*t(ones(ng)) %*% inv(diag(ng)) %*% ones(ng)) )^(-1)
      SS <- diag(ng) - lambda_t[, i]*W[[g]]
      YY <- SS %*% Y[[g]]
      XX <- ZZZ[[g]]
      
      # E  <- (diag(ng) - lambda_t[, (i-1)] * W[[g]]) %*% Y[[g]] - ZZ %*% beta_t[, i] - alpha_t[g, (i-1)]
      # 
      # Alpha_g <- Rg %*% (inv(A.0) %*% Alpha.0 + sige_t[, i]^(-1) * t(ones(ng)) %*% E)
      # 
      
      # alpha_t[g, i] <- rmvnorm(1, mean = Alpha_g, sigma = Rg)
      
      # update aplph_g
      alpha_t[g, i] <- sige_t[i]^(-1) * dd %*% t(ones(ng)) %*% inv(diag(ng)) %*% (YY - XX%*%beta_t[, i])+
        rnorm(1)*sqrt(dd)
      
      
    }
    
    
  }
  
  
  
}









