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
  N <- 1000   # total num of MCMC
  
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
    
    if(i %% 100 == 0) t <- Sys.time()
    
    # =============================       posterior of beta     =============================
    
    ZVX <- 0
    ZVY <- 0
    V   <- sige_t[, i-1]
    
    for (g in 1:G) {
      
      ng <- nrow(W[[g]])
      YY <- (diag(ng) - lambda_t[, (i-1)]*W[[g]]) %*% Y[[g]] - alpha_t[g, (i-1)]
      ZZ <- ZZZ[[g]]
      
      ZVX <- ZVX + (t(ZZ) %*% ZZ) * V^(-1)
      ZVY <- ZVY + (t(ZZ) %*% YY) * V^(-1)
      
    }
    
    B_ <- inv(solve(B.0) + ZVX)
    beta_ <- B_ %*% (solve(B.0) %*% beta.0 + ZVY)
    
    # update beta_t
    beta_t[, i] <- rmvnorm(1, mean = beta_, sigma = B_)
    
    
    # =============================       posterior of sigmae     =============================
    
    ee <- 0
    total_n <- 0
    
    for(g in 1:76){
      
      ng <- nrow(W[[g]])
      ZZ <-  ZZZ[[g]]
      e <- (diag(ng) - lambda_t[, (i-1)] * W[[g]]) %*% Y[[g]] - ZZ %*% beta_t[, i] - alpha_t[g, (i-1)]*ones(ng)
      
      ee <- ee + t(e) %*% e
      total_n <- total_n + nrow(W[[g]])
    }
    
    # update sigma_e2_t
    sige_t[, i] <- (eta.0 + ee) / rchisq(1, df = rho.0 + total_n)
    
    
    
    
    # =============================       posterior of alpha     =============================
    
    for (g in 1:76){
      
      ng <- nrow(W[[g]])
      
      Rg <- solve(solve(A.0) + sige_t[, i]^(-1) * t(ones(ng)) %*% ones(ng))
      
      ZZ <- ZZZ[[g]]
      
      E  <- (diag(ng) - lambda_t[, (i-1)] * W[[g]]) %*% Y[[g]] - ZZ %*% beta_t[, i] - alpha_t[g, (i-1)]
      
      Alpha_g <- Rg %*% (inv(A.0) %*% Alpha.0 + sige_t[, i]^(-1) * t(ones(ng)) %*% E)
      
      
      alpha_t[g, i] <- rmvnorm(1, mean = Alpha_g, sigma = Rg)
      
    }
    
    
    # =============================       posterior of lambda     =============================
    
    # propose a lambda_
    
    if(i < 3){
      
      lambda_ <- rnorm(1, mean = lambda_t[, (i-1)], sd = 0.1)
      
    }else{
      
      sig_propose <- 0.95*( cov(matrix(lambda_t[1:(i-1)]))*2.83^2 ) + 0.05*0.1^2
      
      lambda_ <- rmvnorm(1, mean = lambda_t[, (i-1)], sigma =  sig_propose )[1]
      
    }
    
    
    
    ppl <- 1
    
    for(g in 1:76){
      
      ng <- nrow(W[[g]])
      ZZ <- cbind(X[[g]], W[[g]] %*% X[[g]])
      
      S1 <- ( diag(ng) - lambda_t[, (i-1)] * W[[g]] )
      S2 <- ( diag(ng) - lambda_ * W[[g]] )
      
      e1 <- S1 %*% Y[[g]] - ZZ %*% beta_t[, i] - alpha_t[g, i]
      e2 <- S2 %*% Y[[g]] - ZZ %*% beta_t[, i] - alpha_t[g, i]
      
      likelihood_1 <- det(S1) * exp(-.5 * (t(e1) %*% e1) * sige_t[, i]^(-1))
      likelihood_2 <- det(S2) * exp(-.5 * (t(e2) %*% e2) * sige_t[, i]^(-1))
      
      ppl = ppl * (likelihood_1/likelihood_2)
      # print(ppl)
    }
    
    ppl <- min(ppl, 1)
    
    if(runif(1) <= ppl){
      # accept propose
      lambda_t[, i] <- lambda_
      
    }else{
      # reject
      lambda_t[, i] <- lambda_t[, (i-1)]
    }
    
    if(i %% 100 ==0){
      elapsed.time <- Sys.time() - t
      message("[",i,"] \t using : \t", elapsed.time, "sec")
    }
  }
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
}









