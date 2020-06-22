# library(rstanarm)
library(readstata13)
library(data.table)
library(MASS)

data <- read.dta13('./hw3/5368_data and programs_0/CilibertoTamerEconometrica.dta')
setDT(data)


# question 1 --------------------------------------------------------------

probit_data <- data[, totalN := rowSums(.SD), .SDcols = names(data) %like% 'airline[A-Z]+'
                    ][, lnN := log(totalN)][, lnN := ifelse(lnN <= 0, 0, lnN)][, N := as.ordered(totalN)
                      ][, .(market, marketdistance, marketsize, N, lnN)][, m2 := marketdistance**2]

o.probit <- polr("N ~ marketdistance + marketsize + lnN", data = probit_data, 
                 method = 'probit', start = c(o.probit$coefficients, 90, o.probit$zeta))









# question 2 --------------------------------------------------------------

# simulation
# draw ui0 and uik 
# params <- rep(.5, 6)
foo <- function(params){
  
  # hyperparameters
  beta  <- params[1:2]
  alpha <- params[3:4]
  delta <- params[5]
  rho   <- params[6]
  
  # use N_hat_sum to record all N_hat generated
  N_hat_sum <- rep(0, 2742)
  
  # repetition for 100 times
  for (i in 1:100) {
    
    # draw ui0 for each market, in total 2742 markets.
    ui0 <- rnorm(2742, 0, 1)
    
    
    # first we construct simulated profit pi_hat for each firm k
    
    # k=1
    pi_hat_AA <- cbind(data$marketdistance, data$marketsize) %*% beta + 
      cbind(data$mindistancefromhubAA, data$marketpresenceAA) %*% alpha -
      delta * data$lnN + rho * ui0 + sqrt((1-rho**2)) * rnorm(1)
    
    # k=2
    pi_hat_DL <- cbind(data$marketdistance, data$marketsize) %*% beta + 
      cbind(data$mindistancefromhubDL, data$marketpresenceDL) %*% alpha -
      delta * data$lnN + rho * ui0 + sqrt((1-rho**2)) * rnorm(1)
    
    # k=3
    pi_hat_UA <- cbind(data$marketdistance, data$marketsize) %*% beta + 
      cbind(data$mindistancefromhubUA, data$marketpresenceUA) %*% alpha -
      delta * data$lnN + rho * ui0 + sqrt((1-rho**2)) * rnorm(1)
    
    # k=4
    pi_hat_AL <- cbind(data$marketdistance, data$marketsize) %*% beta + 
      cbind(data$mindistancefromhubAL, data$marketpresenceAL) %*% alpha -
      delta * data$lnN + rho * ui0 + sqrt((1-rho**2)) * rnorm(1)
    
    # k=5
    pi_hat_LCC <- cbind(data$marketdistance, data$marketsize) %*% beta + 
      cbind(data$mindistancefromhubLCC, data$marketpresenceLCC) %*% alpha -
      delta * data$lnN + rho * ui0 + sqrt((1-rho**2)) * rnorm(1)
    
    # k=6
    pi_hat_WN <- cbind(data$marketdistance, data$marketsize) %*% beta + 
      cbind(data$mindistancefromhubWN, data$marketpresenceWN) %*% alpha -
      delta * data$lnN + rho * ui0 + sqrt((1-rho**2)) * rnorm(1)
    
    # calculate whether the simulated pi (pi_hat) is greater than 0
    # if pi_hat is > 0, then the firm will enter
    # sum all firms' decision
    N_hat_sum <- N_hat_sum + (pi_hat_AA>0)+(pi_hat_DL>0)+(pi_hat_UA>0)+(pi_hat_AL>0)+(pi_hat_LCC>0)+(pi_hat_WN>0)
    
  }
  # calculate mean predicted entrants
  N_hat <- N_hat_sum / 100
  
  # calculate predicted error
  pred_error <- data$totalN - N_hat
  
  # the obj. function will yield the prediction error sum of squared
  res <- t(data[, 8:27]) %*% (pred_error)
  ans <- t(res) %*% res
  return(-ans)
}

foo(params = rep(.5, 6))
library(maxLik)

# try optimization
t <- proc.time()
result <- maxLik(foo, start = rep(.5, 6), method = 'BFGSR')
proc.time() - t

initial <- result$estimate
# set the repetition to 100 times
# record each result in all_result
reps <- 100
all_result <- list()
initial <- rep(.5, 6)
for (s in 1:reps) {
  
  message("------Iteration: ", s, " processing...")
  t <- Sys.time()
  result <- maxLik(foo, start = initial, method = 'BFGSR')
  message("Time Elapsed: ", Sys.time() - t)
  
  all_result[[s]] <- result
  initial <- result$estimate
  
}

# use purrr::map to extract each coefficient
# save coefs. to each variables, e.g. beta1s, beta2s.
library(purrr)

params <- map(all_result, 'estimate')

beta1s <- map_dbl(params, 1)
beta2s <- map_dbl(params, 2)
alpha1s <- map_dbl(params, 3)
alpha2s <- map_dbl(params, 4)
deltas <- map_dbl(params, 5)
rhos <- map_dbl(params, 6)

# combine six parameters estimation to one data.table
dk <- cbind(beta1s, beta2s, alpha1s, alpha2s, deltas, rhos) %>% as.data.table()

# summary statistic for the estimation
# showing mean, st. dev ...,etc
library(stargazer)
stargazer(dk, type = 'text', 
          covariate.labels = c('mkt. dist.','mkt. size','mkt.hub dist.','mkt. presense','delta','rho'))













