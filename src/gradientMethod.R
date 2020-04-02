library(maxLik)

loglik <- function(param){
    
    mu <- param[1]
    sigma <- param[2]
    ll <- -0.5*N*log(2*pi) - N*log(sigma) - sum(0.5*(x - mu)^2/sigma^2)

    return(ll)

}

x <- rnorm(1000, 1, 2) # use mean=1, stdd=2
N <- length(x)

res <- maxNR(loglik, start=c(0,1)) # use Newton-Raphson 
summary(res)

res <- maxLik(loglik, start=c(0, 1))
summary(res)
