library(ggplot2)

target <- function(theta){
  
  dnorm(theta, mean = 20, sd = 1)
  
}


n = 50000


sigma = 1
reject = 0
pi <- vector("numeric", (n+1))

for(i in 1:(n+1)){
  
  x <- rnorm(1, mean = pi[i], sd = sigma)
  
  alpha <- min(1, (target(x))/(target(pi[i-1])))
  
  u <- runif(1)
  
  if(u < alpha){
    pi[i+1] <- x
  }else{
    pi[i+1] <- pi[i]
    reject = reject + 1
  }
  
}


reject/n

qplot(pi, geom = 'density')


plot(density(pi))
curve(dnorm(x, 20, 1), -30, 30, add = T, col = 'red')

plot(pi[40000:50000], type = 'l')
