
probit <- function(x){
  
  (1/sqrt(2*pi))*exp((-.5)*x**2)
  
}

x <- revd(400000) - revd(400000)
plot(density(x), ylim = c(0, 1), xlab = 'X', main = 'PDF')
curve(exp(x)/(1+exp(x)**2), from = -6, to = 8, add = T, col = 'blue')# logit 
curve(probit, from = -6, to = 8, add = T, col = 'red')# probit

probit <- function(x){
  
  (1/sqrt(2*pi))*exp((-.5)*x**2)
  
}

