# a test file for vectorization by cppcodes.cpp
Rcpp::sourceCpp(file = './hw1/cppcodes.cpp', verbose = T)


tester_llk <- function(beta1, beta2){
  
  
  res <- C_F(x1, x2, y, beta1, beta2)
  
  return(res)
}


from <- -5
to   <- 5
step <- 0.01
grids <- expand.grid(seq(from, to, step), seq(from, to, step))


t <- proc.time()
grids$llk <- mapply(tester_llk, grids$Var1, grids$Var2)
proc.time() - t



t <- proc.time()
grids$llk  <- CppLoglikelihood(grids$Var1, grids$Var2)
proc.time() - t


CppLoglikelihood(-5:5, -5:5)

t <- proc.time()
a <- mapply(tester_llk, -5000:5000, -5000:5000)
proc.time() - t


t <- proc.time()
a <- CppLoglikelihood(-5000:5000, -5000:5000)
proc.time() - t


