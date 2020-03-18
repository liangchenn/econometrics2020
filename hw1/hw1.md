hw1
================

# Problem Set 1

## Preprocessing

Load necessary packages:

``` r
library(extRemes)# for generating type-I random variable
library(maxLik)# for mle 
library(dplyr)# for pipeline operator `%>%`
library(purrr)# for map() 
```

First we set up parameters:

``` r
set.seed(123)
n <- 400
beta1 <- 1.0
beta2 <- -0.5

x1 <- rnorm(n, 0, 1)
x2 <- rchisq(n, df = 1)

u1 <- revd(n)
u2 <- revd(n)
```

## (1) data generating

``` r
y = ifelse(x1*beta1 + u1 > x2*beta2 + u2, 1, 0)
```

## (2) log-likelihood function

First, we define the log likelihood function,

taking `beta` = `{beta1, beta2}` as arguments and ouputing a numerical
vector of log-likelihood.

(Adding an option `type` to choose the output to be either vector or the
sum of log-likelihood value.)

``` r
# first, build a distribution function called `cdf`.
# using CDf of logit distribution, since subtraction of 2 type-I distribution is logit-distributed.
cdf <- function(x1, x2, beta1, beta2){
  
  u <- x1*beta1 - x2*beta2
  res <- exp(u) / (1+exp(u))
  
  return(res)
  
}


# define a function called `log_likelihood` using beta = {beta1, beta2} as arguments.
log_likelihood <- function(beta1, beta2, type = "vector"){
  
  res <- (y*log(cdf(x1, x2, beta1, beta2)) + (1-y)*log(1 - cdf(x1, x2, beta1, beta2)))
  
  if(type == 'scalar') return(sum(res)) else return(res)
  
}
```

## (3) grid search

generating a grid of \(\beta = \{\beta_1, \beta_2\}\),

here we define a function to expand the grid of parameters and report
the maximum likelihood realisation row.

``` r
get_mle_betas <- function(from = -5, to = 5, step = 0.1){
  # expand grid with step given
  grids <- expand.grid(seq(from, to, step), seq(from, to, step))
  
  # calculate log-likelihood(llk) ratio for each set of beta
  # use mapply function to map log_likelihood function rowwise.
  grids$llk <- mapply(log_likelihood, grids$Var1, grids$Var2, 'scalar')
  
  # get index of the maximum value of log likelihood ratio
  index <- which.max(grids$llk)  
  betas <- grids[index, ]
  # returning betas
  return(betas)
}
```

If choosing step = 0.1.

``` r
get_mle_betas()
```

    ##      Var1 Var2       llk
    ## 4607  1.1 -0.5 -224.4289

Narrowing down with step = 0.01

``` r
get_mle_betas(from = -1, to = 1, step = 0.01)
```

    ##       Var1 Var2       llk
    ## 10251    1 -0.5 -224.4786

## (4) other gradient/ non-gradient method

Here we use `BHHH` and `Nelder Mead` methods provided by `maxLik`
package.

\_(Ref : <https://github.com/cran/maxLik)_>

### BHHH

first, we rewirte log likelihood function to meet the requirement of
`maxLik` function.

``` r
log_likelihood <- function(beta){
  beta1 <- beta[1]
  beta2 <- beta[2]
  
  res <- (y*log(cdf(x1, x2, beta1, beta2)) + (1-y)*log(1 - cdf(x1, x2, beta1, beta2)))

  return(res)
}
```

``` r
BHHH <- maxLik(log_likelihood, start = c(0, 0), method = 'BHHH')
BHHH$estimate
```

    ## [1]  1.0609732 -0.5057706

Generating R = 100, N = 400 sample.

``` r
# R = 100, N = 400
t <- proc.time()
betas <- list()
for (i in 1:100) {
  # setups
  n <- 400
  beta1 <- 1.0
  beta2 <- -0.5
  x1 <- rnorm(n, 0, 1)
  x2 <- rchisq(n, df = 1)
  u1 <- revd(n)
  u2 <- revd(n)
  y <- ifelse(x1*beta1 + u1 > x2*beta2 + u2, 1, 0)
  
  # BHHH estimation
  BHHH <- maxLik(log_likelihood, start = c(0, 0), method = 'BHHH')
  
  betas[[i]] <- BHHH$estimate
  # message(i, '-th simulation...')
}
proc.time() - t
```

    ##    user  system elapsed 
    ##   0.548   0.034   0.586

``` r
# the computation time: 1.087s
```

Define a function to show the final result, including estimated beta and
its standard deviation.

``` r
show_results <- function(betas) {
  beta1s <- map_dbl(.x = betas, .f = 1)
  beta2s <- map_dbl(betas, 2)
  
  message("The mean of beta 1 is: ", mean(beta1s), " , and the sd is: ", sd(beta1s))
  message("The mean of beta 2 is: ", mean(beta2s), " , and the sd is: ", sd(beta2s))
}
```

``` r
show_results(betas)
```

    ## The mean of beta 1 is: 0.99400702590802 , and the sd is: 0.128792152916141

    ## The mean of beta 2 is: -0.50031288370689 , and the sd is: 0.091142855808259

### Nelder Mead Method

Same for NM :

``` r
# R = 100, N = 400
t <- proc.time()
betas_nm <- list()
for (i in 1:100) {
  # setups
  n <- 400
  beta1 <- 1.0
  beta2 <- -0.5
  x1 <- rnorm(n, 0, 1)
  x2 <- rchisq(n, df = 1)
  u1 <- revd(n)
  u2 <- revd(n)
  y <- ifelse(x1*beta1 + u1 > x2*beta2 + u2, 1, 0)
  
  # BHHH estimation
  NM <- maxLik(log_likelihood, start = c(0, 0), method = 'NM')
  
  betas_nm[[i]] <- NM$estimate
  # message(i, '-th simulation...')
}
proc.time() - t
```

    ##    user  system elapsed 
    ##   0.609   0.036   0.651

``` r
show_results(betas_nm)
```

    ## The mean of beta 1 is: 1.00602680797551 , and the sd is: 0.13462494456679

    ## The mean of beta 2 is: -0.509983180077846 , and the sd is: 0.0952679195463455
