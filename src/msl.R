# maximum simulated likelihood

data <- read.table('./demo/mma12p2mslmsm.asc')
u <- data[, 1]
e <- data[, 2]
y <- data[, 3]
obs <- length(u)
simulations <- 1000



