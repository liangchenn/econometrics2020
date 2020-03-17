


x <- revd(40000) - revd(40000)
plot(density(x), ylim = c(0, 1))
curve(exp(x)/(1+exp(x)**2), from = -6, to = 8, add = T)
