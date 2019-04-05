set.seed(123)

y20 <- rpois(20, 0.5)
y100 <- rpois(100, 0.5)
y1000 <- rpois(1000, 0.5)

a <- c(1, 1, 1, 2)
b <- c(.5, 2, 10, 2)

# 1.
gamma_prior <- function(x, k){return(dgamma(x, a[k], b[k]))}

jeffreys_prior <- function(x, theta=0.5){return(1/sqrt(x))}

for(k in 1:4){
	curve(gamma_prior(x, k), col=k, ylim=c(0,5), add=TRUE)
}
curve(jeffreys_prior(x), add=TRUE, col='yellow')

# 2.

gamma_posterior <- function(x, y, k){
	return(dgamma(x, sum(y)+a[k], length(y)+b[k]))
}

jeffreys_posterior <-  function(x, theta=0.5){
	return(dgamma(x, sum(y)+0.5, length(y)))
}

data <- c(y20, y100, y1000)

plot_data <- function(y){
  for(k in 1:4){
    curve(gamma_posterior(x, y20, k), col=k, ylim=c(0,5), add=TRUE)
  }
	curve(jeffreys_posterior(x, y20), add=TRUE, col='yellow')
}




# 3.




# 4.

# 5.


