set.seed(123)

# exercise 1

#####################
# 1.a generate normal rvs from cauchy

# find the maximum c
f = function(x){return(dnorm(x)/dcauchy(x))}

# perché ottimizzo proprio su questo intervallo?
c = optimize(f, c(-4,4), maximum=T)$objective

n = 100000
u = runif(n, 0 , 1)
y = rcauchy(n)
accept = ifelse(u < dnorm(y)/(c*dcauchy(y)), TRUE, FALSE)
x = y[accept]

# comparing the distributions
curve(dnorm(x), xlim=c(-4,4), col='blue', ylab="")
curve(dcauchy(x), add=T)


plot.new()
hist(x, freq=F)
curve(dnorm(x), add=T, col='red')

#################
# 1.b generate gamma random variables

# find the maximum c
f = function(x){return(dgamma(x, 4.3, 6.2)/dgamma(x, 4, 7))}

# perché ottimizzo proprio su questo intervallo?
c = optimize(f, c(0,5), maximum=T)$objective

u = runif(n, 0 , 1)
y = rgamma(n, 4, 7)
accept = ifelse(u < dgamma(y, 4.3, 6.2)/(c*dgamma(y, 4, 7)), TRUE, FALSE)
x = y[accept]

# comparing the distributions
curve(dgamma(x, 4.3, 6.2), xlim=c(0,5), col='blue', ylab="")
curve(dgamma(x, 4, 7), add=T)

plot.new()
hist(x, freq=F)
curve(dgamma(x, 4.3, 6.2), add=T, col='red')

####################
# exercise 2

# 2.a

metropolis <- function(cand, S, df){
  chain = c()
  chain[1] = cand
  R = 0
  for (i in 2:S+1){
    cand = chain[i]
    R = (dt(cand, df=df)*dnorm(chain[i]))/(dt(chain[i], df=df)*dnorm(cand))
    accept_prob = min(R, 1)
    if (runif(1) < accept_prob){
      chain[i+1] = cand
      R = R+1
    }else{
      chain[i+1] = chain[i]
    }
  }
  print(chain)
  return(chain)
}

metropolis(cand=rnorm(1,0,1), S=100, df=4)



