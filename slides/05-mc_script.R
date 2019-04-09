###############################
# AR method
###############################

## simulation from beta density

# 1)

#number of simulations
Nsims=2500
#beta parameters
a=2.7; b=6.3
#find c via optimization
c=optimise(f=function(x) {dbeta(x,a,b)},
           interval=c(0,1), maximum=TRUE)$objective
u=runif(Nsims, max=c)
theta_star=runif(Nsims)
theta=theta_star[u<dbeta(theta_star,a,b)]
#accept probability
1/c

#plots
pdf(file="ar_beta.pdf", width=11, height =6)
par(mfrow=c(1,2), mar =c(5,5,2,1))
hist(theta, prob=T, xlab=expression(theta), ylim=c(0,3),
     main ="Simulation from Beta(2.7,6.3)", cex.lab =1.9)
curve(dbeta(x,a,b), col="black",lwd=4, add=T)
curve(dbeta(x,a,b), col="black",lwd=4, xlab=expression(theta), 
      ylab="Density", ylim=c(0,3), main ="Simulation from Beta(2.7,6.3)", 
      cex.lab=1.9)
points(theta, u[u<dbeta(theta_star,a,b)], col="darkgreen")
points(theta_star[u>=dbeta(theta_star,a,b)],
       u[u>=dbeta(theta_star,a,b)], col="red")
dev.off()

# 2)

Nsims=2500
#beta parameters
a=2; b=3
#find c via optimization
c=optimise(f=function(x) {dbeta(x,a,b)},
           interval=c(0,1), maximum=TRUE)$objective
u=runif(Nsims, max=c)
theta_star=runif(Nsims)
theta=theta_star[u<dbeta(theta_star,a,b)]
#accept probability
1/c

#plots
pdf(file="ar_beta2.pdf", width=11, height =6)
par(mfrow=c(1,2), mar =c(5,5,2,1))
hist(theta, prob=T, xlab=expression(theta), ylim=c(0,3),
     main ="Simulation from Beta(2,3)", cex.lab =1.9)
curve(dbeta(x,a,b), col="black",lwd=4, add=T)
curve(dbeta(x,a,b), col="black",lwd=4, xlab=expression(theta), 
      ylab="Density", ylim=c(0,3), main ="Simulation from Beta(2,3)", 
      cex.lab=1.9)
points(theta, u[u<dbeta(theta_star,a,b)], col="darkgreen")
points(theta_star[u>=dbeta(theta_star,a,b)],
       u[u>=dbeta(theta_star,a,b)], col="red")
dev.off()

#############################
## CLASSICAL MONTE CARLO
############################


## Normal mean with cauchy prior

set.seed(12345)
theta = rnorm(5000, 10, 1)
I = sum(theta/(1 + theta^2))/sum(1/(1 + theta^2))
I


par(mar=c(5,6,4,1))
pdf(file="mc_normal.pdf", height =7, width = 10)
plot(cumsum(theta/(1 + theta^2))/cumsum(1/(1 +theta^2)), type = "l",
      ylab = "Posterior mean",
     xlab = "iterations", cex.lab =1.6)
abline(h = I, col="blue", lwd=3)
#abline(h=I2)
dev.off()

#############################
## IMPORTANCE SAMPLING
############################


## Location for a t9
y.t = rt(n = 9, df = 3)
t.lik = function(theta, data) prod((3 + (data - theta)^2)^(-(3 + 1)/2))
pdf(file ="is_t.pdf", width =8, height =7)
curve(sapply(x, function(x) t.lik(x, data = y.t)),
        from = -3, to = 3, xlab = expression(theta), ylab = "", 
      cex.lab = 2, lwd=3, main = "Posterior for the location of student t",
      cex.main =1.4)
dev.off()


t.medpost = function(nsim, data, l) {
   sim <- data[l] + rt(nsim, 3)
   n <- length(data)
   s <- c(1:n)[-l]
   num <- cumsum(sim * sapply(sim, function(theta) t.lik(theta, data[s])))
   den <- cumsum(sapply(sim, function(theta) t.lik(theta, data[s])))
   num/den
   }
 media.post = t.medpost(nsim = 1500, data = y.t,
                         l = which(y.t == median(y.t)))
 media.post[1500]

 
 pdf(file ="is_t_conv.pdf", width =8, height =7)
 plot(media.post, type = "l", ylim = c(-1, 1), xlab ="Iterations", cex.lab =1.5,
      ylab ="Posterior mean")
 dev.off()

 
 pdf(file="is_t_other_g.pdf", width =10, height =7)
  par(mfrow = c(1, 2))
  plot(c(0, 0), xlim = c(0, 1000), 
       ylim = c(-0.75,0.75), type = "n", ylab = "Posterior mean",
       xlab="Iterations", main = bquote(g[1](theta)), cex.lab=1.8, cex.main =2.4)
  for (i in 1:10) {
    lines(x = c(1:1000), y = t.medpost(nsim = 1000, 
          data = y.t, l = which(y.t == median(y.t))),
           col = 3)
    }
  plot(c(0, 0), xlim = c(0, 1000), ylim = c(-0.75,  0.75), type = "n", ylab = "Posterior mean",
       xlab ="Iterations", main = bquote(g[2](theta)),cex.lab=1.8, cex.main =2.4)
  for (i in 1:10) {
    lines(x = c(1:1000), y = t.medpost(nsim = 1000,
                        data = y.t, l = which(y.t == max(y.t))),
           col = 3)
  }
  dev.off()



