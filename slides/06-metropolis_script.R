#############################
## METROPOLIS-HASTINGS
############################

## Gamma distribution: rw

mh.gamma <-
  function(n.sims=10^5,
           start=1, burnin=0.3*n.sims, cand.sd=2,
           shape=1, rate=1) {
    # draws values from a gamma density
    # with shape and rate
    # sampling from a normal with mean
    # equal to the last vaue and sd = sd
    theta.cur <- start
    draws <- c()
    theta.update <- function(theta.cur, shape, rate) {
      theta.can <- rnorm(1, mean = theta.cur, sd = cand.sd)
      accept.prob <- dgamma(theta.can, shape = shape, rate = rate)/
        dgamma(theta.cur, shape = shape, rate = rate)
      if (runif(1) <= accept.prob)
        theta.can
      else theta.cur
    }
    for (i in 1:n.sims) {
      draws[i] <- theta.cur <- theta.update(theta.cur,
                                            shape = shape, rate = rate)}
    samp<-draws[(burnin+1):n.sims]
    par(mfrow=c(1,2), mar=c(5,5,4,1))
    hist(samp, nclass=50, prob=T,
         main="Posterior density", 
         col="grey", xlab=expression(theta), cex.lab=1.6, cex.main=2)
    curve(dgamma(x, shape=shape, rate=rate), add=T, lwd =3, 
          col ="red")
    plot(samp, type="l", ylab=expression(theta),
         xlab="Iterations", cex.lab=1.6, cex.main=2, main="Markov chain")
    return(samp)}
pdf(file="gamma_rw.pdf", width=12, height=6.5)
a=mh.gamma(shape=4.4, rate=1.7)
dev.off()


## Beta distribution (2.7, 6.3): independent

mh.beta <-
  function(n.sims=10^5,
           start=0.2, burnin=0.3*n.sims,
           lower =0,
           upper =3,
           shape1=2.7, shape2=6.3) {
    # draws values from a gamma density
    # with shape and rate
    # sampling from a normal with mean
    # equal to the last vaue and sd = sd
    theta.cur <- start
    draws <- c()
    draws[1]<- start
    cont<-0
    
      
      for (i in 2:n.sims) {
        theta.can <- runif(1, lower, upper)
      accept.prob <- dbeta(theta.can, shape1 = shape1, shape2 = shape2)*dunif(draws[i-1], lower, upper)/
        dbeta(draws[i-1], shape1 = shape1, shape2 = shape2)*dunif(theta.can, lower, upper )
      if (runif(1) <= accept.prob){
        draws[i] <- theta.can
      cont <- cont +1
      }else{ draws[i] <- draws[i-1]
      }
    }
      
      accept_rate <- cont/n.sims
    
    
    samp<-draws[(burnin+1):n.sims]
    par(mfrow=c(1,2), mar=c(5,5,4,1))
    hist(samp, nclass=50, prob=T,
         main="Posterior density", 
         col="grey", xlab=expression(theta), cex.lab=1.6, cex.main=2)
    curve(dbeta(x, shape1=shape1, shape2=shape2), add=T, lwd =3, 
          col ="red")
    plot(samp, type="l", ylab=expression(theta),
         xlab="Iterations", cex.lab=1.6, cex.main=2, main="Markov chain")
    return(c(samp, cont))}

pdf(file="beta_ind2.pdf",  width=12, height=6.5)
b=mh.beta(shape1=2.7, shape2=6.3)
dev.off()


## Shuttle data

temp <- c(66, 70, 69, 68, 67, 72, 73, 70, 57, 63, 70, 78, 67, 53,
           67, 75, 70, 81, 76, 79, 75, 76, 58)
failure <- c(0, 1, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 1, 0, 0, 0,
                0, 0, 0, 1, 0, 1)

#log-verosimiglianza


loglik <- function(beta) {
  sum(failure * (beta[1] + beta[2] * temp)) -
    sum(7 * log(1 + exp(beta[1] + beta[2] *
                          temp)))
}

#Grafico della verosimiglianza

beta0=seq(-5,15,length=200)
beta1=seq(-0.25,0.05,length=200)
l.val=apply(expand.grid(beta0,beta1),FUN=loglik,MAR=1)
l.val=matrix(l.val,ncol=length(beta1),nrow=length(beta0))
par(mfrow=c(1,1), mar=c(5,2,2,1))
pdf("shuttle_loglik.pdf", width =10, height=8)
contour(beta0,beta1,exp(l.val-max(l.val)),
        levels=c(0.01,0.1,0.5,0.9),xlab=expression(beta[0]),
        ylab=expression(beta[1]),cex.lab=1,main="relative likelihood", cex.main =2,
        cex.lab=2)
dev.off()


# 1) random walk MH with univariate proposals

mh <- function(nsim, s1, s2, b.init) {
  mh.out <- matrix(ncol = 2, nrow = nsim)
   b <- b.init
   for (i in 1:nsim) {
     b.p <- c(rnorm(1, b[1], s1), rnorm(1,b[2], s2))
    if (runif(1) < exp(loglik(b.p) - loglik(b)))
      b <- b.p
      mh.out[i, ] <- b
       }
  mh.out}

par(mfrow = c(2, 2))
pdf("shuttle1.pdf", width=8, height =7)
mh.out1 <- mh(nsim = 5000, s1 = 1, s2 = 1, b.init = c(0,0))
par(mfrow = c(2, 2), mar=c(5,5,2,1))
plot(mh.out1[, 1],type="l", ylab=expression(beta[0]), xlab="Iterations", cex.lab=2)
plot(mh.out1[, 2],type="l", ylab=expression(beta[1]), xlab="Iterations", cex.lab =2)
contour(beta0,beta1,exp(l.val-max(l.val)),levels=c(0.01,0.1,0.5,0.9), 
        xlab=expression(beta[0]), 
        ylab=expression(beta[1]), xlim=c(-10,15), 
        ylim=c(-0.3,0.1), cex.lab=2)
points(mh.out1,col=2)
dev.off()

# 2) random wal with univariate proposals, different variances

 par(mfrow = c(2, 2))
 mh.out2=mh(nsim = 5000, s1 = 2, s2 = 0.1, b.init = c(0,0))
 par(mfrow = c(2, 2))
 pdf("shuttle2.pdf", width=8, height =7)
 par(mfrow = c(2, 2), mar=c(5,5,2,1))
 plot(mh.out2[, 1],type="l", ylab=expression(beta[0]), xlab="Iterations", cex.lab=2)
 plot(mh.out2[, 2],type="l", ylab=expression(beta[1]), xlab="Iterations", cex.lab =2)
 contour(beta0,beta1,exp(l.val-max(l.val)),levels=c(0.01,0.1,0.5,0.9), 
         xlab=expression(beta[0]), 
         ylab=expression(beta[1]), xlim=c(-10,15), 
         ylim=c(-0.3,0.1), cex.lab=2)
 points(mh.out2,col=2)
 dev.off()
 
 # 3) random wal with univariate proposals, estimated variances
 
 par(mfrow = c(2, 2))
 mh.out3=mh(nsim = 5000, s1 = sd(mh.out2[,1]), s2 =sd(mh.out2[,2]), b.init = c(0,0))
 par(mfrow = c(2, 2))
 pdf("shuttle3.pdf", width=8, height =7)
 par(mfrow = c(2, 2), mar=c(5,5,2,1))
 plot(mh.out3[, 1],type="l", ylab=expression(beta[0]), xlab ="Iterations", cex.lab=2)
 plot(mh.out3[, 2],type="l", ylab=expression(beta[1]), xlab ="Iterations", cex.lab=2)
 contour(beta0,beta1,exp(l.val-max(l.val)),levels=c(0.01,0.1,0.5,0.9),
         cex.lab=2, xlab=expression(beta[0]),
         ylab=expression(beta[1]))
 points(mh.out3,col=2)
 dev.off()
 

# 4) random walk with bivariate normal and correlation structure

 library(mvtnorm)
 mhd <- function(nsim, V, b.init) {
   mh.out <- matrix(ncol = 2, nrow = nsim)
   b <- b.init
   for (i in 1:nsim) {
     b.p <- rmvnorm(n = 1, mean = b, sigma = V)
     if (runif(1) < exp(loglik(b.p) - loglik(b)))
       b <- b.p
       mh.out[i, ] <- b
       }
   mh.out}

par(mfrow = c(2, 2), oma=c(0,0,0,0))
v <- var(mh.out3)
par(mfrow=c(2,2))
 mh.out4 <- mhd(nsim = 10000, V = v, b.init = c(0,0))
 pdf("shuttle4.pdf", width=8, height =7)
 par(mfrow=c(2,2), mar=c(5,5,2,1))
 plot(mh.out4[, 1],type="l", ylab =expression(beta[0]), xlab = "Iterations", cex.lab =2)
 plot(mh.out4[, 2],type="l", ylab=expression(beta[1]), xlab="Iterations", cex.lab =2)
 contour(beta0,beta1,exp(l.val-max(l.val)),levels=c(0.01,0.1,0.5,0.9),
         xlab =expression(beta[0]), ylab=expression(beta[1]), cex.lab=2)
 points(mh.out4,col=2)
dev.off()


post_mean <- apply(mh.out4,2,mean)
pdf(file="shuttle_est.pdf", width=12, height=5)
par(mfrow = c(1, 2), mar=c(5,5,2,1))
 plot(density(mh.out4[, 1], bw = 0.5), xlab = expression(beta[0]),main="", cex.lab =2)
 plot(density(mh.out4[, 2], bw = 0.01), xlab = expression(beta[1]),main="", cex.lab=2)
 dev.off()

 apply(mh.out4, 2, mean)
 
 
 ############################
 ## GIBBS SAMPLING
 ###########################
 
 ## bivariate normal
 
 par(mfrow=c(1,1))
 
 bn<-function(init=c(0,0),mu=c(0,0),
              rho=0, nsim=100){
   theta=matrix(nrow=nsim, ncol=2)
   theta[1,]=init
   for (i in 2:nsim){
     theta[i,1]<-rnorm(1,mean= mu[1]+ rho*(theta[i-1,2]-mu[2]), sd=sqrt(1-rho^2))
     theta[i,2]<-rnorm(1,mean= mu[2]+ rho*(theta[i,1]- mu[1]), sd=sqrt(1-rho^2))
   }
   #plot(theta, type="l", xlab=expression(theta[1]),
        #ylab=expression(theta[2]), cex.lab =2, cex.main=2, main= bquote(rho==.(rho)))
   
   
   # 2-dimensional plot
   xy <- theta
   nbins <- 20
   x.bin <- seq(floor(min(xy[,1])),
                ceiling(max(xy[,1])), length=nbins)
   y.bin <- seq(floor(min(xy[,2])),
                ceiling(max(xy[,2])), length=nbins)
   freq <-  as.data.frame(table(findInterval(xy[,1],
                                             x.bin),findInterval(xy[,2], y.bin)))
   freq[,1] <- as.numeric(freq[,1])
   freq[,2] <- as.numeric(freq[,2])
   freq2D <- diag(nbins)*0
   freq2D[cbind(freq[,1], freq[,2])] <- freq[,3]
   res <- persp(x.bin, y.bin,
                freq2D, theta=30, phi=30, xlab = expression(theta[1]),
                ylab = expression(theta[2]), zlab="Density",
                expand=0.5, ltheta=120,
                col = "lightblue",
                shade = 0.1, ticktype = "detailed",
                main= bquote(rho==.(rho)), cex.main=1.5, cex.lab=1.5)
   
   
   }
 
 pdf(file="gibbs_bivn3.pdf", width=10, height =7.5)
 par(mfrow=c(2,2), mar =c(5,5,2,1))
 my.draws<-bn(rho=0.2)
 my.draws<-bn(rho=0.5)
 my.draws<-bn(rho=0.7)
 my.draws<-bn(rho=0.9)
 dev.off()
 
 pdf(file="gibbs_bivn4.pdf", width=10, height =7.5)
 par(mfrow=c(2,2), mar =c(5,5,2,1))
 my.draws<-bn(rho=0.2, nsim =10000)
 my.draws<-bn(rho=0.5, nsim=10000)
 my.draws<-bn(rho=0.7, nsim=10000)
 my.draws<-bn(rho=0.9, nsim=10000)
 dev.off()
 
 ## Nuclear pump
 
  y <- c(5, 1, 5, 14, 3, 19, 1, 1, 4, 22)
  t <- c(94, 16, 63, 126, 5, 31, 1, 1, 2, 10)
  rbind(y, t)
  
  
  
  
  ###########################
  ## MCMC DIAGNOSTICS
  ##########################
  
  
  
bm<-function(init,mu=c(0,0), rho, nsim=100){  
    #inizializzo funzione gibbs
    theta=matrix(nrow=nsim, ncol=2)
    theta[1,]=init  # inizializzo theta
    
    for (i in 2:nsim){
      theta[i,1]=rnorm(1,mu[1]+
                         rho*(theta[i-1,2]-mu[2]),sqrt(1-(rho)^2))
      theta[i,2]=rnorm(1, mu[2]+
                         rho*(theta[i,1]-mu[1]), sqrt(1-(rho)^2))
      
      
    }
    
    
    return(theta[(.5*nsim):nsim,])
  }
  
  
  ## fake example of graphical inspection: just for producing the
  ## plots 
  
  
  fake_chains1 <- bm(init=c(0,0), rho=0.7, nsim =3000)
  fake_chains2 <- bm(init=c(0,0), rho=0.7, nsim =3000)
  
  pdf(file="fake_conv.pdf", width=11.5, height =7)
  par(mfrow=c(1,2), mar=c(5,5,2,1))
  plot(fake_chains1[,1], type="l", col="blue", ylim=c(-8,8),
       ylab =expression(mu[1]), main = "Scenario 1", cex.main =2, cex.lab=2,
       xlab="Iterations")
  lines(fake_chains1[,2]+4, col="red")
  legend(1000,-6, col=c("red", "blue"), lty=1, 
         lwd =3, c("chain 1", "chain 2"))
  plot(fake_chains1[,1]+log(seq(1:1501)), type="l", col="blue", ylim=c(-2,12),
       ylab =expression(mu[1]), main = "Scenario 2", cex.main =2, cex.lab=2,
       xlab="Iterations")
  lines(fake_chains1[,2]+10-log(seq(1:1501)) , col="red")
  legend(1000,-0.5, col=c("red", "blue"), lty=1, 
         lwd =3, c("chain 1", "chain 2"))
  dev.off()
  
  pdf(file="fake_conv_hist.pdf", width=11.5, height =7)
  par(mfrow=c(1,2), mar=c(5,5,2,1))
  # hist(c(fake_chains1[,1], fake_chains1[,2]+4), cex.lab=2, 
  #      xlab =expression(mu[1]), probability =TRUE, breaks=30, 
  #      main ="Scenario 1", xlim=c(-5,6))
  hist(c(fake_chains1[,1]), cex.lab=2, col="blue",
       xlab =expression(mu[1]), probability =TRUE, breaks=30, 
       main ="Scenario 1", xlim=c(-5,6), cex.main=2)
  hist(fake_chains1[,2]+4, breaks=30, add=TRUE, col="red",
       probability=TRUE)
  hist(c(fake_chains1[,1]+log(seq(1:1501)) ), 
       cex.lab=2, cex.main=2,
       xlab =expression(mu[1]), probability =TRUE, breaks=40, 
       main ="Scenario 2", col="blue")
 hist(fake_chains1[,2]+10-log(seq(1:1501)), cex.lab=2, 
      xlab =expression(mu[1]), probability =TRUE, breaks=40, add =TRUE, col="red")
  dev.off()
  
  
  
  
  
  
  ## multiple chains
  
  # 100 iterations
  
  my.draws<-bm(init=c(0,0),rho=0.7)
  my.draws2<-bm(init=c(5,5), rho=0.7)
  my.draws3<-bm(init=c(10,10), rho=0.7)
  
  
     # convergence to stationarity
  
  par(mfrow=c(2,2))
  pdf("multiplechains_bivn.pdf", width=8, height =7)
  par(mfrow=c(2,2), mar=c(5,5,2,1))
  plot(my.draws, type="l", xlab=expression(theta[1]), 
       ylab=expression(theta[2]), cex.lab=2)
  lines(my.draws2, type="l", col="red", 
        xlab=expression(theta[1]), ylab=expression(theta[2]))
  lines(my.draws3, type="l", col="green", 
        xlab=expression(theta[1]), ylab=expression(theta[2]))
  
  plot(my.draws[,1], col="black",type="l", 
       ylab=expression(theta[1]), xlab="Iterations", ylim=c(-3,3), cex.lab=2)
  lines(1:length(my.draws[,1]),my.draws2[,1], col="red",type="l")
  lines(1:length(my.draws[,1]),my.draws3[,1], col="green",type="l")
  
  plot(my.draws[,2], col="black",type="l", 
       ylab=expression(theta[2]), xlab="Iterations", ylim=c(-3,3),
       cex.lab=2)
  lines(1:length(my.draws[,2]),my.draws2[,2], col="red",type="l")
  lines(1:length(my.draws[,2]),my.draws3[,2], col="green",type="l")
  dev.off()
  
  
      # iid samples
  
  
  pdf(file="multiplechains_bivn_acf.pdf", width=8, height =7)
  par(mfrow=c(3,2), mar=c(5,4,5,1))
  acf(my.draws[,1], main =(bquote(theta[1]~": Chain 1")), cex.main = 3)
  acf(my.draws[,2], main =(bquote(theta[2]~": Chain 1")), cex.main = 3)
  acf(my.draws2[,1],main =(bquote(theta[1]~": Chain 2")), cex.main = 3)
  acf(my.draws2[,2],main =(bquote(theta[2]~": Chain 2")), cex.main = 3)
  acf(my.draws3[,1],main =(bquote(theta[1]~": Chain 3")), cex.main = 3)
  acf(my.draws3[,2],main =(bquote(theta[2]~": Chain 3")), cex.main = 3)
  dev.off()
  
  # 10000 iterations
  
  my.draws<-bm(init=c(0,0),rho=0.7, nsim=10000)
  my.draws2<-bm(init=c(5,5), rho=0.7, nsim=10000)
  my.draws3<-bm(init=c(10,10), rho=0.7, nsim=10000)
  
  
        # convergence
  
  par(mfrow=c(2,2))
  pdf("multiplechains_bivn2.pdf", width=8, height =7)
  par(mfrow=c(2,2), mar=c(5,5,2,1))
  plot(my.draws, type="l", xlab=expression(theta[1]), 
       ylab=expression(theta[2]), cex.lab=2)
  lines(my.draws2, type="l", col="red", 
        xlab=expression(theta[1]), ylab=expression(theta[2]))
  lines(my.draws3, type="l", col="green", 
        xlab=expression(theta[1]), ylab=expression(theta[2]))
  
  plot(my.draws[,1], col="black",type="l", 
       ylab=expression(theta[1]), xlab="Iterations", ylim=c(-3,3), cex.lab=2)
  lines(1:length(my.draws[,1]),my.draws2[,1], col="red",type="l")
  lines(1:length(my.draws[,1]),my.draws3[,1], col="green",type="l")
  
  plot(my.draws[,2], col="black",type="l", 
       ylab=expression(theta[2]), xlab="Iterations", ylim=c(-3,3),
       cex.lab=2)
  lines(1:length(my.draws[,2]),my.draws2[,2], col="red",type="l")
  lines(1:length(my.draws[,2]),my.draws3[,2], col="green",type="l")
  dev.off()
  
  
          # iid samples
  
  pdf(file="multiplechains_bivn2_acf.pdf", width=8, height =7)
  par(mfrow=c(3,2), mar=c(5,4,5,1))
  acf(my.draws[,1], main =(bquote(theta[1]~": Chain 1")), cex.main = 3)
  acf(my.draws[,2], main =(bquote(theta[2]~": Chain 1")), cex.main = 3)
  acf(my.draws2[,1],main =(bquote(theta[1]~": Chain 2")), cex.main = 3)
  acf(my.draws2[,2],main =(bquote(theta[2]~": Chain 2")), cex.main = 3)
  acf(my.draws3[,1],main =(bquote(theta[1]~": Chain 3")), cex.main = 3)
  acf(my.draws3[,2],main =(bquote(theta[2]~": Chain 3")), cex.main = 3)
  dev.off()
  
  #Diagnostics: Gelman-Rubin statistic (coda package)
  
  par(mfrow=c(1,1))
  library(coda)
  chain1 <- mcmc(my.draws)
  chain2 <- mcmc(my.draws2)
  chain3 <- mcmc(my.draws3)
  combinechains=mcmc.list(chain1, chain2, chain3)
  
      # gelman-rubin
 
  gelman.diag(combinechains)
  pdf(file = "multiplechains_biv2_gr.pdf", width=10, height=7)
  par(mar=c(5,5,2,1))
  gelman.plot(combinechains)
  dev.off()
  
    # geweke
  
  geweke.diag(chain1)
  geweke.diag(chain2)
  geweke.diag(chain3)
  