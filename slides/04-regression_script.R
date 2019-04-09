## data

head(cars)
n <- dim(cars)[1]
y <- cars$dist
x <- cars$speed
pdf(file="cars.pdf", height=6, width=7)
plot(x, y, xlab = "Speed", ylab="Distance", main ="cars data")
lm_est <- lm(y~x)
lines(x, lm_est$fitted.values , col="red", lwd =2)
dev.off()

p <- 1
X <- matrix(NA, n,2)
X[,1] <- rep(1,n)
X[,2] <- cars$speed
beta_hat <- solve(t(X)%*%X)%*%t(X)%*%y
s2 <- ( n/(n-p)*1 /(n-p))*t(y-X%*%beta_hat)%*%(y-X%*%beta_hat )


###########################
## NONINFORMATIVE ANALYSIS
###########################


#install.packages("invgamma")
library(invgamma)
nsims <- 10^3

# marginal for sigma2 


sigma2_sim <- rinvgamma(nsims, (n-p)/2, (n-p)*s2/2)
pdf(file="marg_post.pdf", width =7, height =6)
hist(sigma2_sim, probability = TRUE, 
     xlab =expression(sigma^2), main =bquote("Marginal posterior for" ~ sigma^2  ), 
     breaks= 30, ylim=c(0,0.01), cex.main =1.5)
curve(dinvgamma(x, (n-p)/2, (n-p)*s2/2 ), col="blue", lwd=3, add=TRUE)
abline(v = s2, col ="red", lwd =3, lty=2)
text(s2+10, 0.009, expression(s^2), col="red", cex =1.3)
legend(400, 0.008, lty=c(1,2), lwd=c(3,3), col=c("blue", "red"), c("Posterior", "MLE est."))
dev.off()

# conditional for beta

library(mvtnorm)
 beta <- matrix(NA, nsims,2)
for (i in 1:nsims){
 beta[i,] <- rmvnorm(1, beta_hat, solve(t(X)%*%X)*sigma2_sim[i])
 }

pdf(file="cond_post_alpha.pdf", width =7, height =6)
hist(beta[,1], probability = TRUE, 
     xlab =expression(alpha), main =bquote("Conditional posterior for" ~ alpha  ), 
     breaks= 30, ylim= c(0, 0.07),  cex.main =1.5, cex.lab =2)
curve(dnorm(x, beta_hat[1], sqrt(solve(t(X)%*%X)[1,1]*as.double(s2))), col="blue", lwd=3, add=TRUE)
abline(v = beta_hat[1], col ="red", lwd =3, lty=2)
text(beta_hat[1]+2, 0.065, expression(hat(alpha)), col="red", cex=1.3)
legend(-5, 0.06, lty=c(1,2), lwd=c(3,3), col=c("blue", "red"), c("Posterior", "MLE est."))
dev.off()

pdf(file="cond_post_beta.pdf", width =7, height =6)
hist(beta[,2], probability = TRUE, 
     xlab =expression(gamma), main =bquote("Conditional posterior for" ~ gamma  ), breaks= 30, ylim=c(0,1.1),
     cex.main =1.5, cex.lab=2)
curve(dnorm(x, beta_hat[2], sqrt(solve(t(X)%*%X)[2,2]*as.double(s2))), col="blue", lwd=3, add=TRUE)
abline(v = beta_hat[2], col ="red", lwd =3, lty=2)
text(beta_hat[2]+0.1, 1.05, expression(hat(gamma)), col="red", cex =1.3)
legend(4.6, 0.8, lty=c(1,2), lwd=c(3,3), col=c("blue", "red"), c("Posterior", "MLE est."))
dev.off()


# Predictions

n_pred <- 9
x_new <- sample(5:max(x), n_pred, replace=FALSE)
X_new <- matrix(c( rep(1,n_pred), x_new), n_pred,2, byrow = FALSE)
X_pred <- rbind(X, X_new)

y_pred <- matrix(NA, nsims, n_pred)


for (n in 1:nsims){
  y_pred[n,] <- rnorm(n_pred, X_new%*%beta[n,], sqrt(sigma2_sim[n]) )
}
  

y_pred_est <- apply(y_pred,2,mean)
y_pred_50 <- apply(y_pred,2,median)
y_pred_025 <- apply(y_pred,2, function(x) quantile(x, c(0.025)))
y_pred_975 <- apply(y_pred,2, function(x) quantile(x, c(0.975)))

pdf(file="cars_prev.pdf", width =7, height =6)
par(mfrow=c(1,1))
plot(x, y, xlab = "Speed", ylab="Distance", main ="cars data", xlim=c(1,25), ylim =c(-30, 120))
points(x_new, y_pred_50, col="red", pch = 16)
     for (i in 1:n_pred){
       segments(x_new[i], y_pred_975[i], x_new[i], y_pred_025[i], col="red", lty=2, lwd=2)
     }
legend(3, 110, col="red", lty=2, lwd=2, "Credible 95% intervals")
dev.off()


library(metRology)

pdf(file="hist_prev.pdf", width=12, height =7)
par(mfrow=c(3,3), mar=c(6,4,2,1))
for (i in 1:n_pred){
  sub<-substitute(a, list(a=i))
  hist(y_pred[,i], probability=TRUE, 
       main =bquote("Posterior predictive for" ~ tilde(y)[.(sub)]  ), 
       breaks= 30, 
       cex.main =1.5, xlab =expression(tilde(y)), cex.lab = 1.5,
       ylim=c(0,0.035))
  
  mean <- (X_new%*%beta_hat)[sub]
  sd <-   sqrt(as.double(s2)*(diag(n_pred) + X_new%*%solve(t(X)%*%X)%*%t(X_new))[sub,sub])
  
  curve(dt.scaled(x, n-p, mean, sd), add =TRUE, col="blue", lwd =3)
  legend(min(y_pred[,i]), 0.035, lwd=3, col="blue", expression(t[n-p]), cex=1.4)
}
dev.off()


# residuals


stand_res <- (y-c(X%*%beta_hat))/sqrt(s2)
pdf(file="cars_res.pdf", width=8.5, height =5.5)
par(mfrow=c(1,1), mar=c(5,5,2,1), xaxt="n", yaxt="n")
plot(y, stand_res, xlab="Distance", ylab="Standardized residual", cex.lab=1.8)
abline(h=1.96, col="blue", lty=2, lwd=2)
abline(h=-1.96, col="blue", lty=2, lwd=2)
par(xaxt="s", yaxt="s")
axis(1, cex.axis=1.6)
axis(2, cex.axis=1.6)
dev.off()

# ppcheck and residuals

n<- dim(cars)[1]
y_rep <- matrix(NA, nsims, n)
for (i in 1:nsims){
  y_rep[i,]<- rnorm(n, X%*%beta[i,], sqrt(sigma2_sim[i]))
}

res_matrix <- matrix(NA, nsims, n)
stand_matrix <- matrix(NA, nsims, n)
for (i in 1:nsims){
  lin.mod <- lm(y_rep[i,] ~ x)
  res_matrix[i,] <- lin.mod$residual
  stand_matrix[i,] <- (y_rep[i,]- as.vector(X%*%as.vector(lin.mod$coefficients[1:2])))/sqrt(s2)
}

stand_res_quant <- apply(stand_matrix, 2, function(x) quantile(x, c(0.025, 0.5, 0.975)))

pdf(file="cars_res_bay.pdf", width=8.5, height =5.5)
par(mfrow=c(1,1), mar=c(5,5,2,1), xaxt="n", yaxt="n")
plot(y, stand_res_quant[2,], xlab="Distance", ylab="Standardized residual", cex.lab=1.8,
     ylim=c(-3,3))
for (i in 1:n){
  segments(y[i], stand_res_quant[1,i], y[i], stand_res_quant[3,i],
           lty=2, col="red", lwd=2)
}

abline(h=1.96, col="blue", lty=2, lwd=2)
abline(h=-1.96, col="blue", lty=2, lwd=2)
par(xaxt="s", yaxt="s")
axis(1, cex.axis=1.6)
axis(2, cex.axis=1.6)
dev.off()

# proportion of outliers


prop_outliers <- apply(stand_matrix, 1, function(x) sum(x>1.96|x< -1.96))

outliers <- quantile(prop_outliers/nsims, c(0.025, 0.5, 0.975))
library(xtable)
xtable(outliers)



###############################
## INFORMATIVE PRIOR ANALYSIS
###############################


a <- 1
b <- 1
mu_beta <- rep(0, 2)
V_beta <- matrix(c(1,0,0,1),2,2, byrow = TRUE)
mu_star <- solve(solve(V_beta)+t(X)%*%X)%*%(solve(V_beta)%*%mu_beta+t(X)%*%y)
V_star <- solve(solve(V_beta)+t(X)%*%X)
a_star <- a+n/2
b_star <- b+(1/2)*(t(mu_beta)%*%solve(V_beta)%*%mu_beta+
                     t(y)%*%(y)-t(mu_star)%*%solve(V_star)%*%mu_star)
sigma2_sim_inf <- c()

sigma2_sim_inf <- rinvgamma(nsims, a_star, b_star)
beta_inf <- matrix(NA, nsims, 2)

for (i in 1:nsims){
  beta_inf[i,] <- rmvnorm(1, mu_star, V_star*sigma2_sim_inf[i])
}


pdf(file="marg_post_inf.pdf", width =7, height =6)
hist(sigma2_sim_inf, probability = TRUE, 
     xlab =expression(sigma^2), main =bquote("Marginal posterior for" ~ sigma^2  ), 
     breaks= 30, ylim=c(0,0.012), cex.main =1.5)
curve(dinvgamma(x, a_star, b_star ), col="blue", lwd=3, add=TRUE)
abline(v = s2, col ="red", lwd =3, lty=2)
text(s2+25, 0.009, expression(s^2), col="red", cex =1.3)
legend(400, 0.008, lty=c(1,2), lwd=c(3,3), col=c("blue", "red"), c("Posterior", "MLE est."))
dev.off()


pdf(file="cond_post_alpha_inf.pdf", width =7, height =6)
hist(beta_inf[,1], probability = TRUE, 
     xlab =expression(alpha), main =bquote("Conditional posterior for" ~ alpha  ), 
     breaks= 30, 
     ylim=c(0,0.1), 
     cex.main =1.5, cex.lab=2)
curve(dnorm(x, as.double(mu_star[1]), sqrt(V_star[1,1]*s2) ), col="blue", lwd=3, add=TRUE)
abline(v = beta_hat[1], col ="red", lwd =3, lty=2)
text(beta_hat[1]+1, 0.08, expression(hat(alpha)), col="red", cex=1.3)
legend(-12, 0.08, lty=c(1,2), lwd=c(3,3), col=c("blue", "red"), c("Posterior", "MLE est."))
dev.off()

pdf(file="cond_post_beta_inf.pdf", width =7, height =6)
hist(beta_inf[,2], probability = TRUE, 
     xlab =expression(gamma), main =bquote("Conditional posterior for" ~ gamma  ), 
     breaks= 30,xlim=c(2.5,6),  cex.main =1.5, cex.lab=2)
curve(dnorm(x, as.double(mu_star[2]), sqrt(V_star[2,2]*s2) ), col="blue", lwd=3, add=TRUE)
abline(v = beta_hat[2], col ="red", lwd =3, lty=2)
text(beta_hat[2]+0.2, 1, expression(hat(gamma)), col="red", cex=1.3)
legend(4.4, 1, lty=c(1,2), lwd=c(3,3), col=c("blue", "red"), c("Posterior", "MLE est."))
dev.off()

# extension: heteroschedasticity

  # cars data: perturbation

n <- 200
new_x <- seq(1,n)
new_y <- c()

for (i in 1:length(new_x)){
 new_y[i] <- rnorm(1, new_x[i], sd = new_x[i])
}

p <- 1
X <- matrix(NA, n,2)
X[,1] <- rep(1,n)
X[,2] <- new_x
beta_hat <- solve(t(X)%*%X)%*%t(X)%*%new_y
s2 <- ( n/(n-p)*1 /(n-p))*t(new_y-X%*%beta_hat)%*%(new_y-X%*%beta_hat )


pdf(file="cars_new.pdf", height =6, width =7)
plot(new_x, new_y, xlab = "x", ylab="y", main ="heteroschedastic data")
lm_est_new <- lm(new_y~new_x)
lines(new_x, lm_est_new$fitted.values , col="red", lwd =2)
dev.off()


