library(arm)
library(foreign)
library(R2WinBUGS)
library(rstanarm)
library(rstan)


#########################
# 1) LOGISTIC REGRESSION
########################


## 1992 polls

source("4.7_Fitting a series of regressions.R") # where data was cleaned; set the directory to be where this file is

yr <- 1992
ok <- nes.year==yr & presvote<3
vote <- presvote[ok] - 1
income <- data$income[ok]

# Estimation with glm
fit.1 <- glm (vote ~ income, family=binomial(link="logit"))
display(fit.1)

# Graph figure 5.1 (a) e 5.2 (b) dal Gelman e Hill (2007)
pdf(file="1992plots.pdf", width=9, height =6)
par(mfrow=c(1,2))
curve (invlogit(fit.1$coef[1] + fit.1$coef[2]*x), 1, 5, ylim=c(-.01,1.01),
       xlim=c(-2,8), xaxt="n", xaxs="i", mgp=c(2,.5,0),
       ylab="Pr (Republican vote)", xlab="Income", lwd=4,
       main ="Fitted logistic vs income")
curve (invlogit(fit.1$coef[1] + fit.1$coef[2]*x), -2, 8, lwd=.5, add=T)
axis (1, 1:5, mgp=c(2,.5,0))
mtext ("(poor)", 1, 1.5, at=1, adj=.5)
mtext ("(rich)", 1, 1.5, at=5, adj=.5)
points (jitter (income, .5), jitter (vote, .08), pch=20, cex=.1)

curve (invlogit(fit.1$coef[1] + fit.1$coef[2]*x), -15, 15, ylim=c(0,1.01),
       xlim=c(-15,15), xaxt="n", xaxs="i", mgp=c(2,.5,0),
       ylab="Pr (Republican vote)", xlab="Income", lwd=2,
       main = expression(logit^-1*(-1.4+0.44*x)))
segments(4.3, -0.2, 4.2, invlogit(fit.1$coef[1] + fit.1$coef[2]*4.2),
         lty=3)
segments(-15, invlogit(fit.1$coef[1] + fit.1$coef[2]*4.2), 4.2, 
         invlogit(fit.1$coef[1] + fit.1$coef[2]*4.2), lty =3)
text(4.3, invlogit(fit.1$coef[1] + fit.1$coef[2]*4.2)+0.05, "Halfway point: 1.44/0.33 = 4.2")

dev.off()

# stan model (rstan package)

data_1992 <- list(N = 1179,
                  vote = vote[is.na(vote)==FALSE],
                  income = income[is.na(income)==FALSE])
fit.rstan <- stan('1992polls.stan',
              data = data_1992, 
              iter =2000,
              chains =4)
print(fit.rstan)

# stan_glm noninformative
fit.2 <-  stan_glm (vote ~ income, 
                    family=binomial(link="logit"),
                    prior=normal(0, 100),
                    prior_intercept=normal(0,100))
print(fit.2)

# stan_glm weakly informative
fit.3 <-  stan_glm (vote ~ income, 
                    family=binomial(link="logit"),
                    prior=normal(0, 2.5),
                    prior_intercept=normal(0,10))
print(fit.3)

# stan_glm probit
fit.4 <- stan_glm (vote ~ income, 
                   family=binomial(link="probit"),
                   prior=normal(0, 2.5),
                   prior_intercept=normal(0,10))
print(fit.4)



## Separation

n <- 100
x1 <- rnorm (n)
x2 <- rbinom (n, 1, .5)
b0 <- 1
b1 <- 1.5
b2 <- 2
y <- rbinom (n, 1, invlogit(b0+b1*x1+b2*x2))

M1 <- glm (y ~ x1 + x2, family=binomial(link="logit"))
display (M1)

M2 <- stan_glm (y ~ x1 + x2, 
                family=binomial(link="logit"),
                prior=normal(0,100),
                prior_intercept = normal(0,100))
print(M2)  # just a test:  this should be very close to classical logit

M3 <- stan_glm (y ~ x1 + x2, 
                family=binomial(link="logit"))
print(M3)
# default normal prior with scale 2.5



logit <- function(x, b0, b1, b2, ind){
  1/(1+exp(-(b0+ b1*x + b2*ind)))
}



pdf(file="regr_log_1.pdf", width=10, height =7)
par(mfrow=c(1,1),xaxt="n", yaxt="n", mar=c(5,5,3,1))
curve(logit(x, b0 = M1$coefficients[1],
               b1 = M1$coefficients[2], 
               b2 = M1$coefficients[3],
               ind = 0), xlim=c(min(x1)-3, max(x1)+1), ylim =c(0,1), ylab ="logit", xlab ="x1",
               lwd = 3, main= "Dataset 1", cex.lab =2, cex.main =2)

curve(logit(x, b0 = M2$coefficients[1],
            b1 = M2$coefficients[2], 
            b2 = M2$coefficients[3],
            ind = 0), ylim =c(0,1), ylab ="logit", xlab ="x1", lwd = 3, col="red", add = TRUE)

curve(logit(x, b0 = M3$coefficients[1],
            b1 = M3$coefficients[2], 
            b2 = M3$coefficients[3],
            ind = 0), ylim =c(0,1), ylab ="logit", xlab ="x1", lwd = 3, col="blue", add = TRUE)
curve(logit(x, b0 = M1$coefficients[1],
            b1 = M1$coefficients[2], 
            b2 = M1$coefficients[3],
            ind = 1), xlim=c(min(x1)-2, max(x1)+1), ylim =c(0,1), ylab ="logit", xlab ="x1",
      lwd = 4, main= "x2=1", cex.lab =2, cex.main =2, lty =2, add =TRUE)

curve(logit(x, b0 = M2$coefficients[1],
            b1 = M2$coefficients[2], 
            b2 = M2$coefficients[3],
            ind = 1), ylim =c(0,1), ylab ="logit", xlab ="x1", lwd = 4, col="red", add = TRUE, lty =2)

curve(logit(x, b0 = M3$coefficients[1],
            b1 = M3$coefficients[2], 
            b2 = M3$coefficients[3],
            ind = 1), ylim =c(0,1), ylab ="logit", xlab ="x1", lwd = 4, col="blue", add = TRUE, lty =2)
par(xaxt="s", yaxt="s")
axis(1, cex.axis =1.7)
axis(2, cex.axis=1.7)
legend(0, 0.4, lty =1, lwd=3, col=c("black", "red", "blue"), c("classical", 
                                                               "Normal(0,100^2)", 
                                                               "Normal(0,2.5^2)"),
       cex =1.9)
text(-2, 0.87, "x2=1", cex= 1.9)
text(-0.7, 0.2, "x2=0", cex=1.9)
dev.off()





# Create separation:  set y=1 whenever x2=1
# Now it should blow up without the prior!

y <- ifelse (x2==1, 1, y)

M1 <- glm (y ~ x1 + x2, family=binomial(link="logit"))
display (M1)

M2 <- stan_glm (y ~ x1 + x2, 
                family=binomial(link="logit"),
                prior=normal(0,100), 
                prior_intercept=normal(0,100)) # Same as M1
print(M2)

M3 <- stan_glm (y ~ x1 + x2, 
                family=binomial(link="logit"))
print(M3)

M4 <- stan_glm (y ~ x1 + x2, 
                family=binomial(link="logit"),
                prior=normal(0,2.5), 
                prior_intercept=normal(0,10))  # Same as M3
print(M4)

M5 <- stan_glm (y ~ x1 + x2, 
                family=binomial(link="logit"),
                prior=cauchy(0,2.5),
                prior_intercept = cauchy(0,10))
print(M5)




logit <- function(x, b0, b1, b2, ind){
  1/(1+exp(-(b0+ b1*x + b2*ind)))
}



pdf(file="regr_log_2.pdf", width=10, height =7)
par(xaxt="n", yaxt="n", mar=c(5,5,3,1))
curve(logit(x, b0 = M1$coefficients[1],
            b1 = M1$coefficients[2], 
            b2 = M1$coefficients[3],
            ind = 0), xlim=c(min(x1)-9, max(x1)+3.2), ylim =c(0,1), ylab ="logit", xlab ="x1",
      lwd = 3, main= "Dataset 2", cex.lab =2, cex.main =2)

curve(logit(x, b0 = M2$coefficients[1],
            b1 = M2$coefficients[2], 
            b2 = M2$coefficients[3],
            ind = 0), ylim =c(0,1), ylab ="logit", xlab ="x1", lwd = 3, col="red", add = TRUE)

curve(logit(x, b0 = M3$coefficients[1],
            b1 = M3$coefficients[2], 
            b2 = M3$coefficients[3],
            ind = 0), ylim =c(0,1), ylab ="logit", xlab ="x1", lwd = 3, col="blue", add = TRUE)
curve(logit(x, b0 = M5$coefficients[1],
            b1 = M5$coefficients[2], 
            b2 = M5$coefficients[3],
            ind = 0), ylim =c(0,1), ylab ="logit", xlab ="x1", lwd = 3, col="darkgreen", add = TRUE)

curve(logit(x, b0 = M1$coefficients[1],
            b1 = M1$coefficients[2], 
            b2 = M1$coefficients[3],
            ind = 1), xlim=c(min(x1)-2, max(x1)+1), ylim =c(0,1), ylab ="logit", xlab ="x1",
      lwd = 4, main= "x2=1", cex.lab =2, cex.main =2, lty =2, add =TRUE)

curve(logit(x, b0 = M2$coefficients[1],
            b1 = M2$coefficients[2], 
            b2 = M2$coefficients[3],
            ind = 1), ylim =c(0,1), ylab ="logit", xlab ="x1", lwd = 4, col="red", add = TRUE, lty =2)

curve(logit(x, b0 = M3$coefficients[1],
            b1 = M3$coefficients[2], 
            b2 = M3$coefficients[3],
            ind = 1), ylim =c(0,1), ylab ="logit", xlab ="x1", lwd = 4, col="blue", add = TRUE, lty =2)
curve(logit(x, b0 = M5$coefficients[1],
            b1 = M5$coefficients[2], 
            b2 = M5$coefficients[3],
            ind = 1), ylim =c(0,1), ylab ="logit", xlab ="x1", lwd = 4, col="darkgreen", add = TRUE, lty=2)
par(xaxt="s", yaxt="s")
axis(1, cex.axis =1.7)
axis(2, cex.axis=1.7)
legend(-1, 0.4, lty =1, lwd=3, 
       col=c("black", "red", "blue", "darkgreen"), 
       c("classical", 
         "Normal(0,100^2)",
         "Normal(0,2.5^2)",
         "Cauchy(0,2.5)"), cex =1.9)
text(-7, 0.89, "x2=1", cex= 1.9)
text(0.5, 0.6, "x2=0", cex=1.9)
dev.off()



# bayesglm with gaussian family (bayes lm)
sigma <- 5
y2 <- rnorm (n, b0+b1*x1+b2*x2, sigma)
M7 <- bayesglm (y2 ~ x1 + x2, prior.scale=Inf, prior.df=Inf)
display (M7)


# bayesglm with categorical variables
z1 <- trunc(runif(n, 4, 9))
levels(factor(z1))
z2 <- trunc(runif(n, 15, 19))
levels(factor(z2))

## drop the base level (R default)
M8 <- bayesglm (y ~ x1 + factor(z1) + factor(z2),
                family=binomial(link="logit"), prior.scale=2.5, prior.df=Inf)
display (M8)

## keep all levels with the intercept, keep the variable order
M9 <- bayesglm (y ~ x1 + x1:x2 + factor(z1) + x2 + factor(z2),
                family=binomial(link="logit"),
                prior.mean=rep(0,12),
                prior.scale=rep(2.5,12),
                prior.df=rep(Inf,12),
                prior.mean.for.intercept=0,
                prior.scale.for.intercept=10,
                prior.df.for.intercept=1,
                drop.baseline=FALSE, keep.order=TRUE)
display (M9)

## keep all levels without the intercept
M10 <- bayesglm (y ~ x1 + factor(z1) + x1:x2 + factor(z2)-1,
                 family=binomial(link="logit"),
                 prior.mean=rep(0,11),
                 prior.scale=rep(2.5,11),
                 prior.df=rep(Inf,11),
                 drop.baseline=FALSE)
display (M10)



###########################
# 2) Probit regression
###########################

#comparison logit probit

logistic <- function(x){
  return(dlogis(x, 0,1))
}
probit <- function(x){
  return(dnorm(x, sd =1.6))
}

pdf(file="logit_probit.pdf", width=9.1, height =6)
par(mar=c(5,5,2,1))
curve(logistic(x), -6,6, col="red", lwd =2, cex.lab=2)
curve(probit(x), -6,6, col="blue", add =TRUE, lwd =2)
legend(2.2, 0.2, lwd=3, col =c("red", "blue"), c("Logit", "Probit (sd =1.6)"),
       cex =1.3)
dev.off()


#######################################
# 3) Poisson regression (Jonah example)
######################################

library(ggplot2)
library(rstanarm)
library(bayesplot)

pest_data <- readRDS('pest_data.RDS')
str(pest_data)


N_buildings <- length(unique(pest_data$building_id))
N_buildings

#? preliminary plots
ggplot(pest_data, aes(x = complaints)) + 
  geom_bar()
ggsave(file="hist_pest.pdf", width=8,height=6)

ggplot(pest_data, aes(x = traps, y = complaints, color = live_in_super == TRUE)) + 
  geom_jitter()

y <- pest_data$complaints
x <- pest_data$traps
M_pois <- stan_glm(y~x, family=poisson(link="log"))
print(M_pois)


library(bayesplot)

## posterior plots


mcmc_hist(as.matrix(M_pois))+
  xaxis_text(on =TRUE, size=rel(1.9))+
  yaxis_text(on =TRUE, size=rel(1.9))+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"))+
  legend_text(size=rel(12))
ggsave(file="hist_pest_post.pdf", width=8,height=6)

mcmc_scatter(as.matrix(M_pois), alpha = 0.2)+
  xaxis_text(on =TRUE, size=rel(1.9))+
  yaxis_text(on =TRUE, size=rel(1.9))+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=19,face="bold"))
ggsave(file="scatter_pest_post.pdf", width =8, height =6)


# extending the model (offset)

exposure <- log(pest_data$total_sq_foot/1e4)
M_pois_exposure <- stan_glm(y~x+offset(exposure), 
                            family=poisson(link="log"))
print(M_pois_exposure)


# overdispersion

pdf("poiss_negbin.pdf", width=6, height =5)
k<-0.5;
sub<-substitute(x, list(x=k))
plot(0:10, dnbinom(0:10, size=k,mu=2), xlab="Y",
     main = bquote(psi==.(sub)),
     ylab="f(y)",pch=21, bg=1, type="h", ylim =c(0,0.5), lwd=3, 
     cex.main =2, cex.lab=1.7)
points(0.2:10.2, dpois(0:10,2), col="red", pch=21, bg=1, type="h", lwd=3)
legend(4.3, 0.4, lty=1, lwd=3,col=c("black", "red"), c( "Neg. Binomial", "Poisson"),
       cex=1.4)
dev.off()

pdf("poiss_negbin2.pdf", width=6, height =5)
k<-10;
sub<-substitute(x, list(x=k))
plot(0:10, dnbinom(0:10, size=k,mu=2), xlab="Y",
     main = bquote(psi==.(sub)),
     ylab="f(y)",pch=21, bg=1, type="h", ylim =c(0,0.5), lwd=3, 
     cex.main =2, cex.lab=1.7)
points(0.2:10.2, dpois(0:10,2), col="red", pch=21, bg=1, type="h", lwd=3)
legend(4.3, 0.4, lty=1, lwd=3,col=c("black", "red"), c( "Neg. Binomial", "Poisson"),
       cex=1.4)
dev.off()




M_negbin <- stan_glm(y ~ x, 
                     family =neg_binomial_2(link="log"))
print(M_negbin)

M_negbin_exposure <- stan_glm(y ~ x, 
                              family =neg_binomial_2(link="log"),
                              offset=exposure)


mcmc_hist(as.matrix(M_negbin_exposure))+
  xaxis_text(on =TRUE, size=rel(1.9))+
  yaxis_text(on =TRUE, size=rel(1.9))+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"))+
  legend_text(size=rel(12))
ggsave(file="hist_pest_post_negbin.pdf", width=11.3,height=7)

mcmc_scatter(as.matrix(M_pois_rstan), alpha = 0.2)+
  xaxis_text(on =TRUE, size=rel(1.9))+
  yaxis_text(on =TRUE, size=rel(1.9))+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=19,face="bold"))
ggsave(file="scatter_pest_post_negbin.pdf", width =8, height =6)

print(M_negbin_exposure)


# residuals

nsims <- 10^2
y_rep <- matrix(NA, nsims, length(y))

lambda <- exp(as.double(M_negbin_exposure$coefficients[1])+
                as.double( M_negbin_exposure$coefficients[2]  )*x+
                exposure) 
inv_phi <- M_negbin_exposure$stan_summary[3,1]
inv_phi_vec <- c()
for (n in 1:nsims){
  y_rep[n,] <- rnbinom(length(y), mu=lambda, size= 1/inv_phi)
  M_temp <- stan_glm(y_rep[n,] ~ x, family =neg_binomial_2(link="log"),
                                   offset=exposure)
  inv_phi_vec[n] <- as.double(M_temp$stan_summary[3,1])
}

mean_inv_phi <- mean(inv_phi_vec)
y_rep <- colMeans(y_rep)
std_resid <- (y - y_rep) / sqrt(y_rep + y_rep^2*mean_inv_phi)
qplot(y_rep, std_resid) + hline_at(2) + hline_at(-2)+
  xaxis_text(on =TRUE, size=rel(1.9))+
  yaxis_text(on =TRUE, size=rel(1.9))+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=19,face="bold"))
ggsave(file="stand_res_pest.pdf", width=9, height=6)



