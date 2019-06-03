## ----global_options, include=FALSE---------------------------------------
knitr::opts_chunk$set(fig.align = 'center', warning=FALSE, message=FALSE, fig.asp=0.625, dev='png', global.par = TRUE, dev.args=list(pointsize=10), fig.path = 'figs/',  eval=TRUE)
library(MASS)

## ----setup, include=FALSE------------------------------------------------
library(knitr)
local({
    hook_plot = knit_hooks$get('plot')
    knit_hooks$set(plot = function(x, options) {
        paste0('\n\n----\n\n', hook_plot(x, options))
    })
})


## ----eval=TRUE,echo=TRUE, warning=FALSE, results='hide',message=FALSE----
library(rstan)
library(dplyr)
library(lubridate)
library(ggplot2)
library(bayesplot)

theme_set(bayesplot::theme_default())

# con questo stan durante la compilazione salva un file rds
# che contiene il modello compilato e ogni volta che lo 
# chiamiamo fa un check e in caso ricompila
rstan_options(auto_write=TRUE)

## ----eval=TRUE,echo=TRUE-------------------------------------------------
pest_data <- readRDS('lab4/pest_data.RDS')
str(pest_data)
summary(pest_data)

##number of buildings
N_buildings <- length(unique(pest_data$building_id))
N_buildings


## ----stan-data-----------------------------------------------------------
## arrange data into a list
stan_dat_simple <- list(
  N = nrow(pest_data), 
  complaints = pest_data$complaints,
  traps = pest_data$traps,
  log_sq_foot = log(pest_data$total_sq_foot/1e4),
  live_in_super = pest_data$live_in_super
)
str(stan_dat_simple)


## ---- cache=TRUE, results="hide", message=FALSE--------------------------
comp_model_NB <- stan_model('lab4/multiple_NB_regression.stan')


## ----runNB---------------------------------------------------------------
fitted_model_NB <- sampling(comp_model_NB, data = stan_dat_simple)
samps_NB <- rstan::extract(fitted_model_NB)


## ----ppc-group_means-----------------------------------------------------
y_rep <- samps_NB$y_rep

ppc_stat_grouped(
  y = stan_dat_simple$complaints, 
  yrep = y_rep, 
  group = pest_data$building_id, 
  stat = 'mean',
  binwidth = 0.2
)


## ----prep-data-----------------------------------------------------------
N_months <- length(unique(pest_data$date))

# Add some IDs for building and month
pest_data <- pest_data %>%
  mutate(
    building_fac = factor(building_id, levels = unique(building_id)),
    building_idx = as.integer(building_fac),
    ids = rep(1:N_months, N_buildings),
    mo_idx = lubridate::month(date)
  )

# Center and rescale the building specific data
building_data <- pest_data %>%
    select(
      building_idx,
      live_in_super,
      age_of_building,
      total_sq_foot,
      average_tenant_age,
      monthly_average_rent
    ) %>%
    unique() %>%
    arrange(building_idx) %>%
    select(-building_idx) %>%
    scale(scale=FALSE) %>%
    as.data.frame() %>%
    mutate( # scale by constants
      age_of_building = age_of_building / 10,
      total_sq_foot = total_sq_foot / 10000,
      average_tenant_age = average_tenant_age / 10,
      monthly_average_rent = monthly_average_rent / 1000
    ) %>%
    as.matrix()

str(building_data)

stan_dat_hier <-
  with(pest_data,
        list(complaints = complaints,
             traps = traps,
             N = length(traps),
             J = N_buildings,
             log_sq_foot = log(pest_data$total_sq_foot/1e4),
             building_data = building_data[,-3],
             mo_idx = as.integer(as.factor(date)),
             K = 4,
             building_idx = building_idx
             )
       )
str(stan_dat_hier)


## ----comp-NB-hier, cache=TRUE, results="hide", message=FALSE-------------
comp_model_NB_hier <- stan_model('lab5/hier_NB_regression.stan')

## ----run-NB-hier---------------------------------------------------------
fitted_model_NB_hier <- sampling(comp_model_NB_hier, data = stan_dat_hier,
  chains = 4, cores = 4, iter = 4000)


## ------------------------------------------------------------------------
samps_hier_NB <- rstan::extract(fitted_model_NB_hier)


## ----print-NB-hier-------------------------------------------------------
print(fitted_model_NB_hier, pars = c('sigma_mu','beta','alpha','phi','mu'))

## ------------------------------------------------------------------------
mcmc_trace(
  as.array(fitted_model_NB_hier,pars = 'sigma_mu'),
  np = nuts_params(fitted_model_NB_hier),
  window = c(500,1000)
)

## ------------------------------------------------------------------------
# assign to object so we can compare to another plot later
scatter_with_divs <- mcmc_scatter(
  as.array(fitted_model_NB_hier),
  pars = c("mu[4]", 'sigma_mu'),
  transform = list('sigma_mu' = "log"),
  np = nuts_params(fitted_model_NB_hier)
)
scatter_with_divs


## ------------------------------------------------------------------------
N_sims <- 1000
log_sigma <- rep(NA, N_sims)
theta <- rep(NA, N_sims)
for (j in 1:N_sims) {
  log_sigma[j] <- rnorm(1, mean = 0, sd = 1)
  theta[j] <- rnorm(1, mean = 0, sd = exp(log_sigma[j]))
}
draws <- cbind("mu" = theta, "log(sigma_mu)" = log_sigma)
mcmc_scatter(draws)


## ------------------------------------------------------------------------
parcoord_with_divs <- mcmc_parcoord(
  as.array(fitted_model_NB_hier, pars = c("sigma_mu", "mu")),
  np = nuts_params(fitted_model_NB_hier)
)
parcoord_with_divs


## ----comp-NB-hier-ncp, cache=TRUE----------------------------------------
comp_model_NB_hier_ncp <- stan_model('lab5/hier_NB_regression_ncp.stan')


## ----run-NB-hier-ncp-----------------------------------------------------
fitted_model_NB_hier_ncp <- sampling(comp_model_NB_hier_ncp, data = stan_dat_hier, chains = 4, cores = 4)


## ----n-eff-NB-hier-ncp-check---------------------------------------------
print(fitted_model_NB_hier_ncp, pars = c('sigma_mu','beta','alpha','phi','mu'))


## ------------------------------------------------------------------------
scatter_no_divs <- mcmc_scatter(
  as.array(fitted_model_NB_hier_ncp),
  pars = c("mu[4]", 'sigma_mu'),
  transform = list('sigma_mu' = "log"),
  np = nuts_params(fitted_model_NB_hier_ncp)
)
bayesplot_grid(scatter_with_divs, scatter_no_divs,
               grid_args = list(ncol = 2), ylim = c(-11, 1))


## ------------------------------------------------------------------------
parcoord_no_divs <- mcmc_parcoord(
  as.array(fitted_model_NB_hier_ncp, pars = c("sigma_mu", "mu")),
  np = nuts_params(fitted_model_NB_hier_ncp)
)
bayesplot_grid(parcoord_with_divs, parcoord_no_divs,
               ylim = c(-3, 3))


## ----samps-full-hier-----------------------------------------------------
samps_NB_hier_ncp <- rstan::extract(fitted_model_NB_hier_ncp, pars = c('y_rep','inv_phi'))


## ----ppc-full-hier-------------------------------------------------------
y_rep <- as.matrix(fitted_model_NB_hier_ncp, pars = "y_rep")
ppc_dens_overlay(stan_dat_hier$complaints, y_rep[1:200,])

# qua ci sono 200 linee azzurre, una per ogni simulazione
# quella blu invece corrisponde ai dati

## ----ppc-group_means-hier------------------------------------------------
ppc_stat_grouped(
  y = stan_dat_hier$complaints,
  yrep = y_rep,
  group = pest_data$building_id,
  stat = 'mean',
  binwidth = 0.5
)


## ------------------------------------------------------------------------
ppc_stat_grouped(
  y = stan_dat_hier$complaints,
  yrep = y_rep,
  group = pest_data$building_id,
  stat = 'sd',
  binwidth = 0.5
)


## ------------------------------------------------------------------------
ppc_intervals(
  y = stan_dat_hier$complaints,
  yrep = y_rep,
  x = stan_dat_hier$traps
) +
  labs(x = "Number of traps", y = "Number of complaints")

## ------------------------------------------------------------------------
mean_y_rep <- colMeans(y_rep)
mean_inv_phi <- mean(as.matrix(fitted_model_NB_hier_ncp, pars = "inv_phi"))
std_resid <- (stan_dat_hier$complaints - mean_y_rep) / sqrt(mean_y_rep + mean_y_rep^2*mean_inv_phi)
qplot(mean_y_rep, std_resid) + hline_at(2) + hline_at(-2)


#################### 24 maggio

## modello gerarchico con pendenza variabile
stan_dat_hier <- readRDS('lab5/pest_data_longer_stan_dat.RDS')

## ----comp-NB-hier-slopes, cache=TRUE, results="hide", message=FALSE------
comp_model_NB_hier_slopes <- stan_model('lab5/hier_NB_regression_ncp_slopes_mod.stan')


## ----run-NB-hier-slopes--------------------------------------------------
fitted_model_NB_hier_slopes <-
  sampling(
    comp_model_NB_hier_slopes,
    data = stan_dat_hier,
    chains = 4, cores = 4,
    control = list(adapt_delta = 0.95)
  )


## ------------------------------------------------------------------------
mcmc_hist(
  as.matrix(fitted_model_NB_hier_slopes, pars = "sigma_kappa"),
  binwidth = 0.005
)


## ------------------------------------------------------------------------
print(fitted_model_NB_hier_slopes, pars = c('kappa','beta','alpha','phi','sigma_mu','sigma_kappa','mu'))


## ------------------------------------------------------------------------
mcmc_hist(
  as.matrix(fitted_model_NB_hier_slopes, pars = "beta"),
  binwidth = 0.005
)


## ----ppc-full-hier-slopes------------------------------------------------
y_rep <- as.matrix(fitted_model_NB_hier_slopes, pars = "y_rep")
ppc_dens_overlay(
  y = stan_dat_hier$complaints,
  yrep = y_rep[1:200,]
)


# potrebbe esserci una differenza tra un mese e l'altro


## ----ppc-group_max-hier-slopes-mean-by-mo--------------------------------
select_vec <- which(stan_dat_hier$mo_idx %in% 1:12)
ppc_stat_grouped(
  y = stan_dat_hier$complaints[select_vec],
  yrep = y_rep[,select_vec],
  group = stan_dat_hier$mo_idx[select_vec],
  stat = 'mean'
) + xlim(0, 11)


## ----comp-NB-hier-mos, cache=TRUE, results="hide", message=FALSE---------
comp_model_NB_hier_mos <- stan_model('lab5/hier_NB_regression_ncp_slopes_mod_mos.stan')


## ----run-NB-hier-slopes-mos----------------------------------------------
fitted_model_NB_hier_mos <- sampling(comp_model_NB_hier_mos, data = stan_dat_hier, chains = 4, cores = 4, control = list(adapt_delta = 0.9))


## ----ppc-full-hier-mos---------------------------------------------------
y_rep <- as.matrix(fitted_model_NB_hier_mos, pars = "y_rep")
ppc_dens_overlay(
  y = stan_dat_hier$complaints,
  yrep = y_rep[1:200,]
)


## ------------------------------------------------------------------------
select_vec <- which(stan_dat_hier$mo_idx %in% 1:12)
ppc_stat_grouped(
  y = stan_dat_hier$complaints[select_vec],
  yrep = y_rep[,select_vec],
  group = stan_dat_hier$mo_idx[select_vec],
  stat = 'mean'
)


## ------------------------------------------------------------------------
# 1) compare draws from prior and draws from posterior
rho_draws <- cbind(
  2 * rbeta(4000, 10, 5) - 1, # draw from prior
  as.matrix(fitted_model_NB_hier_mos, pars = "rho")
)
colnames(rho_draws) <- c("prior", "posterior")
mcmc_hist(rho_draws, freq = FALSE, binwidth = 0.025,
          facet_args = list(nrow = 2)) + xlim(-1, 1)


# 2) overlay prior density curve on posterior draws
gen_rho_prior <- function(x) {
  alpha <- 10; beta <- 5
  a <- -1; c <- 1
  lp <- (alpha - 1) * log(x - a) +
    (beta - 1) * log(c - x) -
    (alpha + beta - 1) * log(c - a) -
    lbeta(alpha, beta)
  return(exp(lp))
}
mcmc_hist(as.matrix(fitted_model_NB_hier_mos, pars = "rho"),
          freq = FALSE, binwidth = 0.01) +
  overlay_function(fun = gen_rho_prior) +
  xlim(-1,1)

## ------------------------------------------------------------------------
print(fitted_model_NB_hier_mos, pars = c('rho','sigma_mu','sigma_kappa','gamma'))


## ------------------------------------------------------------------------
ppc_intervals(
  y = stan_dat_hier$complaints,
  yrep = y_rep,
  x = stan_dat_hier$traps
) +
  labs(x = "Number of traps", y = "Number of complaints")

## LOOCV

## ------------------------------------------------------------------------
comp_model_NB_hier_slopes_llik <- stan_model('lab5/hier_NB_regression_ncp_slopes_mod_mos_loglik.stan')

fitted_model_NB_hier_slopes_llik <-
  sampling(
    comp_model_NB_hier_slopes_llik,
    data = stan_dat_hier,
    chains = 4, cores = 4,
    control = list(adapt_delta = 0.95)
  )


## ------------------------------------------------------------------------
comp_model_NB_hier_mos_llik <- stan_model('hier_NB_regression_ncp_slopes_mod_mos_loglik.stan')

fitted_model_NB_hier_mos_llik <-
  sampling(
    comp_model_NB_hier_mos_llik,
    data = stan_dat_hier,
    chains = 4, cores = 4,
    control = list(adapt_delta = 0.95)
  )



## ----run-NB-hier-slopes-mos-loo------------------------------------------
library(loo)
log_lik_slopes <- extract_log_lik(fitted_model_NB_hier_slopes_llik)
loo_slopes <- loo(log_lik_slopes)
loo_slopes
log_lik_mos <- extract_log_lik(fitted_model_NB_hier_mos_llik)
loo_mos <- loo(log_lik_mos)
loo_mos

compare(loo_slopes, loo_mos)


