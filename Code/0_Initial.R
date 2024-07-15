
n.param.gp <- n.sites
n.param.beta <- 1

load(file=paste0('./init/init_vals_', species.no, '.RData'))

alpha_n <- (0.1)
lambda_n <- (0.7)

###### Initializing matrices
x_pp <- x_pp_2 <- x_mu <- x_sigma <-  x_bin <- matrix(nrow= B, ncol=(n.param.gp ))

hyper_gp_latent_pp <- matrix(nrow= B, ncol=(2 ))
hyper_gp_latent_pp2 <- matrix(nrow= B, ncol=(2 ))
hyper_gp_latent_mu <- matrix(nrow= B, ncol=(2 ))
hyper_gp_latent_sigma <- matrix(nrow= B, ncol=(2 ))
hyper_gp_latent_re <- matrix(nrow= B, ncol=(2 ))
hyper_gp_latent_bin <- matrix(nrow= B, ncol=(2 ))

beta_pp <- beta_pp_2 <- beta_bin <-
  beta_mu <- beta_sigma <- beta_xi <- matrix(nrow= B, ncol=(n.param.beta ))

fe_bin <- fe_mu <-matrix(nrow= B, ncol=(1 ))
fe_pp2 <- NULL

beta_share <- matrix(nrow= B, ncol=5)

# acceptance probability tracking
ind_intercept <- matrix(NA, nrow= B, ncol=10)
ind_fe <- matrix(NA, nrow= B, ncol=10)
ind_share <- matrix(NA, nrow= B, ncol=5)
ind_hyper_pp <- ind_hyper_pp2 <- ind_hyper_mu <-
  ind_hyper_sigma <- ind_hyper_re <- ind_hyper_bin <- matrix(NA, nrow= B, ncol=2)
ind_latent_gp <- matrix(NA, nrow= B, ncol=5)
ind <- matrix(NA, nrow= 2, ncol=2)

data_test <- cbind(apply(dat.mat, 2, function(x) mean(x, na.rm=TRUE)),
                   coords[,1:2])
alpha_n_time <- 0.1
lambda_n_time <- 0.7

#### SET initial values

##### for hyperpriors
hyper_gp_latent_pp[1,] <- hyper_gp_latent_pp2[1,] <- hyper_gp_latent_mu[1,] <-
  hyper_gp_latent_sigma[1,] <- hyper_gp_latent_bin[1,] <- c(alpha_n, lambda_n)
hyper_gp_latent_re[1,] <- c(alpha_n_time, lambda_n_time)

mu <- 0
sigma <- 100

#### Starting values for latent

set.seed(2023)
x_pp[1,] <- init_vals$x_pp
x_pp_2[1,] <- rnorm(n.sites, 0, 0.1)
x_mu[1,] <- rnorm(n.sites, 0, 0.1)
#x_sigma[1,] <- rnorm(n.sites, 0, 0.1)
x_sigma[1,] <- rnorm(n.sites, 0, 0.1)
x_bin[1,] <- 0 ### this is anyway fixed to 0 in MCMC

#sum to zero constraint
x_pp[1,]<- x_pp[1,]-mean(x_pp[1,])
x_pp_2[1,]<- x_pp_2[1,]-mean(x_pp_2[1,])
x_mu[1,]<- x_mu[1,]-mean(x_mu[1,])
x_sigma[1,]<- x_sigma[1,]-mean(x_sigma[1,])

#### Starting values for sharing/fixed effects

beta_xi[1,] <- -0.2
beta_sigma[1,] <- log(sd(dat.mat, na.rm=TRUE))
beta_mu[1,] <- log(mean(dat.mat, na.rm=TRUE))
beta_pp[1,] <- init_vals$beta_pp
beta_pp_2[1,] <- init_vals$beta_pp_2 #
beta_bin[1,] <- init_vals$beta_bin

1- exp( -exp( beta_bin[1,]) )


#size of spatial pixel
size_A <- 20*20


# sharing
beta_share[1,1] <- 0 # pref sharing
beta_share[1,2] <- 2 # intercept effort
beta_share[1,4] <- 1 # this is fixed in the MCMC
beta_share[1,3] <- 0 # GEV-niche
beta_share[1,5] <- 0 # activity sharing

fe_mu[1,] <- 0
fe_bin[1,] <- -1 # should start from negative values

fe_pp <- matrix(nrow= B, ncol=(n.year ))
fe_pp[1,] <- init_vals$fe_pp
hyper_gp_latent_re[1,] <- c(alpha_n, lambda_n)


## this should be 0, no infinite loglik to start with...
print(sum(test_likelihood(1)$loglik))

NNarray <- find_ordered_nn(as.matrix(coords),neigh)
NNarray_time <- find_ordered_nn(as.matrix(1:21, ncol=1),neigh)

log_posteriors <- log_posteriors_2 <- log_posteriors_3 <- log_posteriors_4 <- numeric()

