

#######################################################
#### 0. Cluster inputs ################################
#######################################################


input1 <- 1 # for species number
input2 <- 1 # for naming purpose only

species.no <- as.numeric(input1)
mcmc.no <- as.numeric(input2)

####################################################
#### 1. Source base functions, data, libraries #####
####################################################
source(paste0("./functions_all_new.R") )
source(paste0("./0_Data_prep.R") )
weight_route <- (t(df_routes[,-c(1,2)])/route_length)

library(GpGp)
# B is number of MCMC samples
B <- 20001

sum <- function(x){
  base::sum(x, na.rm=TRUE)
}

### For PC priors/ GP
sigma.sq_long=1
sigma.sq_medium=1
sigma.sq_short=1
# scale (range)
phi_long=100
phi_medium=50
phi_short=10
neigh <- 5

# For Vecchia approximation:
ord <- order_maxmin(coords) # calculate an ordering
locsord <- coords[ord,] # reorder locations
ord_time <- order_maxmin(matrix(1:21, ncol=1)) # calculate an ordering
locsord_time <- matrix(1:21, ncol=1)[ord_time,] # reorder locations

## g_func


g_func <- function(theta_1, theta_3,x,theta_2){
  exp(theta_2)*( 1/(1+exp(-(theta_3+x)*(theta_1) )) )
}

g_d_func <- function(theta_1, theta_3,x,theta_2){
  exp(theta_2) * (exp(-(theta_3 + x) * (theta_1)) * (theta_1)/(1 +
                                                                 exp(-(theta_3 + x) * (theta_1)))^2)
}

g_d_theta2_func <- function(theta_1, theta_3,x,theta_2){
  exp(theta_2) * (1/(1 + exp(-(theta_3 + x) * (theta_1))))
}


## PRIORS..
#checking out priors
shape_alpha_n <- 1
rate_alpha_n <- 2

shape_lambda_n <- 3
rate_lambda_n <- 1


#############################
##### 2. Initialising #######
#############################

source('./0_Initial.R')
##############
## stepsizes
##############
hess=1
# MALA: effort, niche, mu, bin (unused), sigma
delta=c(0.001, 0.01, 0.005, 0.01, 0.005)
# hyperparam1,2,
stepsize <- c(0.15, 0.3)
# ###                     (mu,   sigma,    xi,  pp1, bin ,pp2)
stepsize_intercept <- c(0.001, 0.001, 0.002, 0.001 , 0.001,  0.01   )
# ###             fe_pp, fe_bin, fe_mu
stepsize_fe <- c(0.001, 0.1, 0.001)
stepsize_share <- c(0.05, 0.05, 0.005, NA, 0.05)

## Adaptive step sizes? Set to true in beginning of chain

## Adaptive step sizes? Set to true in beginning of chain
# change this after burn-in!
if (mcmc.no>3){
  adaptive=FALSE
} else {
  adaptive=TRUE
}

## load from previous 20k run if mcmc.no>2

if (mcmc.no>=2){

  load(paste0('./Output/', mcmc.no-1, '/mcmc_data_', mcmc.no-1, '_', species.no, '.RData') )

  hyper_gp_latent_pp[1,] <- function_list_store$hyper_gp_latent_pp[sim,]
  hyper_gp_latent_pp2[1,] <- function_list_store$hyper_gp_latent_pp2[sim,]
  hyper_gp_latent_mu[1,] <- function_list_store$hyper_gp_latent_mu[sim,]
  hyper_gp_latent_sigma[1,] <- function_list_store$hyper_gp_latent_sigma[sim,]
  hyper_gp_latent_re[1,] <- function_list_store$hyper_gp_latent_re[sim,]
  hyper_gp_latent_bin[1,] <- function_list_store$hyper_gp_latent_bin[sim,]

  # latent fields
  x_mu[1,] <- function_list_store$x_mu[sim,]
  x_pp[1,] <- function_list_store$x_pp[sim,]
  x_pp_2[1,] <- function_list_store$x_pp_2[sim,]
  x_sigma[1,] <- function_list_store$x_sigma[sim,]
  x_bin[1,] <- function_list_store$x_bin[sim,]

  # intercepts
  beta_pp[1,] <- function_list_store$beta_pp[sim,]
  beta_pp_2[1,] <- function_list_store$beta_pp_2[sim,]
  beta_mu[1,] <- function_list_store$beta_mu[sim,]
  beta_sigma[1,] <- function_list_store$beta_sigma[sim,]
  beta_xi[1,] <- function_list_store$beta_xi[sim,]
  beta_bin[1,] <- function_list_store$beta_bin[sim,]

  #sharing params
  beta_share[1,] <- function_list_store$beta_share[sim,]
  #fixed effects
  fe_mu[1,] <- function_list_store$fe_mu[sim,]
  fe_pp[1,] <- function_list_store$fe_pp[sim,]
  fe_pp2[1,] <- function_list_store$fe_pp2[sim,]
  fe_bin[1,] <- function_list_store$fe_bin[sim,]

  # list of acceptances
  ind <- function_list_store$ind
  # one of the log posteriors to check if it is increasing
  log_posteriors[1] <- function_list_store$log_posteriors[sim]
  log_posteriors_2[1] <- function_list_store$log_posteriors_2[sim]
  log_posteriors_3[1] <- function_list_store$log_posteriors_3[sim]
  log_posteriors_4[1] <- function_list_store$log_posteriors_4[sim]

}


start_time <- Sys.time()

for (sim in 1:(B-1)){

  if (adaptive==TRUE){
    if (sim>2000 & sim%%500==0){

      which.low <- which(apply(ind_latent_gp[(sim-450):(sim-1),], 2, mean)<0.3)
      which.high <- which(apply(ind_latent_gp[(sim-450):(sim-1),], 2, mean)>0.95)
      delta[which.low] <- 0.8*delta[which.low]
      delta[which.high] <- 1.2*delta[which.high]

      if (mean(ind_fe[(sim-450):(sim-1),1])<0.3){
        stepsize_fe[1] <- 0.8*stepsize_fe[1]
      } else if (mean(ind_fe[(sim-450):(sim-1),1])>0.95){
        stepsize_fe[1] <- 1.2*stepsize_fe[1]
      }

      which.low <- which(apply(ind_intercept[(sim-450):(sim-1),], 2, mean)<0.15)
      which.high <- which(apply(ind_intercept[(sim-450):(sim-1),], 2, mean)>0.65)
      stepsize_intercept[which.low] <- 0.8*stepsize_intercept[which.low]
      stepsize_intercept[which.high] <- 1.2*stepsize_intercept[which.high]

      which.low <- which(apply(ind_fe[(sim-450):(sim-1),-1], 2, mean)<0.15)
      which.high <- which(apply(ind_fe[(sim-450):(sim-1),-1], 2, mean)>0.65)
      stepsize_fe[c(2:length(stepsize_fe))[which.low] ] <- 0.8*stepsize_fe[c(2:length(stepsize_fe))[which.low] ]
      stepsize_fe[c(2:length(stepsize_fe))[which.high] ] <- 1.2*stepsize_fe[c(2:length(stepsize_fe))[which.high] ]

      which.low <- which(apply(ind_share[(sim-450):(sim-1),], 2, mean)<0.15)
      which.high <- which(apply(ind_share[(sim-450):(sim-1),], 2, mean)>0.65)
      stepsize_share[which.low] <- 0.8*stepsize_share[which.low]
      stepsize_intercept[which.high] <- 1.2*stepsize_share[which.high]

    }
  }



  function_list <- MCMC_function_GP(sim, hess, delta, stepsize, log_posteriors[c(sim,sim+1)],
                                    hyper_gp_latent_pp[c(sim,sim+1),], hyper_gp_latent_pp2[c(sim,sim+1),],
                                    hyper_gp_latent_mu[c(sim,sim+1),], hyper_gp_latent_sigma[c(sim,sim+1),],
                                    x_pp[c(sim,sim+1),], x_pp_2[c(sim,sim+1),], x_mu[c(sim,sim+1),], x_sigma[c(sim,sim+1),],
                                    matrix(beta_pp[c(sim,sim+1),]), matrix(beta_mu[c(sim,sim+1),]),
                                    matrix(beta_sigma[c(sim,sim+1),]), matrix(beta_xi[c(sim,sim+1),]),
                                    beta_share[c(sim,sim+1),], matrix(beta_pp_2[c(sim,sim+1),]),
                                    matrix(beta_bin[c(sim,sim+1),]),
                                    matrix(fe_mu[c(sim,sim+1),]),fe_pp[c(sim,sim+1),],
                                    fe_pp2[c(sim,sim+1),], matrix(fe_bin[c(sim,sim+1),]),
                                    log_posteriors_2[c(sim,sim+1)], log_posteriors_3[c(sim,sim+1)], log_posteriors_4[c(sim,sim+1)],
                                    hyper_gp_latent_bin[c(sim,sim+1),], x_bin[c(sim,sim+1),],
                                    hyper_gp_latent_re[c(sim,sim+1),])




  hyper_gp_latent_pp[c(sim,sim+1),] <- function_list$hyper_gp_latent_pp
  hyper_gp_latent_pp2[c(sim,sim+1),] <- function_list$hyper_gp_latent_pp2
  hyper_gp_latent_mu[c(sim,sim+1),] <- function_list$hyper_gp_latent_mu
  hyper_gp_latent_sigma[c(sim,sim+1),] <- function_list$hyper_gp_latent_sigma
  hyper_gp_latent_re[c(sim,sim+1),] <- function_list$hyper_gp_latent_re
  hyper_gp_latent_bin[c(sim,sim+1),] <- function_list$hyper_gp_latent_bin

  # latent fields
  x_mu[c(sim,sim+1),] <- function_list$x_mu
  x_pp[c(sim,sim+1),] <- function_list$x_pp
  x_pp_2[c(sim,sim+1),] <- function_list$x_pp_2
  x_sigma[c(sim,sim+1),] <- function_list$x_sigma
  x_bin[c(sim,sim+1),] <- function_list$x_bin

  # intercepts
  beta_pp[c(sim,sim+1),] <- function_list$beta_pp
  beta_pp_2[c(sim,sim+1),] <- function_list$beta_pp_2
  beta_mu[c(sim,sim+1),] <- function_list$beta_mu
  beta_sigma[c(sim,sim+1),] <- function_list$beta_sigma
  beta_xi[c(sim,sim+1),] <- function_list$beta_xi
  beta_bin[c(sim,sim+1),] <- function_list$beta_bin

  #sharing params
  beta_share[c(sim,sim+1),] <- function_list$beta_share
  #fixed effects
  fe_mu[c(sim,sim+1),] <- function_list$fe_mu
  fe_pp[c(sim,sim+1),] <- function_list$fe_pp
  fe_pp2[c(sim,sim+1),] <- function_list$fe_pp2
  fe_bin[c(sim,sim+1),] <- function_list$fe_bin

  # list of acceptances
  ind <- function_list$ind
  # one of the log posteriors to check if it is increasing
  log_posteriors[c(sim,sim+1)] <- function_list$log_posteriors
  log_posteriors_2[c(sim,sim+1)] <- function_list$log_posteriors_2
  log_posteriors_3[c(sim,sim+1)] <- function_list$log_posteriors_3
  log_posteriors_4[c(sim,sim+1)] <- function_list$log_posteriors_4

  print(sim)

  # save output first 100, 500, and 1000, then every 2000 iterations

  if (sim%%2000==0 | sim==100 | sim==500 | sim==1000){
    #### store
    function_list_store <- list()

    #hyperparameters
    function_list_store$hyper_gp_latent_pp <- hyper_gp_latent_pp
    function_list_store$hyper_gp_latent_pp2 <- hyper_gp_latent_pp2
    function_list_store$hyper_gp_latent_mu <- hyper_gp_latent_mu
    function_list_store$hyper_gp_latent_sigma <- hyper_gp_latent_sigma
    function_list_store$hyper_gp_latent_re <- hyper_gp_latent_re
    function_list_store$hyper_gp_latent_bin <- hyper_gp_latent_bin

    # latent fields
    function_list_store$x_mu <- x_mu
    function_list_store$x_pp <- x_pp
    function_list_store$x_pp_2 <- x_pp_2
    function_list_store$x_sigma <- x_sigma
    function_list_store$x_bin <- x_bin

    # intercepts
    function_list_store$beta_pp <- beta_pp
    function_list_store$beta_mu <- beta_mu
    function_list_store$beta_sigma <- beta_sigma
    function_list_store$beta_xi <- beta_xi
    function_list_store$beta_pp_2 <- beta_pp_2
    function_list_store$beta_bin <- beta_bin


    #sharing params
    function_list_store$beta_share <- beta_share

    #fixed effects
    function_list_store$fe_mu <- fe_mu
    function_list_store$fe_pp <- fe_pp
    function_list_store$fe_pp2 <- fe_pp2
    function_list_store$fe_bin <- fe_bin

    # list of acceptances
    function_list_store$ind <- ind

    # one of the log posteriors to check if it is increasing
    function_list_store$log_posteriors <- log_posteriors
    function_list_store$log_posteriors_2 <- log_posteriors_2
    function_list_store$log_posteriors_3 <- log_posteriors_3
    function_list_store$log_posteriors_4 <- log_posteriors_4


    # save(function_list_store, sim, ind_latent_gp, ind_intercept, ind_fe, ind_share, coords, delta,
    #      stepsize, stepsize_intercept, stepsize_fe, stepsize_share, ind_hyper_pp, ind_hyper_pp2, ind_hyper_mu,
    #        ind_hyper_sigma, ind_hyper_re, ind_hyper_bin,
    #      file=paste0('./Output/', mcmc.no, '/mcmc_data_', mcmc.no, '_', species.no, '.RData'))
  }
}

end_time <- Sys.time()

end_time - start_time

# ~6.28 sec/iteration !


########################
#### Some basic plots #######
########################


burn_in <- 1

plot(beta_share[burn_in:sim,1], type='l')
plot(beta_share[burn_in:sim,2], type='l')
plot(beta_share[burn_in:sim,3], type='l')
abline(h=0)
plot(beta_share[burn_in:sim,5], type='l')
abline(h=0)

x <- x_pp_2
plot(x[burn_in:sim,120], type='l')
plot(x[burn_in:sim,50], type='l')

apply(ind_latent_gp[burn_in:(sim-1),], 2, mean)
apply(ind_intercept[burn_in:(sim-1),], 2, mean)
apply(ind_fe[burn_in:(sim-1),], 2, mean)
apply(ind_share[burn_in:(sim-1),], 2, mean)

plot(beta_bin[burn_in:sim,], type='l')
plot(beta_mu[burn_in:sim,], type='l')
plot(beta_xi[burn_in:sim,], type='l')
plot(beta_sigma[burn_in:sim,], type='l')
plot(beta_pp_2[burn_in:sim,], type='l')
plot(beta_pp[burn_in:sim,], type='l')

plot(log_posteriors[burn_in:sim], type='l')
abline(v=which.max(log_posteriors[burn_in:sim]))
plot(log_posteriors_2[burn_in:sim], type='l')
plot(log_posteriors_3[burn_in:sim], type='l')
plot(log_posteriors_4[burn_in:sim], type='l')

