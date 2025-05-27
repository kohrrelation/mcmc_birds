

######################################################
###### 1. Functions for GEV components
######################################################

## with constant to prevent numerical underflow/overflow

gev_like_fn <- function(x,mu,sigma,xi, constant=0){
  #A <- ( (1+xi*(x-mu)/sigma)^(-1/xi) )
  if (any( (1+xi*(x-mu)/sigma)  < 0,na.rm=TRUE)==TRUE){
    logL <- -Inf
  } else {
    logL <- ( -constant + log( 1/sigma * ( (1+xi*(x-mu)/sigma)^(-1/xi) )^(xi+1) * exp(constant-( (1+xi*(x-mu)/sigma)^(-1/xi) )) ) )
  }
  #return(A)
  if (xi==0){
    logL <- ( -constant+ log( 1/sigma * ( exp(-(x-mu)/sigma)  )^(xi+1) * exp(constant-( exp(-(x-mu)/sigma) ) ) ) )
  }
  return(logL)
}

# fExtremes::dgev(1.5,xi=1,mu=1,beta=1.5, log=TRUE)
#
# evd::dgev(x=1.5,loc=0.6,scale=0.1,shape=0.01, log=TRUE)
# gev_like_fn(x=1.5,mu=0.6,sigma=0.1,xi=0.01)
#
# gev_like_fn(x=1.5,mu=0.6,sigma=0.1,xi=-0.5)
# mev::gev.score(par=c(0.6,0.1,1), dat=1.5)
# gev_like_d_fn(x=1.5,mu=0.6,sigma=0.1,xi=1)
# gev_like_dsigma_fn(x=1.5,mu=0.6,sigma=0.1,xi=1)

gev_like_d_fn <- function(x,mu,sigma,xi){
  #A <- ( (1+xi*(x-mu)/sigma)^(-1/xi) )
  # if (( (1+xi*(x-mu)/sigma) ) < 0){
  #   logL_d <- -Inf
  # } else {
  logL_d <-  (1/sigma * ((1 + xi * (x - mu)/sigma)^(-1/xi))^(xi + 1) * (exp(-((1 +
                                                                                 xi * (x - mu)/sigma)^(-1/xi))) * ((1 + xi * (x - mu)/sigma)^((-1/xi) -
                                                                                                                                                1) * ((-1/xi) * (xi/sigma)))) - 1/sigma * (((1 + xi * (x -
                                                                                                                                                                                                         mu)/sigma)^(-1/xi))^((xi + 1) - 1) * ((xi + 1) * ((1 + xi *
                                                                                                                                                                                                                                                              (x - mu)/sigma)^((-1/xi) - 1) * ((-1/xi) * (xi/sigma))))) *
                exp(-((1 + xi * (x - mu)/sigma)^(-1/xi))))/(1/sigma * ((1 +
                                                                          xi * (x - mu)/sigma)^(-1/xi))^(xi + 1) * exp(-((1 + xi *
                                                                                                                            (x - mu)/sigma)^(-1/xi))))
  #   }
  #return(A)
  if (xi==0){
    logL_d <- ( (1/sigma * ((exp(-(x - mu)/sigma))^((xi + 1) - 1) * ((xi + 1) *
                                                                       (exp(-(x - mu)/sigma) * (1/sigma)))) * exp(-(exp(-(x - mu)/sigma))) -
                   1/sigma * (exp(-(x - mu)/sigma))^(xi + 1) * (exp(-(exp(-(x -
                                                                              mu)/sigma))) * (exp(-(x - mu)/sigma) * (1/sigma))))/(1/sigma *
                                                                                                                                     (exp(-(x - mu)/sigma))^(xi + 1) * exp(-(exp(-(x - mu)/sigma)))) )
  }
  return(logL_d)

}


gev_like_d_fn <- function(x,mu,sigma,xi, constant=0){
  logL_d <-  (1/sigma * ((1 + xi * (x - mu)/sigma)^(-1/xi))^(xi + 1) * (exp(constant -
                                                                              ((1 + xi * (x - mu)/sigma)^(-1/xi))) * ((1 + xi * (x - mu)/sigma)^((-1/xi) -
                                                                                                                                                   1) * ((-1/xi) * (xi/sigma)))) - 1/sigma * (((1 + xi * (x -
                                                                                                                                                                                                            mu)/sigma)^(-1/xi))^((xi + 1) - 1) * ((xi + 1) * ((1 + xi *
                                                                                                                                                                                                                                                                 (x - mu)/sigma)^((-1/xi) - 1) * ((-1/xi) * (xi/sigma))))) *
                exp(constant - ((1 + xi * (x - mu)/sigma)^(-1/xi))))/(1/sigma *
                                                                        ((1 + xi * (x - mu)/sigma)^(-1/xi))^(xi + 1) * exp(constant -
                                                                                                                             ((1 + xi * (x - mu)/sigma)^(-1/xi))))
  #   }
  #return(A)
  if (xi==0){
    logL_d <- (1/sigma * ((exp(-(x - mu)/sigma))^((xi + 1) - 1) * ((xi + 1) *
                                                                     (exp(-(x - mu)/sigma) * (1/sigma)))) * exp(constant - (exp(-(x -
                                                                                                                                    mu)/sigma))) - 1/sigma * (exp(-(x - mu)/sigma))^(xi + 1) *
                 (exp(constant - (exp(-(x - mu)/sigma))) * (exp(-(x - mu)/sigma) *
                                                              (1/sigma))))/(1/sigma * (exp(-(x - mu)/sigma))^(xi +
                                                                                                                1) * exp(constant - (exp(-(x - mu)/sigma))))
  }
  return(logL_d)

}


## FOR SIGMA

gev_like_dsigma_fn <- function(x,mu,sigma,xi, constant=0){
  logL_d <-  (1/sigma * ((1 + xi * (x - mu)/sigma)^(-1/xi))^(xi + 1) * (exp(constant -
                                                                              ((1 + xi * (x - mu)/sigma)^(-1/xi))) * ((1 + xi * (x - mu)/sigma)^((-1/xi) -
                                                                                                                                                   1) * ((-1/xi) * (xi * (x - mu)/sigma^2)))) - (1/sigma * (((1 +
                                                                                                                                                                                                                xi * (x - mu)/sigma)^(-1/xi))^((xi + 1) - 1) * ((xi + 1) *
                                                                                                                                                                                                                                                                  ((1 + xi * (x - mu)/sigma)^((-1/xi) - 1) * ((-1/xi) * (xi *
                                                                                                                                                                                                                                                                                                                           (x - mu)/sigma^2))))) + 1/sigma^2 * ((1 + xi * (x - mu)/sigma)^(-1/xi))^(xi +
                                                                                                                                                                                                                                                                                                                                                                                                      1)) * exp(constant - ((1 + xi * (x - mu)/sigma)^(-1/xi))))/(1/sigma *
                                                                                                                                                                                                                                                                                                                                                                                                                                                                    ((1 + xi * (x - mu)/sigma)^(-1/xi))^(xi + 1) * exp(constant -
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         ((1 + xi * (x - mu)/sigma)^(-1/xi))))

  if (xi==0){
    logL_d <- ((1/sigma * ((exp(-(x - mu)/sigma))^((xi + 1) - 1) * ((xi + 1) *
                                                                      (exp(-(x - mu)/sigma) * ((x - mu)/sigma^2)))) - 1/sigma^2 *
                  (exp(-(x - mu)/sigma))^(xi + 1)) * exp(constant - (exp(-(x -
                                                                             mu)/sigma))) - 1/sigma * (exp(-(x - mu)/sigma))^(xi + 1) *
                 (exp(constant - (exp(-(x - mu)/sigma))) * (exp(-(x - mu)/sigma) *
                                                              ((x - mu)/sigma^2))))/(1/sigma * (exp(-(x - mu)/sigma))^(xi +
                                                                                                                         1) * exp(constant - (exp(-(x - mu)/sigma))))
  }
  return(logL_d)

}


gev_like_d <- function(data,mu_func,sigma,xi,constant=0, centroidassignment=centroidassignment,fe_mu=NULL){

  # if (!is.null(fe_mu)){
  #   covariate <- as.matrix( data[2:nrow(data),] )
  #   mu <- mu_func(covariate%*%fe_mu)
  # } else {
  #   mu <- mu_func(0)
  # }

  mu <- mu_func

  logL_d_test <- gev_like_d_fn(x=data[1,],mu=mu[centroidassignment],sigma=sigma[centroidassignment],xi=xi,constant)
  logL_test <- gev_like_fn(data[1,],mu[centroidassignment],sigma[centroidassignment],xi,constant)

  logL_d <- data.frame( val=rep(0, length(sigma)) , Category=1:length(sigma))
  logL <- data.frame( val=rep(0, length(sigma)) , Category=1:length(sigma))

  test_d <- rbind(logL_d,data.frame(val=logL_d_test,Category=centroidassignment) )
  test <- rbind(logL,data.frame(val=logL_test,Category=centroidassignment) )

  log_L <- aggregate(test$val, by=list(Category=test$Category), FUN=sum)
  logL_d <- aggregate(test_d$val, by=list(Category=test_d$Category), FUN=sum)

  terms <- list()
  terms$loglik_d <- logL_d[,2]
  terms$loglik <- log_L[,2]

  return(terms)
}



gev_like_d_list <- function(data_list,mu_func,sigma,xi,constant=0,fe_mu=NULL, g_d_func){
  ## use mapply maybe?
  mat_loglik <- matrix(NA,nrow=length(data_list$data), ncol=length(sigma))
  mat_loglik_d <- matrix(NA,nrow=length(data_list$data), ncol=length(sigma))

  for (j in 1:length(data_list$data)){

    if (!is.null(fe_mu)){
      covariate <- as.matrix( data_list$data[[j]][2:nrow(data_list$data[[j]]),] )
      mu <- mu_func(covariate%*%fe_mu)[[j]]
    } else {
      mu <- mu_func(0)[[j]]
    }

    loglike_all <- gev_like_d(data=data_list$data[[j]],mu_func=mu, sigma=sigma, xi=xi, constant=constant,
                              centroidassignment=data_list$centroidassignment_list[[j]],fe_mu)
    # chain rule
    mat_loglik[j,] <- loglike_all$loglik
    mat_loglik_d[j,] <- loglike_all$loglik_d*g_d_func[[j]]
  }

  terms <- list()
  terms$loglik_d <- apply(mat_loglik_d, 2, sum)
  terms$loglik <- apply(mat_loglik, 2, sum)

  return(terms)
}



gev_like_d_list_time <- function(data_list,mu_func,sigma,xi,constant=0,fe_mu=NULL, g_d_func){
  ## use mapply maybe?
  mat_loglik <- matrix(NA,nrow=length(data_list$data), ncol=length(sigma))
  mat_loglik_d <- matrix(NA,nrow=length(data_list$data), ncol=length(sigma))

  for (j in 1:length(data_list$data)){

    if (!is.null(fe_mu)){
      covariate <- as.matrix( data_list$data[[j]][2:nrow(data_list$data[[j]]),] )
      mu <- mu_func(covariate%*%fe_mu)[[j]]
    } else {
      mu <- mu_func(0)[[j]]
    }

    loglike_all <- gev_like_d(data=data_list$data[[j]],mu_func=mu, sigma=sigma, xi=xi, constant=constant,
                              centroidassignment=data_list$centroidassignment_list[[j]],fe_mu)
    # chain rule
    mat_loglik[j,] <- loglike_all$loglik
    mat_loglik_d[j,] <- loglike_all$loglik_d*g_d_func[[j]]
  }

  terms <- list()
  terms$loglik_d <- apply(mat_loglik_d, 1, sum)
  terms$loglik <- apply(mat_loglik, 1, sum)

  return(terms)
}



gev_like_d_theta2_list <- function(data_list,mu_func,sigma,xi,constant=0,fe_mu=NULL, g_d_theta2_func){
  ## use mapply maybe?
  mat_loglik <- matrix(NA,nrow=length(data_list$data), ncol=length(sigma))
  mat_loglik_d <- matrix(NA,nrow=length(data_list$data), ncol=length(sigma))

  for (j in 1:length(data_list$data)){

    if (!is.null(fe_mu)){
      covariate <- as.matrix( data_list$data[[j]][2:nrow(data_list$data[[j]]),] )
      mu <- mu_func(covariate%*%fe_mu)[[j]]
    } else {
      mu <- mu_func(0)[[j]]
    }

    loglike_all <- gev_like_d(data=data_list$data[[j]],mu_func=mu, sigma=sigma, xi=xi, constant=constant,
                              centroidassignment=data_list$centroidassignment_list[[j]],fe_mu)
    # chain rule
    mat_loglik[j,] <- loglike_all$loglik
    mat_loglik_d[j,] <- loglike_all$loglik_d*g_d_theta2_func[[j]]
  }

  terms <- list()
  terms$loglik_d <- apply(mat_loglik_d, 2, sum)
  terms$loglik <- apply(mat_loglik, 2, sum)

  return(terms)
}



gev_like_dsigma <- function(data,mu_func,sigma,xi,constant=0, centroidassignment=centroidassignment,fe_mu=NULL){

  # if (!is.null(fe_mu)){
  #   covariate <- as.matrix( data[2:nrow(data_list$data[[j]]),] )
  #   mu <- mu_func(covariate%*%fe_mu)
  # } else {
  #   mu <- mu_func(0)
  # }
  mu <- mu_func

  logL_d_test <- gev_like_dsigma_fn(data[1,],mu[centroidassignment],sigma[centroidassignment],xi,constant)
  logL_test <- gev_like_fn(data[1,],mu[centroidassignment],sigma[centroidassignment],xi,constant)

  logL_d <- data.frame( val=rep(0, length(sigma)) , Category=1:length(sigma))
  logL <- data.frame( val=rep(0, length(sigma)) , Category=1:length(sigma))

  test_d <- rbind(logL_d,data.frame(val=logL_d_test,Category=centroidassignment) )
  test <- rbind(logL,data.frame(val=logL_test,Category=centroidassignment) )

  log_L <- aggregate(test$val, by=list(Category=test$Category), FUN=sum)
  logL_d <- aggregate(test_d$val, by=list(Category=test_d$Category), FUN=sum)

  terms <- list()
  terms$loglik_d <- logL_d[,2]
  terms$loglik <- log_L[,2]

  return(terms)
}


gev_like_dsigma_list <- function(data_list,mu_func,sigma,xi,constant=0,fe_mu=NULL){
  ## use mapply maybe?
  mat_loglik <- matrix(NA,nrow=length(data_list$data), ncol=length(sigma))
  mat_loglik_d <- matrix(NA,nrow=length(data_list$data), ncol=length(sigma))

  for (j in 1:length(data_list$data)){

    if (!is.null(fe_mu)){
      covariate <- as.matrix( data_list$data[[j]][2:nrow(data_list$data[[j]]),] )
      mu <- mu_func(covariate%*%fe_mu)[[j]]
    } else {
      mu <- mu_func(0)[[j]]
    }

    loglike_all <- gev_like_dsigma(data=data_list$data[[j]],mu_func=mu, sigma=sigma, xi=xi, constant=constant,
                                   centroidassignment=data_list$centroidassignment_list[[j]],fe_mu)
    mat_loglik[j,] <- loglike_all$loglik
    mat_loglik_d[j,] <- loglike_all$loglik_d
  }


  terms <- list()
  terms$loglik_d <- apply(mat_loglik_d, 2, sum)
  terms$loglik <- apply(mat_loglik, 2, sum)

  return(terms)
}



gev_like_only <- function(data,mu_func,sigma,xi,constant=0, centroidassignment=centroidassignment, fe_mu=NULL){

  # if (!is.null(fe_mu)){
  #   covariate <- as.matrix( data[2:nrow(data_list$data[[j]]),] )
  #   mu <- mu_func(covariate%*%fe_mu)
  # } else {
  #   mu <- mu_func(0)
  # }
  mu <- mu_func

  logL_test <- gev_like_fn(data[1,],mu[centroidassignment],sigma[centroidassignment],xi,constant)

  logL <- data.frame( val=rep(0, length(sigma)) , Category=1:length(sigma))

  test <- rbind(logL,data.frame(val=logL_test,Category=centroidassignment) )

  log_L <- aggregate(test$val, by=list(Category=test$Category), FUN=sum)

  terms <- list()
  terms$loglik <- log_L[,2]

  return(terms)
}


gev_like_list <- function(data_list,mu_func,sigma,xi,constant=0,fe_mu=NULL){
  ## use mapply maybe?
  mat_loglik <- matrix(NA,nrow=length(data_list$data), ncol=length(sigma))

  for (j in 1:length(data_list$data)){

    if (!is.null(fe_mu)){
      covariate <- as.matrix( data_list$data[[j]][2:nrow(data_list$data[[j]]),] )
      mu <- mu_func(covariate%*%fe_mu)[[j]]
    } else {
      mu <- mu_func(0)[[j]]
    }

    loglike_all <- gev_like_only(data=data_list$data[[j]],mu_func=mu, sigma=sigma, xi=xi, constant=constant,
                                 centroidassignment=data_list$centroidassignment_list[[j]], fe_mu)
    mat_loglik[j,] <- loglike_all$loglik
  }


  terms <- apply(mat_loglik, 2, sum)

  return(terms)
}





######################################################
###### 2. Functions for PP components
######################################################
# LOG LIKE #


hpp_likelihood <- function( lambda, size_A, centroids_count, fe_pp=NULL){
  n_i <- centroids_count[,3]
  if (!is.null(fe_pp)){
    covariate <- as.matrix( centroids_count[,4:ncol(centroids_count)] )
    lambda <- matrix(lambda, nrow=nrow(centroids_count))+
      fe_pp[covariate]
  } else {
    lambda <- matrix(lambda, nrow=nrow(centroids_count))
  }

  loglik <-  sum(( lambda )* n_i - size_A* ( exp(lambda) ))
  # this term here is with the intergral approximated
  terms <- list()
  terms$loglik <- loglik
  terms$lambda <- lambda
  return( terms )
}

## THIS IS NOT WORKING. TO FIX OR DON't USE
hpp_likelihood_list <- function(lambda, size_A, centroids_count_list,fe_pp=NULL){
  term <- lapply(centroids_count_list, function(x) hpp_likelihood( lambda, size_A, x, fe_pp) )

  terms <- list()
  terms$loglik <- term$loglik
  terms$lambda <- lapply( term, function(x) x$lambda)
  return(terms)
}


hpp_likelihood_grad <- function( lambda, size_A, centroids_count, fe_pp=NULL){
  n_i <- centroids_count[,3]

  if (!is.null(fe_pp)){
    covariate <- as.matrix( centroids_count[,4:ncol(centroids_count)] )
    lambda <- matrix(lambda, nrow=nrow(centroids_count))+
      fe_pp[covariate]
  } else {
    lambda <- matrix(lambda, nrow=nrow(centroids_count))
  }

  loglik <-  ( lambda*n_i - size_A* ( exp(lambda) ) )

  loglik_d <- (n_i - size_A*exp(lambda)*1)

  terms <- list()
  terms$loglik_d <- loglik_d
  terms$loglik <- loglik
  terms$lambda <- lambda

  return(terms)
}


hpp_likelihood_d_list <- function(lambda, size_A, centroids_count_list, fe_pp=NULL){
  term <- lapply(centroids_count_list, function(x) hpp_likelihood_grad( lambda, size_A, x,
                                                                        fe_pp) )
  loglik <- apply ( matrix( unlist( lapply( term, function(x) x$loglik) ),
                            nrow=length(lambda)) , 1, sum)
  loglik_d <- apply ( matrix( unlist( lapply( term, function(x) x$loglik_d) ),
                              nrow=length(lambda)) , 1, sum)
  terms <- list()
  terms$loglik_d <- loglik_d
  terms$loglik <- loglik
  terms$lambda <- lapply( term, function(x) x$lambda)
  return(terms)
}



hpp_likelihood_d_list_time <- function(lambda, size_A, centroids_count_list, fe_pp=NULL){
  term <- lapply(centroids_count_list, function(x) hpp_likelihood_grad( lambda, size_A, x,
                                                                        fe_pp) )
  loglik <- apply ( matrix( unlist( lapply( term, function(x) x$loglik) ),
                            nrow=length(lambda)) , 2, sum)
  loglik_d <- apply ( matrix( unlist( lapply( term, function(x) x$loglik_d) ),
                              nrow=length(lambda)) , 2, sum)
  terms <- list()
  terms$loglik_d <- loglik_d
  terms$loglik <- loglik
  terms$lambda <- lapply( term, function(x) x$lambda)
  return(terms)
}


get_intensity <- function( lambda, size_A, centroids_count, fe_pp=NULL){
  n_i <- centroids_count[,3]
  if (!is.null(fe_pp)){
    covariate <- as.matrix( centroids_count[,4:ncol(centroids_count)] )
    lambda <- matrix(lambda, nrow=nrow(centroids_count))+
      fe_pp[covariate]
  } else {
    lambda <- matrix(lambda, nrow=nrow(centroids_count))
  }
  terms <- list()
  terms$lambda <- lambda
  return(terms)
}

get_intensity_list <- function(lambda, size_A, centroids_count_list, fe_pp=NULL){
  term <- lapply(centroids_count_list, function(x) get_intensity( lambda, size_A, x,
                                                                        fe_pp) )
  return(lapply( term, function(x) x$lambda) )
}




######################################################
###### 3. Functions for prior
######################################################
# log-likelihood in the second (latent) layer.

hpp_likelihood_2 <- function( alpha_l, lambda_l, lambda, neigh,intercept=0){
  covparms <- c( (alpha_l),(lambda_l), 0)
  loglik <- vecchia_meanzero_loglik( covparms, "exponential_isotropic",
                                     lambda[ord]-intercept,
                                     (locsord), NNarray )
  return(loglik$loglik)
}

# for RE
hpp_likelihood_2_time <- function( alpha_l, lambda_l, lambda, neigh,intercept=0){
  covparms <- c( (alpha_l),(lambda_l), 0)
  loglik <- vecchia_meanzero_loglik( covparms, "exponential_isotropic",
                                     lambda[ord_time]-intercept,
                                     matrix(c(1:21)[ord_time], ncol=1), NNarray_time )
  return(loglik$loglik)
}

#
# prior_igamma <- function(x, shape=1, rate=12){
#   return(MCMCpack::dinvgamma((x), shape=shape, scale = rate) )
# }
#
# prior_gamma <- function(x, shape=5, rate=3){
#   return(dgamma(x=(x), shape=shape, rate=rate, log=TRUE) )
# }


# prior_normal <- function(x, mu, sigma){
#   # x is column vector of observations
#   # mu is row vector of parameters
#   # sigma is matrix
#   x_re <- as.numeric(x)
#   return(dmvnorm(x=x_re, mean=mu, sigma=sigma, log=TRUE) )
# }


prior_normal <- function(x, mu, sigma){
  # x is column vector of observations
  # mu is row vector of parameters
  # sigma is matrix
  x_re <- as.numeric(x)
  return(dnorm(x=x_re, mean=mu, sd=sigma, log=TRUE) )
}




pc_prior <- function(sigma,rho, sigma_0, rho_0, alpha_1=0.5, alpha_2=0.5){
  lambda_1 <- -log(alpha_1)*rho_0
  lambda_2 <- -log(alpha_2)/sigma_0

  lambda_1*lambda_2*rho^(-2)*exp(-lambda_1*rho^(-1) - lambda_2*sigma)
}





### Poisson likelihood


poisson_likelihood_list <- function(lambda, centroids_count_list,fe_pp=NULL){
  lapply(centroids_count_list, function(x) poisson_likelihood( lambda, x, fe_pp) )
}

poisson_likelihood_d_list <- function(lambda, centroids_count_list, fe_pp=NULL){
  term <- lapply(centroids_count_list, function(x) poisson_likelihood_grad( lambda, x,
                                                                            fe_pp) )
  loglik <- apply ( matrix( unlist( lapply( term, function(x) x$loglik) ),
                            nrow=nrow(centroids_count_list[[2]])) , 1, sum)
  loglik_d <- apply ( matrix( unlist( lapply( term, function(x) x$loglik_d) ),
                              nrow=nrow(centroids_count_list[[2]])) , 1, sum)
  terms <- list()
  terms$loglik_d <- loglik_d
  terms$loglik <- loglik
  return(terms)
}

## new one where just remove n_stop=0 data

poisson_likelihood <- function( lambda, centroids_count, fe_pp=NULL){
  n_i <- centroids_count[,2]
  n_stop <- centroids_count[,3]

  # if covariate is on route level. if it is on grid level then need to adapt. or other way...
  if (!is.null(fe_pp)){
    covariate <- as.matrix( centroids_count[,3:ncol(centroids_count)] )
    lambda <- matrix(lambda, nrow=n.loc)+
      (covariate)%*%fe_pp
  } else {
    lambda <- matrix(lambda, nrow=n.loc)
  }

  lambda_route <-  weight_route%*%exp(lambda)
  loglik <-  (log(n_stop*lambda_route )* n_i - ( n_stop*lambda_route ))
  loglik[which(n_stop==0)] <- 0
  return( loglik )
}


poisson_likelihood_grad <- function( lambda, centroids_count, fe_pp=NULL){
  n_i <- centroids_count[,2]
  n_stop <- centroids_count[,3]
  if (!is.null(fe_pp)){
    covariate <- as.matrix( centroids_count[,3:ncol(centroids_count)] )
    lambda <- matrix(lambda, nrow=n.loc)+
      (covariate)%*%fe_pp
  } else {
    lambda <- matrix(lambda, nrow=n.loc)
  }

  lambda_route <-  weight_route%*%exp(lambda)

  loglik <-  log(n_stop*lambda_route )* n_i - ( n_stop*lambda_route )
  loglik_d <- n_stop/(n_stop * lambda_route) * n_i - n_stop

  loglik[which(n_stop==0)] <- 0
  loglik_d[which(n_stop==0)] <- 0

  terms <- list()
  terms$loglik_d <- loglik_d
  terms$loglik <- loglik

  return(terms)
}




### Bionmial likelihood


bin_likelihood <- function( lambda, centroids_count, fe_bin=NULL){
  n_i <- centroids_count[,3]
  n <- centroids_count[,4]
  if (!is.null(fe_bin)){
    covariate <- as.matrix( centroids_count[,5:ncol(centroids_count)] )
    lambda <- matrix(lambda, nrow=nrow(centroids_count))+
      (1/  covariate )%*%fe_bin 
  } else {
    lambda <- matrix(lambda, nrow=nrow(centroids_count))
  }
  p <- 1- exp( -exp( lambda) )

  loglik <-  (n_i*log(p) + (n-n_i)*log(1-p)   )
  loglik[which(n_i==0)] <- 0

  loglik[which(n<n_i)] <- 0

  # this term here is with the intergral approximated
  return( loglik )
}


bin_likelihood_list <- function(lambda, centroids_count_list,fe_bin=NULL){
  lapply(centroids_count_list, function(x) bin_likelihood( lambda, x, fe_bin) )
}


bin_likelihood_grad <- function( lambda, centroids_count, fe_bin=NULL){
  n_i <- centroids_count[,3]
  n <- centroids_count[,4]
  if (!is.null(fe_bin)){
    covariate <- as.matrix( centroids_count[,5:ncol(centroids_count)] )
    lambda <- matrix(lambda, nrow=nrow(centroids_count))+
      (1/  covariate )%*%fe_bin 
  } else {
    lambda <- matrix(lambda, nrow=nrow(centroids_count))
  }
  p <- 1- exp( -exp( lambda) )
  #chain rule: D(expression(1- exp( -exp( lambda) )), 'lambda')
  loglik <-  ( (n_i*log(p) + (n-n_i)*log(1-p))   )
  loglik_d <-  ( n_i/p - (n-n_i)/(1-p)  )*exp(-exp(lambda)) * exp(lambda)

  # n_i * (exp(-exp(lambda)) * exp(lambda)/(1 - exp(-exp(lambda)))) +
  #   (n - n_i) * (exp(-exp(lambda)) * exp(lambda)/(1 - 1 - exp(-exp(lambda))))

  loglik[which(n==0)] <- 0
  loglik_d[which(n==0)] <- 0

  loglik[which(n<n_i)] <- 0
  loglik_d[which(n<n_i)] <- 0

  terms <- list()
  terms$loglik_d <- loglik_d
  terms$loglik <- loglik
  terms$effort <- 1/( (covariate) )

  # this term here is with the intergral approximated
  return( terms )
}



bin_likelihood_d_list <- function(lambda, centroids_count_list, fe_bin=NULL){
  term <- lapply(centroids_count_list, function(x) bin_likelihood_grad( lambda, x,
                                                                            fe_bin) )
  loglik <- apply ( matrix( unlist( lapply( term, function(x) x$loglik) ),
                            nrow=length(lambda)) , 1, sum)
  loglik_d <- apply ( matrix( unlist( lapply( term, function(x) x$loglik_d) ),
                              nrow=length(lambda)) , 1, sum)
  terms <- list()
  terms$loglik_d <- loglik_d
  terms$loglik <- loglik
  terms$effort <- lapply( term, function(x) x$effort)
  return(terms)
}



get_effort <- function( lambda, centroids_count, fe_bin=NULL){

  if (!is.null(fe_bin)){
    covariate <- as.matrix( centroids_count[,5:ncol(centroids_count)] )
    # lambda <- matrix(lambda, nrow=nrow(centroids_count))+
    #   1/ ( (covariate)%*%fe_bin )
  } else {
    # lambda <- matrix(lambda, nrow=nrow(centroids_count))
  }

  terms <- list()
  terms$effort <- 1/( (covariate) )

  return( terms )
}


get_effort_list <- function(lambda, centroids_count_list, fe_bin=NULL){
  term <- lapply(centroids_count_list, function(x) get_effort( lambda, x,
                                                                        fe_bin) )
  terms <- list()
  terms$effort <- lapply( term, function(x) x$effort)
  return(terms)
}

### Bernouli


### Bionmial likelihood


# poisson_likelihood <- function( lambda, centroids_count, fe_bin=NULL){
#   n_i <- centroids_count[,2]
#   n_stop <- centroids_count[,3]
#   n <- 1
#   if (!is.null(fe_bin)){
#     lambda <- matrix(lambda, nrow=nrow(centroids_count))
#   } else {
#     lambda <- matrix(lambda, nrow=nrow(centroids_count))
#   }
#   p <- n_stop*(1- exp( -exp( lambda) ))
#
#   loglik <-  (n_i*log(p) + (n-n_i)*log(1-p)   )
#   loglik[which(n_i==0)] <- 0
#
#   loglik[which(n<n_i)] <- 0
#
#   # this term here is with the intergral approximated
#   return( loglik )
# }




# poisson_likelihood_grad <- function( lambda, centroids_count, fe_bin=NULL){
#   n_i <- centroids_count[,2]
#   n_stop <- centroids_count[,3]
#   n <- 1
#   if (!is.null(fe_bin)){
#     lambda <- matrix(lambda, nrow=nrow(centroids_count))
#   } else {
#     lambda <- matrix(lambda, nrow=nrow(centroids_count))
#   }
#   p <- n_stop*(1- exp( -exp( lambda) ))
#   #chain rule: D(expression(1- exp( -exp( lambda) )), 'lambda')
#   #chain rule: D(expression(n_stop*(1- exp( -exp( lambda) ))), 'lambda')
#   loglik <-  ( (n_i*log(p) + (n-n_i)*log(1-p))   )
#   loglik_d <-  ( n_i/p - (n-n_i)/(1-p)  )*n_stop * (exp(-exp(lambda)) * exp(lambda))
#
#   loglik[which(n==0)] <- 0
#   loglik_d[which(n==0)] <- 0
#
#   loglik[which(n<n_i)] <- 0
#   loglik_d[which(n<n_i)] <- 0
#
#   terms <- list()
#   terms$loglik_d <- loglik_d
#   terms$loglik <- loglik
#   # this term here is with the intergral approximated
#   return( terms )
# }

### Bionmial likelihood


# poisson_likelihood <- function( lambda, centroids_count, fe_bin=NULL){
#   n_i <- centroids_count[,2]
#   n_stop <- centroids_count[,3]
#   n <- 1
#   if (!is.null(fe_bin)){
#     lambda <- matrix(lambda, nrow=nrow(centroids_count))
#   } else {
#     lambda <- matrix(lambda, nrow=nrow(centroids_count))
#   }
#   p <- (1- exp( -exp( lambda) ))
#
#   loglik <-  (n_i*log(p) + (n-n_i)*log(1-p)   )
#   loglik[which(n_i==0)] <- 0
#   loglik[which(n<n_i)] <- 0
#   loglik[which(n_stop==0)] <- 0
#
#   # this term here is with the intergral approximated
#   return( loglik )
# }
#
#
# poisson_likelihood_list <- function(lambda, centroids_count_list,fe_bin=NULL){
#   lapply(centroids_count_list, function(x) poisson_likelihood( lambda, x, fe_bin) )
# }
#
#
# poisson_likelihood_grad <- function( lambda, centroids_count, fe_bin=NULL){
#   n_i <- centroids_count[,2]
#   n_stop <- centroids_count[,3]
#   n <- 1
#   if (!is.null(fe_bin)){
#     lambda <- matrix(lambda, nrow=nrow(centroids_count))
#   } else {
#     lambda <- matrix(lambda, nrow=nrow(centroids_count))
#   }
#   p <- (1- exp( -exp( lambda) ))
#   #chain rule: D(expression(1- exp( -exp( lambda) )), 'lambda')
#   #chain rule: D(expression(n_stop*(1- exp( -exp( lambda) ))), 'lambda')
#   loglik <-  ( (n_i*log(p) + (n-n_i)*log(1-p))   )
#   loglik_d <-  ( n_i/p - (n-n_i)/(1-p)  )
#
#   loglik[which(n==0)] <- 0
#   loglik_d[which(n==0)] <- 0
#
#   loglik[which(n<n_i)] <- 0
#   loglik_d[which(n<n_i)] <- 0
#
#   loglik[which(n_stop==0)] <- 0
#   loglik_d[which(n_stop==0)] <- 0
#
#   terms <- list()
#   terms$loglik_d <- loglik_d
#   terms$loglik <- loglik
#   # this term here is with the intergral approximated
#   return( terms )
# }
#
#
#
# poisson_likelihood_d_list <- function(lambda, centroids_count_list, fe_bin=NULL){
#   term <- lapply(centroids_count_list, function(x) poisson_likelihood_grad( lambda, x,
#                                                                         fe_bin) )
#   loglik <- apply ( matrix( unlist( lapply( term, function(x) x$loglik) ),
#                             nrow=length(lambda)) , 1, sum)
#   loglik_d <- apply ( matrix( unlist( lapply( term, function(x) x$loglik_d) ),
#                               nrow=length(lambda)) , 1, sum)
#   terms <- list()
#   terms$loglik_d <- loglik_d
#   terms$loglik <- loglik
#   return(terms)
# }



######################################################
###### 4. MAIN MCMC CODE
######################################################

# stepsize for hyperparams
# delta for latent GPs, ind_latent_gp
# stepsize_intercept for model intercepts, ind_intercept
# stepsize_share for sharing params, ind_share
# stepsize_fe for FE/trends, ind_fe

MCMC_function_GP <- function(sim_ind, hess, delta, stepsize, log_posteriors,
                             hyper_gp_latent_pp, hyper_gp_latent_pp2, hyper_gp_latent_mu, hyper_gp_latent_sigma,
                             x_pp, x_pp_2, x_mu, x_sigma,
                             beta_pp, beta_mu, beta_sigma, beta_xi, beta_share, beta_pp_2, beta_bin,
                             fe_mu,fe_pp,fe_pp2, fe_bin, log_posteriors_2, log_posteriors_3, log_posteriors_4,
                             hyper_gp_latent_bin, x_bin,hyper_gp_latent_re){

  sim <- 1

  ## IND !!
  # update hyperparameters first

  # hyper_gp_latent=hyper_gp_latent_sigma
  # x=x_sigma
  # long_range='short'

  update_hyper <- function(hyper_gp_latent, x, intercept=0, long_range='long', type_re='space'){

    hyper_gp_latent[sim+1,1] <- exp(rnorm(1, log(hyper_gp_latent[sim,1]), stepsize[2]))

    if (long_range=='long'){
      prior_terms <- log(pc_prior(sigma=sqrt(hyper_gp_latent[sim+1,1]),
                                  rho=hyper_gp_latent[sim,2]/2, sigma_0=sqrt(sigma.sq_long), rho_0=phi_long/2 , alpha_1=0.5, alpha_2=0.5))-
        log(pc_prior(sigma=sqrt(hyper_gp_latent[sim,1]),
                     rho=hyper_gp_latent[sim,2]/2, sigma_0=sqrt(sigma.sq_long), rho_0=phi_long/2 , alpha_1=0.5, alpha_2=0.5))
    } else if (long_range=='medium'){
      prior_terms <- log(pc_prior(sigma=sqrt(hyper_gp_latent[sim+1,1]),
                                  rho=hyper_gp_latent[sim,2]/2, sigma_0=sqrt(sigma.sq_medium), rho_0=phi_medium/2 , alpha_1=0.5, alpha_2=0.5)) -
        log(pc_prior(sigma=sqrt(hyper_gp_latent[sim,1]),
                     rho=hyper_gp_latent[sim,2]/2, sigma_0=sqrt(sigma.sq_medium), rho_0=phi_medium/2 , alpha_1=0.5, alpha_2=0.5))
    } else {
      prior_terms <- log(pc_prior(sigma=sqrt(hyper_gp_latent[sim+1,1]),
                                  rho=hyper_gp_latent[sim,2]/2, sigma_0=sqrt(sigma.sq_short), rho_0=phi_short/2 , alpha_1=0.5, alpha_2=0.5)) -
        log(pc_prior(sigma=sqrt(hyper_gp_latent[sim,1]),
                     rho=hyper_gp_latent[sim,2]/2, sigma_0=sqrt(sigma.sq_short), rho_0=phi_short/2 , alpha_1=0.5, alpha_2=0.5))
    }


    if (type_re=='space'){
      alpha <- +
        (dnorm(log(hyper_gp_latent[sim,1]), log(hyper_gp_latent[sim+1,1]), stepsize[1], log = T) + -log(hyper_gp_latent[sim,1])) -
        (dnorm(log(hyper_gp_latent[sim+1,1]), log(hyper_gp_latent[sim,1]), stepsize[1], log = T) + -log(hyper_gp_latent[sim+1,1]))+
        -(hpp_likelihood_2((hyper_gp_latent[sim,1]), (hyper_gp_latent[sim,2]), x[sim,], neigh, intercept=intercept) ) +
        (hpp_likelihood_2((hyper_gp_latent[sim+1,1]), (hyper_gp_latent[sim,2]), x[sim,], neigh, intercept=intercept) ) +
        prior_terms
    } else {
      alpha <- +
        (dnorm(log(hyper_gp_latent[sim,1]), log(hyper_gp_latent[sim+1,1]), stepsize[1], log = T) + -log(hyper_gp_latent[sim,1])) -
        (dnorm(log(hyper_gp_latent[sim+1,1]), log(hyper_gp_latent[sim,1]), stepsize[1], log = T) + -log(hyper_gp_latent[sim+1,1]))+
        -(hpp_likelihood_2_time((hyper_gp_latent[sim,1]), (hyper_gp_latent[sim,2]), x[sim,], neigh, intercept=intercept) ) +
        (hpp_likelihood_2_time((hyper_gp_latent[sim+1,1]), (hyper_gp_latent[sim,2]), x[sim,], neigh, intercept=intercept) ) +
        prior_terms
    }

    u <- runif(1)
    alpha <- log(min(1,exp(alpha)))

    if (log( u ) >alpha){
      ind[sim,1] <<- 0
      hyper_gp_latent[sim+1,1] <- hyper_gp_latent[sim,1] } else {ind[sim,1] <<- 1}

    # second hyperparameter

    hyper_gp_latent[sim+1,2] <- exp(rnorm(1, log(hyper_gp_latent[sim,2]), stepsize[2]))


    if (long_range=='long'){
      prior_terms <- log(pc_prior(sigma=sqrt(hyper_gp_latent[sim+1,1]),
                                  rho=hyper_gp_latent[sim+1,2]/2, sigma_0=sqrt(sigma.sq_long), rho_0=phi_long/2 , alpha_1=0.5, alpha_2=0.5))-
        log(pc_prior(sigma=sqrt(hyper_gp_latent[sim+1,1]),
                     rho=hyper_gp_latent[sim,2]/2, sigma_0=sqrt(sigma.sq_long), rho_0=phi_long/2 , alpha_1=0.5, alpha_2=0.5))
    } else if (long_range=='medium') {
      prior_terms <- log(pc_prior(sigma=sqrt(hyper_gp_latent[sim+1,1]),
                                  rho=hyper_gp_latent[sim+1,2]/2, sigma_0=sqrt(sigma.sq_medium), rho_0=phi_medium/2 , alpha_1=0.5, alpha_2=0.5)) -
        log(pc_prior(sigma=sqrt(hyper_gp_latent[sim+1,1]),
                     rho=hyper_gp_latent[sim,2]/2, sigma_0=sqrt(sigma.sq_medium), rho_0=phi_medium/2 , alpha_1=0.5, alpha_2=0.5))
    } else {
      prior_terms <- log(pc_prior(sigma=sqrt(hyper_gp_latent[sim+1,1]),
                                  rho=hyper_gp_latent[sim+1,2]/2, sigma_0=sqrt(sigma.sq_short), rho_0=phi_short/2 , alpha_1=0.5, alpha_2=0.5)) -
        log(pc_prior(sigma=sqrt(hyper_gp_latent[sim+1,1]),
                     rho=hyper_gp_latent[sim,2]/2, sigma_0=sqrt(sigma.sq_short), rho_0=phi_short/2 , alpha_1=0.5, alpha_2=0.5))
    }

    if (type_re=='space'){
      alpha <- (dnorm(log(hyper_gp_latent[sim,2]), log(hyper_gp_latent[sim+1,2]), stepsize[2], log = T) + -log(hyper_gp_latent[sim,2])) -
        (dnorm(log(hyper_gp_latent[sim+1,2]), log(hyper_gp_latent[sim,2]), stepsize[2], log = T) + -log(hyper_gp_latent[sim+1,2]))+
        (hpp_likelihood_2((hyper_gp_latent[sim+1,1]), (hyper_gp_latent[sim+1,2]), x[sim,], neigh, intercept=intercept) ) -
        (hpp_likelihood_2((hyper_gp_latent[sim+1,1]), (hyper_gp_latent[sim,2]), x[sim,], neigh, intercept=intercept) ) +
        prior_terms
    } else{
      alpha <- (dnorm(log(hyper_gp_latent[sim,2]), log(hyper_gp_latent[sim+1,2]), stepsize[2], log = T) + -log(hyper_gp_latent[sim,2])) -
        (dnorm(log(hyper_gp_latent[sim+1,2]), log(hyper_gp_latent[sim,2]), stepsize[2], log = T) + -log(hyper_gp_latent[sim+1,2]))+
        (hpp_likelihood_2_time((hyper_gp_latent[sim+1,1]), (hyper_gp_latent[sim+1,2]), x[sim,], neigh, intercept=intercept) ) -
        (hpp_likelihood_2_time((hyper_gp_latent[sim+1,1]), (hyper_gp_latent[sim,2]), x[sim,], neigh, intercept=intercept) ) +
        prior_terms
    }

    u <- runif(1)
    alpha <- log(min(1,exp(alpha)))


    if (log( u ) >alpha){
      ind[sim,2] <<- 0
      hyper_gp_latent[sim+1,2] <- hyper_gp_latent[sim,2] } else {ind[sim,2] <<- 1}

    return(hyper_gp_latent)
  }

  # effort
  hyper_gp_latent_pp <- update_hyper(hyper_gp_latent_pp, x_pp, long_range='long')
  ind_hyper_pp[sim_ind,] <<- ind[sim,]
  #hyper_gp_latent_pp[sim+1,] <- c(sigma.sq_long, phi_long)
  # niche
  hyper_gp_latent_pp2 <- update_hyper(hyper_gp_latent_pp2, x_pp_2, long_range='long')
  ind_hyper_pp2[sim_ind,] <<- ind[sim,]
  #hyper_gp_latent_pp2[sim+1,] <- c(sigma.sq_long, phi_long)
  #phenology
  hyper_gp_latent_mu <- update_hyper(hyper_gp_latent_mu, x_mu, long_range='medium')
  ind_hyper_mu[sim_ind,] <<- ind[sim,]
  #volatility
  hyper_gp_latent_sigma <- update_hyper(hyper_gp_latent_sigma, x_sigma, long_range='short')
  ind_hyper_sigma[sim_ind,] <<- ind[sim,]

  #time RE
  hyper_gp_latent_re <- update_hyper(hyper_gp_latent_re, fe_pp, long_range='short', type_re = 'time')
  ind_hyper_re[sim_ind,] <<- ind[sim,]

  #x_bin
  # hyper_gp_latent_bin <- update_hyper(hyper_gp_latent_bin, x_bin, long_range='medium')
  # ind_hyper_bin[sim_ind,] <<- ind[sim,]





  ################################################
  ### latent field - NICHE - x_bin
  ################################################


  # covparms <- c(hyper_gp_latent_bin[sim+1,1], hyper_gp_latent_bin[sim+1,2], 0)
  #
  # # bin likelihood
  # test <- bin_likelihood_d_list(beta_bin[sim,]+beta_share[sim,4]*x_pp_2[sim,]+x_bin[sim,], centroids_count_list_bin, fe_bin=fe_bin[sim,])
  #
  # L_inv <- vecchia_Linv(covparms, "exponential_isotropic", coords, NNarray)
  #
  # grad_prior_1 <- -Linv_mult( L_inv, Linv_t_mult(L_inv, x_bin[sim,], NNarray) , NNarray)
  #
  #
  # for (j in 1:n.sites) {
  #   x_bin[sim+1,j] <- rnorm(1, mean=x_bin[sim,j]+delta[4]^2/2*(test$loglik_d[j]+
  #                                                                  grad_prior_1[j])*hess[1], sd=delta[4]*sqrt(hess[1]))
  # }
  #
  # # bin likelihood
  # test_2 <- bin_likelihood_d_list(beta_bin[sim,]+beta_share[sim,4]*x_pp_2[sim,]+x_bin[sim+1,], centroids_count_list_bin, fe_bin=fe_bin[sim,])
  #
  # grad_prior_2 <- -Linv_mult( L_inv, Linv_t_mult(L_inv, x_bin[sim+1,], NNarray) , NNarray)
  #
  # prop_density_1 <- sum(dnorm(x_bin[sim,], mean=x_bin[sim+1,]+delta[4]^2/2*(test_2$loglik_d+grad_prior_2)*hess, sd=delta[4]*sqrt(hess), log=TRUE))
  # prop_density_2 <- sum(dnorm(x_bin[sim+1,], mean=x_bin[sim,]+delta[4]^2/2*(test$loglik_d+grad_prior_1)*hess, sd=delta[4]*sqrt(hess), log=TRUE))
  #
  # alpha <-
  #   # PROPOSALS
  #   prop_density_1 - prop_density_2 +
  #   sum(test_2$loglik)+
  #   hpp_likelihood_2(hyper_gp_latent_bin[sim+1,1], hyper_gp_latent_bin[sim+1,2],
  #                    x_bin[sim+1,], neigh, intercept=0) -
  #   # LIKELIHOOD
  #   sum(test$loglik) -
  #   # PRIORS
  #   hpp_likelihood_2(hyper_gp_latent_bin[sim+1,1], hyper_gp_latent_bin[sim+1,2],
  #                    x_bin[sim,], neigh, intercept=0 )
  #
  # # for -Inf likelihood. just stay away from that space
  # if (is.nan(alpha)){alpha=-10000}
  #
  # u <- runif(1)
  #
  #
  # alpha <- log(min(1,exp(alpha)))
  #
  # if (log( u ) >alpha){x_bin[sim+1,] <- x_bin[sim,]
  # ind_latent_gp[sim_ind,4] <<- 0
  # } else {ind_latent_gp[sim_ind,4] <<- 1}
  #
  # ################
  #
  # ###############
  # ## SUM TO ZERO CONSTRAINT
  # ###############
  #
  # x_bin[sim+1,]<- x_bin[sim+1,]-mean(x_bin[sim+1,]) # make sure that sum of x_pp is zero


  x_bin[sim+1,]<- 0


  ################################################
  ### latent field - SHARED PP-GEV - x_pp
  ################################################

  covparms <- c(hyper_gp_latent_pp[sim+1,1], hyper_gp_latent_pp[sim+1,2], 0)

  test <- hpp_likelihood_d_list(beta_pp[sim,]+x_pp[sim,],size_A, centroids_count_list,
                                fe_pp[sim,])

  share_intensity <- test$lambda

  test_bin <- get_effort_list(beta_bin[sim,]+beta_share[sim,4]*x_pp_2[sim,]+x_bin[sim+1,], centroids_count_list_bin,
                              fe_bin=fe_bin[sim,])

  share_intensity_1 <- list()
  for (j in 1:length(share_intensity)){
    share_intensity_1[[j]] <- beta_share[sim,1]*share_intensity[[j]]+beta_share[sim,5]*test_bin$effort[[j]]
  }

  #mu
  mu_func <- function(y){
    lapply(share_intensity_1,
           function(x) g_func(1,beta_share[sim,2],x,x_mu[sim,]+beta_mu[sim,]+beta_share[sim,3]*x_pp_2[sim,]+y)+
             0 )
  }

  g_d_func_list <- lapply(share_intensity_1,
                          function(x) beta_share[sim,1]*g_d_func(1,beta_share[sim,2],x,x_mu[sim,]+beta_mu[sim,]+beta_share[sim,3]*x_pp_2[sim,]) )


  test_gev <- gev_like_d_list(data_list=data_list,
                              mu_func=mu_func,
                              sigma=exp(beta_sigma[sim,]+x_sigma[sim,]), xi=beta_xi[sim,], constant=0,
                              fe_mu[sim,], g_d_func_list)

  #sum(test_gev$loglik)
  L_inv <- vecchia_Linv(covparms, "exponential_isotropic", locsord, NNarray)

  grad_prior_1 <- -Linv_mult( L_inv, Linv_t_mult(L_inv, x_pp[sim,ord], NNarray) , NNarray)


  ## DONT FORGET THE CHAIN RULE
  for (j in 1:n.sites) {
    x_pp[sim+1,j] <- rnorm(1, mean=x_pp[sim,j]+delta[1]^2/2*
                             (test_gev$loglik_d[j]+test$loglik_d[j]+grad_prior_1[j])*hess[1], sd=delta[1]*sqrt(hess[1]))
  }

  test_2 <- hpp_likelihood_d_list(beta_pp[sim,]+x_pp[sim+1,],size_A, centroids_count_list,
                                  fe_pp[sim,])
  share_intensity <- test_2$lambda


  share_intensity_2 <- list()
  for (j in 1:length(share_intensity)){
    share_intensity_2[[j]] <- beta_share[sim,1]*share_intensity[[j]]+beta_share[sim,5]*test_bin$effort[[j]]
  }


  mu_func <- function(y){
    lapply(share_intensity_2,
           function(x) g_func(1,beta_share[sim,2],x,x_mu[sim,]+beta_mu[sim,]+beta_share[sim,3]*x_pp_2[sim,]+y)+
             0 )
  }

  g_d_func_list <- lapply(share_intensity_2,
                          function(x) beta_share[sim,1]*g_d_func(1,beta_share[sim,2],x,x_mu[sim,]+beta_mu[sim,]+beta_share[sim,3]*x_pp_2[sim,]) )


  test_gev_2 <- gev_like_d_list(data=data_list,
                                mu_func=mu_func,
                                sigma=exp(beta_sigma[sim,]+x_sigma[sim,]), xi=beta_xi[sim,], constant=0,
                                fe_mu[sim,], g_d_func_list)

  grad_prior_2 <- -Linv_mult( L_inv, Linv_t_mult(L_inv, x_pp[sim+1,ord], NNarray) , NNarray)


  prop_density_1 <- sum(dnorm(x_pp[sim,],
                              mean=x_pp[sim+1,]+delta[1]^2/2*
                                (test_gev_2$loglik_d+test_2$loglik_d+grad_prior_2)*hess,
                              sd=delta[1]*sqrt(hess), log=TRUE))
  prop_density_2 <- sum(dnorm(x_pp[sim+1,],
                              mean=x_pp[sim,]+delta[1]^2/2*
                                (test_gev$loglik_d+test$loglik_d+grad_prior_1)*hess,
                              sd=delta[1]*sqrt(hess), log=TRUE))


  log_posteriors[sim] <-
    sum(test_2$loglik)+ sum(test_gev_2$loglik) +
    hpp_likelihood_2(hyper_gp_latent_pp[sim+1,1], hyper_gp_latent_pp[sim+1,2],
                     x_pp[sim+1,], neigh, intercept=0)

  log_posteriors_2[sim] <- sum(test_gev_2$loglik)

  alpha <-
    # PROPOSALS
    prop_density_1 - prop_density_2 +
    log_posteriors[sim] -
    # LIKELIHOOD
    sum(test$loglik) - sum(test_gev$loglik) -
    # PRIORS
    hpp_likelihood_2(hyper_gp_latent_pp[sim+1,1], hyper_gp_latent_pp[sim+1,2],
                     x_pp[sim,], neigh, intercept=0 )

  # for -Inf likelihood. just stay away from that space
  if (is.nan(alpha)){alpha=-10000}

  u <- runif(1)

  alpha <- log(min(1,exp(alpha)))


  if (log( u ) >alpha){x_pp[sim+1,] <- x_pp[sim,]
  ind_latent_gp[sim_ind,1] <<- 0
  } else {ind_latent_gp[sim_ind,1] <<- 1}

  ###############
  ## SUM TO ZERO CONSTRAINT
  ###############

  x_pp[sim+1,]<- x_pp[sim+1,]-mean(x_pp[sim+1,]) # make sure that sum of x_pp is zero


  ################################################
  ### latent field - NICHE - x_pp_2
  ################################################


  covparms <- c(hyper_gp_latent_pp2[sim+1,1], hyper_gp_latent_pp2[sim+1,2], 0)

  # bin likelihood
  test <- bin_likelihood_d_list(beta_bin[sim,]+beta_share[sim,4]*x_pp_2[sim,]+x_bin[sim+1,], centroids_count_list_bin, fe_bin=fe_bin[sim,])
  #sum(test$loglik)
  #bbs
  test_pp <- poisson_likelihood_d_list(beta_pp_2[sim,]+ x_pp_2[sim,] , centroids_count_list_2,
                                   fe_pp2[sim,])

  #mu
  share_intensity_1 <- get_intensity_list(beta_pp[sim,]+x_pp[sim+1,],size_A, centroids_count_list,
                                          fe_pp[sim,])
  share_intensity <- list()
  for (j in 1:length(share_intensity_1)){
    share_intensity[[j]] <- beta_share[sim,1]*share_intensity_1 [[j]]+beta_share[sim,5]*test$effort[[j]]
  }

  mu_func <- function(y){
    lapply(share_intensity,
           function(x) g_func(1,beta_share[sim,2],x,x_mu[sim,]+beta_mu[sim,]+beta_share[sim,3]*x_pp_2[sim,]+y)+
             0 )
  }

  g_d_func_list <- lapply(share_intensity,
                          function(x) beta_share[sim,3]*g_d_theta2_func(1,beta_share[sim,2],x,x_mu[sim,]+beta_mu[sim,]+beta_share[sim,3]*x_pp_2[sim,]) )

  test_gev <- gev_like_d_theta2_list(data_list=data_list,
                                     mu_func=mu_func,
                                     sigma=exp(beta_sigma[sim,]+x_sigma[sim,]), xi=beta_xi[sim,], constant=0,
                                     fe_mu[sim,], g_d_func_list)


  L_inv <- vecchia_Linv(covparms, "exponential_isotropic", locsord, NNarray)

  grad_prior_1 <- -Linv_mult( L_inv, Linv_t_mult(L_inv, x_pp_2[sim,ord], NNarray) , NNarray)


  like_loop <- (test_pp$loglik_d)
  #chainrule
  chain <- t(weight_route)*exp(beta_pp_2[sim,]+ x_pp_2[sim,] )

  log_d_pp <- numeric()
  for (j in 1:n.sites) {
    log_d_pp[j] <- sum(like_loop * chain[j,] )
    x_pp_2[sim+1,j] <- rnorm(1, mean=x_pp_2[sim,j]+delta[2]^2/2*(beta_share[sim,4]*test$loglik_d[j]+
                                                                   log_d_pp[j]+ test_gev$loglik_d[j]+
                                                                   grad_prior_1[j])*hess[1], sd=delta[2]*sqrt(hess[1]))
  }

  # bin likelihood
  test_2 <- bin_likelihood_d_list(beta_bin[sim,]+beta_share[sim,4]*x_pp_2[sim+1,]+x_bin[sim+1,], centroids_count_list_bin, fe_bin=fe_bin[sim,])
  #test_2$loglik
  #bbs
  test_pp_2 <- poisson_likelihood_d_list(beta_pp_2[sim,]+ x_pp_2[sim+1,], centroids_count_list_2,
                                       fe_pp2[sim,])

  ##mu
  share_intensity <- list()
  for (j in 1:length(share_intensity_1)){
    share_intensity[[j]] <- beta_share[sim,1]*share_intensity_1[[j]]+beta_share[sim,5]*test_2$effort[[j]]
  }


  mu_func <- function(y){
    lapply(share_intensity,
           function(x) g_func(1,beta_share[sim,2],x,x_mu[sim,]+beta_mu[sim,]+beta_share[sim,3]*x_pp_2[sim+1,]+y)+
             0 )
  }

  g_d_func_list <- lapply(share_intensity,
                          function(x) beta_share[sim,3]*g_d_theta2_func(1,beta_share[sim,2],x,x_mu[sim,]+beta_mu[sim,]+beta_share[sim,3]*x_pp_2[sim+1,]) )

  test_gev_2 <- gev_like_d_theta2_list(data_list=data_list,
                                       mu_func=mu_func,
                                       sigma=exp(beta_sigma[sim,]+x_sigma[sim,]), xi=beta_xi[sim,], constant=0,
                                       fe_mu[sim,], g_d_func_list)


  log_d_pp2 <- numeric()
  like_loop <- (test_pp_2$loglik_d)
  #chainrule
  chain <- t(weight_route)*exp(beta_pp_2[sim,]+ x_pp_2[sim+1,] )


  for (j in 1:n.sites) {
    log_d_pp2[j] <- sum( like_loop* chain[j,])
  }

  grad_prior_2 <- -Linv_mult( L_inv, Linv_t_mult(L_inv, x_pp_2[sim+1,ord], NNarray) , NNarray)

  prop_density_1 <- sum(dnorm(x_pp_2[sim,], mean=x_pp_2[sim+1,]+delta[2]^2/2*(beta_share[sim,4]*test_2$loglik_d+
                                                                                log_d_pp2+test_gev_2$loglik_d+grad_prior_2)*hess, sd=delta[2]*sqrt(hess), log=TRUE))
  prop_density_2 <- sum(dnorm(x_pp_2[sim+1,], mean=x_pp_2[sim,]+delta[2]^2/2*(beta_share[sim,4]*test$loglik_d+
                                                                                log_d_pp+test_gev$loglik_d+grad_prior_1)*hess, sd=delta[2]*sqrt(hess), log=TRUE))

  alpha <-
    # PROPOSALS
    prop_density_1 - prop_density_2 +
    sum(test_2$loglik)+ sum(test_pp_2$loglik) + sum(test_gev_2$loglik) +
    hpp_likelihood_2(hyper_gp_latent_pp2[sim+1,1], hyper_gp_latent_pp2[sim+1,2],
                     x_pp_2[sim+1,], neigh, intercept=0) -
    # LIKELIHOOD
    sum(test$loglik) - sum(test_pp$loglik) - sum(test_gev$loglik) -
    # PRIORS
    hpp_likelihood_2(hyper_gp_latent_pp2[sim+1,1], hyper_gp_latent_pp2[sim+1,2],
                     x_pp_2[sim,], neigh, intercept=0 )

  # for -Inf likelihood. just stay away from that space
  if (is.nan(alpha)){alpha=-10000}

  u <- runif(1)


  alpha <- log(min(1,exp(alpha)))

  if (log( u ) >alpha){x_pp_2[sim+1,] <- x_pp_2[sim,]
  ind_latent_gp[sim_ind,2] <<- 0
  } else {ind_latent_gp[sim_ind,2] <<- 1}

  ################

  ###############
  ## SUM TO ZERO CONSTRAINT
  ###############

  x_pp_2[sim+1,]<- x_pp_2[sim+1,]-mean(x_pp_2[sim+1,]) # make sure that sum of x_pp is zero

  ################################################
  ### latent field - GEV - x_mu
  ################################################

  covparms <- c(hyper_gp_latent_mu[sim+1,1], hyper_gp_latent_mu[sim+1,2], 0)

  #mu
  share_intensity_test <- get_intensity_list(beta_pp[sim,]+x_pp[sim+1,],size_A, centroids_count_list,
                                        fe_pp[sim,])
  test_bin <- get_effort_list(beta_bin[sim,]+beta_share[sim,4]*x_pp_2[sim+1,]+x_bin[sim+1,], centroids_count_list_bin,
                              fe_bin=fe_bin[sim,])

  share_intensity <- list()
  for (j in 1:length(share_intensity_test)){
    share_intensity[[j]] <- beta_share[sim,1]*share_intensity_test[[j]]+beta_share[sim,5]*test_bin$effort[[j]]
  }


  mu_func <- function(y){
    lapply(share_intensity,
           function(x) g_func(1,beta_share[sim,2],x,x_mu[sim,]+beta_mu[sim,]+beta_share[sim,3]*x_pp_2[sim+1,]+y)+
             0 )
  }

  g_d_func_list <- lapply(share_intensity,
                          function(x) g_d_theta2_func(1,beta_share[sim,2],x,x_mu[sim,]+beta_mu[sim,]+beta_share[sim,3]*x_pp_2[sim+1,]) )

  test_gev <- gev_like_d_theta2_list(data_list=data_list,
                                     mu_func=mu_func,
                                     sigma=exp(beta_sigma[sim,]+x_sigma[sim,]), xi=beta_xi[sim,], constant=0,
                                     fe_mu[sim,], g_d_func_list)

  L_inv <- vecchia_Linv(covparms, "exponential_isotropic", locsord, NNarray)

  grad_prior_1 <- -Linv_mult( L_inv, Linv_t_mult(L_inv, x_mu[sim,ord], NNarray) , NNarray)

  for (j in 1:n.sites) {
    x_mu[sim+1,j] <- rnorm(1, mean=x_mu[sim,j]+delta[3]^2/2*( test_gev$loglik_d[j]+grad_prior_1[j])*hess[1], sd=delta[3]*sqrt(hess[1]))
  }

  if (any(is.nan(x_mu[sim+1,]) )){

    x_mu[sim+1,] <- x_mu[sim,]
    ind_latent_gp[sim_ind,3] <<- 0
  } else {

    mu_func <- function(y){
      lapply(share_intensity,
             function(x) g_func(1,beta_share[sim,2],x,x_mu[sim+1,]+beta_mu[sim,]+beta_share[sim,3]*x_pp_2[sim+1,]+y)+
               0 )
    }

    g_d_func_list <- lapply(share_intensity,
                            function(x) g_d_theta2_func(1,beta_share[sim,2],x,x_mu[sim+1,]+beta_mu[sim,]+beta_share[sim,3]*x_pp_2[sim+1,]) )

    test_gev_2 <- gev_like_d_theta2_list(data_list=data_list,
                                       mu_func=mu_func,
                                       sigma=exp(beta_sigma[sim,]+x_sigma[sim,]), xi=beta_xi[sim,], constant=0,
                                       fe_mu[sim,], g_d_func_list)

    grad_prior_2 <- -Linv_mult( L_inv, Linv_t_mult(L_inv, x_mu[sim+1,ord], NNarray) , NNarray)


    prop_density_1 <- sum(dnorm(x_mu[sim,], mean=x_mu[sim+1,]+delta[3]^2/2*(test_gev_2$loglik_d+grad_prior_2)*hess, sd=delta[3]*sqrt(hess), log=TRUE))
    prop_density_2 <- sum(dnorm(x_mu[sim+1,], mean=x_mu[sim,]+delta[3]^2/2*(test_gev$loglik_d+grad_prior_1)*hess, sd=delta[3]*sqrt(hess), log=TRUE))



    alpha <-
      # PROPOSALS
      prop_density_1 - prop_density_2 +
      sum(test_gev_2$loglik) +
      hpp_likelihood_2(hyper_gp_latent_mu[sim+1,1], hyper_gp_latent_mu[sim+1,2],
                       x_mu[sim+1,], neigh, intercept=0) -
      # LIKELIHOOD
      sum(test_gev$loglik) -
      # PRIORS
      hpp_likelihood_2(hyper_gp_latent_mu[sim+1,1], hyper_gp_latent_mu[sim+1,2],
                       x_mu[sim,], neigh, intercept=0 )

    # for -Inf likelihood. just stay away from that space
    if (is.nan(alpha)){alpha=-10000}

    u <- runif(1)

    alpha <- log(min(1,exp(alpha)))


    if (log( u ) >alpha){x_mu[sim+1,] <- x_mu[sim,]
    ind_latent_gp[sim_ind,3] <<- 0
    } else {ind_latent_gp[sim_ind,3] <<- 1}

  }
  ###############
  ## SUM TO ZERO CONSTRAINT
  ###############

  x_mu[sim+1,]<- x_mu[sim+1,]-mean(x_mu[sim+1,]) # make sure that sum of x_pp is zero

  ################################################
  ### latent field - GEV - x_sigma
  ################################################

  covparms <- c(hyper_gp_latent_sigma[sim+1,1], hyper_gp_latent_sigma[sim+1,2], 0)

  mu_func <- function(y){
    lapply(share_intensity,
           function(x) g_func(1,beta_share[sim,2],x,x_mu[sim+1,]+beta_mu[sim,]+beta_share[sim,3]*x_pp_2[sim+1,]+y)+
             0 )
  }

  test_gev <- gev_like_dsigma_list(data=data_list,
                                   mu_func=mu_func,
                                   sigma=exp(beta_sigma[sim,]+x_sigma[sim,]), xi=beta_xi[sim,], constant=0,
                                   fe_mu[sim,])

  L_inv <- vecchia_Linv(covparms, "exponential_isotropic", locsord, NNarray)

  grad_prior_1 <- -Linv_mult( L_inv, Linv_t_mult(L_inv, x_sigma[sim,ord], NNarray) , NNarray)

  ##CHAIN RULE dont forget
  for (j in 1:n.sites) {
    x_sigma[sim+1,j] <- rnorm(1, mean=x_sigma[sim,j]+delta[5]^2/2*(exp(beta_sigma[sim,]+x_sigma[sim,j])*test_gev$loglik_d[j]+
                                                                     grad_prior_1[j])*hess[1], sd=delta[5]*sqrt(hess[1]))
  }

  if (any(is.nan(x_sigma[sim+1,]) )){

    x_sigma[sim+1,] <- x_sigma[sim,]
    ind_latent_gp[sim_ind,5] <<- 0
  } else {


    mu_func <- function(y){
      lapply(share_intensity,
             function(x) g_func(1,beta_share[sim,2],x,x_mu[sim+1,]+beta_mu[sim,]+beta_share[sim,3]*x_pp_2[sim+1,]+y)+
               0 )
    }

    test_gev_2 <- gev_like_dsigma_list(data=data_list,
                                     mu_func=mu_func,
                                     sigma=exp(beta_sigma[sim,]+x_sigma[sim+1,]), xi=beta_xi[sim,], constant=0,
                                     fe_mu[sim,])

    grad_prior_2 <- -Linv_mult( L_inv, Linv_t_mult(L_inv, x_sigma[sim+1,ord], NNarray) , NNarray)


    prop_density_1 <- sum(dnorm(x_sigma[sim,], mean=x_sigma[sim+1,]+delta[5]^2/2*(exp(beta_sigma[sim,]+x_sigma[sim+1,])*test_gev_2$loglik_d+grad_prior_2)*hess, sd=delta[5]*sqrt(hess), log=TRUE))
    prop_density_2 <- sum(dnorm(x_sigma[sim+1,], mean=x_sigma[sim,]+delta[5]^2/2*(exp(beta_sigma[sim,]+x_sigma[sim,])*test_gev$loglik_d+grad_prior_1)*hess, sd=delta[5]*sqrt(hess), log=TRUE))



    alpha <-
      # PROPOSALS
      prop_density_1 - prop_density_2 +
      sum(test_gev_2$loglik) +
      hpp_likelihood_2(hyper_gp_latent_sigma[sim+1,1], hyper_gp_latent_sigma[sim+1,2],
                       x_sigma[sim+1,], neigh, intercept=0) -
      # LIKELIHOOD
      sum(test_gev$loglik) -
      # PRIORS
      hpp_likelihood_2(hyper_gp_latent_sigma[sim+1,1], hyper_gp_latent_sigma[sim+1,2],
                       x_sigma[sim,], neigh, intercept=0 )

    # for -Inf likelihood. just stay away from that space
    if (is.nan(alpha)){alpha=-10000}

    u <- runif(1)

    alpha <- log(min(1,exp(alpha)))


    if (log( u ) >alpha){x_sigma[sim+1,] <- x_sigma[sim,]
    ind_latent_gp[sim_ind,5] <<- 0
    } else {ind_latent_gp[sim_ind,5] <<- 1}

  }
  ###############
  ## SUM TO ZERO CONSTRAINT
  ###############

  x_sigma[sim+1,]<- x_sigma[sim+1,]-mean(x_sigma[sim+1,]) # make sure that sum of x_pp is zero


  #x_sigma[sim+1,]<- 0

  ################################################
  ### intercept mu
  ################################################

  beta_mu[sim+1,] <- rnorm(n=1, mean=beta_mu[sim,], sd=stepsize_intercept[1])

  mu_func <- function(y){
    lapply(share_intensity,
           function(x) g_func(1,beta_share[sim,2],x,x_mu[sim+1,]+beta_mu[sim,]+beta_share[sim,3]*x_pp_2[sim+1,]+y)+
             0 )
  }

  test_gev_intercept <- gev_like_list(data=data_list,
                                        mu_func=mu_func,
                                        sigma=exp(beta_sigma[sim,]+x_sigma[sim+1,]),
                                        xi=beta_xi[sim,], constant=0,
                                        fe_mu[sim,])


  mu_func <- function(y){
    lapply(share_intensity,
           function(x) g_func(1,beta_share[sim,2],x,x_mu[sim+1,]+beta_mu[sim+1,]+beta_share[sim,3]*x_pp_2[sim+1,]+y)+
             0 )
  }

  test_gev_2_intercept <- gev_like_list(data=data_list,
                                        mu_func=mu_func,
                                        sigma=exp(beta_sigma[sim,]+x_sigma[sim+1,]),
                                        xi=beta_xi[sim,], constant=0,
                                        fe_mu[sim,])


  alpha <-
    sum(test_gev_2_intercept) -
    sum(test_gev_intercept) +
    prior_normal(x=beta_mu[sim+1,], mu=mu, sigma=sigma ) -
    prior_normal(x=beta_mu[sim,], mu=mu, sigma=sigma)

  # for -Inf likelihood. just stay away from that space
  if (is.nan(alpha)){alpha=-10000}
  u <- runif(1)
  alpha <- log(min(1,exp(alpha)))

  if (log( u ) >alpha){
    ind_intercept[sim_ind,1] <<- 0
    beta_mu[sim+1,] <- beta_mu[sim,] } else {ind_intercept[sim_ind,1] <<- 1}


  # sigma

  beta_sigma[sim+1,] <- rnorm(n=1, mean=beta_sigma[sim,], sd=stepsize_intercept[2])


  mu_func <- function(y){
    lapply(share_intensity,
           function(x) g_func(1,beta_share[sim,2],x,x_mu[sim+1,]+beta_mu[sim+1,]+beta_share[sim,3]*x_pp_2[sim+1,]+y)+
             0 )
  }

  test_gev_intercept <- gev_like_list(data=data_list,
                                      mu_func=mu_func,
                                      sigma=exp(beta_sigma[sim,]+x_sigma[sim+1,]),
                                      xi=beta_xi[sim,], constant=0,
                                      fe_mu[sim,])

  test_gev_2_intercept <- gev_like_list(data=data_list,
                                        mu_func=mu_func,
                                        sigma=exp(beta_sigma[sim+1,]+x_sigma[sim+1,]),
                                        xi=beta_xi[sim,], constant=0,
                                        fe_mu[sim,])


  alpha <-
    sum(test_gev_2_intercept) -
    sum(test_gev_intercept) +
    prior_normal(x=beta_sigma[sim+1,], mu=mu, sigma=sigma) -
    prior_normal(x=beta_sigma[sim,], mu=mu, sigma=sigma)

  # for -Inf likelihood. just stay away from that space
  if (is.nan(alpha)){alpha=-10000}

  u <- runif(1)
  alpha <- log(min(1,exp(alpha)))

  if (log( u ) >alpha){
    ind_intercept[sim_ind,2] <<- 0
    beta_sigma[sim+1,] <- beta_sigma[sim,] } else {ind_intercept[sim_ind,2] <<- 1}




  # xi

  beta_xi[sim+1,] <- rnorm(n=1, mean=beta_xi[sim,], sd=stepsize_intercept[3])

  test_gev_intercept <- gev_like_list(data=data_list,
                                      mu_func=mu_func,
                                      sigma=exp(beta_sigma[sim,]+x_sigma[sim+1,]),
                                      xi=beta_xi[sim,], constant=0,
                                      fe_mu[sim,])

  test_gev_2_intercept <- gev_like_list(data=data_list,
                                        mu_func=mu_func,
                                        sigma=exp(beta_sigma[sim+1,]+x_sigma[sim+1,]),
                                        xi=beta_xi[sim+1,], constant=0,
                                        fe_mu[sim,])


  alpha <-
    sum(test_gev_2_intercept) -
    sum(test_gev_intercept) +
    prior_normal(x=beta_xi[sim+1,], mu=mu, sigma=sigma) -
    prior_normal(x=beta_xi[sim,], mu=mu, sigma=sigma)

  # for -Inf likelihood. just stay away from that space
  if (is.nan(alpha)){alpha=-10000}

  u <- runif(1)
  alpha <- log(min(1,exp(alpha)))

  if (log( u ) >alpha){
    ind_intercept[sim_ind,3] <<- 0
    beta_xi[sim+1,] <- beta_xi[sim,] } else {ind_intercept[sim_ind,3] <<- 1}



  ### SHARING 1

  beta_share[sim+1,1] <- rnorm(1,beta_share[sim,1], stepsize_share[1])

  mu_func <- function(y){
    lapply(share_intensity,
           function(x) g_func(1,beta_share[sim,2],x,x_mu[sim+1,]+beta_mu[sim+1,]+beta_share[sim,3]*x_pp_2[sim+1,]+y)+
             0 )
  }

  test_gev_intercept <- gev_like_list(data=data_list,
                                      mu_func=mu_func,
                                      sigma=exp(beta_sigma[sim+1,]+x_sigma[sim+1,]),
                                      xi=beta_xi[sim+1,], constant=0,
                                      fe_mu[sim,])


  share_intensity <- list()
  for (j in 1:length(share_intensity_test)){
    share_intensity[[j]] <- beta_share[sim+1,1]*share_intensity_test[[j]]+beta_share[sim,5]*test_bin$effort[[j]]
  }


  mu_func <- function(y){
    lapply(share_intensity,
           function(x) g_func(1,beta_share[sim,2],x,x_mu[sim+1,]+beta_mu[sim+1,]+beta_share[sim,3]*x_pp_2[sim+1,]+y)+
             0 )
  }

  test_gev_2_intercept <- gev_like_list(data=data_list,
                                        mu_func=mu_func,
                                        sigma=exp(beta_sigma[sim+1,]+x_sigma[sim+1,]),
                                        xi=beta_xi[sim+1,], constant=0,
                                        fe_mu[sim,])

  alpha <- dnorm(beta_share[sim+1,1], 0, 100, log = T) -
    dnorm(beta_share[sim,1], 0, 100, log = T) +
    sum(test_gev_2_intercept) - sum(test_gev_intercept)

  # for -Inf likelihood. just stay away from that space
  if (is.nan(alpha)){alpha=-10000}

  u <- runif(1)

  alpha <- log(min(1,exp(alpha)))

  if (log( u ) >alpha){
    ind_share[sim_ind,1] <<- 0
    beta_share[sim+1,1] <- beta_share[sim,1] } else {ind_share[sim_ind,1] <<- 1}



  ### SHARING 2

  
  beta_share[sim+1,2] <- rnorm(1,beta_share[sim,2], stepsize_share[2])
  
  share_intensity <- list()
  for (j in 1:length(share_intensity_test)){
    share_intensity[[j]] <- beta_share[sim+1,1]*share_intensity_test[[j]]+beta_share[sim,5]*test_bin$effort[[j]]
  }
  
  mu_func <- function(y){
    lapply(share_intensity,
           function(x) g_func(1,beta_share[sim,2],x,x_mu[sim+1,]+beta_mu[sim+1,]+beta_share[sim,3]*x_pp_2[sim+1,]+y)+
             0 )
  }
  
  test_gev_intercept <- gev_like_list(data=data_list,
                                      mu_func=mu_func,
                                      sigma=exp(beta_sigma[sim+1,]+x_sigma[sim+1,]),
                                      xi=beta_xi[sim+1,], constant=0,
                                      fe_mu[sim,])
  
  
  mu_func <- function(y){
    lapply(share_intensity,
           function(x) g_func(1,beta_share[sim+1,2],x,x_mu[sim+1,]+beta_mu[sim+1,]+beta_share[sim,3]*x_pp_2[sim+1,]+y)+
             0 )
  }
  
  test_gev_2_intercept <- gev_like_list(data=data_list,
                                        mu_func=mu_func,
                                        sigma=exp(beta_sigma[sim+1,]+x_sigma[sim+1,]),
                                        xi=beta_xi[sim+1,], constant=0,
                                        fe_mu[sim,])
  
  alpha <- dnorm(beta_share[sim+1,2], 0, 100, log = T) -
    dnorm(beta_share[sim,2], 0, 100, log = T) +
    sum(test_gev_2_intercept) - sum(test_gev_intercept)
  
  # for -Inf likelihood. just stay away from that space
  if (is.nan(alpha)){alpha=-10000}
  
  u <- runif(1)
  
  alpha <- log(min(1,exp(alpha)))
  
  if (log( u ) >alpha){
    ind_share[sim_ind,2] <<- 0
    beta_share[sim+1,2] <- beta_share[sim,2] } else {ind_share[sim_ind,2] <<- 1}
  


  ### intercept COX

  # COX

  beta_pp[sim+1,] <- rnorm(n=1, mean=beta_pp[sim,], sd=stepsize_intercept[4])

  test_pp_intercept <- hpp_likelihood_d_list(beta_pp[sim,]+x_pp[sim+1,],size_A, centroids_count_list,
                                fe_pp[sim,])

  share_intensity_1 <- test_pp_intercept$lambda

  test_bin <- get_effort_list(beta_bin[sim,]+beta_share[sim,4]*x_pp_2[sim+1,]+x_bin[sim+1,], centroids_count_list_bin,
                              fe_bin=fe_bin[sim,])
  share_intensity <- list()
  for (j in 1:length(share_intensity_1)){
    share_intensity[[j]] <- beta_share[sim+1,1]*share_intensity_1[[j]]+beta_share[sim,5]*test_bin$effort[[j]]
  }

  #mu
  mu_func <- function(y){
    lapply(share_intensity,
           function(x) g_func(1,beta_share[sim+1,2],x,x_mu[sim+1,]+beta_mu[sim+1,]+
                                beta_share[sim,3]*x_pp_2[sim+1,]+y)+
             0 )
  }

  test_gev <- gev_like_list(data=data_list,
                  mu_func=mu_func,
                  sigma=exp(beta_sigma[sim+1,]+x_sigma[sim+1,]),
                  xi=beta_xi[sim+1,], constant=0,
                  fe_mu[sim,])


  test_pp_2_intercept <- hpp_likelihood_d_list(beta_pp[sim+1,]+x_pp[sim+1,],size_A, centroids_count_list,
                                           fe_pp[sim,])

  share_intensity_2 <- test_pp_2_intercept$lambda

  share_intensity <- list()
  for (j in 1:length(share_intensity_2)){
    share_intensity[[j]] <- beta_share[sim+1,1]*share_intensity_2[[j]]+beta_share[sim,5]*test_bin$effort[[j]]
  }

  #mu
  mu_func <- function(y){
    lapply(share_intensity,
           function(x) g_func(1,beta_share[sim+1,2],x,x_mu[sim+1,]+beta_mu[sim+1,]+
                                beta_share[sim,3]*x_pp_2[sim+1,]+y)+
             0 )
  }

  test_gev_2 <- gev_like_list(data=data_list,
                            mu_func=mu_func,
                            sigma=exp(beta_sigma[sim+1,]+x_sigma[sim+1,]),
                            xi=beta_xi[sim+1,], constant=0,
                            fe_mu[sim,])

  alpha <-
    sum(unlist(test_pp_2_intercept$loglik)) +sum(test_gev_2) -
    sum(unlist(test_pp_intercept$loglik)) - sum(test_gev) +
    prior_normal(x=beta_pp[sim+1,], mu=mu, sigma=sigma) -
    prior_normal(x=beta_pp[sim,], mu=mu, sigma=sigma)

  # for -Inf likelihood. just stay away from that space
  if (is.nan(alpha)){alpha=-10000}

  u <- runif(1)
  alpha <- log(min(1,exp(alpha)))

  if (log( u ) >alpha){
    ind_intercept[sim_ind,4] <<- 0
    beta_pp[sim+1,] <- beta_pp[sim,] } else {ind_intercept[sim_ind,4] <<- 1}




  ### intercept BIN

  # affects both the BIN and LOC components

  beta_bin[sim+1,] <- rnorm(n=1, mean=beta_bin[sim,], sd=stepsize_intercept[5])

  test_bin <- bin_likelihood_d_list(beta_bin[sim,]+beta_share[sim,4]*x_pp_2[sim+1,]+x_bin[sim+1,], centroids_count_list_bin, fe_bin=fe_bin[sim,])

  # share_intensity <- get_intensity_list(beta_pp[sim+1,]+x_pp[sim+1,],size_A, centroids_count_list,
  #                                       fe_pp[sim,])
  # combined_intensity <- list()
  # for (j in 1:length(share_intensity)){
  #   combined_intensity[[j]] <- beta_share[sim+1,1]*share_intensity[[j]]+beta_share[sim,5]*test_bin$effort[[j]]
  # }
  #
  # #mu
  # mu_func <- function(y){
  #   lapply(combined_intensity,
  #          function(x) g_func(1,beta_share[sim+1,2],x,x_mu[sim+1,]+beta_mu[sim+1,]+
  #                               beta_share[sim,3]*x_pp_2[sim+1,]+y)+
  #            0 )
  # }
  #
  # test_gev <- gev_like_list(data=data_list,
  #                           mu_func=mu_func,
  #                           sigma=exp(beta_sigma[sim+1,]+x_sigma[sim+1,]),
  #                           xi=beta_xi[sim+1,], constant=0,
  #                           fe_mu[sim,])


  test_bin_2 <- bin_likelihood_d_list(beta_bin[sim+1,]+beta_share[sim,4]*x_pp_2[sim+1,]+x_bin[sim+1,], centroids_count_list_bin, fe_bin=fe_bin[sim,])

  # combined_intensity <- list()
  # for (j in 1:length(share_intensity)){
  #   combined_intensity[[j]] <- beta_share[sim+1,1]*share_intensity[[j]]+beta_share[sim,5]*test_bin_2$effort[[j]]
  # }
  #
  # #mu
  # mu_func <- function(y){
  #   lapply(combined_intensity,
  #          function(x) g_func(1,beta_share[sim+1,2],x,x_mu[sim+1,]+beta_mu[sim+1,]+
  #                                 beta_share[sim,3]*x_pp_2[sim+1,]+y)+
  #            0 )
  # }
  #
  # test_gev_2 <- gev_like_list(data=data_list,
  #                           mu_func=mu_func,
  #                           sigma=exp(beta_sigma[sim+1,]+x_sigma[sim+1,]),
  #                           xi=beta_xi[sim+1,], constant=0,
  #                           fe_mu[sim,])

  log_posteriors_3[sim] <- sum(test_bin$loglik)


  alpha <-
    sum(test_bin_2$loglik) -#+sum(test_gev_2) -
    sum(test_bin$loglik) + #- sum(test_gev) +
    prior_normal(x=beta_bin[sim+1,], mu=mu, sigma=sigma) -
    prior_normal(x=beta_bin[sim,], mu=mu, sigma=sigma)

  # for -Inf likelihood. just stay away from that space
  if (is.nan(alpha)){alpha=-10000}

  u <- runif(1)
  alpha <- log(min(1,exp(alpha)))

  if (log( u ) >alpha){
    ind_intercept[sim_ind,5] <<- 0
    beta_bin[sim+1,] <- beta_bin[sim,] } else {ind_intercept[sim_ind,5] <<- 1}





  ### SHARING 3

  #beta_share[sim+1,3] <- 0

  beta_share[sim+1,3] <- rnorm(1,beta_share[sim,3], stepsize_share[3])

  share_intensity_test <- get_intensity_list(beta_pp[sim+1,]+x_pp[sim+1,],size_A, centroids_count_list,
                                        fe_pp[sim,])
  test_bin <- bin_likelihood_d_list(beta_bin[sim+1,]+x_pp_2[sim+1,]+x_bin[sim+1,], centroids_count_list_bin, fe_bin=fe_bin[sim,])

  share_intensity <- list()
  for (j in 1:length(share_intensity_test)){
    share_intensity[[j]] <- beta_share[sim+1,1]*share_intensity_test[[j]]+beta_share[sim,5]*test_bin$effort[[j]]
  }

  #mu
  mu_func <- function(y){
    lapply(share_intensity,
           function(x) g_func(1,beta_share[sim+1,2],x,x_mu[sim+1,]+beta_mu[sim+1,]+
                                beta_share[sim,3]*x_pp_2[sim+1,]+y)+
             0 )
  }

  test_gev <- gev_like_list(data=data_list,
                              mu_func=mu_func,
                              sigma=exp(beta_sigma[sim+1,]+x_sigma[sim+1,]),
                              xi=beta_xi[sim+1,], constant=0,
                              fe_mu[sim,])


  mu_func <- function(y){
    lapply(share_intensity,
           function(x) g_func(1,beta_share[sim+1,2],x,x_mu[sim+1,]+beta_mu[sim+1,]+
                                beta_share[sim+1,3]*x_pp_2[sim+1,]+y)+
             0 )
  }

  test_gev_2 <- gev_like_list(data=data_list,
                            mu_func=mu_func,
                            sigma=exp(beta_sigma[sim+1,]+x_sigma[sim+1,]),
                            xi=beta_xi[sim+1,], constant=0,
                            fe_mu[sim,])

  alpha <-
    sum(test_gev_2) -
   sum(test_gev) + dnorm(beta_share[sim+1,3], 0, 100, log = T) -
    dnorm(beta_share[sim,3], 0, 100, log = T)

    # for -Inf likelihood. just stay away from that space
    if (is.nan(alpha)){alpha=-10000}

  u <- runif(1)
  alpha <- log(min(1,exp(alpha)))

  if (log( u ) >alpha){
    ind_share[sim_ind,3] <<- 0
    beta_share[sim+1,3] <- beta_share[sim,3] } else {ind_share[sim_ind,3] <<- 1}



  ##########
  ## COX 2


  beta_pp_2[sim+1,] <- rnorm(n=1, mean=beta_pp_2[sim,], sd=stepsize_intercept[6])

  test_pp_intercept <- poisson_likelihood_list(beta_pp_2[sim,]+ x_pp_2[sim+1,],
                                           centroids_count_list_2,fe_pp2[sim,])
  test_pp_2_intercept <- poisson_likelihood_list(beta_pp_2[sim+1,]+ x_pp_2[sim+1,] ,
                                             centroids_count_list_2,fe_pp2[sim,])

  log_posteriors_4[sim] <- sum(unlist(test_pp_intercept))


  alpha <-
    sum(unlist(test_pp_2_intercept)) -
    sum(unlist(test_pp_intercept)) +
    prior_normal(x=beta_pp_2[sim+1,], mu=mu, sigma=sigma) -
    prior_normal(x=beta_pp_2[sim,], mu=mu, sigma=sigma)

  # for -Inf likelihood. just stay away from that space
  if (is.nan(alpha)){alpha=-10000}

  u <- runif(1)
  alpha <- log(min(1,exp(alpha)))

  if (log( u ) >alpha){
    ind_intercept[sim_ind,6] <<- 0
    beta_pp_2[sim+1,] <- beta_pp_2[sim,] } else {ind_intercept[sim_ind,6] <<- 1}



  ## fixed effects

  ## pp

  #if (!is.null(fe_pp)){
    covparms <- c(hyper_gp_latent_re[sim+1,1], hyper_gp_latent_re[sim+1,2], 0)

    test <- hpp_likelihood_d_list_time(beta_pp[sim+1,]+x_pp[sim+1,],size_A, centroids_count_list,
                                  fe_pp[sim,])

    share_intensity <- test$lambda

    test_bin <- get_effort_list(beta_bin[sim+1,]+beta_share[sim,4]*x_pp_2[sim+1,]+x_bin[sim+1,], centroids_count_list_bin,
                                fe_bin=fe_bin[sim,])

    share_intensity_1 <- list()
    for (j in 1:length(share_intensity)){
      share_intensity_1[[j]] <- beta_share[sim+1,1]*share_intensity[[j]]+beta_share[sim,5]*test_bin$effort[[j]]
    }

    #mu
    mu_func <- function(y){
      lapply(share_intensity_1,
             function(x) g_func(1,beta_share[sim+1,2],x,x_mu[sim+1,]+beta_mu[sim+1,]+beta_share[sim+1,3]*x_pp_2[sim+1,]+y)+
               0 )
    }

    g_d_func_list <- lapply(share_intensity_1,
                            function(x) g_d_func(1,beta_share[sim+1,2],x,x_mu[sim+1,]+beta_mu[sim+1,]+beta_share[sim+1,3]*x_pp_2[sim+1,]) )


    test_gev <- gev_like_d_list_time(data_list=data_list,
                                mu_func=mu_func,
                                sigma=exp(beta_sigma[sim+1,]+x_sigma[sim+1,]), xi=beta_xi[sim+1,], constant=0,
                                fe_mu[sim,], g_d_func_list)

    #sum(test_gev$loglik)
    L_inv <- vecchia_Linv(covparms, "exponential_isotropic", matrix(locsord_time,ncol=1), NNarray_time)

    grad_prior_1 <- -Linv_mult( L_inv, Linv_t_mult(L_inv, fe_pp[sim,ord_time], NNarray_time) , NNarray_time)


    ## DONT FORGET THE CHAIN RULE
    for (j in 1:n.year) {
      fe_pp[sim+1,j] <- rnorm(1, mean=fe_pp[sim,j]+stepsize_fe[1]^2/2*
                               (test_gev$loglik_d[j]+test$loglik_d[j]+grad_prior_1[j])*hess[1], sd=stepsize_fe[1]*sqrt(hess[1]))
    }

    test_2 <- hpp_likelihood_d_list_time(beta_pp[sim+1,]+x_pp[sim+1,],size_A, centroids_count_list,
                                    fe_pp[sim+1,])
    share_intensity <- test_2$lambda


    share_intensity_2 <- list()
    for (j in 1:length(share_intensity)){
      share_intensity_2[[j]] <- beta_share[sim+1,1]*share_intensity[[j]]+beta_share[sim,5]*test_bin$effort[[j]]
    }


    mu_func <- function(y){
      lapply(share_intensity_2,
             function(x) g_func(1,beta_share[sim+1,2],x,x_mu[sim+1,]+beta_mu[sim+1,]+
                                  beta_share[sim+1,3]*x_pp_2[sim+1,]+y)+
               0 )
    }

    g_d_func_list <- lapply(share_intensity_2,
                            function(x) g_d_func(1,beta_share[sim+1,2],x,x_mu[sim+1,]+beta_mu[sim+1,]+beta_share[sim+1,3]*x_pp_2[sim+1,]) )


    test_gev_2 <- gev_like_d_list_time(data=data_list,
                                  mu_func=mu_func,
                                  sigma=exp(beta_sigma[sim+1,]+x_sigma[sim+1,]), xi=beta_xi[sim+1,], constant=0,
                                  fe_mu[sim,], g_d_func_list)

    grad_prior_2 <- -Linv_mult( L_inv, Linv_t_mult(L_inv, fe_pp[sim+1,ord_time], NNarray_time) , NNarray_time)


    prop_density_1 <- sum(dnorm(fe_pp[sim,],
                                mean=fe_pp[sim+1,]+stepsize_fe[1]^2/2*
                                  (test_gev_2$loglik_d+test_2$loglik_d+grad_prior_2)*hess,
                                sd=stepsize_fe[1]*sqrt(hess), log=TRUE))
    prop_density_2 <- sum(dnorm(fe_pp[sim+1,],
                                mean=fe_pp[sim,]+stepsize_fe[1]^2/2*
                                  (test_gev$loglik_d+test$loglik_d+grad_prior_1)*hess,
                                sd=stepsize_fe[1]*sqrt(hess), log=TRUE))

    alpha <-
      # PROPOSALS
      prop_density_1 - prop_density_2 +
      sum(test_2$loglik)+ sum(test_gev_2$loglik) +
      hpp_likelihood_2_time(hyper_gp_latent_re[sim+1,1], hyper_gp_latent_re[sim+1,2],
                       fe_pp[sim+1,], neigh, intercept=0) -
      # LIKELIHOOD
      sum(test$loglik) - sum(test_gev$loglik) -
      # PRIORS
      hpp_likelihood_2_time(hyper_gp_latent_re[sim+1,1], hyper_gp_latent_re[sim+1,2],
                       fe_pp[sim,], neigh, intercept=0 )

    # for -Inf likelihood. just stay away from that space
    if (is.nan(alpha)){alpha=-10000}

    u <- runif(1)

    alpha <- log(min(1,exp(alpha)))

    if (log( u ) >alpha){
      ind_fe[sim_ind,1] <<- 0
      fe_pp[sim+1,] <- fe_pp[sim,] } else {ind_fe[sim_ind,1] <<- 1}


    ###############
    ## SUM TO ZERO CONSTRAINT
    ###############

    fe_pp[sim+1,]<- fe_pp[sim+1,]-mean(fe_pp[sim+1,]) # make sure that sum of x_pp is zero


 # }


  #effort for duration, which affects both bin and loc


  ## bin

  if (!is.null(fe_bin)){
    fe_bin[sim+1,] <- rnorm(n=1, mean=fe_bin[sim,], sd=stepsize_fe[2])

    test_bin <- bin_likelihood_d_list(beta_bin[sim+1,]+beta_share[sim,4]*x_pp_2[sim+1,]+x_bin[sim+1,], centroids_count_list_bin,
                                      fe_bin=fe_bin[sim,])

    test_bin_2 <- bin_likelihood_d_list(beta_bin[sim+1,]+beta_share[sim,4]*x_pp_2[sim+1,]+x_bin[sim+1,], centroids_count_list_bin,
                                        fe_bin=fe_bin[sim+1,])

    alpha <-
      sum(test_bin_2$loglik) -
      sum(test_bin$loglik) +
      prior_normal(x=fe_bin[sim+1,], mu=mu, sigma=sigma) -
      prior_normal(x=fe_bin[sim,], mu=mu, sigma=sigma)

    # for -Inf likelihood. just stay away from that space
    if (is.nan(alpha)){alpha=-10000}

    u <- runif(1)
    alpha <- log(min(1,exp(alpha)))

    if (log( u ) >alpha){
      ind_fe[sim_ind,2] <<- 0
      fe_bin[sim+1,] <- fe_bin[sim,] } else {ind_fe[sim_ind,2] <<- 1}
  }


  #### FE_MU


  if (!is.null(fe_mu)){
    fe_mu[sim+1,] <- rnorm(n=1, mean=fe_mu[sim,], sd=stepsize_fe[3])

    test_bin <- get_effort_list(beta_bin[sim+1,]+beta_share[sim,4]*x_pp_2[sim+1,]+x_bin[sim+1,], centroids_count_list_bin,
                                      fe_bin=fe_bin[sim+1,])

    combined_intensity <- list()
    for (j in 1:length(share_intensity)){
      combined_intensity[[j]] <- beta_share[sim+1,1]*share_intensity[[j]]+beta_share[sim,5]*test_bin$effort[[j]]
    }

    #mu
    mu_func <- function(y){
      lapply(combined_intensity,
             function(x) g_func(1,beta_share[sim+1,2],x,x_mu[sim+1,]+beta_mu[sim+1,]+
                                  beta_share[sim+1,3]*x_pp_2[sim+1,]+y)+
               0 )
    }

    test_gev <- gev_like_list(data=data_list,
                              mu_func=mu_func,
                              sigma=exp(beta_sigma[sim+1,]+x_sigma[sim+1,]),
                              xi=beta_xi[sim+1,], constant=0,
                              fe_mu[sim,])

    test_gev_2 <- gev_like_list(data=data_list,
                              mu_func=mu_func,
                              sigma=exp(beta_sigma[sim+1,]+x_sigma[sim+1,]),
                              xi=beta_xi[sim+1,], constant=0,
                              fe_mu[sim+1,])

    alpha <-
      sum(test_gev_2) -
     sum(test_gev)  +
      prior_normal(x=fe_mu[sim+1,], mu=mu, sigma=sigma) -
      prior_normal(x=fe_mu[sim,], mu=mu, sigma=sigma)

    # for -Inf likelihood. just stay away from that space
    if (is.nan(alpha)){alpha=-10000}

    u <- runif(1)
    alpha <- log(min(1,exp(alpha)))

    if (log( u ) >alpha){
      ind_fe[sim_ind,3] <<- 0
      fe_mu[sim+1,] <- fe_mu[sim,] } else {ind_fe[sim_ind,3] <<- 1}
  }


    ##########
    ### Beta 4
    ##########
    beta_share[sim+1,4] <- 1

    # beta_share[sim+1,4] <- rnorm(n=1, mean=beta_share[sim,4], sd=stepsize_share[4])
    #
    # test_pp_intercept <- bin_likelihood_list(beta_bin[sim+1,]+beta_share[sim,4]*x_pp_2[sim+1,]+x_bin[sim+1,], centroids_count_list_bin,
    #                                            fe_bin=fe_bin[sim+1,])
    # test_pp_2_intercept <- bin_likelihood_list(beta_bin[sim+1,]+beta_share[sim+1,4]*x_pp_2[sim+1,]+x_bin[sim+1,], centroids_count_list_bin,
    #                                              fe_bin=fe_bin[sim+1,])
    #
    # alpha <-
    #   sum(unlist(test_pp_2_intercept)) -
    #   sum(unlist(test_pp_intercept)) +
    #   prior_normal(x=beta_share[sim+1,4], mu=mu, sigma=sigma) -
    #   prior_normal(x=beta_share[sim,4], mu=mu, sigma=sigma)
    #
    # # for -Inf likelihood. just stay away from that space
    # if (is.nan(alpha)){alpha=-10000}
    #
    # u <- runif(1)
    # alpha <- log(min(1,exp(alpha)))
    #
    # if (log( u ) >alpha){
    #   ind_share[sim_ind,4] <<- 0
    #   beta_share[sim+1,4] <- beta_share[sim,4] } else {ind_share[sim_ind,4] <<- 1}

    ##########
    ### Beta 5
    ##########

    beta_share[sim+1,5] <- rnorm(n=1, mean=beta_share[sim,5], sd=stepsize_share[5])

    test_bin <- get_effort_list(beta_bin[sim+1,]+beta_share[sim,4]*x_pp_2[sim+1,]+x_bin[sim+1,], centroids_count_list_bin,
                                fe_bin=fe_bin[sim+1,])

    combined_intensity <- list()
    for (j in 1:length(share_intensity)){
      combined_intensity[[j]] <- beta_share[sim+1,1]*share_intensity[[j]]+beta_share[sim,5]*test_bin$effort[[j]]
    }

    #mu
    mu_func <- function(y){
      lapply(combined_intensity,
             function(x) g_func(1,beta_share[sim+1,2],x,x_mu[sim+1,]+beta_mu[sim+1,]+
                                  beta_share[sim+1,3]*x_pp_2[sim+1,]+y)+
               0 )
    }

    combined_intensity_2 <- list()
    for (j in 1:length(share_intensity)){
      combined_intensity_2[[j]] <- beta_share[sim+1,1]*share_intensity[[j]]+beta_share[sim+1,5]*test_bin$effort[[j]]
    }

    #mu
    mu_func_2 <- function(y){
      lapply(combined_intensity_2,
             function(x) g_func(1,beta_share[sim+1,2],x,x_mu[sim+1,]+beta_mu[sim+1,]+
                                  beta_share[sim+1,3]*x_pp_2[sim+1,]+y)+
               0 )
    }

    test_gev <- gev_like_list(data=data_list,
                              mu_func=mu_func,
                              sigma=exp(beta_sigma[sim+1,]+x_sigma[sim+1,]),
                              xi=beta_xi[sim+1,], constant=0,
                              fe_mu[sim+1,])

    test_gev_2 <- gev_like_list(data=data_list,
                                mu_func=mu_func_2,
                                sigma=exp(beta_sigma[sim+1,]+x_sigma[sim+1,]),
                                xi=beta_xi[sim+1,], constant=0,
                                fe_mu[sim+1,])

    alpha <-
      sum(test_gev_2) -
      sum(test_gev)  +
      dnorm(beta_share[sim+1,5], 0, 100, log = T) -
      dnorm(beta_share[sim,5], 0, 100, log = T)


    # for -Inf likelihood. just stay away from that space
    if (is.nan(alpha)){alpha=-10000}

    u <- runif(1)
    alpha <- log(min(1,exp(alpha)))

    if (log( u ) >alpha){
      ind_share[sim_ind,5] <<- 0
      beta_share[sim+1,5] <- beta_share[sim,5] } else {ind_share[sim_ind,5] <<- 1}






  ##########

  function_list <- list()

  #hyperparameters
  function_list$hyper_gp_latent_pp <- hyper_gp_latent_pp
  function_list$hyper_gp_latent_pp2 <- hyper_gp_latent_pp2
  function_list$hyper_gp_latent_mu <- hyper_gp_latent_mu
  function_list$hyper_gp_latent_sigma <- hyper_gp_latent_sigma
  function_list$hyper_gp_latent_re <- hyper_gp_latent_re
  function_list$hyper_gp_latent_bin <- hyper_gp_latent_bin

  # latent fields
  function_list$x_mu <- x_mu
  function_list$x_pp <- x_pp
  function_list$x_pp_2 <- x_pp_2
  function_list$x_sigma <- x_sigma
  function_list$x_bin <- x_bin

  # intercepts
  function_list$beta_pp <- beta_pp
  function_list$beta_mu <- beta_mu
  function_list$beta_sigma <- beta_sigma
  function_list$beta_xi <- beta_xi
  function_list$beta_pp_2 <- beta_pp_2
  function_list$beta_bin <- beta_bin


  #sharing params
  function_list$beta_share <- beta_share

  #fixed effects
  function_list$fe_mu <- fe_mu
  function_list$fe_pp <- fe_pp
  function_list$fe_pp2 <- fe_pp2
  function_list$fe_bin <- fe_bin

  # list of acceptances
  function_list$ind <- ind

  # one of the log posteriors to check if it is increasing
  function_list$log_posteriors <- log_posteriors
  function_list$log_posteriors_2 <- log_posteriors_2
  function_list$log_posteriors_3 <- log_posteriors_3
  function_list$log_posteriors_4 <- log_posteriors_4

  return(function_list)

}








## test likelihood to see if initial position is viable

test_likelihood <- function(sim){
  test <- hpp_likelihood_d_list(beta_pp[sim,]+x_pp[sim,],size_A, centroids_count_list,
                                fe_pp[sim,])

  share_intensity <- test$lambda

  test_bin <- get_effort_list(beta_bin[sim,]+beta_share[sim,4]*x_pp_2[sim,], centroids_count_list_bin,
                              fe_bin=fe_bin[sim,])

  share_intensity_1 <- list()
  for (j in 1:length(share_intensity)){
    share_intensity_1[[j]] <- beta_share[sim,1]*share_intensity[[j]]+beta_share[sim,5]*test_bin$effort[[j]]
  }

  #mu
  mu_func <- function(y){
    lapply(share_intensity_1,
           function(x) g_func(1,beta_share[sim,2],x,x_mu[sim,]+beta_mu[sim,]+beta_share[sim,3]*x_pp_2[sim,]+y)+
             0 )
  }

  g_d_func_list <- lapply(share_intensity_1,
                          function(x) g_func(1,beta_share[sim,2],x,x_mu[sim,]+beta_mu[sim,]+beta_share[sim,3]*x_pp_2[sim,]) )


  test_gev <- gev_like_d_list(data_list=data_list,
                              mu_func=mu_func,
                              sigma=exp(beta_sigma[sim,]+x_sigma[sim,]), xi=beta_xi[sim,], constant=0,
                              fe_mu[sim,], g_d_func_list)
  return(test_gev)
}




