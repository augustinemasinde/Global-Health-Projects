library(dplyr)
library(binom)
library(ggplot2)
library(MASS)
library(data.table)

#=============================================================================================================================#
#                           MCMC function (single dataset)                                                                    #
#=============================================================================================================================#

#============================================================#
#       Simple FoI model (single dataset)                    #
#============================================================#

# # prior distributions for each parameter #
prior <- function(par) {
  
  lambda <- par[1]; se<- par[2]; sp <- par[3]
  
  lambda = lambda
  
  lambda_prior = dunif(lambda, min = 0.000001, max = 12, log = T)
  
  se = se
  se_prior = dbeta(se, alpha_se, beta_se, log = T)
  
  sp = sp
  sp_prior = dbeta(sp, alpha_sp, beta_sp, log = T)
  
  return(sum(c(lambda_prior, se_prior, sp_prior)))
  
}


# likelihood calculation given each parameter set & data
loglike_simple <- function(data, par){
  predicted_seroprevalence = predicted_prev_func(age = data$age, par)
  sum(dbinom(data$pos, data$n, predicted_seroprevalence, log=T))
}


# Posterior Function
posterior_function_simple <- function(data, par){
  loglike_simple(data, par) + prior(par)
}

# proposal function for MCMC
proposal_function_simple <- function(par, cov) {
  #lambda <- par[1]; se<- par[2]; sp <- par[3]
  
  ## draw propopsals all from a multivariate normal
  repeat {
    proposed <- mvrnorm(1, par, cov)
    if(all(proposed[2]>0 & all(proposed[2]<1))& all(proposed[3]>0 & all(proposed[3]<1)) & all(proposed[1]>0 & all(proposed[1]<12))){break}
  }
  
  
  return(proposed)
  
}  


# Run MCMC Function 
MCMC_simple_model <- function(inits,  number_of_iterations, cov, fitting) {
  
  # Storage for Output
  MCMC_output <- matrix(nrow = number_of_iterations + 1, ncol=length(inits)) # ncol is number of parameters
  MCMC_output[1,] <- inits
  Acceptances <- vector(length = number_of_iterations + 1)
  LogLikelihood_storage <- vector(length = number_of_iterations + 1)
  Logposterior_storage <- vector(length = number_of_iterations + 1)
  
  
  # Running the Actual MCMC
  for (i in 1:number_of_iterations){
    
    if(fitting == "single dataset"){
      proposed_parameter_value <- proposal_function_simple(MCMC_output[i,], cov)  #new proposed paramater value(s) with a given s.d. (step size)
    }
    
    if(fitting == "multiple datasets"){
      proposed_parameter_value <- proposal_function_simple_multidata(MCMC_output[i,], cov)  #new proposed paramater value(s) with a given s.d. (step size)
    }
    
    if(fitting == "single dataset"){
      current_likelihood <- loglike_simple(data, MCMC_output[i,]) # likelihood 
    }
    
    if(fitting == "multiple datasets"){
      current_likelihood <- loglike_simple_multidata(data, MCMC_output[i,]) # likelihood 
    }
    
    if(fitting == "single dataset"){
      current_posterior <- posterior_function_simple(data, MCMC_output[i,]) # current posterior likelihood from MCMC
    }
    
    if(fitting == "multiple datasets"){
      current_posterior <- posterior_function_simple_multidata(data, MCMC_output[i,]) # current posterior likelihood from MCMC
    }
    
    if(fitting == "single dataset"){
      proposed_posterior <- posterior_function_simple(data, proposed_parameter_value) # proposed posterior likelihood with new proposed par value
    }
    
    if(fitting == "multiple datasets"){
      proposed_posterior <- posterior_function_simple_multidata(data, proposed_parameter_value) # proposed posterior likelihood with new proposed par value
    }
    
    likelihood_ratio = exp(proposed_posterior - current_posterior);
    
    if(i %% (number_of_iterations/20) == 0){
      message(round(i/number_of_iterations*100), ' % completed')
    }
    
    if(runif(1) < likelihood_ratio) {
      
      MCMC_output[i + 1,] <- proposed_parameter_value
      Acceptances[i] <- 1
      LogLikelihood_storage[i + 1] <- current_likelihood
      Logposterior_storage[i + 1] <- proposed_posterior
      
    } else{
      
      MCMC_output[i + 1,] <- MCMC_output[i,]
      Acceptances[i] <- 0
      LogLikelihood_storage[i + 1] <- current_likelihood
      Logposterior_storage[i + 1] <- current_posterior
      
    }
    
  } # likelihood ratio comparison step (exponentiated because on log scale) - accept (1) if posterior proposed improved LL on posterior current (runif- generates random variable between 0 -1)
  
  list <- list()
  list[["MCMC_Output"]] <- MCMC_output
  list[["Acceptances"]] <- Acceptances
  list[["Likelihood_Output"]] <- LogLikelihood_storage
  list[["Posterior_Output"]] <- Logposterior_storage
  return(list)
  
}

#===========================================================#
#       Reversible FoI model (single dataset)               #
#===========================================================#

# prior distributions for each parameter #
# specify lambda median from simple model to inform the lognromal prior
prior_function_reversible <- function(par, simple_lambda_median) {
  
  lambda <- par[1]; se<- par[2]; sp <- par[3]; rho <- par[4]
  
  lambda = lambda
  lambda_prior = dlnorm(lambda, log(simple_lambda_median), 1, log=TRUE)
  
  se = se
  se_prior = dbeta(se, alpha_se, beta_se, log = T)
  
  sp = sp
  sp_prior = dbeta(sp, alpha_sp, beta_sp, log = T)
  
  rho = rho
  rho_prior = dunif(rho, min = 0.000001, max = 12, log= T)
  
  return(sum(c(lambda_prior, se_prior, sp_prior, rho_prior)))
  
}

# when calculating postrior after MCMC fitting (i.e. replace simple median lambda with new fitted lambda)
# prior_function_reversible2 <- function(par, simple_lambda_median) {
#   
#   lambda <- par[1]; se<- par[2]; sp <- par[3]; rho <- par[4]
#   
#   lambda = lambda
#   lambda_prior = dlnorm(lambda, log(simple_lambda_median), 1, log=TRUE)
#   
#   se = se
#   se_prior = dbeta(se, alpha_se, beta_se, log = T)
#   
#   sp = sp
#   sp_prior = dbeta(sp, alpha_sp, beta_sp, log = T)
#   
#   rho = rho
#   rho_prior = dunif(rho, min = 0.000001, max = 12, log= T)
#   
#   return(sum(c(lambda_prior, se_prior, sp_prior, rho_prior)))
#   
# }

# likelihood calculation given each parameter set & data
log_lik_func_reversible <- function(data, par){
  
  predicted_seroprevalence = predicted_prev_reversible_func(data$age, par)
  
  sum(dbinom(data$pos, data$n, predicted_seroprevalence, log=T))
}

# posterior calculation (for MCMC fitting)
posterior_function_reversible<- function(data, par, simple_lambda_median){
  
  lambda <- par[1]; se<- par[2]; sp <- par[3]; rho <- par[4]
  
  log_lik_func_reversible(data, par) + prior_function_reversible(par, simple_lambda_median)
}

# posterior calculation (post mCMC fitting when have fitted lambda)
# posterior_function_reversible2 <- function(data, par){
#   
#   lambda <- par[1]; se<- par[2]; sp <- par[3]; rho <- par[4]
#   
#   log_lik_func_reversible(data, par) + prior_function_reversible2(par)
# }

# proposal function for MCMC
proposal_function_reversible <- function(par, cov) {
  
  ## draw propopsals all from a multivariate normal
  repeat {
    proposed <- mvrnorm(1, par, cov)
    if(all(proposed[2]>0 & all(proposed[2]<1)) & all(proposed[3]>0 & all(proposed[3]<1)) & all(proposed[1]>0 & all(proposed[1]<12)) & all(proposed[4]>0 & all(proposed[4]<12))){break}
  }
  
  return(proposed)
  
}  

MCMC_reversible_model <- function(inits,  number_of_iterations, cov, simple_lambda_median, fitting) {
  
  # Storage for Output
  MCMC_output <- matrix(nrow = number_of_iterations + 1, ncol=length(inits)) # ncol is number of parameters
  MCMC_output[1,] <- inits
  Acceptances <- vector(length = number_of_iterations + 1)
  LogLikelihood_storage <- vector(length = number_of_iterations + 1)
  Logposterior_storage <- vector(length = number_of_iterations + 1)
  
  
  # Running the Actual MCMC
  for (i in 1:number_of_iterations){
    
    
    if(fitting == "single dataset"){
      proposed_parameter_value <- proposal_function_reversible(MCMC_output[i,], cov)  #new proposed paramater value(s) with a given s.d. (step size)
    }
    
    if(fitting == "multiple datasets"){
      proposed_parameter_value <- proposal_function_reversible_multidata(MCMC_output[i,], cov)  #new proposed paramater value(s) with a given s.d. (step size)
    }
    
    if(fitting == "single dataset"){
      current_likelihood <- log_lik_func_reversible(data, MCMC_output[i,]) # likelihood 
    }
    
    if(fitting == "multiple datasets"){
      current_likelihood <- loglike_reversible_multidata(data, MCMC_output[i,]) # likelihood 
    }
    
    if(fitting == "single dataset"){
      current_posterior <- posterior_function_reversible(data, MCMC_output[i,], simple_lambda_median) # current posterior likelihood from MCMC
    }
    
    if(fitting == "multiple datasets"){
      current_posterior <- posterior_function_reversible_multidata(data, MCMC_output[i,], simple_lambda_median) # current posterior likelihood from MCMC
    }
    
    if(fitting == "single dataset"){
      proposed_posterior <- posterior_function_reversible(data, proposed_parameter_value, simple_lambda_median) # proposed posterior likelihood with new proposed par value
    }
    
    if(fitting == "multiple datasets"){
      proposed_posterior <- posterior_function_reversible_multidata(data, proposed_parameter_value, simple_lambda_median) # proposed posterior likelihood with new proposed par value
    }
    
    likelihood_ratio = exp(proposed_posterior - current_posterior);
    
    if(i %% (number_of_iterations/20) == 0){
      message(round(i/number_of_iterations*100), ' % completed')
    }
    
    if(runif(1) < likelihood_ratio) {
      
      MCMC_output[i + 1,] <- proposed_parameter_value
      Acceptances[i] <- 1
      LogLikelihood_storage[i + 1] <- current_likelihood
      Logposterior_storage[i + 1] <- proposed_posterior
      
      
    } else{
      
      MCMC_output[i + 1,] <- MCMC_output[i,]
      Acceptances[i] <- 0
      LogLikelihood_storage[i + 1] <- current_likelihood
      Logposterior_storage[i + 1] <- current_posterior
      
      
    }
    
  } # likelihood ratio comparison step (exponentiated because on log scale) - accept (1) if posterior proposed improved LL on posterior current (runif- generates random variable between 0 -1)
  
  list <- list()
  list[["MCMC_Output"]] <- MCMC_output
  list[["Acceptances"]] <- Acceptances
  list[["Likelihood_Output"]] <- LogLikelihood_storage
  list[["Posterior_Output"]] <- Logposterior_storage
  return(list)
}


#=============================================================================================================================#
#                                Estimating diagnostic (prior) parameters                                                     #
#=============================================================================================================================#

estimate_alpha_beta_par_diagnostic <-
  function(input_alpha, input_beta, target_l, target_u) {
    
    # Define targets we want to fit to
    #test_alpha <- 84.5 # original
    #test_beta <- 15.5 # original
    
    # mean (of beta prior distribution given by alpha and beta shape parameter as:)
    target_mean <- input_alpha / (input_alpha + input_beta)
    
    # Lower CI
    #target_l <- qbeta(0.025, input_alpha, input_beta)
    # Upper CI
    #target_u <- qbeta(0.975, input_alpha, input_beta)
    
    # Actual upper and low targets to calibrate on for below #
    target_l <- target_l
    target_u <- target_u
    
    # Find parameters
    # Vector of possible alphas
    alpha <- seq(0, 100, 0.1)
    # Vector of corresponding betas, such that alpha / (alpha + beta) = m
    beta <- (alpha - (alpha * target_mean)) / target_mean
    # Lower CI | alpha, beta
    l <- qbeta(0.025, alpha, beta)
    # Upper CI | alpha, beta
    u <- qbeta(0.975, alpha, beta)
    # Best fit
    best <- which.min(abs(l - target_l) + abs(u - target_u))
    
    # Recover input parameters
    alpha[best]
    beta[best]
    
    actual_alpha <- alpha[best] # updated, best fit
    actual_beta <- beta[best] # updated, best fit
    
    actual_mean <- actual_alpha / (actual_alpha + actual_beta) # mean
    actual_mean
    actual_lower <- qbeta(0.025, actual_alpha, actual_beta)
    actual_lower
    actual_upper <- qbeta(0.975, actual_alpha, actual_beta)
    actual_upper
    
    # 0.6592385 (l) and  0.9939759 (u) # given by 11 for alpha and 1.2612613 for beta
    # p = seq(0, 1, length = 100)
    # #par(mfrow = c(1, 1)) #
    # 
    # a <- plot(
    #   p,
    #   dbeta(p, actual_alpha, actual_beta),
    #   ylab = "density",
    #   type = "l",
    #   col = 4
    # ) #
    
    
    
    return(list(actual_alpha, actual_beta))
  }



#==============================================================================================================================#
#                                       MASTER SCRIPT  - single dataset fitting                                                #
#==============================================================================================================================#

rm(list = ls())

#=============================================#
#       Load data files                       #
#=============================================#

#=======================#                                                                                             
# Initiatie sub-scripts #                                                                                             
source('libraries.R')
source('plot_data_modelfit_functions.R')
source('diagnostic_parameter_functions.R')
source('predicted_prevalence_functions.R')
source('MCMC_functions_singledataset.R')
source('plot_MCMC_output_functions.R')
source('process_MCMC_output_functions.R')

#======================#
#    load data         #
#======================#

test_data_singledataset <- read.csv("~/human_tsol_FoI_modelling/data/test_data_singledataset.csv")

data <- test_data_singledataset

#data <- data_frame_Theis

names(data) <- c("age", "pos", "n","prev","lower","upper")

plot_ageprev_func(data) # plot age-prevalence data

#===============================================#
#  optimise parameters for diagnostic priors    #

# sensitivity #
sensitivity_parameters <- estimate_alpha_beta_par_diagnostic(input_alpha = 96, input_beta = 4, 
                                                             target_l = 0.93, target_u = 0.99)

alpha_se <- sensitivity_parameters[[1]]
alpha_se 
beta_se <- sensitivity_parameters[[2]]
beta_se

# specifcity #
specificity_parameters <- estimate_alpha_beta_par_diagnostic(input_alpha = 98, input_beta = 2, 
                                                             target_l = 0.96, target_u = 1)

alpha_sp <- specificity_parameters[[1]]
alpha_sp
beta_sp <- specificity_parameters[[2]]
beta_sp

#===============================================#
#  run MCMC (single dataset; simple FoI model)  #

inits1 <- c(0.1, 0.93, 0.96)   # Initial parameter values to initiate chains
inits2 <- c(0.0001, 0.99, 0.999) # e.g. Ab−EITB, rT24H (se: 0.96 (0.93-0.99), sp: 0.98 (0.96-1), Noh et al. 2014) 
sd <- 0.001 # set standard deviation of proposal distribution; aim for 0.25 acceptance
cov <- diag(sd^2, 3) # covariance
niter <- 100000 # number of iterations
burnin <- 50000 # burnin (chains to discard before convergence)

# run MCMC (chain 1)
set.seed(123) # for reproducibility
simple_out_chain1 <- MCMC_simple_model(inits1, niter, cov, fitting = "single dataset")  # initiate the MCMC

# run MCMC (chain 1)
set.seed(123)
simple_out_chain2 <- MCMC_simple_model(inits2, niter, cov, fitting = "single dataset")  # initiate the MCMC

# whats the acceptance ratio (aiming for 0.25)
sum(simple_out_chain1$Acceptances)/niter
sum(simple_out_chain2$Acceptances)/niter

# plot MCMC outputs #
chains_plot <- chains_plot_func(inits1 = inits1, chain1 = simple_out_chain1, chain2 = simple_out_chain2,
                                model = "simple")

histograms_plot <- histogram_plot_func(inits1 = inits1, burnin = burnin, niter = niter, 
                                       chain1 = simple_out_chain1, chain2 = simple_out_chain2, model = "simple")

#=============================================# 
#         process MCMC outputs                # 

chains1_output <- simple_out_chain1$MCMC_Output
chains2_output <- simple_out_chain2$MCMC_Output

# remove burnin and proceed with reducing autocorrelation (thinning by sub-sampling)
PC_simple <- Process_chains(chains1_output, chains2_output, burnin = burnin, sample = 50) # set burnin to 0 if already

# View the process chains (from the autocorrelation plots, sampling every 20th value seems appropriate)
plot_chains(PC_simple[[1]], PC_simple[[2]])

# check autocorrelation of chains for each parameter (to inform sub-sampling)
check_autocorr <- determine_autocorrelation_func1(processed_chain = PC_simple, number_datasets = 1)

check_autocorr[1] # autocorrelation significance parameter 1 (e.g. lambda)
check_autocorr[2] # autocorrelation significance parameter 2 (e.g. se)
check_autocorr[3] # autocorrelation significance parameter 1 (e.g. sp)

# plot loglikelihood
loglikchains_plot_func(chain1 = simple_out_chain1, chain2 = simple_out_chain2)

#==============================================================================#
# Obtain parameter values (median & credible) & plot posterior distributions   #

simple_model_parameters <- obtain_parameter_values_func(processed_chains = PC_simple, model = "simple", number_datasets = 1)
simple_model_parameters

plot_posterior_distrib_func(processed_chains = PC_simple, model = "simple", number_datasets = 1)

#============================================================#
# calculate Deviance Information Criterion (DIC) - model fit #

DIC_result <- calculate_DIC_func1(chain = simple_out_chain1, burnin = burnin, subsample = 50, 
                                  parameters = simple_model_parameters, number_datasets = 1)
DIC_result  # 1) D bar model1, 2) modal posterior likelihood, 3) modal posterior deviance, 4) DIC

#====================================================================================================================#
# calculate (with posterior parameter estimates) predicted (sero)prevalence and unceetainty intervals (for plotting) #

predicted_prev_output <- calculate_predicted_prevalence_function(max_age_toplot = 90, data = data, 
                                                                 pars = simple_model_parameters,
                                                                 processed_chains = PC_simple, model = "simple")

predicted_prev_output[[1]] # predicted prevalence plot (ylim 0-100% prev)
predicted_prev_output[[2]] # predicted prevalence plot (ylim 0-50% prev)
predicted_prev_output[[3]] # predicted prevalence plot (ylim 0-25% prev)

#=========================================================================================================================#
#===================================================#
#  run MCMC (single dataset; reversible FoI model)  #

inits1 <- c(0.004, 0.93, 0.96, 0.001)   # Initial parameter values to initiate chains
inits2 <- c(0.0001, 0.99, 0.999, 0.01) # e.g. Ab−EITB, rT24H (se: 0.96 (0.93-0.99), sp: 0.98 (0.96-1), Noh et al. 2014) 

sd <- 0.004 # set standard deviation of proposal distribution; aim for 0.25 acceptance
cov <- diag(sd^2, 4)# covariance
niter <- 100000 # number of iterations
burnin <- 50000 # burnin (chains to discard before convergence)

# run MCMC (chain 1)
set.seed(123) # for reproducibility
reversible_out_chain1 <- MCMC_reversible_model(inits1, niter, cov, simple_lambda_median = 0.00042,
                                               fitting = "single dataset")  # initiate the MCMC (& specify lambda median from simple model to inform the lognromal prior)

# run MCMC (chain 1)
set.seed(123)
reversible_out_chain2 <- MCMC_reversible_model(inits2, niter, cov, simple_lambda_median = 0.00042,
                                               fitting = "single dataset")   # initiate the MCMC (& specify lambda median from simple model to inform the lognromal prior)

# whats the acceptance ratio (aiming for 0.25)
sum(reversible_out_chain1$Acceptances)/niter
sum(reversible_out_chain2$Acceptances)/niter

# plot MCMC outputs #
chains_plot <- chains_plot_func(inits1 = inits1, chain1 = reversible_out_chain1, chain2 = reversible_out_chain2, 
                                model = "reversible")

histograms_plot <- histogram_plot_func(inits1 = inits1, burnin = burnin, niter = niter, 
                                       chain1 = reversible_out_chain1, chain2 = reversible_out_chain2,
                                       model = "reversible")

#=============================================# 
#         process MCMC outputs                # 

chains1_output <- reversible_out_chain1$MCMC_Output
chains2_output <- reversible_out_chain2$MCMC_Output

# remove burnin and proceed with reducing autocorrelation (thinning by sub-sampling)
PC_reversible <-Process_chains(chains1_output, chains2_output, burnin = burnin, sample = 50) # set burnin to 0 if already

# View the process chains (from the autocorrelation plots, sampling every 20th value seems appropriate)
plot_chains(PC_reversible[[1]], PC_reversible[[2]])

# check autocorrelation of chains for each parameter (to inform sub-sampling)
check_autocorr <- determine_autocorrelation_func2(processed_chain = PC_reversible, number_datasets = 1)

check_autocorr[1] # autocorrelation significance parameter 1 (e.g. lambda)
check_autocorr[2] # autocorrelation significance parameter 2 (e.g. se)
check_autocorr[3] # autocorrelation significance parameter 1 (e.g. sp)
check_autocorr[4] # autocorrelation significance parameter 1 (e.g. sp)

# plot loglikelihood
loglikchains_plot_func(chain1 = reversible_out_chain1, chain2 = reversible_out_chain2)

#==============================================================================#
# Obtain parameter values (median & credible) & plot posterior distributions   #

reversible_model_parameters <- obtain_parameter_values_func(processed_chains = PC_reversible, model = "reversible", number_datasets = 1)
reversible_model_parameters

plot_posterior_distrib_func(processed_chains = PC_reversible, model = "reversible", number_datasets = 1)

#============================================================#
# calculate Deviance Information Criterion (DIC) - model fit #

DIC_result <- calculate_DIC_func2(chain = reversible_out_chain1, burnin = burnin, subsample = 50, 
                                  parameters = reversible_model_parameters, number_datasets = 1,
                                  simple_lambda_median = 0.00042)
DIC_result  # 1) D bar model1, 2) modal posterior likelihood, 3) modal posterior deviance, 4) DIC

#====================================================================================================================#
# calculate (with posterior parameter estimates) predicted (sero)prevalence and unceetainty intervals (for plotting) #

predicted_prev_output <- calculate_predicted_prevalence_function(max_age_toplot = 90, data = data, 
                                                                 pars = reversible_model_parameters,
                                                                 processed_chains = PC_reversible, model = "reversible")

predicted_prev_output[[1]] # predicted prevalence plot (ylim 0-100% prev)
predicted_prev_output[[2]] # predicted prevalence plot (ylim 0-50% prev)
predicted_prev_output[[3]] # predicted prevalence plot (ylim 0-25% prev)


#=============================================================================================================================#
#                                       Plot MCMC outputs                                                                     #
#=============================================================================================================================#

#======================================#
#       Single datasets                #

#### chains plot (single dataset) - not processed #####

chains_plot_func <- function(inits1, chain1, chain2, model) {
  
  if(model == "simple"){
  par(mfrow = (c(1, length(inits1))))
  for (i in 1:length(inits1)) {
    if (i == 1) {
      ylab = "lambda"
    } else if (i == 2) {
      ylab = "se"
    } else {
      ylab = "sp"
    }
    plot(
      simple_out_chain1$MCMC_Output[, i],
      type = "l",
      ylab = ylab,
      xlab = "iter"
    )
    lines(simple_out_chain2$MCMC_Output[, i], col = "red")
    
  }}
  
  if(model == "reversible"){
    par(mfrow=(c(1,length(inits1))))
    for (i in 1:length(inits1)) {
      if (i==1) {
        ylab="lambda"
      } else if (i==2) {
        ylab="se"
      } else if (i==3) {
        ylab="sp"
      } else {
        ylab="rho"
      }
      plot(reversible_out_chain1$MCMC_Output[,i], type = "l", ylab=ylab, xlab="iter", col = "black")
      lines(reversible_out_chain2$MCMC_Output[,i], col="red")
    }}
  
}


#### posterior histogram plots - not processed #####

histogram_plot_func <- function(inits1, burnin, niter, chain1, chain2, model) {
    
  if(model == "simple"){
  par(mfrow = (c(1, length(inits1))))
    for (i in 1:length(inits1)) {
      if (i == 1) {
        ylab = "lambda"
      } else if (i == 2) {
        ylab = "se"
      } else {
        ylab = "sp"
      }
      
      hist(
        c(simple_out_chain1$MCMC_Output[burnin:niter, i], simple_out_chain2$MCMC_Output[burnin:niter, i]),
        xlab = ylab,
        main = ""
      )
      
    }}
  
  if(model == "reversible"){
    par(mfrow=(c(1,length(inits1))))
    for (i in 1:length(inits1)) {
      if (i==1) {
        ylab="lambda"
      } else if (i==2) {
        ylab="se"
      } else if (i==3) {
        ylab="sp"
      } else {
        ylab="rho"
      }
      
      hist(c(reversible_out_chain1$MCMC_Output[burnin:niter,i],reversible_out_chain2$MCMC_Output[burnin:niter,i]), 
           xlab = ylab, main="")
      
    }}
  
  }

##### plot log likelihood #####

loglikchains_plot_func <- function(chain1, chain2){
  par(mfrow=c(1,1))
  plot(chain1$Likelihood_Output, t='l', ylab='Loglikelihood', xlab='iteration')
  lines(chain2$Likelihood_Output, col='red')
}

#### function to plot chains - post-processing ####
plot_chains <- function(run1, run2){
  par(mfrow=c(ncol(run1),1))
  
  for(i in 1:ncol(run1)){
    plot(run1[,i], t='l', col='deeppink',
         ylim=c(min(c(run1[,i], run2[,i])),max(c(run1[,i], run2[,i]))),
         xlab='', ylab=paste('Parameter', i, sep=' '))
    lines(run2[,i], col='dodgerblue')
  }
  
}

plot_chains_multidatasets <- function(run1, run2, inits, number_datasets, model){
  par(mfrow=c(2,length(inits)/2))
  
  if(number_datasets == 5 && model == "simple"){
    for(i in 1:ncol(run1)){
      if (i == 1) {
        ylab = "sp"
      } else if (i == 2) {
        ylab = "se"
      }  else if (i ==3) {
        ylab="lambda (site 1)"
      } else if (i==4) {
        ylab="lambda (site 2)"
      } else if (i==5) {
        ylab="lambda (site 3)"
      } else if (i==6) {
        ylab="lambda (site 4)"
      } else {
        ylab="lambda (site 5)"
      }
      plot(run1[,i], t='l', col='deeppink',
           ylim=c(min(c(run1[,i], run2[,i])),max(c(run1[,i], run2[,i]))),
           xlab='', ylab = ylab)
      lines(run2[,i], col='dodgerblue')
    }}
  
  if(number_datasets == 5 && model == "reversible"){
    for(i in 1:ncol(run1)){
      if (i == 1) {
        ylab = "sp"
      } else if (i == 2) {
        ylab = "se"
      }  else if (i ==3) {
        ylab="lambda (site 1)"
      } else if (i==4) {
        ylab="lambda (site 2)"
      } else if (i==5) {
        ylab="lambda (site 3)"
      } else if (i==6) {
        ylab="lambda (site 4)"
      } else if (i==7) {
        ylab="lambda (site 5)"
      }  else if (i==8) {
        ylab="rho (site 1)"
      } else if (i==9) {
        ylab="rho (site 2)"
      } else if (i==10) {
        ylab="rho (site 3)"
      } else if (i==11) {
        ylab="rho (site 4)"
      } else {
        ylab="rho (site 5)"
      }
      
      plot(run1[,i], t='l', col='deeppink',
           ylim=c(min(c(run1[,i], run2[,i])),max(c(run1[,i], run2[,i]))),
           xlab='', ylab = ylab)
      lines(run2[,i], col='dodgerblue')
    }}
}


#### plot posterior distributions (after processing) ###

plot_posterior_distrib_func <- function(processed_chains, model, number_datasets){
  
  if(model == "simple" && number_datasets == 1){
  par(mfrow=c(1,3))
  
  hist(c(processed_chains[[1]][,1], processed_chains[[2]][,1]), breaks=30, xlab='Lambda', main="")  # Parameter 1 - lambda 
  hist(c(processed_chains[[1]][,2], processed_chains[[2]][,2]), breaks=30, xlab='sensitivity', main="")      # Parameter 2 - se
  hist(c(processed_chains[[1]][,3], processed_chains[[2]][,3]), breaks=30, xlab='specificity', main="")      # Parameter 3 - sp
  }
  
  if(model == "simple" && number_datasets == 5){
    par(mfrow=c(2,ncol(PC_simple[[1]])/2))
    
    hist(c(processed_chains[[1]][,1], processed_chains[[2]][,1]), breaks=30, xlab='specifcity', main="")  # Parameter 1  
    hist(c(processed_chains[[1]][,2], processed_chains[[2]][,2]), breaks=30, xlab='sensitivity', main="")      # Parameter 2 
    hist(c(processed_chains[[1]][,3], processed_chains[[2]][,3]), breaks=30, xlab='lambda (site 1)', main="")      # Parameter 3 
    hist(c(processed_chains[[1]][,4], processed_chains[[2]][,4]), breaks=30, xlab='lambda (site 2)', main="")  
    hist(c(processed_chains[[1]][,5], processed_chains[[2]][,5]), breaks=30, xlab='lambda (site 3)', main="")     
    hist(c(processed_chains[[1]][,6], processed_chains[[2]][,6]), breaks=30, xlab='lambda (site 4)', main="")      
    hist(c(processed_chains[[1]][,7], processed_chains[[2]][,7]), breaks=30, xlab='lambda (site 5)', main="")      
  }
  
  if(model == "reversible" && number_datasets == 1){
    par(mfrow=c(1,4))
    hist(c(processed_chains[[1]][,1], processed_chains[[2]][,1]), breaks=30, xlab='Lambda', main="")  # Parameter 1 - lambda 
    hist(c(processed_chains[[1]][,2], processed_chains[[2]][,2]), breaks=30, xlab='sensitivity', main="")      # Parameter 2 - se
    hist(c(processed_chains[[1]][,3], processed_chains[[2]][,3]), breaks=30, xlab='specificity', main="")      # Parameter 3 - sp
    hist(c(processed_chains[[1]][,4], processed_chains[[2]][,4]), breaks=30, xlab='rho', main="")      # Parameter 4 - rho
  }
  
  if(model == "reversible" && number_datasets == 5){
    par(mfrow=c(2,ncol(PC_reversible[[1]])/2))
    
    hist(c(processed_chains[[1]][,1], processed_chains[[2]][,1]), breaks=30, xlab='specifcity', main="")   
    hist(c(processed_chains[[1]][,2], processed_chains[[2]][,2]), breaks=30, xlab='sensitivity', main="")      
    hist(c(processed_chains[[1]][,3], processed_chains[[2]][,3]), breaks=30, xlab='lambda (site 1)', main="")      
    hist(c(processed_chains[[1]][,4], processed_chains[[2]][,4]), breaks=30, xlab='lambda (site 2)', main="")  
    hist(c(processed_chains[[1]][,5], processed_chains[[2]][,5]), breaks=30, xlab='lambda (site 3)', main="")     
    hist(c(processed_chains[[1]][,6], processed_chains[[2]][,6]), breaks=30, xlab='lambda (site 4)', main="")
    hist(c(processed_chains[[1]][,7], processed_chains[[2]][,7]), breaks=30, xlab='lambda (site 5)', main="")
    hist(c(processed_chains[[1]][,8], processed_chains[[2]][,8]), breaks=30, xlab='rho (site 1)', main="")
    hist(c(processed_chains[[1]][,9], processed_chains[[2]][,9]), breaks=30, xlab='rho (site 2)', main="")      
    hist(c(processed_chains[[1]][,10], processed_chains[[2]][,10]), breaks=30, xlab='rho (site 3)', main="")  
    hist(c(processed_chains[[1]][,11], processed_chains[[2]][,11]), breaks=30, xlab='rho (site 4)', main="")     
    hist(c(processed_chains[[1]][,12], processed_chains[[2]][,12]), breaks=30, xlab='rho (site 5)', main="")      

  }
  
  
}




#=============================================================================================================================#
#                           plot data & model fit functions                                                                   #
#=============================================================================================================================#


plot_ageprev_func <- function(data) {
  
  p <- ggplot() +    
    geom_point(data=data, aes(x=age, y=prev))+
    #geom_errorbar(data=predicted_simple,aes(x=age, y=prev, ymin=lower, ymax=upper, width=wd))+
    geom_errorbar(data=data,aes(x=age, y=prev, ymin=lower, ymax=upper))+
    #facet_wrap(~ref, scales = "free")+
    #facet_wrap(~ref)+
    ylim(0,1)+
    labs(x="Age (months) of human host", y="(Sero)prevalence (0 -1)")+
    theme_bw() +
    theme(strip.background=element_rect(fill=NA, color=NA),
          strip.text=element_text(size=11, face= "bold"),
          axis.text.x = element_text(size = 18),
          axis.text.y = element_text(size = 18),
          plot.title = element_text(size = 18, hjust=0.5),
          axis.title.x = element_text(size = 18, face= "bold"),
          axis.title.y = element_text(size = 16, angle = 90, face= "bold"),
          legend.position = c(0.83, 0.15),
          legend.title=element_text(size=20), 
          legend.text=element_text(size=16))
  
  return(p)
}




#=============================================================================================================================#
#                                   Predicted prevalence functions                                                            #
#=============================================================================================================================#

#===============================#
# Simple Model (single dataset) #

# for MCMC 
predicted_prev_func <- function(age, par){
  
  lambda <- par[1]; se<- par[2]; sp <- par[3]
  
  tp <-  1 - exp(-lambda * (data$age))   # true prevalence
  op <- (1-sp) + (se+sp-1)*tp      # observed prevalence Diggle et al 2011 Epi Res Int
  op
}

# for single & multiple datasets (once fitted w/ parameter estimates to obtain prevalence curves for each datasets)

predicted_prev_func2 <- function(age, par){
  
  lambda <- par[1]; se<- par[2]; sp <- par[3]
  
  tp <-  1 - exp(-lambda * (age))   # true prevalence
  op <- (1-sp) + (se+sp-1)*tp      # observed prevalence Diggle et al 2011 Epi Res Int
  op
}


#===================================#
# reversible model (single dataset) #

predicted_prev_reversible_func <- function(age, par){
  
  lambda <- par[1]; se<- par[2]; sp <- par[3]; rho <- par[4]
  
  tp <-  (lambda/(lambda + rho)) * (1 - exp(-(lambda + rho) *(data$age)))   # true prevalence
  op <- (1-sp) + (se+sp-1) * tp      # observed prevalence Diggle et al 2011 Epi Res Int
  op
}

predicted_prev_reversible_func2 <- function(age, par){
  
  lambda <- par[1]; se<- par[2]; sp <- par[3]; rho <- par[4]
  
  tp <-  (lambda/(lambda + rho)) * (1 - exp(-(lambda + rho) *(age)))   # true prevalence
  op <- (1-sp) + (se+sp-1) * tp      # observed prevalence Diggle et al 2011 Epi Res Int
  op
}

#================================================================#
# produce predicted prevalence curves (fitted to single dataset) #

calculate_predicted_prevalence_function <- function (max_age_toplot, data, pars, processed_chains, model) {
  
  if(model == "simple"){
    # specify posterior median parameters to enable predicted prevalence calculation
    lambda.median <- pars[[1]]
    se.median <- pars[[2]]
    sp.median <- pars[[3]]
    
    # set up dataframe and calculate (median) predicted prevalence 
    age_dum <- seq(from=0, to = max_age_toplot, by=0.005)  ## If not already performed this step for simple catalytic model
    fitted_curve_df <- as.data.frame(age_dum) ## make sequence of numbers (mean ages) for predicted variable
    names(fitted_curve_df)[names(fitted_curve_df)=="age_dum"] <- "age"
    
    predicted_median_curve <- full_join(fitted_curve_df, data) 
    predicted_median_curve$predicted <- sapply(1:nrow(predicted_median_curve), 
                                               function(i) predicted_prev_func2(age = predicted_median_curve$age[i], 
                                                                                c(lambda.median, se.median, 
                                                                                  sp.median)))
    
    # create uncertainty (credible interval) of model run resulting from posterior #
    subsampled_model_outputs <- matrix(NA, nrow = length(processed_chains[[1]][,1]), ncol = length(seq(0, 90, 0.005)))
    
    for (i in 1:length(processed_chains[[1]][,1])){
      
      single_model_output <- predicted_prev_func2(seq(0, max_age_toplot, 0.005),c(processed_chains[[1]][i,1],
                                                                                  processed_chains[[1]][i,2],
                                                                                  processed_chains[[1]][i,3]))
      subsampled_model_outputs[i, ] <- single_model_output
      
    }
    
    lower_credible_interval_processed <- apply(subsampled_model_outputs, MARGIN = 2, quantile, prob = 0.025)
    upper_credible_interval_processed <- apply(subsampled_model_outputs, MARGIN = 2, quantile, prob = 0.975)
    
    Lower_obs <- as.data.frame(lower_credible_interval_processed)
    Upper_obs <- as.data.frame(upper_credible_interval_processed)
    
    predicted_CrI <- cbind(Lower_obs, Upper_obs)
    predicted_CrI <- as.data.frame(predicted_CrI)
    predicted_CrI$age <- fitted_curve_df$age
    
  }
  
  if(model == "reversible"){
    # specify posterior median parameters to enable predicted prevalence calculation
    lambda.median <- pars[[1]]
    se.median <- pars[[2]]
    sp.median <- pars[[3]]
    rho.median <- pars[[4]]
    
    # set up dataframe and calculate (median) predicted prevalence 
    age_dum <- seq(from=0, to = max_age_toplot, by=0.005)  ## If not already performed this step for simple catalytic model
    fitted_curve_df <- as.data.frame(age_dum) ## make sequence of numbers (mean ages) for predicted variable
    names(fitted_curve_df)[names(fitted_curve_df)=="age_dum"] <- "age"
    
    predicted_median_curve <- full_join(fitted_curve_df, data) 
    predicted_median_curve$predicted <- sapply(1:nrow(predicted_median_curve), 
                                               function(i) predicted_prev_reversible_func2(
                                                 age = predicted_median_curve$age[i], 
                                                 c(lambda.median, se.median, sp.median, rho.median)))
    
    # create uncertainty (credible interval) of model run resulting from posterior #
    subsampled_model_outputs <- matrix(NA, nrow = length(processed_chains[[1]][,1]), ncol = length(seq(0, 90, 0.005)))
    
    for (i in 1:length(processed_chains[[1]][,1])){
      
      single_model_output <- predicted_prev_reversible_func2(seq(0, max_age_toplot, 0.005),c(processed_chains[[1]][i,1],
                                                                                             processed_chains[[1]][i,2],
                                                                                             processed_chains[[1]][i,3],
                                                                                             processed_chains[[1]][i,4]))
      subsampled_model_outputs[i, ] <- single_model_output
      
    }
    
    lower_credible_interval_processed <- apply(subsampled_model_outputs, MARGIN = 2, quantile, prob = 0.025)
    upper_credible_interval_processed <- apply(subsampled_model_outputs, MARGIN = 2, quantile, prob = 0.975)
    
    Lower_obs <- as.data.frame(lower_credible_interval_processed)
    Upper_obs <- as.data.frame(upper_credible_interval_processed)
    
    predicted_CrI <- cbind(Lower_obs, Upper_obs)
    predicted_CrI <- as.data.frame(predicted_CrI)
    predicted_CrI$age <- fitted_curve_df$age
    
  }
  
  # plot predicted prevalence curve & uncertainty band vs data
  p <- ggplot() +   
    theme(legend.position = 'none') +   
    geom_point(data=predicted_median_curve, aes(x=age, y=prev))+
    geom_errorbar(data=predicted_median_curve,aes(x=age, y=prev, ymin=lower, ymax=upper), width=0.8)+
    geom_line(data=predicted_median_curve,aes(x=age, y=predicted), size= 1.1, colour='purple')+
    geom_ribbon(data=predicted_CrI,aes(x=age, ymin=lower_credible_interval_processed,
                                       ymax=upper_credible_interval_processed), fill="purple", alpha=0.1)+
    ylim(0,1.0)+
    labs(x="Age (months) of human host", y="(Sero)prevalence (0 -1)")+
    theme_bw() +
    theme(strip.background=element_rect(fill=NA, color=NA),
          strip.text=element_text(size=11, face= "bold"),
          axis.text.x = element_text(size = 18),
          axis.text.y = element_text(size = 18),
          plot.title = element_text(size = 18, hjust=0.5),
          axis.title.x = element_text(size = 18, face= "bold"),
          axis.title.y = element_text(size = 16, angle = 90, face= "bold"),
          legend.position = c(0.83, 0.15),
          legend.title=element_text(size=20), 
          legend.text=element_text(size=16))
  
  p2 <- ggplot() +   
    theme(legend.position = 'none') +   
    geom_point(data=predicted_median_curve, aes(x=age, y=prev))+
    geom_errorbar(data=predicted_median_curve,aes(x=age, y=prev, ymin=lower, ymax=upper), width=0.8)+
    geom_line(data=predicted_median_curve,aes(x=age, y=predicted), size= 1.1, colour='purple')+
    geom_ribbon(data=predicted_CrI,aes(x=age, ymin=lower_credible_interval_processed,
                                       ymax=upper_credible_interval_processed), fill="purple", alpha=0.1)+
    ylim(0,0.5)+
    labs(x="Age (months) of human host", y="(Sero)prevalence (0 -1)")+
    theme_bw() +
    theme(strip.background=element_rect(fill=NA, color=NA),
          strip.text=element_text(size=11, face= "bold"),
          axis.text.x = element_text(size = 18),
          axis.text.y = element_text(size = 18),
          plot.title = element_text(size = 18, hjust=0.5),
          axis.title.x = element_text(size = 18, face= "bold"),
          axis.title.y = element_text(size = 16, angle = 90, face= "bold"),
          legend.position = c(0.83, 0.15),
          legend.title=element_text(size=20), 
          legend.text=element_text(size=16))
  
  
  p3 <- ggplot() +   
    theme(legend.position = 'none') +   
    geom_point(data=predicted_median_curve, aes(x=age, y=prev))+
    geom_errorbar(data=predicted_median_curve,aes(x=age, y=prev, ymin=lower, ymax=upper), width=0.8)+
    geom_line(data=predicted_median_curve,aes(x=age, y=predicted), size= 1.1, colour='purple')+
    geom_ribbon(data=predicted_CrI,aes(x=age, ymin=lower_credible_interval_processed,
                                       ymax=upper_credible_interval_processed), fill="purple", alpha=0.1)+
    ylim(0,0.25)+
    labs(x="Age (months) of human host", y="(Sero)prevalence (0 -1)")+
    theme_bw() +
    theme(strip.background=element_rect(fill=NA, color=NA),
          strip.text=element_text(size=11, face= "bold"),
          axis.text.x = element_text(size = 18),
          axis.text.y = element_text(size = 18),
          plot.title = element_text(size = 18, hjust=0.5),
          axis.title.x = element_text(size = 18, face= "bold"),
          axis.title.y = element_text(size = 16, angle = 90, face= "bold"),
          legend.position = c(0.83, 0.15),
          legend.title=element_text(size=20), 
          legend.text=element_text(size=16))
  
  
  return(list(p, p2, p3, predicted_median_curve, predicted_CrI))
  
}

#===============================================================================#
# produce predicted prevalence curves (simple model - fitted to single dataset) #


calculate_predicted_prevalence_multipledatasets_function1 <- function (max_age_toplot, data, pars, processed_chains, 
                                                                       model, number_datasets) {
  ## simple model ##
  if(model == "simple"){
    if(number_datasets == 5){
      
      # specify posterior median parameters to enable predicted prevalence calculation
      sp.median <- pars[[1]]
      se.median <- pars[[2]]
      lambda1.median <- pars[[3]]
      lambda2.median <- pars[[4]]
      lambda3.median <- pars[[5]]
      lambda4.median <- pars[[6]]
      lambda5.median <- pars[[7]]
      
      # set up dataframes for predicted prevalence calculations (n = 5 datasets)
      n <- length(unique(data$dataset))
      eval(parse(text = paste0("data", seq(1:n), " <- ", split(data, data$dataset))))
      eval(parse(text = paste0("data", seq(1:n), " <-  as.data.frame(data", seq(1:number_datasets), ")"))) # change seq to 1: number of individual datasets
      
      age_dum <- seq(from=0, to=max_age_toplot, by=0.005)  
      
      fitted_curve_df <- as.data.frame(age_dum) # make sequence of numbers (mean ages) for predicted variable
      
      names(fitted_curve_df)[names(fitted_curve_df)=="age_dum"] <- "age"
      
      predicted_median_curve1 <- full_join(fitted_curve_df, data1) 
      predicted_median_curve2 <- full_join(fitted_curve_df, data2) 
      predicted_median_curve3 <- full_join(fitted_curve_df, data3)
      predicted_median_curve4 <- full_join(fitted_curve_df, data4)
      predicted_median_curve5 <- full_join(fitted_curve_df, data5)
      
      # calculate predicted prevalence for each dataset
      predicted_median_curve1$predicted <- sapply(1:nrow(predicted_median_curve1), function(i) predicted_prev_multidataset_func2(age = predicted_median_curve1$age[i], 
                                                                                                                                 c(sp.median, se.median, lambda1.median)))
      predicted_median_curve2$predicted <- sapply(1:nrow(predicted_median_curve2), function(i) predicted_prev_multidataset_func2(age = predicted_median_curve2$age[i], 
                                                                                                                                 c(sp.median, se.median, lambda2.median)))
      predicted_median_curve3$predicted <- sapply(1:nrow(predicted_median_curve3), function(i) predicted_prev_multidataset_func2(age = predicted_median_curve3$age[i], 
                                                                                                                                 c(sp.median, se.median, lambda3.median)))
      predicted_median_curve4$predicted <- sapply(1:nrow(predicted_median_curve4), function(i) predicted_prev_multidataset_func2(age = predicted_median_curve4$age[i], 
                                                                                                                                 c(sp.median, se.median, lambda4.median)))
      predicted_median_curve5$predicted <- sapply(1:nrow(predicted_median_curve5), function(i) predicted_prev_multidataset_func2(age = predicted_median_curve5$age[i], 
                                                                                                                                 c(sp.median, se.median, lambda5.median)))
      # make master dataframe combining all datasets
      predicted_median_curve1$dataset_name <- rep(as.factor("data1"))
      predicted_median_curve2$dataset_name <- rep(as.factor("data2"))
      predicted_median_curve3$dataset_name <- rep(as.factor("data3"))
      predicted_median_curve4$dataset_name <- rep(as.factor("data4"))
      predicted_median_curve5$dataset_name <- rep(as.factor("data5"))
      
      
      predicted_median_curve <- rbind(predicted_median_curve1, predicted_median_curve2, predicted_median_curve3, predicted_median_curve4, predicted_median_curve5)
      
      # chains for each parameter to sample uncertainty from
      sp.chain <- c(PC_simple[[1]][,1], PC_simple[[2]][,1])
      se.chain <- c(PC_simple[[1]][,2], PC_simple[[2]][,2])
      lambda1.chain <- c(PC_simple[[1]][,3], PC_simple[[2]][,3])
      lambda2.chain <- c(PC_simple[[1]][,4], PC_simple[[2]][,4])
      lambda3.chain <- c(PC_simple[[1]][,5], PC_simple[[2]][,5])
      lambda4.chain <- c(PC_simple[[1]][,6], PC_simple[[2]][,6])
      lambda5.chain <- c(PC_simple[[1]][,7], PC_simple[[2]][,7])
      
      # calculate uncertainty (credible interval) of model run resulting from posterior #
      # dataset 1 #
      subsampled_model_outputs_dat1 <- matrix(NA, nrow = length(processed_chains[[1]][,1]), ncol = length(seq(0, 90, 0.005)))
      
      for (i in 1:length(processed_chains[[1]][,1])){
        
        single_model_output <- predicted_prev_multidataset_func2(seq(0, max_age_toplot, 0.005),c(sp.chain[i], se.chain[i], lambda1.chain[i]))
        subsampled_model_outputs_dat1[i, ] <- single_model_output
        
      }
      
      lower_credible_interval_processed <- apply(subsampled_model_outputs_dat1, MARGIN = 2, quantile, prob = 0.025)
      upper_credible_interval_processed <- apply(subsampled_model_outputs_dat1, MARGIN = 2, quantile, prob = 0.975)
      
      Lower_obs <- as.data.frame(lower_credible_interval_processed)
      Upper_obs <- as.data.frame(upper_credible_interval_processed)
      
      predicted_CrI_dat1 <- cbind(Lower_obs, Upper_obs)
      predicted_CrI_dat1 <- as.data.frame(predicted_CrI_dat1)
      predicted_CrI_dat1$age <- fitted_curve_df$age
      predicted_CrI_dat1$dataset_name <- rep(as.factor("data1"))
      names(predicted_CrI_dat1) <- c("lower", "upper", "age", "dataset_name")
      
      # dataset 2 #
      subsampled_model_outputs_dat2 <- matrix(NA, nrow = length(processed_chains[[1]][,1]), ncol = length(seq(0, 90, 0.005)))
      
      for (i in 1:length(processed_chains[[1]][,1])){
        
        single_model_output <- predicted_prev_multidataset_func2(seq(0, max_age_toplot, 0.005),c(sp.chain[i], se.chain[i], lambda2.chain[i]))
        subsampled_model_outputs_dat2[i, ] <- single_model_output
        
      }
      
      lower_credible_interval_processed <- apply(subsampled_model_outputs_dat2, MARGIN = 2, quantile, prob = 0.025)
      upper_credible_interval_processed <- apply(subsampled_model_outputs_dat2, MARGIN = 2, quantile, prob = 0.975)
      
      Lower_obs <- as.data.frame(lower_credible_interval_processed)
      Upper_obs <- as.data.frame(upper_credible_interval_processed)
      
      predicted_CrI_dat2 <- cbind(Lower_obs, Upper_obs)
      predicted_CrI_dat2 <- as.data.frame(predicted_CrI_dat2)
      predicted_CrI_dat2$age <- fitted_curve_df$age
      predicted_CrI_dat2$dataset_name <- rep(as.factor("data2"))
      names(predicted_CrI_dat2) <- c("lower", "upper", "age", "dataset_name")
      
      # dataset 3 #
      subsampled_model_outputs_dat3 <- matrix(NA, nrow = length(processed_chains[[1]][,1]), ncol = length(seq(0, 90, 0.005)))
      
      for (i in 1:length(processed_chains[[1]][,1])){
        
        single_model_output <- predicted_prev_multidataset_func2(seq(0, max_age_toplot, 0.005),c(sp.chain[i], se.chain[i], lambda3.chain[i]))
        subsampled_model_outputs_dat3[i, ] <- single_model_output
        
      }
      
      lower_credible_interval_processed <- apply(subsampled_model_outputs_dat3, MARGIN = 2, quantile, prob = 0.025)
      upper_credible_interval_processed <- apply(subsampled_model_outputs_dat3, MARGIN = 2, quantile, prob = 0.975)
      
      Lower_obs <- as.data.frame(lower_credible_interval_processed)
      Upper_obs <- as.data.frame(upper_credible_interval_processed)
      
      predicted_CrI_dat3 <- cbind(Lower_obs, Upper_obs)
      predicted_CrI_dat3 <- as.data.frame(predicted_CrI_dat3)
      predicted_CrI_dat3$age <- fitted_curve_df$age
      predicted_CrI_dat3$dataset_name <- rep(as.factor("data3"))
      names(predicted_CrI_dat3) <- c("lower", "upper", "age", "dataset_name")
      
      # dataset 4 #
      subsampled_model_outputs_dat4 <- matrix(NA, nrow = length(processed_chains[[1]][,1]), ncol = length(seq(0, 90, 0.005)))
      
      for (i in 1:length(processed_chains[[1]][,1])){
        
        single_model_output <- predicted_prev_multidataset_func2(seq(0, max_age_toplot, 0.005),c(sp.chain[i], se.chain[i], lambda4.chain[i]))
        subsampled_model_outputs_dat4[i, ] <- single_model_output
        
      }
      
      lower_credible_interval_processed <- apply(subsampled_model_outputs_dat4, MARGIN = 2, quantile, prob = 0.025)
      upper_credible_interval_processed <- apply(subsampled_model_outputs_dat4, MARGIN = 2, quantile, prob = 0.975)
      
      Lower_obs <- as.data.frame(lower_credible_interval_processed)
      Upper_obs <- as.data.frame(upper_credible_interval_processed)
      
      predicted_CrI_dat4 <- cbind(Lower_obs, Upper_obs)
      predicted_CrI_dat4 <- as.data.frame(predicted_CrI_dat4)
      predicted_CrI_dat4$age <- fitted_curve_df$age
      predicted_CrI_dat4$dataset_name <- rep(as.factor("data4"))
      names(predicted_CrI_dat4) <- c("lower", "upper", "age", "dataset_name")
      
      # dataset 5 #
      subsampled_model_outputs_dat5 <- matrix(NA, nrow = length(processed_chains[[1]][,1]), ncol = length(seq(0, 90, 0.005)))
      
      for (i in 1:length(processed_chains[[1]][,1])){
        
        single_model_output <- predicted_prev_multidataset_func2(seq(0, max_age_toplot, 0.005),c(sp.chain[i], se.chain[i], lambda5.chain[i]))
        subsampled_model_outputs_dat5[i, ] <- single_model_output
        
      }
      
      lower_credible_interval_processed <- apply(subsampled_model_outputs_dat5, MARGIN = 2, quantile, prob = 0.025)
      upper_credible_interval_processed <- apply(subsampled_model_outputs_dat5, MARGIN = 2, quantile, prob = 0.975)
      
      Lower_obs <- as.data.frame(lower_credible_interval_processed)
      Upper_obs <- as.data.frame(upper_credible_interval_processed)
      
      predicted_CrI_dat5 <- cbind(Lower_obs, Upper_obs)
      predicted_CrI_dat5 <- as.data.frame(predicted_CrI_dat5)
      predicted_CrI_dat5$age <- fitted_curve_df$age
      predicted_CrI_dat5$dataset_name <- rep(as.factor("data5"))
      names(predicted_CrI_dat5) <- c("lower", "upper", "age", "dataset_name")
      
      # combined uncertainty dataframes across datasets
      predicted_CrI <- rbind(predicted_CrI_dat1, predicted_CrI_dat2, predicted_CrI_dat3, predicted_CrI_dat4, predicted_CrI_dat5)
    }
  }
  
  ## reversible model ##
  if(model == "reversible"){
    if(number_datasets == 5){
      
      # specify posterior median parameters to enable predicted prevalence calculation
      sp.median <- pars[[1]]
      se.median <- pars[[2]]
      lambda1.median <- pars[[3]]
      lambda2.median <- pars[[4]]
      lambda3.median <- pars[[5]]
      lambda4.median <- pars[[6]]
      lambda5.median <- pars[[7]]
      rho1.median <- pars[[8]]
      rho2.median <- pars[[9]]
      rho3.median <- pars[[10]]
      rho4.median <- pars[[11]]
      rho5.median <- pars[[12]]
      
      # set up dataframes for predicted prevalence calculations (n = 5 datasets)
      n <- length(unique(data$dataset))
      eval(parse(text = paste0("data", seq(1:n), " <- ", split(data, data$dataset))))
      eval(parse(text = paste0("data", seq(1:n), " <-  as.data.frame(data", seq(1:number_datasets), ")"))) # change seq to 1: number of individual datasets
      
      age_dum <- seq(from=0, to=max_age_toplot, by=0.005)  
      
      fitted_curve_df <- as.data.frame(age_dum) # make sequence of numbers (mean ages) for predicted variable
      
      names(fitted_curve_df)[names(fitted_curve_df)=="age_dum"] <- "age"
      
      predicted_median_curve1 <- full_join(fitted_curve_df, data1) 
      predicted_median_curve2 <- full_join(fitted_curve_df, data2) 
      predicted_median_curve3 <- full_join(fitted_curve_df, data3)
      predicted_median_curve4 <- full_join(fitted_curve_df, data4)
      predicted_median_curve5 <- full_join(fitted_curve_df, data5)
      
      # calculate predicted prevalence for each dataset
      predicted_median_curve1$predicted <- sapply(1:nrow(predicted_median_curve1), function(i) predicted_prev_reversible_func2_multidataset(age = predicted_median_curve1$age[i], 
                                                                                                                                            c(sp.median, se.median, 
                                                                                                                                              lambda1.median, rho1.median)))
      predicted_median_curve2$predicted <- sapply(1:nrow(predicted_median_curve2), function(i) predicted_prev_reversible_func2_multidataset(age = predicted_median_curve2$age[i], 
                                                                                                                                            c(sp.median, se.median, 
                                                                                                                                              lambda2.median, rho2.median)))
      predicted_median_curve3$predicted <- sapply(1:nrow(predicted_median_curve3), function(i) predicted_prev_reversible_func2_multidataset(age = predicted_median_curve3$age[i], 
                                                                                                                                            c(sp.median, se.median, 
                                                                                                                                              lambda3.median, rho3.median)))
      predicted_median_curve4$predicted <- sapply(1:nrow(predicted_median_curve4), function(i) predicted_prev_reversible_func2_multidataset(age = predicted_median_curve4$age[i], 
                                                                                                                                            c(sp.median, se.median, 
                                                                                                                                              lambda4.median, rho4.median)))
      predicted_median_curve5$predicted <- sapply(1:nrow(predicted_median_curve5), function(i) predicted_prev_reversible_func2_multidataset(age = predicted_median_curve5$age[i], 
                                                                                                                                            c(sp.median, se.median, 
                                                                                                                                              lambda5.median, rho5.median)))
      # make master dataframe combining all datasets
      predicted_median_curve1$dataset_name <- rep(as.factor("data1"))
      predicted_median_curve2$dataset_name <- rep(as.factor("data2"))
      predicted_median_curve3$dataset_name <- rep(as.factor("data3"))
      predicted_median_curve4$dataset_name <- rep(as.factor("data4"))
      predicted_median_curve5$dataset_name <- rep(as.factor("data5"))
      
      
      predicted_median_curve <- rbind(predicted_median_curve1, predicted_median_curve2, predicted_median_curve3, predicted_median_curve4, predicted_median_curve5)
      
      # chains for each parameter to sample uncertainty from
      sp.chain <- c(PC_reversible[[1]][,1], PC_reversible[[2]][,1])
      se.chain <- c(PC_reversible[[1]][,2], PC_reversible[[2]][,2])
      lambda1.chain <- c(PC_reversible[[1]][,3], PC_reversible[[2]][,3])
      lambda2.chain <- c(PC_reversible[[1]][,4], PC_reversible[[2]][,4])
      lambda3.chain <- c(PC_reversible[[1]][,5], PC_reversible[[2]][,5])
      lambda4.chain <- c(PC_reversible[[1]][,6], PC_reversible[[2]][,6])
      lambda5.chain <- c(PC_reversible[[1]][,7], PC_reversible[[2]][,7])
      rho1.chain <- c(PC_reversible[[1]][,8], PC_reversible[[2]][,8])
      rho2.chain <- c(PC_reversible[[1]][,9], PC_reversible[[2]][,9])
      rho3.chain <- c(PC_reversible[[1]][,10], PC_reversible[[2]][,10])
      rho4.chain <- c(PC_reversible[[1]][,11], PC_reversible[[2]][,11])
      rho5.chain <- c(PC_reversible[[1]][,12], PC_reversible[[2]][,12])
      
      # calculate uncertainty (credible interval) of model run resulting from posterior #
      # dataset 1 #
      subsampled_model_outputs_dat1 <- matrix(NA, nrow = length(processed_chains[[1]][,1]), ncol = length(seq(0, 90, 0.005)))
      
      for (i in 1:length(processed_chains[[1]][,1])){
        
        single_model_output <- predicted_prev_reversible_func2_multidataset(seq(0, max_age_toplot, 0.005),c(sp.chain[i], se.chain[i], 
                                                                                                            lambda1.chain[i], rho1.chain[i]))
        subsampled_model_outputs_dat1[i, ] <- single_model_output
        
      }
      
      lower_credible_interval_processed <- apply(subsampled_model_outputs_dat1, MARGIN = 2, quantile, prob = 0.025)
      upper_credible_interval_processed <- apply(subsampled_model_outputs_dat1, MARGIN = 2, quantile, prob = 0.975)
      
      Lower_obs <- as.data.frame(lower_credible_interval_processed)
      Upper_obs <- as.data.frame(upper_credible_interval_processed)
      
      predicted_CrI_dat1 <- cbind(Lower_obs, Upper_obs)
      predicted_CrI_dat1 <- as.data.frame(predicted_CrI_dat1)
      predicted_CrI_dat1$age <- fitted_curve_df$age
      predicted_CrI_dat1$dataset_name <- rep(as.factor("data1"))
      names(predicted_CrI_dat1) <- c("lower", "upper", "age", "dataset_name")
      
      # dataset 2 #
      subsampled_model_outputs_dat2 <- matrix(NA, nrow = length(processed_chains[[1]][,1]), ncol = length(seq(0, 90, 0.005)))
      
      for (i in 1:length(processed_chains[[1]][,1])){
        
        single_model_output <- predicted_prev_reversible_func2_multidataset(seq(0, max_age_toplot, 0.005),c(sp.chain[i], se.chain[i], 
                                                                                                            lambda2.chain[i], rho2.chain[i]))
        subsampled_model_outputs_dat2[i, ] <- single_model_output
        
      }
      
      lower_credible_interval_processed <- apply(subsampled_model_outputs_dat2, MARGIN = 2, quantile, prob = 0.025)
      upper_credible_interval_processed <- apply(subsampled_model_outputs_dat2, MARGIN = 2, quantile, prob = 0.975)
      
      Lower_obs <- as.data.frame(lower_credible_interval_processed)
      Upper_obs <- as.data.frame(upper_credible_interval_processed)
      
      predicted_CrI_dat2 <- cbind(Lower_obs, Upper_obs)
      predicted_CrI_dat2 <- as.data.frame(predicted_CrI_dat2)
      predicted_CrI_dat2$age <- fitted_curve_df$age
      predicted_CrI_dat2$dataset_name <- rep(as.factor("data2"))
      names(predicted_CrI_dat2) <- c("lower", "upper", "age", "dataset_name")
      
      # dataset 3 #
      subsampled_model_outputs_dat3 <- matrix(NA, nrow = length(processed_chains[[1]][,1]), ncol = length(seq(0, 90, 0.005)))
      
      for (i in 1:length(processed_chains[[1]][,1])){
        
        single_model_output <- predicted_prev_reversible_func2_multidataset(seq(0, max_age_toplot, 0.005),c(sp.chain[i], se.chain[i], 
                                                                                                            lambda3.chain[i], rho3.chain[i]))
        subsampled_model_outputs_dat3[i, ] <- single_model_output
        
      }
      
      lower_credible_interval_processed <- apply(subsampled_model_outputs_dat3, MARGIN = 2, quantile, prob = 0.025)
      upper_credible_interval_processed <- apply(subsampled_model_outputs_dat3, MARGIN = 2, quantile, prob = 0.975)
      
      Lower_obs <- as.data.frame(lower_credible_interval_processed)
      Upper_obs <- as.data.frame(upper_credible_interval_processed)
      
      predicted_CrI_dat3 <- cbind(Lower_obs, Upper_obs)
      predicted_CrI_dat3 <- as.data.frame(predicted_CrI_dat3)
      predicted_CrI_dat3$age <- fitted_curve_df$age
      predicted_CrI_dat3$dataset_name <- rep(as.factor("data3"))
      names(predicted_CrI_dat3) <- c("lower", "upper", "age", "dataset_name")
      
      # dataset 4 #
      subsampled_model_outputs_dat4 <- matrix(NA, nrow = length(processed_chains[[1]][,1]), ncol = length(seq(0, 90, 0.005)))
      
      for (i in 1:length(processed_chains[[1]][,1])){
        
        single_model_output <- predicted_prev_reversible_func2_multidataset(seq(0, max_age_toplot, 0.005),c(sp.chain[i], se.chain[i], 
                                                                                                            lambda4.chain[i], rho4.chain[i]))
        subsampled_model_outputs_dat4[i, ] <- single_model_output
        
      }
      
      lower_credible_interval_processed <- apply(subsampled_model_outputs_dat4, MARGIN = 2, quantile, prob = 0.025)
      upper_credible_interval_processed <- apply(subsampled_model_outputs_dat4, MARGIN = 2, quantile, prob = 0.975)
      
      Lower_obs <- as.data.frame(lower_credible_interval_processed)
      Upper_obs <- as.data.frame(upper_credible_interval_processed)
      
      predicted_CrI_dat4 <- cbind(Lower_obs, Upper_obs)
      predicted_CrI_dat4 <- as.data.frame(predicted_CrI_dat4)
      predicted_CrI_dat4$age <- fitted_curve_df$age
      predicted_CrI_dat4$dataset_name <- rep(as.factor("data4"))
      names(predicted_CrI_dat4) <- c("lower", "upper", "age", "dataset_name")
      
      # dataset 5 #
      subsampled_model_outputs_dat5 <- matrix(NA, nrow = length(processed_chains[[1]][,1]), ncol = length(seq(0, 90, 0.005)))
      
      for (i in 1:length(processed_chains[[1]][,1])){
        
        single_model_output <- predicted_prev_reversible_func2_multidataset(seq(0, max_age_toplot, 0.005),c(sp.chain[i], se.chain[i], 
                                                                                                            lambda5.chain[i], rho5.chain[i]))
        subsampled_model_outputs_dat5[i, ] <- single_model_output
        
      }
      
      lower_credible_interval_processed <- apply(subsampled_model_outputs_dat5, MARGIN = 2, quantile, prob = 0.025)
      upper_credible_interval_processed <- apply(subsampled_model_outputs_dat5, MARGIN = 2, quantile, prob = 0.975)
      
      Lower_obs <- as.data.frame(lower_credible_interval_processed)
      Upper_obs <- as.data.frame(upper_credible_interval_processed)
      
      predicted_CrI_dat5 <- cbind(Lower_obs, Upper_obs)
      predicted_CrI_dat5 <- as.data.frame(predicted_CrI_dat5)
      predicted_CrI_dat5$age <- fitted_curve_df$age
      predicted_CrI_dat5$dataset_name <- rep(as.factor("data5"))
      names(predicted_CrI_dat5) <- c("lower", "upper", "age", "dataset_name")
      
      # combined uncertainty dataframes across datasets
      predicted_CrI <- rbind(predicted_CrI_dat1, predicted_CrI_dat2, predicted_CrI_dat3, predicted_CrI_dat4, predicted_CrI_dat5)
    }
  }
  
  
  # plot predicted prevalence curve & uncertainty band vs data
  p <- ggplot() +   
    theme(legend.position = 'none') +   
    geom_point(data=predicted_median_curve, aes(x=age, y=prev))+
    geom_errorbar(data=predicted_median_curve,aes(x=age, y=prev, ymin=lower, ymax=upper), width=0.8)+
    geom_line(data=predicted_median_curve,aes(x=age, y=predicted), size= 1.1, colour='purple')+
    geom_ribbon(data=predicted_CrI,aes(x=age, ymin=lower, ymax=upper), fill="purple", alpha=0.1)+
    ylim(0,1.0)+
    labs(x="Age (months) of human host", y="(Sero)prevalence (0 -1)")+
    facet_wrap(~dataset_name, scales = "free")+
    theme_bw() +
    theme(strip.background=element_rect(fill=NA, color=NA),
          strip.text=element_text(size=11, face= "bold"),
          axis.text.x = element_text(size = 18),
          axis.text.y = element_text(size = 18),
          plot.title = element_text(size = 18, hjust=0.5),
          axis.title.x = element_text(size = 18, face= "bold"),
          axis.title.y = element_text(size = 16, angle = 90, face= "bold"),
          legend.position = c(0.83, 0.15),
          legend.title=element_text(size=20), 
          legend.text=element_text(size=16))
  
  p2 <- ggplot() +   
    theme(legend.position = 'none') +   
    geom_point(data=predicted_median_curve, aes(x=age, y=prev))+
    geom_errorbar(data=predicted_median_curve,aes(x=age, y=prev, ymin=lower, ymax=upper), width=0.8)+
    geom_line(data=predicted_median_curve,aes(x=age, y=predicted), size= 1.1, colour='purple')+
    geom_ribbon(data=predicted_CrI,aes(x=age, ymin=lower, ymax=upper), fill="purple", alpha=0.1)+
    ylim(0,0.5)+
    labs(x="Age (months) of human host", y="(Sero)prevalence (0 -1)")+
    facet_wrap(~dataset_name, scales = "free")+
    theme_bw() +
    theme(strip.background=element_rect(fill=NA, color=NA),
          strip.text=element_text(size=11, face= "bold"),
          axis.text.x = element_text(size = 18),
          axis.text.y = element_text(size = 18),
          plot.title = element_text(size = 18, hjust=0.5),
          axis.title.x = element_text(size = 18, face= "bold"),
          axis.title.y = element_text(size = 16, angle = 90, face= "bold"),
          legend.position = c(0.83, 0.15),
          legend.title=element_text(size=20), 
          legend.text=element_text(size=16))
  
  
  p3 <- ggplot() +   
    theme(legend.position = 'none') +   
    geom_point(data=predicted_median_curve, aes(x=age, y=prev))+
    geom_errorbar(data=predicted_median_curve,aes(x=age, y=prev, ymin=lower, ymax=upper), width=0.8)+
    geom_line(data=predicted_median_curve,aes(x=age, y=predicted), size= 1.1, colour='purple')+
    geom_ribbon(data=predicted_CrI,aes(x=age, ymin=lower, ymax=upper), fill="purple", alpha=0.1)+
    ylim(0,0.25)+
    labs(x="Age (months) of human host", y="(Sero)prevalence (0 -1)")+
    facet_wrap(~dataset_name, scales = "free")+
    theme_bw() +
    theme(strip.background=element_rect(fill=NA, color=NA),
          strip.text=element_text(size=11, face= "bold"),
          axis.text.x = element_text(size = 18),
          axis.text.y = element_text(size = 18),
          plot.title = element_text(size = 18, hjust=0.5),
          axis.title.x = element_text(size = 18, face= "bold"),
          axis.title.y = element_text(size = 16, angle = 90, face= "bold"),
          legend.position = c(0.83, 0.15),
          legend.title=element_text(size=20), 
          legend.text=element_text(size=16))
  
  
  return(list(p, p2, p3, predicted_median_curve, predicted_CrI))
  
}

#=============================================================================================================================#
#                                      Process MCMC outputs functions                                                         #
#=============================================================================================================================#

# function to remove burnin
Burn <- function(chains, burnin){
  chains[-(1:burnin),]
}

# function to sub-sample
Downsample <- function(chains, sample){
  chains[seq(1, nrow(chains), sample),]
}

# function calling burnin and sub-sampling
Process_chains <- function(run1, run2, burnin, sample){
  C1 <- Burn(run1, burnin)
  C2<- Burn(run2, burnin)
  
  C1 <- Downsample(C1, sample)
  C2 <- Downsample(C2, sample)
  
  
  return(list(C1, C2))
}

# function to calculate autocorrelation for each parameter (simple model)
determine_autocorrelation_func1 <- function(processed_chain, number_datasets){
  
  if(number_datasets == 1){
    # Repplot auto/-correlation post processing #
    autocor.par1 <- cor(processed_chain[[1]][,1][-1], processed_chain[[1]][,1][-length(processed_chain[[1]][,1])])
    autocor.par2 <- cor(processed_chain[[1]][,2][-1], processed_chain[[1]][,2][-length(processed_chain[[1]][,2])])
    autocor.par3 <- cor(processed_chain[[1]][,3][-1], processed_chain[[1]][,3][-length(processed_chain[[1]][,3])])
    
    
    a <- signif(autocor.par1,3)
    b <- signif(autocor.par2,3)
    c <- signif(autocor.par3,3)
    
    return(list(a, b, c))
  }
  
  if(number_datasets == 5){
    # Repplot auto/-correlation post processing #
    autocor.par1 <- cor(processed_chain[[1]][,1][-1], processed_chain[[1]][,1][-length(processed_chain[[1]][,1])])
    autocor.par2 <- cor(processed_chain[[1]][,2][-1], processed_chain[[1]][,2][-length(processed_chain[[1]][,2])])
    autocor.par3 <- cor(processed_chain[[1]][,3][-1], processed_chain[[1]][,3][-length(processed_chain[[1]][,3])])
    autocor.par4 <- cor(processed_chain[[1]][,4][-1], processed_chain[[1]][,4][-length(processed_chain[[1]][,4])])
    autocor.par5 <- cor(processed_chain[[1]][,5][-1], processed_chain[[1]][,5][-length(processed_chain[[1]][,5])])
    autocor.par6 <- cor(processed_chain[[1]][,6][-1], processed_chain[[1]][,6][-length(processed_chain[[1]][,6])])
    autocor.par7 <- cor(processed_chain[[1]][,7][-1], processed_chain[[1]][,7][-length(processed_chain[[1]][,7])])
    
    a <- signif(autocor.par1,3)
    b <- signif(autocor.par2,3)
    c <- signif(autocor.par3,3)
    d <- signif(autocor.par4,3)
    e <- signif(autocor.par5,3)
    f <- signif(autocor.par6,3)
    g <- signif(autocor.par7,3)
    
    return(list(a, b, c, d, e, f, g))
  }
  
}

# function to calculate autocorrelation for each parameter (reversible model)
determine_autocorrelation_func2 <- function(processed_chain, number_datasets){
  
  if(number_datasets == 1){
    # Repplot auto/-correlation post processing #
    autocor.par1 <- cor(processed_chain[[1]][,1][-1], processed_chain[[1]][,1][-length(processed_chain[[1]][,1])])
    autocor.par2 <- cor(processed_chain[[1]][,2][-1], processed_chain[[1]][,2][-length(processed_chain[[1]][,2])])
    autocor.par3 <- cor(processed_chain[[1]][,3][-1], processed_chain[[1]][,3][-length(processed_chain[[1]][,3])])
    autocor.par4 <- cor(processed_chain[[1]][,4][-1], processed_chain[[1]][,4][-length(processed_chain[[1]][,4])])
    
    a <- signif(autocor.par1,3)
    b <- signif(autocor.par2,3)
    c <- signif(autocor.par3,3)
    d <- signif(autocor.par3,4)
    
    return(list(a, b, c, d))
  }
  
  if(number_datasets == 5){
    # Repplot auto/-correlation post processing #
    autocor.par1 <- cor(processed_chain[[1]][,1][-1], processed_chain[[1]][,1][-length(processed_chain[[1]][,1])])
    autocor.par2 <- cor(processed_chain[[1]][,2][-1], processed_chain[[1]][,2][-length(processed_chain[[1]][,2])])
    autocor.par3 <- cor(processed_chain[[1]][,3][-1], processed_chain[[1]][,3][-length(processed_chain[[1]][,3])])
    autocor.par4 <- cor(processed_chain[[1]][,4][-1], processed_chain[[1]][,4][-length(processed_chain[[1]][,4])])
    autocor.par5 <- cor(processed_chain[[1]][,5][-1], processed_chain[[1]][,5][-length(processed_chain[[1]][,5])])
    autocor.par6 <- cor(processed_chain[[1]][,6][-1], processed_chain[[1]][,6][-length(processed_chain[[1]][,6])])
    autocor.par7 <- cor(processed_chain[[1]][,7][-1], processed_chain[[1]][,7][-length(processed_chain[[1]][,7])])
    autocor.par8 <- cor(processed_chain[[1]][,8][-1], processed_chain[[1]][,8][-length(processed_chain[[1]][,8])])
    autocor.par9 <- cor(processed_chain[[1]][,9][-1], processed_chain[[1]][,9][-length(processed_chain[[1]][,9])])
    autocor.par10 <- cor(processed_chain[[1]][,10][-1], processed_chain[[1]][,10][-length(processed_chain[[1]][,10])])
    autocor.par11 <- cor(processed_chain[[1]][,11][-1], processed_chain[[1]][,11][-length(processed_chain[[1]][,11])])
    autocor.par12 <- cor(processed_chain[[1]][,12][-1], processed_chain[[1]][,12][-length(processed_chain[[1]][,12])])
    
    
    a <- signif(autocor.par1,3)
    b <- signif(autocor.par2,3)
    c <- signif(autocor.par3,3)
    d <- signif(autocor.par4,3)
    e <- signif(autocor.par5,3)
    f <- signif(autocor.par6,3)
    g <- signif(autocor.par7,3)
    h <- signif(autocor.par8,3)
    i <- signif(autocor.par9,3)
    j <- signif(autocor.par10,3)
    k <- signif(autocor.par11,3)
    l <- signif(autocor.par12,3)
    
    return(list(a, b, c, d, e, f, g, h, i, j, k, l))
  }
}


# function to obtain median and credible intervals for each parameter

obtain_parameter_values_func <- function(processed_chains, model, number_datasets){
  
  if(number_datasets == 1){
    lambda_simple<-quantile(c(processed_chains[[1]][,1], processed_chains[[2]][,1]), c(0.025,0.5,0.975))
    se_simple<-quantile(c(processed_chains[[1]][,2], processed_chains[[2]][,2]), c(0.025,0.5,0.975))
    sp_simple<-quantile(c(processed_chains[[1]][,3], processed_chains[[2]][,3]), c(0.025,0.5,0.975))
    
    lambda.median<-quantile(c(processed_chains[[1]][,1], processed_chains[[2]][,1]), c(0.5))
    se.median<-quantile(c(processed_chains[[1]][,2], processed_chains[[2]][,2]), c(0.5))
    sp.median<-quantile(c(processed_chains[[1]][,3], processed_chains[[2]][,3]), c(0.5))
    
    lambda.credible<-quantile(c(processed_chains[[1]][,1], processed_chains[[2]][,1]), c(0.025, 0.975))
    se.credible<-quantile(c(processed_chains[[1]][,2], processed_chains[[2]][,2]), c(0.025, 0.975))
    sp.credible<-quantile(c(processed_chains[[1]][,3], processed_chains[[2]][,3]), c(0.025, 0.975)) 
  }
  
  if(model == "simple" && number_datasets == 1){
    return(list(lambda.median, se.median, sp.median,
                lambda_simple, se_simple, sp_simple,
                lambda.credible, se.credible, sp.credible))
  }
  
  if(model == "reversible" && number_datasets == 1){
    rho_reversible <- quantile(c(processed_chains[[1]][,4], processed_chains[[2]][,4]), c(0.025,0.5,0.975))
    rho.median <- quantile(c(processed_chains[[1]][,4], processed_chains[[2]][,4]), c(0.5))
    rho.credible <- quantile(c(processed_chains[[1]][,4], processed_chains[[2]][,4]), c(0.025, 0.975))
    
    return(list(lambda.median, se.median, sp.median, rho.median,
                lambda_simple, se_simple, sp_simple, rho_reversible,
                lambda.credible, se.credible, sp.credible, rho.credible))
  }
  
  if(number_datasets == 5){
    
    sp_simple<-quantile(c(processed_chains[[1]][,1], processed_chains[[2]][,1]), c(0.025,0.5,0.975))
    se_simple<-quantile(c(processed_chains[[1]][,2], processed_chains[[2]][,2]), c(0.025,0.5,0.975))
    lambda1_simple<-quantile(c(processed_chains[[1]][,3], processed_chains[[2]][,3]), c(0.025,0.5,0.975))
    lambda2_simple<-quantile(c(processed_chains[[1]][,4], processed_chains[[2]][,4]), c(0.025,0.5,0.975))
    lambda3_simple<-quantile(c(processed_chains[[1]][,5], processed_chains[[2]][,5]), c(0.025,0.5,0.975))
    lambda4_simple<-quantile(c(processed_chains[[1]][,6], processed_chains[[2]][,6]), c(0.025,0.5,0.975))
    lambda5_simple<-quantile(c(processed_chains[[1]][,7], processed_chains[[2]][,7]), c(0.025,0.5,0.975))
    
    sp.median <-quantile(c(processed_chains[[1]][,1], processed_chains[[2]][,1]), c(0.5))
    se.median <-quantile(c(processed_chains[[1]][,2], processed_chains[[2]][,2]), c(0.5))
    lambda1.median <-quantile(c(processed_chains[[1]][,3], processed_chains[[2]][,3]), c(0.5))
    lambda2.median <-quantile(c(processed_chains[[1]][,4], processed_chains[[2]][,4]), c(0.5))
    lambda3.median <-quantile(c(processed_chains[[1]][,5], processed_chains[[2]][,5]), c(0.5))
    lambda4.median <-quantile(c(processed_chains[[1]][,6], processed_chains[[2]][,6]), c(0.5))
    lambda5.median <-quantile(c(processed_chains[[1]][,7], processed_chains[[2]][,7]), c(0.5))
    
    sp.credible <-quantile(c(processed_chains[[1]][,1], processed_chains[[2]][,1]), c(0.025,0.975))
    se.credible <-quantile(c(processed_chains[[1]][,2], processed_chains[[2]][,2]), c(0.025,0.975))
    lambda1.credible <-quantile(c(processed_chains[[1]][,3], processed_chains[[2]][,3]), c(0.025,0.975))
    lambda2.credible <-quantile(c(processed_chains[[1]][,4], processed_chains[[2]][,4]), c(0.025,0.975))
    lambda3.credible <-quantile(c(processed_chains[[1]][,5], processed_chains[[2]][,5]), c(0.025,0.975))
    lambda4.credible <-quantile(c(processed_chains[[1]][,6], processed_chains[[2]][,6]), c(0.025,0.975))
    lambda5.credible <-quantile(c(processed_chains[[1]][,7], processed_chains[[2]][,7]), c(0.025,0.975))
  }
  
  if(model == "simple" ){
    return(list(sp.median, se.median, lambda1.median, lambda2.median, lambda3.median, lambda4.median, lambda5.median,
                sp_simple, se_simple, lambda1_simple, lambda2_simple, lambda3_simple, lambda4_simple, lambda5_simple,
                sp.credible, se.credible, lambda1.credible, lambda2.credible, lambda3.credible, lambda4.credible, lambda5.credible))
  }
  
  if(model == "reversible" && number_datasets == 5){
    
    rho1_reversible <- quantile(c(processed_chains[[1]][,8], processed_chains[[2]][,8]), c(0.025,0.5,0.975))
    rho2_reversible <- quantile(c(processed_chains[[1]][,9], processed_chains[[2]][,9]), c(0.025,0.5,0.975))
    rho3_reversible <- quantile(c(processed_chains[[1]][,10], processed_chains[[2]][,10]), c(0.025,0.5,0.975))
    rho4_reversible <- quantile(c(processed_chains[[1]][,11], processed_chains[[2]][,11]), c(0.025,0.5,0.975))
    rho5_reversible <- quantile(c(processed_chains[[1]][,12], processed_chains[[2]][,12]), c(0.025,0.5,0.975))
    
    rho1.median <- quantile(c(processed_chains[[1]][,8], processed_chains[[2]][,8]), c(0.5))
    rho2.median <- quantile(c(processed_chains[[1]][,9], processed_chains[[2]][,9]), c(0.5))
    rho3.median <- quantile(c(processed_chains[[1]][,10], processed_chains[[2]][,10]), c(0.5))
    rho4.median <- quantile(c(processed_chains[[1]][,11], processed_chains[[2]][,11]), c(0.5))
    rho5.median <- quantile(c(processed_chains[[1]][,12], processed_chains[[2]][,12]), c(0.5))
    
    rho1.credible <- quantile(c(processed_chains[[1]][,8], processed_chains[[2]][,8]), c(0.025, 0.975))
    rho2.credible <- quantile(c(processed_chains[[1]][,9], processed_chains[[2]][,9]), c(0.025, 0.975))
    rho3.credible <- quantile(c(processed_chains[[1]][,10], processed_chains[[2]][,10]), c(0.025, 0.975))
    rho4.credible <- quantile(c(processed_chains[[1]][,11], processed_chains[[2]][,11]), c(0.025, 0.975))
    rho5.credible <- quantile(c(processed_chains[[1]][,12], processed_chains[[2]][,12]), c(0.025, 0.975))
    
    return(list(sp.median, se.median, lambda1.median, lambda2.median, lambda3.median, lambda4.median, lambda5.median, rho1.median,
                rho2.median, rho3.median, rho4.median, rho5.median,
                sp_simple, se_simple, lambda1_simple, lambda2_simple, lambda3_simple, lambda4_simple, lambda5_simple, rho1_reversible,
                rho2_reversible, rho3_reversible, rho4_reversible, rho5_reversible,
                sp.credible, se.credible, lambda1.credible, lambda2.credible, lambda3.credible, lambda4.credible, lambda5.credible,
                rho1.credible, rho2.credible, rho3.credible, rho4.credible, rho5.credible))
    
  }
}

# calculate model fit : deviance information criterion (DIC) statistic #
# for simple model
calculate_DIC_func1 <- function(chain, burnin, subsample, parameters, number_datasets){
  
  #remove burnin from loglik and logposterior chains
  loglike.chain.no.burnin_model1 <- chain$Likelihood[-(1:burnin)]
  logpost.chain.no.burnin_model1 <- chain$Posterior_Output[-(1:burnin)]
  
  ## thinning 
  loglike.chain.sub_model1 <- loglike.chain.no.burnin_model1[seq(1,length(loglike.chain.no.burnin_model1),subsample)]
  logpost.chain.sub_model1 <- logpost.chain.no.burnin_model1[seq(1,length(logpost.chain.no.burnin_model1),subsample)]
  
  ## computing goodness of fit (mean deviance) - deviance for given theta --> want to calculate mean of these 
  D_bar_model1 <- mean(-2 * loglike.chain.sub_model1) # mean deviance
  
  
  if(number_datasets == 1){
    # Posterior log likelihood
    modal_posterior_likelihood <- posterior_function_simple(data=data, c(parameters[[1]], 
                                                                         parameters[[2]], 
                                                                         parameters[[3]]))
  }
  
  if(number_datasets == 5){
    # Posterior log likelihood
    modal_posterior_likelihood <- posterior_function_simple_multidata(data=data, c(parameters[[1]], parameters[[2]], 
                                                                                   parameters[[3]], parameters[[4]], parameters[[5]],
                                                                                   parameters[[6]], parameters[[7]]))
  }
  
  modal_posterior_deviance <- -2 * modal_posterior_likelihood
  
  # calculating the DIC 
  DIC <- 2 * D_bar_model1 - modal_posterior_deviance
  
  return(list(D_bar_model1, modal_posterior_likelihood, modal_posterior_deviance, DIC))
  
}

# for reversible model
calculate_DIC_func2 <- function(chain, burnin, subsample, parameters, number_datasets, simple_lambda_median){
  
  #remove burnin from loglik and logposterior chains
  loglike.chain.no.burnin_model1 <- chain$Likelihood[-(1:burnin)]
  logpost.chain.no.burnin_model1 <- chain$Posterior_Output[-(1:burnin)]
  
  ## thinning 
  loglike.chain.sub_model1 <- loglike.chain.no.burnin_model1[seq(1,length(loglike.chain.no.burnin_model1),subsample)]
  logpost.chain.sub_model1 <- logpost.chain.no.burnin_model1[seq(1,length(logpost.chain.no.burnin_model1),subsample)]
  
  ## computing goodness of fit (mean deviance) - deviance for given theta --> want to calculate mean of these 
  D_bar_model1 <- mean(-2 * loglike.chain.sub_model1) # mean deviance
  
  
  if(number_datasets == 1){
    # Posterior log likelihood
    modal_posterior_likelihood <- posterior_function_reversible(data=data, c(parameters[[1]], 
                                                                             parameters[[2]], 
                                                                             parameters[[3]],
                                                                             parameters[[4]]), simple_lambda_median)
  }
  
  if(number_datasets == 5){
    # Posterior log likelihood
    modal_posterior_likelihood <- posterior_function_reversible_multidata(data=data, c(parameters[[1]], parameters[[2]], parameters[[3]], 
                                                                                       parameters[[4]], parameters[[5]], parameters[[6]], 
                                                                                       parameters[[7]], parameters[[8]], parameters[[9]], 
                                                                                       parameters[[10]], parameters[[11]], parameters[[12]]),
                                                                          simple_lambda_median)
  }
  
  
  modal_posterior_deviance <- -2 * modal_posterior_likelihood
  
  # calculating the DIC 
  DIC <- 2 * D_bar_model1 - modal_posterior_deviance
  
  return(list(D_bar_model1, modal_posterior_likelihood, modal_posterior_deviance, DIC))
  
}


