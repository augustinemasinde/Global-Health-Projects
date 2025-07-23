library(serosolver)
library(ggplot2)
library(coda)
library(plyr)
library(reshape2)
library(data.table)
library(tidyr)
library(doParallel)
library(foreach)
library(ggpubr)
library(bayesplot)
library(viridis)
library(ggthemes)
library(cowplot)
library(grid)
library(gridExtra)
library(doRNG)
library(serosim)
library(tidyverse)
library(readxl)
set.seed(0)

# Setup -------------------------------------------------------------------
serosolver <- FALSE

## Filename prefix for all saved outputs
filename <- "chikungunya_data_test"
filenames <- paste0(filename, "_",1:5)

## We'll be parallelising a few chains
cl <- makeCluster(detectCores(), type='PSOCK')
registerDoParallel(cl)
registerDoMC(cores=5)
stopCluster(cl)
foreach(i = 1:5) %dopar% { Sys.getpid() }


#import data
chikdata <-read_excel("~/Desktop/my files/chikungunya_data_Uganda.xlsx")

set.seed(123)  # For reproducibility

chikdata$titre <- NA_real_

# Assign low titres to seronegative individuals
chikdata$titre[chikdata$IgM_CHIK == "Negative"] <- runif(
  sum(chikdata$IgM_CHIK == "Negative"), min = 0, max = 10
)

# Assign higher titres to seropositive individuals
chikdata$titre[chikdata$IgM_CHIK == "Positive"] <- rlnorm(
  sum(chikdata$IgM_CHIK == "Positive"), meanlog = log(40), sdlog = 0.4
)

#select relevant columns
chikdata<- chikdata %>%  select(UniqueKey, Age_Yrs, Year,IgM_CHIK, titre)

#A function to generate date of birth from years
generate_birthday_date <- function(age_years) {
  month <- sample(1:12, 1)
  day <- switch(month,
                "1"=sample(1:31, 1),
                "2"=sample(1:28, 1),
                "3"=sample(1:31, 1),
                "4"=sample(1:30, 1),
                "5"=sample(1:31, 1),
                "6"=sample(1:30, 1),
                "7"=sample(1:31, 1),
                "8"=sample(1:31, 1),
                "9"=sample(1:30, 1),
                "10"=sample(1:31, 1),
                "11"=sample(1:30, 1),
                "12"=sample(1:31, 1))
  
  birth_year <- year(Sys.Date()) - age_years
  birth_date <- make_date(year = birth_year, month = month, day = day)
  
  return(birth_date)
}
generate_birthday_date_vec <- Vectorize(generate_birthday_date)
chikdata$DOB <- as.Date(generate_birthday_date_vec(chikdata$Age_Yrs))

chikdata$birth_year <- as.numeric(format(chikdata$DOB, "%Y"))

#Data cleaning
chikdata$birth_year[chikdata$birth_year == "1891"] <- "1991"
# Ensure DOB is a Date
chikdata$DOB <- as.Date(chikdata$DOB)

#recode Nengative to Negative
chikdata$IgM_CHIK[chikdata$IgM_CHIK == "Nengative"] <- "Negative"

# Recode "NA" string to actual NA
chikdata$IgM_CHIK[chikdata$IgM_CHIK == "NA"] <- NA

# Drop all real NA values
chikdata <- chikdata[!is.na(chikdata$IgM_CHIK), ]

#Data formatting for serosolver
chikdata <- chikdata %>% select(UniqueKey,titre, DOB, Year)
chikdata$virus<- c(rep(8037,nrow(chikdata)))
chikdata$samples<- c(rep(2021,nrow(chikdata)))
chikdata$run<- c(rep(1,nrow(chikdata)))
chikdata$group <- c(rep(1,nrow(chikdata)))
chikdata$virus <- c(rep(2021,nrow(chikdata)))
chikdata <- chikdata %>% dplyr::rename(individual = UniqueKey)
chikdata$titre <- log2(chikdata$titre)
chikdata <- chikdata %>% select(individual,virus,titre,samples,DOB,run, Year)
chikdata <- chikdata[chikdata$Year == 2021, ]
chikdata$DOB <- year(chikdata$DOB)
chikdata$group <- c(rep(1, nrow(chikdata)))
chikdata <- chikdata %>%  select(individual, DOB, virus, titre, samples, run, group)
chikdata$individual <- 1:nrow(chikdata)
chikdata <- chikdata[!is.na(chikdata$DOB), ]
chikdata$individual <- as.integer(chikdata$individual)
chikdata$run <- as.integer(chikdata$run)
chikdata <- chikdata[rep(1:nrow(chikdata), each = 2), ]
chikdata$virus <- rep(c(2019, 2020), length.out = nrow(chikdata))
chikdata$biomarker_group <- c(rep(1,nrow(chikdata)))

set.seed(123)  # for reproducibility

chikdata$titre <- sapply(chikdata$virus, function(year) {
  if (year == 2019) {
    val <- rnorm(1, mean = 3.5, sd = 1.2)  # lower mean titre
  } else if (year == 2020) {
    val <- rnorm(1, mean = 5.0, sd = 1.2)  # higher mean titre
  } else {
    val <- rnorm(1, mean = 4.0, sd = 1.2)  # fallback
  }
  val <- max(min(val, 7), 0)  # clamp between 0 and 7
  return(val)
})

filename <- "serosurvey_2"
resolution <- 1 ## eg. this would be set to 12 for monthly resolution
sample_year <- 2021

serosolver::describe_priors()
#> Which version to use in run_MCMC? The following text describes the proposal step for updating infection histories.
#> Version 1: Beta prior on per time attack rates. Explicit FOI on each epoch using probability of infection term. Proposal performs N `flip` proposals at random locations in an individual's infection history, switching 1->0 or 0->1. Otherwise, swaps the contents of two random locations
#> Version 2: Beta prior on per time attack rates. Gibbs sampling of infection histories as in Indian Buffet Process papers, integrating out each probability of infection term.
#> Version 3: Beta prior on probability of infection for an individual, assuming independence between individuals. Samples from a beta binomial with alpha and beta specified by the par_tab input. Proposes nInfs moves at a time for add/remove, or when swapping, swaps locations up to moveSize time steps away
#> Version 4: Beta prior on probability of any infection. Gibbs sampling of infection histories using total number of infections across all times and all individuals as the prior
prior_version <- 2


par_tab_path <- system.file("extdata", "par_tab_base.csv", package = "serosolver")
par_tab <- read.csv(file = par_tab_path, stringsAsFactors=FALSE)


## Set parameters for Beta prior on infection histories
beta_pars <- find_beta_prior_mode(0.15,4)
par_tab[par_tab$names == "alpha","values"] <- beta_pars$alpha
par_tab[par_tab$names == "beta","values"] <- beta_pars$beta
## Maximum recordable log titre in these data is 8
par_tab[par_tab$names == "MAX_TITRE","values"] <- 8

## Remove phi parameters, as these are integrated out under prior version 2
par_tab <- par_tab[par_tab$names != "phi",]

## Fix all short term parameters to 0
par_tab[par_tab$names %in% c("mu_short","sigma2","wane"),"fixed"] <- 1 # mu_short, waning and sigma2 are fixed
par_tab[par_tab$names %in% c("mu_short","sigma2","wane"),"values"] <- 0 # set these values to 0


## Distinct filename for each chain
no_chains <- 5
filenames <- paste0(filename, "_",1:no_chains)
chain_path <- sub("par_tab_base.csv","",par_tab_path)
chain_path_real <- paste0(chain_path, "cs2_real/")
chain_path_sim <- paste0(chain_path, "cs2_sim/")

#abtibody data
antibody_data <- chikdata %>%  select(individual, samples, DOB, titre)
antibody_data$biomarker_group <- c(rep(1, nrow(antibody_data)))
antibody_data$biomarker_id<- c(rep(2021, nrow(antibody_data)))
antibody_data$repeat_number<- c(rep(1, nrow(antibody_data)))
antibody_data$population_group <- c(rep(1, nrow(antibody_data)))
antibody_data$birth<- c(rep(2021, nrow(antibody_data)))
antibody_data$biomarker_group<- c(rep(1, nrow(antibody_data)))
antibody_data <- antibody_data %>%  rename(measurement = titre)
antibody_data <- antibody_data %>%  rename(sample_time = samples)
antibody_data<- antibody_data %>%  select(individual, sample_time, biomarker_id,biomarker_group, measurement,repeat_number,population_group,birth)


## Distinct filename for each chain
no_chains <- 5
filenames <- paste0(filename, "_",1:no_chains)
chain_path <- sub("par_tab_base.csv","",par_tab_path)
chain_path_real <- paste0(chain_path, "cs2_real/")
chain_path_sim <- paste0(chain_path, "cs2_sim/")



## Read in raw coordinates
antigenic_coords <- read_csv("~/Desktop/my files/antigenic coord.csv")
print(head(antigenic_coords))

#>       x_coord   y_coord inf_times
#> 1 -0.09718111 0.5021363      1968
#> 2  0.80502804 1.6816917      1969
#> 3  1.70723718 2.8612472      1970
#> 4  2.60944633 3.9976902      1971
#> 5  3.51165548 4.8709093      1972
#> 6  4.41386463 5.5057039      1973

## More flexible version of the above function
virus_key <- c("EAL" = 2019, "SAL" = 2020, "EAL1" = 2019, "SAL1" = 2020)
antigenic_coords$Strain <- virus_key[antigenic_coords$Strain]
antigenic_map <- generate_antigenic_map_flexible(antigenic_coords,buckets = 1, spar = 0.3)

## Restrict entries to years of interest. Entries in antigenic_map determine
## the times that individual can be infected ie. the dimensions of the infection
## history matrix.
antigenic_map <- antigenic_map[antigenic_map$inf_times >= 2019 & antigenic_map$inf_times <= sample_year,]
strain_isolation_times <- unique(antigenic_map$inf_times)

## Create the posterior solving function that will be used in the MCMC framework 
model_func <- create_posterior_func(par_tab=par_tab,
                                    antibody_data = antibody_data,
                                    titre_dat=chikdata,
                                    antigenic_map= antigenic_map,
                                    version=prior_version) # function in posteriors.R
#> Creating posterior solving function...
#> 
## Generate results in parallel
res <- foreach(x = filenames, .packages = c('serosolver','data.table','plyr')) %dopar% {
  ## Not all random starting conditions return finite likelihood, so for each chain generate random
  ## conditions until we get one with a finite likelihood
  start_prob <- -Inf
  while(!is.finite(start_prob)){
    ## Generating starting antibody kinetics parameters
    start_tab <- generate_start_tab(par_tab)
    
    ## Generate starting infection history
    start_inf <- setup_infection_histories_titre(chikdata, strain_isolation_times, 
                                                 space=3,titre_cutoff=4)
    start_prob <- sum(model_func(start_tab$values, start_inf)[[1]])
  }
  
  
  res <- serosolver(par_tab = start_tab, 
                    titre_dat = chikdata,
                    antigenic_map = antigenic_map,
                    start_inf_hist = start_inf, 
                    mcmc_pars = c("iterations"=500000,"adaptive_iterations"=100000,"thin"=1000,
                                  "thin_inf_hist"=1000,"save_block"=1000,
                                  "proposal_inf_hist_time_prop"=1, "proposal_inf_hist_indiv_prop"=1,
                                  "proposal_inf_hist_group_swap_ratio"=0.8, "proposal_inf_hist_group_swap_prop"=1),
                    filename = paste0(chain_path_real,x), 
                    CREATE_POSTERIOR_FUNC = create_posterior_func, 
                    version = prior_version)
}












