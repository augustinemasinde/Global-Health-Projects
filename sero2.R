# -------------------------------
# 1. Load libraries
# -------------------------------
library(devtools)    # if needed to install serosolver
library(serosolver)
library(tidyverse)
library(lubridate)
library(data.table)
library(plyr)
library(foreach)
library(doParallel)
library(readxl)

# -------------------------------
# 2. Constants & seed
# -------------------------------
set.seed(123)
no_chains <- 5
prior_version <- 2  # version 2 integrates over infection histories

# -------------------------------
# 3. Load & clean data
# -------------------------------
chikdata <- read_excel("~/Downloads/chikungunya_data_Uganda.xlsx")

# Fix typos and missing values in IgM_CHIK column
chikdata$IgM_CHIK[chikdata$IgM_CHIK == "Nengative"] <- "Negative"
chikdata$IgM_CHIK[chikdata$IgM_CHIK == "NA"] <- NA
chikdata <- chikdata[!is.na(chikdata$IgM_CHIK), ]

# Simulate titres based on IgM_CHIK result
chikdata$titre <- NA_real_
chikdata$titre[chikdata$IgM_CHIK == "Negative"] <- runif(sum(chikdata$IgM_CHIK == "Negative"), 0, 10)
chikdata$titre[chikdata$IgM_CHIK == "Positive"] <- rlnorm(sum(chikdata$IgM_CHIK == "Positive"), meanlog = log(40), sdlog = 0.4)

# Clean gender
chikdata$Gender[chikdata$Gender == "2-Female"] <- "Female"
chikdata$Gender[chikdata$Gender == "1-Male"] <- "Male"
chikdata <- chikdata[!is.na(chikdata$Gender), ]

# Generate synthetic birth dates from age
generate_birthday_date <- function(age_years) {
  month <- sample(1:12, 1)
  day <- sample(1:ifelse(month %in% c(4,6,9,11), 30, ifelse(month == 2, 28, 31)), 1)
  birth_year <- year(Sys.Date()) - age_years
  make_date(year = birth_year, month = month, day = day)
}
chikdata$birth_date <- sapply(chikdata$Age_Yrs, generate_birthday_date) %>% as.Date()

# Convert date_case_reported to Date and calculate sample_time as days since first sample
chikdata$date_case_reported <- as.Date(chikdata$date_case_reported)
origin_date <- min(chikdata$date_case_reported, na.rm = TRUE)
chikdata$sample_time <- as.numeric(chikdata$date_case_reported - origin_date)

# -------------------------------
# 4. Format titre data for serosolver
# -------------------------------
titre_dat <- chikdata %>%
  mutate(
    individual = UniqueKey,
    biomarker_id = 1,
    repeat_number = 1
  ) %>%
  select(individual, sample_time, biomarker_id, repeat_number, titre) %>%
  as.data.frame()

# antibody_data is the same structure here
antibody_data <- titre_dat

# -------------------------------
# 5. Define strain isolation times in days
# -------------------------------
resolution <- 1
strain_isolation_times <- seq(
  from = 0,
  to = max(titre_dat$sample_time, na.rm = TRUE),
  by = resolution
)

# -------------------------------
# 6. Load and prepare parameter table
# -------------------------------
par_tab_path <- system.file("extdata", "par_tab_base.csv", package = "serosolver")
par_tab <- read.csv(par_tab_path, stringsAsFactors = FALSE)

# Set beta prior on attack rate parameters alpha and beta
par_tab[par_tab$names %in% c("alpha", "beta"), "values"] <- c(1/3, 1/3)

# Set maximum titre to 9 based on simulated data scale
par_tab[par_tab$names == "MAX_TITRE", "values"] <- 9

# Remove phi if using integrated infection histories (prior_version = 2)
par_tab <- par_tab[par_tab$names != "phi", ]

# Fix cross-reactivity and antigenic seniority parameters (set fixed and values to 0)
par_tab[par_tab$names %in% c("tau", "sigma1", "sigma2"), "fixed"] <- 1
par_tab[par_tab$names %in% c("tau", "sigma1", "sigma2"), "values"] <- 0

# Ensure mu_short has a valid lower bound
par_tab[par_tab$names == "mu_short", "lower_bound"] <- 1

# -------------------------------
# 7. Create posterior function
# -------------------------------
model_func <- create_posterior_func(
  par_tab = par_tab,
  titre_dat = titre_dat,
  antibody_data = antibody_data,
  strain_isolation_times = strain_isolation_times,
  version = prior_version
)

# -------------------------------
# 8. Run MCMC chains in parallel
# -------------------------------
registerDoParallel(cores = no_chains)

chain_path_real <- "chik_output/cs1_real/"
dir.create(chain_path_real, recursive = TRUE, showWarnings = FALSE)

filenames <- paste0("chik_chain_", 1:no_chains)

res <- foreach(x = filenames, .packages = c("serosolver", "data.table", "plyr")) %dopar% {
  
  start_prob <- -Inf
  
  while (!is.finite(start_prob)) {
    
    start_tab <- generate_start_tab(par_tab)
    
    start_inf <- setup_infection_histories(
      titre_dat = titre_dat,
      strain_isolation_times = strain_isolation_times,
      space = 3,
      titre_cutoff = 4
    )
    
    start_prob <- sum(model_func(start_tab$values, start_inf)[[1]])
  }
  
  serosolver(
    par_tab = start_tab,
    titre_dat = titre_dat,
    antigenic_map = NULL,
    strain_isolation_times = strain_isolation_times,
    start_inf_hist = start_inf,
    mcmc_pars = c(
      iterations = 2000000,
      target_acceptance_rate_theta = 0.44,
      target_acceptance_rate_inf_hist = 0.44,
      adaptive_frequency = 1000,
      thin = 1000,
      adaptive_iterations = 500000,
      save_block = 1000,
      thin_inf_hist = 100,
      proposal_inf_hist_indiv_prop = 1,
      proposal_ratio = 2,
      burnin = 0,
      proposal_inf_hist_time_prop = 0.5,
      proposal_inf_hist_distance = 3,
      proposal_inf_hist_adaptive = 1,
      proposal_inf_hist_indiv_swap_ratio = 0.5,
      proposal_inf_hist_group_swap_ratio = 0.5,
      proposal_inf_hist_group_swap_prop = 1
    ),
    filename = paste0(chain_path_real, x),
    CREATE_POSTERIOR_FUNC = create_posterior_func,
    version = prior_version
  )
}

# -------------------------------
# 9. End
# -------------------------------
cat("MCMC chains started.\n")
