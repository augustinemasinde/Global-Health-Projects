library(serosim)
library(tidyverse)
library(data.table)
library(ggplot2)
library(patchwork)
library(reshape2)

##import the data
library(readxl)
chikdata <- read_excel("~/Desktop/my files/chikungunya_data_Uganda.xlsx")



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

chikdata<- chikdata %>%  select(UniqueKey,Gender, Age_Yrs, Year,IgM_CHIK, titre)


###data cleaning
chikdata$IgM_CHIK[chikdata$IgM_CHIK == "Nengative"] <- "Negative"

# Recode "NA" string to actual NA
chikdata$IgM_CHIK[chikdata$IgM_CHIK == "NA"] <- NA

# Drop all real NA values
chikdata <- chikdata[!is.na(chikdata$IgM_CHIK), ]


# Recode gender values
chikdata$Gender[chikdata$Gender == "2-Female"] <- "Female"
chikdata$Gender[chikdata$Gender == "1-Male"] <- "Male"

# Drop rows with NA in Gender
chikdata <- chikdata[!is.na(chikdata$Gender), ]

# Confirm result
unique(chikdata$Gender)


library(lubridate)

set.seed(123)  # For reproducibility

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
chikdata$birth_date <- as.Date(generate_birthday_date_vec(chikdata$Age_Yrs))


# Load required packages
library(tibble)
library(dplyr)

# STEP 1: Filter your chikdata
chik_filtered <- chikdata %>%
  filter(!is.na(birth_date), !is.na(Gender))

# STEP 2: Set N and times
N <- nrow(chik_filtered)
times <- 0:100

# STEP 3: Convert birth_date to simulation time units
# We'll assume simulation ends at max birth year
birth_years <- as.numeric(format(chik_filtered$birth_date, "%Y"))
sim_end_year <- max(birth_years, na.rm = TRUE)
birth_times <- sim_end_year - birth_years  # Now on simulation time scale

# STEP 4: Define aux (e.g., Gender distribution)
gender_table <- table(chik_filtered$Gender)
gender_props <- prop.table(gender_table)
gender_values <- names(gender_table)
gender_probs <- as.numeric(gender_props)

aux <- list(
  list(
    var = "Gender",
    values = gender_values,
    probs = gender_probs
  )
)

# STEP 5: Patched version of generate_pop_demography
generate_pop_demography <- function(
    N,
    times,
    birth_times = NULL,
    age_min = 0,
    removal_min = min(times),
    removal_max = max(times),
    prob_removal,
    aux = NULL
) {
  # Generate default birth_times if not provided
  if (is.null(birth_times)) {
    birth_times <- sample(times[1:(length(times) - age_min)], N, replace = TRUE)
  }
  
  # Removal time logic
  removal <- if (is.na(prob_removal)) {
    rep(max(times) + 1, N)
  } else {
    removal <- rep(max(times) + 1, N)
    remove_i <- sample(1:N, size = floor(prob_removal * N), replace = FALSE)
    removal_ages <- sample(removal_min:removal_max, length(remove_i), replace = TRUE)
    removal[remove_i] <- birth_times[remove_i] + removal_ages
    removal
  }
  
  # Create basic tibble
  demog <- tibble(
    i = seq_len(N),
    birth = birth_times,
    removal = removal
  )
  
  # Handle aux variables safely
  if (!is.null(aux)) {
    for (j in seq_along(aux)) {
      values <- aux[[j]]$values
      probs <- aux[[j]]$probs
      
      # Safety checks
      if (length(values) == 0 || length(probs) == 0) {
        stop(paste("Aux variable", aux[[j]]$var, "has invalid values or probs"))
      }
      if (length(values) != length(probs)) {
        stop(paste("Length mismatch in aux variable", aux[[j]]$var))
      }
      
      probs <- probs / sum(probs)  # Normalize
      demog[[aux[[j]]$var]] <- sample(values, N, replace = TRUE, prob = probs)
    }
  }
  
  return(demog)
}

# STEP 6: Run it
demography <- generate_pop_demography(
  N = N,
  times = times,
  birth_times = birth_times,
  age_min = 0,
  removal_min = 0,
  removal_max = max(times),
  prob_removal = 0.2,
  aux = aux
)

library(ggplot2)
library(dplyr)

# Calculate age at removal or current age if not removed
demography <- demography %>%
  mutate(
    age_at_removal = removal - birth,
    age_at_removal = ifelse(age_at_removal > max(times), NA, age_at_removal)
  )

# 1. Age structure (distribution of birth times / ages)
ggplot(demography, aes(x = birth)) +
  geom_histogram(binwidth = 5, fill = "steelblue", color = "black") +
  labs(title = "Age Structure (Birth Times)", x = "Birth Time (simulation units)", y = "Count")

# 2. Removal timing distribution
ggplot(demography, aes(x = age_at_removal)) +
  geom_histogram(binwidth = 5, fill = "tomato", color = "black", na.rm = TRUE) +
  labs(title = "Removal Timing Distribution", x = "Age at Removal", y = "Count")

# 3. Gender distribution
ggplot(demography, aes(x = Gender, fill = Gender)) +
  geom_bar() +
  labs(title = "Gender Distribution", y = "Count") +
  theme(legend.position = "none")




#catalytic models-constant
library(serofoi)
data("chagas2012")
serosurvey<- chagas2012



# how the disease circulates
foi_df <- data.frame(
  year = seq(2000, 2049, 1),
  foi = rep(0.02, 50)
)

# specify 5-year age bins for observations
survey_features <- data.frame(
  age_min = seq(1, 50, 5),
  age_max = seq(5, 50, 5),
  n_sample = rep(25, 10)
)

serosurvey_constant <- simulate_serosurvey(
  "time",
  foi_df,
  survey_features
) |>
  mutate(survey_year = 2050)

seromodel_constant <- fit_seromodel(
  serosurvey = serosurvey_constant,
  model_type = "constant",
  iter = 800
)
plot_constant<-plot_seromodel(
  seromodel_constant,
  serosurvey = serosurvey_constant,
  foi_df = foi_df,
  size_text = 6
)


#time varying models- slow time varying

foi_df <- data.frame(
  year = seq(2000, 2049, 1),
  foi = c(
    rep(0.2, 25),
    rep(0.1, 10),
    rep(0.00001, 15)
  )
)

survey_features <- data.frame(
  age_min = seq(1, 50, 5),
  age_max = seq(5, 50, 5),
  n_sample = rep(25, 10)
)

serosurvey_sw_dec <- simulate_serosurvey(
  "time",
  foi_df,
  survey_features
) |>
  mutate(survey_year = 2050)


foi_index <- data.frame(
  year = seq(2000, 2049),
  foi_index = rep(c(1, 2, 3), c(25, 10, 15))
)
seromodel_time_normal <- fit_seromodel(
  serosurvey = serosurvey_sw_dec,
  model_type = "time",
  foi_index = foi_index,
  iter = 1500
)
plot_time_normal<-plot_seromodel(
  seromodel_time_normal,
  serosurvey = serosurvey_sw_dec,
  foi_df = foi_df,
  size_text = 6
)

#time varying models-fast time

foi_df <- data.frame(
  year = seq(2000, 2049, 1),
  foi = c(
    rep(0, 30),
    rep(0.7, 3),
    rep(0, 17)
  )
)

survey_features <- data.frame(
  age_min = seq(1, 50, 5),
  age_max = seq(5, 50, 5),
  n_sample = rep(25, 10)
)

serosurvey_large_epi <- simulate_serosurvey(
  survey_features = survey_features,
  foi_df,
  model = "time"
) |>
  mutate(survey_year = 2050)

foi_index <- data.frame(
  year = seq(2000, 2049),
  foi_index = rep(c(1, 2, 3), c(30, 3, 17))
)
seromodel_log_time_normal <- fit_seromodel(
  serosurvey = serosurvey_large_epi,
  model_type = "time",
  is_log_foi = TRUE,
  foi_index = foi_index,
  iter = 2000
)

plot_log_time_normal <- plot_seromodel(
  seromodel_log_time_normal,
  serosurvey = serosurvey_large_epi,
  foi_df = foi_df,
  size_text = 5,
  foi_max = 0.7
)
plot(plot_log_time_normal)


#model comparison
cowplot::plot_grid(
  plot_constant, plot_time_normal, plot_log_time_normal,
  nrow = 1, ncol = 3, labels = "AUTO"
)











