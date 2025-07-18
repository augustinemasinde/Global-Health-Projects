##import the data
library(readxl)
chikdata <- read_excel("~/Downloads/chikungunya_data_Uganda.xlsx")

# how the disease circulates
 foi <- chik_filtered %>% select(Year) %>%  mutate(foi = rep(0.02, 1405)) %>% 
   rename(year="Year")

 # Define age breaks and labels
 breaks <- c(1, 12, 18, 50, 65, 134)
 labels <- c("1-11", "12-18", "18-49", "50-64", "65+")
 
 # Create age_cat and n_sample columns
 survey_features <- chik_filtered %>%
   mutate(age_cat = cut(Age_Yrs, 
                        breaks = breaks, 
                        labels = labels, 
                        right = FALSE, 
                        include.lowest = TRUE)) %>%
   group_by(age_cat) %>%
   mutate(n_sample = n()) %>%
   ungroup() %>% select(age_cat, n_sample) %>%  distinct()
 
 
serosurvey_constant <- simulate_serosurvey(
  "time",
  foi,
  survey_features
) |>
  mutate(survey_year = 2070)
 
  