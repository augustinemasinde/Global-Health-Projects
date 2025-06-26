install.packages("pacman")
pacman::p_load(
  rio,        # importing data  
  here,       # relative file pathways  
  janitor,    # data cleaning and tables
  lubridate,  # working with dates
  matchmaker, # dictionary-based cleaning
  epikit,     # age_categories() function
  tidyverse   # data management and visualization
)
linelist_raw<- linelist_raw
skimr::skim(linelist_raw)
linelist <- linelist %>% rename(genz = 2)
