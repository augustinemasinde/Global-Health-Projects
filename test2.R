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
library(janitor)
library(tidyverse)
linelist_raw<- linelist_raw
skimr::skim(linelist_raw)
linelist <- linelist %>% rename(genz = 2)
colnames(linelist)
linelist %>% tabyl(gender, )
linelist %>% 
  count(outcome, age_unit, name = "Number", na.rm = T)

table(fct_explicit_na(linelist$gender), fct_explicit_na(linelist$outcome)) %>% 
  addmargins() %>% 
  as.data.frame.matrix() %>% 
  tibble::rownames_to_column(var = "Age Category") %>% 
  flextable::flextable()

