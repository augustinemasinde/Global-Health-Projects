library(devtools)
## Install and load serosim 
library(serosim)
library(serosolver)
## Load additional packages required 
library(tidyverse)
library(data.table)
library(ggplot2)
library(patchwork)
library(reshape2)

## Specify the number of time steps in the simulation
times <- seq(1,120,by=1)
## Generate the population demography tibble
demography <- generate_pop_demography(N=100, times=times, prob_removal=0)
#> Joining with `by = join_by(i)`


## Create biomarker map
biomarker_map_original <- tibble(exposure_id=c("ifxn","vacc"),biomarker_id=c("IgG","IgG"))

## Reformat biomarker_map for runserosim
biomarker_map <-reformat_biomarker_map(biomarker_map_original)


## Create an empty array to store the force of exposure for all exposure types
foe_pars <- array(0, dim=c(1,max(times),n_distinct(biomarker_map$exposure_id)))
## Specify the force of exposure for exposure ID 1 which represents natural infection
foe_pars[,,1] <- 0.01
## Specify the force of exposure for exposure ID 2 which represents vaccination
foe_pars[,,2] <- 0.1


## Bring in the antibody parameters needed for the antibody model
## Note that the observation error parameter needed for the observation model (Section 1.7) is defined here too.
model_pars_path <- system.file("extdata", "model_pars_README.csv", package = "serosim")
model_pars_original <- read.csv(file = model_pars_path, header = TRUE)

## Reformat model_pars for runserosim so that exposure_id and biomarker_id are numeric and match the exposure to biomarker map
model_pars<-reformat_biomarker_map(model_pars_original)


## Run the core simulation and save outputs in "res"
res<- runserosim(
  simulation_settings=list("t_start"=1,"t_end"=max(times)),
  demography,
  observation_times=tibble(i=1:max(demography$i),t=120, b=1),
  foe_pars,
  biomarker_map,
  model_pars,
  exposure_model=exposure_model_simple_FOE,
  immunity_model=immunity_model_vacc_ifxn_simple,
  antibody_model=antibody_model_monophasic,
  observation_model=observation_model_continuous_noise,
  draw_parameters=draw_parameters_random_fx,
  
  ## Other arguments needed
  max_events=c(1,1),
  vacc_exposures=2,
  vacc_age=c(NA,9),
  sensitivity=0.85,
  specificity=0.9
)

## Plot biomarker kinetics and immune histories for 10 individuals 
plot_subset_individuals_history(res$biomarker_states, res$immune_histories_long, subset=10, demography)


## Plot the serosurvey results (observed biomarker quantities)
plot_obs_biomarkers_one_sample(res$observed_biomarker_states)

library(serosolver)
library(serosolver)
library(plyr)
library(data.table)
library(ggplot2)

#data format
titre_dat <- data.frame(individual=c(rep(1,4),rep(2,4)),
                        samples=c(8039,8040,8044,8047,8039,8041,8045,8048),
                        virus=c(rep(8036,8)),
                        titre=c(0,0,7,7,0,5,6,5),
                        run=c(rep(1,8)),
                        DOB=c(rep(8036,8)),
                        group=c(rep(1,8))
)
knitr::kable(head(titre_dat))

data(example_par_tab)

# generate starting parameter values
start_tab <- example_par_tab
for(i in 1:nrow(start_tab)){
  if(start_tab[i,"fixed"] == 0){
    start_tab[i,"values"] <- runif(1,start_tab[i,"lower_start"], 
                                   start_tab[i,"upper_start"])
  }
}


##case study 1

# Required to run serosolver
devtools::install_github("seroanalytics/serosolver")
library(serosolver)
library(plyr)
library(data.table)

## Required for this analysis
library(reshape2)
library(foreach)
library(doParallel)
library(bayesplot)
library(coda)
library(ggplot2)
library(viridis)
library(ggpubr)

# set up cluster
set.seed(0)
cl <- makeCluster(5)

## Note that this vignette was generated on a Windows machine,
## and the setup for parallelisation is different on a Linux or Mac:

if(Sys.info()[["sysname"]]=="Darwin" | Sys.info()[["sysname"]]=="Linux"){
  library(doMC)
  library(doRNG)
  registerDoMC(cores=5)
}else{
  registerDoParallel(cl)
}








