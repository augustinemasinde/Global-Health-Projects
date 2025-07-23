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

# set up cluster
set.seed(1234)
cl <- makeCluster(5)
registerDoParallel(cl)

## Note that this vignette was generated on a Windows machine,
## and the setup for parallelisation is different on a Linux machine
## for Linux machine:
library(doMC)
library(doRNG)
registerDoMC(cores=5)

filename <- "case_study_2"
resolution <- 1 ## eg. this would be set to 12 for monthly resolution
sample_year <- 2009

serosolver::describe_priors()
#> Which version to use in run_MCMC? The following text describes the proposal step for updating infection histories.
#> Version 1: Beta prior on per time attack rates. Explicit FOI on each epoch using probability of infection term. Proposal performs N `flip` proposals at random locations in an individual's infection history, switching 1->0 or 0->1. Otherwise, swaps the contents of two random locations
#> Version 2: Beta prior on per time attack rates. Gibbs sampling of infection histories as in Indian Buffet Process papers, integrating out each probability of infection term.
#> Version 3: Beta prior on probability of infection for an individual, assuming independence between individuals. Samples from a beta binomial with alpha and beta specified by the par_tab input. Proposes nInfs moves at a time for add/remove, or when swapping, swaps locations up to moveSize time steps away
#> Version 4: Beta prior on probability of any infection. Gibbs sampling of infection histories using total number of infections across all times and all individuals as the prior
prior_version <- 2

# Read in data
raw_dat_path <- system.file("extdata", "Fluscape_HI_data.csv", package = "serosolver")
raw_dat <- read.csv(file = raw_dat_path, stringsAsFactors = FALSE)
print(head(raw_dat))
#>   Age HI.H3N2.1968 HI.H3N2.1975 HI.H3N2.1979 HI.H3N2.1989 HI.H3N2.1995
#> 1  75           80           40           40           80          160
#> 2  35           20           80          160           40           80
#> 3  71           80           40           20           20           40
#> 4  65           80           40           40           20           40
#> 5  64          160           80           40           10           40
#> 6  33           40           20          160           80           80
#>   HI.H3N2.2002 HI.H3N2.2003 HI.H3N2.2005 HI.H3N2.2008
#> 1          160           40           80           40
#> 2           20           10            0            0
#> 3           80           20           10            0
#> 4           20            0            0            0
#> 5           40            0           20           20
#> 6          160           40           40           20

## Add indexing column for each individual
raw_dat$individual <- 1:nrow(raw_dat)

## Convert data to long format
melted_dat <- reshape2::melt(raw_dat, id.vars=c("individual","Age"),stringsAsFactors=FALSE)

## Modify column names to meet serosolver's expectations
colnames(melted_dat) <- c("individual","DOB","virus","titre")
melted_dat$virus <- as.character(melted_dat$virus)

## Extract circulation years for each virus code, which will be used 
## by serosolver as the circulation time
melted_dat$virus <- as.numeric(sapply(melted_dat$virus, function(x) strsplit(x,split = "HI.H3N2.")[[1]][2]))

## Clean and log transform the data
melted_dat <- melted_dat[complete.cases(melted_dat),]
melted_dat[melted_dat$titre == 0,"titre"] <- 5
melted_dat$titre <- log2(melted_dat$titre/5)

## Convert ages to DOB
melted_dat$DOB <- sample_year - melted_dat$DOB

## All samples taken at the same time
melted_dat$samples <- sample_year

## Add column for titre repeats, enumerating for each measurement for the same virus/sample/individual
melted_dat <- plyr::ddply(melted_dat,.(individual,virus,samples),function(x) cbind(x,"run"=1:nrow(x),"group"=1))

## Rename to data expected by serosolver
titre_dat <- melted_dat
print(head(titre_dat))
#>   individual  DOB virus titre samples run group
#> 1          1 1934  1968     4    2009   1     1
#> 2          1 1934  1975     3    2009   1     1
#> 3          1 1934  1979     3    2009   1     1
#> 4          1 1934  1989     4    2009   1     1
#> 5          1 1934  1995     5    2009   1     1
#> 6          1 1934  2002     5    2009   1     1
#> 
#> 

antigenic_coords_path <- system.file("extdata", "fonville_map_approx.csv", package = "serosolver")
antigenic_coords <- read.csv(file = antigenic_coords_path, stringsAsFactors=FALSE)
print(head(antigenic_coords))
#>   Strain    X    Y
#> 1   HK68  1.8  2.4
#> 2   EN72  2.7  4.9
#> 3   VI75  7.6  6.3
#> 4   TX77  7.9  8.8
#> 5   BK79  9.6 11.0
#> 6   SI87 15.0  6.8

## Convert to form expected by serosolver
antigenic_map <- generate_antigenic_map(antigenic_coords, resolution)
print(head(antigenic_map))
#>       x_coord   y_coord inf_times
#> 1 -0.09718111 0.5021363      1968
#> 2  0.80502804 1.6816917      1969
#> 3  1.70723718 2.8612472      1970
#> 4  2.60944633 3.9976902      1971
#> 5  3.51165548 4.8709093      1972
#> 6  4.41386463 5.5057039      1973

## More flexible version of the above function
virus_key <- c(
  "HK68" = 1968, "EN72" = 1972, "VI75" = 1975, "TX77" = 1977, 
  "BK79" = 1979, "SI87" = 1987, "BE89" = 1989, "BJ89" = 1989,
  "BE92" = 1992, "WU95" = 1995, "SY97" = 1997, "FU02" = 2002, 
  "CA04" = 2004, "WI05" = 2005, "PE06" = 2006
)
antigenic_coords$Strain <- virus_key[antigenic_coords$Strain]
antigenic_map <- generate_antigenic_map_flexible(antigenic_coords)

## Restrict entries to years of interest. Entries in antigenic_map determine
## the times that individual can be infected ie. the dimensions of the infection
## history matrix.
antigenic_map <- antigenic_map[antigenic_map$inf_times >= 1968 & antigenic_map$inf_times <= sample_year,]
strain_isolation_times <- unique(antigenic_map$inf_times)





