
#.libPaths("M:/My Documents/R/win-library/4.2") # Adjust .libPaths if needed
library(gridExtra)
library(ggplot2)

# Adjust to your working directory
# setwd('C:/Users/linge006...')

### Read data of interest (make sure one '#' is removed)
#dataset <- "Gpulex"
#dataset <- "ChRiparius"
#dataset <- "Fcandida_Cu01"
dataset <- "Fcandida_Cu02"
#dataset <- "Fcandida_Cd01"
#dataset <- "Fcandida_Cd02"
#dataset <- "Fcandida_Cd03"
#dataset <- "simul"

# Read data set of choice along with nls model fit starting values, toxicant name and time at which elimination phase starts
if (dataset=="Gpulex") {t_scale <- "_h"; toxkin_dat <- read.csv(paste0("Gpulex_propanolol_vs",t_scale,".csv"), sep=";", header=T)
  start_nls <- list(C_0=30, k_1=0.618, k_2=0.037); toxic <- "Propanolol"; t_d <- 48}
if (dataset=="ChRiparius") {toxkin_dat <- read.csv("Ch_riparius_Ag_Silva2023.txt", sep=",", header=T)
  start_nls <- list(C_0=30, k_1=0.618, k_2=0.037); toxic <- "Ag"; t_scale <- "_h"; t_d <- 48}
if (dataset=="Fcandida_Cu01") {toxkin_dat <- read.csv("Fcandida_Cu_Ardestani2013_01.txt", sep=",", header=T)
  start_nls <- list(C_0=44.08, k_1=0.1308, k_2=0.0406); toxic <- "Cu"; t_scale <- "_d"; t_d <- 14}
if (dataset=="Fcandida_Cu02") {toxkin_dat <- read.csv("Fcandida_Cu_Ardestani2013_02.txt", sep=",", header=T)
  start_nls <- list(C_0=30, k_1=0.618, k_2=0.037); toxic <- "Cu"; t_scale <- "_d"; t_d <- 14}
if (dataset=="Fcandida_Cd01") {toxkin_dat <- read.csv("Fcandida_Cd_Ardestani2013_01.txt", sep=",", header=T)
  start_nls <- list(C_0=1, k_1=0.618, k_2=0.037); toxic <- "Cd"; t_scale <- "_d"; t_d <- 21}
if (dataset=="Fcandida_Cd02") {toxkin_dat <- read.csv("Fcandida_Cd_Ardestani2013_02.txt", sep=",", header=T)
  start_nls <- list(C_0=1, k_1=0.618, k_2=0.037); toxic <- "Cd"; t_scale <- "_d"; t_d <- 21}
if (dataset=="Fcandida_Cd03") {toxkin_dat <- read.csv("Fcandida_Cd_Ardestani2013_03.txt", sep=",", header=T)
  start_nls <- list(C_0=0.07, k_1=0.408, k_2=0.22); toxic <- "Cd"; t_scale <- "_d"; t_d <- 21}


if (dataset=="simul") {
  t_d <- 14; start_nls <- list(C_0=1, k_1=0.618, k_2=0.037); toxic <- "Toxicant"; t_scale <- "_d";
  n_row <- 28; #set.seed(95616)
  #set.seed(222222)
  toxkin_dat <- data.frame(Time=0:n_row, k_1=runif(n_row+1,0.5,1.5), k_2=runif(n_row+1,0.05,0.15), 
                         C_exp=rep(1,n_row+1), t_d=rep(t_d,n_row+1), error=runif(n_row+1,-2,2))
  #toxkin_dat <- data.frame(toxkin_dat, toxicant=with(toxkin_dat, one_compartment(1, k_1, k_2, C_exp, Time, t_d) + error)) 
  toxkin_dat <- data.frame(toxkin_dat, pred_toxicant=with(toxkin_dat, one_compartment(2, 5, 1, 1, Time, t_d) )) # C_0=1, k_1=1, k_2=0.1, C_exp=1
  toxkin_dat <- data.frame(toxkin_dat, toxicant=with(toxkin_dat, one_compartment(2, 5, 1, 1, Time, t_d) + 
                                                       rnorm(nrow(toxkin_dat), 0, pred_toxicant*0.2) ))
} # End IFs()


# Plot data
plot(toxicant~Time, data=toxkin_dat, pch=20, ylim=c(0,max(toxkin_dat$toxicant)), xlim=c(0,max(toxkin_dat$Time)))


# Assign observation to uptake or elimination for estimating two residual variance parameters
toxkin_dat <- data.frame(toxkin_dat,phase=NA) # add one column to toxkin_dat object to assign a phase to the observations
toxkin_dat$phase[toxkin_dat$Time <= t_d] <- "Uptake" # assign observations to Uptake based on sampling time
toxkin_dat$phase[toxkin_dat$Time > t_d] <- "Elimination" # assign observations to Elimination based on sampling time

C_expsr <- mean(toxkin_dat$C_exp[toxkin_dat$Time <= t_d]) # Compute mean toxicant concentration in exposure medium
toxkin_dat$C_exp[toxkin_dat$Time > t_d] <- C_expsr # Assign mean exposure concentration 

head(toxkin_dat) # Show head of data object

# Fitting of a nonlinear least squares (nls) model to the data
# the response variable toxicant is modeled using the one_compartment function with parameters C_0, k_1, k_2 and variables C_exp and Time.
# start_nls is a list that provides the initial guesses for the C_0, k_1, k_2 parameters in the model. 
# toxkin_dat is the data.frame object containing the data for fitting the model.txkin_dat contains the toxicant, Time, and t_e variables used in the one_compartment function.

# Initialize the one_compartment() and compute_BAF() functions and load the nlme package
source("ToxKin_gnls_fun.R") 

toxicant.nls <- nls(toxicant~one_compartment(C_0,k_1,k_2,C_exp,Time), 
                    start=start_nls, data=toxkin_dat)


### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ####
### gnls code for fitting frequentist (i.e. non-Bayesian) models ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### #### 

# Fit generalized nonlinear least squares (gnls) models to the data
source("ToxKin_gnls_run.R")

# Note that the message: 
# "Error in gnls(toxicant ~ one_compartment(C_0, k_1, k_2, C_exp, Time), : step halving factor reduced below minimum in NLS step"
# appears if inaccurate starting values for C_0, k_1 and k_2 are provided and may be solved by providing more accurate starting values


### ### ### ### ### ### ### ### ### ### ### ### 
### RSTAN code for fitting Bayesian models ####
### ### ### ### ### ### ### ### ### ### ### ### 

# Run and initialize the RSTAN models and additional R functions
source("ToxKin_rstan_fun.R")

### Consideration for RSTAN simulation: plot(dgamma(seq(0,10,0.01),3,1)~seq(0,10,0.01))
#   for (j in 1:N) y[j] ~ normal(m[j], sqrt(y[j])*sigma_e[no_sigma]);

### y ~ gamma(alpha, beta)
#   Increment target log probability density with gamma_lupdf(y | alpha, beta).

### real gamma_lpdf(reals y | reals alpha, reals beta)
#   The log of the gamma density of y given shape alpha and inverse scale beta


# Fit Bayesian generalized nonlinear least squares models to the data
source("ToxKin_rstan_run.R")
