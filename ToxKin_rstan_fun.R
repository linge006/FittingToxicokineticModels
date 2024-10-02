
# Note that Rtools compatible with your R version should be installed for running the rstan package. 
# It is also recommended to update R and Rtools before installing rstan
# RSTAN may be installed through running: install.packages("rstan") or install.packages("rstan", repos = c('https://stan-dev.r-universe.dev', getOption("repos")))
# Python users may like to explore PyStan, the python implementation of STAN

library(rstan)
## For execution on a local, multicore CPU with excess RAM we recommend calling
## options(mc.cores = parallel::detectCores()).
## To avoid recompilation of unchanged Stan programs, we recommend calling
## rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores()) #options(mc.cores = NUM_CORES)
rstan_options(auto_write = TRUE)
packageVersion("rstan")

# This loo package provides Bayesian model diagnostics such as the approximate leave-one-out cross-validation (LOO-CV). 
# It primarily focuses on providing methods to evaluate the predictive accuracy and performance of Bayesian models.
library(loo) 

# The bayesplot package is used for creating visualizations (such as plotting pairs using mcmc_pairs) to evaluate and understand Bayesian models. 
# It provides a flexible and comprehensive set of plotting functions specifically designed for Bayesian analysis
library(bayesplot)

print('Initializing RSTAN non-linear one-compartmental model. Note that this may take a few minutes...')

#

mod_OneComp_d <- "
// The data block contains all data that is read in from R
data {
  int<lower=0> N; // Number of observations
  real x[N]; 
  real y[N]; 
  real z; int t_d; 
  int grp[N]; int no_sigma;
  int<lower=N> P; real pred_t[P]; int s_grp[P];
  real mean_C0; real sd_C0;
  real mean_k1; real sd_k1;
  real mean_k2; real sd_k2;
  int varPow; int varExp; int varFix;
} // End of data block

// The parameters block contains all model parameters that need to be estimated
parameters {
  real<lower=0> C_0; 
  real<lower=0> k[2+varPow+varExp]; // Contains k1 and k2. It might also contain a third k representing the delta parameter in case of exponential and power residual models.
  real<lower=0> sigma_e[no_sigma]; // residual sd
} // End of parameters block

// The transformed parameters block contains the one-compartmental function for toxicokinetic modeling
transformed parameters {
  real m[N];
  
  for (i in 1:N){
    m[i] = C_0 + (x[i] <= t_d ? k[1]/k[2]*z*(1-exp(-k[2]*x[i])) : k[1]/k[2]*z*(exp(-k[2]*(x[i]-t_d))-exp(-k[2]*x[i])));
  } // End of for loop i
  
} // End of transformed parameters block

// The model block contains the priors for all parameters and the likelihood of the models 
model {
  // priors
  C_0 ~ normal(mean_C0, sd_C0); // Assign a normal prior to C_0
  k[1] ~ normal(mean_k1, sd_k1); k[2] ~ normal(mean_k2, sd_k2); // Assign normal priors to k1 and k2 
  if (varExp+varPow==1) k[3] ~ normal(0,1000); // Assign a normal prior to k3 (i.e. delta); these priors may be chosen differently
  sigma_e ~ cauchy(0,5); // Assign a half-cauch prior to residual variance parameter sigma
  
  // likelihoods
  if (no_sigma == 1 && (varFix+varExp+varPow)==0) for (j in 1:N) y[j] ~ normal(m[j], sigma_e[no_sigma]); // Homoscedastic model
  if (no_sigma > 1) for (j in 1:N) y[j] ~ normal(m[j], sigma_e[grp[j]]); // varIdent model 
  if (varFix==1) for (j in 1:N) y[j] ~ normal(m[j], sqrt(y[j])*sigma_e[no_sigma]); // varFixed model
  if (varExp==1) for (j in 1:N) y[j] ~ normal(m[j], exp(k[3]*y[j])*sigma_e[no_sigma]); // varExp model
  if (varPow==1) for (j in 1:N) y[j] ~ normal(m[j], pow(y[j],k[3])*sigma_e[no_sigma]); // varPower model
} // End of model block

// The generated quantities block contains all quantities to be returned in addition to the estimated parameters
generated quantities {
  // Initialize (additional) quantities to be generated
  vector[N] log_lik; // to use loo package and compute waic - the logLikelihood MUST be called log_lik!
  real pred[N]; // Quantities for fitted values (not including noise/error estimate)
  real res[N]; // Residuals
  real cpred[P]; real ypred[P]; // Predicted values (including noise/error estimate)
  real BAF; // Bioaccumulation factor
  
  // Assign values to the initialized quantities to be returned
  BAF = k[1]/k[2];

  for (n in 1:N){
    pred[n] = C_0 + (x[n] <= t_d ? k[1]/k[2]*z*(1-exp(-k[2]*x[n])) : k[1]/k[2]*z*(exp(-k[2]*(x[n]-t_d))-exp(-k[2]*x[n])));
    res[n] = y[n] - pred[n];
    
    if (no_sigma > 1) {log_lik[n] = normal_lpdf(y[n] | pred[n], sigma_e[grp[n]]);} // varIdent model
    if (no_sigma == 1 && (varFix+varExp+varPow)==0) {log_lik[n] = normal_lpdf(y[n] | pred[n], sigma_e[no_sigma]);} // Homoscedastic model
    if (varFix==1) {log_lik[n] = normal_lpdf(y[n] | pred[n], sqrt(y[n])*sigma_e[no_sigma]);} // varFixed model
    if (varExp==1) {log_lik[n] = normal_lpdf(y[n] | pred[n], exp(k[3]*y[n])*sigma_e[no_sigma]);} // varExp model
    if (varPow==1) {log_lik[n] = normal_lpdf(y[n] | pred[n], pow(y[n],k[3])*sigma_e[no_sigma]);} // varPower model
  
    } // End of for loop n
  
  for (p in 1:P){
    cpred[p] = C_0 + (pred_t[p] <= t_d ? k[1]/k[2]*z*(1-exp(-k[2]*pred_t[p])) : k[1]/k[2]*z*(exp(-k[2]*(pred_t[p]-t_d))-exp(-k[2]*pred_t[p])));

    if (no_sigma > 1)               {ypred[p] = normal_rng(cpred[p], sigma_e[s_grp[p]]);} // varId model
    if (no_sigma == 1 && (varFix+varExp+varPow)==0) {ypred[p] = normal_rng(cpred[p], sigma_e[no_sigma]);} // Homosc model
    if (varFix==1) {ypred[p] = normal_rng(cpred[p], sqrt(cpred[p])*sigma_e[no_sigma]);} // varFixed model
    if (varExp==1) {ypred[p] = normal_rng(cpred[p], exp(k[3]*cpred[p])*sigma_e[no_sigma]);} // varExp model
    if (varPow==1) {ypred[p] = normal_rng(cpred[p], pow(cpred[p],k[3])*sigma_e[no_sigma]);} // varPower model
    
  } // End of for loop p

} // End of generated quantities block

"

#

print("Compiling one-compartmental stan model. Note that this may take a couple of minutes...!")
stan_mod_Onecomp_d <- stan_model(model_code=mod_OneComp_d, model_name="mod_OneComp_d") # Compile stan model
print('One-compartmental RSTAN model has been compiled.')

#

# The ToxKin_data4rstan() function prepares a list with data for fitting a Bayesian model using Stan, tailored for toxicokinetic data
# The function contains a variety of arguments. These arguments are:
# 1) stan_mod: The Stan model object; 2) datafile: The dataset to be used for fitting the model.
# 3) var_y: The variable name for the response variable (y); 4) var_z: The variable name for a specific variable (z) to compute its mean.
# 5) var_t: The variable name for another predictor variable (x); 6) parms_out: Parameters to be monitored and outputted from the Stan model.
# 7) priors: A list of prior distributions for the model parameters; 8) t_e: A time-related variable, defaults to the value of t_d.
# 9) var_grp (optional): The grouping variable, default is "phase_no"; 10) burnin: The number of warmup iterations for MCMC sampling
# 11) niter: The total number of iterations for MCMC sampling; 12) n_chains: The number of MCMC chains to run
# 13) n_sigma: Number of sigma parameters; 14) varFixed: Binary variable that should take 1 for a varFixed model and 0 otherwise.
# 15 varPower: Binary variable that should take 1 for a varPower model and 0 otherwise. 16) varExponent: Binary variable that should take 1 for a varPower model and 0 otherwise.

Run_RSTAN_ToxKin_mod <- function(datafile,var_y,var_z,var_t,parms_out,priors,stan_mod=stan_mod_Onecomp_d,t_e=t_d,var_grp="phase_no",
                                burnin=3e3,niter=53e3,n_chains=4,n_sigma=1, varFixed=0,varPower=0,varExponent=0){
  cat('Number of sigmas:',n_sigma,', varFixed is:',varFixed,', varExp is:',varExponent,', varPower is:',varPower,'\n')
  t_rng <- seq(min(datafile$Time),t_e,0.1); rng_grp <- rep(1,length(t_rng)); t_rng <- c(t_rng,t_e+t_rng); rng_grp <- c(rng_grp,rng_grp+1)
  print(priors) # Print the C_0, k_1 and k_2 prior means and standard deviations
  rs.dat <- list(N=nrow(datafile), z=mean(datafile[,var_z]), t_d=t_e, y=datafile[,var_y], x=datafile[,var_t], grp=datafile[,var_grp], 
                 no_sigma=n_sigma, varFix=varFixed, varPow=varPower, varExp=varExponent,
                 P=length(t_rng), pred_t=t_rng, s_grp=rng_grp,
                 mean_C0=priors[["C_0","Value"]], sd_C0=priors["C_0","Std.Error"], 
                 mean_k1=priors[["k_1","Value"]], sd_k1=priors["k_1","Std.Error"], 
                 mean_k2=priors[["k_2","Value"]], sd_k2=priors["k_2","Std.Error"]) # Prepare data object for stan model, which should be a list
  print(rs.dat$y) # Print toxicant concentration values
  print(paste0(var_y," data from experiment initialized, RSTAN runs a simulation that started at: ",Sys.time()," ...")) # Notify user that RSTAN simulation starts
  return(sampling(stan_mod, data=rs.dat, iter=niter,
                  pars=parms_out, warmup=burnin, thin=25, chains=n_chains)) # Run and return the RSTAN Bayesian model fit object
} # End Run_RSTAN_ToxKin_mod()


# The plot_RSTAN_res() function is computes and visualizes the standardized residuals from a Bayesian model fitted using Stan
# These standardized residuals are plotted against against fitted values. 
# 1) x_rstan: The Stan fit object; 2) header: The title of the plot.
# 3) hsk: Specifies the residual variance function. Default is FALSE for Homoscedastic models.
# 4) y_lim: The y-axis limit for the plot. Default is 3; 5) x_lim): The x-axis limits for the plot. Default is c(0,100).
# 6) y_ij: Observed data points used to compute residuals and leverage; 7) y_lab: The label for the y-axis. Default is 'Residuals'.
# 8) x_lab: The label for the x-axis. Default is an empty string (no label).
# 9) phase: A variable used to define different phases for computing residual sigma when hsk is set to 'varId'.

plot_RSTAN_res <- function(x_rstan, header, hsk=FALSE, y_lim=3, x_lim=c(0,100), y_ij, y_lab='Residuals', x_lab="", phase=NULL){
  res <- summary(x_rstan,pars=c("res"))$summary[,'mean'] # Extract residuals from RSTAN object
  n_y <- length(y_ij) # Extract number of observations
  
  # Compute residual scaling factor(s)
  if (hsk=='varId'){
    sigma <- summary(x_rstan,pars=c("sigma_e[1]","sigma_e[2]"))$summary[,'mean'] # Extract sigma parameters
    sigma <- c(rep(sigma[1],table(phase)[1]), rep(sigma[2],table(phase)[2])) # Compute residual SD per observation
  } # if() statement for computing residual SD for varIdent model
  if (hsk=='varFix'){
    sigma <- summary(x_rstan,pars=c("sigma_e[1]"))$summary[,'mean'] # Extract sigma parameter
    sigma <- sqrt(y_ij)*sigma # Compute estimated model variance
  } # if() statement for computing residual SD for varFixed model
  if (hsk=='varExp'){
    sigma <- summary(x_rstan,pars=c("sigma_e[1]"))$summary[,'mean'] # Extract sigma parameter
    delta <- summary(x_rstan,pars=c("delta"))$summary[,'mean'] # Extract delta parameter
    sigma <- exp(delta*y_ij)*sigma # Compute estimated model variance
  } # if() statement for computing residual SD for varExponential model
  if (hsk=='varPow'){
    sigma <- summary(x_rstan,pars=c("sigma_e[1]"))$summary[,'mean'] # Extract sigma parameter
    delta <- summary(x_rstan,pars=c("delta"))$summary[,'mean'] # Extract delta parameter
    sigma <- y_ij^delta*sigma # Compute estimated model variance
  } # if() statement for computing residual SD for varPower model
  if (hsk==FALSE) {
    sigma <- summary(x_rstan,pars=c("sigma_e"))$summary[,'mean'] # Extract sigma parameter
  } # if() statement for computing residual SD for Homoscedastic model
  
  leverage <- 1/n_y*(y_ij-mean(y_ij))**2/(sum(y_ij-mean(y_ij)**2)) # Compute leverage
  res <- res/(sqrt((sigma**2)*(1-leverage))) # standardize residuals
  
  if (max(abs(res)) > 3) y_lim <- max(abs(res)) # Check if standardized residuals are more than 3*SD away from zero and extend y-axis range accordingly
  # Plot standardized residuals vs fitted values
  plot(res ~ summary(x_rstan,pars=c("pred"))$summary[,'mean'],
       pch=16, cex.axis=1.5, cex.lab=1.5, ylab=y_lab, xlab=x_lab, main=header, cex.main=1.75,
       ylim=y_lim*c(-1,1), xlim=c(0.8,1.1)*x_lim); abline(h=0, lty=2)
  
} # End of plot_RSTAN_res function


# The Plot_PostPredDistr() function generates 2 plots containing the posterior predictive distributions from two Bayesian models (The 95% 
# credible interval as a shaded area), the raw data points, and the model's predicted values as a curve. The first model assumes homoscedastic 
# residual variance, the second heteroscedastic. This function is useful for visually assessing the fit of the two models and comparing how 
# well they capture the variability in the data. The one_compartment function mentioned in the code is assumed to be defined. 
# The function has various arguments: 1) summ_mod01: A data.frame containing parameter and generated quantities percentiles for the homoscedastic model.
# 2) summ_mod02: A data.frame containing parameter and generated quantities percentiles for heteroscedastic model.
# 3) raw_data: The raw data used for the model fit; 4) var_x: The variable name for the x-axis, which is time here.
# 5) var_y: The variable name for the y-axis (the response variable); 6) title: The title for the combined plot.
# 7) y_lim: The y-axis limits for the plots. Default is c(-15,345).
# 8) main01: The title for the first model plot. Default is "Homoscedastic".
# 9) main02: The title for the second model plot. Default is "Heteroscedastic".

Plot_PostPredDistr <- function(summ_mod01, summ_mod02, raw_data, var_x, var_y, title, y_lim=c(-15,345), main01="Homoscedastic", main02="Heteroscedastic"){
  pars <- as.numeric(summ_mod01[1:3,1]); pars_hsk <- as.numeric(summ_mod02[1:3,1]) # Extract C_0, k_1 and k_2 for both models
  x_max <- max(raw_data[,var_x]) # Extract maximum value of raw data toxicant concentration in the organism
  C_exp <- mean(raw_data$C_exp[raw_data$Time <= 0.5*x_max]) # Compute mean of the toxicant concentration in the exposure medium
  
  pred_rows <- grepl("ypred",rownames(summ_mod01)) # Select Predicted values (including noise/error estimate) for homoscedastic model
  pred_hsk_rows <- grepl("ypred",rownames(summ_mod02)) # Select Predicted values (including noise/error estimate) from heteroscedastic model

  t_rng <- seq(min(raw_data$Time),0.5*x_max,0.1); t_rng <- c(t_rng,0.5*x_max+t_rng) # Generate time points at which error margin will be computed
  
  c_y2 <- summ_mod01[pred_rows,"X97.5."]; c_y1 <- summ_mod01[pred_rows,"X2.5."] # Extract 2.5 and 97.5 percentiles of ypred for homoscedastic model
  c_y4 <- summ_mod02[pred_hsk_rows,"X97.5."]; c_y3 <- summ_mod02[pred_hsk_rows,"X2.5."] # Extract 2.5 and 97.5 percentiles of ypred for heteroscedastic model
  
  # Model 1 plotting
  curve(one_compartment(pars[1],pars[2],pars[3],C_exp,x,0.5*x_max), from=0, to=x_max, ylim=y_lim, las=1,
        cex.lab=1.25, cex.axis=1.25, ylab="", xlab="", main=paste0(main01," model"), xaxt='n') # Plot model's predicted values as a curve
  axis(1,seq(0,x_max,0.25*x_max),cex.axis=1.25) # Plot x-axis tick marks
  polygon(c(t_rng, rev(t_rng)), c(c_y2, rev(c_y1)),
          col="lightgrey", lty=0) # Plot 95% credible interval of posterior predictive distribution as a grey shaded area
  points(raw_data[,var_x], raw_data[,var_y], pch=20) # Plot raw data points
  curve(one_compartment(pars[1],pars[2],pars[3],C_exp,x,0.5*x_max), add=T) # Plot model's predicted values as a curve again
  
  # Model 2 plotting
  curve(one_compartment(pars_hsk[1],pars_hsk[2],pars_hsk[3],C_exp,x,0.5*x_max), from=0, to=x_max, ylim=y_lim, 
        cex.lab=1.25, cex.axis=1.25, ylab="", xlab="", main=paste0(main02," model"), xaxt='n', yaxt='n') # Plot model's predicted values as a curve
  axis(1,seq(0,x_max,0.25*x_max),cex.axis=1.25) # Plot x-axis tick marks
  polygon(c(t_rng, rev(t_rng)), c(c_y4, rev(c_y3)),
          col="lightgrey", lty=0) # Plot 95% credible interval of posterior predictive distribution as a grey shaded area
  points(raw_data[,var_x], raw_data[,var_y], pch=20) # Plot raw data points
  curve(one_compartment(pars_hsk[1],pars_hsk[2],pars_hsk[3],C_exp,x,0.5*x_max), add=T) # Plot model's predicted values as a curve again
  
  # Paste outer margin axis labels
  mtext(list("Time",title), 1:2, cex=1.25, outer=T, line=c(0.25,1.0)) 
} # End plot Plot_PostPredDistr()


# The compute_loo() function computes model evaluation metrics for a Bayesian model object fitted using RSTAN.
# It calculates the Leave-One-Out Information Criterion (LOOIC) for a RSTAN model fit, it extracts the log-likelihood from the RSTAN
# object, calculates the WAIC, converts the log-likelihood to an array format and computes the relative efficiency, then
# computes the LOOIC using the extracted log-likelihood and relative efficiency and returns the LOOIC value.
# This function is useful for model selection and comparison in Bayesian analysis, where lower LOOIC values indicate a better model fit.
# The function has only one argument: a Bayesian model object fitted using RSTAN. 

compute_loo <- function(stan_out.rs){
  loglik_mod.rs <- extract_log_lik(stan_out.rs, parameter_name="log_lik", merge_chains=T) # Extract log likelihood; merge_chains=T indicates that log-likelihoods should be merged across different chains if applicable
  rs_waic <-  suppressWarnings(waic(loglik_mod.rs)) # Extract WAIC (Widely Applicable Information Criterion)
  print(rs_waic) # Print WAIC
  LL <- as.array(stan_out.rs, pars="log_lik") # Convert log Likelihood into an array
  r_eff <- relative_eff(exp(LL), cores=1) # Compute likelihood rather than log-likelihood and the relative efficiency of the MCMC chains based on these likelihoods. Calculation should be done using a single core.
  ltest <- loo(loglik_mod.rs, r_eff=r_eff, cores=1) # Compute the LOOIC using extracted log-likelihood and r_eff using a single core (cores=1).
  return(ltest) # Return LOOIC
} # End of compute_loo()

print('Functions for generating RSTAN objects, extracting residuals from RSTAN objects and computing logLik and LOOIC have been compiled.')
