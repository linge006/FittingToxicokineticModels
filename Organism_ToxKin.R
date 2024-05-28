
.libPaths("M:/My Documents/R/win-library/4.2") # Adjust if needed
library(nlme)
library(gridExtra)
library(ggplot2)

setwd('...') # Set your working directory

### Read data
dataset <- "Gpulex"
#dataset <- "ChRiparius"
#dataset <- "Fcandida_Cu01"
#dataset <- "Fcandida_Cu02"
#dataset <- "Fcandida_Cd01"
#dataset <- "Fcandida_Cd02"
#dataset <- "Fcandida_Cd03"
#dataset <- "simul"

if (dataset=="Gpulex") {t_scale <- "_h"; toxkin_dat <- read.csv(paste0("Gpulex_propanolol_vs",t_scale,".csv"), sep=";", header=T)
  start_Gpulex <- list(C_0=30, k_1=0.618, k_2=0.037); toxic <- "Propanolol"; t_d <- 48}
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
  start_nls <- list(C_0=1, k_1=0.618, k_2=0.037); toxic <- "Cd"; t_scale <- "_d"; t_d <- 21}


# One compartmental model function

one_compartment <- function(C_0,k_1,k_2,C_exposure,time,t_e=t_d){
  ifelse(time <= t_e, 
         C_0+k_1/k_2*C_exposure*(1-exp(-k_2*time)), 
         C_0+k_1/k_2*C_exposure*(exp(-k_2*(time-t_e))-exp(-k_2*time)))
} # End one_compartment()

if (dataset=="simul") {
  t_d <- 14; start_nls <- list(C_0=1, k_1=0.618, k_2=0.037); toxic <- "Toxicant"; t_scale <- "_d";
  n_row <- 28; set.seed(6708)
  toxkin_dat <- data.frame(Time=0:n_row, k_1=runif(n_row+1,0.5,1.5), k_2=runif(n_row+1,0.05,0.15), 
                         C_exp=rep(1,n_row+1), t_d=rep(t_d,n_row+1), sigma=runif(n_row+1,-1,1))
  #toxkin_dat <- data.frame(toxkin_dat, toxicant=with(toxkin_dat, one_compartment(1, k_1, k_2, C_exp, Time, t_d) + sigma)) 
  toxkin_dat <- data.frame(toxkin_dat, toxicant=with(toxkin_dat, one_compartment(1, 1, 0.1, 1, Time, t_d) + sigma)) 
} # End IFs()


# Plot data
plot(toxicant~Time, data=toxkin_dat, pch=20, ylim=c(0,max(toxkin_dat$toxicant)), xlim=c(0,max(toxkin_dat$Time)))


# Assign observation to uptake or elimination for estimating two residual variance parameters
toxkin_dat <- data.frame(toxkin_dat,phase=NA)
toxkin_dat$phase[toxkin_dat$Time <= t_d] <- "Uptake"
toxkin_dat$phase[toxkin_dat$Time > t_d] <- "Elimination"

C_expsr <- mean(toxkin_dat$C_exp[toxkin_dat$Time <= t_d])
toxkin_dat$C_exp[toxkin_dat$Time > t_d] <- C_expsr

head(toxkin_dat)

# non-linear (one compartmental) gnls() model fitting 

toxicant.nls <- nls(toxicant~one_compartment(C_0,k_1,k_2,C_exp,Time), start=start_nls, 
     data=toxkin_dat)

if (Frequentist==TRUE){
  toxicant.nls <- gnls(toxicant~one_compartment(C_0,k_1,k_2,C_exp,Time), start=as.list(coef(toxicant.nls)), 
                       weights=NULL, 
                       data=toxkin_dat); summary(toxicant.nls) # only one sigma^2 is estimated
  toxicant.gnls_Id <- gnls(toxicant~one_compartment(C_0,k_1,k_2,C_exp,Time), start=as.list(coef(toxicant.nls)), 
                           weights=varIdent(form = ~ 1|phase), 
                           data=toxkin_dat); summary(toxicant.gnls_Id) # two sigma^2 are estimated
  toxicant.gnls_Fix <- gnls(toxicant~one_compartment(C_0,k_1,k_2,C_exp,Time), start=as.list(coef(toxicant.nls)), 
                            weights=varFixed(~ toxicant),
                            data=toxkin_dat); summary(toxicant.gnls_Fix) # sigma^2 scaled to y_i
  toxicant.gnls_Exp <- gnls(toxicant~one_compartment(C_0,k_1,k_2,C_exp,Time), start=as.list(coef(toxicant.nls)), 
                            weights=varExp(form=~toxicant), 
                            data=toxkin_dat); summary(toxicant.gnls_Exp)
  toxicant.gnls_Pow <- gnls(toxicant~one_compartment(C_0,k_1,k_2,C_exp,Time), start=as.list(coef(toxicant.gnls_Fix)), 
                            weights=varPower(form=~toxicant),
                            data=toxkin_dat); summary(toxicant.gnls_Pow) # Doesn't work for Cu01
  
  anova(toxicant.nls,toxicant.gnls_Fix,toxicant.gnls_Id,toxicant.gnls_Exp,toxicant.gnls_Pow)
  p_model <- anova(toxicant.nls,toxicant.gnls_Id,toxicant.gnls_Fix,toxicant.gnls_Exp,toxicant.gnls_Pow)
  write.csv(p_model,paste0(dataset,"_AIC_BIC_logLik.csv"))

  ### Compute BAF + propagated error
  
  # Extract coefficient means and standard errors
  
  compute_BAF <- function(summ_toxic.gnls){
    coef.gnls <- coef(summ_toxic.gnls)
    sigma_k1 <- coef.gnls['k_1',2]; sigma_k2 <- coef.gnls['k_2',2]
    
    # Compute k1 and k2 covariance
    sigma_k1k2 <- summ_toxic.gnls$corBeta['k_1','k_2']*sigma_k1*sigma_k2
    
    BAF <- coef.gnls['k_1',1]/coef.gnls['k_2',1] # Compute BAF
    sigma_BAF <- BAF * sqrt((sigma_k1/coef.gnls['k_1',1])^2 + 
                              (sigma_k2/coef.gnls['k_2',1])^2 -
                              2*(sigma_k1k2/prod(coef.gnls[c('k_1','k_2'),1]))^2) # Compute error propagated SE of BAF
    BAF <- data.frame(mean=BAF, se=sigma_BAF, row.names='BAF') # Write mean and se of BAF to BAF object
    return(BAF)
    #write.table(BAF, paste0(dataset,"_gnls.csv"), append=T, sep=",") # Write mean and se of BAF to file
  } # End compute_BAF()
  
  summ_toxic.gnls <- summary(list(toxicant.nls,toxicant.gnls_Id,toxicant.gnls_Fix,toxicant.gnls_Exp,toxicant.gnls_Pow)[[which.min(p_model$BIC)]])
  write.csv(summary(toxicant.nls)$tTable, paste0(dataset,"_gnls.csv")) # Write output homoscedastic model to file 
  sapply(list(compute_BAF(summary(toxicant.nls)), c(sigma=toxicant.nls$sigma)), 
         function(x) write.table(x, paste0(dataset,"_gnls.csv"), append=T, sep=",", col.names=F)) # Write BAF+sigma homoscedastic model to file 
  write.table(summ_toxic.gnls$tTable, paste0(dataset,"_gnls.csv"), append=T, sep=",") # Write output heteroscedastic model to file 
  sapply(list(compute_BAF(summ_toxic.gnls), c(sigma=summ_toxic.gnls$sigma, delta=coef(summ_toxic.gnls$modelStruct$varStruct, allCoef=TRUE))), 
         function(x) write.table(x, paste0(dataset,"_gnls.csv"), append=T, sep=",", col.names=F)) # Write BAF+sigma heteroscedastic model to file 
  
  #
  
  y_range <- c(-3.1,3.1); x_range <- range(toxkin_dat$toxicant)
  p01 <- plot(toxicant.nls, main="Homoscedastic", lty=0, ylim=y_range, xlim=x_range, pch=20)
  p02 <- plot(toxicant.gnls_Fix, main="varFixed", ylab="", ylim=y_range, xlim=x_range, pch=20)
  p03 <- plot(toxicant.gnls_Id, main="varIdent", ylab="", ylim=y_range, xlim=x_range, pch=20)
  p04 <- plot(toxicant.gnls_Exp, main="varExp", ylab="", ylim=y_range, xlim=x_range, pch=20)
  p05 <- plot(toxicant.gnls_Pow, main="varPower", ylab="", ylim=y_range, xlim=x_range, pch=20)
  
  pdf(paste0(dataset,"_ModRes_scedastic.pdf"), height=6.5, width=10)
  grid.arrange(p01, p02, p03, p04, p05,
               ncol=3, bottom="")
  dev.off()

  # ggplot fitted models
  
  x_range <- 0:max(toxkin_dat$Time)
  y_vs_x <- data.frame(var_x=x_range,
                       homosk = one_compartment(coef(toxicant.nls)[1],coef(toxicant.nls)[2],coef(toxicant.nls)[3], 
                                                mean(toxkin_dat$C_exp[toxkin_dat$Time <= t_d]),x_range),
                       var_Id = one_compartment(coef(toxicant.gnls_Id)[1],coef(toxicant.gnls_Id)[2],coef(toxicant.gnls_Id)[3], 
                                                mean(toxkin_dat$C_exp[toxkin_dat$Time <= t_d]),x_range),
                       var_Fix = one_compartment(coef(toxicant.gnls_Fix)[1],coef(toxicant.gnls_Fix)[2],coef(toxicant.gnls_Fix)[3], 
                                                 mean(toxkin_dat$C_exp[toxkin_dat$Time <= t_d]),x_range),
                       var_Exp = one_compartment(coef(toxicant.gnls_Exp)[1],coef(toxicant.gnls_Exp)[2],coef(toxicant.gnls_Exp)[3], 
                                                 mean(toxkin_dat$C_exp[toxkin_dat$Time <= t_d]),x_range),
                       var_Pow = one_compartment(coef(toxicant.gnls_Pow)[1],coef(toxicant.gnls_Pow)[2],coef(toxicant.gnls_Pow)[3], 
                                                 mean(toxkin_dat$C_exp[toxkin_dat$Time <= t_d]),x_range))
  
  p1 <- ggplot(data=y_vs_x, aes(x=x_range)) + 
    geom_line(aes(y = homosk, colour = "Homoscedastic"), linewidth=1) + 
    geom_line(aes(y = var_Id, colour = "var Id"), linewidth=1) +
    geom_line(aes(y = var_Fix, colour = "var Fixed"), linewidth=1) +
    geom_line(aes(y = var_Exp, colour = "var Exp"), linewidth=1) +
    geom_line(aes(y = var_Pow, colour = "var Pow"), linewidth=1) +
    theme(legend.position=c(0.775,0.9), legend.text=element_text(size=11), legend.title=element_text("",size=0), legend.margin=margin(0.1,0.1,0.1,0.1,'cm')) +
    ggtitle(paste0(toxic," exposure")) + ylab(expression("Content in"~paste(mu,g,"/g"))) + xlab("") #
  
  pdf(paste0(dataset,"_gnls.pdf"), height=6, width=7)
  p1 + geom_point(data=toxkin_dat, aes(y=toxicant,x=Time), color='black') 
  dev.off()
} # End if(Frequentist==TRUE)

library(rstan)
## Loading required package: StanHeaders
## rstan (Version 2.18.2, GitRev: 2e1f913d3ca3)
## For execution on a local, multicore CPU with excess RAM we recommend calling
## options(mc.cores = parallel::detectCores()).
## To avoid recompilation of unchanged Stan programs, we recommend calling
## rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores()) #options(mc.cores = NUM_CORES)
rstan_options(auto_write = TRUE)

library(loo)
library(bayesplot)

print('Initializing RSTAN non-linear one-compartmental model...')

mod_OneComp_no_d <- "
data {
  int<lower=0> N; 
  real x[N]; 
  real y[N]; 
  real z; int t_d;
  int grp[N]; int no_sigma;
  real mean_C0; real sd_C0;
  real mean_k1; real sd_k1;
  real mean_k2; real sd_k2;
  int varPow; int varExp; int varFix;
} 

parameters {
  real<lower=0> C_0; 
  real<lower=0> k_1;  
  real<lower=0> k_2;
  real<lower=0> sigma_e[no_sigma]; // residual sd
} 

transformed parameters {
  real m[N];
  
  for (i in 1:N){
    m[i] = C_0 + (x[i] < t_d ? k_1/k_2*z*(1-exp(-k_2*x[i])) : k_1/k_2*z*(exp(-k_2*(x[i]-t_d))-exp(-k_2*x[i])));
  } // End of for loop i
  
} // End tr parms

model {
  // priors
  C_0 ~ normal(mean_C0, sd_C0); 
  k_1 ~ normal(mean_k1, sd_k1); 
  k_2 ~ normal(mean_k2, sd_k2); 
  sigma_e ~ cauchy(0,5);
  
  // likelihood
  if (no_sigma == 1 && varFix==0) for (j in 1:N) y[j] ~ normal(m[j], sigma_e[no_sigma]); // varIdent
  if (no_sigma > 1) for (j in 1:N) y[j] ~ normal(m[j], sigma_e[grp[j]]); // varIdent
  if (varFix==1) for (j in 1:N) y[j] ~ normal(m[j], sqrt(y[j])*sigma_e[no_sigma]); // varFixed
}

generated quantities {
  vector[N] log_lik;
  real ypred[N]; real res[N]; real pred[N];
  real BAF;
  
  BAF = k_1/k_2;
  
  // to use loo package and compute waic - MUST be called log_lik
  for (n in 1:N){
    pred[n] = C_0 + (x[n] < t_d ? k_1/k_2*z*(1-exp(-k_2*x[n])) : k_1/k_2*z*(exp(-k_2*(x[n]-t_d))-exp(-k_2*x[n])));
    res[n] = y[n] - pred[n];
    
    if (no_sigma > 1) {ypred[n] = normal_rng(pred[n], sigma_e[grp[n]]);
      log_lik[n] = normal_lpdf(y[n] | pred[n], sigma_e[grp[n]]);} // varIdent
    if (no_sigma == 1 && varFix==0) {ypred[n] = normal_rng(pred[n], sigma_e[no_sigma]);
      log_lik[n] = normal_lpdf(y[n] | pred[n], sigma_e[no_sigma]);} // Homoscedasticity
    if (varFix==1) {ypred[n] = normal_rng(pred[n], sqrt(y[n])*sigma_e[no_sigma]); 
      log_lik[n] = normal_lpdf(y[n] | pred[n], sqrt(y[n])*sigma_e[no_sigma]);} // varFixed

    } // End of for loop

}

"
#
mod_OneComp_d <- "
data {
  int<lower=0> N; 
  real x[N]; 
  real y[N]; 
  real z; int t_d; 
  int grp[N]; int no_sigma;
  real mean_C0; real sd_C0;
  real mean_k1; real sd_k1;
  real mean_k2; real sd_k2;
  int varPow; int varExp; int varFix;
} 

parameters {
  real<lower=0> C_0; 
  real<lower=0> k_1;  
  real<lower=0> k_2;
  real<lower=0> sigma_e[no_sigma]; // residual sd
  real<lower=0> delta;
} 

transformed parameters {
  real m[N];
  
  for (i in 1:N){
    m[i] = C_0 + (x[i] < t_d ? k_1/k_2*z*(1-exp(-k_2*x[i])) : k_1/k_2*z*(exp(-k_2*(x[i]-t_d))-exp(-k_2*x[i])));
  } // End of for loop i
  
} // End tr parms

model {
  // priors
  C_0 ~ normal(mean_C0, sd_C0); 
  k_1 ~ normal(mean_k1, sd_k1); 
  k_2 ~ normal(mean_k2, sd_k2); 
  delta ~  normal(0,1000);
  sigma_e ~ cauchy(0,5);
  
  // likelihood
  if (varExp==1) for (j in 1:N) y[j] ~ normal(m[j], exp(delta*y[j])*sigma_e[no_sigma]); // varExp
  if (varPow==1) for (j in 1:N) y[j] ~ normal(m[j], pow(y[j],delta)*sigma_e[no_sigma]); // varPower
}

generated quantities {
  vector[N] log_lik;
  real ypred[N]; real res[N]; real pred[N];
  real BAF;
  
  BAF = k_1/k_2;
  
  // to use loo package and compute waic - MUST be called log_lik
  for (n in 1:N){
    pred[n] = C_0 + (x[n] < t_d ? k_1/k_2*z*(1-exp(-k_2*x[n])) : k_1/k_2*z*(exp(-k_2*(x[n]-t_d))-exp(-k_2*x[n])));
    res[n] = y[n] - pred[n];
    
    if (varExp==1) {ypred[n] = normal_rng(pred[n], exp(delta*y[n])*sigma_e[no_sigma]);
      log_lik[n] = normal_lpdf(y[n] | pred[n], exp(delta*y[n])*sigma_e[no_sigma]);} // varExp
    if (varPow==1) {ypred[n] = normal_rng(pred[n], pow(y[n],delta)*sigma_e[no_sigma]);
      log_lik[n] = normal_lpdf(y[n] | pred[n], pow(y[n],delta)*sigma_e[no_sigma]);} // varPower
  
    } // End of for loop

}

"

stan_mod_Onecomp_no_d <- stan_model(model_code=mod_OneComp_no_d, model_name="mod_OneComp_no_d")
stan_mod_Onecomp_d <- stan_model(model_code=mod_OneComp_d, model_name="mod_OneComp_d")

print('One-compartmental model has been compiled...')

#

plot_res <- function(x_rstan, titel, hsk=FALSE, y_lim=2, x_lim=c(60,80), y_ij, phase=NULL){
  res <- summary(x_rstan,pars=c("res"))$summary[,'mean'] # Extract residuals
  n_y <- length(y_ij)
  
  # Compute residual scaling factor(s)
  if (hsk=='varId'){
    sigma <- summary(x_rstan,pars=c("sigma_e[1]","sigma_e[2]"))$summary[,'mean']
    sigma <- c(rep(sigma[1],table(phase)[1]), rep(sigma[2],table(phase)[2]))
  } 
  if (hsk=='varFix'){
    sigma <- summary(x_rstan,pars=c("sigma_e[1]"))$summary[,'mean']
    sigma <- sqrt(y_ij)*sigma
  } 
  if (hsk=='varExp'){
    sigma <- summary(x_rstan,pars=c("sigma_e[1]"))$summary[,'mean']
    delta <- summary(x_rstan,pars=c("delta"))$summary[,'mean'] # Extract residuals
    sigma <- exp(delta*y_ij)*sigma
  } 
  if (hsk=='varPow'){
    sigma <- summary(x_rstan,pars=c("sigma_e[1]"))$summary[,'mean']
    delta <- summary(x_rstan,pars=c("delta"))$summary[,'mean'] # Extract residuals
    sigma <- y_ij^delta*sigma
  } 
  if (hsk==FALSE) {
    sigma <- summary(x_rstan,pars=c("sigma_e"))$summary[,'mean']
  } # End if()s
  
  leverage <- 1/n_y*(y_ij-mean(y_ij))**2/(sum(y_ij-mean(y_ij)**2))
  res <- res/(sqrt((sigma**2)*(1-leverage))) # standardize residuals
  
  if (max(abs(res)) > 3) y_lim <- max(abs(res))
  plot(res ~ summary(x_rstan,pars=c("pred"))$summary[,'mean'],
       pch=16, cex.axis=1.25, cex.lab=1.25, ylab='Residuals', xlab='Fitted', main=titel, 
       ylim=y_lim*c(-1,1), xlim=1.1*x_lim); abline(h=0)
  
} # End of plot_res function


### RSTAN sampling/simulation function ###

toxkin_dat <- data.frame(toxkin_dat, phase_no=NA)
toxkin_dat$phase_no[toxkin_dat$phase == "Uptake"] <- 1
toxkin_dat$phase_no[toxkin_dat$phase == "Elimination"] <- 2

organism_data4rstan <- function(stan_mod,datafile,var_y,var_z,var_t,parms_out,priors,t_e=t_d,
                              burnin=3e3,niter=53e3,n_chains=4,n_sigma=1, varFixed=0,varPower=0,varExponent=0){
  cat('Number of sigmas:',n_sigma,', varFixed is:',varFixed,', varExp is:',varExponent,', varPower is:',varPower,'\n')
  print(priors)
  rs.dat <- list(N=nrow(datafile), z=mean(datafile[,var_z]), t_d=t_e, y=datafile[,var_y], x=datafile[,var_t], grp=datafile$phase_no, no_sigma=n_sigma,
                 varFix=varFixed, varPow=varPower, varExp=varExponent,
                 mean_C0=priors[["C_0","Value"]], sd_C0=priors["C_0","Std.Error"], 
                 mean_k1=priors[["k_1","Value"]], sd_k1=priors["k_1","Std.Error"], 
                 mean_k2=priors[["k_2","Value"]], sd_k2=priors["k_2","Std.Error"])
  print(rs.dat$y)
  print(paste0(var_y," data from experiment initialized, RSTAN runs..."))
  return(sampling(stan_mod, data=rs.dat, iter=niter,
                  pars=parms_out, warmup=burnin, thin=25, chains=n_chains))
} # End organism_data4rstan()

#

pars_no_d <- c("C_0","k_1","k_2","BAF","sigma_e","pred","res","ypred","log_lik")
pars_d <- c("C_0","k_1","k_2","BAF","delta","sigma_e","pred","res","ypred","log_lik")

#nls_prior <- summary(toxicant.gnls_Exp)$tTable
nls_prior <- summary(toxicant.nls)$parameters; colnames(nls_prior)[1:2] <- c("Value","Std.Error")
#if (dataset=="Fcandida_Cu01") {nls_prior <- summary(toxicant.nls)$parameters; colnames(nls_prior)[1:2] <- c("Value","Std.Error")}

stan_toxic_out.rs <- organism_data4rstan(stan_mod_Onecomp_no_d, toxkin_dat,"toxicant","C_exp","Time",pars_no_d,nls_prior)
stan_toxic_Id_out.rs <- organism_data4rstan(stan_mod_Onecomp_no_d, toxkin_dat,"toxicant","C_exp","Time",pars_no_d,nls_prior,n_sigma=2)
stan_toxic_Fix_out.rs <- organism_data4rstan(stan_mod_Onecomp_no_d, toxkin_dat,"toxicant","C_exp","Time",pars_no_d,nls_prior, varFixed=1) # Fit S
stan_toxic_Exp_out.rs <- organism_data4rstan(stan_mod_Onecomp_d, toxkin_dat,"toxicant","C_exp","Time",pars_d,nls_prior, varExponent=1) # Fit STAN model
stan_toxic_Pow_out.rs <- organism_data4rstan(stan_mod_Onecomp_d, toxkin_dat,"toxicant","C_exp","Time",pars_d,nls_prior, varPower=1) # Fit S

summ_toxic_rs <- summary(stan_toxic_out.rs,pars=pars_no_d[1:5])$summary # Extract model parameter summary
summ_toxic_Id_rs <- summary(stan_toxic_Id_out.rs,pars=pars_no_d[1:5])$summary # Extract model parameter summary
summ_toxic_Fix_rs <- summary(stan_toxic_Fix_out.rs,pars=pars_no_d[1:5])$summary # Extract model parameter summary
summ_toxic_Exp_rs <- summary(stan_toxic_Exp_out.rs,pars=pars_d[1:6])$summary # Extract model parameter summary
summ_toxic_Pow_rs <- summary(stan_toxic_Pow_out.rs,pars=pars_d[1:6])$summary # Extract model parameter summary

print(summ_toxic_rs); print(summ_toxic_Id_rs); print(summ_toxic_Fix_rs); print(summ_toxic_Exp_rs); print(summ_toxic_Pow_rs)

# Write model parameter summaries
write.csv(summary(stan_toxic_out.rs)$summary,paste0(dataset,'_toxkin_STAN.csv')) # For toxicant 
write.csv(summary(stan_toxic_Id_out.rs)$summary,paste0(dataset,'_toxkin_Id_STAN.csv')) # For toxicant 
write.csv(summary(stan_toxic_Fix_out.rs)$summary,paste0(dataset,'_toxkin_Fix_STAN.csv')) # For toxicant 
write.csv(summary(stan_toxic_Exp_out.rs)$summary,paste0(dataset,'_toxkin_Exp_STAN.csv')) # For toxicant 
write.csv(summary(stan_toxic_Pow_out.rs)$summary,paste0(dataset,'_toxkin_Pow_STAN.csv')) # For toxicant 


# Write chains to files
dev.off()
pdf(paste0('Chains_',toxic,'_',dataset,'.pdf'), height=7, width=10) # ,units='in',res=600,compression='lzw'
traceplot(stan_toxic_out.rs, pars=c("C_0","k_1","k_2","BAF","sigma_e"))
traceplot(stan_toxic_Id_out.rs, pars=c("C_0","k_1","k_2","BAF","sigma_e"))
traceplot(stan_toxic_Fix_out.rs, pars=c("C_0","k_1","k_2","BAF","sigma_e"))
traceplot(stan_toxic_Exp_out.rs, pars=c("C_0","k_1","k_2","BAF","delta","sigma_e"))
traceplot(stan_toxic_Pow_out.rs, pars=c("C_0","k_1","k_2","BAF","delta","sigma_e"))
dev.off()


# Write residual figures to file
x_range <- c(0,max(toxkin_dat$toxicant))
pdf(paste0('Res_STAN_',dataset,'.pdf'), height=7, width=8)
par(mfrow=c(2,3))
plot_res(stan_toxic_out.rs, hsk=FALSE, "Homoscedastic", y_lim=3, x_lim=x_range, y_ij=toxkin_dat$toxicant)
plot_res(stan_toxic_Id_out.rs, hsk="varId", "varId", y_lim=3, x_lim=x_range, y_ij=toxkin_dat$toxicant, toxkin_dat$phase)
plot_res(stan_toxic_Fix_out.rs, hsk="varFix", "varFix", y_lim=3, x_lim=x_range, y_ij=toxkin_dat$toxicant)
plot_res(stan_toxic_Exp_out.rs, hsk="varExp", "varExp", y_lim=3, x_lim=x_range, y_ij=toxkin_dat$toxicant)
plot_res(stan_toxic_Pow_out.rs, hsk="varPow", "varPow", y_lim=3, x_lim=x_range, y_ij=toxkin_dat$toxicant)
dev.off()

###

# Generate fitted model and data Figure
x_range <- 0:max(toxkin_dat$Time)
y_vs_x <- data.frame(var_x=x_range,
                     Homosk = one_compartment(summ_toxic_rs['C_0','mean'],summ_toxic_rs['k_1','mean'],summ_toxic_rs['k_2','mean'],C_expsr,x_range),
                     var_Id = one_compartment(summ_toxic_Id_rs['C_0','mean'],summ_toxic_Id_rs['k_1','mean'],summ_toxic_Id_rs['k_2','mean'],C_expsr,x_range),
                     var_Fix = one_compartment(summ_toxic_Fix_rs['C_0','mean'],summ_toxic_Fix_rs['k_1','mean'],summ_toxic_Fix_rs['k_2','mean'],C_expsr,x_range),
                     var_Exp = one_compartment(summ_toxic_Exp_rs['C_0','mean'],summ_toxic_Exp_rs['k_1','mean'],summ_toxic_Exp_rs['k_2','mean'],C_expsr,x_range),
                     var_Pow = one_compartment(summ_toxic_Pow_rs['C_0','mean'],summ_toxic_Pow_rs['k_1','mean'],summ_toxic_Pow_rs['k_2','mean'],C_expsr,x_range))

p4 <- ggplot(data=toxkin_dat, aes(y=toxicant,x=Time)) + geom_point(color='black') +
  geom_line(data=y_vs_x, aes(var_x, Homosk, color = "black"), lwd=1) + theme(legend.position="none") + ggtitle("Homoscedastic variance") + 
  ylab(expression("Toxicant content in"~paste("ng/g"))) + xlab("")

p5 <- ggplot(data=toxkin_dat, aes(y=toxicant,x=Time)) + geom_point(color='black') +
  geom_line(data=y_vs_x, aes(var_x, var_Exp, color = "black"), lwd=1) + theme(legend.position="none") + ggtitle("Exponential variance") + 
  ylab(expression("Toxicant content in"~paste("ng/g"))) + xlab("")

p6 <- ggplot(data=toxkin_dat, aes(y=toxicant,x=Time)) + geom_point(color='black') +
  geom_line(data=y_vs_x, aes(var_x, var_Pow, color = "black"), lwd=1) + theme(legend.position="none") + ggtitle("Power variance") + 
  ylab(expression("Toxicant content in"~paste("ng/g"))) + xlab("")

grid.arrange(p4, p5, p6,
             ncol=3, bottom="")
dev.off()

p1 <- ggplot(data=y_vs_x, aes(x=x_range)) + 
  geom_line(aes(y = Homosk, colour = "Homoscedastic"), linewidth=1) + 
  geom_line(aes(y = var_Id, colour = "var Id"), linewidth=1) +
  geom_line(aes(y = var_Fix, colour = "var Fixed"), linewidth=1) +
  geom_line(aes(y = var_Exp, colour = "var Exp"), linewidth=1) +
  geom_line(aes(y = var_Pow, colour = "var Pow"), linewidth=1) +
  theme(legend.position=c(0.775,0.9), legend.text=element_text(size=11), legend.title=element_text("",size=0), legend.margin=margin(0.1,0.1,0.1,0.1,'cm')) +
  ggtitle(paste0(toxic," exposure")) + ylab(expression("Content in"~paste(mu,g,"/g"))) + xlab("") #

pdf(paste0(toxic,'_',dataset,t_scale,'_RSTAN.pdf'), height=6, width=7)
p1 + geom_point(data=toxkin_dat, aes(y=toxicant,x=Time), color='black') 
dev.off()

### ### ###

### LOOIC ###

compute_loo <- function(stan_out.rs){
  loglik_mod.rs <- extract_log_lik(stan_out.rs, parameter_name="log_lik", merge_chains=T) # Extract loglikelihood
  rs_waic <-  suppressWarnings(waic(loglik_mod.rs)) # Extract WAIC
  print(rs_waic)
  LL <- as.array(stan_out.rs, pars="log_lik")
  r_eff <- relative_eff(exp(LL), cores=1)
  ltest <- loo(loglik_mod.rs, r_eff=r_eff, cores=1) # Extract log-likelihood, r_eff and looic
  return(ltest)
} # End of compute_loo()


loo_homosk <- compute_loo(stan_toxic_out.rs)
loo_varId <- compute_loo(stan_toxic_Id_out.rs)
loo_varFix <- compute_loo(stan_toxic_Fix_out.rs)
loo_varExp <- compute_loo(stan_toxic_Exp_out.rs)
loo_varPow <- compute_loo(stan_toxic_Pow_out.rs)

# Compare models using logLik
loo_compare(loo_homosk,loo_varId,loo_varFix,loo_varExp,loo_varPow)

file_name <- paste0("loo_",dataset,"_RSTAN.csv")
write.csv(loo_homosk$estimates, file_name)
write.table("varId", file_name, col.names=F, row.names=F, append=T)
write.table(loo_varId$estimates, file_name, col.names=F, append=T, sep=",")
write.table("varFixed", file_name, col.names=F, row.names=F, append=T)
write.table(loo_varFix$estimates, file_name, col.names=F, append=T, sep=",")
write.table("varExponential", file_name, col.names=F, row.names=F, append=T)
write.table(loo_varExp$estimates, file_name, col.names=F, append=T, sep=",")
write.table("varPower", file_name, col.names=F, row.names=F, append=T)
write.table(loo_varPow$estimates, file_name, col.names=F, append=T, sep=",")

##### #####
