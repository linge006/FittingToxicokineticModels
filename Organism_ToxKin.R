
.libPaths("M:/My Documents/R/win-library/4.2") # Adjust .libPaths if needed
library(gridExtra)
library(ggplot2)

# Adjust to your working directory
setwd('C:/Users/linge006/OneDrive - Wageningen University & Research/Documenten/WUR_HOME/Diagonal/Del 4.8/Experiment ZnOMn & Mn/Chemical analysis')

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


#dataset <- 'Algae_Zn_4n4'
#dataset <- 'Algae_Zn_4n6'
#dataset <- 'Algae_Mn_4n6'
#if (dataset=="Algae_Zn_4n4") {toxkin_dat <- read.csv("toxicodynamics4n4v2.txt", sep="\t", header=T)
#start_nls <- list(C_0=1, k_1=0.618, k_2=0.037); toxic <- "Zn"; t_scale <- "_d"; t_d <- 48}
#if (dataset=="Algae_Zn_4n6") {toxkin_dat <- read.csv("toxicodynamics4n6vZn.txt", sep="\t", header=T)
#start_nls <- list(C_0=1, k_1=0.618, k_2=0.037); toxic <- "Zn"; t_scale <- "_d"; t_d <- 48}
#if (dataset=="Algae_Mn_4n6") {toxkin_dat <- read.csv("toxicodynamics4n6vZnMn.txt", sep="\t", header=T)
#start_nls <- list(C_0=1, k_1=0.618, k_2=0.037); toxic <- "Mn"; t_scale <- "_d"; t_d <- 48}

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

# Initialize the one_compartment() and compute_BAF() functions
source("ToxKin_gnls_fun.R")

# Fitting of a nonlinear least squares (nls) model to the data
# the response variable toxicant is modeled using the one_compartment function with parameters C_0, k_1, k_2 and variables C_exp and Time.
# start_nls is a list that provides the initial guesses for the C_0, k_1, k_2 parameters in the model. 
# toxkin_dat is the data.frame object containing the data for fitting the model.txkin_dat contains the toxicant, Time, and t_e variables used in the one_compartment function.

toxicant.nls <- nls(toxicant~one_compartment(C_0,k_1,k_2,C_exp,Time), 
                    start=start_nls, data=toxkin_dat)

# Fitting of a generalized nonlinear least squares (gnls) models to the data
# variables and parameters are as for the nls model fit
# In addition, the weights argument determines the residual variance function. NULL indicates homoscedastic residual variance, 
  # whereas the varIdent(), varFixed(), varExp() and varPower() functions are the considered heteroscedastic variance assumptions 
# The estimated parameters returned from the nls() model are used as initial values for the gnls() model fits


# If frequentist==T, then residual variance models will be fitted using frequentist (non-Bayesian) statistics
# Otherwise, the user may proceed with fitting models in a Bayesian framework
Frequentist <- FALSE 

if (Frequentist==TRUE){
  toxicant.gnls <- gnls(toxicant~one_compartment(C_0,k_1,k_2,C_exp,Time), start=as.list(coef(toxicant.nls)), 
                       weights=NULL, 
                       data=toxkin_dat); summary(toxicant.gnls) # only one sigma^2 is estimated
  toxicant.gnls_Id <- gnls(toxicant~one_compartment(C_0,k_1,k_2,C_exp,Time), start=as.list(coef(toxicant.nls)), 
                           weights=varIdent(form = ~ 1|phase), 
                           data=toxkin_dat); summary(toxicant.gnls_Id) # two sigma^2 are estimated
  toxicant.gnls_Fix <- gnls(toxicant~one_compartment(C_0,k_1,k_2,C_exp,Time), start=as.list(coef(toxicant.nls)), 
                            weights=varFixed(~ toxicant),
                            data=toxkin_dat); summary(toxicant.gnls_Fix) # sigma^2 scaled to y_i
  toxicant.gnls_Exp <- gnls(toxicant~one_compartment(C_0,k_1,k_2,C_exp,Time), start=as.list(coef(toxicant.nls)), 
                            weights=varExp(form=~toxicant), 
                            data=toxkin_dat); summary(toxicant.gnls_Exp)
  toxicant.gnls_Pow <- gnls(toxicant~one_compartment(C_0,k_1,k_2,C_exp,Time), start=as.list(coef(toxicant.nls)), 
                            weights=varPower(form=~toxicant),
                            data=toxkin_dat); summary(toxicant.gnls_Pow) # Doesn't work for Cu01

    
  # The anova() function performs an analysis of variance (ANOVA) to compare the model fits of the considered gnls() models to the data
  # The anova() function returns a table showing statistics such as degrees of freedom, log likelihood, goodness of fit relative to 
    # model complexity (AIC and BIC), and the likelihood ratio and the associated P-value of nested models
  # These statistics help to determine if more complex models (with more parameters) provide a significantly better fit than simpler models.

  anova(toxicant.gnls,toxicant.gnls_Fix,toxicant.gnls_Id,toxicant.gnls_Exp,toxicant.gnls_Pow)
  p_model <- anova(toxicant.gnls,toxicant.gnls_Id,toxicant.gnls_Fix,toxicant.gnls_Exp,toxicant.gnls_Pow)
  write.csv(p_model,paste0(dataset,"_AIC_BIC_logLik.csv")) # Write output of anova() function for the considered models to file

  
  # Select the model with the lowest BIC (best goodness of relative to number of model parameters) and take summary of that model
  summ_toxic.gnls <- summary(list(toxicant.gnls,toxicant.gnls_Id,toxicant.gnls_Fix,toxicant.gnls_Exp,toxicant.gnls_Pow)[[which.min(p_model$BIC)]])
  write.csv(summary(toxicant.gnls)$tTable, paste0(dataset,"_gnls.csv")) # Write output homoscedastic model to file 
  sapply(list(compute_BAF(summary(toxicant.gnls)), c(sigma=toxicant.gnls$sigma)), 
         function(x) write.table(x, paste0(dataset,"_gnls.csv"), append=T, sep=",", col.names=F)) # Write BAF+sigma homoscedastic model to file 
  write.table(summ_toxic.gnls$tTable, paste0(dataset,"_gnls.csv"), append=T, sep=",") # Write output heteroscedastic model with lowest BIC to file 
  sapply(list(compute_BAF(summ_toxic.gnls), c(sigma=summ_toxic.gnls$sigma, delta=coef(summ_toxic.gnls$modelStruct$varStruct, allCoef=TRUE))), 
         function(x) write.table(x, paste0(dataset,"_gnls.csv"), append=T, sep=",", col.names=F)) # Write BAF+sigma heteroscedastic model with lowest BIC to file 
  
  
  ### Generate residuals plots (standardized residuals vs fitted values) for all 5 gnls() models
  # A range is set for the y-axis (residuals) and the x-axis (fitted values) using the values of the toxicant
  y_range <- c(-3.1,3.1); x_range <- range(toxkin_dat$toxicant)
  
  # The plot() function is used for generating residual plots, with the first argument being the gnls model fit. 
    # main determines the title of the plot, ylim and xlim the ranges of the y- and x-axes, respectively
    # pch determines the symbol for the data points, ylab="" and xlab="" omit axis labels
  p01 <- plot(toxicant.gnls, main="Homoscedastic", xlab="", ylim=y_range, xlim=x_range, pch=20)
  p02 <- plot(toxicant.gnls_Fix, main="varFixed", xlab="", ylab="", ylim=y_range, xlim=x_range, pch=20)
  p03 <- plot(toxicant.gnls_Id, main="varIdent", xlab="", ylab="", ylim=y_range, xlim=x_range, pch=20)
  p04 <- plot(toxicant.gnls_Exp, main="varExp", xlab="", ylim=y_range, xlim=x_range, pch=20)
  p05 <- plot(toxicant.gnls_Pow, main="varPower", ylab="", ylim=y_range, xlim=x_range, pch=20)
  
  # The residuals are saved as a .pdf file. The pdf() function takes the file name and then the height and width of the pdf in inches
  # The grid.arrange() function takes the five residual plots as figure panels and orders due to the numbers of columns given to ncol
  # bottom="" omits any outer label at the x-axis
  # Running dev.off() closes the pdf file and ensures that no more plots are saved in the pdf file
  pdf(paste0(dataset,"_ModRes_scedastic.pdf"), height=6.5, width=10)
  grid.arrange(p01, p02, p03, p04, p05,
               ncol=3, bottom="")
  dev.off()

  
  ### Show fitted models using ggplot
  
  # seq() generates a sequence of numbers that represent the time values in the data set
  # The sequence starts at the minimum value, with the floor() function rounding it down to the nearest integer
  # The sequence ends at the maximum value, with the ceiling() function rounding it up to the nearest integer
  t_range <- seq(floor(min(toxkin_dat$Time)), ceiling(max(toxkin_dat$Time)), 1)
  
  # Create a data.frame named y_vs_t for the 5 different models, toxicant values are generated using the one_compartment() model function and 
    # the estimated gnls() model coefficients are inputs along with the generated sequence of numbers and C_expsr (medium exposure concentration)
  y_vs_t <- data.frame(var_t=t_range,
                       homosk = one_compartment(coef(toxicant.gnls)[1],coef(toxicant.gnls)[2],coef(toxicant.gnls)[3], 
                                                C_expsr,t_range),
                       var_Id = one_compartment(coef(toxicant.gnls_Id)[1],coef(toxicant.gnls_Id)[2],coef(toxicant.gnls_Id)[3], 
                                                C_expsr,t_range),
                       var_Fix = one_compartment(coef(toxicant.gnls_Fix)[1],coef(toxicant.gnls_Fix)[2],coef(toxicant.gnls_Fix)[3], 
                                                 C_expsr,t_range),
                       var_Exp = one_compartment(coef(toxicant.gnls_Exp)[1],coef(toxicant.gnls_Exp)[2],coef(toxicant.gnls_Exp)[3], 
                                                 C_expsr,t_range),
                       var_Pow = one_compartment(coef(toxicant.gnls_Pow)[1],coef(toxicant.gnls_Pow)[2],coef(toxicant.gnls_Pow)[3], 
                                                 C_expsr,t_range))
  
  # Generate a plot with multiple lines, each representing the predictions from the 5 different models assigned to object p1.
  # The ggplot() function uses the y_vs_t data.frame and sets the aesthetic mapping for the x-axis to t_range using aes(x = t_range).
  # Then, to plot a line for each model, predicted values are taken from columns of the y_vs_t data.frame, the line color 
    # is set to the variance functions, and the linewidth argument sets the line width
  # theme(.) customizes the theme using: legend.position = c(0.775, 0.9) to places the legend at the specified coordinates within the plot area, 
    # legend.text = element_text(size=11) sets the font size of the legend text to 11, legend.title = element_text("",size=0) sets the 
    # legend title to an empty string with a size of 0 (effectively removing it), legend.margin=margin(0.1,0.1,0.1,0.1,'cm') sets the margin around the legend.
  # ggtitle(paste0(toxic, " exposure")) sets the title of the plot to the string concatenated from the toxic variable and " exposure", 
    # ylab(expression("Content in" ~ paste(mu, g, "/g"))) sets the y-axis label to "Content in µg/g" using the expression function to properly format the Greek letter mu (and the unit "g/g".
  p1 <- ggplot(data=y_vs_t, aes(x=t_range)) + 
    geom_line(aes(y = homosk, colour = "Homoscedastic"), linewidth=1) + 
    geom_line(aes(y = var_Id, colour = "var Id"), linewidth=1) +
    geom_line(aes(y = var_Fix, colour = "var Fixed"), linewidth=1) +
    geom_line(aes(y = var_Exp, colour = "var Exp"), linewidth=1) +
    geom_line(aes(y = var_Pow, colour = "var Pow"), linewidth=1) +
    theme(legend.position=c(0.775,0.9), legend.text=element_text(size=11), legend.title=element_text("",size=0), legend.margin=margin(0.1,0.1,0.1,0.1,'cm')) +
    ggtitle(paste0(toxic," exposure")) + ylab(expression("Content in"~paste(mu,g,"/g"))) + xlab("") #
  
  # Generate a pdf file named according to the dataset variable with the suffix "_gnls.pdf". The pdf contains the plot p1 with additional black 
    # points representing data from the toxkin_dat data.frame, where the toxicant variable is plotted against the Time variable. 
    # The height and widht of the pdf file are set to 6 and 7 inches, respectively. The dev.off() function closes and saves the pdf file.
  pdf(paste0(dataset,"_gnls.pdf"), height=6, width=7)
  p1 + geom_point(data=toxkin_dat, aes(y=toxicant,x=Time), color='black') 
  dev.off()
} # End if(Frequentist==TRUE)


### ### ### ### ### ### ### ### ### ### ### ### 
### RSTAN code for fitting Bayesian models ####
### ### ### ### ### ### ### ### ### ### ### ### 

# Run and initialize the RSTAN models and additional R functions
source("ToxKin_rstan_fun.R")

toxkin_dat <- data.frame(toxkin_dat, phase_no=NA) # Create extra column in data.frame for toxicokinetic phases. These should be integers in STAN
toxkin_dat$phase_no[toxkin_dat$phase == "Uptake"] <- 1 # Assign uptake phase to 1
toxkin_dat$phase_no[toxkin_dat$phase == "Elimination"] <- 2 # Assign elimination phase to 2

#

pars_no_d <- c("C_0","k_1","k_2","BAF","sigma_e","pred","res","ypred","log_lik") # Initialize parameters and generated quantities to be returned for Homoscedastic, varId and varFixed models
pars_d <- c("C_0","k_1","k_2","BAF","delta","sigma_e","pred","res","ypred","log_lik") # Initialize parameters and generated quantities to be returned for varExp and varPower models

# Take means and standard errors from nls model as priors for mean and standard deviation for C_0, k_1 and k_2 parameters of RSTAN model fit 
nls_prior <- abs(summary(toxicant.nls)$parameters[,1:2]); colnames(nls_prior) <- c("Value","Std.Error") 

# Fit Homoscedastic, varId, varFix, varExp and varPow models using RSTAN
stan_toxic_out.rs <- Run_RSTAN_ToxKin_mod(stan_mod_Onecomp_no_d, toxkin_dat,"toxicant","C_exp","Time",pars_no_d,nls_prior) 
stan_toxic_Id_out.rs <- Run_RSTAN_ToxKin_mod(stan_mod_Onecomp_no_d, toxkin_dat,"toxicant","C_exp","Time",pars_no_d,nls_prior, n_sigma=2)
stan_toxic_Fix_out.rs <- Run_RSTAN_ToxKin_mod(stan_mod_Onecomp_no_d, toxkin_dat,"toxicant","C_exp","Time",pars_no_d,nls_prior, varFixed=1)
stan_toxic_Exp_out.rs <- Run_RSTAN_ToxKin_mod(stan_mod_Onecomp_d, toxkin_dat,"toxicant","C_exp","Time",pars_d,nls_prior, varExponent=1)
stan_toxic_Pow_out.rs <- Run_RSTAN_ToxKin_mod(stan_mod_Onecomp_d, toxkin_dat,"toxicant","C_exp","Time",pars_d,nls_prior, varPower=1)

# Extract model parameter summaries
summ_toxic_rs <- summary(stan_toxic_out.rs,pars=pars_no_d[1:5])$summary 
summ_toxic_Id_rs <- summary(stan_toxic_Id_out.rs,pars=pars_no_d[1:5])$summary 
summ_toxic_Fix_rs <- summary(stan_toxic_Fix_out.rs,pars=pars_no_d[1:5])$summary 
summ_toxic_Exp_rs <- summary(stan_toxic_Exp_out.rs,pars=pars_d[1:6])$summary 
summ_toxic_Pow_rs <- summary(stan_toxic_Pow_out.rs,pars=pars_d[1:6])$summary 

# Print RSTAN fit model posteriors from summaries
print(summ_toxic_rs); print(summ_toxic_Id_rs); print(summ_toxic_Fix_rs); print(summ_toxic_Exp_rs); print(summ_toxic_Pow_rs)

# Write model parameter summaries to files
write.csv(summary(stan_toxic_out.rs)$summary,paste0(dataset,'_toxkin_STAN.csv')) 
write.csv(summary(stan_toxic_Id_out.rs)$summary,paste0(dataset,'_toxkin_Id_STAN.csv'))
write.csv(summary(stan_toxic_Fix_out.rs)$summary,paste0(dataset,'_toxkin_Fix_STAN.csv'))
write.csv(summary(stan_toxic_Exp_out.rs)$summary,paste0(dataset,'_toxkin_Exp_STAN.csv'))
write.csv(summary(stan_toxic_Pow_out.rs)$summary,paste0(dataset,'_toxkin_Pow_STAN.csv'))


# Write residual plots to file
t_range <- range(predict(toxicant.nls)); y_range <- 3 # Set ranges for x-axis (fitted value) based on nls model predicted values and y-axis (residuals) 
if (grepl("Fcandida",dataset)) t_range[1] <- 0  # Adjust x-axis lower limit to 0 for Fcandida datasets

# Save Figure with 2x3 panels to file
pdf(paste0('Res_STAN_',dataset,'.pdf'), height=6.5, width=10)
par(mfrow=c(2,3))
plot_RSTAN_res(stan_toxic_out.rs, hsk=FALSE, "Homoscedastic", y_lim=y_range, x_lim=t_range, y_ij=toxkin_dat$toxicant)
plot_RSTAN_res(stan_toxic_Id_out.rs, hsk="varId", "varId", y_lim=y_range, x_lim=t_range, y_ij=toxkin_dat$toxicant, y_lab="", phase=toxkin_dat$phase)
plot_RSTAN_res(stan_toxic_Fix_out.rs, hsk="varFix", "varFix", y_lim=y_range, x_lim=t_range, y_ij=toxkin_dat$toxicant)
plot_RSTAN_res(stan_toxic_Exp_out.rs, hsk="varExp", "varExp", y_lim=y_range, x_lim=t_range, y_lab="", x_lab="Fitted values", y_ij=toxkin_dat$toxicant)
plot_RSTAN_res(stan_toxic_Pow_out.rs, hsk="varPow", "varPow", y_lim=y_range, x_lim=t_range, y_lab="", y_ij=toxkin_dat$toxicant)
dev.off()


# Save parameter correlations to png files. Note that using pdf files instead may results in very large file sizes
png(paste0(dataset,'_pairs_Homosc.png'), height=8.5, width=13, res=300, units="in")
mcmc_pairs(as.array(stan_toxic_out.rs), pars = c("C_0","k_1","k_2","BAF","sigma_e[1]"),
           off_diag_args = list(size = 0.75)); dev.off()
png(paste0(dataset,'_pairs_varIdc.png'), height=8.5, width=13, res=300, units="in")
mcmc_pairs(as.array(stan_toxic_Id_out.rs), pars = c("C_0","k_1","k_2","BAF","sigma_e[1]","sigma_e[2]"),
           off_diag_args = list(size = 0.75)); dev.off()
png(paste0(dataset,'_pairs_varFix.png'), height=8.5, width=13, res=300, units="in")
mcmc_pairs(as.array(stan_toxic_Fix_out.rs), pars = c("C_0","k_1","k_2","BAF","sigma_e[1]"),
           off_diag_args = list(size = 0.75)); dev.off()
png(paste0(dataset,'_pairs_varExp.png'), height=8.5, width=13, res=300, units="in")
mcmc_pairs(as.array(stan_toxic_Exp_out.rs), pars = c("C_0","k_1","k_2","BAF","delta","sigma_e[1]"),
           off_diag_args = list(size = 0.75)); dev.off()
png(paste0(dataset,'_pairs_varPow.png'), height=8.5, width=13, res=300, units="in")
mcmc_pairs(as.array(stan_toxic_Pow_out.rs), pars = c("C_0","k_1","k_2","BAF","delta","sigma_e[1]"),
           off_diag_args = list(size = 0.75)); dev.off()


# Save MCMC chains per model to file
pdf(paste0('Chains_',toxic,'_',dataset,'.pdf'), height=7, width=10) # ,units='in',res=600,compression='lzw'
traceplot(stan_toxic_out.rs, pars=c("C_0","k_1","k_2","BAF","sigma_e"))
traceplot(stan_toxic_Id_out.rs, pars=c("C_0","k_1","k_2","BAF","sigma_e"))
traceplot(stan_toxic_Fix_out.rs, pars=c("C_0","k_1","k_2","BAF","sigma_e"))
traceplot(stan_toxic_Exp_out.rs, pars=c("C_0","k_1","k_2","BAF","delta","sigma_e"))
traceplot(stan_toxic_Pow_out.rs, pars=c("C_0","k_1","k_2","BAF","delta","sigma_e"))
dev.off()


###
# generates a sequence of numbers that represent the time values in the data set
# The sequence starts at 0
# The sequence ends at the maximum value
t_range <- 0:max(toxkin_dat$Time)

# Create a data.frame named y_vs_t for the 5 different models, toxicant values are generated using the one_compartment() model function and 
  # the estimated RSTAN model coefficients are inputs along with the generated sequence of numbers and C_expsr (medium exposure concentration)
y_vs_t <- data.frame(var_t=t_range,
                     Homosk = one_compartment(summ_toxic_rs['C_0','mean'],summ_toxic_rs['k_1','mean'],summ_toxic_rs['k_2','mean'],C_expsr,t_range),
                     var_Id = one_compartment(summ_toxic_Id_rs['C_0','mean'],summ_toxic_Id_rs['k_1','mean'],summ_toxic_Id_rs['k_2','mean'],C_expsr,t_range),
                     var_Fix = one_compartment(summ_toxic_Fix_rs['C_0','mean'],summ_toxic_Fix_rs['k_1','mean'],summ_toxic_Fix_rs['k_2','mean'],C_expsr,t_range),
                     var_Exp = one_compartment(summ_toxic_Exp_rs['C_0','mean'],summ_toxic_Exp_rs['k_1','mean'],summ_toxic_Exp_rs['k_2','mean'],C_expsr,t_range),
                     var_Pow = one_compartment(summ_toxic_Pow_rs['C_0','mean'],summ_toxic_Pow_rs['k_1','mean'],summ_toxic_Pow_rs['k_2','mean'],C_expsr,t_range))

# Generate ggplot2 for homoscedastic model
p4 <- ggplot(data=toxkin_dat, aes(y=toxicant,x=Time)) + geom_point(color='black') +
  geom_line(data=y_vs_t, aes(var_t, Homosk, color = "black"), lwd=1) + theme(legend.position="none") + ggtitle("Homoscedastic variance") + 
  ylab(expression("Toxicant content in"~paste("ng/g"))) + xlab("")

# Generate ggplot2 for varId model
p8 <- ggplot(data=toxkin_dat, aes(y=toxicant,x=Time)) + geom_point(color='black') +
  geom_line(data=y_vs_t, aes(var_t, var_Id, color = "black"), lwd=1) + theme(legend.position="none") + ggtitle("Id variance") + 
  ylab(expression("Toxicant content in"~paste("ng/g"))) + xlab("")

# Generate ggplot2 for varFixed model
p5 <- ggplot(data=toxkin_dat, aes(y=toxicant,x=Time)) + geom_point(color='black') +
  geom_line(data=y_vs_t, aes(var_t, var_Fix, color = "black"), lwd=1) + theme(legend.position="none") + ggtitle("Fixed variance") + 
  ylab(expression("Toxicant content in"~paste("ng/g"))) + xlab("")

# Generate ggplot2 for varExp model
p6 <- ggplot(data=toxkin_dat, aes(y=toxicant,x=Time)) + geom_point(color='black') +
  geom_line(data=y_vs_t, aes(var_t, var_Exp, color = "black"), lwd=1) + theme(legend.position="none") + ggtitle("Exponential variance") + 
  ylab(expression("Toxicant content in"~paste("ng/g"))) + xlab("")

# Generate ggplot2 for varPower model
p7 <- ggplot(data=toxkin_dat, aes(y=toxicant,x=Time)) + geom_point(color='black') +
  geom_line(data=y_vs_t, aes(var_t, var_Pow, color = "black"), lwd=1) + theme(legend.position="none") + ggtitle("Power variance") + 
  ylab(expression("Toxicant content in"~paste("ng/g"))) + xlab("")

grid.arrange(p4, p8, p5, p6, p7,
             ncol=3, bottom="")
dev.off()

# Generate a plot with multiple lines, each representing the predictions from the 5 different models assigned to object p1.
# The ggplot() function uses the y_vs_t data.frame and sets the aesthetic mapping for the x-axis to t_range using aes(x = t_range).
# Then, to plot a line for each model, predicted values are taken from columns of the y_vs_t data.frame, the line color 
# is set to the variance functions, and the linewidth argument sets the line width
# theme(.) customizes the theme using: legend.position = c(0.775, 0.9) to places the legend at the specified coordinates within the plot area, 
# legend.text = element_text(size=11) sets the font size of the legend text to 11, legend.title = element_text("",size=0) sets the 
# legend title to an empty string with a size of 0 (effectively removing it), legend.margin=margin(0.1,0.1,0.1,0.1,'cm') sets the margin around the legend.
# ggtitle(paste0(toxic, " exposure")) sets the title of the plot to the string concatenated from the toxic variable and " exposure", 
# ylab(expression("Content in" ~ paste(mu, g, "/g"))) sets the y-axis label to "Content in Âµg/g" using the expression function to properly format the Greek letter mu (and the unit "g/g".

p1 <- ggplot(data=y_vs_t, aes(x=t_range)) + 
  geom_line(aes(y = Homosk, colour = "Homoscedastic"), linewidth=1) + 
  geom_line(aes(y = var_Id, colour = "var Id"), linewidth=1) +
  geom_line(aes(y = var_Fix, colour = "var Fixed"), linewidth=1) +
  geom_line(aes(y = var_Exp, colour = "var Exp"), linewidth=1) +
  geom_line(aes(y = var_Pow, colour = "var Pow"), linewidth=1) +
  theme(legend.position=c(0.775,0.9), legend.text=element_text(size=11), legend.title=element_text("",size=0), legend.margin=margin(0.1,0.1,0.1,0.1,'cm')) +
  ggtitle(paste0(toxic," exposure")) + ylab(expression("Content in"~paste(mu,g,"/g"))) + xlab("") #

# Generate a pdf file named according to the dataset variable with the suffix "_RSTAN.pdf". The pdf contains the plot p1 with additional black 
# points representing data from the toxkin_dat data.frame, where the toxicant variable is plotted against the Time variable. 
# The height and width of the pdf file are set to 6 and 7 inches, respectively. The dev.off() function closes and saves the pdf file.
pdf(paste0(toxic,'_',dataset,t_scale,'_RSTAN.pdf'), height=6, width=7)
p1 + geom_point(data=toxkin_dat, aes(y=toxicant,x=Time), color='black') 
dev.off()


# Compute LOOIC and log likelihood for the 5 different models
loo_homosk <- compute_loo(stan_toxic_out.rs)
loo_varId <- compute_loo(stan_toxic_Id_out.rs)
loo_varFix <- compute_loo(stan_toxic_Fix_out.rs)
loo_varExp <- compute_loo(stan_toxic_Exp_out.rs)
loo_varPow <- compute_loo(stan_toxic_Pow_out.rs)

# Compare models using logLik
loo_compare(loo_homosk,loo_varId)
loo_compare(loo_homosk,loo_varId,loo_varFix,loo_varExp,loo_varPow)

# Write logLik and LOOIC per model to file
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


# Red posteriors and posterior predictive distributions from files
mod_homsc <- read.csv(paste0(dataset,'_toxkin_STAN.csv'),row.names=1)
mod_hetsc <- read.csv(paste0(dataset,'_toxkin_Id_STAN.csv'),row.names=1)
mod_hetsc02 <- read.csv(paste0(dataset,'_toxkin_Exp_STAN.csv'),row.names=1)
mod_hetsc03 <- read.csv(paste0(dataset,'_toxkin_Pow_STAN.csv'),row.names=1)

# Set y-axis limit from 2.5% percentile posterior predictive distribution to 1.25 times maximum observed toxicant concentration in the organism
y_range <- c(min(mod_homsc[grepl("ypred",rownames(mod_homsc)),"X2.5."]),
             1.25*max(toxkin_dat$toxicant))
# Alternatively one may set the lower bound at zero using: y_range <- c(0, 1.25*max(toxkin_dat$toxicant))


# Save plot with posterior predictive distribution, predicted value and raw data points to file for the various models
pdf(paste0(dataset,"_PredInt_RSTAN.pdf"), width=11, height=5)
par(mfrow=c(1,2), mar=c(2.0,2.0,1.5,0), oma=c(1.5,2.25,0,0.5))
Plot_PostPredDistr(mod_homsc, mod_hetsc, raw_data=toxkin_dat, "Time","toxicant", paste0(toxic," content"), y_lim=y_range, main02="var Id")
Plot_PostPredDistr(mod_hetsc02, mod_hetsc03, raw_data=toxkin_dat, "Time","toxicant", paste0(toxic," content"), y_lim=y_range, main01="var het02", main02="var het03")
dev.off()
