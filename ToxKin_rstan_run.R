
# Prepare data for RSTAN
toxkin_dat <- data.frame(toxkin_dat, phase_no=NA) # Create extra column in data.frame for toxicokinetic phases. These should be integers in STAN
toxkin_dat$phase_no[toxkin_dat$phase == "Uptake"] <- 1 # Assign uptake phase to 1
toxkin_dat$phase_no[toxkin_dat$phase == "Elimination"] <- 2 # Assign elimination phase to 2

#

pars_no_d <- c("C_0","k","BAF","sigma_e","pred","res","ypred","log_lik") # Initialize parameters and generated quantities to be returned for Homoscedastic, varId and varFixed models
pars_d <- c("C_0","k","BAF","sigma_e","pred","res","ypred","log_lik") # Initialize parameters and generated quantities to be returned for varExp and varPower models

# Take means and standard errors from nls model as priors for mean and standard deviation for C_0, k_1 and k_2 parameters of RSTAN model fit 
nls_prior <- abs(summary(toxicant.nls)$coefficients[,1:2]); colnames(nls_prior) <- c("Value","Std.Error") 


# Fit Homoscedastic, varId, varFix, varExp and varPow models using RSTAN
stan_toxic_out.rs <- Run_RSTAN_ToxKin_mod(toxkin_dat,"toxicant","C_exp","Time",pars_no_d,nls_prior) 
stan_toxic_Id_out.rs <- Run_RSTAN_ToxKin_mod(toxkin_dat,"toxicant","C_exp","Time",pars_no_d,nls_prior, n_sigma=2)
stan_toxic_Fix_out.rs <- Run_RSTAN_ToxKin_mod(toxkin_dat,"toxicant","C_exp","Time",pars_no_d,nls_prior, varFixed=1)
stan_toxic_Exp_out.rs <- Run_RSTAN_ToxKin_mod(toxkin_dat,"toxicant","C_exp","Time",pars_d,nls_prior, varExponent=1); names(stan_toxic_Exp_out.rs)[4] <- "delta" # Change k[3] into delta
stan_toxic_Pow_out.rs <- Run_RSTAN_ToxKin_mod(toxkin_dat,"toxicant","C_exp","Time",pars_d,nls_prior, varPower=1); names(stan_toxic_Pow_out.rs)[4] <- "delta" # Change k[3] into delta

# If messages such as 
# "here are whatever error messages were returned
# [[1]]
# Stan model 'mod_OneComp_d' does not contain samples.
# [[2]]
# Stan model 'mod_OneComp_d' does not contain samples."
# or
# "Warning message:
# In .local(object, ...) :
# some chains had errors; consider specifying chains = 1 to debug"
# appear, then run a simulation with only 1 chain using:
# stan_toxic_Exp_out.rs <- Run_RSTAN_ToxKin_mod(toxkin_dat, "toxicant","C_exp","Time", pars_d, nls_prior, varExponent=1, n_chains=1); names(stan_toxic_Exp_out.rs)[4] <- "delta" # Change k[3] into delta


# Extract model parameter summaries

summ_toxic_rs <- summary(stan_toxic_out.rs,pars=pars_no_d[1:4])$summary 
summ_toxic_Id_rs <- summary(stan_toxic_Id_out.rs,pars=pars_no_d[1:4])$summary 
summ_toxic_Fix_rs <- summary(stan_toxic_Fix_out.rs,pars=pars_no_d[1:4])$summary 
summ_toxic_Exp_rs <- summary(stan_toxic_Exp_out.rs,pars=pars_d[1:4])$summary 
summ_toxic_Pow_rs <- summary(stan_toxic_Pow_out.rs,pars=pars_d[1:4])$summary 

# Print RSTAN fit model posteriors from summaries
print(summ_toxic_rs); print(summ_toxic_Id_rs); print(summ_toxic_Fix_rs); print(summ_toxic_Exp_rs); print(summ_toxic_Pow_rs)

# Write model parameter summaries to files
write.csv(summary(stan_toxic_out.rs)$summary,paste0(dataset,'_toxkin_STAN.csv')) 
write.csv(summary(stan_toxic_Id_out.rs)$summary,paste0(dataset,'_toxkin_Id_STAN.csv'))
write.csv(summary(stan_toxic_Fix_out.rs)$summary,paste0(dataset,'_toxkin_Fix_STAN.csv'))
write.csv(summary(stan_toxic_Exp_out.rs)$summary,paste0(dataset,'_toxkin_Exp_STAN.csv'))
write.csv(summary(stan_toxic_Pow_out.rs)$summary,paste0(dataset,'_toxkin_Pow_STAN.csv'))


# Write residual plots to file
x_range <- range(predict(toxicant.nls)); y_range <- 3 # Set ranges for x-axis (fitted value) based on nls model predicted values and y-axis (residuals) 
if (grepl("Fcandida",dataset)) x_range[1] <- 0  # Adjust x-axis lower limit to 0 for Fcandida datasets

# Save Figure with 2x3 panels to file
pdf(paste0('Res_STAN_',dataset,'.pdf'), height=6.5, width=10)
par(mfrow=c(2,3))
plot_RSTAN_res(stan_toxic_out.rs, hsk=FALSE, "Homoscedastic", y_lim=y_range, x_lim=x_range, y_ij=toxkin_dat$toxicant); print(Sys.time())
plot_RSTAN_res(stan_toxic_Id_out.rs, hsk="varId", "varId", y_lim=y_range, x_lim=x_range, y_ij=toxkin_dat$toxicant, y_lab="", phase=toxkin_dat$phase)
plot_RSTAN_res(stan_toxic_Fix_out.rs, hsk="varFix", "varFix", y_lim=y_range, x_lim=x_range, y_ij=toxkin_dat$toxicant)
plot_RSTAN_res(stan_toxic_Exp_out.rs, hsk="varExp", "varExp", y_lim=y_range, x_lim=x_range, y_lab="", x_lab="Fitted values", y_ij=toxkin_dat$toxicant)
plot_RSTAN_res(stan_toxic_Pow_out.rs, hsk="varPow", "varPow", y_lim=y_range, x_lim=x_range, y_lab="", y_ij=toxkin_dat$toxicant)
dev.off()


# Save parameter correlations to png files. Note that using pdf files instead may results in very large file sizes
png(paste0(dataset,'_pairs_Homosc.png'), height=8.5, width=13, res=300, units="in")
mcmc_pairs(as.array(stan_toxic_out.rs), pars = c("C_0","k[1]","k[2]","BAF","sigma_e[1]"),
           off_diag_args = list(size = 0.75)); dev.off()
png(paste0(dataset,'_pairs_varId.png'), height=8.5, width=13, res=300, units="in")
mcmc_pairs(as.array(stan_toxic_Id_out.rs), pars = c("C_0","k[1]","k[2]","BAF","sigma_e[1]","sigma_e[2]"),
           off_diag_args = list(size = 0.75)); dev.off()
png(paste0(dataset,'_pairs_varFix.png'), height=8.5, width=13, res=300, units="in")
mcmc_pairs(as.array(stan_toxic_Fix_out.rs), pars = c("C_0","k[1]","k[2]","BAF","sigma_e[1]"),
           off_diag_args = list(size = 0.75)); dev.off()
png(paste0(dataset,'_pairs_varExp.png'), height=8.5, width=13, res=300, units="in")
mcmc_pairs(as.array(stan_toxic_Exp_out.rs), pars = c("C_0","k[1]","k[2]","BAF","delta","sigma_e[1]"),
           off_diag_args = list(size = 0.75)); dev.off()
png(paste0(dataset,'_pairs_varPow.png'), height=8.5, width=13, res=300, units="in")
mcmc_pairs(as.array(stan_toxic_Pow_out.rs), pars = c("C_0","k[1]","k[2]","BAF","delta","sigma_e[1]"),
           off_diag_args = list(size = 0.75)); dev.off()


# Save MCMC chains per model to file
pdf(paste0('Chains_',toxic,'_',dataset,'.pdf'), height=7, width=10) # ,units='in',res=600,compression='lzw'
traceplot(stan_toxic_out.rs, pars=c("C_0","k[1]","k[2]","BAF","sigma_e"))
traceplot(stan_toxic_Id_out.rs, pars=c("C_0","k[1]","k[2]","BAF","sigma_e"))
traceplot(stan_toxic_Fix_out.rs, pars=c("C_0","k[1]","k[2]","BAF","sigma_e"))
traceplot(stan_toxic_Exp_out.rs, pars=c("C_0","k[1]","k[2]","BAF","delta","sigma_e"))
traceplot(stan_toxic_Pow_out.rs, pars=c("C_0","k[1]","k[2]","BAF","delta","sigma_e"))
dev.off()


###
# generates a sequence of numbers that represent the time values in the data set
# The sequence starts at 0
# The sequence ends at the maximum value
t_range <- 0:max(toxkin_dat$Time)

# Create a data.frame named y_vs_t for the 5 different models, toxicant values are generated using the one_compartment() model function and 
# the estimated RSTAN model coefficients are inputs along with the generated sequence of numbers and C_expsr (medium exposure concentration)
y_vs_t <- data.frame(var_t=t_range,
                     Homosk = one_compartment(summ_toxic_rs['C_0','mean'],summ_toxic_rs['k[1]','mean'],summ_toxic_rs['k[2]','mean'],C_expsr,t_range),
                     var_Id = one_compartment(summ_toxic_Id_rs['C_0','mean'],summ_toxic_Id_rs['k[1]','mean'],summ_toxic_Id_rs['k[2]','mean'],C_expsr,t_range),
                     var_Fix = one_compartment(summ_toxic_Fix_rs['C_0','mean'],summ_toxic_Fix_rs['k[1]','mean'],summ_toxic_Fix_rs['k[2]','mean'],C_expsr,t_range),
                     var_Exp = one_compartment(summ_toxic_Exp_rs['C_0','mean'],summ_toxic_Exp_rs['k[1]','mean'],summ_toxic_Exp_rs['k[2]','mean'],C_expsr,t_range),
                     var_Pow = one_compartment(summ_toxic_Pow_rs['C_0','mean'],summ_toxic_Pow_rs['k[1]','mean'],summ_toxic_Pow_rs['k[2]','mean'],C_expsr,t_range))

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

# Note that the message:
# "Some Pareto k diagnostic values are too high. See help('pareto-k-diagnostic') for details."
# is a warning related to the computation of the WAIC, not an error. 
# If this warning message is shown, then only use the LOOIC

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


# Read posteriors and posterior predictive distributions from files
mod_homsc <- read.csv(paste0(dataset,'_toxkin_STAN.csv'),row.names=1)
mod_hetsc02 <- read.csv(paste0(dataset,'_toxkin_Id_STAN.csv'),row.names=1)
mod_hetsc03 <- read.csv(paste0(dataset,'_toxkin_Fix_STAN.csv'),row.names=1)
mod_hetsc04 <- read.csv(paste0(dataset,'_toxkin_Exp_STAN.csv'),row.names=1)
mod_hetsc05 <- read.csv(paste0(dataset,'_toxkin_Pow_STAN.csv'),row.names=1)

# Set y-axis limit from 2.5% percentile posterior predictive distribution to 1.25 times maximum observed toxicant concentration in the organism
y_range <- c(min(mod_homsc[grepl("ypred",rownames(mod_homsc)),"X2.5."]),
             1.25*max(toxkin_dat$toxicant))
# Alternatively one may set the lower bound at zero using: y_range <- c(0, 1.25*max(toxkin_dat$toxicant))


# Save plot with posterior predictive distribution, predicted value and raw data points to file for the various models
pdf(paste0(dataset,"_PredInt_RSTAN.pdf"), width=11, height=5)
par(mfrow=c(1,2), mar=c(2.0,2.0,1.5,0), oma=c(1.5,2.25,0,0.5))
Plot_PostPredDistr(mod_homsc, mod_hetsc02, raw_data=toxkin_dat, "Time","toxicant", paste0(toxic," content"), y_lim=y_range, main02="var Id")
Plot_PostPredDistr(mod_hetsc04, mod_hetsc05, raw_data=toxkin_dat, "Time","toxicant", paste0(toxic," content"), y_lim=y_range, main01="var Exp", main02="var Power")
dev.off()
