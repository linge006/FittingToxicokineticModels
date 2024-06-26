
# Load nlme library
library(nlme)

### one_compartment() function ###
# This function is used for fitting the one-compartmental non-linear toxicokinetic model
# The function predicts the toxicant concentration in the organism that is studied
# The function requires starting values for the C_0, k_1 and k_2 parameters to be estimated
# The function requires C_exposure (concentration of the toxicant in the medium) 
# The function requires a value for the t_e parameter representing the timepoint that the organism were transferred from exposed to clean media
# The function requires time values associated with the study design on which toxicant concentrations in organism and exposure medium were measured

one_compartment <- function(C_0,k_1,k_2,C_exposure,time,t_e=t_d){

  ifelse(time <= t_e, 
         C_0+k_1/k_2*C_exposure*(1-exp(-k_2*time)), 
         C_0+k_1/k_2*C_exposure*(exp(-k_2*(time-t_e))-exp(-k_2*time)))
  
} # End one_compartment()


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


# Resurces for computing prediction intervals:
  # https://stackoverflow.com/questions/14358811/extract-prediction-band-from-lme-fit
  # https://rpubs.com/Daniel_He/1041063
