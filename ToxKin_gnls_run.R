
# Fitting of a generalized nonlinear least squares (gnls) models to the data
# variables and parameters are as for the nls model fit
# In addition, the weights argument determines the residual variance function. NULL indicates homoscedastic residual variance, 
# whereas the varIdent(), varFixed(), varExp() and varPower() functions are the considered heteroscedastic variance assumptions 
# The estimated parameters returned from the nls() model are used as initial values for the gnls() model fits

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
