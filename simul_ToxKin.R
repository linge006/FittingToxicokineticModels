
t_d <- 14; 

### ### ### ### ### ### ### ### ### ### ### ### ###
### Case 1 - No plateaus and Heteroscedasticity ###
### ### ### ### ### ### ### ### ### ### ### ### ### 

t_end <- 28; repl <- 4

toxkin_dat <- data.frame(Time=rep(c(0:t_end), each=repl), k_1=runif(repl+t_end*repl,0.5,1.5), k_2=runif(repl+t_end*repl,0.05,0.15), 
                         C_exp=rep(1,repl+t_end*repl), t_d=rep(t_d,repl+t_end*repl), error=runif(repl+t_end*repl,-2,2))

#toxkin_dat <- data.frame(toxkin_dat, toxicant=with(toxkin_dat, one_compartment(1, k_1, k_2, C_exp, Time, t_d) + error)) # Randomize/Simulate "fixed-effects" parameters as well

toxkin_dat <- data.frame(toxkin_dat, pred_toxicant=with(toxkin_dat, one_compartment(2, 1, 0.1, 1, Time, t_d) )) # C_0=2, k_1=1, k_2=0.1, C_exp=1
toxkin_dat <- data.frame(toxkin_dat, toxicant=with(toxkin_dat, one_compartment(2, 1, 0.1, 1, Time, t_d) + 
                                                     rnorm(nrow(toxkin_dat), 0, pred_toxicant*0.2) ))

y_max <- max(toxkin_dat$toxicant)
plot(toxicant~Time, data=toxkin_dat, pch=20, ylim=c(0,y_max))

var_names <- c("Time","toxicant","C_exp","t_d")
red_dat <- c(1,2,7,14:16,21,28)

write.csv(toxkin_dat[,var_names], "simul_toxkin_1_hete_full.csv", row.names=F)
write.csv(toxkin_dat[toxkin_dat$Time%in%red_dat,var_names], "simul_toxkin_1_hete_red.csv", row.names=F)

### ### ### ### ### ### ### ### ### ### ### ### ###
### Case 1 - No plateaus and Homoscedasticity   ###
### ### ### ### ### ### ### ### ### ### ### ### ### 

toxkin_dat$toxicant <- toxkin_dat$pred_toxicant+toxkin_dat$error

plot(toxicant~Time, data=toxkin_dat, pch=20, ylim=c(0,y_max), col="magenta")

write.csv(toxkin_dat[,var_names], "simul_toxkin_1_homo_full.csv", row.names=F)
write.csv(toxkin_dat[toxkin_dat$Time%in%red_dat,var_names], "simul_toxkin_1_homo_red.csv", row.names=F)


### ### ### ### ### ### ### ### ### ### ### ### 
### Case 2 - Plateau and Heteroscedasticty  ###
### ### ### ### ### ### ### ### ### ### ### ### 

toxkin_dat <- data.frame(Time=rep(c(0:t_end), each=repl), k_1=runif(repl+t_end*repl,4,6), k_2=runif(repl+t_end*repl,0.5,1.5), 
                         C_exp=rep(1,repl+t_end*repl), t_d=rep(t_d,repl+t_end*repl), error=runif(repl+t_end*repl,-2,2))
toxkin_dat <- data.frame(toxkin_dat, pred_toxicant=with(toxkin_dat, one_compartment(2, 5, 1, 1, Time, t_d) )) # C_0=1, k_1=1, k_2=0.1, C_exp=1
toxkin_dat <- data.frame(toxkin_dat, toxicant=with(toxkin_dat, one_compartment(2, 5, 1, 1, Time, t_d) + 
                                                     rnorm(nrow(toxkin_dat), 0, pred_toxicant*0.2) ))

y_max <- max(toxkin_dat$toxicant)
plot(toxicant~Time, data=toxkin_dat, pch=20, ylim=c(0,y_max))

write.csv(toxkin_dat[,var_names], "simul_toxkin_2_hete_full.csv", row.names=F)
write.csv(toxkin_dat[toxkin_dat$Time%in%red_dat,var_names], "simul_toxkin_2_hete_red.csv", row.names=F)

### ### ### ### ### ### ### ### ### ### ### ### 
### Case 2 - Plateau and Homoscedasticty    ###
### ### ### ### ### ### ### ### ### ### ### ### 

toxkin_dat$toxicant <- toxkin_dat$pred_toxicant+toxkin_dat$error

plot(toxicant~Time, data=toxkin_dat, pch=20, ylim=c(0,y_max), col="magenta")

write.csv(toxkin_dat[,var_names], "simul_toxkin_2_homo_full.csv", row.names=F)
write.csv(toxkin_dat[toxkin_dat$Time%in%red_dat,var_names], "simul_toxkin_2_homo_red.csv", row.names=F)


### ### ### ### ### ### ### ### ### ### ### ### 
### Case 3 - Plateau02 and Heteroscedasticty  ###
### ### ### ### ### ### ### ### ### ### ### ### 

toxkin_dat <- data.frame(Time=rep(c(0:t_end),each=repl), k_1=runif(repl+t_end*repl,8,12), k_2=runif(repl+t_end*repl,0.5,1.5), 
                         C_exp=rep(1,repl+t_end*repl), t_d=rep(t_d,repl+t_end*repl), error=runif(repl+t_end*repl,-3,3))
toxkin_dat <- data.frame(toxkin_dat, pred_toxicant=with(toxkin_dat, one_compartment(5, 10, 1, 1, Time, t_d) )) # C_0=1, k_1=1, k_2=0.1, C_exp=1
toxkin_dat <- data.frame(toxkin_dat, toxicant=with(toxkin_dat, one_compartment(5, 10, 1, 1, Time, t_d) + 
                                                     rnorm(nrow(toxkin_dat), 0, pred_toxicant*0.2) ))

y_max <- max(toxkin_dat$toxicant)
plot(toxicant~Time, data=toxkin_dat, pch=20, ylim=c(0,y_max))

### ### ### ### ### ### ### ### ### ### ### ### 
### Case 3 - Plateau02 and Homoscedasticty  ###
### ### ### ### ### ### ### ### ### ### ### ### 

toxkin_dat$toxicant <- toxkin_dat$pred_toxicant+toxkin_dat$error

plot(toxicant~Time, data=toxkin_dat, pch=20, ylim=c(0,y_max), col='magenta')
