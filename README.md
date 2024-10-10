Download the 6 .R files and start running from "Organisms_ToxKin.R". When working yourself through this R file, you read in toxicokinetic data for _Eisenia fetida_, _Gammarus Pulex_ and _Folsomia candida_ as well as the simulated datasets from either "Efetida_Zn.txt", "Gpulex_propanolol_vs_h.csv", "Fcandida_Cu_Ardestani2013_02.txt" or "simul_toxkin_2_---_---.csv". Running "simul_ToxKin.R" would perform a new data simulation. 

After reading in data, running "ToxKin_gnls_fun.R" that contains various functions and operations prepares for fitting frequentist models (i.e. non-Bayesian models) that assume homoscedastic or heteroscedastic residual variance. Then the models are fitted by running "ToxKin_gnls_run.R". People interested in Bayesian models should first run "ToxKin_rstan_fun.R" to initialize various functions that prepare for fitting. Then the models are fitted by running "ToxKin_rstan_run.R". 

Furthermore, note that the Apache 2.0 license applies to the use of the provided files. 
