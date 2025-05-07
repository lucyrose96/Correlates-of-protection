### Run parametric CoP model ###

#Note these simulations were run using parallel processing, with multiple files run simultaneously on different cores within a cluster. Each model also ran a separate chain on a different core (4 chains, 4 cores for each simulation). Parallel processing can also be used with the foreach() instead of the for() funtion.

#1. Setup 
#   Install packages and register cores for parallel processing
library(rstan)
library(foreach)
library(doParallel)
library(tidyverse)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
numCores <- 4
registerDoParallel(numCores)

#2. Define model,  data, and parameters (one dataset for each simulation)

##Model
model <- stan_model(file = 'Parametric_model.stan')

##Data
args <- commandArgs(trailingOnly=TRUE)
index=as.numeric(args[1])
data_file<-print(paste("c_",index,".csv", collapse = NULL, sep =""))
c<-read.csv(data_file)

##Parameters to take from model
pars=c("ps", "foi", "kk_in", "kk_pr", "n_50_in", "n_50_pr", "scale_in", "scale_pr", "VE_direct_in", "VE_direct_pr", "VE_indirect_in_mean", "VE_indirect_pr_mean", "p_direct_in", "p_direct_pr", "VE_in_mean", "VE_pr_mean", "VE_sym_mean", "VE_asym_mean")

##Parameters for extracting results values
lloq<- -5  
zero<-0
perc50<-c$perc50[1]
perc60<-c$perc60[1]
perc70<-c$perc70[1]
perc80<-c$perc80[1]
perc90<-c$perc90[1]
perc95<-c$perc95[1]
cop_ests<-c(lloq, perc50, perc60, zero, perc70, perc80, perc90, perc95)

#3. Run model
trial_dat_c <- list(N=nrow(c), 
                        Total=c$Total,
                        C=c$C,
                        A=c$A,
                        Z=c$s,
                        min_Z=min(c$s),
                        max_Z=max(c$s),
                        vaccinated=c$arm_01
                         
    )
    
fit<-sampling(object = model, data = trial_dat_c, chains = 4, cores = 4)
print(fit, pars = pars) 

#4. Extract estimates

##Parameter estimates
est_params<-data.frame(rstan::extract(fit, pars=pars))
 
##VE values at given CoP values     
cop_dat<-data.frame(matrix(NA, nrow = nrow(est_params)*length(cop_ests), ncol = 6))
colnames(cop_dat)<-c("cop_ests", "iter", "ve_in", "ve_pr", "ve_sym", "ve_asym")
cop_dat$cop_ests<-rep(cop_ests, nrow(est_params))
cop_dat$iter<-rep(seq(1:nrow(est_params)), each=length(cop_ests))

      
for(i in 1:nrow(cop_dat)){
        cop_dat$ve_in[i]<- est_params[cop_dat$iter[i], "p_direct_in"]*est_params[cop_dat$iter[i], "scale_in"] + ((1-est_params[cop_dat$iter[i], "p_direct_in"])*est_params[cop_dat$iter[i], "scale_in"])*(1/(1+exp(-est_params[cop_dat$iter[i], "kk_in"]*(cop_dat$cop_ests[i]-est_params[cop_dat$iter[i], "n_50_in"]))))
        cop_dat$ve_pr[i]<- est_params[cop_dat$iter[i], "p_direct_pr"]*est_params[cop_dat$iter[i], "scale_pr"] + ((1-est_params[cop_dat$iter[i], "p_direct_pr"])*est_params[cop_dat$iter[i], "scale_pr"])*(1/(1+exp(-est_params[cop_dat$iter[i], "kk_pr"]*(cop_dat$cop_ests[i]-est_params[cop_dat$iter[i], "n_50_pr"]))))
        cop_dat$ve_sym<- 1-(1-cop_dat$ve_in)*(1-cop_dat$ve_pr)
        cop_dat$ve_asym<- 1-((1-est_params[cop_dat$iter[i], "ps"])*(1-cop_dat$ve_in) + est_params[cop_dat$iter[i], "ps"]*(1-cop_dat$ve_in)*cop_dat$ve_pr) / (1-est_params[cop_dat$iter[i], "ps"])
        
      }
      
cop_dat_plot<-cop_dat %>%
        group_by(cop_ests)%>%
        summarise(
          mean_ve_in=mean(ve_in),
          mean_ve_pr=mean(ve_pr),
          mean_ve_sym=mean(ve_sym),
          mean_ve_asym=mean(ve_asym),
          
          lci_ve_in=quantile(ve_in, 0.025),
          lci_ve_pr=quantile(ve_pr, 0.025),
          lci_ve_sym=quantile(ve_sym, 0.025),
          lci_ve_asym=quantile(ve_asym, 0.025),
          
          uci_ve_in=quantile(ve_in, 0.975),
          uci_ve_pr=quantile(ve_pr, 0.975),
          uci_ve_sym=quantile(ve_sym, 0.975),
          uci_ve_asym=quantile(ve_asym, 0.975)
          
        )
 
#5. Output results
output_file=print(paste("median/results_para.", index, ".csv", collapse = NULL, sep = ""))
write.csv(cop_dat_plot, output_file, row.names = FALSE)

