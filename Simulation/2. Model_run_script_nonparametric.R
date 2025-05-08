### Run non-parametric CoP model ###

#Note these simulations were run using parallel processing, with multiple files run simultaneously on different cores within a cluster. Each model also ran a separate chain on a different core (4 chains, 4 cores for each simulation). Parallel processing can also be used with the foreach() instead of the for() funtion.

#1. Setup 
#   Install packages and register cores for parallel processing
library(rstan)
library(foreach)
library(doParallel)
library(tidyverse)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = FALSE)
numCores <- 4
registerDoParallel(numCores)

#2. Define model,  data, and parameters (one dataset for each simulation)

##Model
model <- stan_model(file = 'Nonparametric_model.stan')

##Data
args <- commandArgs(trailingOnly=TRUE)
index=as.numeric(args[1])
data_file<-print(paste("c_",index,".csv", collapse = NULL, sep =""))#loading one simulation dataset for each model run 
c<-read.csv(data_file)

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

#neuts = CoP value sequence (e.g. neutralising antibodies, at intervals of 0.1), Z_range = range of spline positions (Z gives the position along the X axis), dX = value of intervals along the spline, Z_pos = spline position.
neuts<-seq(min(c$s[c$arm=="V"]), max(c$s[c$arm=="V"]),0.1)
Z_range<-seq(min(neuts),max(neuts),0.1)
dX<-0.1
Z_pos<-floor((c$s-Z_range[1])/dX)
Z_pos[Z_pos<1]<-1
Z_pos[Z_pos>length(Z_range)]<-length(Z_range)
num_subjects<-nrow(c)

trial_dat <- list(N=num_subjects, 
                  C=c$C,
                  A=c$A,
                  vaccinated=c$arm_01,
                  Z=c$s,
                  num_knots=5,
                  knots=unname(quantile(c$s[c$arm_01==1],probs=seq(from=0, to=1, length.out = 5))), 
                  spline_degree=3,
                  X=Z_range, 
                  spline_length=length(Z_range),
                  dX=dX,
                  Z_index=Z_pos
)

fit <- sampling(object = model, data = trial_dat, chains = 4, cores = 4)


#4. Extract VE estimates

##VEin
cVE_in_hat<-rstan::extract(fit)[["cVE_hat"]]
cVE_in_hat_mean<-apply(cVE_in_hat,2,mean)
cVE_in_hat_interval<-apply(cVE_in_hat,2,quantile,c(0.025,0.975))
##VEpr
cVE_pr_hat<-rstan::extract(fit)[["cVE_pr_hat"]]
cVE_pr_hat_mean<-apply(cVE_pr_hat,2,mean)
cVE_pr_hat_interval<-apply(cVE_pr_hat,2,quantile,c(0.025,0.975))
##VEsym
cVE_sym_hat<-1-(1-cVE_in_hat)*(1-cVE_pr_hat)
cVE_sym_hat_mean<-apply(cVE_sym_hat,2,mean)
cVE_sym_hat_interval<-apply(cVE_sym_hat,2,quantile,c(0.025,0.975))
##VEasym
ps<-rstan::extract(fit)[["ps"]]
ps_mean<-mean(ps)
cVE_asym_hat<-   1-((1-ps_mean)*(1-cVE_in_hat) + ps_mean*(1-cVE_in_hat)*cVE_pr_hat) / (1-ps_mean)
cVE_asym_hat_mean<-apply(cVE_asym_hat,2,mean)
cVE_asym_hat_interval<-apply(cVE_asym_hat,2,quantile,c(0.025,0.975))

##Combine to dataset
plot_dat<-data.frame(cbind(neuts,cVE_in_hat_mean,t(cVE_in_hat_interval), cVE_pr_hat_mean,t(cVE_pr_hat_interval), cVE_sym_hat_mean,t(cVE_sym_hat_interval), cVE_asym_hat_mean,t(cVE_asym_hat_interval)))
names(plot_dat)<-c("cop_ests","mean_ve_in","lci_ve_in","uci_ve_in","mean_ve_pr","lci_ve_pr","uci_ve_pr","mean_ve_sym","lci_ve_sym","uci_ve_sym","mean_ve_asym","lci_ve_asym","uci_ve_asym")


#5. Extract estimates - VE values at given CoP values   

##The non-parametric model estimates VE at CoP value intervals of 0.1. The closest one to each of the key parameters is taken here.
plot_dat_lloq<-plot_dat[which(plot_dat$cop_ests >= lloq-0.05 & plot_dat$cop_ests <=  lloq+0.05),]
plot_dat_lloq<-plot_dat_lloq[1,]
plot_dat_lloq$cop<-"LLOQ"

plot_dat_zero<-plot_dat[which(plot_dat$cop_ests >= zero-0.05 & plot_dat$cop_ests <=  zero+0.05),]
plot_dat_zero<-plot_dat_zero[1,]
plot_dat_zero$cop<-"zero"

plot_dat_50perc<-plot_dat[which(plot_dat$cop_ests >= perc50-0.05 & plot_dat$cop_ests <=  perc50+0.05),]
plot_dat_50perc<-plot_dat_50perc[1,]
plot_dat_50perc$cop<-"50 percentile"

plot_dat_60perc<-plot_dat[which(plot_dat$cop_ests >= perc60-0.05 & plot_dat$cop_ests <=  perc60+0.05),]
plot_dat_60perc<-plot_dat_60perc[1,]
plot_dat_60perc$cop<-"60 percentile"

plot_dat_70perc<-plot_dat[which(plot_dat$cop_ests >= perc70-0.05 & plot_dat$cop_ests <=  perc70+0.05),]
plot_dat_70perc<-plot_dat_70perc[1,]
plot_dat_70perc$cop<-"70 percentile"

plot_dat_80perc<-plot_dat[which(plot_dat$cop_ests >= perc80-0.05 & plot_dat$cop_ests <=  perc80+0.05),]
plot_dat_80perc<-plot_dat_80perc[1,]
plot_dat_80perc$cop<-"80 percentile"

plot_dat_90perc<-plot_dat[which(plot_dat$cop_ests >= perc90-0.05 & plot_dat$cop_ests <=  perc90+0.05),]
plot_dat_90perc<-plot_dat_90perc[1,]
plot_dat_90perc$cop<-"90 percentile"

plot_dat_95perc<-plot_dat[which(plot_dat$cop_ests >= perc95-0.05 & plot_dat$cop_ests <=  perc95+0.05),]
plot_dat_95perc<-plot_dat_95perc[1,]
plot_dat_95perc$cop<-"95 percentile"

cop_dat_plot<-rbind(plot_dat_lloq, plot_dat_zero, plot_dat_50perc, plot_dat_60perc, plot_dat_70perc, plot_dat_80perc, plot_dat_90perc, plot_dat_95perc)
cop_dat_plot$index<-index
cop_dat_plot

#6. Output results
output_file=print(paste("median/results_nonpara.", index, ".csv", collapse = NULL, sep = ""))
write.csv(cop_dat_plot, output_file, row.names = FALSE)

