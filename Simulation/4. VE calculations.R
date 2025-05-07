# VE calculations - VE at given CoP values

#First run "1. Simulation_setup" to generate the dataset "c" and all the true parameter estimates 

#1. Setup
library(tidyverse)
lci<-function(x){quantile(x, 0.025)}
uci<-function(x){quantile(x, 0.975)}

#2. Download results (parametric or non-parametric)
#cop_dat_plot<-read.csv("results_p.csv")
cop_dat_plot<-read.csv("results_np.csv")


#3. Set CoP values to estimate at 
lloq<- c$lloq[1]  
zero<- c$zero[1]
perc50<- c$perc50[1] #percentiles across all simulations
perc60<- c$perc60[1]
perc70<- c$perc70[1]
perc80<- c$perc80[1]
perc90<- c$perc90[1]
perc95<- c$perc95[1]

#4. Get the true VE estimates at these CoP values
##Extract
(true_vein_lloq<-c$VE_direct_in[1])
(true_vepr_lloq<-c$VE_direct_pr[1])
true_vein_zero<-  (p_direct_in*scale_vein+((1-p_direct_in)*scale_vein)*(((1/(1+exp(-kk_in*(zero-n_50_in)))))))
true_vepr_zero<-  (p_direct_pr*scale_vepr+((1-p_direct_pr)*scale_vepr)*(((1/(1+exp(-kk_pr*(zero-n_50_pr)))))))
true_vein_50perc<-  (p_direct_in*scale_vein+((1-p_direct_in)*scale_vein)*(((1/(1+exp(-kk_in*(perc50-n_50_in)))))))
true_vepr_50perc<-  (p_direct_pr*scale_vepr+((1-p_direct_pr)*scale_vepr)*(((1/(1+exp(-kk_pr*(perc50-n_50_pr)))))))
true_vein_60perc<-  (p_direct_in*scale_vein+((1-p_direct_in)*scale_vein)*(((1/(1+exp(-kk_in*(perc60-n_50_in)))))))
true_vepr_60perc<-  (p_direct_pr*scale_vepr+((1-p_direct_pr)*scale_vepr)*(((1/(1+exp(-kk_pr*(perc60-n_50_pr)))))))
true_vein_70perc<-  (p_direct_in*scale_vein+((1-p_direct_in)*scale_vein)*(((1/(1+exp(-kk_in*(perc70-n_50_in)))))))
true_vepr_70perc<-  (p_direct_pr*scale_vepr+((1-p_direct_pr)*scale_vepr)*(((1/(1+exp(-kk_pr*(perc70-n_50_pr)))))))
true_vein_80perc<-  (p_direct_in*scale_vein+((1-p_direct_in)*scale_vein)*(((1/(1+exp(-kk_in*(perc80-n_50_in)))))))
true_vepr_80perc<-  (p_direct_pr*scale_vepr+((1-p_direct_pr)*scale_vepr)*(((1/(1+exp(-kk_pr*(perc80-n_50_pr)))))))
true_vein_90perc<-  (p_direct_in*scale_vein+((1-p_direct_in)*scale_vein)*(((1/(1+exp(-kk_in*(perc90-n_50_in)))))))
true_vepr_90perc<-  (p_direct_pr*scale_vepr+((1-p_direct_pr)*scale_vepr)*(((1/(1+exp(-kk_pr*(perc90-n_50_pr)))))))
true_vein_95perc<-  (p_direct_in*scale_vein+((1-p_direct_in)*scale_vein)*(((1/(1+exp(-kk_in*(perc95-n_50_in)))))))
true_vepr_95perc<-  (p_direct_pr*scale_vepr+((1-p_direct_pr)*scale_vepr)*(((1/(1+exp(-kk_pr*(perc95-n_50_pr)))))))

##Round
decimals=2

true_vein_lloq<-round(true_vein_lloq, digits = decimals)
true_vein_zero<-round(true_vein_zero, digits = decimals)
true_vein_50perc<-round(true_vein_50perc, digits = decimals)
true_vein_60perc<-round(true_vein_60perc, digits = decimals)
true_vein_70perc<-round(true_vein_70perc, digits = decimals)
true_vein_80perc<-round(true_vein_80perc, digits = decimals)
true_vein_90perc<-round(true_vein_90perc, digits = decimals)
true_vein_95perc<-round(true_vein_95perc, digits = decimals)

true_vepr_lloq<-round(true_vepr_lloq, digits = decimals)
true_vepr_zero<-round(true_vepr_zero, digits = decimals)
true_vepr_50perc<-round(true_vepr_50perc, digits = decimals)
true_vepr_60perc<-round(true_vepr_60perc, digits = decimals)
true_vepr_70perc<-round(true_vepr_70perc, digits = decimals)
true_vepr_80perc<-round(true_vepr_80perc, digits = decimals)
true_vepr_90perc<-round(true_vepr_90perc, digits = decimals)
true_vepr_95perc<-round(true_vepr_95perc, digits = decimals)


cop_dat_plot$lci_ve_in<-round(cop_dat_plot$lci_ve_in, digits=decimals)
cop_dat_plot$lci_ve_pr<-round(cop_dat_plot$lci_ve_pr, digits=decimals)
cop_dat_plot$lci_ve_sym<-round(cop_dat_plot$lci_ve_sym, digits=decimals)
cop_dat_plot$lci_ve_asym<-round(cop_dat_plot$lci_ve_asym, digits=decimals)
cop_dat_plot$uci_ve_in<-round(cop_dat_plot$uci_ve_in, digits=decimals)
cop_dat_plot$uci_ve_pr<-round(cop_dat_plot$uci_ve_pr, digits=decimals)
cop_dat_plot$uci_ve_sym<-round(cop_dat_plot$uci_ve_sym, digits=decimals)
cop_dat_plot$uci_ve_asym<-round(cop_dat_plot$uci_ve_asym, digits=decimals)


##Set names (8 parameters to estimate)
cop_dat_plot$cop<-"LLOQ"
cop_dat_plot$cop[seq(2,nrow(cop_dat_plot)-(8-2), 8)]<-"median"
cop_dat_plot$cop[seq(3,nrow(cop_dat_plot)-(8-3), 8)]<-"60 percentile"
cop_dat_plot$cop[seq(4,nrow(cop_dat_plot)-(8-4), 8)]<-"zero"
cop_dat_plot$cop[seq(5,nrow(cop_dat_plot)-(8-5), 8)]<-"70 percentile"
cop_dat_plot$cop[seq(6,nrow(cop_dat_plot)-(8-6), 8)]<-"80 percentile"
cop_dat_plot$cop[seq(7,nrow(cop_dat_plot)-(8-7), 8)]<-"90 percentile"
cop_dat_plot$cop[seq(8,nrow(cop_dat_plot)-(8-8), 8)]<-"95 percentile"
table(cop_dat_plot$cop) 

#5. Summarise VE at each CoP value
##LLOQ
(cop_dat_plot_overall_lloq<-cop_dat_plot %>%
    subset(cop== "LLOQ")%>%
    summarise(
      
      mean_vein=mean(mean_ve_in),
      mean_vepr=mean(mean_ve_pr),
      
      lci_vein=lci(mean_ve_in),
      lci_vepr=lci(mean_ve_pr),
      
      uci_vein=uci(mean_ve_in),
      uci_vepr=uci(mean_ve_pr),
      
      cov_in=sum(true_vein_lloq >= lci_ve_in & true_vein_lloq <= uci_ve_in )/n(),
      cov_pr=sum(true_vepr_lloq >= lci_ve_pr & true_vepr_lloq <= uci_ve_pr )/n(),
      
      RMSE_in= sqrt(sum((true_vein_lloq - mean_ve_in)^2)/n()),
      RMSE_pr= sqrt(sum((true_vepr_lloq - mean_ve_pr)^2)/n()),
      
      RRMSE_in= sqrt(sum(((true_vein_lloq - mean_ve_in)/true_vein_lloq)^2)/n()),
      RRMSE_pr= sqrt(sum(((true_vepr_lloq - mean_ve_pr)/true_vepr_lloq)^2)/n()),
      
    ))


##Zero
(cop_dat_plot_overall_zero<-cop_dat_plot %>%
    subset(cop== "zero")%>%
    summarise(
      
      mean_vein=mean(mean_ve_in),
      mean_vepr=mean(mean_ve_pr),
      
      lci_vein=lci(mean_ve_in),
      lci_vepr=lci(mean_ve_pr),
      
      uci_vein=uci(mean_ve_in),
      uci_vepr=uci(mean_ve_pr),
      
      cov_in=sum(true_vein_zero >= lci_ve_in & true_vein_zero <= uci_ve_in )/n(),
      cov_pr=sum(true_vepr_zero >= lci_ve_pr & true_vepr_zero <= uci_ve_pr )/n(),
      
      RMSE_in= sqrt(sum((true_vein_zero - mean_ve_in)^2)/n()),
      RMSE_pr= sqrt(sum((true_vepr_zero - mean_ve_pr)^2)/n()),
      
      RRMSE_in= sqrt(sum(((true_vein_zero - mean_ve_in)/true_vein_zero)^2)/n()),
      RRMSE_pr= sqrt(sum(((true_vepr_zero - mean_ve_pr)/true_vepr_zero)^2)/n()),
      
    ))

##50th percentile
(cop_dat_plot_overall_50perc<-cop_dat_plot %>%
    subset(cop== "median")%>%
    summarise(
      
      mean_vein=mean(mean_ve_in),
      mean_vepr=mean(mean_ve_pr),
      
      lci_vein=lci(mean_ve_in),
      lci_vepr=lci(mean_ve_pr),
      
      uci_vein=uci(mean_ve_in),
      uci_vepr=uci(mean_ve_pr),
      
      
      cov_in=sum(true_vein_50perc >= lci_ve_in & true_vein_50perc <= uci_ve_in )/n(),
      cov_pr=sum(true_vepr_50perc >= lci_ve_pr & true_vepr_50perc <= uci_ve_pr )/n(),
      
      RMSE_in= sqrt(sum((true_vein_50perc - mean_ve_in)^2)/n()),
      RMSE_pr= sqrt(sum((true_vepr_50perc - mean_ve_pr)^2)/n()),
      
      RRMSE_in= sqrt(sum(((true_vein_50perc - mean_ve_in)/true_vein_50perc)^2)/n()),
      RRMSE_pr= sqrt(sum(((true_vepr_50perc - mean_ve_pr)/true_vepr_50perc)^2)/n()),
      
    ))

##60th percentile
(cop_dat_plot_overall_60perc<-cop_dat_plot %>%
    subset(cop== "60 percentile")%>%
    summarise(
      
      mean_vein=mean(mean_ve_in),
      mean_vepr=mean(mean_ve_pr),
      
      lci_vein=lci(mean_ve_in),
      lci_vepr=lci(mean_ve_pr),
      
      uci_vein=uci(mean_ve_in),
      uci_vepr=uci(mean_ve_pr),
      
      
      cov_in=sum(true_vein_60perc >= lci_ve_in & true_vein_60perc <= uci_ve_in )/n(),
      cov_pr=sum(true_vepr_60perc >= lci_ve_pr & true_vepr_60perc <= uci_ve_pr )/n(),
      
      RMSE_in= sqrt(sum((true_vein_60perc - mean_ve_in)^2)/n()),
      RMSE_pr= sqrt(sum((true_vepr_60perc - mean_ve_pr)^2)/n()),
      
      RRMSE_in= sqrt(sum(((true_vein_60perc - mean_ve_in)/true_vein_60perc)^2)/n()),
      RRMSE_pr= sqrt(sum(((true_vepr_60perc - mean_ve_pr)/true_vepr_60perc)^2)/n()),
      
    ))

##70th percentile
(cop_dat_plot_overall_70perc<-cop_dat_plot %>%
    subset(cop== "70 percentile")%>%
    summarise(
      
      
      mean_vein=mean(mean_ve_in),
      mean_vepr=mean(mean_ve_pr),
      
      lci_vein=lci(mean_ve_in),
      lci_vepr=lci(mean_ve_pr),
      
      uci_vein=uci(mean_ve_in),
      uci_vepr=uci(mean_ve_pr),
      
      
      cov_in=sum(true_vein_70perc >= lci_ve_in & true_vein_70perc <= uci_ve_in )/n(),
      cov_pr=sum(true_vepr_70perc >= lci_ve_pr & true_vepr_70perc <= uci_ve_pr )/n(),
      
      RMSE_in= sqrt(sum((true_vein_70perc - mean_ve_in)^2)/n()),
      RMSE_pr= sqrt(sum((true_vepr_70perc - mean_ve_pr)^2)/n()),
      
      RRMSE_in= sqrt(sum(((true_vein_70perc - mean_ve_in)/true_vein_70perc)^2)/n()),
      RRMSE_pr= sqrt(sum(((true_vepr_70perc - mean_ve_pr)/true_vepr_70perc)^2)/n()),
      
    ))

##80th percentile
(cop_dat_plot_overall_80perc<-cop_dat_plot %>%
    subset(cop== "80 percentile")%>%
    summarise(
      
      mean_vein=mean(mean_ve_in),
      mean_vepr=mean(mean_ve_pr),
      
      lci_vein=lci(mean_ve_in),
      lci_vepr=lci(mean_ve_pr),
      
      uci_vein=uci(mean_ve_in),
      uci_vepr=uci(mean_ve_pr),
      
      
      cov_in=sum(true_vein_80perc >= lci_ve_in & true_vein_80perc <= uci_ve_in )/n(),
      cov_pr=sum(true_vepr_80perc >= lci_ve_pr & true_vepr_80perc <= uci_ve_pr )/n(),
      
      RMSE_in= sqrt(sum((true_vein_80perc - mean_ve_in)^2)/n()),
      RMSE_pr= sqrt(sum((true_vepr_80perc - mean_ve_pr)^2)/n()),
      
      RRMSE_in= sqrt(sum(((true_vein_80perc - mean_ve_in)/true_vein_80perc)^2)/n()),
      RRMSE_pr= sqrt(sum(((true_vepr_80perc - mean_ve_pr)/true_vepr_80perc)^2)/n()),
      
    ))

##90th percentile
(cop_dat_plot_overall_90perc<-cop_dat_plot %>%
    subset(cop== "90 percentile")%>%
    summarise(
      
      mean_vein=mean(mean_ve_in),
      mean_vepr=mean(mean_ve_pr),
      
      lci_vein=lci(mean_ve_in),
      lci_vepr=lci(mean_ve_pr),
      
      uci_vein=uci(mean_ve_in),
      uci_vepr=uci(mean_ve_pr),
      
      
      cov_in=sum(true_vein_90perc >= lci_ve_in & true_vein_90perc <= uci_ve_in )/n(),
      cov_pr=sum(true_vepr_90perc >= lci_ve_pr & true_vepr_90perc <= uci_ve_pr )/n(),
      
      RMSE_in= sqrt(sum((true_vein_90perc - mean_ve_in)^2)/n()),
      RMSE_pr= sqrt(sum((true_vepr_90perc - mean_ve_pr)^2)/n()),
      
      RRMSE_in= sqrt(sum(((true_vein_90perc - mean_ve_in)/true_vein_90perc)^2)/n()),
      RRMSE_pr= sqrt(sum(((true_vepr_90perc - mean_ve_pr)/true_vepr_90perc)^2)/n()),
      
    ))


##95th percentile
(cop_dat_plot_overall_95perc<-cop_dat_plot %>%
    subset(cop== "95 percentile")%>%
    summarise(
      
      mean_vein=mean(mean_ve_in),
      mean_vepr=mean(mean_ve_pr),
      
      lci_vein=lci(mean_ve_in),
      lci_vepr=lci(mean_ve_pr),
      
      uci_vein=uci(mean_ve_in),
      uci_vepr=uci(mean_ve_pr),
      
      
      cov_in=sum(true_vein_95perc >= lci_ve_in & true_vein_95perc <= uci_ve_in )/n(),
      cov_pr=sum(true_vepr_95perc >= lci_ve_pr & true_vepr_95perc <= uci_ve_pr )/n(),
      
      RMSE_in= sqrt(sum((true_vein_95perc - mean_ve_in)^2)/n()),
      RMSE_pr= sqrt(sum((true_vepr_95perc - mean_ve_pr)^2)/n()),
      
      RRMSE_in= sqrt(sum(((true_vein_95perc - mean_ve_in)/true_vein_95perc)^2)/n()),
      RRMSE_pr= sqrt(sum(((true_vepr_95perc - mean_ve_pr)/true_vepr_95perc)^2)/n()),
      
    )%>%
    as.data.frame(.))

##Add true values to datasets
cop_dat_plot_overall_lloq$true_vein<-true_vein_lloq
cop_dat_plot_overall_lloq$true_vepr<-true_vepr_lloq
cop_dat_plot_overall_95perc$true_vein<-true_vein_95perc
cop_dat_plot_overall_95perc$true_vepr<-true_vepr_95perc


#6. Save parametric and non-parametric results

#write.csv(cop_dat_plot_overall_95perc, "ests_para_2dp_95perc.csv", row.names = FALSE)
#write.csv(cop_dat_plot_overall_zero, "ests_para_2dp_zero.csv", row.names = FALSE)
#write.csv(cop_dat_plot_overall_lloq, "ests_para_2dp_lloq.csv", row.names = FALSE)

write.csv(cop_dat_plot_overall_95perc, "ests_nonpara_2dp_95perc.csv", row.names = FALSE)
write.csv(cop_dat_plot_overall_zero, "ests_nonpara_2dp_zero.csv", row.names = FALSE)
write.csv(cop_dat_plot_overall_lloq, "ests_nonpara_2dp_lloq.csv", row.names = FALSE)

#7. Plot
cop_dat_plot %>%
  subset(cop== "LLOQ")%>%
  ggplot(., aes(x=lci_ve_in))+
  geom_histogram(fill = "darkblue")+
  geom_histogram(fill = "pink", aes(x=uci_ve_in))+
  geom_vline(xintercept = true_vein_lloq, color="red")


cop_dat_plot %>%
  subset(cop== "zero")%>%
  ggplot(., aes(x=lci_ve_in))+
  geom_histogram(fill = "darkblue")+
  geom_histogram(fill = "pink", aes(x=uci_ve_in))+
  geom_vline(xintercept = true_vein_zero, color="red")

cop_dat_plot %>%
  subset(cop== "95 percentile")%>%
  ggplot(., aes(x=lci_ve_in))+
  geom_histogram(fill = "darkblue")+
  geom_histogram(fill = "pink", aes(x=uci_ve_in))+
  geom_vline(xintercept = true_vein_95perc, color="red")

cop_dat_plot %>%
  subset(cop== "LLOQ")%>%
  ggplot(., aes(x=lci_ve_pr))+
  geom_histogram(fill = "darkblue")+
  geom_histogram(fill = "pink", aes(x=uci_ve_pr))+
  geom_vline(xintercept = true_vepr_lloq, color="red")

cop_dat_plot %>%
  subset(cop== "zero")%>%
  ggplot(., aes(x=lci_ve_pr))+
  geom_histogram(fill = "darkblue")+
  geom_histogram(fill = "pink", aes(x=uci_ve_pr))+
  geom_vline(xintercept = true_vepr_zero, color="red")


cop_dat_plot %>%
  subset(cop== "95 percentile")%>%
  ggplot(., aes(x=lci_ve_pr))+
  geom_histogram(fill = "darkblue")+
  geom_histogram(fill = "pink", aes(x=uci_ve_pr))+
  geom_vline(xintercept = true_vepr_95perc, color="red")

