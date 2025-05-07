// Stan model for joint analysis of asymptomatic and symptomatic cases in COVID-19 vaccine trial data (parametric)

data {
  int<lower=0> N;
  int<lower=0> C[N];
  int<lower=0> A[N];
  int <lower=0, upper=1> vaccinated[N];
  real Z[N];
  //real<lower=0> pers_yrs_at_risk[N]; //removed for the simulations, added for COV002 data
  real min_Z;
  real max_Z;
}

parameters {
  real alpha;
  real gam;
 // real <lower=-1, upper=1> VE_direct_in;
  //real <lower=0, upper=1> VE_direct_pr; //assuming >0
  real <lower=0, upper=1> p_direct_in;
  real <lower=0, upper=1> p_direct_pr;
  real <lower=min_Z, upper=max_Z> n_50_in;
  real <lower=min_Z, upper=max_Z> n_50_pr;
  real <lower=0, upper =4> log_kk_in; 
  real <lower=0, upper =4> log_kk_pr;
  real <lower=0, upper=1> scale_in;
  real <lower=0, upper=1> scale_pr;
}

transformed parameters {
  real <lower=0> mu_c[N];
  real <lower=0> mu_a[N];
  real <lower=0, upper=1> s[N];
  real <lower=0> lam[N];
  real <lower=-1, upper=1> VE_in[N];
  real <lower=0, upper=1> VE_pr[N];  
  real kk_in;
  real kk_pr;
  real <lower=0, upper=1> VE_indirect_in[N];
  real <lower=0, upper=1> VE_indirect_pr[N]; 
 
  kk_in=exp(log_kk_in)-1;
  kk_pr=exp(log_kk_pr)-1;
  
  for (i in 1:N){
    s[i]=inv_logit(gam);
    lam[i]=exp(alpha);//*pers_yrs_at_risk[i];
    
    VE_indirect_in[i]=vaccinated[i]*((1-p_direct_in)*scale_in)*(inv_logit(kk_in*(Z[i]-n_50_in))); 
    VE_indirect_pr[i]=vaccinated[i]*((1-p_direct_pr)*scale_pr)*(inv_logit(kk_pr*(Z[i]-n_50_pr)));
   
    VE_in[i]=vaccinated[i]*(p_direct_in*scale_in + VE_indirect_in[i]);
    VE_pr[i]=vaccinated[i]*(p_direct_pr*scale_pr + VE_indirect_pr[i]);
   
    
    if(vaccinated[i]==1){
      mu_c[i] =(1-VE_in[i])*(1-VE_pr[i])*lam[i]*s[i];
      mu_a[i] =(1-VE_in[i])*lam[i]*((1-s[i])+s[i]*VE_pr[i]);
    }else if(vaccinated[i]==0){
      mu_c[i] =lam[i]*s[i];
      mu_a[i] =lam[i]*(1-s[i]);
    }
  }
}

model {
  //likelihood
  target+= poisson_lpmf(C | mu_c);
  target+= poisson_lpmf(A | mu_a);
}

generated quantities{
  real VE_sym[N];
  real VE_asym[N];
  real VE_sym_mean;
  real VE_asym_mean;
  real ps = inv_logit(gam);
  real foi = exp(alpha);
  vector[N] log_lik;
  real VE_in_mean = sum(VE_in)/sum(vaccinated);
  real VE_pr_mean = sum(VE_pr)/sum(vaccinated);
  real <lower=0, upper=1> VE_direct_in;
  real <lower=0, upper=1> VE_direct_pr; 
  real <lower=0, upper=1> VE_indirect_in_max;
  real <lower=0, upper=1> VE_indirect_pr_max; 
  real <lower=0, upper=1> VE_indirect_in_mean;
  real <lower=0, upper=1> VE_indirect_pr_mean; 
  
  VE_direct_in=scale_in*p_direct_in;
  VE_direct_pr=scale_pr*p_direct_pr;
  VE_indirect_in_max=scale_in*(1-p_direct_in);
  VE_indirect_pr_max=scale_pr*(1-p_direct_pr);
    
  for (i in 1:N){
  VE_sym[i]= 1-(1-VE_in[i])*(1-VE_pr[i]);
  VE_asym[i]= 1-((1-VE_in[i])-s[i]*(1-VE_in[i])*(1-VE_pr[i]))/(1-s[i]);

  }
  VE_sym_mean = sum(VE_sym)/sum(vaccinated);
  VE_asym_mean = sum(VE_asym)/sum(vaccinated);
  VE_indirect_in_mean = sum(VE_indirect_in)/sum(vaccinated);
  VE_indirect_pr_mean = sum(VE_indirect_pr)/sum(vaccinated);
  
  // LOO cross validation: define the log likelihood as a vector
  for (i in 1:N){
  log_lik[i] = poisson_lpmf(C[i] | mu_c[i]) + poisson_lpmf(A[i] | mu_a[i]);
  }
  
}
