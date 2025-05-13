data {
  int<lower=0> N;
  int<lower=0> C[N];
  int<lower=0> A[N];
  int <lower=0, upper=1> vaccinated[N];
  real Z[N];
  real<lower=0> pers_yrs_at_risk[N]; 
  real min_Z;
  real max_Z;

  real cov1[N]; //HCW 0 covid patients
  real cov2[N]; //HCW 1+ covid patients
  real cov3[N]; // age (years)
  real cov4[N]; // obese (y/n)
  real cov5[N]; // white (y/n)
  real cov6[N]; // comorbidity (y/n)
  real <lower=0> weight[N];
}

parameters {
  real alpha_ref;
  real alpha1;
  real alpha2;
  real alpha3;
  real alpha4;
  real alpha5;
  real alpha6;
  real gam;
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
    lam[i]=exp(alpha_ref + alpha1*cov1[i] + alpha2*cov2[i] + alpha3*cov3[i] + alpha4*cov4[i] + alpha5*cov5[i] + alpha6*cov6[i] )*pers_yrs_at_risk[i];
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
  for(i in 1:N){ 
    target += poisson_lpmf(C[i] | mu_c[i]) * weight[i];
    target += poisson_lpmf(A[i] | mu_a[i]) * weight[i];
  }
    }

generated quantities{
  real VE_sym[N];
  real VE_asym[N];
  real VE_sym_mean;
  real VE_asym_mean;
  real ps = inv_logit(gam);
  real foi_ref = exp(alpha_ref);
  real foi_cov1 = exp(alpha_ref + alpha1);
  real foi_cov2 = exp(alpha_ref + alpha2);
  real foi_cov3 = exp(alpha_ref + alpha3);
  real foi_cov4 = exp(alpha_ref + alpha4);
  real foi_cov5 = exp(alpha_ref + alpha5);
  real foi_cov6 = exp(alpha_ref + alpha6);
  vector[N] log_lik;
  real log_lik_sum;
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
  log_lik[i] = (poisson_lpmf(C[i] | mu_c[i]) + poisson_lpmf(A[i] | mu_a[i]))*weight[i];
  }

log_lik_sum = sum(log_lik);
  
}
