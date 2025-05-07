/// Stan model for joint analysis of asymptomatic and symptomatic cases in COVID-19 vaccine trial data (non-parametric)

//splines code adapted from https://mc-stan.org/users/documentation/case-studies/splines_in_stan.html

functions {
  //spline functions
  vector build_b_spline(real[] t, real[] ext_knots, int ind, int order);
  vector build_b_spline(real[] t, real[] ext_knots, int ind, int order) {
    // INPUTS:
    //    t:          the points at which the b_spline is calculated
    //    ext_knots:  the set of extended knots
    //    ind:        the index of the b_spline
    //    order:      the order of the b-spline
    vector[size(t)] b_spline;
    vector[size(t)] w1 = rep_vector(0, size(t));
    vector[size(t)] w2 = rep_vector(0, size(t));
    if (order==1)
      for (i in 1:size(t)) // B-splines of order 1 are piece-wise constant
        b_spline[i] = (ext_knots[ind] <= t[i]) && (t[i] < ext_knots[ind+1]);
    else {
      if (ext_knots[ind] != ext_knots[ind+order-1])
        w1 = (to_vector(t) - rep_vector(ext_knots[ind], size(t))) /
             (ext_knots[ind+order-1] - ext_knots[ind]);
      if (ext_knots[ind+1] != ext_knots[ind+order])
        w2 = 1 - (to_vector(t) - rep_vector(ext_knots[ind+1], size(t))) /
                 (ext_knots[ind+order] - ext_knots[ind+1]);
      // Calculating the B-spline recursively as linear interpolation of two lower-order splines
      b_spline = w1 .* build_b_spline(t, ext_knots, ind, order-1) +
                 w2 .* build_b_spline(t, ext_knots, ind+1, order-1);
    }
    return b_spline;
  }
}

data {
  int<lower=0> N;
  int<lower=0> C[N];
  int<lower=0> A[N];
  int <lower=0, upper=1> vaccinated[N];
  real Z[N];
  int num_knots;            // num of knots
  vector[num_knots] knots;  // the sequence of knots
  int spline_degree;        // the degree of spline (is equal to order - 1)
  int spline_length;        //number of discrete points in the spline
  real X[spline_length];    //immune CoP values for spline
  real dX;                  //step size of immune CoP
  int Z_index[N];           //value of Z as location in spline (x position)
}

transformed data {
  int num_basis = num_knots + spline_degree - 1; // total number of B-splines
  matrix[num_basis, spline_length] B;  // matrix of B-splines for VE_in
  matrix[num_basis, spline_length] B_pr;  // matrix of B-splines for VE_pr
  vector[spline_degree + num_knots] ext_knots_temp;
  vector[2*spline_degree + num_knots] ext_knots; // set of extended knots
  ext_knots_temp = append_row(rep_vector(knots[1], spline_degree), knots);
  ext_knots = append_row(ext_knots_temp, rep_vector(knots[num_knots], spline_degree));
  for (ind in 1:num_basis)
    B[ind,:] = to_row_vector(build_b_spline(X, to_array_1d(ext_knots), ind, spline_degree + 1));
  B[num_knots + spline_degree - 1, spline_length] = 1;
  for (ind in 1:num_basis)
    B_pr[ind,:] = to_row_vector(build_b_spline(X, to_array_1d(ext_knots), ind, spline_degree + 1));
  B_pr[num_knots + spline_degree - 1, spline_length] = 1;
}

parameters {
  real alpha;
  real gam;
  
  //controlled VE_in spline parameters
  row_vector[num_basis] a_raw;
  real a0;  
  
  //controlled VE_pr spline parameters
  row_vector[num_basis] a_raw_pr;
  real a0_pr;  
}

transformed parameters {
  real <lower=0> mu_c[N];
  real <lower=0> mu_a[N];
  real <lower=0, upper=1> s[N];
  real <lower=0> lam[N];
  real <lower=0, upper=1> VE_in[N];
  real <lower=0, upper=1> VE_pr[N];

  row_vector[num_basis] a; // spline coefficients
  vector[spline_length] cVE;
  row_vector[num_basis] a_pr; // spline coefficients
  vector[spline_length] cVE_pr;
  
  a = a_raw;
  cVE = a0 + to_vector(a*B);
  
  a_pr = a_raw_pr;
  cVE_pr = a0_pr + to_vector(a_pr*B_pr);
  
  for (i in 1:N){
    s[i]=inv_logit(gam);
    lam[i]=exp(alpha);
    VE_in[i]=vaccinated[i]*inv_logit(cVE[Z_index[i]]);
    VE_pr[i]=vaccinated[i]*inv_logit(cVE_pr[Z_index[i]]);
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
  //priors
  a_raw ~ normal(0, 1);
  a0 ~ normal(-5, 2);
 
  a_raw_pr ~ normal(0, 1);
  a0_pr ~ normal(-5, 2);
 
  //likelihood
  target+= poisson_lpmf(C | mu_c);
  target+= poisson_lpmf(A | mu_a);
}

generated quantities {
  vector[spline_length] cVE_hat  = inv_logit(cVE);
  vector[spline_length] cVE_pr_hat  = inv_logit(cVE_pr);
  
  real VE_sym[N];
  real VE_asym[N];
  real ps = inv_logit(gam);
  for (i in 1:N){
  VE_sym[i]= 1-(1-VE_in[i])*(1-VE_pr[i]);
  VE_asym[i]= 1-((1-VE_in[i])-s[i]*(1-VE_in[i])*(1-VE_pr[i]))/(1-s[i]);
  }
}


