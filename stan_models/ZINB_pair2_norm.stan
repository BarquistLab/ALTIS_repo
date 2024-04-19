// ZINB model for the analysis of TraDIS data
// uses pre-determined dropout rates fitted to neutral/pseudo genes
// center logFC around zero, can be set to zero
// models abundance dependent dispersion
functions {
  int num_zeros(array[] int y) {
    int value = 0;
    for (n in 1:size(y)) value += (y[n] == 0);
    return value;
  }
  
  int num_neg(array[] int y) {
    int value = 0;
    for (n in 1:size(y)) value += (y[n] < 0);
    return value;
  }
  
  array[] int y_nonzero(array[] int it_g, array[] int y_in, array[] int y_out, int N_nonzero) {
    array[3*N_nonzero] int nonzero;
    int counter = 1;
    for (n in 1:size(y_in)) {
      if (y_in[n] > 0) {
        nonzero[counter] = it_g[n]; // gene index
        nonzero[N_nonzero+counter] = y_in[n];
        nonzero[2*N_nonzero+counter] = y_out[n];
        counter += 1;
      }
    }
    return nonzero;
  }
  
  vector zinb_model(vector param, vector theta, array[] real X, array[] int y){
    
    int N_ti = y[1]; // number of non-zero input TIs
    int N_g = y[2];
    array[N_ti] int gene = y[3:(N_ti+2)]; // gene index
    array[N_ti] int Y_out = y[(N_ti+3):(2*N_ti+2)]; // read counts out
    
    vector[N_ti] alpha_ti = to_vector(X[1:N_ti]); // read counts in
    real alpha_out_mean = X[1+N_ti];
    real alpha_max = X[2+N_ti];
    
    vector[N_g] inv_rho_g = param[1:N_g];
    vector[2] a = param[(N_g+2):(N_g+3)];
    real gamma_ps = X[3+N_ti];
    real logN_rel = X[4+N_ti];
    real b = param[N_g+1];
    real norm_fac = theta[1];
    real ll=0;
    
    for (it_ti in 1:N_ti) {
      int it_g = gene[it_ti];
      // expected dropout rate for a given insertion site (depends on the expected count)
      real gamma_ti = gamma_ps + a[1] + a[2]*alpha_ti[it_ti];
      // abundance dependent dispersion
      real rho_ti = inv(inv_rho_g[it_g]*(1+b*alpha_ti[it_ti]/alpha_max)).*inv(inv_rho_g[it_g]*(1+b*alpha_ti[it_ti]/alpha_max));
      if (Y_out[it_ti]==0) {
        ll += log_sum_exp(bernoulli_logit_lpmf(1 | gamma_ti),
                          bernoulli_logit_lpmf(0 | gamma_ti) + neg_binomial_2_log_lpmf(0 | alpha_out_mean + alpha_ti[it_ti] + logN_rel + norm_fac, rho_ti)
              );
      } else {
        ll += bernoulli_logit_lpmf(0 | gamma_ti);
        ll += neg_binomial_2_log_lpmf(Y_out[it_ti] | alpha_out_mean + alpha_ti[it_ti] + logN_rel + norm_fac, rho_ti);
      }
    }
    return [ll]';
  }
}

data{
  int<lower=1> N_ti; // total number of transposon insertions
  int<lower=1> N_ti_max; // max. number of non-zero TIs per input sample
  int<lower=1> N_g; // number of genes
  int<lower=1> N_exp; // number of input-output pairs
  vector [N_exp] Z_ps; // fraction of pseudo genes with zero count in output (sample-wise)
  array [N_ti] int<lower=1> it_g; // gene index
  array [N_exp, N_ti] int raw_counts_in; // raw count matrix input
  array [N_exp, N_ti] int raw_counts_out; // raw count matrix output, assumed to be same order as input
}

transformed data{
  array [N_exp, 4+N_ti_max] real X;
  array [N_exp, 2+2*N_ti_max] int y_out;
  vector [N_exp] alpha_out_mean;
  vector [N_exp] alpha_max;
  vector [N_exp] logN_rel;
  
  for (it_exp in 1:N_exp) {
    array [N_ti] int counts_in = raw_counts_in[it_exp,];
    array [N_ti] int counts_out = raw_counts_out[it_exp,];
    logN_rel[it_exp] = log(sum(to_vector(counts_out))/sum(to_vector(counts_in)));
    int N_g0_in = num_zeros(counts_in);
    int N_g_in = N_ti - N_g0_in;
    //remove TIs with zero count in input (if there are any)
    array [3*N_g_in] int nonzero_counts = y_nonzero(it_g, counts_in, counts_out, N_g_in);
    vector [N_g_in] alpha_ti = log(to_vector(nonzero_counts[(N_g_in+1):(2*N_g_in)]));
    alpha_out_mean[it_exp] = mean(alpha_ti);
    alpha_max[it_exp] = max(abs(alpha_ti - alpha_out_mean[it_exp]));
    
    y_out[it_exp, 1] = N_g_in;
    y_out[it_exp, 2] = N_g;
    y_out[it_exp, 3:(N_g_in+2)] = nonzero_counts[1:N_g_in];
    y_out[it_exp, (N_g_in+3):(2*N_g_in+2)] = nonzero_counts[(2*N_g_in+1):(3*N_g_in)];
    
    X[it_exp, 1:N_g_in] = to_array_1d(alpha_ti - alpha_out_mean[it_exp]);
    X[it_exp, 1+N_g_in] = alpha_out_mean[it_exp];
    X[it_exp, 2+N_g_in] = alpha_max[it_exp];
    X[it_exp, 3+N_g_in] = logit(Z_ps[it_exp]);
    X[it_exp, 4+N_g_in] = logN_rel[it_exp];
  }
}

parameters{
  real<lower=0> sigma_rho;
  real<lower=0.2> mu_rho;
  real mu;
  real<lower=0> sigma;
  
  vector<lower=0>[N_g] inv_rho_g; // actually 1/sqrt(rho_g)
  
  real<lower=-1, upper=1> b;
  vector[2] a;
  array[N_exp] real norm_fac_std;
}

transformed parameters{
}

model{
  // priors
  sigma_rho ~ cauchy(0,1);
  mu_rho ~ cauchy(1,0.5);

  sigma ~ cauchy(0,1);
  mu ~ cauchy(0,1);
  
  b ~ std_normal();
  a ~ std_normal();
  
  norm_fac_std ~ normal(mu, sigma);
  inv_rho_g ~ normal(mu_rho, sigma_rho);
  
  {
    vector[N_g+3] param;
    array[N_exp] vector[1] theta;
    
    param[1:N_g] = inv_rho_g;
    param[N_g+1] = b;
    param[(N_g+2):(N_g+3)] = a;
    theta[, 1] = norm_fac_std;
    target += sum(map_rect(zinb_model, param, theta, X, y_out));
  }
}

generated quantities{
  vector[N_exp] gamma_sa = logit(Z_ps);
  vector[N_exp] alpha_sa_mean = alpha_out_mean;
  vector[N_exp] alpha_sa_max = alpha_max;
  vector[N_exp] norm_fac = to_vector(norm_fac_std) + logN_rel;
}
