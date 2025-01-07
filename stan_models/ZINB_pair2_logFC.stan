// ZINB model to extract gene-fitness scores from TraDIS infection screens
// uses pre-determined stochastic loss and normalization factors
// fitted to neutral genes (e.g. pseudogenes) (using ZINB_pair2_norm.stan)
functions {
  // number of zero entries in array y
  int num_zeros(array[] int y) {
    int value = 0;
    for (n in 1:size(y)) value += (y[n] == 0);
    return value;
  }
  
  // number of negative entries in array y
  int num_neg(array[] int y) {
    int value = 0;
    for (n in 1:size(y)) value += (y[n] < 0);
    return value;
  }
  
  // returns only input-output pairs where the input
  // transposon insertion site (TIS) is non-zero
  // arguments:
  // it_g: gene index of TIS in y_in and y_out
  // y_in: input counts of TIS
  // y_out: output counts of TIS
  // N_nonzero: number of non-zero input counts, determined by
  // running num_zeros on y_in before.
  array[] int y_nonzero(array[] int it_g, array[] int y_in, array[] int y_out, int N_nonzero) {
    // return value has 3 times the length of the number of non-zero entries in y_in
    // 1:N_nonzero : gene index
    // (N_nonzero+1):(2*N_nonzero) : non-zero input count
    // (2*N_nonzero+1):(3*N_nonzero) : output count, may be zero
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
  
  // zero-inflated negative binomial model which is run in parallel for every pair of
  // input-output sequencing libraries
  vector zinb_model(vector param, vector theta, array[] real X, array[] int y){
    
    int N_ti = y[1]; // number of non-zero input TIs
    int N_g = y[2]; // number of genes
    array[N_ti] int gene = y[3:(N_ti+2)]; // gene index
    array[N_ti] int Y_out = y[(N_ti+3):(2*N_ti+2)]; // output read counts per TIS
    
    vector[N_ti] alpha_ti = to_vector(X[1:N_ti]); // normalized input log-counts
    real alpha_out_mean = X[1+N_ti];// maximum difference from library-wise input log-count
    real alpha_max = X[2+N_ti]; // library-wise maximum of normalized input log-counts
    real b = X[4+N_ti]; // abundance correction of dispersion
    array[2] real a = X[(5+N_ti):(6+N_ti)];
    
    vector[N_g] logFC_g = param[1:N_g]; // geme-wise log-fold changes between input and output condition, fitness scores
    vector[N_g] inv_rho_g = param[(N_g+1):(2*N_g)]; // gene-wise coefficient of variation (CV) for TIS with log-count alpha_max
    real gamma_ps = X[3+N_ti]; // logit of sample-wise stochastic loss
    real ll=0; // log-likelihood
    
    for (it_ti in 1:N_ti) {
      int it_g = gene[it_ti];
      // abundance dependent dispersion
      real rho_ti = inv(inv_rho_g[it_g]*(1+b*alpha_ti[it_ti]/alpha_max)).*inv(inv_rho_g[it_g]*(1+b*alpha_ti[it_ti]/alpha_max));
      // logit of expected stochastic loss for a given insertion site (depends on the expected count)
      real gamma_ti = gamma_ps + a[1] + a[2]*alpha_ti[it_ti];
      // contribution of zero output counts to log-likelihood
      if (Y_out[it_ti]==0) {
        ll += log_sum_exp(bernoulli_logit_lpmf(1 | gamma_ti),
                          bernoulli_logit_lpmf(0 | gamma_ti) + neg_binomial_2_log_lpmf(0 | alpha_out_mean + alpha_ti[it_ti] + logFC_g[it_g], rho_ti)
              );
      // contribution on non-zero output counts
      } else {
        ll += bernoulli_logit_lpmf(0 | gamma_ti);
        ll += neg_binomial_2_log_lpmf(Y_out[it_ti] | alpha_out_mean + alpha_ti[it_ti] + logFC_g[it_g], rho_ti);
      }
    }
    return [ll]';
  }
}

data{
  int<lower=1> N_ti; // total number of transposon insertion sites (TIS)
  int<lower=1> N_ti_max; // max. number of non-zero TIS per input sample
  int<lower=1> N_g; // number of genes
  int<lower=1> N_exp; // number of input-output pairs of sequencing libraries
  array [N_exp] real norm_fac; // normalization factors obtained from running NB_pair_pseudo
  vector [N_exp] Z_ps; // fraction of pseudogenes with zero count in output (library-wise)
  array [N_ti] int<lower=1> it_g; // gene index
  array [N_exp, N_ti] int raw_counts_in; // raw count matrix input
  array [N_exp, N_ti] int raw_counts_out; // raw count matrix output, assumed to be same order as input
  real b; // abundance dependence of dispersion, fitted during normalization
  array[2] real a_in; // abundance dependence of stochastic loss, fitted during normalization
}

transformed data{
  // we bring the data into a library-wise (input-output pair) structure for parallelization
  array [N_exp, 6+N_ti_max] real X; // non-integer library-wise variables
  array [N_exp, 2+2*N_ti_max] int y_out; // mostly library-wise output count matrix
  vector [N_exp] alpha_out_mean; // library-wise mean of input log-counts
  vector [N_exp] alpha_max; // maximum difference from library-wise input log-count
  
  // format data for every input-output pair of sequencing libraries
  for (it_exp in 1:N_exp) {
    array [N_ti] int counts_in = raw_counts_in[it_exp,];
    array [N_ti] int counts_out = raw_counts_out[it_exp,];
    int N_g0_in = num_zeros(counts_in); // number of zeros in input library
    int N_g_in = N_ti - N_g0_in; // number of non-zero TIS in input library
    //remove TIS with zero count in input (if there are any) from input AND output library
    array [3*N_g_in] int nonzero_counts = y_nonzero(it_g, counts_in, counts_out, N_g_in);
    vector [N_g_in] alpha_ti = log(to_vector(nonzero_counts[(N_g_in+1):(2*N_g_in)])) + norm_fac[it_exp];  //normalized log-count of input pair
    alpha_out_mean[it_exp] = mean(alpha_ti);
    alpha_max[it_exp] = max(abs(alpha_ti - alpha_out_mean[it_exp]));
    
    y_out[it_exp, 1] = N_g_in; // number of non-zero TIS in input library
    y_out[it_exp, 2] = N_g; // total number of genes
    y_out[it_exp, 3:(N_g_in+2)] = nonzero_counts[1:N_g_in]; // gene indices for output counts
    y_out[it_exp, (N_g_in+3):(2*N_g_in+2)] = nonzero_counts[(2*N_g_in+1):(3*N_g_in)]; // output count matrix for non-zero input counts
    
    // subtract library-wise mean from log input count
    X[it_exp, 1:N_g_in] = to_array_1d(alpha_ti - alpha_out_mean[it_exp]);
    X[it_exp, 1+N_g_in] = alpha_out_mean[it_exp];
    X[it_exp, 2+N_g_in] = alpha_max[it_exp];
    X[it_exp, 3+N_g_in] = logit(Z_ps[it_exp]);
    X[it_exp, 4+N_g_in] = b;
    X[it_exp, (5+N_g_in):(6+N_g_in)] = a_in;
  }
}

parameters{
  // hyperparameters of hierarchical priors
  // Coefficient of variation (CV)
  real<lower=0> sigma_rho;
  real<lower=0.2> mu_rho;
  // mean log-fold change
  real mu;
  
  // gene-wise log-fold change
  vector[N_g] logFC_g;
  // maximal gene-wise coefficient of variation (CV) of NB distribution
  vector<lower=0>[N_g] inv_rho_g; // actually 1/sqrt(rho_g)
}

transformed parameters{
}

model{
  // priors
  sigma_rho ~ cauchy(0,1);
  mu_rho ~ cauchy(1,0.5);
  mu ~ cauchy(0,1);
  
  // hierarchical priors
  // we keep the width of the prior on the log-fold changes fixed to avoid shrinkage of the
  // effect size
  logFC_g ~ normal(mu, 5);
  inv_rho_g ~ normal(mu_rho, sigma_rho);
  
  {
    vector[2*N_g] param; // global model parameters
    array[N_exp] vector[0] theta; // no library-wise model parameters
    
    param[1:N_g] = logFC_g;
    param[(N_g+1):(2*N_g)] = inv_rho_g;
    // add log-likelihood for every input-output pair of sequencing libraries
    target += sum(map_rect(zinb_model, param, theta, X, y_out));
  }
}

generated quantities{
  // library-wise stochastic loss between input and output
  vector[N_exp] gamma_sa = logit(Z_ps);
  // library-wise mean input log-count
  vector[N_exp] alpha_sa_mean = alpha_out_mean;
  // maximum difference from library-wise input log-count
  vector[N_exp] alpha_sa_max = alpha_max;
  // abundance dependence of stochastic dropout fitted during normalization
  vector[2] a = to_vector(a_in);
}
