data {
  int<lower = 0> N;
  vector[N] lower1;
  vector[N] upper1;
  vector[N] lower2;
  vector[N] upper2;
}

parameters {
  real<lower=0, upper=1> p1;
  real<lower=0, upper=1> p2;
  // Upgrade about mod2: impose contingency tab on the c vector,
  // s.t. c_{1,j} and c_{2,j} are not independent.
  real<lower=-1, upper=1> deltaOrig; 
  ordered[2] Beta1;
  ordered[2] Beta2;
  vector<lower=0>[2] Sigma;
  corr_matrix[2] R;
  vector<lower=0, upper=1>[N] latent_raw1;
  vector<lower=0, upper=1>[N] latent_raw2;
}

transformed parameters {
  vector[N] latent1 = lower1 + (upper1-lower1) .* latent_raw1;
  vector[N] latent2 = lower2 + (upper2-lower2) .* latent_raw2;
  matrix[N, 2] latenty = append_col(latent1, latent2); 
  vector[2] Beta_nn = [Beta1[1], Beta2[1]]';
  vector[2] Beta_rn = [Beta1[2], Beta2[1]]';
  vector[2] Beta_nr = [Beta1[1], Beta2[2]]';
  vector[2] Beta_rr = [Beta1[2], Beta2[2]]';
  // Upgrade about mod3: scale deltaOrig to [-min(p1p2, (1-p1)(1-p2)), min(p1, p2)-p1p2]
  real<lower=-1, upper=0> deltaLower = -fmin(p1*p2, (1-p1)*(1-p2));
  real<lower=0, upper=1> deltaUpper = fmin(p1,p2)-p1*p2;
  real delta = (deltaUpper-deltaLower)/2*deltaOrig + (deltaUpper+deltaLower)/2;
  real<lower=0, upper=1> p_nn = (1-p1) * (1-p2) + delta;
  real<lower=0, upper=1> p_rn = p1 * (1-p2) - delta;
  real<lower=0, upper=1> p_nr = (1-p1) * p2 - delta;
  real<lower=0, upper=1> p_rr = p1 * p2 + delta;
}

model {
  real ps[4]; 
  p1 ~ beta(0.5, 0.5);
  p2 ~ beta(0.5, 0.5);
  Beta1 ~ normal(0, 100);
  Beta2 ~ normal(0, 100);
  Sigma ~ cauchy(0, 2);
  R ~ lkj_corr(1);
  
  for (n in 1:N) {
    ps[1] = log(p_nn) + multi_normal_lpdf(latenty[n] | Beta_nn, quad_form_diag(R, Sigma));
    ps[2] = log(p_rn) + multi_normal_lpdf(latenty[n] | Beta_rn, quad_form_diag(R, Sigma));
    ps[3] = log(p_nr) + multi_normal_lpdf(latenty[n] | Beta_nr, quad_form_diag(R, Sigma));
    ps[4] = log(p_rr) + multi_normal_lpdf(latenty[n] | Beta_rr, quad_form_diag(R, Sigma));
    target += log_sum_exp(ps);
  }
}

generated quantities {
  real<lower=-1, upper=1> cor_bin = delta / (sqrt(p1 * (1 - p1)) * sqrt(p2 * (1 - p2)));
}



