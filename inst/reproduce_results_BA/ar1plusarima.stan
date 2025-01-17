data{
  int<lower=1> N;
  int<lower=1> N_obs;
  int<lower=1> N_missing;
  //int<lower=1> Nval;
  int<lower=1, upper=N_obs+N_missing>  which_missing[N_missing];
  int<lower=1, upper=N_obs+N_missing>  which_obs[N_obs];
  int<lower=1> K;         // Number of predictors (columns in the design matrix)

  matrix[N, K] X;         // Design matrix (N x K)
  vector[N] dy;            // Response variable (dy)
  vector[N_obs] y_obs;

  //vector[Nval] Qv;
  //vector[N+1] Qu;
}

parameters{
  vector[K] beta;
  real theta_phi;  // AR(1) coefficient
  real theta_phi_y;  // AR(1) coefficient
  real<lower=0> kappa_dy;
  real<lower=0> kappa_y;

  vector[N_missing] y_missing;
}

transformed parameters{
  real<lower=0> innov_sigma_dy;
  //real<lower=0> innov_sigma_y;
  real<lower=-1, upper=1> phi;
  real<lower=-1, upper=1> phi_y;

  vector[N] X_mu;
  vector[N] y;

  //vector[N] diag_Q;
  //vector[N-1] off_diag_Q;
  //vector[Nval] Qw;
  matrix[N,N] Q;


  X_mu = X * beta;
  phi = 2*exp(theta_phi)/(1+exp(theta_phi))-1;
  phi_y = 2*exp(theta_phi_y)/(1+exp(theta_phi_y))-1;

  innov_sigma_dy = 1/sqrt(kappa_dy)*sqrt(1-phi^2);
  //innov_sigma_y = sigma_y*sqrt(1-phi_y^2);

  y[which_missing] = y_missing;
  y[which_obs] = y_obs;

  // diag_Q[1] = kappa_y;  // First element
  //for (i in 2:N - 1) {
//    diag_Q[i] = (1 + phi_y^2) *kappa_y;
  //}
  //diag_Q[N] = kappa_y;  // Last element

  //for (i in 1:N - 1) {
    //off_diag_Q[i] = -phi_y *kappa_y;
  //}

/*
//noise = sigma^2*(1-rho^2)
Qw[1] = kappa_y/(1-square(phi_y))*1;
Qw[2] = -kappa_y/(1-square(phi_y))*square(1+phi_y);
Qw[2] = kappa_y/(1-square(phi_y))*phi_y;
Qw[4] = -kappa_y/(1-square(phi_y))*square(1+phi_y);
Qw[5] = kappa_y/(1-square(phi_y))*2*(1+phi_y+square(phi_y));
Qw[6] = -kappa_y/(1-square(phi_y))*square(1+phi_y);
Qw[7] = kappa_y/(1-square(phi_y))*phi_y;
for(i in 1:(Nval-14)){
  if( i%5 == 1 || i%5 == 5) Qw[i+7] = kappa_y/(1-square(phi_y))*phi_y;
  if( i%5 == 2 || i%5 == 4) Qw[i+7] = -kappa_y/(1-square(phi_y))*square(1+phi_y);
  if( i%5 == 3) Qw[i+7] = kappa_y/(1-square(phi_y))*2*(1+phi_y+square(phi_y));
}
Qw[Nval-6] = kappa_y/(1-square(phi_y))*phi_y;
Qw[Nval-5] = -kappa_y/(1-square(phi_y))*square(1+phi_y);
Qw[Nval-4] = kappa_y/(1-square(phi_y))*(2+2*phi_y+square(phi_y));
Qw[Nval-3] = -kappa_y/(1-square(phi_y))*(1+phi_y);
Qw[Nval-2] = kappa_y/(1-square(phi_y))*phi_y;
Qw[Nval-1] = -kappa_y/(1-square(phi_y))*(1+phi_y);
Qw[Nval] =kappa_y/(1-square(phi_y));
*/
Q = rep_matrix(0.0, N, N); // Initialize all elements to zero
Q[1,1] = kappa_y/(1-square(phi_y))*(2+2*phi_y+square(phi_y));
Q[1,2] = -kappa_y/(1-square(phi_y))*square(1+phi_y);
Q[1,3] = kappa_y/(1-square(phi_y))*phi_y;
Q[2,1] = -kappa_y/(1-square(phi_y))*square(1+phi_y);
Q[2,3] = -kappa_y/(1-square(phi_y))*square(1+phi_y);
Q[2,4] = kappa_y/(1-square(phi_y))*phi_y;
for(i in 2:(N-2)){
  Q[i,i] = kappa_y/(1-square(phi_y))*2*(1+phi_y+square(phi_y));
}
for(i in 2:(N-1)){
  Q[i,i-1] = -kappa_y/(1-square(phi_y))*square(1+phi_y);
  Q[i-1,i] = -kappa_y/(1-square(phi_y))*square(1+phi_y);
}
for(i in 3:N){
  Q[i,i-2] = kappa_y/(1-square(phi_y))*phi_y;
  Q[i-2,i] = kappa_y/(1-square(phi_y))*phi_y;
}
Q[N-1,N-1] = kappa_y/(1-square(phi_y))*(2+2*phi_y+square(phi_y));
Q[N-1,N] = -kappa_y/(1-square(phi_y))*(1+phi_y);
Q[N,N-1] = -kappa_y/(1-square(phi_y))*(1+phi_y);
Q[N,N] = kappa_y/(1-square(phi_y));

}

model {

  //vector[N] wx;
  vector[N] cumsum_dy;
  //matrix[N,N] Q;


  1/kappa_dy ~ gamma(1,5e-05);
  1/kappa_y ~ gamma(1,5e-05);

  beta ~ normal(0, sqrt(100));         // Prior for regression coefficients
  theta_phi ~ normal(0, 1);  // Prior for AR(1) coefficient
  theta_phi_y ~ normal(0, 1);  // Prior for AR(1) coefficient

  // Likelihood for layer increments: dy
  dy[1] ~ normal(X_mu[1], 1/sqrt(kappa_dy));  // First observation
  dy[2:N] ~ normal(X_mu[2:N] + phi * (dy[1:(N-1)] - X_mu[1:(N-1)]), innov_sigma_dy  );


  cumsum_dy = cumulative_sum(dy);
  //wx = csr_matrix_times_vector(N, N, Rvalues, col_indices, row_pointers,  y);
  //Q = csr_matrix_t(N, N, Qw, Qv, Qu); // might need to switch Qv and Qu


  //y ~ normal(cumsum_dy - (wx ./ sds), sigma_y/sqrt(sds));
  y ~ multi_normal_prec(cumsum_dy, Q);

}

