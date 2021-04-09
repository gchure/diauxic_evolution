data {
    // Dimensional parameters
    int<lower=1> J; // Number of biological replicates
    int<lower=1> N; // Total number of measurements
    int<lower=1, upper=J> idx[N]; // ID vector

    // Observed data
    vector<lower=0>[N] OD; // optical density measurements
    vector<lower=0>[N] time; // Elapsed time
}

parameters {
    // Level-1 parameters
    vector<lower=0>[J] OD_init;
    vector<lower=0>[J] lam;

    // Hyperparameters
    real<lower=0> OD_init_mu;
    real<lower=0> OD_init_sigma;
    real<lower=0> lam_mu;
    real<lower=0> lam_sigma; 

    // Homoscedastic error
    real<lower=0> sigma;
}

model {
    // Hyperpriors
    sigma ~ std_normal();
    lam_mu ~ std_normal();
    lam_sigma ~ std_normal();
    OD_init_mu ~ std_normal();
    OD_init_sigma ~ std_normal();

    // Level-1 priors
    lam ~ normal(lam_mu, lam_sigma);
    OD_init ~ normal(OD_init_mu, OD_init_sigma);

    // Gaussian likelihood
    OD ~ normal(OD_init[idx] .* exp(time .* lam[idx]), sigma);
}

generated quantities {
    real<lower=0> od_ppc[N];
    for (i in 1:N) {
        od_ppc[i] = normal_rng(OD_init[idx[i]] * exp(time[i] * lam[idx[i]]), sigma);
    }
}