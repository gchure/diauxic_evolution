data {
    int<lower=3> N; // Number of timepoints
    vector<lower=0.04, upper=0.4>[N] absorbance; // Absorbance measurements 
    vector<lower=0>[N] time; // Time 
}

parameters {
    real<lower=0> lambda; // Growth rate in units of inverse time
    real<lower=0> absorbance_0; // Absorbance intercept
    real<lower=0> sigma; // Homoscedastic error
}

model {
    // Prior definition
    lambda ~ normal(0, 1);
    absorbance_0 ~ normal(0, 0.1);
    sigma ~ normal(0, 0.1);

    // Likelihood definition
    absorbance ~ normal(absorbance_0 * exp(lambda * time), sigma);
}