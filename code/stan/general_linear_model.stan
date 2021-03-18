data {
    int<lower=1> N; // Number of measurement
    vector[N] x;
    vector[N] y;
}

transformed data {
    vector[N] centered_x = (x - mean(x)) / sd(x);
    vector[N] centered_y = (y - mean(y)) / sd(y);
}

parameters {
    real slope_std;
    real intercept_std;
    real<lower=0> sigma_std;
}

model {
    slope_std ~ std_normal();
    intercept_std ~ std_normal();
    sigma_std ~ std_normal();
    centered_y ~ normal(centered_x * slope_std + intercept_std, sigma_std);
} 

generated quantities {
    real slope = slope_std * sd(y)/sd(x);
    real intercept = sd(y) * (intercept_std - slope_std * mean(x) / sd(x)) + mean(y);
    real sigma = sd(y)  * sigma_std;
}
