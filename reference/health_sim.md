# Simulated health data

Simulated health data associated with ln(NO2) concentration on Jan 10,
2012. For details, see `health_sim.R`.

## Usage

``` r
data(health_sim)
```

## Format

A data frame with n = 2000 rows and 6 variables:

- Y:

  simulated continuous health outcome

- Ybinary:

  simulated binary health outcome

- lon:

  simulated health subject longitude

- lat:

  simulated health subject latitude

- Z:

  simulated covariate (p=1) that is not subject to measurement error

- X_true:

  true ln(NO2) exposure used for simulating health outcome
