# Bayesian generalized linear models with spatial exposure measurement error.

This function fits a Bayesian generalized linear model in the presence
of spatial exposure measurement error for covariate(s) \\X\\. One of the
most important features of this function is that it allows a sparse
matrix input for the prior precision matrix of \\X\\ for scalable
computation. As of version 1.0.0, only the Bayesian logistic regression
model is supported among GLMs, and function `bglm_me()` runs a Gibbs
sampler to carry out posterior inference using Polya-Gamma augmentation
(Polson et al., 2013). See the "Details" section below for the model
description and Lee et al. (2024) for an application example in
environmental epidemiology.

## Usage

``` r
bglm_me(
  Y,
  X_mean,
  X_prec,
  Z,
  family = binomial(link = "logit"),
  nburn = 5000,
  nsave = 5000,
  nthin = 1,
  prior = NULL,
  saveX = FALSE
)
```

## Arguments

- Y:

  *vector\<int\>*, n by 1 binary response vector.

- X_mean:

  *vector\<num\>*, n by 1 prior mean vector \\\mu_X\\. When there are q
  multiple exposures subject to measurement error, it can be a length q
  list of n by 1 vectors.

- X_prec:

  *matrix\<num\>*, n by n prior precision matrix \\Q_X\\, which allows
  sparse format from
  [Matrix](https://rdrr.io/pkg/Matrix/man/Matrix.html) package. When
  there are q multiple exposures subject to measurement error, it can be
  a length q list of n by n matrices.

- Z:

  *matrix\<num\>*, n by p matrix containing p covariates that are not
  subject to measurement error.

- family:

  *class [family](https://rdrr.io/r/stats/family.html)*, a description
  of the error distribution and the link function to be used in the
  model. Currently, it only supports binomial(link = "logit").

- nburn:

  *integer*, number of burn-in iterations (default=5000).

- nsave:

  *integer*, number of posterior samples (default=5000). Total number of
  MCMC iteration is `nburn + nsave * nthin`.

- nthin:

  *integer*, thin-in rate (default=1).

- prior:

  *list*, list of prior parameters of the regression model. Default is
  `list(var_beta = 100)`.

- saveX:

  *logical*, default FALSE, whether save posterior samples of X
  (exposure).

## Value

List of the following:

- posterior:

  `nsave` by (q+p) matrix of posterior samples of \\\beta_x\\ and
  \\\beta_z\\ as a coda::[mcmc](https://rdrr.io/pkg/coda/man/mcmc.html)
  object.

- time:

  time taken for running MCMC (in seconds)

- X_save:

  (if `saveX = TRUE`) posterior samples of X

## Details

Let \\Y_i\\ be a binary response, \\X_i\\ be a \\q\times 1\\ covariate
vector that is subject to spatial exposure measurement error, and
\\Z_i\\ be a \\p\times 1\\ covariate vector without measurement error.
Consider a logistic regression model, independently for each
\\i=1,\dots,n,\\ \$\$\log(Pr(Y_i = 1)/Pr(Y_i = 0)) = \beta_0 + X_i^\top
\beta_X + Z_i^\top \beta_Z.\$\$ Spatial exposure measurement error of
\\X_i\\ (for \\i=1,\dots,n\\) is incorporated into the model using a
multivariate normal prior. For example, when \\q=1\\, we have an
\\n-\\dimensional multivariate normal prior on \\X =
(X_1,\dots,X_n)^\top\\, \$\$(X_1,\dots,X_n)\sim N_n(\mu_X,
Q_X^{-1}).\$\$ Most importantly, it allows a sparse matrix input for the
prior precision matrix \\Q_X\\ for scalable computation, which can be
obtained by Vecchia approximation. When \\q\>1\\, \\q\\ independent
\\n-\\dimensional multivariate normal priors are assumed.

We consider normal priors for regression coefficients, \$\$\beta_0 \sim
N(0, V\_\beta), \quad \beta\_{X,j} \stackrel{iid}{\sim} N(0, V\_\beta),
\quad \beta\_{Z,k} \stackrel{iid}{\sim} N(0, V\_\beta)\$\$ where
`var_beta` corresponds to \\V\_\beta\\.

## References

Polson, N. G., Scott, J. G., & Windle, J. (2013). Bayesian inference for
logistic models using Pólya–Gamma latent variables. Journal of the
American statistical Association, 108(504), 1339-1349.

Lee, C. J., Symanski, E., Rammah, A., Kang, D. H., Hopke, P. K., & Park,
E. S. (2024). A scalable two-stage Bayesian approach accounting for
exposure measurement error in environmental epidemiology. arXiv preprint
arXiv:2401.00634.

## Examples

``` r
if (FALSE) { # \dontrun{
library(bspme)
data(NO2_Jan2012)
data(health_sim)
library(fields)
library(maps)
# Obtain the predicted exposure mean and covariance at simulated health subject locations
# based on NO2 data obtained on Jan 10, 2012
# using a Gaussian process prior with mean zero and exponential covariance kernel
# with a fixed range 8 (in km) and standard deviation 1.

# exposure data
data_jan10 = NO2_Jan2012[NO2_Jan2012$date == as.POSIXct("2012-01-10"),]
coords_monitor = cbind(data_jan10$lon, data_jan10$lat)

# health data
coords_health = cbind(health_sim$lon, health_sim$lat)

distmat_xx <- rdist.earth(coords_monitor, miles = F)
distmat_xy <- rdist.earth(coords_monitor, coords_health, miles = F)
distmat_yy <- rdist.earth(coords_health, miles = F)

a = 8; sigma = 1; # assume known

Sigmaxx = fields::Matern(distmat_xx, smoothness = 0.5, range = a, phi = sigma^2)
Sigmaxy = fields::Matern(distmat_xy, smoothness = 0.5, range = a, phi = sigma^2)
Sigmayy = fields::Matern(distmat_yy, smoothness = 0.5, range = a, phi = sigma^2)

# posterior predictive mean and covariance of exposure at health subject locations
X_mean <- t(Sigmaxy) %*% solve(Sigmaxx, data_jan10$lnNO2)
X_cov <- Sigmayy - t(Sigmaxy) %*% solve(Sigmaxx,Sigmaxy) # n_y by n_y

# visualize
# monitoring station exposure data
quilt.plot(cbind(data_jan10$lon, data_jan10$lat),
           data_jan10$lnNO2, main = "NO2 exposures (in log) at 21 monitoring stations",
           xlab = "longitude", ylab= "latitude", xlim = c(-96.5, -94.5), ylim = c(29, 30.5))
maps::map("county", "Texas", add = T)

# posterior predictive mean of exposure at health subject locations
quilt.plot(cbind(health_sim$lon, health_sim$lat),
           X_mean, main = "posterior predictive mean of exposure at health subject locations",
           xlab = "longitude", ylab= "latitude", xlim = c(-96.5, -94.5), ylim = c(29, 30.5))
maps::map("county", "Texas", add = T)

# posterior predictive sd of exposure at health subject locations
quilt.plot(cbind(health_sim$lon, health_sim$lat),
           sqrt(diag(X_cov)), main = "posterior predictive sd of exposure at health subject locations",
           xlab = "longitude", ylab= "latitude", xlim = c(-96.5, -94.5), ylim = c(29, 30.5))
maps::map("county", "Texas", add = T)


# vecchia approximation
run_vecchia = vecchia_cov(X_cov, coords = cbind(health_sim$lon, health_sim$lat),
                          n.neighbors = 10)
Q_sparse = run_vecchia$Q
run_vecchia$cputime

# fit the model, binary response
fit = bglm_me(Y = health_sim$Ybinary,
                X_mean = X_mean,
                X_prec = Q_sparse, # sparse precision matrix
                Z = health_sim$Z,
                family = binomial(link = "logit"),
                nburn = 5000,
                nsave = 5000,
                nthin = 1)
fit$cputime
summary(fit$posterior)
library(bayesplot)
bayesplot::mcmc_trace(fit$posterior)
} # }
```
