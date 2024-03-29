% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/blm_me.R
\name{blm_me}
\alias{blm_me}
\title{Bayesian linear regression models with spatial exposure measurement error.}
\usage{
blm_me(
  Y,
  X_mean,
  X_prec,
  Z,
  nburn = 5000,
  nsave = 5000,
  nthin = 1,
  prior = NULL,
  saveX = FALSE
)
}
\arguments{
\item{Y}{\emph{vector<int>}, n by 1 continuous response vector.}

\item{X_mean}{\emph{vector<num>}, n by 1 prior mean vector \eqn{\mu_X}. When there are q multiple exposures subject to measurement error, it can be a length q list of n by 1 vectors.}

\item{X_prec}{\emph{matrix<num>}, n by n prior precision matrix \eqn{Q_X}, which allows sparse format from \link{Matrix} package. When there are q multiple exposures subject to measurement error, it can be a length q list of n by n matrices.}

\item{Z}{\emph{matrix<num>}, n by p matrix containing p covariates that are not subject to measurement error.}

\item{nburn}{\emph{integer}, number of burn-in iterations (default=5000).}

\item{nsave}{\emph{integer}, number of posterior samples (default=5000). Total number of MCMC iteration is \code{nburn + nsave * nthin}.}

\item{nthin}{\emph{integer}, thin-in rate (default=1).}

\item{prior}{\emph{list}, list of prior parameters of the regression model. Default is \code{list(var_beta = 100, a_Y = 0.01, b_Y = 0.01)}.}

\item{saveX}{\emph{logical}, default FALSE, whether save posterior samples of X (exposure).}
}
\value{
list of the following:
\describe{
\item{posterior}{\code{nsave} by (q + p + 1) matrix of posterior samples of \eqn{\beta_X}, \eqn{\beta_Z}, \eqn{\sigma_Y^2} as a coda::\link[coda]{mcmc} object.}
\item{time}{time taken for running MCMC (in seconds)}
\item{X_save}{ (if \code{saveX = TRUE}) posterior samples of X}
}
}
\description{
This function fits a Bayesian linear regression model in the presence of spatial exposure measurement error for covariate(s) \eqn{X}.
One of the most important features of this function is that it allows a sparse matrix input for the prior precision matrix of \eqn{X} for scalable computation.
Function \code{blm_me()} runs a Gibbs sampler to carry out posterior inference;  see the "Details" section below for the model description, and Lee et al. (2024) for an application example in environmental epidemiology.
}
\details{
Let \eqn{Y_i} be a continuous response, \eqn{X_i} be a \eqn{q\times 1} covariate vector that is subject to spatial exposure measurement error,
and \eqn{Z_i} be a \eqn{p\times 1} covariate vector without measurement error.
Consider a normal linear regression model,
\deqn{Y_i = \beta_0 + X_i^\top \beta_X +  Z_i^\top \beta_Z + \epsilon_i,\quad  \epsilon_i \stackrel{iid}{\sim} N(0, \sigma^2_Y), \quad i=1,\dots,n.}
Spatial exposure measurement error of \eqn{X_i} for \eqn{i=1,\dots,n} is incorporated into the model using a multivariate normal prior.
For example when \eqn{q=1}, we have an \eqn{n-}dimensional multivariate normal prior on \eqn{X = (X_1,\dots,X_n)^\top},
\deqn{(X_1,\dots,X_n)\sim N_n(\mu_X, Q_X^{-1}).}
Most importantly, it allows a sparse matrix input for the prior precision matrix \eqn{Q_X} for scalable computation, which can be obtained by Vecchia approximation.
When \eqn{q>1}, \eqn{q} independent \eqn{n-}dimensional multivariate normal priors are assumed.

We consider semiconjugate priors for regression coefficients and error variance,
\deqn{\beta_0 \sim N(0, V_\beta), \quad \beta_{X,j} \stackrel{iid}{\sim} N(0, V_\beta), \quad \beta_{Z,k} \stackrel{iid}{\sim} N(0, V_\beta), \quad \sigma_Y^2 \sim IG(a_Y, b_Y).}
where \code{var_beta} corresponds to \eqn{V_\beta}, and \code{a_Y} and \code{b_Y} correspond to hyperparameters of an inverse gamma prior for \eqn{\sigma^2_Y}.
}
\examples{
\dontrun{
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
X_mean <- t(Sigmaxy) \%*\% solve(Sigmaxx, data_jan10$lnNO2)
X_cov <- Sigmayy - t(Sigmaxy) \%*\% solve(Sigmaxx,Sigmaxy) # n_y by n_y

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

# fit the model, continuous response
fit = blm_me(Y = health_sim$Y,
                X_mean = X_mean,
                X_prec = Q_sparse, # sparse precision matrix
                Z = health_sim$Z,
                nburn = 5000,
                nsave = 5000,
                nthin = 1)
fit$cputime
summary(fit$posterior)
library(bayesplot)
bayesplot::mcmc_trace(fit$posterior)
}

}
\references{
Lee, C. J., Symanski, E., Rammah, A., Kang, D. H., Hopke, P. K., & Park, E. S. (2024). A scalable two-stage Bayesian approach accounting for exposure measurement error in environmental epidemiology. arXiv preprint arXiv:2401.00634.
}
