
<!-- README.md is generated from README.Rmd. Please edit that file -->

# bspme

<!-- badges: start -->
<!-- badges: end -->

Go to the package website:
[\[link\]](https://changwoo-lee.github.io/bspme/)

See a vignette with NO2 exposure data:
[\[link\]](https://changwoo-lee.github.io/bspme/articles/no2-exposure-and-health-data-analysis.html)

See bspme_1.0.0.pdf for the pdf file of the package manual.

**bspme** is an R package that provides fast, scalable inference tools
for **B**ayesian **sp**atial exposure **m**easurement **e**rror models,
namely, the Bayesian linear and generalized linear models with the
presence of spatial exposure measurement error of covariate(s). These
models typically arise from a two-stage Bayesian analysis of
environmental exposures and health outcomes. From a first-stage model,
predictions of the covariate of interest (“exposure”) and their
uncertainty information (typically contained in MCMC samples) are
obtained and used to form a multivariate normal prior distribution
$X\sim N(\mu, \Sigma)$ for exposure in a second-stage regression model.
Naive, non-sparse choices of the precision matrix $Q = \Sigma^{-1}$ of
the multivariate normal (such as a sample precision matrix) lead to the
MCMC posterior inference algorithm being infeasible to run for a large
number of subjects $n$ because of the cubic computational cost
associated with the $n$-dimensional MVN prior. With a sparse precision
matrix $Q$ obtained from the Vecchia approximation, the **bspme**
package offers fast, scalable algorithms to conduct posterior inference
for large health datasets, with the number of subjects $n$ possibly
reaching tens of thousands. For more details, please see the following
paper:

> Lee, C. J., Symanski, E., Rammah, A., Kang, D. H., Hopke, P. K., &
> Park, E. S. (2024). A scalable two-stage Bayesian approach accounting
> for exposure measurement error in environmental epidemiology. arXiv
> preprint arXiv:2401.00634. <https://arxiv.org/abs/2401.00634> \##
> Installation

You can install the development version of bspme with the following
code:

``` r
# install.packages("devtools")
devtools::install_github("changwoo-lee/bspme", build_vignettes = T)
```

To browse and see vignettes, run

``` r
browseVignettes("bspme")
```

## Functionality

| Function        | Description                                                                 |
|-----------------|-----------------------------------------------------------------------------|
| `blm_me()`      | Bayesian linear regression models with spatial exposure measurement error.  |
| `bglm_me()`     | Bayesian generalized linear models with spatial exposure measurement error. |
| `vecchia_cov()` | Run Vecchia approximation given a covariance matrix.                        |

To see function description in R environment, run the following lines:

``` r
?blm_me
?bglm_me
?vecchia_cov
```

## datasets

| Dataset call          | Description                                                                       |
|-----------------------|-----------------------------------------------------------------------------------|
| `data("NO2_Jan2012")` | Daily average NO2 concentrations in and around Harris County, Texas, in Jan 2012. |
| `data("health_sim")`  | Simulated health data associated with ln(NO2) concentration on Jan 10, 2012.      |

## Examples

Please see the vignette “NO2-exposure-and-health-data-analysis”.
