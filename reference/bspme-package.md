# bspme: Bayesian Spatial Measurement Error Models

Scalable methods for fitting Bayesian linear and generalized linear
models in the presence of spatial exposure measurement error. These
models typically arise from a two-stage Bayesian analysis of
environmental exposures and health outcomes. From a first-stage model,
predictions of the covariate of interest (”exposure”) and their
uncertainty information (typically contained in MCMC samples) are used
to form a multivariate normal prior distribution for exposure in a
second-stage regression model. This package also provides implementation
of the methods used in Lee et al. (2024)
<https://arxiv.org/abs/2401.00634>.

## See also

Useful links:

- <https://changwoo-lee.github.io/bspme/>

- Report bugs at <https://github.com/changwoo-lee/bspme/issues>

## Author

**Maintainer**: Changwoo Lee <c.lee@stat.tamu.edu>

Authors:

- Eun Sug Park
