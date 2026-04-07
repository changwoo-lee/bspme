# Run Vecchia approximation given a covariance matrix

Given a multivariate normal (MVN) distribution with covariance matrix
\\\Sigma\\, this function finds a sparse precision matrix (inverse
covariance) \\Q\\ based on the Vecchia approximation (Vecchia 1988,
Katzfuss and Guinness 2021), where \\N(\mu, Q^{-1})\\ is the sparse MVN
that approximates the original MVN \\N(\mu, \Sigma)\\. The algorithm is
based on the pseudocode 2 of Finley et al. (2019). Nearest neighbor
finding algorithm is based on GpGp::find_ordered_nn.

## Usage

``` r
vecchia_cov(
  Sigma,
  coords,
  n.neighbors,
  ord = NULL,
  KLdiv = FALSE,
  lonlat = FALSE
)
```

## Arguments

- Sigma:

  *matrix\<num\>*, n by n covariance matrix

- coords:

  *matrix\<num\>*, n by 2 coordinate matrix for nearest neighborhood
  search

- n.neighbors:

  *integer*, the number of nearest neighbors (k) to determine
  conditioning set of Vecchia approximation

- ord:

  *vector\<int\>*, length n vector, ordering of data. If NULL, ordering
  based on the first coordinate will be used.

- KLdiv:

  *logical*, If TRUE, return KL divergence \\D\_{KL}(p \|\| \tilde{p})\\
  where \\p\\ is multivariate normal with original covariance matrix and
  \\\tilde{p}\\ is the approximated multivariate normal with sparse
  precision matrix.

- lonlat:

  *logical*, (default FALSE) If TRUE, the coordinates are in longitude
  and latitude. If FALSE, the coordinates are in Euclidean space.

## Value

list of the following:

- Q:

  n by n sparse precision matrix in
  [Matrix](https://rdrr.io/pkg/Matrix/man/Matrix.html) format

- ord:

  ordering used for Vecchia approximation

- cputime:

  time taken to run Vecchia approximation

- KLdiv:

  (if `KLdiv = TRUE`) KL divergence \\D\_{KL}(p \|\| \tilde{p})\\ where
  \\p\\ is the multivariate normal with original covariance matrix and
  \\\tilde{p}\\ is the approximated multivariate normal with a sparse
  precision matrix.

- Linv:

  n by n sparse reverse Cholesky factor of Q, a lower triangular matrix
  such that Q = t(Linv)%\*%Linv (before ordering changes). In other
  words, Linv = chol(Qn:1,n:1)n:1,n:1 (before ordering changes).

## References

Vecchia, A. V. (1988). Estimation and model identification for
continuous spatial processes. Journal of the Royal Statistical Society
Series B: Statistical Methodology, 50(2), 297-312.

Katzfuss, M., & Guinness, J. (2021). A General Framework for Vecchia
Approximations of Gaussian Processes. Statistical Science, 36(1).

Finley, A. O., Datta, A., Cook, B. D., Morton, D. C., Andersen, H. E., &
Banerjee, S. (2019). Efficient algorithms for Bayesian nearest neighbor
Gaussian processes. Journal of Computational and Graphical Statistics,
28(2), 401-414.

## Examples

``` r
library(bspme)
n = 1000
coords = cbind(runif(n), runif(n))
library(GpGp)
ord <- GpGp::order_maxmin(coords) # from GpGp
Sigma = fields::Exp.cov(coords, aRange = 1)
fit5 = vecchia_cov(Sigma, coords, n.neighbors = 5, ord = ord, KLdiv = TRUE)
fit5$KLdiv
#> [1] 15.75398
fit10 = vecchia_cov(Sigma, coords, n.neighbors = 10, ord = ord, KLdiv = TRUE)
fit10$KLdiv
#> [1] 3.531344

# Check Linv
Q = fit5$Q
Q_reordered = Q[ord,ord]
all.equal(fit5$Linv, chol(Q_reordered[n:1,n:1])[n:1,n:1])
#> [1] TRUE
all.equal(Q_reordered, t(fit5$Linv) %*% fit5$Linv)
#> [1] TRUE
```
