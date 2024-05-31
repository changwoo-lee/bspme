#' Run Vecchia approximation given a covariance matrix
#'
#' Given a multivariate normal (MVN) distribution with covariance matrix \eqn{\Sigma},
#' this function finds a sparse precision matrix (inverse covariance) \eqn{Q} based on the Vecchia approximation (Vecchia 1988, Katzfuss and Guinness 2021),
#' where \eqn{N(\mu, Q^{-1})} is the sparse MVN that approximates the original MVN \eqn{N(\mu, \Sigma)}.
#' The algorithm is based on the pseudocode 2 of Finley et al. (2019).
#' Nearest neighbor finding algorithm is based on GpGp::find_ordered_nn.
#'
#' @param Sigma *matrix&lt;num&gt;*, n by n covariance matrix
#' @param coords *matrix&lt;num&gt;*, n by 2 coordinate matrix for nearest neighborhood search
#' @param n.neighbors *integer*, the number of nearest neighbors (k) to determine conditioning set of Vecchia approximation
#' @param ord *vector&lt;int&gt;*, length n vector, ordering of data. If NULL, ordering based on the first coordinate will be used.
#' @param KLdiv *logical*, If TRUE, return KL divergence \eqn{D_{KL}(p || \tilde{p})} where \eqn{p} is multivariate normal with original covariance matrix and \eqn{\tilde{p}} is the approximated multivariate normal with sparse precision matrix.
#' @param lonlat *logical*, (default FALSE) If TRUE, the coordinates are in longitude and latitude. If FALSE, the coordinates are in Euclidean space.
#'
#' @return list of the following:
#' \describe{
#'   \item{Q}{n by n sparse precision matrix in \link{Matrix} format}
#'   \item{ord}{ordering used for Vecchia approximation}
#'   \item{cputime}{time taken to run Vecchia approximation}
#'   \item{KLdiv}{ (if \code{KLdiv = TRUE}) KL divergence \eqn{D_{KL}(p || \tilde{p})} where \eqn{p} is the multivariate normal with original covariance matrix and \eqn{\tilde{p}} is the approximated multivariate normal with a sparse precision matrix.}
#'   \item{Linv}{ n by n sparse reverse Cholesky factor of Q, a lower triangular matrix such that Q = t(Linv)%*%Linv (before ordering changes). In other words, Linv = chol(Q[n:1,n:1])[n:1,n:1] (before ordering changes).}
#' }
#'
#' @references Vecchia, A. V. (1988). Estimation and model identification for continuous spatial processes. Journal of the Royal Statistical Society Series B: Statistical Methodology, 50(2), 297-312.
#'
#' Katzfuss, M., & Guinness, J. (2021). A General Framework for Vecchia Approximations of Gaussian Processes. Statistical Science, 36(1).
#'
#' Finley, A. O., Datta, A., Cook, B. D., Morton, D. C., Andersen, H. E., & Banerjee, S. (2019). Efficient algorithms for Bayesian nearest neighbor Gaussian processes. Journal of Computational and Graphical Statistics, 28(2), 401-414.
#'
#' @importFrom Matrix Diagonal
#'
#' @export
#'
#' @examples
#'
#' library(bspme)
#' n = 1000
#' coords = cbind(runif(n), runif(n))
#' library(GpGp)
#' ord <- GpGp::order_maxmin(coords) # from GpGp
#' Sigma = fields::Exp.cov(coords, aRange = 1)
#' fit5 = vecchia_cov(Sigma, coords, n.neighbors = 5, ord = ord, KLdiv = TRUE)
#' fit5$KLdiv
#' fit10 = vecchia_cov(Sigma, coords, n.neighbors = 10, ord = ord, KLdiv = TRUE)
#' fit10$KLdiv
#'
#' # Check Linv
#' Q = fit5$Q
#' Q_reordered = Q[ord,ord]
#' all.equal(fit5$Linv, chol(Q_reordered[n:1,n:1])[n:1,n:1])
#' all.equal(Q_reordered, t(fit5$Linv) %*% fit5$Linv)
#'
vecchia_cov = function(Sigma, coords, n.neighbors, ord = NULL, KLdiv = FALSE, lonlat = FALSE){
  start.time = Sys.time()

  n = ncol(Sigma)
  if(nrow(Sigma) != n) stop("Sigma must be square matrix")
  if(nrow(coords) != n) stop("number of coordinates are not same as dimension of Sigma")
  if (missing(ord)) {
    ord <- order(coords[, 1])
  }
  else {
    if (length(ord) != n) {
      stop("error: supplied order vector ord must be of length n")
    }
  }
  if(n.neighbors == 0){
    Q = Matrix::Diagonal(n, x=1/diag(Sigma))
    Q = as(Q, "dgCMatrix")
    out = list(Q = Q, ord = ord, cputime =  difftime(Sys.time(),start.time, units = "secs"))
    if(KLdiv){
      out$KLdiv = as.numeric(0.5*(- Matrix::determinant(Q)$modulus - determinant(Sigma)$modulus))
    }
    return(out)
  }

  Sigma = Sigma[ord,ord] # Sigma is now reordered, overwrite to reduce memory usage

  # GpGp::find_ordered_nn allows lon,lat input;
  NNarray <- GpGp::find_ordered_nn(coords[ord,],n.neighbors, lonlat = lonlat)
  NNarray = NNarray[,-1] # remove first column
  Nlist = list()
  Nlist[[1]] = c()
  for(i in 1:(n-1)){
    Nlist[[i+1]] <- NNarray[i+1,!is.na(NNarray[i+1,])]
  }

  # sparse matrix
  A = Matrix::Matrix(0, n, n);
  D = numeric(n) #D = Matrix(0, n, n);
  D[1] = 1
  # pseudocode 2 of Finley et al.
  #Nlist[[i]] be the set of indices j  <i such that A[i,j] not qe 0
  for(i in 1:(n-1)) {
    temp =  solve(Sigma[Nlist[[i+1]], Nlist[[i+1]]],Sigma[Nlist[[i+1]], i+1])
    A[i+1, Nlist[[i+1]]]= temp
    D[i+1]= Sigma[i+1, i+1] - sum(Sigma[i+1, Nlist[[i+1]]]*temp)
  }
  Dmat = Matrix::Diagonal(n, x = D)
  Q_reordered = Matrix::t((Matrix::Diagonal(n) - A))%*%Matrix::solve(Dmat)%*%(Matrix::Diagonal(n) - A)
  Linv =  Matrix::Diagonal(n, x = 1/sqrt(D))%*%(Matrix::Diagonal(n) - A)

  Q = Q_reordered[order(ord),order(ord)]
  Q = as(Q, "dgCMatrix")

  end.time = Sys.time()

  out = list()
  out$Q = Q
  out$ord = ord
  out$cputime = difftime(end.time,start.time, units = "secs")
  if(KLdiv){
    out$KLdiv = as.numeric(0.5*(sum(Q*Sigma[order(ord),order(ord)]) - n - Matrix::determinant(Q)$modulus - determinant(Sigma)$modulus))
  }
  out$Linv = Linv
  return(out)
}


