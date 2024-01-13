#' Run Vecchia approximation given a covariance matrix
#'
#' Given a multivariate normal (MVN) distribution with covariance matrix \eqn{\Sigma},
#' this function finds a sparse precision matrix (inverse covariance) \eqn{Q} based on the Vecchia approximation (Vecchia 1988, Katzfuss and Guinness 2021),
#' where \eqn{N(\mu, Q^{-1})} is the sparse MVN that approximates the original MVN \eqn{N(\mu, \Sigma)}.
#' The algorithm is based on the pseudocode 2 of Finley et al. (2019).
#'
#' @param Sigma *matrix&lt;num&gt;*, n by n covariance matrix
#' @param coords *matrix&lt;num&gt;*, n by 2 coordinate matrix for nearest neighborhood search
#' @param n.neighbors *integer*, the number of nearest neighbors (k) to determine conditioning set of Vecchia approximation
#' @param ord *vector&lt;int&gt;*, length n vector, ordering of data. If NULL, ordering based on the first coordinate will be used.
#' @param KLdiv *logical*, If TRUE, return KL divergence \eqn{D_{KL}(p || \tilde{p})} where \eqn{p} is multivariate normal with original covariance matrix and \eqn{\tilde{p}} is the approximated multivariate normal with sparse precision matrix.
#'
#' @return list of the following:
#' \describe{
#'   \item{Q}{n by n sparse precision matrix in \link{Matrix} format}
#'   \item{ord}{ordering used for Vecchia approximation}
#'   \item{cputime}{time taken to run Vecchia approximation}
#'   \item{KLdiv}{ (if \code{KLdiv = TRUE}) KL divergence \eqn{D_{KL}(p || \tilde{p})} where \eqn{p} is the multivariate normal with original covariance matrix and \eqn{\tilde{p}} is the approximated multivariate normal with a sparse precision matrix.}
#' }
#'
#' @references Vecchia, A. V. (1988). Estimation and model identification for continuous spatial processes. Journal of the Royal Statistical Society Series B: Statistical Methodology, 50(2), 297-312.
#'
#' Katzfuss, M., & Guinness, J. (2021). A General Framework for Vecchia Approximations of Gaussian Processes. Statistical Science, 36(1).
#'
#' Finley, A. O., Datta, A., Cook, B. D., Morton, D. C., Andersen, H. E., & Banerjee, S. (2019). Efficient algorithms for Bayesian nearest neighbor Gaussian processes. Journal of Computational and Graphical Statistics, 28(2), 401-414.
#'
#' Zhang, L., (2020), public Github repository <https://github.com/LuZhangstat/NNGP_STAN>.
#'
#' @importFrom Matrix Diagonal
#'
#' @export
#'
#' @examples
#'
#' n = 1000
#' coords = cbind(runif(n), runif(n))
#' Sigma = fields::Exp.cov(coords, aRange = 1)
#' fit5 = vecchia_cov(Sigma, coords, n.neighbors = 5, KLdiv = TRUE)
#' fit5$KLdiv
#' fit10 = vecchia_cov(Sigma, coords, n.neighbors = 10, KLdiv = TRUE)
#' fit10$KLdiv
#'
vecchia_cov = function(Sigma, coords, n.neighbors, ord = NULL, KLdiv = FALSE){
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

  indx <- spNNGP:::mkNNIndxCB(coords, n.neighbors, n.omp.threads = 1) # search.type == "cb"
  nn.indx <- indx$nnIndx
  n.indx = spNNGP:::mk.n.indx.list(nn.indx, n, n.neighbors)

  if(n.neighbors > 1){
    NN_ind <- t(sapply(1: (n - 1), get_NN_ind, n.indx[-1], n.neighbors))
  }else if(n.neighbors == 1){ # n.neighbors == 1
    NN_ind <- as.matrix(sapply(1: (n - 1), get_NN_ind, n.indx[-1], n.neighbors))
  }else{ # n.neighbors == 0
    Q = Matrix::Diagonal(n, x=1/diag(Sigma))
    Q = as(Q, "dgCMatrix")
    out = list(Q = Q, ord = ord, cputime =  difftime(Sys.time(),start.time, units = "secs"))
    if(KLdiv){
      out$KLdiv = as.numeric(0.5*(- Matrix::determinant(Q)$modulus - determinant(Sigma)$modulus))
    }
    return(out)
  }
  Sigma = Sigma[ord,ord] # Sigma is now reordered, overwrite to reduce memory usage

  Nlist = list()
  Nlist[[1]] = c()
  for(i in 1:(n-1)){
    if(i<=n.neighbors){
      Nlist[[i+1]] <- NN_ind[i,1:i]
    }else{
      Nlist[[i+1]] <- NN_ind[i,]
    }
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
  D = Matrix::Diagonal(n, x = D)
  Q_reordered = Matrix::t((Matrix::Diagonal(n) - A))%*%Matrix::solve(D)%*%(Matrix::Diagonal(n) - A)
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
  return(out)
}






#' helper function, get nearest neighbor indices
#' source: https://github.com/LuZhangstat/NNGP_STAN/blob/master/NNMatrix.R
#' @noRd
get_NN_ind <- function (ind, ind_distM_i, M) {
  if (ind < M ){l = ind } else {l = M};
  D_i <- rep(0, M);
  D_i[1:l] <- c(ind_distM_i[[ind]])[1:l]
  return(D_i)
}



