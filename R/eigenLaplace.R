#' Generate eigenvalues of discrete Laplace matrix.
#'
#' The eigenvalues of a discrete Laplace matrix with dimension (\code{mm}, \code{nn}) are
#' calculated.
#'
#' @param mm Number of rows of the discrete Laplace matrix.
#' @param nn Number of columns of the discrete Laplace matrix.
#' @return A row vector containing the eigenvalues of the discrete
#'     laplace matrix.
#' @export
#' @examples
#' eigval <- eigenLaplace(5, 5)
#'
eigenLaplace <- function(mm, nn) {

  N <- mm * nn
  Mu <- rep(0, N)
  In <- seq(0, mm - 1, 1)
  Im <- seq(0, nn - 1, 1)
  In <- t(In)
  Im <- t(Im)

  mu <- 2 - 2 * cos(pi * Im / nn)
  lambda <- 2 - 2 * cos(pi * In / mm)

  for (k in 1:nn) {
    for (s in 1:mm) {
      Mu[(k - 1) * mm + s] <- lambda[s] + mu[k]
    }
  }

  return(Mu)
}

