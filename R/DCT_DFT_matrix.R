#' Create a n-by-n discrete cosine transform matrix.
#'
#' The discrete cosine transform (DCT) matrix for a given dimension n is
#' calculated.
#'
#' The function can be used for 1D- or 2D-DCT transforms of data.
#' \itemize{
#' \item \strong{1D:} Let \code{Q} be a m-by-n matrix with some data. \code{D} is a
#' m-by-m DCT matrix created by \code{dctMatrix(m)}. Then \code{D \%*\% Q} returns the
#' discrete cosine transform of the columns of Q. \code{t(D) \%*\% Q} returns the
#' inverse DCT of the columns of Q. As D is orthogonal, \code{solve(D) = t(D)}.
#' \item \strong{2D:} Let \code{Q} be a m-by-n matrix with some data. \code{D_m} is a
#' m-by-m DCT matrix created by \code{dctMatrix(m)}, \code{D_n} a n-by-n DCT matrix
#' created by \code{dctMatrix(n)}. \code{D_m \%*\% Q \%*\% t(D_n)} computes the 2D-DCT
#' of Q. The inverse 2D-DCT of Q can be computed via \cr \code{t(D_mm) \%*\% DCT_Q \%*\% D_n}.
#' D_m transforms along columns, D_n along rows. Since D is orthogonal, \code{solve(D) = t(D)}.
#' }
#' It can be faster to use \code{dctMatrix} than using a direct transformation,
#' especially when calculating several DCT's.
#'
#' @param n Dimension for the DCT matrix.
#' @return The n-by-n DCT matrix.
#' @export
#' @examples
#' D <- dctMatrix(5)
dctMatrix <- function(n){
  dctMat <- matrix(0, nrow = n, ncol = n)
  for (i in 1:n) {
    for (j in 1:n) {
      if (i == 1) {
        dctMat[i, j] <- 1 / sqrt(n)
      }
      if (i > 1) {
        dctMat[i, j] <- sqrt(2 / n) * cos((2 * (j - 1) + 1) * (i - 1) * pi / (2 * n))
      }
    }
  }
  return(dctMat)
}




#' Create a n-by-n discrete Fourier transform matrix.
#' 
#' The discrete Fourier transform (DFT) matrix for a given dimension n is
#' calculated.
#' 
#' The DFT matrix can be used for computing the discrete Fourier transform of 
#' a matrix or vector. \code{dftMatrix(n) \%*\% testMatrix} is the same as
#' \code{apply(testMatrix, MARGIN = 2, FUN = fft)}.
#' 
#' @param n Dimension for the DFT matrix.
#' @return The n-by-n DFT matrix.
#' @export
#' @examples
#' set.seed(987)
#' testMatrix <- matrix(sample(1:10, size = 25, replace = TRUE), nrow = 5)
#' D <- dftMatrix(5)
#' 
#' # Discrete Fourier transform with matrix multiplication:
#' D %*% testMatrix 
#' 
#' # Discrete Fourier transform with function fft: 
#' apply(testMatrix, MARGIN = 2, FUN = fft)
#' 
dftMatrix <- function(n) {
  omega <- exp(-2 * pi * (0 + 1i) / n)
  sq <- rep(0:(n - 1), n) * rep(0:(n - 1), each = n)
  AA <- sapply(sq, function(nn, omega) omega^nn, omega)
  dft.mat <- matrix(AA, n, n)
  return(dft.mat)
}

