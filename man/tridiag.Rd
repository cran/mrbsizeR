% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/HelpFunctions.R
\name{tridiag}
\alias{tridiag}
\title{Generate a tridiagonal matrix.}
\usage{
tridiag(mainDiag, upperDiag, lowerDiag)
}
\arguments{
\item{mainDiag}{Diagonal of tridiagonal matrix.}

\item{upperDiag}{Superdiagonal of tridiagonal matrix. Must have length \code{length(mainDiag) - 1}.}

\item{lowerDiag}{Subdiagonal of tridiagonal matrix. Must have length \code{length(mainDiag) - 1}.}
}
\value{
Tridiagonal matrix.
}
\description{
Generate a tridiagonal matrix with \code{upperDiag} as superdiagonal, 
\code{lowerDiag} as subdiagonal and \code{mainDiag} as diagonal.
}
\examples{
set.seed(987)
mainDiag <- sample(100:110, size = 6, replace = TRUE)
upperDiag <- sample(10:20, size = 5, replace = TRUE)
lowerDiag <- sample(1:10, size = 5, replace = TRUE)

tridiag(mainDiag, upperDiag, lowerDiag)  
 
}
