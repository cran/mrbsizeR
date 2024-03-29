% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot.HPWmapGrid.R
\name{plot.HPWmapGrid}
\alias{plot.HPWmapGrid}
\title{Plotting of pointwise and highest pointwise probabilities.}
\usage{
\method{plot}{HPWmapGrid}(
  x,
  plotWhich = "Both",
  color = c("firebrick1", "gainsboro", "dodgerblue3"),
  turnOut = TRUE,
  title,
  aspRatio = 1,
  ...
)
}
\arguments{
\item{x}{List containing the pointwise (PW) and highest pointwise (HPW)
probabilities of all differences of smooths.}

\item{plotWhich}{Which probabilities shall be plotted? \code{HPW}, \code{PW}
or \code{Both}?}

\item{color}{Vector of length 3 containing the colors to be used in the
credibility maps. The first color represents the credibly negative pixels,
the second color the pixels that are not credibly different from zero
and the third color the credibly positive pixels.}

\item{turnOut}{Logical. Should the output images be turned 90 degrees
counter-clockwise?}

\item{title}{Vector containing one string per plot. The required
number of titles is equal to \code{length(mrbOut$hpout)}. If no \code{title}
is passed, defaults are used.}

\item{aspRatio}{Adjust the aspect ratio of the plots. The default \code{aspRatio = 1}
produces square plots.}

\item{...}{Further graphical parameters can be passed.}
}
\value{
Plots of pointwise and/or highest pointwise probabilities for all
    differences of smooths are created.
}
\description{
Maps with pointwise (PW) probabilities and/or highest pointwise (HPW)
probabilities of all differences of smooths at neighboring scales are plotted.
}
\details{
The default colors of the maps have the following meaning:
\itemize{
\item \strong{Blue}: Credibly positive pixels.
\item \strong{Red}: Credibly negative pixels.
\item \strong{Grey}: Pixels that are not credibly different from zero.
}
\code{x} corresponds to the \code{hpout}-part of the
output of \code{\link{mrbsizeRgrid}}.
}
\examples{
# Artificial sample data
set.seed(987)
sampleData <- matrix(stats::rnorm(100), nrow = 10)
sampleData[4:6, 6:8] <- sampleData[4:6, 6:8] + 5

# Generate samples from multivariate t-distribution
tSamp <- rmvtDCT(object = sampleData, lambda = 0.2, sigma = 6, nu0 = 15,
                   ns = 1000)

# mrbsizeRgrid analysis
mrbOut <- mrbsizeRgrid(posteriorFile = tSamp$sample, mm = 10, nn = 10,
                       lambdaSmoother = c(1, 1000), prob = 0.95)

# Posterior mean of the differences of smooths
plot(x = mrbOut$smMean, turnOut = TRUE)

# Credibility analysis using pointwise (PW) maps
plot(x = mrbOut$hpout, plotWhich = "PW", turnOut = TRUE)

# Credibility analysis using highest pointwise probability (HPW) maps
plot(x = mrbOut$hpout, plotWhich = "HPW", turnOut = TRUE)

}
