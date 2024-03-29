# mrbsizeR

## Installation

The preferred way to install mrbsizeR from CRAN is:
``` r
install.packages("mrbsizeR")
```

You can install the mrbsizeR devel version from github with:
``` r
# install.packages("devtools")
devtools::install_github("romanflury/mrbsizeR")
```


## Example

``` r
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
plot(x = mrbOut$smMean, turn_out = TRUE)

# Credibility analysis using simultaneous credible intervals
plot(x = mrbOut$ciout, turn_out = TRUE) 
```
[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/mrbsizeR)](https://cran.r-project.org/package=mrbsizeR)
[![](https://cranlogs.r-pkg.org/badges/mrbsizeR)](https://cran.r-project.org/package=mrbsizeR)
