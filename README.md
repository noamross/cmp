# An R Package for the Conway-Maxwell-Poisson Distribution

[![Project Status: Wip - Initial development is in progress, but there has not yet been a stable, usable release suitable for the public.](http://www.repostatus.org/badges/0.1.0/wip.svg)](http://www.repostatus.org/#wip)
[![Build Status](https://travis-ci.org/noamross/compoisson.png?branch=master)](https://travis-ci.org/noamross/compoisson)
[![Coverage Status](https://coveralls.io/repos/noamross/compoisson/badge.png?style=flat)](https://coveralls.io/r/noamross/compoisson)

This in-development package implements distribution, likelihood, and fitting functions of the
[Conway-Maxwell-Poisson (CMP) Distribution](http://en.wikipedia.org/wiki/Conway%E2%80%93Maxwell%E2%80%93Poisson_distribution),
which is suitable for both under- and over-dispersed count data.

For an explanation of the distribution see the [Shiny App](https://noamross.shinyapps.io/compois/Plotting-compoisson-distributions.Rmd).

Several other packages implement various approaches to the CMP distribution. This
package aims to improve both speed and numerical accuracy by

-   Using C++ rather than R for all essential functions
-   Parallelizing as possible using RcppParallel
-   Working in log-space to prevent numerical overflow or underflow
-   Using forking algorithms that take different strategies based on the
    parameters
-   Calculating stopping criteria using bounded errors

Many of these strategies are described in:

*Shmueli, G., Minka, T. P., Kadane, J. B., Borle, S. and Boatwright, P. (2005) A
useful distribution for fitting discrete data: Revival of the 
Conway-Maxwell-Poisson distribution.  J. Royal Statist. Soc. vol. 54, pp. 127-142, 2005.*

The package also owes much to Jeffery Dunn's [compoisson](http://cran.r-project.org/web/packages/compoisson/)
package, which implements some of these strategies in pure R.


## Install

To install in R, run

```
library(devtools) # If you don't have it, run install.packages('devtools')
install_github('noamross/cmp')
```


