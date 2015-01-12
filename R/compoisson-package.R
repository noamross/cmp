

#' Operations in Log-space
#' 
#' Computes the difference of two values in log-space.
#' 
#' \code{com.log.difference} computes the difference of two values in
#' log-space, \eqn{log( e^x - e^y )}, without significant chance of overflow or
#' underflow.
#' 
#' \code{com.log.sum} computes the sum of two values in log-space, \eqn{log(
#' e^x + e^y )}, without significant change of overflow or underflow.
#' 
#' \code{com.log.factorial} computes \eqn{log(x!)} which is equivalent to a
#' summation.
#' 
#' @aliases com.log.difference com.log.sum com.log.factorial
#' @param x first value
#' @param y second value
#' @return The requested computation in log-space.
#' @author Jeffrey Dunn
#' @keywords manip
#' @examples
#' 
#' 	a = exp(com.log.difference(log(100), log(20))); # a = 80
#' 	b = exp(com.log.sum(log(100), log(20))); # b = 120
#' 	c = exp(com.log.factorial(4)); # c = 24
#' 
NULL





#' Insurance Count Datasets
#' 
#' Two auto insurance datasets compiled from published works. The Lemaire
#' dataset contains published aggregate claim numbers for automobile
#' third-party liability insurance of a Belgian insurance company in the early
#' 1990's. The Buhlmann dataset originates from aggregate accident claims in
#' 1961 for a class of auto insurance in Switzerland.
#' 
#' 
#' @name compoisson-data
#' @aliases insurance Lemaire Buhlmann
#' @docType data
#' @format Each dataset is a matrix with two columns. The first column contains
#' the levels and the second contains the number of customers who submitted the
#' corresponding level of claims.
#' @source Lemaire, Jean. \dQuote{Bonus-Malus Systems for Automobile
#' Insurance}. Kluwer Academic Publishers, 1995.
#' 
#' Panjer, Harry. \dQuote{Actuarial Mathematics (Proceedings of Symposia in
#' Applied Mathematics)}. Providence: American Mathematical Society, 1986.
#' @keywords datasets
#' @examples
#' 
#' 	data(insurance)
#' 	Lemaire
#' 	Buhlmann
#' 
NULL





#' Conway-Maxwell Poisson Distribution
#' 
#' Provides routines for computing the density of the Conway-Maxwell Poisson
#' distribution and fitting parameters to data.
#' 
#' \tabular{ll}{ Package: \tab compoisson\cr Type: \tab Package\cr Version:
#' \tab 0.2\cr Date: \tab 2008-04-21\cr License: \tab BSD\cr }
#' 
#' @name compoisson-package
#' @aliases compoisson-package compoisson
#' @docType package
#' @author Jeffrey Dunn
#' 
#' Maintainer: Jeffrey Dunn <jsd115@@gmail.com>
#' @seealso See \code{\link{dcom}} for calculating the pmf of the distribution,
#' see \code{\link{com.fit}} for fitting parameters.
#' @references Shmueli, G., Minka, T. P., Kadane, J. B., Borle, S. and
#' Boatwright, P., \dQuote{A useful distribution for fitting discrete data:
#' Revival of the Conway-Maxwell-Poisson distribution,} J. Royal Statist. Soc.,
#' v54, pp. 127-142, 2005.
#' @keywords package models regression
NULL



