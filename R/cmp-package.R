#' Conway-Maxwell Poisson Distribution
#'
#' Provides routines for computing the density of the Conway-Maxwell Poisson
#' distribution and fitting parameters to data.
#'
#' \tabular{ll}{ Package: \tab cmp\ cr Type: \tab Package\cr Version:
#' \tab 0.2\cr Date: \tab 2008-04-21\cr License: \tab BSD\cr }
#'
#' @name cmp-package
#' @aliases cmp-package cmp
#' @docType package
#' @author Noam Ross
#' @references Shmueli, G., Minka, T. P., Kadane, J. B., Borle, S. and
#' Boatwright, P., \dQuote{A useful distribution for fitting discrete data:
#' Revival of the Conway-Maxwell-Poisson distribution,} J. Royal Statist. Soc.,
#' v54, pp. 127-142, 2005.
#' @keywords package models regression
NULL

#' Insurance Count Datasets
#'
#' Two auto insurance datasets compiled from published works. The Lemaire
#' dataset contains published aggregate claim numbers for automobile
#' third-party liability insurance of a Belgian insurance company in the early
#' 1990's. The Buhlmann dataset originates from aggregate accident claims in
#' 1961 for a class of auto insurance in Switzerland.
#'
#' @name cmp-data
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
#'  data(insurance)
#'  Lemaire
#'  Buhlmann
#'
NULL

#' @useDynLib cmp
#' @importFrom Rcpp sourceCpp
#' @import RcppParallel
#' @import RcppArmadillo

NULL
