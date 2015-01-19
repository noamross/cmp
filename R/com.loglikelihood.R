#' Computes Log-Likelihood of COM-Poisson
#'
#' Given a set of data, computes the log-likelihood of the data under the
#' COM-Poisson distribution for values of the parameters.
#'
#' The argument x should consist of a matrix where the first column is the
#' level and the second column is the count for the corresponding level.
#'
#' @param x matrix of count data
#' @param lambda value of lambda parameter
#' @param nu value of nu parameter
#' @return The log-likelihood of the data.
#' @author Jeffrey Dunn
#' @seealso \code{\link{com.fit}}, \code{\link{dcom}}
#' @references Shmueli, G., Minka, T. P., Kadane, J. B., Borle, S. and
#' Boatwright, P., \dQuote{A useful distribution for fitting discrete data:
#' Revival of the Conway-Maxwell-Poisson distribution,} J. Royal Statist. Soc.,
#' v54, pp. 127-142, 2005.
#' @keywords models
#' @export com.loglikelihood
com.loglikelihood = function(x, lambda, nu, ...) {
  if (lambda < 0 || nu < 0)
      return (-Inf)
  log.z = com.compute.log.z(lambda, nu, ...)
  return (x[,2] %*% ( x[,1] * log(lambda) - nu * lfactorial(x[,1]) - log.z ))
}
