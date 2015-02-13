#' Computes CMP Regression
#'
#' Computes the maximum likelihood estimates of the CMP model for given
#' count data.
#'
#' The argument x should consist of a matrix where the first column is the
#' level and the second column is the count for the corresponding level.
#'
#' @param x matrix of count data
#' @param par initial guesses for \deqn{lambda} and \deqn{nu}
#' @return Returns an object containing four fields: \item{lambda }{ Estimate
#' of the lambda parameter} \item{nu }{ Estimate of the nu parameter} \item{z
#' }{ Normalizing constant} \item{fitted.values }{ Estimated counts at given
#' levels }
#' @author Jeffrey Dunn
#' @seealso \code{\link{com_compute_z}}, \code{\link{com.loglikelihood}}
#' @references Shmueli, G., Minka, T. P., Kadane, J. B., Borle, S. and
#' Boatwright, P., \dQuote{A useful distribution for fitting discrete data:
#' Revival of the Conway-Maxwell-Poisson distribution,} J. Royal Statist. Soc.,
#' v54, pp. 127-142, 2005.
#' @keywords models regression
#' @examples
#'
#'   data(insurance)
#'   cmp_fit(Lemaire)
#'
#' @importFrom stats optim
#' @export cmp_fit
cmp_fit = function(x, par=NULL, ...) {

  if (is.null(par)) {
    xbar = (x[,1] %*% x[,2]) / sum(x[,2])
    par = c(xbar, 1)
  }

  options(warn = -1)
  result = optim(log(par),
                 function(p) {return (-cmp_loglik(x, exp(p[1]), exp(p[2])))},
                 )
  options(warn = 0)

  lambda = exp(result$par[1])
  nu = exp(result$par[2])
  fit = list(lambda = lambda,
             nu = nu,
             z = compute_z(lambda, nu),
             fitted.values = sum(x[,2]) * dcmp(x[,1], lambda, nu),
             log.likelihood = cmp_loglik(x, lambda, nu))
  return (fit)
}

#' @importFrom stats optim
#' @export
pois_fit = function(x, par=NULL) {

  if (is.null(par)) {
    par = (x[,1] %*% x[,2]) / sum(x[,2])
  }

  if (is.nan(par) || par == 0) {
    lambda = 0
  } else {
    options(warn = -1)
    result = optim(log(par),
                   function(p) {return (-pois_loglik(x, exp(p[1])))},
    )
    options(warn = 0)
    
    lambda = exp(result$par[1])
  }
  
  fit = list(lambda = lambda,
             fitted.values = sum(x[,2]) * dpois(x[,1], lambda),
             log.likelihood = pois_loglik(x, lambda))
  return (fit)
}

#' @importFrom stats optim
#' @export nb_fit
nb_fit = function(x, par=NULL) {

  if (is.null(par)) {
    xbar = (x[,1] %*% x[,2]) / sum(x[,2])
    par = c(xbar, 1e6)
  }
  
  if (is.nan(par[1]) || par[1] == 0) {
    mu = 0
    size = 1
  } else {
    
    options(warn = -1)
    result = optim(log(par),
                   function(p) {return (-nb_loglik(x, exp(p[1]), exp(p[2])))},
    )
    options(warn = 0)
    
    mu = exp(result$par[1])
    size = exp(result$par[2])
  }
  fit = list(mu = mu,
             size = size,
             fitted.values = sum(x[,2]) * dnbinom(x = x[,1], mu=mu, size=size),
             log.likelihood = nb_loglik(x, mu, size))
  return (fit)
}

