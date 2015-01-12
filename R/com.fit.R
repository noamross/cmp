#' Computes COM-Poisson Regression
#' 
#' Computes the maximum likelihood estimates of the COM-Poisson model for given
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
#' @seealso \code{\link{com.compute.z}}, \code{\link{com.loglikelihood}}
#' @references Shmueli, G., Minka, T. P., Kadane, J. B., Borle, S. and
#' Boatwright, P., \dQuote{A useful distribution for fitting discrete data:
#' Revival of the Conway-Maxwell-Poisson distribution,} J. Royal Statist. Soc.,
#' v54, pp. 127-142, 2005.
#' @keywords models regression
#' @examples
#' 
#' 	data(insurance)
#' 	com.fit(Lemaire);
#' 
#' @export com.fit
com.fit = function(x, par=NULL)
{
	if (is.null(par)) {
		xbar = (x[,1] %*% x[,2]) / sum(x[,2]);
		par = c(xbar, 1)
	}
#	options(warn = -1);
	result = optim(par, function(p) {return (-com.loglikelihood(x, p[1], p[2]));},
		method="L-BFGS-B", lower=c(1e-10,0));
#	options(warn = 0);
	
	lambda = result$par[1];
	nu = result$par[2];
	fit = list( lambda = lambda,
	            nu = nu,
	            z = com.compute.z(lambda, nu),
	            fitted.values = sum(x[,2]) * dcom(x[,1], lambda, nu),
				log.likelihood = com.loglikelihood(x, lambda, nu) );

	return (fit);
}
