#' Computes Expectation of a Function of a COM-Poisson Random Variable
#' 
#' Computes an expectation of a function of a COM-Poisson random variable.
#' 
#' Computes the expectation \eqn{E[f(X)]}{E[f(X)]} where X is a COM-Poisson
#' random variable.
#' 
#' @param f function taking as a single argument the value of x
#' @param lambda value of lambda parameter
#' @param nu value of nu parameter
#' @param log.error precision in the log of the expectation
#' @return The expectation as a real number.
#' @author Jeffrey Dunn
#' @seealso \code{\link{com.mean}}, \code{\link{com.var}},
#' \code{\link{com.fit}}
#' @references Shmueli, G., Minka, T. P., Kadane, J. B., Borle, S. and
#' Boatwright, P., \dQuote{A useful distribution for fitting discrete data:
#' Revival of the Conway-Maxwell-Poisson distribution,} J. Royal Statist. Soc.,
#' v54, pp. 127-142, 2005.
#' @keywords models
#' @importFrom matrixStats logSumExp
#' @export com.expectation
com.expectation = function(f, lambda, nu, log.error = 0.001)
{
	log.z = com.compute.log.z(lambda, nu);

	# Initialize variables
	ex = -.Machine$double.xmax;
	ex.last = 0;
	j = 0;

	# Continue until we have reached specified precision
	while ((ex == -.Machine$double.xmax && ex.last == -.Machine$double.xmax) || abs(ex - ex.last) > log.error)
	{
		ex.last = ex;
		ex = logSumExp(c(ex, log(f(j)) + com.log.density(j, lambda, nu, log.z)));
		j = j + 1;
	}
	return (exp(ex));
}





#' Computes Mean of the COM-Poisson Distribution
#' 
#' Computes the mean of the COM-Poisson distribution for given values of the
#' parameters.
#' 
#' Uses \code{\link{com.expectation}} to compute the first moment of the
#' distribution.
#' 
#' @param lambda value of lambda parameter
#' @param nu value of the nu parameter
#' @return The mean of the distribution.
#' @author Jeffrey Dunn
#' @seealso \code{\link{com.expectation}}, \code{\link{com.var}}
#' @references Shmueli, G., Minka, T. P., Kadane, J. B., Borle, S. and
#' Boatwright, P., \dQuote{A useful distribution for fitting discrete data:
#' Revival of the Conway-Maxwell-Poisson distribution,} J. Royal Statist. Soc.,
#' v54, pp. 127-142, 2005.
#' @keywords models
#' @examples
#' 
#' 	data(insurance)
#' 	model = com.fit(Lemaire);
#' 	com.mean(model$lambda, model$nu);
#' 
#' @export com.mean
com.mean = function(lambda, nu)
{
	return ( com.expectation(function (x) {x}, lambda, nu) );
}





#' Computes Variance of the COM-Poisson Distribution
#' 
#' Computes the variance of the COM-Poisson distribution for given values of
#' the parameters.
#' 
#' Uses \code{\link{com.expectation}} to compute the second moment of the
#' distribution and subtracts the squared mean, computed using
#' \code{\link{com.mean}}.
#' 
#' @param lambda value of lambda parameter
#' @param nu value of the nu parameter
#' @return The variance of the distribution.
#' @author Jeffrey Dunn
#' @seealso \code{\link{com.expectation}}, \code{\link{com.mean}}
#' @references Shmueli, G., Minka, T. P., Kadane, J. B., Borle, S. and
#' Boatwright, P., \dQuote{A useful distribution for fitting discrete data:
#' Revival of the Conway-Maxwell-Poisson distribution,} J. Royal Statist. Soc.,
#' v54, pp. 127-142, 2005.
#' @keywords models
#' @examples
#' 
#' 	data(insurance)
#' 	model = com.fit(Lemaire);
#' 	com.var(model$lambda, model$nu);
#' 
#' @export com.var
com.var = function(lambda, nu)
{
	return ( com.expectation(function(x) {x^2}, lambda, nu) - (com.mean(lambda,nu))^2 );
}


