#' Compute COM-Poisson Normalizing Constant
#' 
#' Computes the normalizing constant in the COM-Poisson model for given values
#' of the parameters.
#' 
#' \code{com.compute.z} computes the COM-Poisson normalizing constant \deqn{
#' }{z = Sum (lambda^j)/(j!^nu) }\deqn{ z = \sum_{i=0}^\infty
#' \frac{\lambda^j}{(j!)^\nu} }{z = Sum (lambda^j)/(j!^nu) }\deqn{ }{z = Sum
#' (lambda^j)/(j!^nu) } to the specified precision. If no precision is
#' specified, then the package default is used.
#' 
#' \code{com.compute.log.z} is equivalent to \code{log(com.compute.z(lambda,
#' nu))} but provudes additional precision.
#' 
#' @aliases com.compute.z com.compute.log.z
#' @param lambda Lambda value in COM-Poisson distribution
#' @param nu Nu value in COM-Poisson distribution
#' @param log.error Precision in the log of the normalizing constant
#' @return The normalizing constant as a real number with specified precision.
#' @author Jeffrey Dunn
#' @seealso \code{\link{com.fit}}
#' @references Shmueli, G., Minka, T. P., Kadane, J. B., Borle, S. and
#' Boatwright, P., \dQuote{A useful distribution for fitting discrete data:
#' Revival of the Conway-Maxwell-Poisson distribution,} J. Royal Statist. Soc.,
#' v54, pp. 127-142, 2005.
#' @keywords models
#' @examples
#' 
#' 	data(insurance);
#' 	fit = com.fit(Lemaire);
#' 	z = com.compute.z(fit$lambda, fit$nu);
#' 
#' @export com.compute.z
com.compute.z = function(lambda, nu, log.error = 0.001)
{
	return (exp(com.compute.log.z(lambda,nu,log.error)));
}

#' @importFrom matrixStats logSumExp
com.compute.log.z = function(lambda, nu, log.error = 0.001)
{
	# Perform argument checking
	if (lambda < 0 || nu < 0)
		stop("Invalid arguments, only defined for lambda >= 0, nu >= 0");
	
	# Initialize values
	z = -Inf;
	z.last = 0;
	j = 0;

	# Continue until we have reached specified precision
	while (abs(z - z.last) > log.error)
	{
		z.last = z;
		z = logSumExp(c(z, j * log(lambda) - nu * lfactorial(j)));

		j = j + 1;
	}
	return (z);
}