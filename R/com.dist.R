#' The COM-Poisson Distribution
#' 
#' Probability mass function and random generation for the COM-Poisson
#' distribution for given values of the parameters.
#' 
#' Computes the probability mass function of the COM-Poisson distribution
#' \deqn{ }{f(x) = (1/Z) (lambda^x)/(x!^nu).}\deqn{ f(x) =
#' \frac{1}{Z(\lambda,\nu)}\frac{\lambda^x}{(x!)^\nu}. }{f(x) = (1/Z)
#' (lambda^x)/(x!^nu).}\deqn{ }{f(x) = (1/Z) (lambda^x)/(x!^nu).}
#' 
#' @aliases dcom rcom
#' @param x level to evaluate the PMF at
#' @param lambda value of lambda parameter
#' @param nu value of nu parameter
#' @param z normalizing constant, computed if not specified
#' @param n number of random values to return
#' @param log.z natural log of z
#' @return \code{dcom} gives the probability that a random COM-Poisson variable
#' X takes value x.  \code{rcom} gives a vector of \code{n} random values
#' sampled from the COM-Poisson distribution.
#' @author Jeffrey Dunn
#' @seealso \code{\link{com.loglikelihood}}, \code{\link{com.log.density}}
#' @references Shmueli, G., Minka, T. P., Kadane, J. B., Borle, S. and
#' Boatwright, P., \dQuote{A useful distribution for fitting discrete data:
#' Revival of the Conway-Maxwell-Poisson distribution,} J. Royal Statist. Soc.,
#' v54, pp. 127-142, 2005.
#' @keywords models
#' @examples
#' 
#' 	data(insurance);
#' 	fit = com.fit(Lemaire);
#' 	dcom(0, fit$lambda, fit$nu, fit$z);
#' 	r = rcom(10, fit$lambda, fit$nu);
#' 
#' @export dcom
dcom = function(x, lambda, nu, z = NULL)
{
	# Perform argument checking
	if (lambda < 0 || nu < 0)
		stop("Invalid arguments, only defined for lambda >= 0, nu >= 0");
	if (x < 0 || x != floor(x))
		return (0);
	if (is.null(z) || z <= 0)
		z = com.compute.z(lambda, nu);
	
	# Return pmf
	return ((lambda ^ x) * ((factorial(x)) ^ -nu) / z);
}

#' @export
rcom = function(n, lambda, nu, log.z = NULL)
{
  # Check arguments
	if (lambda < 0 || nu < 0)
		stop("Invalid arguments, only defined for lambda >= 0, nu >= 0");
	if (is.null(log.z))
		log.z = com.compute.log.z(lambda, nu);

	r = NULL;	# Vector of random values

	for (i in 1:n)
	{
		# Get a uniform random variable and find the smallest value x such that
		# the cdf of x is greater than the random value
		log.prob = log(runif(1));
		j = 0;
		while (1)
		{
			new.log.prob = com.log.difference( log.prob, com.log.density(j, lambda, nu, log.z) );
			if (is.nan(new.log.prob))
				break;

			log.prob = new.log.prob;
			j = j + 1;
		}

		r = c(r, j);
	}

	return (r);
}


#' Computes the Log PMF of the COM-Poisson Distribution
#' 
#' Computes the log probability mass function of the COM-Poisson distribution
#' for given values of the parameters.
#' 
#' Computes the log probability mass function of the COM-Poisson distribution
#' \deqn{ }{log f(x) = x * log(lambda) - log(Z) - nu * log(x!). }\deqn{ \log
#' f(x) = x \log \lambda - \log(Z(\lambda,\nu)) - \nu \sum_{i=1}^x x. }{log
#' f(x) = x * log(lambda) - log(Z) - nu * log(x!). }\deqn{ }{log f(x) = x *
#' log(lambda) - log(Z) - nu * log(x!). }
#' 
#' @param x level to evaulate the log PMF at
#' @param lambda value of the lambda parameter
#' @param nu value of the nu parameter
#' @param log.z log of the normalizing constant, computed if not specified
#' @return The log probability that a random COM-Poisson variable X takes value
#' x.
#' @author Jeffrey Dunn
#' @seealso \code{\link{com.loglikelihood}}, \code{\link{dcom}}
#' @references Shmueli, G., Minka, T. P., Kadane, J. B., Borle, S. and
#' Boatwright, P., \dQuote{A useful distribution for fitting discrete data:
#' Revival of the Conway-Maxwell-Poisson distribution,} J. Royal Statist. Soc.,
#' v54, pp. 127-142, 2005.
#' @keywords models
#' @examples
#' 
#' 	data(insurance);
#' 	fit = com.fit(Lemaire);
#' 	com.log.density(0, fit$lambda, fit$nu, fit$z);
#' 
#' @export com.log.density
com.log.density = function(x, lambda, nu, log.z = NULL)
{
	# Perform argument checking
	if (lambda < 0 || nu < 0)
		stop("Invalid arguments, only defined for lambda >= 0, nu >= 0");
	if (x < 0 || x != floor(x))
		return (0);
	if (is.null(log.z)) { log.z = com.compute.log.z(lambda, nu); }
	
	# Return log pmf
	return ((x * log(lambda) - nu * com.log.factorial(x)) - log.z);
}