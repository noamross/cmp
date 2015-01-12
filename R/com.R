#' Computes COM-Poisson Regression
#' 
#' Computes the maximum likelihood estimates of the COM-Poisson model for given
#' count data.
#' 
#' The argument x should consist of a matrix where the first column is the
#' level and the second column is the count for the corresponding level.
#' 
#' @param x matrix of count data
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
com.fit = function(x)
{
	xbar = (x[,1] %*% x[,2]) / sum(x[,2]);
	options(warn = -1);
	result = optim(c(xbar,1), function(p) {return (-com.loglikelihood(x, p[1], p[2]));},
		method="L-BFGS-B", lower=c(1e-10,0));
	options(warn = 0);
	
	lambda = result$par[1];
	nu = result$par[2];
	fit = list( lambda = lambda,
	            nu = nu,
	            z = com.compute.z(lambda, nu),
	            fitted.values = sum(x[,2]) * dcom(x[,1], lambda, nu),
				log.likelihood = com.loglikelihood(x, lambda, nu) );

	return (fit);
}

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
		z = com.log.sum(z, j * log(lambda) - nu * com.log.factorial(j));

		j = j + 1;
	}
	return (z);
}





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
com.loglikelihood = function(x, lambda, nu)
{
	# Perform argument checking
	if (lambda < 0 || nu < 0)
		return (-Inf);

	log.z = com.compute.log.z(lambda, nu);
	return (x[,2] %*% ( x[,1] * log(lambda) - nu * com.log.factorial(x[,1]) - log.z ));
}





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
#' @export com.expectation
com.expectation = function(f, lambda, nu, log.error = 0.001)
{
	log.z = com.compute.log.z(lambda, nu);

	# Initialize variables
	ex = -Inf;
	ex.last = 0;
	j = 0;

	# Continue until we have reached specified precision
	while ((ex == -Inf && ex.last == -Inf) || abs(ex - ex.last) > log.error)
	{
		ex.last = ex;
		ex = com.log.sum(ex, log(f(j)) + com.log.density(j, lambda, nu, log.z));

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






#' Computes a confidence interval for parameter estimates of the COM-Poisson
#' Distribution
#' 
#' Computes a pivotal bootstrap confidence interval for maximum likelihood
#' parameter estimates.
#' 
#' Uses a standard pivotal confidence interval from a bootstrap sample.
#' 
#' @param data the matrix of data to fit
#' @param level the level of the confidence interval
#' @param B number of repetitions of the bootstrap
#' @param n number of data points in each bootstrap sample
#' @return A matrix containing the confidence intervals for each parameter
#' @author Akshaya Jha, Jeffrey Dunn
#' @seealso \code{\link{com.fit}}
#' @references Wasserman, L. (2005). "All of Statistics: A Concise Course in
#' Statistical Inference," Springer Texts in Statistics.
#' @keywords models
#' @export com.confint
com.confint = function(data, level=0.95, B=1000, n=1000)
{
	# B = Number of repetitions of bootstrap
	# n = number of obs. in each sample used to find bootstrap mean
	# data = the original data
	# level - confidence interval level


	# Check arguments
	if (level <= 0 || level >= 1)
		stop("Invalid arguments, 0 < level < 1");
	if (n <= 0)
		stop("Invalid arguments, n > 0");
	if (B <= 0)
		stop("Invalid arguments, B > 0");


	boot.lambda = matrix(0,B,1)
	boot.nu = matrix(0,B,1)

	# Sample statistic
	COMobject = com.fit(data)
	lambda.mle = COMobject$lambda 
	nu.mle = COMobject$nu

	# Changing data into form we can sample from
	fulldata = matrix(0,sum(data[,2]),1)
	index = 0
	for (i in 1:length(data[,1]))
	{
		index2 = index + data[i,2]
		fulldata[(index+1):index2] = data[i,1]
	}

	# Creates a vector of means (the bootstrap vector)
	for (i in 1:B)
	{
		samplewR = sample(fulldata, n, replace=TRUE);
		sample = data.matrix( as.data.frame(table(samplewR)) );
		sample = cbind( sample[,1] - 1, sample[,2] );
		COMsampleobject = com.fit(sample);
		boot.lambda[i] = COMsampleobject$lambda;
		boot.nu[i] = COMsampleobject$nu;
		remove(samplewR);
	}

	# Pivotal method calculation
	boot.lambda = sort(boot.lambda)
	boot.nu = sort(boot.nu)

	lower.bound = (1 - level) / 2;
	upper.bound = 1 - lower.bound;

	lower.index = floor(lower.bound * B)+1;
	upper.index = floor(upper.bound * B);

	pivotal.ci = matrix(0,2,2);
	pivotal.ci[1,1] = max(0, 2*lambda.mle - boot.lambda[upper.index]);
	pivotal.ci[1,2] = 2*lambda.mle - boot.lambda[lower.index];

	pivotal.ci[2,1] = max(0, 2*nu.mle - boot.nu[upper.index]);
	pivotal.ci[2,2] = 2*nu.mle - boot.nu[lower.index];

	rownames(pivotal.ci) = c("lambda", "nu");
	colnames(pivotal.ci) = c(lower.bound, upper.bound);

	return (pivotal.ci);
}

