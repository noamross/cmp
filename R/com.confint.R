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

