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
#'  @export
com.log.difference = function(x,y) { # log.difference(x,y) = log(exp(x)- xp(y))
	if (x == -Inf) {
		return (NaN)
	} else if (y == -Inf) {
		return (x)
	} else if (x > y) {
	  return (x + log( 1 - exp(y - x) ) )
	} else {
		return (NaN)
	}
}

#' @export
com.log.factorial = base::lfactorial
  
#   function(x)	# log(factorial(x))
# {
# 	if (is.vector(x) && length(x) > 1)
# 	{
# 		for (i in 1:length(x))
# 			x[i] = com.log.factorial(x[i]);
# 		return (x);
# 	}
# 	else if (is.numeric(x))
# 	{
# 		if (x == 0) { x = 1; }
# 		return (sum(log(seq(from = 1, to = x, by = 1))));
# 	}
# 	else { stop("x must be a vector or number."); }
# }

#' @export
com.log.sum = function(x,y)  	# log.sum(x,y) = log( exp(x) + exp(y) )
  
{
	if (x == -Inf)
		{ return (y); }
	else if (y == -Inf)
		{ return (x); }
	else if (x > y)
		{ return (x + log( 1 + exp(y - x) ) ); }
	else
		{ return (y + log( 1 + exp(x - y) ) ); }
}

