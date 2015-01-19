#' Operations in Log-space
#'
#' Computes the difference of two values in log-space.
#'
#' \code{com.log.difference} computes the difference of two values in
#' log-space, \eqn{log( e^x - e^y )}, without significant chance of overflow or
#' underflow.
#'
#' @param x first value
#' @param y second value
#' @return The requested computation in log-space.
#' @author Jeffrey Dunn
#' @keywords manip
#' @examples
#'
#'  a = exp(com.log.difference(log(100), log(20))) # a = 80
#'  b = exp(logsumexp(log(100), log(20))) # b = 120
#'  c = exp(lfactorial(4)) # c = 24
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
