#Renamed functions provided for the sake of backwards compatibility

#Note Rcpp cannot export functions with "." in them, 
#' @rdname com_compute_z
#' @export
com.compute.z = com_compute_z

#' @rdname com_compute_z
#' @export
com.compute.log.z = com_compute_log_z

#' @rdname dcom
#' @export
com.log.density = function(x, lambda, nu, log.z = NULL) {
  return(dcom(x = x, lambda = lambda, nu = nu, log = TRUE,
              z = ifelse(!is.null(log.z), exp(z), NA_real_)))
}
