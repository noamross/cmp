master_odes <- function(t, y, parms) {
  #parms = relist(parm_vector)
  derivs = rep(NA, length(y))
  is = 0:(parms["max_i"])
  pop = sum(y)
  
  Lambda = parms["lambda"] * sum(is * y) + parms["lambda_ex"]
  derivs[1] = parms["r"]*pop*(1 - pop / parms["K"]) - y[1]*parms["d"] - y[1]*Lambda
  derivs[-c(1, length(y))] = Lambda*y[-c(length(y), length(y) - 1)] -
                               (Lambda + parms["d"])*y[-c(1, length(y))] - 
                               parms["alpha"]*((is[-c(1, length(y))])^parms["alpha_power"])*y[-c(1, length(y))] -
                               parms["mu"]*is[-c(1, length(y))]*y[-c(1, length(y))] +
                               parms["mu"]*is[-c(1, 2)]*y[-c(1, 2)]
  derivs[length(y)] = Lambda*y[length(y) - 1] - (Lambda + parms["d"])*y[length(y)] - 
                parms["alpha"]*((is[length(y-1)])^parms["alpha_power"])*y[length(y-1)] -
                parms["mu"]*is[length(y)]*y[length(y)]
  
  return(list(derivs))
}


dfe_init = function(parms, init = NULL) {
  if(is.null(init)) {
    init = c(parms$N0, rep(0, parms$max_i))
  } else {
    init = c(sum(init), rep(0, parms$max_i))
  }
  dfe = rootSolve::stode(y = init, time = 0, func = master_odes,
                         parms = unlist(within(parms, { lambda_ex = 0 })))
  return(dfe)
}


run_sodp = function(parms, init = NULL, matrix = FALSE) {
  if(is.null(init)) {
    init = c(parms$N0 - parms$P0, parms$P0, rep(0, parms$max_i - 1))
  }
  times = seq(0, parms$time_max, parms$step)
  parm_vector = unlist(parms)
  out = deSolve::lsoda(y = init, times = times, func = master_odes,
                       parms = parm_vector)
  if (matrix) return(out)
  out = data.frame(0, 0, out[,1], 0, out[,-1])
  names(out) = c("run", "start", "time", "control",
                     as.character(0:(ncol(out)-5)))
  return(dplyr::tbl_df(out))
}

#run_sodp <- R.cache::addMemoization(run_sodp.reg)
#dfe_init <- R.cache::addMemoization(dfe_init.reg)