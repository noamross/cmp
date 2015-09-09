## ----odesims------------------------------------------------------------

ode_outputs = plyr::alply(parms_grid, 1, function(chg_parms) {
  parms = within(baseparms, {
    r = chg_parms$r
    alpha_power = chg_parms$alpha_power
  })
  ode_output = run_sodp(parms,
                        init = parms$N0*dpois(0:parms$max_i, parms$P0/parms$N0))
})


## ----odeprocess--------------------------------------------

ode_processed = plyr::llply(ode_outputs, function(ode_output) {
  
  ode_output_tidy = ode_output %>% 
    arrange(run, time) %>% 
    select(-control, -start) %>% 
    gather("infections", "population", -run, -time) %>% 
    arrange(run, time, infections) %>% 
    mutate(run=as.factor(run), infections = as.integer(as.character(infections))) 
  
  ode_stats = ode_output_tidy %>% 
    group_by(run, time) %>% 
    summarize(N = sum(population),
              P = sum(population*infections),
              mean_inf = P/N,
              var_inf = sum(population * ((infections - mean_inf)^2))/sum(population),
              var_mean_ratio = var_inf/mean_inf,
              log_vmr = log(var_mean_ratio)
    ) %>% 
    group_by()
  
  ode_stats2 = ode_output_tidy %>% 
    group_by(run, time) %>% 
    do({
      dat = cbind(.$infections, .$population)
     # cmp = cmp_fit(dat)
      cmp_k = cmp_fit_kld(dat)
    #  pois = pois_fit(dat)
      pois_k = pois_fit_kld(dat)
     # nb = nb_fit(dat)
      nb_k = nb_fit_kld(dat)
      data_frame(#cmp_lambda = cmp$lambda, cmp_nu = cmp$nu, cmp_ll = cmp$log.likelihood,
                 #pois_lambda = pois$lambda, pois_ll = pois$log.likelihood,
                 #nb_mu = nb$mu, nb_size = nb$size, nb_ll = nb$log.likelihood,
                 cmp_kld_lambda = cmp_k$lambda, cmp_kld_nu = cmp_k$nu, cmp_kld = cmp_k$kld,
                 pois_kld_lambda = pois_k$lambda, pois_kld = pois_k$kld,
                 nb_kld_my = nb_k$mu, nb_kld_size = nb_k$size, nb_kld = nb_k$kld)
    }) %>%
 #   mutate(cmp_AIC = 4 - 2*cmp_ll, pois_AIC = 2 - 2*pois_ll, nb_AIC = 4 - 2*nb_ll) %>% 
 #   mutate(rel_AIC = pois_AIC - cmp_AIC) %>% 
    mutate(cmp_kld2 = pois_kld - cmp_kld, nb_kld2 = pois_kld - nb_kld, rel_kld = nb_kld - cmp_kld)
  
  ode_stats = ode_stats %>% 
    inner_join(ode_stats2, by = c("run", "time"))
  
  ode_stats_tidy = ode_stats %>% 
    gather("variable", "value", -run, -time)
  
  return(list(
    ode_output_tidy = ode_output_tidy,
    ode_stats = ode_stats,
    ode_stats_tidy = ode_stats_tidy
  ))

}, .progress = ifelse(interactive(), "time", "none"))

