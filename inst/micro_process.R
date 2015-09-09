
## ----microprocess--------------------------------------------

micro_processed = plyr::llply(outputs, function(micro_output) {
  micro_stats = micro_output %>% 
    group_by(run, time) %>% 
    summarize(N = sum(population),
              P = sum(population*infections),
              mean_inf = P/N,
              var_inf = sum(population * ((infections - mean_inf)^2))/sum(population),
              var_mean_ratio = var_inf/mean_inf,
              log_vmr = log(var_mean_ratio)
    ) %>% 
    group_by()
  
  micro_stats2 = micro_output %>% 
    group_by(run, time) %>% 
    do({
      dat = cbind(.$infections, .$population)
   #   cmp = cmp_fit(dat)
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
#    mutate(cmp_AIC = 4 - 2*cmp_ll, pois_AIC = 2 - 2*pois_ll, nb_AIC = 4 - 2*nb_ll) %>% 
#    mutate(rel_AIC = pois_AIC - cmp_AIC) %>% 
    mutate(cmp_kld2 = pois_kld - cmp_kld, nb_kld2 = pois_kld - nb_kld, rel_kld = nb_kld - cmp_kld)
  
  micro_stats = micro_stats %>%
    inner_join(micro_stats2, by = c("run", "time")) %>% 
    group_by(run) %>% 
    filter(all(P[1:25] !=0)) %>% 
    group_by()
  
  micro_stats_tidy = micro_stats %>% 
    gather("variable", "value", -run, -time)
  
  micro_stats_mean = micro_stats_tidy %>%
    group_by(time, variable) %>% 
    summarise(mean=mean(value)) %>% 
    group_by() %>% 
    arrange(variable)
  
  return(list(
    micro_stats = micro_stats,
    micro_stats_tidy = micro_stats_tidy,
    micro_stats_mean = micro_stats_mean
  ))
    
}, .progress = ifelse(interactive(), "time", "none"))


