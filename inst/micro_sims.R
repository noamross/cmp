## ----microsims--------------------------------------------


file.create("micro.txt", showWarnings = FALSE)

outputs = plyr::alply(parms_grid, 1, function(chg_parms) {
  
  parms = within(baseparms, {
    r = chg_parms$r
    alpha_power = chg_parms$alpha_power
  })
  
  micro_state = c(parms$N0 - parms$P0, parms$P0, rep(0, parms$max_i - 1))
  
  parms$micro_record = file("micro.txt", open="w+")
  
  for(run in 1:parms$n_comp_sims) {
    micro_state_c.stepto(micro_state, parms = parms, time = 0,
                         timeto=parms$time_max, record = parms$micro_record,
                         run = run, control = 0)
    if(interactive()) cat(run, "\r")
  }
  
  close(parms$micro_record); parms$micro_record = NULL
  
  output = read_delim("micro.txt", " ", col_names = FALSE)
  
  colnames(output) = c("run", "start", "time", "control",
                       as.character(0:(ncol(output)-5)))
  used_inf_classes = as.logical(sapply(output[-(1:4)], sum))
  used_inf_classes[is.na(used_inf_classes)] = FALSE
  output = output[, c(TRUE, TRUE, TRUE, TRUE, used_inf_classes)]
  
  micro_output_ints = output %>%
    group_by(run) %>%
    filter(!duplicated(plyr::round_any(time, parms$step, floor), fromLast=TRUE) |
             time == 0) %>%
    group_by() %>%
    mutate(time = plyr::round_any(time, parms$step, f=ceiling)) %>%
    mutate_each(funs(na_to_0), -run, -start, -time, -control) %>%
    mutate(run = factor(run)) %>% 
    filter(time <= parms$time_max) %>% 
    group_by() 
  
  micro_output_full = micro_output_ints %>% 
    full_join(expand.grid(time=0:parms$time_max, run=as.factor(1:parms$n_comp_sims))) %>% 
    arrange(run, time) %>% 
    group_by(run) %>% 
    mutate_each(funs(na.locf), -run, -time)
  
  micro_output = micro_output_full %>% 
    arrange(run, time) %>%
    gather("infections", "population", -run, -start, -time, -control) %>%
    arrange(run, time, infections) %>%
    mutate(run=as.factor(run), infections = as.integer(as.character(infections))) 
  
  return(micro_output)
  
}, .progress = ifelse(interactive(), "time", "none"))
