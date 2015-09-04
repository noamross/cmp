#This chunk sets the base parameters for all of the simulations.

## ----loadpkgs------------------------------------------------------------

library(spore)
library(cmp)
library(dplyr)
library(tidyr)
library(magrittr)
library(ggplot2)
library(scales)
library(zoo)
library(spore)
library(readr)
library(deSolve)
library(parallel)
library(noamtools)
library(cowplot)
library(nloptr)
library(R.cache)
library(extrafont)

## ----setparms------------------------------------------------------------
parms = list(
  max_i = 100,
  is = 0:100,
  lambda = 0.0004,
  lambda_ex = 0.01,
  alpha = 0.088066,
  alpha_power = 1,
  mu =   0.1,
  r = 0.1,
  d = 0.0,
  K = 1000,
  time_max = 100,
  prevent_inf = 0,
  prevent_ex = 0,
  macro_timestep = 1,
  micro_timestep = 0.1,
  micro_relax_steps = 0,
  delta = 0,
  project = FALSE,
  n_comp_sims = 10,
  n_sims = 1000,
  n_sims_jacob = 10000000,
  control_min = 0,
  control_max = 10,
  v = 4,
  c = 400,
  nocontrol = TRUE,
  progress = TRUE,
  micro_record = 0,
  parallel_cores = 2,
  step = 1,
  N0 = 1000,
  P0 = 1
  
  #micro_record = file("micro.txt", open="w+")
  #  macro_record = file("macro.txt", open="w")
)

source('inst/run_sodp.R')

# ode_output = run_sodp(parms)
ode_output = run_sodp(parms, init = parms$N0*dpois(0:parms$max_i, parms$P0/parms$N0))

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
    cmp = cmp_fit(dat)
    cmp_k = cmp_fit_kld(dat)
    pois = pois_fit(dat)
    pois_k = pois_fit_kld(dat)
    nb = nb_fit(dat)
    nb_k = nb_fit_kld(dat)
    data_frame(cmp_lambda = cmp$lambda, cmp_nu = cmp$nu, cmp_ll = cmp$log.likelihood,
               pois_lambda = pois$lambda, pois_ll = pois$log.likelihood,
               nb_mu = nb$mu, nb_size = nb$size, nb_ll = nb$log.likelihood,
               cmp_kld_lambda = cmp_k$lambda, cmp_kld_nu = cmp_k$nu, cmp_kld = cmp_k$kld,
               pois_kld_lambda = pois_k$lambda, pois_kld = pois_k$kld,
               nb_kld_my = nb_k$mu, nb_kld_size = nb_k$size, nb_kld = nb_k$kld)
  }) %>%
  mutate(cmp_AIC = 4 - 2*cmp_ll, pois_AIC = 2 - 2*pois_ll, nb_AIC = 4 - 2*nb_ll) %>% 
  mutate(rel_AIC = pois_AIC - cmp_AIC) %>% 
  mutate(cmp_kld2 = pois_kld - cmp_kld, nb_kld2 = pois_kld - nb_kld)


ode_stats_tidy = ode_stats %>% 
  inner_join(ode_stats2, by = c("run", "time")) %>% 
  gather("variable", "value", -run, -time)

ggplot() +
  geom_line(data = ode_stats_tidy, mapping = aes(x=time, y=value, group=run), 
            col="red", lwd=1) +
  facet_wrap(~variable, scales="free") +
  theme_nr

ggplot() +
  geom_line(data = filter(ode_stats_tidy, variable %in% c("cmp_kld2", "nb_kld2")),
                          mapping = aes(x=time, y=value, col=variable), lwd=1) + theme_nr

micro_state = c(parms$N0 - parms$P0, parms$P0, rep(0, parms$max_i - 1))
macro_state = restrict.micro_state(micro_state)

closeAllConnections()
parms$micro_record = file("micro.txt", open="w+")
for(run in 1:parms$n_comp_sims) {
  micro_state_c.stepto(micro_state, parms = parms, time = 0,
                       timeto=parms$time_max, record = parms$micro_record,
                       run = run, control = 0)
  if(interactive()) cat(run, "\r")
}
closeAllConnections()

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

micro_output = micro_output_ints %>% 
  arrange(run, time) %>%
  gather("infections", "population", -run, -start, -time, -control) %>%
  arrange(run, time, infections) %>%
  mutate(run=as.factor(run), infections = as.integer(as.character(infections))) 

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
    cmp = cmp_fit(dat)
    cmp_k = cmp_fit_kld(dat)
    pois = pois_fit(dat)
    pois_k = pois_fit_kld(dat)
    nb = nb_fit(dat)
    nb_k = nb_fit_kld(dat)
    data_frame(cmp_lambda = cmp$lambda, cmp_nu = cmp$nu, cmp_ll = cmp$log.likelihood,
               pois_lambda = pois$lambda, pois_ll = pois$log.likelihood,
               nb_mu = nb$mu, nb_size = nb$size, nb_ll = nb$log.likelihood,
               cmp_kld_lambda = cmp_k$lambda, cmp_kld_nu = cmp_k$nu, cmp_kld = cmp_k$kld,
               pois_kld_lambda = pois_k$lambda, pois_kld = pois_k$kld,
               nb_kld_my = nb_k$mu, nb_kld_size = nb_k$size, nb_kld = nb_k$kld)
  }) %>%
  mutate(cmp_AIC = 4 - 2*cmp_ll, pois_AIC = 2 - 2*pois_ll, nb_AIC = 4 - 2*nb_ll) %>% 
  mutate(rel_AIC = pois_AIC - cmp_AIC) %>% 
  mutate(cmp_kld2 = pois_kld - cmp_kld, nb_kld2 = pois_kld - nb_kld)


micro_stats_tidy = micro_stats %>%
  inner_join(micro_stats2, by = c("run", "time")) %>% 
  gather("variable", "value", -run, -time)

micro_stats_mean = micro_stats_tidy %>% 
  group_by(time, variable) %>% 
  summarise(mean=mean(value)) %>% 
  group_by() %>% 
  arrange(variable)

ggplot() +
  geom_line(data = micro_stats_tidy, mapping = aes(x=time, y=value, group=run), 
            col="slateblue", alpha=0.5, lwd=1) +
  geom_line(data = ode_stats_tidy, mapping = aes(x=time, y=value, group=run), 
            col="red", lwd=1) +
  geom_line(data = micro_stats_mean, mapping = aes(x=time, y=mean), col="green", lwd=1) + 
  facet_wrap(~variable, scales="free") +
  theme_nr

ggplot() +
  geom_line(data = filter(micro_stats_tidy, variable %in% c("cmp_kld2", "nb_kld2")),
            mapping = aes(x=time, y=value, group=run), 
            col="slateblue", alpha=0.5, lwd=1) +
  geom_line(data = filter(ode_stats_tidy, variable %in% c("cmp_kld2", "nb_kld2")),
            mapping = aes(x=time, y=value, group=run), 
            col="red", lwd=1) +
  geom_line(data = filter(micro_stats_mean, variable %in% c("cmp_kld2", "nb_kld2")),
            mapping = aes(x=time, y=mean), col="green", lwd=1) + 
  theme_nr


ggplot() +
  geom_line(data = filter(ode_stats_tidy, variable %in% c("cmp_kld2", "nb_kld2")),
                          mapping = aes(x=time, y=value, col=variable), lwd=1) + 
  
  theme_nr


plot_times = round(seq(0, parms$time_max, length.out = 12))
ggplot() +
  geom_line(data = filter(micro_output, time %in% plot_times),
            mapping = aes(x=infections, y=population, group=run),
            col="slateblue", alpha=0.24, lwd=1) +
  geom_line(data = filter(ode_output_tidy, time %in% plot_times,
                          infections <= max(micro_output$infections)),
            mapping = aes(x=infections, y=population, group=run),
            col="red", lwd=1, lty=3) +
  scale_y_log10(ylim(1, 200), expand=c(0,0), breaks=c(1,10,50,200)) +
  facet_wrap(~time, nrow=2)
