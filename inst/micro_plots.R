## ----fig5--------------------------------------------

parsets = as.integer(c(1, 2, 3))
vars = c("N", "P", "var_mean_ratio")

micro_data = plyr::ldply(micro_processed[parsets], function(z) z$micro_stats_tidy) %>% 
    mutate(r = plyr::mapvalues(as.integer(.id), parsets, parms_grid$r[parsets])) %>% 
    filter(variable %in% vars) %>% 
    mutate(variable = plyr::revalue(variable, c("var_mean_ratio" = "Var/Mean Ratio"))) %>% 

    droplevels

micro_data_mean = plyr::ldply(micro_processed[parsets], function(z) z$micro_stats_mean) %>% 
    mutate(r = plyr::mapvalues(as.integer(.id), parsets, parms_grid$r[parsets])) %>% 
    filter(variable %in% vars) %>% 
    mutate(variable = plyr::revalue(variable, c("var_mean_ratio" = "Var/Mean Ratio"))) %>% 

    droplevels

ode_data = plyr::ldply(ode_processed[parsets], function(z) z$ode_stats_tidy) %>% 
  mutate(r = plyr::mapvalues(as.integer(.id), parsets, parms_grid$r[parsets])) %>% 
  filter(variable %in% vars) %>% 
  mutate(variable = plyr::revalue(variable, c("var_mean_ratio" = "Var/Mean Ratio"))) %>% 
  droplevels

Figure5 = ggplot() +
  geom_line(data=micro_data,
            mapping=aes(x = time, y = value, group=run),
            lwd = 1, alpha = 0.1, color = "slateblue") +
  geom_line(data=ode_data,
            mapping=aes(x = time, y = value),
            lwd = 1, color = "orange") +
  geom_line(data=micro_data_mean,
            mapping=aes(x = time, y = mean, col=variable),
            lwd = 1) +
  scale_color_manual(values = c("blue", "red", "darkgreen")) +

#  scale_color_manual(values = c("blue", "red", "darkgreen")) +
  facet_grid(variable ~ r, scales = "free") +
  xlab("Time") + ylab("") +
  theme_nr + 
  theme(legend.position = "none",
        axis.line.y = element_blank(),
        panel.grid.major.y = element_line(colour = "grey", size=1),
        panel.background = element_rect(),
        axis.text = element_text(size=12))

if(interactive()) print(Figure5)

## ----fig6--------------------------------------------

parsets = as.integer(c(1, 2, 3))
vars = c("rel_kld") #"cmp_kld2", "nb_kld2", 

micro_data = plyr::ldply(micro_processed[parsets], function(z) z$micro_stats_tidy) %>% 
    mutate(r = plyr::mapvalues(as.integer(.id), parsets, parms_grid$r[parsets])) %>% 
    filter(variable %in% vars) %>% 
    mutate(variable = plyr::revalue(variable, c("var_mean_ratio" = "Var/Mean Ratio"))) %>% 

    droplevels

micro_data_mean = plyr::ldply(micro_processed[parsets], function(z) z$micro_stats_mean) %>% 
    mutate(r = plyr::mapvalues(as.integer(.id), parsets, parms_grid$r[parsets])) %>% 
    filter(variable %in% vars) %>% 
    mutate(variable = plyr::revalue(variable, c("var_mean_ratio" = "Var/Mean Ratio"))) %>% 

    droplevels

ode_data = plyr::ldply(ode_processed[parsets], function(z) z$ode_stats_tidy) %>% 
  mutate(r = plyr::mapvalues(as.integer(.id), parsets, parms_grid$r[parsets])) %>% 
  filter(variable %in% vars) %>% 
  mutate(variable = plyr::revalue(variable, c("var_mean_ratio" = "Var/Mean Ratio"))) %>% 
  droplevels

Figure6 = ggplot() +
  geom_line(data=micro_data,
            mapping=aes(x = time, y = value, group=run),
            lwd = 1, alpha = 0.1, color = "slateblue") +
  geom_line(data=ode_data,
            mapping=aes(x = time, y = value),
            lwd = 1, color = "orange") +
  geom_line(data=micro_data_mean,
            mapping=aes(x = time, y = mean, col=variable),
            lwd = 1) +
  scale_color_manual(values = c("blue", "red", "darkgreen")) +

#  scale_color_manual(values = c("blue", "red", "darkgreen")) +
  facet_grid(. ~ r, scales = "free") +
  xlab("Time") + ylab("Relative KLD (Neg. Binomial - CMP)") +
  theme_nr + 
  theme(legend.position = "none",
        axis.line.y = element_blank(),
        panel.grid.major.y = element_line(colour = "grey", size=1),
        panel.background = element_rect(),
        axis.text = element_text(size=12))

if(interactive()) print(Figure6)

## ----fig7--------------------------------------------

parsets = as.integer(c(4, 5, 6))
vars = c("N", "P", "var_mean_ratio")

micro_data = plyr::ldply(micro_processed[parsets], function(z) z$micro_stats_tidy) %>% 
    mutate(r = plyr::mapvalues(as.integer(.id), parsets, parms_grid$r[parsets])) %>% 
    filter(variable %in% vars) %>% 
    mutate(variable = plyr::revalue(variable, c("var_mean_ratio" = "Var/Mean Ratio"))) %>% 

    droplevels

micro_data_mean = plyr::ldply(micro_processed[parsets], function(z) z$micro_stats_mean) %>% 
    mutate(r = plyr::mapvalues(as.integer(.id), parsets, parms_grid$r[parsets])) %>% 
    filter(variable %in% vars) %>% 
    mutate(variable = plyr::revalue(variable, c("var_mean_ratio" = "Var/Mean Ratio"))) %>% 

    droplevels

ode_data = plyr::ldply(ode_processed[parsets], function(z) z$ode_stats_tidy) %>% 
  mutate(r = plyr::mapvalues(as.integer(.id), parsets, parms_grid$r[parsets])) %>% 
  filter(variable %in% vars) %>% 
  mutate(variable = plyr::revalue(variable, c("var_mean_ratio" = "Var/Mean Ratio"))) %>% 
  droplevels

Figure7 = ggplot() +
  geom_line(data=micro_data,
            mapping=aes(x = time, y = value, group=run),
            lwd = 1, alpha = 0.1, color = "slateblue") +
  geom_line(data=ode_data,
            mapping=aes(x = time, y = value),
            lwd = 1, color = "orange") +
  geom_line(data=micro_data_mean,
            mapping=aes(x = time, y = mean, col=variable),
            lwd = 1) +
  scale_color_manual(values = c("blue", "red", "darkgreen")) +

#  scale_color_manual(values = c("blue", "red", "darkgreen")) +
  facet_grid(variable ~ r, scales = "free") +
  xlab("Time") + ylab("") +
  theme_nr + 
  theme(legend.position = "none",
        axis.line.y = element_blank(),
        panel.grid.major.y = element_line(colour = "grey", size=1),
        panel.background = element_rect(),
        axis.text = element_text(size=12))

if(interactive()) print(Figure7)

## ----fig8--------------------------------------------

parsets = as.integer(c(4, 5, 6))
vars = c("rel_kld") #"cmp_kld2", "nb_kld2", 

micro_data = plyr::ldply(micro_processed[parsets], function(z) z$micro_stats_tidy) %>% 
    mutate(r = plyr::mapvalues(as.integer(.id), parsets, parms_grid$r[parsets])) %>% 
    filter(variable %in% vars) %>% 
    mutate(variable = plyr::revalue(variable, c("var_mean_ratio" = "Var/Mean Ratio"))) %>% 

    droplevels

micro_data_mean = plyr::ldply(micro_processed[parsets], function(z) z$micro_stats_mean) %>% 
    mutate(r = plyr::mapvalues(as.integer(.id), parsets, parms_grid$r[parsets])) %>% 
    filter(variable %in% vars) %>% 
    mutate(variable = plyr::revalue(variable, c("var_mean_ratio" = "Var/Mean Ratio"))) %>% 

    droplevels

ode_data = plyr::ldply(ode_processed[parsets], function(z) z$ode_stats_tidy) %>% 
  mutate(r = plyr::mapvalues(as.integer(.id), parsets, parms_grid$r[parsets])) %>% 
  filter(variable %in% vars) %>% 
  mutate(variable = plyr::revalue(variable, c("var_mean_ratio" = "Var/Mean Ratio"))) %>% 
  droplevels

Figure8 = ggplot() +
  geom_line(data=micro_data,
            mapping=aes(x = time, y = value, group=run),
            lwd = 1, alpha = 0.1, color = "slateblue") +
  geom_line(data=ode_data,
            mapping=aes(x = time, y = value),
            lwd = 1, color = "orange") +
  geom_line(data=micro_data_mean,
            mapping=aes(x = time, y = mean, col=variable),
            lwd = 1) +
  scale_color_manual(values = c("blue", "red", "darkgreen")) +

#  scale_color_manual(values = c("blue", "red", "darkgreen")) +
  facet_grid(. ~ r, scales = "free") +
  xlab("Time") + ylab("Relative KLD (Neg. Binomial - CMP)") +
  theme_nr + 
  theme(legend.position = "none",
        axis.line.y = element_blank(),
        panel.grid.major.y = element_line(colour = "grey", size=1),
        panel.background = element_rect(),
        axis.text = element_text(size=12))

if(interactive()) print(Figure8)

## ----fig9--------------------------------------------
# 
# ode_stats = ode_processed[[1]]$ode_stats
# index = which.min(ode_stats$var_mean_ratio)
# index_time = ode_stats$time[index]
# sim_output = data.frame(infections=0:baseparms$max_i, cmp=ode_stats$N[index] *
#     dcmp(0:baseparms$max_i, ode_stats$cmp_kld_lambda[index],
#          ode_stats$cmp_kld_nu[index]),
#     nb=ode_stats$N[index] *  dnbinom(0:baseparms$max_i, mu=ode_stats$nb_kld_my[index], size=ode_stats$nb_kld_size[index]))
# 
# plotruns = micro_processed[[1]]$micro_stats %>% group_by(run) %>% summarize(ext = all(P[1:25] != 0))
# plotruns = plotruns$run[plotruns$ext]
# micro_output_plot = outputs[[1]] %>% filter(run %in% plotruns)
# 
# micro_output_plot_mean = micro_output_plot %>% 
#   group_by(time, infections) %>% 
#   summarize(population = mean(population))
# 
# ggplot() +
#   geom_line(data = filter(micro_output_plot, time %in% index_time),
#             mapping = aes(x=infections, y=population, group=run),
#             col="slateblue", alpha=0.24, lwd=1) +
#   geom_line(data = filter(micro_output_plot_mean, time %in% index_time),
#             mapping = aes(x=infections, y=population),
#             col="green", lwd=1) +
#   geom_line(data = filter(ode_processed[[1]]$ode_output_tidy, time %in% index_time,
#                           infections <= max(micro_output_plot$infections)),
#             mapping = aes(x=infections, y=population, group=run),
#             col="red", lwd=1, lty=2) +
#   geom_line(data = filter(sim_output, infections <= max(micro_output_plot$infections)),
#             mapping = aes(x=infections, y = cmp), col="orange", lwd=1) +
#   #  geom_line(data = filter(sim_output, infections <= max(micro_output_plot$infections)),
#    #         mapping = aes(x=infections, y = nb), col="brown", lwd=1) +
#   scale_y_continuous(trans=log1p_trans(), oob=rescale_none) + theme_nr 
# 
# 
