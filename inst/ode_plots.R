## ----fig1------------------------------------------------------------
# equal_label = function (variable, value) paste(variable, value, sep = " = ") 

theme_nr = theme_nr + theme(axis.text.y = element_text(size=11))

parsets = as.integer(c(1, 2, 3))
vars = c("N", "P", "var_mean_ratio")

data = plyr::ldply(ode_processed[parsets], function(z) z$ode_stats_tidy) %>% 
  mutate(r = plyr::mapvalues(as.integer(.id), parsets, parms_grid$r[parsets])) %>% 
  filter(variable %in% vars) %>% 
  mutate(variable = plyr::revalue(variable, c("var_mean_ratio" = "Var/Mean Ratio"))) %>% 
  droplevels

Figure1 = ggplot(data, aes(x = time, y = value, col = variable)) +
  scale_color_manual(values = c("blue", "red", "darkgreen")) +
  geom_line(lwd = 1) +
  facet_grid(variable ~ r, scales = "free") +
  xlab("Time") + ylab("") +
  theme_nr + 
  theme(legend.position = "none",
        axis.line.y = element_blank(),
        panel.grid.major.y = element_line(colour = "grey", size=1),
        panel.background = element_rect(),
        axis.text = element_text(size=12))

if(interactive()) print(Figure1)

## ----fig2--------------------------------------------

parsets = as.integer(c(1, 2, 3))
vars = c("rel_kld") #"cmp_kld2", "nb_kld2", 

data = plyr::ldply(ode_processed[parsets], function(z) z$ode_stats_tidy) %>% 
  mutate(r = plyr::mapvalues(as.integer(.id), parsets, parms_grid$r[parsets])) %>% 
  filter(variable %in% vars) %>% 
  mutate(variable = plyr::revalue(variable, c("var_mean_ratio" = "Var/Mean Ratio"))) %>% 
  droplevels

Figure2 = ggplot(data, aes(x = time, y = value, col = variable)) +
  scale_color_manual(values = c("blue", "red", "darkgreen")) +
  geom_line(lwd = 1) +
  facet_grid(. ~ r, scales = "free") +
  xlab("Time") + ylab("Relative KLD (Neg. Binomial - CMP)") +
  theme_nr + 
  theme(legend.position = "none",
        axis.line.y = element_blank(),
        panel.grid.major.y = element_line(colour = "grey", size=1),
        panel.background = element_rect(),
        axis.text = element_text(size=12))

if(interactive()) print(Figure2)

## ----fig3------------------------------------------------------------

parsets = as.integer(c(4, 5, 6))
vars = c("N", "P", "var_mean_ratio")

data = plyr::ldply(ode_processed[parsets], function(z) z$ode_stats_tidy) %>% 
  mutate(r = plyr::mapvalues(as.integer(.id), parsets, parms_grid$r[parsets])) %>% 
  filter(variable %in% vars) %>% 
  mutate(variable = plyr::revalue(variable, c("var_mean_ratio" = "Var/Mean Ratio"))) %>% 
  droplevels

Figure3 = ggplot(data, aes(x = time, y = value, col = variable)) +
  scale_color_manual(values = c("blue", "red", "darkgreen")) +
  geom_line(lwd = 1) +
  facet_grid(variable ~ r, scales = "free") +
  xlab("Time") + ylab("") +
  theme_nr + 
  theme(legend.position = "none",
        axis.line.y = element_blank(),
        panel.grid.major.y = element_line(colour = "grey", size=1),
        panel.background = element_rect(),
        axis.text = element_text(size=12))

if(interactive()) print(Figure3)

## ----fig4--------------------------------------------

parsets = as.integer(c(4, 5, 6))
vars = c("rel_kld") #"cmp_kld2", "nb_kld2", 

data = plyr::ldply(ode_processed[parsets], function(z) z$ode_stats_tidy) %>% 
  mutate(r = plyr::mapvalues(as.integer(.id), parsets, parms_grid$r[parsets])) %>% 
  filter(variable %in% vars) %>% 
  mutate(variable = plyr::revalue(variable, c("var_mean_ratio" = "Var/Mean Ratio"))) %>% 
  droplevels

Figure4 = ggplot(data, aes(x = time, y = value, col = variable)) +
  scale_color_manual(values = c("blue", "red", "darkgreen")) +
  geom_line(lwd = 1) +
  facet_grid(. ~ r, scales = "free") +
  xlab("Time") + ylab("Relative KLD (Neg. Binomial - CMP)") +
  theme_nr + 
  theme(legend.position = "none",
        axis.line.y = element_blank(),
        panel.grid.major.y = element_line(colour = "grey", size=1),
        panel.background = element_rect(),
        axis.text = element_text(size=12))

if(interactive()) print(Figure4)