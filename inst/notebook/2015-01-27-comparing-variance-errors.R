library(compoisson)
library(dplyr)
library(tidyr)
library(ggplot2)


parmvals = expand.grid(lambda = c(seq(0.1, 3, by=0.1)), nu =exp(seq(log(0.005), log(3), length.out = 50)))
parmvals %<>%
  group_by(lambda, nu) %>%
  mutate(com.var=com.var(lambda, nu),
         com_var=com_var(lambda, nu),
         com_var_approx=com_var_approx(lambda, nu),
         est_z_iter = exp(log(lambda)/nu),
         log_z = com_compute_z(lambda, nu),
         z = exp(log_z)) %>%
  gather(var_measure, variance, com.var, com_var, com_var_approx)

ggplot(parmvals, aes(x=est_z_iter, y=variance, col=var_measure)) +
  geom_line() +
  scale_x_log10(limits=c(1e-3, 1e9), breaks=10^(-3:9)) +
  scale_y_log10(limits=c(1e-2, 1e15), breaks=10^(-2:15))

parmvals %>%
  filter(var_measure == "com.var") %>%
  select(lambda, nu, est_z_iter) %>%
  spread(nu, est_z_iter) %>%
  print.data.frame

parmvals %>%
  filter(var_measure == "com.var") %>%
  select(lambda, nu, log_z) %>%
  spread(nu, log_z) %>%
  print.data.frame

parmvals %>%
  filter(var_measure == "com.var") %>%
  mutate(variance = signif(variance, 2)) %>%
  select(lambda, nu, variance) %>%
  spread(nu, variance) %>%
  print.data.frame()

parmvals %>%
  filter(var_measure == "com_var") %>%
  select(lambda, nu, variance) %>%
  mutate(variance = signif(variance, 2)) %>%
  spread(nu, variance) %>%
  print.data.frame


parmvals %>%
  filter(var_measure == "com_var_approx") %>%
  select(lambda, nu, variance) %>%
  spread(nu, variance) %>%
  print.data.frame

parmvals %>%
  spread(var_measure, variance) %>%
  mutate(ratio1 = com_var_approx - com_var) %>%
  select(lambda, nu, ratio1) %>%
  spread(nu, ratio1) %>%
  print.data.frame
  
parmvals %>%
  filter(est_z_iter > 1e-3 & est_z_iter < 1e7 & log_z != Inf) %>%
  spread(var_measure, variance) %>%
  mutate(ratio1 = (com_var-com_var_approx)/com_var) %>%
  ggplot(aes(x = est_z_iter, y=ratio1)) +
    geom_point() +
    scale_x_log10() + scale_y_log10()

parmvals %>%
  filter(log_z != Inf & est_z_iter != Inf) %>%
ggplot(aes(x = est_z_iter, y = log_z)) + geom_line() +
  scale_x_log10(limits=c(1e-3, 1e4), breaks=10^(-3:4)) + 
  scale_y_log10(limits=c(1e-2, 1e30), breaks=10^(-2:30))