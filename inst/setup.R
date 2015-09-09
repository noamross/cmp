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
library(noamtools)
library(cowplot)
library(R.cache)
library(extrafont)
library(scales)

source('run_sodp.R')

options(dplyr.show_progress = FALSE)

## ----setparms------------------------------------------------------------
baseparms = list(
  max_i = 100,
  is = 0:100,
  lambda = 0.0005,
  lambda_ex = 0.000000,
  alpha = 0.05, #been using 0.05
  alpha_power = 1,  # ONE OR TWO
  mu =   0.1,
  r = 0.0001, # THIS SHOULD VARY FROM 0.2 to 0.  Bifurcation is near 0.01
  d = 0.00000001,
  K = 1000,
  time_max = 100,
  prevent_inf = 0,
  prevent_ex = 0,
  macro_timestep = 1,
  micro_timestep = 0.1,
  micro_relax_steps = 0,
  delta = 0,
  project = FALSE,
  n_comp_sims = 100,
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
)

r_vals = c(0.001, 0.05, 0.2)
alpha_power_vals = c(1, 2)
parms_grid = expand.grid(r = r_vals, alpha_power = alpha_power_vals)
