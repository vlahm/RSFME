library(devtools)
devtools::install_github('https://github.com/cran/EcoHydRology')
library(EcoHydRology)

#load simulated_series from 1_ts_simulation_analysis.R ####
load(here('paper','ts_simulation', 'simulated_series.Rdata'))

# test unaltered synthetic ts ####
obj <- BaseflowSeparation(simulated_series[[1]], filter_parameter = 0.925, passes = 3)
sum(obj$qft)/(sum(obj$bt)+sum(obj$qft))

# test stormflow synthetic ts ####
obj <- BaseflowSeparation(simulated_series[[2]], filter_parameter = 0.925, passes = 3)
sum(obj$qft)/(sum(obj$bt)+sum(obj$qft))

#test unaltered synthetic ts ####
obj <- BaseflowSeparation(simulated_series[[3]], filter_parameter = 0.925, passes = 3)
sum(obj$qft)/(sum(obj$bt)+sum(obj$qft))
