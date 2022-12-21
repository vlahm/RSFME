library(devtools)
devtools::install_github('https://github.com/cran/EcoHydRology')
library(EcoHydRology)



obj <- BaseflowSeparation(simulated_series[[1]], filter_parameter = 0.925, passes = 3)


sum(obj$qft)/(sum(obj$bt)+sum(obj$qft))

obj <- BaseflowSeparation(simulated_series[[2]], filter_parameter = 0.925, passes = 3)


sum(obj$qft)/(sum(obj$bt)+sum(obj$qft))

obj <- BaseflowSeparation(simulated_series[[3]], filter_parameter = 0.925, passes = 3)


sum(obj$qft)/(sum(obj$bt)+sum(obj$qft))
