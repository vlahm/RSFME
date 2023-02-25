library(tidyverse)
library(forecast)
library(feather)
library(xts)
library(imputeTS)
library(here)
library(lfstat)
library(lubridate)
library(ggpubr)
library(patchwork)
library(ggthemes)

set.seed(53045)
source(here('paper/coarsen_plot/coarsen_helpers.R'))

# set watershed attributes #####
area <- 42.4
site_code = 'w3'

##set solute #####
target_solute = 'IS_NO3'

# read in data ####
d <- read_feather('C:/Users/gubbi/desktop/w3_sensor_wdisch.feather') %>%
    mutate(wy = water_year(datetime, origin = 'usgs'))

# subset to 2016 wy ####
target_wy <- 2016
dn <- d %>%
    filter(wy == target_wy) %>%
    select(datetime, con = all_of(target_solute))

# make gradually coarsened chem ###
## iniitalize output and loop
coarse_chem <- list()
loopid = 0

# create vector of nth elements
# go from full 15 minute data to daily by hour
# go from daily to biweekly by day
# set monthly and bimonthly discretely
loop_vec <- c(seq(from = 1, to = 92, by = 4),
              seq(from = 96, to = 672, by = 96),
              3592,
              6512)

# Start coarsening loop ####
reps <- 1
for(i in loop_vec){
    n = i
    print(paste0('i=', n))

    for(j in 1:reps){
        loopid <- loopid+1
        start_pos <- sample(1:n, size = 1) # take a random starting position from inside the interval
        coarse_chem[[loopid]] <- tibble(date =  nth_element(dn$datetime, 1, n = start_pos),
                                        con = nth_element(dn$con, 1, n = start_pos))
        names(coarse_chem)[loopid] <- paste0('sample_',n)
    }
}

# plot coarsened chem ####
full <- tibble(coarse_chem[1]$sample_1) %>%
    rename(`1` = con)
for(i in 2:length(coarse_chem)){
    n <- loop_vec[i]
    loop_dat <- tibble(coarse_chem[i][[1]])
    colnames(loop_dat)[2] <- n
    full <- full_join(full, loop_dat, by = 'date')

}

label_df <- tibble(label = c('15 Minute', 'Hourly', 'Daily','Weekly', 'Bimonthly'),
                   n =  as.character(c(1,5,96, 672,6512)))

plot_dat <- full %>%
    pivot_longer(names_to = 'n',
                 values_to = 'con',
                 cols = -date) %>%
    filter(n %in% c(1, 5, 96, 672, 6512)) %>%
    left_join(., label_df, by = 'n') %>%
    mutate(n = as.integer(n))
plot_dat$label <- factor(plot_dat$label, levels = c('15 Minute', 'Hourly', 'Daily','Weekly', 'Bimonthly'))

ggplot(plot_dat, aes(x = date, y = con))+
    geom_point()+
    geom_line()+
    facet_wrap(~label, ncol = 1) +
    theme_few()+
    theme(legend.position = 'none',
          text=element_text(size=20))+
    labs(y = 'Nitrate (mg/L)', x = '')

ggsave(filename = here('paper','coarsen_plot', 'nitrate_thinning.png'), width = 13, height = 6)
