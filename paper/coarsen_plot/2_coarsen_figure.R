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

set.seed(53045)

source(here('source/flux_methods.R'))

# create calcium figure #####
## set watershed attributes ####
area <- 42.4
site_code = 'w3'
target_solute = 'IS_spCond'


## read in and prep data ####
d <- read_feather('C:/Users/gubbi/desktop/w3_sensor_wdisch.feather') %>%
    mutate(wy = water_year(datetime, origin = 'usgs'))

#### subset to 2016 wy ####
target_wy <- 2016
dn <- d %>%
    filter(wy == target_wy) %>%
    select(date, all_of(target_solute), IS_discharge)
colnames(dn)[2] <- 'con'
#### convert from specific conductivity to calcium ####
dn$con <- dn$con*0.06284158

load(file = here('paper','coarsen_plot', '100reps_annual_Ca.RData'))

## calculate 'truth' ####
chem_df <- dn %>%
    group_by(lubridate::yday(date)) %>%
    summarize(date = date(date),
              con = mean(con)) %>%
    ungroup() %>%
    unique() %>%
    select(date, con) %>%
    mutate(site_code = 'w3', wy = target_wy)

q_df <- dn %>%
    select(date, q_lps = IS_discharge)%>%
    group_by(lubridate::yday(date)) %>%
    summarize(date = date(date),
              q_lps = mean(q_lps)) %>%
    ungroup() %>%
    unique() %>%
    mutate(site_code = 'w3', wy = target_wy)

out_val <- generate_residual_corrected_con(chem_df = chem_df, q_df = q_df, sitecol = 'site_code') %>%
    rename(datetime = date) %>%
    calculate_composite_from_rating_filled_df() %>%
    pull(flux)
truth <- tibble(method = 'truth', estimate = out_val)

## calculate error from truth #####
plot_tbl <- out_tbl %>%
    unique() %>%
    mutate(error = ((estimate-truth$estimate[1])/truth$estimate[1])*100,
           error_abs = abs(error),
           method = factor(method, levels = c('pw', 'beale', 'rating', 'composite')),
           percent_coverage = (nrow(dn)/n)/nrow(dn),
           hours = n/4)

## set names and breaks #####
method_names <- c(
    `pw` = "Linear Interpolation",
    `beale` = "Beale",
    `rating` = "Rating",
    `composite` = "Composite"
)

breaks <- c(1,24,96,192,384,768)

x_labels <- c('Hourly', 'Daily', 'Weekly', 'Biweekly', 'Monthly', 'Bimonthly')

y_min = -20
y_max = 20

## generate plot ####
plot_tbl %>%
    group_by(method, hours) %>%
    mutate(min = min(error), max = max(error), median = median(error)) %>%
    filter(hours <= 899) %>%
    ggplot(., aes(x = hours, y = median))+
    annotate('rect', xmin = -Inf, xmax = Inf,
             ymin = -5, ymax = 5, fill = 'green4', alpha = .25)+
    annotate('rect', xmin = -Inf, xmax = Inf,
             ymin = -20, ymax = -5, fill = 'yellow', alpha = .25)+
    annotate('rect', xmin = -Inf, xmax = Inf,
             ymin = 5, ymax = 20, fill = 'yellow', alpha = .25)+
    # annotate('rect', xmin = -Inf, xmax = Inf,
    #          ymin = 20, ymax = Inf, fill = 'red', alpha = .1)+
    # annotate('rect', xmin = -Inf, xmax = Inf,
    #          ymin = -Inf, ymax = -20, fill = 'red', alpha = .1)+
    geom_hline(yintercept = 0)+
    geom_line(size = 1)+
    geom_line(aes(y = max), linetype = 'dashed', size = .75)+
    geom_line(aes(y = min), linetype = 'dashed', size = .75)+
    #geom_point()+
    #geom_ribbon(aes(ymin = min, ymax = max), alpha = .2 )+
    facet_wrap(vars(method), ncol = 2, labeller = as_labeller(method_names))+
    #scale_y_reverse(limits = c(100,0)) +
    labs(x = 'Frequency',
         y = 'Error (%)'
         # y = 'Estimate (kg/hr/yr)',
         # caption = '15 minute NO3 data from HBEF W3 2016 WY resampled by every nth measurement, compared to truth using every sample and the composite method.
         # \n Vertical bars indicate hourly, daily, weekly, biweekly, monthly, and bimonthly intervals, black line is the median prediction and grey area the range of possible predictions.'
    )+
    theme_classic()+
    scale_x_continuous(breaks = breaks, labels = x_labels,guide = guide_axis(check.overlap = TRUE)
    )+
    scale_y_continuous(limits = c(y_min, y_max))+
    theme(text = element_text(size = 20),
          axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 20),
          panel.spacing = unit(.25,'lines'))+
    geom_vline(xintercept = 1)+ #hourly
    geom_vline(xintercept = 24)+ #daily
    geom_vline(xintercept = 96)+ #weekly
    geom_vline(xintercept = 192)+ #biweekly
    geom_vline(xintercept = 384)+ #monthly
    geom_vline(xintercept = 768)+ #bimonthly
    labs(title = 'Calcium Load Accuracy')
ggsave(filename = here('paper','coarsen_plot', 'ca_annual.png'), width = 13, height = 6)

# create nitrate figure #####
## set watershed attributes ####
area <- 42.4
site_code = 'w3'
target_solute = 'IS_NO3'

## read in and prep data ####
load(file = here('paper','coarsen_plot', '100reps_annual_NO3.RData'))
#### subset to 2016 wy ####
target_wy <- 2016
dn <- d %>%
    filter(wy == target_wy) %>%
    select(date, all_of(target_solute), IS_discharge)
colnames(dn)[2] <- 'con'

## calculate 'truth' ####
chem_df <- dn %>%
    group_by(lubridate::yday(date)) %>%
    summarize(date = date(date),
              con = mean(con)) %>%
    ungroup() %>%
    unique() %>%
    select(date, con) %>%
    mutate(site_code = 'w3', wy = target_wy)

q_df <- dn %>%
    select(date, q_lps = IS_discharge)%>%
    group_by(lubridate::yday(date)) %>%
    summarize(date = date(date),
              q_lps = mean(q_lps)) %>%
    ungroup() %>%
    unique() %>%
    mutate(site_code = 'w3', wy = target_wy)

out_val <- generate_residual_corrected_con(chem_df = chem_df, q_df = q_df, sitecol = 'site_code') %>%
    rename(datetime = date) %>%
    calculate_composite_from_rating_filled_df() %>%
    pull(flux)
truth <- tibble(method = 'truth', estimate = out_val)

## calculate error from truth #####
plot_tbl <- out_tbl %>%
    unique() %>%
    mutate(error = ((estimate-truth$estimate[1])/truth$estimate[1])*100,
           error_abs = abs(error),
           method = factor(method, levels = c('pw', 'beale', 'rating', 'composite')),
           percent_coverage = (nrow(dn)/n)/nrow(dn),
           hours = n/4)

## set names and breaks #####
method_names <- c(
    `pw` = "Linear Interpolation",
    `beale` = "Beale",
    `rating` = "Rating",
    `composite` = "Composite"
)

breaks <- c(1,24,96,192,384,768)

x_labels <- c('Hourly', 'Daily', 'Weekly', 'Biweekly', 'Monthly', 'Bimonthly')

y_min = -100
y_max = 150

## generate plot ####
plot_tbl %>%
    group_by(method, hours) %>%
    mutate(min = min(error), max = max(error), median = median(error)) %>%
    filter(hours <= 899) %>%
    ggplot(., aes(x = hours, y = median))+
    annotate('rect', xmin = -Inf, xmax = Inf,
             ymin = -5, ymax = 5, fill = 'green4', alpha = .25)+
    annotate('rect', xmin = -Inf, xmax = Inf,
             ymin = -20, ymax = -5, fill = 'yellow', alpha = .25)+
    annotate('rect', xmin = -Inf, xmax = Inf,
             ymin = 5, ymax = 20, fill = 'yellow', alpha = .25)+
    # annotate('rect', xmin = -Inf, xmax = Inf,
    #          ymin = 20, ymax = Inf, fill = 'red', alpha = .1)+
    # annotate('rect', xmin = -Inf, xmax = Inf,
    #          ymin = -Inf, ymax = -20, fill = 'red', alpha = .1)+
    geom_hline(yintercept = 0)+
    geom_line(size = 1)+
    geom_line(aes(y = max), linetype = 'dashed', size = .75)+
    geom_line(aes(y = min), linetype = 'dashed', size = .75)+
    #geom_point()+
    #geom_ribbon(aes(ymin = min, ymax = max), alpha = .2 )+
    facet_wrap(vars(method), ncol = 2, labeller = as_labeller(method_names))+
    #scale_y_reverse(limits = c(100,0)) +
    labs(x = 'Frequency',
         y = 'Error (%)'
         # y = 'Estimate (kg/hr/yr)',
         # caption = '15 minute NO3 data from HBEF W3 2016 WY resampled by every nth measurement, compared to truth using every sample and the composite method.
         # \n Vertical bars indicate hourly, daily, weekly, biweekly, monthly, and bimonthly intervals, black line is the median prediction and grey area the range of possible predictions.'
    )+
    theme_classic()+
    scale_x_continuous(breaks = breaks, labels = x_labels,guide = guide_axis(check.overlap = TRUE)
    )+
    scale_y_continuous(limits = c(y_min, y_max))+
    theme(text = element_text(size = 20),
          axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 20),
          panel.spacing = unit(.25,'lines'))+
    geom_vline(xintercept = 1)+ #hourly
    geom_vline(xintercept = 24)+ #daily
    geom_vline(xintercept = 96)+ #weekly
    geom_vline(xintercept = 192)+ #biweekly
    geom_vline(xintercept = 384)+ #monthly
    geom_vline(xintercept = 768)+ #bimonthly
    labs(title = 'Nitrate Load Accuracy')+
    coord_cartesian(ylim = c(-50,50))

ggsave(filename = here('paper','coarsen_plot', 'nitrate_annual.png'), width = 13, height = 6)

