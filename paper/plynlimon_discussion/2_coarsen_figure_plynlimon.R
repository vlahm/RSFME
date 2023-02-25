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

# read in shared data ####
area <- 122
site_code <- 'UHF'
target_wy <- 2008

d <- read_csv(here('paper/plynlimon_discussion/PlynlimonHighFrequencyHydrochemistry.csv')) %>%
    filter(Site == site_code) %>%
    select(date_time, `NO3-N mg/l`, `Ca mg/l`, `water flux mm/hr`) %>%
    mutate(wy = water_year(date_time, origin = 'usgs'),
           q_lps = `water flux mm/hr`*area*(1000/1)*(1/10000)*(1/3600)*(1000/1))



# create calcium figure #####
# set watershed attributes #####
target_solute = 'Ca mg/l'

dn <- d %>%
    filter(wy == target_wy) %>%
    select(date = date_time, con = all_of(target_solute), q_lps)

out_tbl <- read_csv(file = here('paper','plynlimon_discussion', '100reps_annual_Ca.csv'))

## calculate 'truth' ####
chem_df <- dn %>%
    group_by(lubridate::yday(date)) %>%
    summarize(date = date(date),
              con = mean(con)) %>%
    ungroup() %>%
    unique() %>%
    select(date, con) %>%
    mutate(site_code = site_code, wy = target_wy)

q_df <- dn %>%
    select(date, q_lps)%>%
    group_by(lubridate::yday(date)) %>%
    summarize(date = date(date),
              q_lps = mean(q_lps)) %>%
    ungroup() %>%
    unique() %>%
    mutate(site_code = site_code, wy = target_wy)

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
           hours = n*7)

## set names and breaks #####
method_names <- c(
    `pw` = "Linear Interpolation",
    `beale` = "Beale",
    `rating` = "Rating",
    `composite` = "Composite"
)

breaks <- c(24,96,192,384,768)

x_labels <- c('Daily', 'Weekly', 'Biweekly', 'Monthly', 'Bimonthly')

y_min = -30
y_max = 30

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
    #geom_vline(xintercept = 1)+ #hourly
    geom_vline(xintercept = 24)+ #daily
    geom_vline(xintercept = 96)+ #weekly
    geom_vline(xintercept = 192)+ #biweekly
    geom_vline(xintercept = 384)+ #monthly
    geom_vline(xintercept = 768)+ #bimonthly
    labs(title = 'Calcium Load Accuracy')+
    coord_cartesian(ylim = c(-25,25))
ggsave(filename = here('paper','plynlimon_discussion', 'ca_annual.png'), width = 13, height = 6)

# create nitrate figure #####
target_solute = 'NO3-N mg/l'

## subset data #####
dn <- d %>%
    filter(wy == target_wy) %>%
    select(date = date_time, con = all_of(target_solute), q_lps)

out_tbl <- read_csv(file = here('paper','plynlimon_discussion', '100reps_annual_NO3.csv'))

## calculate 'truth' ####
chem_df <- dn %>%
    group_by(lubridate::yday(date)) %>%
    summarize(date = date(date),
              con = mean(con)) %>%
    ungroup() %>%
    unique() %>%
    select(date, con) %>%
    mutate(site_code = site_code, wy = target_wy)

q_df <- dn %>%
    select(date, q_lps)%>%
    group_by(lubridate::yday(date)) %>%
    summarize(date = date(date),
              q_lps = mean(q_lps)) %>%
    ungroup() %>%
    unique() %>%
    mutate(site_code = site_code, wy = target_wy)

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
           hours = n*7)

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
    #geom_vline(xintercept = 1)+ #hourly
    geom_vline(xintercept = 24)+ #daily
    geom_vline(xintercept = 96)+ #weekly
    geom_vline(xintercept = 192)+ #biweekly
    geom_vline(xintercept = 384)+ #monthly
    geom_vline(xintercept = 768)+ #bimonthly
    labs(title = 'Nitrate Load Accuracy')+
    coord_cartesian(ylim = c(-25,25))
ggsave(filename = here('paper','plynlimon_discussion', 'nitrate_annual.png'), width = 13, height = 6)