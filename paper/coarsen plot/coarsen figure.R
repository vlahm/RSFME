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

# set watershed attributes ####
area <- 42.4
site_code = 'w3'
#target_solute = 'IS_NO3'
target_solute = 'IS_spCond'

load(file = here('paper','coarsen plot', '100reps_annual_Ca.RData'))


plot_tbl <- out_tbl %>%
    unique() %>%
    mutate(error = ((estimate-truth$estimate[1])/truth$estimate[1])*100,
           error_abs = abs(error),
           method = factor(method, levels = c('pw', 'beale', 'rating', 'composite')),
           percent_coverage = (nrow(dn)/n)/nrow(dn),
           hours = n/4)

method_names <- c(
    `pw` = "Linear Interpolation",
    `beale` = "Beale",
    `rating` = "Rating",
    `composite` = "Composite"
)

breaks <- c(1,24,96,192,384,768)

x_labels <- c('Hourly', 'Daily', 'Weekly', 'Biweekly', 'Monthly', 'Bimonthly')

if(target_solute == 'IS_NO3'){
    y_min = -100
    y_max = 150
}
if(target_solute == ('IS_spCond')){
    y_min = -20
    y_max = 20
}

plot_tbl %>%
    group_by(method, hours) %>%
    mutate(min = min(error), max = max(error), median = median(error)) %>%
    filter(hours <= 899) %>%
    ggplot(., aes(x = hours, y = median))+
    annotate('rect', xmin = -Inf, xmax = Inf,
             ymin = -5, ymax = 5, fill = 'green', alpha = .1)+
    annotate('rect', xmin = -Inf, xmax = Inf,
             ymin = -20, ymax = -5, fill = 'yellow', alpha = .1)+
    annotate('rect', xmin = -Inf, xmax = Inf,
             ymin = 5, ymax = 20, fill = 'yellow', alpha = .1)+
    annotate('rect', xmin = -Inf, xmax = Inf,
             ymin = 20, ymax = Inf, fill = 'red', alpha = .1)+
    annotate('rect', xmin = -Inf, xmax = Inf,
             ymin = -Inf, ymax = -20, fill = 'red', alpha = .1)+
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
    labs(title = 'Nitrate Load Accuracy')
#ggsave(filename = here('paper','coarsen plot', 'nitrate_annual.png'), width = 14, height = 6)
#ggsave(filename = here('paper','coarsen plot', 'ca_annual.png'), width = 14, height = 6)

plot_tbl %>%
    filter(method == 'pw') %>%
    group_by(hours) %>%
    summarize(error_mean = mean(error),
              sd = sd(error),
              hi = error_mean+(sd*1.96),
              lo = error_mean -(sd*1.96),
              max = max(error),
              min = min(error)) %>%
    ggplot(aes(x = hours, y = error_mean))+
    geom_point()


