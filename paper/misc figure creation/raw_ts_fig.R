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

area <- 42.4
site_code = 'w3'
#target_solute = 'IS_NO3'


# read in data ####
d <- read_feather('C:/Users/gubbi/desktop/w3_sensor_wdisch.feather') %>%
    mutate(wy = water_year(datetime, origin = 'usgs'))
#slice(1:ts_len)

# subset to 2016 wy
target_wy <- 2016
target_solute = 'IS_spCond'

dn <- d %>%
    filter(wy == target_wy) %>%
    mutate(IS_discharge = na.approx(IS_discharge),
           IS_NO3 = na.approx(IS_NO3),
           IS_spCond = na.approx(IS_spCond)*0.06284158) %>%
    select(datetime, IS_spCond, IS_NO3, IS_discharge)


# supplement figures #####

q_breaks = c(1e-2, 1, 1e2)
q_labels = c('.01', '1', '100')

dn %>%
    ggplot(aes(x = IS_discharge, y = IS_spCond))+
    geom_point()+
    scale_y_log10()+
    scale_x_log10(breaks = q_breaks,
                  labels = q_labels) +
    theme_classic() +
    geom_smooth(method = 'lm')+
    labs(x = 'Q (lps)',
         y = 'C (mg/L)',
         title = 'Calcium C:Q Relationship')+
    theme(text = element_text(size = 20))
#ggsave(filename = here('paper','coarsen plot', 'ca_cq.png'), width = 6, height = 6)

dn %>%
    rename(Nitrate = IS_NO3,
           Calcium = IS_spCond) %>%
    pivot_longer(cols = c(Nitrate, Calcium),
                 names_to = 'var',
                 values_to = 'val') %>%
    ggplot(aes(x = datetime, y = val)) +
    geom_line()+
    facet_wrap(~var, ncol = 1, scales = 'free')+
    theme_classic()+
    theme(text = element_text(size = 20))+
    labs(x = '',
         y = 'C (mg/L)',
         title = 'HBEF Watershed 3 - 2016 Water Year')

#ggsave(filename = here('paper','coarsen plot', 'rawchem.png'), width = 12, height = 6)

dn %>%
    ggplot(aes(x = datetime, y = IS_discharge)) +
    geom_line()+
    theme_classic()+
    theme(text = element_text(size = 20))+
    scale_y_log10(breaks = q_breaks,
                  labels = q_labels)+
    labs(x = '',
         y = 'Q (Lps)',
         title = 'Streamflow at HBEF Watershed 3 - 2016 Water Year')

#ggsave(filename = here('paper','coarsen plot', 'rawQ.png'), width = 12, height = 6)
