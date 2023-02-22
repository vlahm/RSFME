library(ggplot2)
library(dplyr)
library(tidyr)
library(feather)
library(macrosheds)
library(RColorBrewer)
library(lfstat)
library(here)
library(ggthemes)
library(patchwork)

source(here('source/flux_methods.R'))
source(here('paper','ts_simulation','calculate_truth_ts.R'))

area <- 42.4
site_code = 'w3'

# HBEF Flux Method Comparison
w3_chem <- read_feather('data/ms/hbef/stream_chemistry/w3.feather') %>%
    mutate(var = ms_drop_var_prefix(var)) %>%
    distinct(datetime, site_code, var, .keep_all = TRUE) %>%
    pivot_wider(names_from = var, values_from = val, id_cols = c('datetime', 'site_code')) %>%
    filter(!is.na(Ca), !is.na(spCond))

w3_flux <- read_feather('data/ms/hbef/stream_flux/w3.feather') %>%
    mutate(var = ms_drop_var_prefix(var)) %>%
    select(-ms_recommended) %>%
    distinct(wy, site_code, method, var, .keep_all = TRUE) %>%
    pivot_wider(names_from = var, values_from = val, id_cols = c('wy', 'site_code', 'method')) %>%
    filter(!is.na(Ca),
           site_code == 'w3')

w3_recc <- read_feather('data/ms/hbef/stream_flux/w3.feather') %>%
    mutate(var = ms_drop_var_prefix(var)) %>%
    filter(ms_recommended == 1,
           var == 'Ca') %>%
    select(-ms_recommended) %>%
    distinct(wy, site_code, method, var, .keep_all = TRUE) %>%
    pivot_wider(names_from = var, values_from = val, id_cols = c('wy', 'site_code', 'method')) %>%
    filter(!is.na(Ca), site_code == 'w3') %>%
    mutate(method = 'recommended') %>%
    select(wy, site_code, method, Ca)

# create multiple linear model
lm_fit <- lm(Ca ~ spCond, data=w3_chem)
summary(lm_fit)

ggplotRegression <- function (fit) {
    ggplot(fit$model, aes_string(x = names(fit$model)[2], y = names(fit$model)[1])) +
        geom_point() +
        stat_smooth(method = "lm", col = "red") +
        labs(title = paste("Adj R2 = ",signif(summary(fit)$adj.r.squared, 5),
                           "Intercept =",signif(fit$coef[[1]],5 ),
                           " Slope =",signif(fit$coef[[2]], 5),
                           " P =",signif(summary(fit)$coef[2,4], 5))) +
        theme(text = element_text(size = 26))
}

ggplotRegression(lm_fit)

# flux methods sf
w3_flux_methods <- w3_flux %>%
    select(wy, method, Ca)

# bring in 'true' flux Ca from sensor
w3_flux_true <- read_feather("data/ms/hbef/true/w3_sensor_wdisch.feather") %>%
    mutate(wy = water_year(date, origin = 'usgs')) %>%
    group_by(wy) %>%
    summarise(Ca = sum(IS_spCond, na.rm = TRUE)*lm_fit$coef[[2]]) %>%
    mutate(site_code = 'w3',
           method = 'true') %>%
    select(-Ca, Ca)

d <- read_feather('C:/Users/gubbi/desktop/w3_sensor_wdisch.feather') %>%
    mutate(wy = water_year(datetime, origin = 'usgs'))

w3_true <- tibble(wy = as.integer(),
                  Ca = as.numeric())
for(i in unique(water_year(d$date, origin = 'usgs'))){
    target_wy = as.integer(i)

    dn <- d %>%
        filter(wy == target_wy)

    w3_chem <- dn$IS_spCond*lm_fit$coef[[2]]

    w3_q <- dn %>%
        group_by(date) %>%
        summarize(q_lps = mean(IS_discharge)) %>%
        select(date, q_lps)

    truth <- calculate_truth(w3_chem, w3_q, period = 'annual')$estimate[1]

    out <- tibble(wy = target_wy,
                  Ca = truth)

    w3_true <- bind_rows(w3_true, out)
}

w3_true <- w3_true %>%
    mutate(method = 'true',
           wy = as.character(wy)) %>%
    select(wy, method, Ca)

# bring in published flux Ca
w3_flux_pub <- read.csv('paper/hbef_comparison_fig/hbef_published_flux/ws3_stream_monthly_flux_gHa.csv') %>%
    mutate(date = paste0(Year, '-', Month, '-', '01'),
           wy = water_year(date, origin = 'usgs')) %>%
    select(wy, Ca_flux) %>%
    group_by(wy) %>%
    summarize(Ca = sum(Ca_flux)/1000) %>%
    filter(Ca > 0) %>%
    mutate(site_code = 'w3',
           method = 'published') %>%
    select(-Ca, Ca)

w3_all <- bind_rows(w3_flux_methods, w3_flux_pub, w3_true, w3_recc)

# look at flux time series
fluxpal <- brewer.pal(n=9, name='Set1')
breaks <- c('average', 'pw', 'beale', 'rating','composite', 'true', 'published', 'recommended')
labels <- c('Average', 'LI', 'Beale', 'Rating','Composite', 'True', 'Published', 'Recommended')

p_ts <- w3_all %>%
    select(-site_code) %>%
    filter(as.integer(wy) > 2012,
           as.integer(wy) < 2018) %>%
ggplot( aes(x = as.integer(wy), y= Ca)) +
    geom_point(aes(color = method), size = 5) +
    geom_line(aes(color = method))+
    theme_few() +
    theme(text = element_text(size = 50))+
    theme(panel.grid.major = element_blank(),
          ## legend.position="none",
          text = element_text(size = 24),
          plot.title = element_text(size = 24, face = "bold")) +
    scale_color_manual(breaks = breaks,
                       values = fluxpal,
                       labels = labels)+
    labs(x = '', color = 'Method',
         y = 'Annual Ca Load (kg/ha/yr)')

p_ts

ggsave(filename = here('paper','hbef_comparison_fig', 'method_ts.png'), width = 8, height = 6)

# make 1:1 line figure
p_comp <- w3_all %>%
    left_join(w3_true, by = 'wy') %>%
    filter(!is.na(Ca.y),
           method.x != 'true',
           method.x != 'wrtds') %>%
    ggplot( aes(x = Ca.y, y= Ca.x)) +
    geom_point(aes(color = method.x), size = 5) +
    geom_abline(slope = 1)+
    theme_few() +
    theme(text = element_text(size = 50))+
    theme(panel.grid.major = element_blank(),
          ## legend.position="none",
          text = element_text(size = 24),
          plot.title = element_text(size = 24, face = "bold")) +
    scale_color_manual(breaks = breaks,
                       values = fluxpal,
                       labels = labels
                       )+
    labs(y = 'Estimated Load', color = 'Method',
         x = 'Sensor Derived Load')

p_comp

ggsave(filename = here('paper','hbef_comparison_fig', 'method_comparison.png'), width = 8, height = 6)
