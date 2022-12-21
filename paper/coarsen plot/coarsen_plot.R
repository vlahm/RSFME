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
target_solute = 'IS_spCond'

# read in data ####
d <- read_feather('C:/Users/gubbi/desktop/w3_sensor_wdisch.feather') %>%
    mutate(wy = water_year(datetime, origin = 'usgs'))
#slice(1:ts_len)

# subset to 2016 wy
target_wy <- 2016
dn <- d %>%
    filter(wy == target_wy) %>%
    mutate(IS_discharge = na.approx(IS_discharge),
           IS_NO3 = na.approx(IS_NO3),
           IS_FDOM = na.approx(IS_FDOM),
           IS_spCond = na.approx(IS_spCond)) %>%
    select(date, all_of(target_solute), IS_discharge)
colnames(dn)[2] <- 'con'

if(target_solute == 'IS_spCond'){
    dn$con <- dn$con*0.06284158
}

# calculate truth ####
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

# make gradually coarsened chem ####
nth_element <- function(vector, starting_position, n) {
    vector[seq(starting_position, length(vector), n)]
}

coarse_chem <- list()
loopid = 0
loop_vec <- c(seq(from = 1, to = 92, by = 4),
              seq(from = 96, to = 672, by = 96),
              3592,
              6512)
reps <- 100
for(i in loop_vec){
n = i
    for(j in 1:reps){
    loopid <- loopid+1
    start_pos <- sample(1:n, size = 1)
    coarse_chem[[loopid]] <- tibble(date =  nth_element(dn$date, 1, n = start_pos),
                                    con = nth_element(dn$con, 1, n = start_pos))
    names(coarse_chem)[loopid] <- paste0('sample_',n)
    }
}




# apply flux methods to each #####
apply_methods_coarse <- function(chem_df, q_df){
    n = nrow(chem_df)

    out <- tibble(method = as.character(), estimate = as.numeric())#, se)
    #pw
    out[1,2] <- calculate_pw(chem_df, q_df)
    # theta <- function(x, chem_df, q_df){calculate_pw(chem_df[x,], q_df, datecol = 'datetime') }
    # out_jack[1,3] <- jackknife(1:n,theta, chem_df, q_df)$jack.se
    #beale
    out[2,2] <- calculate_beale(chem_df, q_df)
    # theta <- function(x, chem_df, q_df){calculate_beale(chem_df[x,], q_df, datecol = 'datetime') }
    # out_jack[2,3] <- jackknife(1:n,theta, chem_df, q_df)$jack.se
    #rating
    out[3,2] <- calculate_rating(chem_df, q_df)
    # theta <- function(x, chem_df, q_df){calculate_rating(chem_df[x,], q_df, datecol = 'datetime') }
    # out_jack[3,3] <- jackknife(1:n,theta, chem_df, q_df)$jack.se
    #comp
    out[4,2] <- generate_residual_corrected_con(chem_df = chem_df, q_df = q_df, sitecol = 'site_code') %>%
        rename(datetime = date) %>%
        calculate_composite_from_rating_filled_df() %>%
        pull(flux)

    # theta <- function(x, chem_df, q_df){generate_residual_corrected_con(chem_df = chem_df[x,], q_df = q_df, sitecol = 'site_code') %>%
    #         rename(datetime = date) %>%
    #         calculate_composite_from_rating_filled_df() %>%
    #         pull(flux) }
    # out_jack[4,3] <- jackknife(1:n,theta, chem_df, q_df)$jack.se

    out$method <- c('pw', 'beale', 'rating', 'composite')
    return(out)
}


out_tbl <- tibble(method = as.character(), estimate = as.numeric(), n = as.integer())
for(i in 2:length(coarse_chem)){

    n <- as.numeric(str_split_fixed(names(coarse_chem[i]), pattern = 'sample_', n = 2)[2])

    chem_df <- coarse_chem[[i]] %>%
        group_by(lubridate::yday(date)) %>%
        summarize(date = date(date),
                  con = mean(con)) %>%
        ungroup() %>%
        unique() %>%
        select(date, con) %>%
        mutate(site_code = 'w3', wy = target_wy)

    out_tbl <- apply_methods_coarse(chem_df, q_df) %>%
        mutate(n = n) %>%
        rbind(., out_tbl)
}
# save/load data from previous runs #####

#load('paper/coarsen plot/100repswithbootstrap_v2.RData')
#save(out_tbl, file = here('paper','coarsen plot', '100reps_annual_Ca.RData'))
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



