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
library(RiverLoad)

set.seed(53045)

source(here('source/flux_methods.R'))
source(here('paper/coarsen_plot/coarsen_helpers.R'))

# set watershed attributes #####
area <- 42.4
site_code = 'w3'

# begin solute loop ####
for(i in c('IS_NO3', 'IS_spCond')){
    ## set solute #####
    target_solute = i

    ## read in data ####
    d <- read_feather('C:/Users/gubbi/desktop/w3_sensor_wdisch.feather') %>%
        mutate(wy = water_year(datetime, origin = 'usgs'))

    ## subset to 2016 wy ####
    target_wy <- 2016
    dn <- d %>%
        filter(wy == target_wy) %>%
        select(date, con = all_of(target_solute), IS_discharge)

    # convert specific conductivity to calcium concentration
    # see the 'hbef ca correlation test' folder for more information
    if(target_solute == 'IS_spCond'){
        dn$con <- dn$con*0.06284158
    }

    ## calculate truth ####
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

    ## make gradually coarsened chem ###
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

    ## Start coarsening loop ####
    reps <- 100
    for(i in loop_vec){
        n = i
        print(paste0('i=', n))

        for(j in 1:reps){
            loopid <- loopid+1
            start_pos <- sample(1:n, size = 1) # take a random starting position from inside the interval
            coarse_chem[[loopid]] <- tibble(date =  nth_element(dn$date, 1, n = start_pos),
                                            con = nth_element(dn$con, 1, n = start_pos))
            names(coarse_chem)[loopid] <- paste0('sample_',n)
        }
        }

    ## Start method application loop ####
    out_tbl <- tibble(method = as.character(), estimate = as.numeric(), n = as.integer())
    for(k in 2:length(coarse_chem)){

        n <- as.numeric(str_split_fixed(names(coarse_chem[k]), pattern = 'sample_', n = 2)[2])

        chem_df <- coarse_chem[[k]] %>%
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

    ## save/load data from previous runs #####
    if(target_solute == 'IS_spCond'){save(out_tbl, file = here('paper','coarsen_plot', '100reps_annual_Ca.RData'))}
    if(target_solute == 'IS_NO3'){save(out_tbl, file = here('paper','coarsen_plot', '100reps_annual_NO3.RData'))}
}
