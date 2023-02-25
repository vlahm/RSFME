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
area <- 122
site_code <- 'UHF'

# begin solute loop ####
for(s in c('NO3-N mg/l', 'Ca mg/l')){
    ## set solute #####
    target_solute = s

    ## read in data ####
    d <- read_csv(here('paper/plynlimon_discussion/PlynlimonHighFrequencyHydrochemistry.csv')) %>%
        filter(Site == site_code) %>%
        select(date_time, `NO3-N mg/l`, `Ca mg/l`, `water flux mm/hr`) %>%
        mutate(wy = water_year(date_time, origin = 'usgs'),
               q_lps = `water flux mm/hr`*area*(1000/1)*(1/10000)*(1/3600)*(1000/1)) #convert from mm to mm, ha to m2, hr to sec, and m3 to L

    ## subset to 2008 wy ####
    target_wy <- 2008
    dn <- d %>%
        filter(wy == target_wy) %>%
        select(date = date_time, con = all_of(target_solute), q_lps)

    ## calculate truth ####
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

    ## make gradually coarsened chem ###
    ## iniitalize output and loop
    coarse_chem <- list()
    loopid = 0

    # create vector of nth elements
    # go from 7 hour data to daily by 7 hours
    # go from daily to biweekly by day
    # set monthly and bimonthly discretely
    loop_vec <- c(seq(from = 1, to = 3, by = 1),
                  seq(from = 4, to = 48, by = 3),
                  96,
              192)

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
    for(k in 1:length(coarse_chem)){
        n <- as.numeric(str_split_fixed(names(coarse_chem[k]), pattern = 'sample_', n = 2)[2])
        chem_df <- coarse_chem[[k]] %>%
            group_by(lubridate::yday(date)) %>%
            summarize(date = date(date),
                      con = mean(con)) %>%
            ungroup() %>%
            unique() %>%
            select(date, con) %>%
            mutate(site_code = site_code, wy = target_wy)
        out_tbl <- apply_methods_coarse(chem_df, q_df) %>%
            mutate(n = n) %>%
            rbind(., out_tbl)
        print(paste0(k, ' DONE'))
    }
    ## save/load data from previous runs #####
    if(target_solute == 'Ca mg/l'){write_csv(out_tbl, file = here('paper','plynlimon_discussion', '100reps_annual_Ca.csv'))}
    if(target_solute == 'NO3-N mg/l'){write_csv(out_tbl, file = here('paper','plynlimon_discussion', '100reps_annual_NO3.csv'))}


}
