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
library(cowplot)
library(zoo)

set.seed(53045)


source(here('source/flux_methods.R'))
# make random sampling function
rtnorm <- function(n, mean, sd, a = 0, b = 5){
    qnorm(runif(n, pnorm(a, mean, sd), pnorm(b, mean, sd)), mean, sd)
}

# Set watershed attributes
area <- 42.4
site_code = 'w3'
# Thinning Frequency Loop Start #####
thin_freqs <- c('weekly','biweekly', 'monthly')
for(n in 1:3){
thin_freq = thin_freqs[n]
#period <- 'month'
period <- 'annual'
reps = 100

# Set data coarsening function
if(thin_freq == 'weekly'){
    coarsen_data <- function(chem_df){
        out <- chem_df %>%
            filter(hour(datetime) %in% c(13:18)) %>%
            filter(lubridate::mday(datetime) %in% c(1, 8, 15, 23)) %>%
            mutate(date = lubridate::date(datetime)) %>%
            distinct(date, .keep_all = T)
        return(out)
    }
}

if(thin_freq == 'biweekly'){
    coarsen_data <- function(chem_df){
        out <- chem_df %>%
            filter(hour(datetime) %in% c(13:18)) %>%
            filter(lubridate::mday(datetime) %in% c(1, 15)) %>%
            mutate(date = lubridate::date(datetime)) %>%
            distinct(date, .keep_all = T)
        return(out)
    }
}

if(thin_freq == 'monthly'){
    coarsen_data <- function(chem_df){
        out <- chem_df %>%
            filter(hour(datetime) %in% c(13:18)) %>%
            filter(lubridate::mday(datetime) %in% c(1)) %>%
            mutate(date = lubridate::date(datetime)) %>%
            distinct(date, .keep_all = T)
        return(out)
    }
}

# Repitition Loop Start #####
# Initialize output by aggregation period
if(period == 'annual'){
loop_out <- tibble(method = as.character(), estimate = as.numeric(),
                  flow = as.character(), cq = as.character(), runid = as.integer())
}
if(period == 'month'){
    loop_out <- tibble(method = as.character(), date = as.character(), estimate = as.numeric(),
                       flow = as.character(), cq = as.character(), runid = as.integer())
}
for(i in 1:reps){
# Read in data
d <- read_feather('C:/Users/gubbi/desktop/w3_sensor_wdisch.feather') %>%
    mutate(wy = water_year(datetime, origin = 'usgs'))

## Subset to 2016 wy
target_wy <- 2016
dn <- d %>%
    filter(wy == target_wy)

## fit ARIMA model to series ####
fit <- auto.arima(xts(dn$IS_discharge, order.by = dn$datetime))

## resample residuals ####
resampled_residuals = sample(fit$residuals,
                             #size = ts_len,
                             replace = TRUE)

# Initialize output of simulated series
simulated_series = list()

# sum total q for hold factor ####
sum_q <- sum(dn$IS_discharge)

## Flow series creation #####
### unaltered #####
reg <- dn$IS_discharge
simulated_series[[1]] = reg + resampled_residuals + 5
simulated_series[[1]][which(simulated_series[[1]] <= 0)] = 0.1
hold_factor <- (sum_q/sum(simulated_series[[1]]))
simulated_series[[1]] <- simulated_series[[1]]*hold_factor

### stormflow dominated ####
storm <- dn$IS_discharge
simulated_series[[2]] = storm^1.5 + resampled_residuals + 5
simulated_series[[2]][which(simulated_series[[2]] <= 0)] = 0.1
hold_factor <- (sum_q/sum(simulated_series[[2]]))
simulated_series[[2]] <- simulated_series[[2]]*hold_factor

###baseflow dominated ####
base <- dn$IS_discharge
simulated_series[[3]] <- base^0.9 + resampled_residuals + 5
simulated_series[[3]] <- rollmean(simulated_series[[3]], k = 10, fill = T)
simulated_series[[3]][which(simulated_series[[3]] <= 0)] = 0.1
hold_factor <- (sum_q/sum(simulated_series[[3]]))
simulated_series[[3]] <- simulated_series[[3]]*hold_factor

## Chemistry series creation  ####
# randomly sample for chemostatic
simulated_series[[4]] <- rtnorm(n = nrow(dn), sd = 0.1, mean = 2)

### no pattern ####
# randomly sample for no pattern
simulated_series[[5]] <- rtnorm(n = nrow(dn), sd = 0.5, mean = 2)

### enriching ####
error_vec <- rnorm(length(simulated_series[[1]]), mean = 1, sd = 0.1)
simulated_series[[6]] <- (10^((log10(simulated_series[[1]])*1)-1))*error_vec*hold_factor
en_unalt <- simulated_series[[6]]

simulated_series[[8]] <- (10^((log10(simulated_series[[2]])*1)-1))*error_vec*hold_factor
en_storm <- simulated_series[[8]]

simulated_series[[9]] <- (10^((log10(simulated_series[[3]])*1)-1))*error_vec*hold_factor
en_base <- simulated_series[[9]]

### simple dilution ####
error_vec <- rnorm(length(simulated_series[[1]]), mean = 1, sd = 0.1)
simulated_series[[7]] <- (10^((log10(simulated_series[[1]])*-1)+1.25))*error_vec*hold_factor
di_unalt <- simulated_series[[7]]
plot(log10(simulated_series[[7]])~log10(simulated_series[[1]]))

simulated_series[[10]] <- (10^((log10(simulated_series[[2]])*-1)+1.25))*error_vec*hold_factor
di_storm <- simulated_series[[10]]

simulated_series[[11]] <- (10^((log10(simulated_series[[3]])*-1)+1.25))*error_vec*hold_factor
di_base <- simulated_series[[11]]

# Estimate flux #####
## daily q aggregation function #####
make_q_daily <- function(q_df){
out <- q_df %>%
    group_by(lubridate::yday(datetime)) %>%
    summarize(date = date(datetime),
              q_lps = mean(q_lps)) %>%
    ungroup() %>%
    unique() %>%
    select(date, q_lps)
}

# truth calculation function
source(here('paper','ts_simulation','calculate_truth_ts.R'))


# apply methods function
apply_methods <- function(chem_df, q_df, period = period, flow_regime = NULL, cq = NULL){
    if(period == 'annual'){
    out <- tibble(method = as.character(), estimate = as.numeric(),
                  flow = as.character(), cq = as.character())
    #pw
    out[1,2] <- calculate_pw(chem_df, q_df)
    #beale
    out[2,2] <- calculate_beale(chem_df, q_df)
    #rating
    out[3,2] <- calculate_rating(chem_df, q_df)
    #comp
    out[4,2] <- generate_residual_corrected_con(chem_df = chem_df, q_df = q_df, sitecol = 'site_code') %>%
        rename(datetime = date) %>%
        calculate_composite_from_rating_filled_df() %>%
        pull(flux)

    out$method <- c('pw', 'beale', 'rating', 'composite')
    out$flow <- flow_regime
    out$cq <- cq
    }
    if(period == 'month'){
    out <- tibble(method = as.character(), date = as_date(NA), estimate = as.numeric(),
                  flow = as.character(), cq = as.character())
    #pw
    out_pw <- calculate_pw(chem_df, q_df, datecol = 'date', period = 'month') %>%
        mutate(method = 'pw') %>%
        select(method, date, estimate = flux)
    #beale
    out_beale <- calculate_beale(chem_df, q_df, datecol = 'date', period = 'month') %>%
        mutate(method = 'beale')%>%
        select(method, date = date, estimate = flux)
    #rating
    out_rating <- calculate_rating(chem_df, q_df, datecol = 'date', period = 'month') %>%
        mutate(method = 'rating')%>%
        select(method, date, estimate = flux)
    #comp
    out_comp <- generate_residual_corrected_con(chem_df = chem_df, q_df = q_df, sitecol = 'site_code') %>%
        select(datetime = date, q_lps, con, con_com, wy) %>%
        calculate_composite_from_rating_filled_df(., period = 'month') %>%
        mutate(method = 'composite',
               date_fixed = paste0(
                   str_split_fixed(as.character(date), '-', n = 3)[,1],
                   '-',
                   str_split_fixed(as.character(date), '-', n = 3)[,2]
                                      )
               ) %>%
        ungroup()%>%
        select(method, date = date_fixed,
               estimate = flux)

    out <- rbind(out_pw, out_beale, out_rating, out_comp) %>%
        mutate(date)

    out$flow <- flow_regime
    out$cq <- cq
    }


    return(out)


}
# initialize output for load estimates
if(period == 'annual'){
run_out <- tibble(method = as.character(), estimate = as.numeric(),
                  flow = as.character(), cq = as.character())
}
if(period == 'month'){
    run_out <- tibble(method = as.character(), date = as.character(), estimate = as.numeric(),
                      flow = as.character(), cq = as.character())
}
# Each set of estimates is generated with the same method
# First the truth is calculated  from the full simulated series
# Then each estimate is calculated from the corasened series
# Each are appended to the running output
### chemostatic ####
q_df <- make_q_daily(tibble(datetime = dn$datetime, q_lps = simulated_series[[1]])) %>%
    mutate(site_code = 'w3', wy = target_wy)

chem_df <- coarsen_data(tibble(datetime = dn$datetime, con = simulated_series[[4]])) %>%
    #filter(max(q_df$date) >= date,
    #       min(q_df$date) <= date) %>%
    mutate(site_code = 'w3', wy = target_wy)

#### unaltered flow ####
run_out <- rbind(
    calculate_truth(raw_chem_list = simulated_series[[4]], q_df, period = period, flow_regime = 'unaltered', cq = 'chemostatic'),
    run_out)
run_out <- rbind(
    apply_methods(chem_df, q_df, period = period, flow_regime = 'unaltered', cq = 'chemostatic'),
    run_out)

##### under storm flow ####
q_df <- make_q_daily(tibble(datetime = dn$datetime, q_lps = simulated_series[[2]])) %>%
    mutate(site_code = 'w3', wy = target_wy)
#truth
run_out <- rbind(
    calculate_truth(raw_chem_list = simulated_series[[4]], q_df, period = period, flow_regime = 'storm', cq = 'chemostatic'),
    run_out)
# apply
run_out <- rbind(
    apply_methods(chem_df, q_df, period = period, flow_regime = 'storm', cq = 'chemostatic'),
    run_out)

##### under base flow ####
q_df <- make_q_daily(tibble(datetime = dn$datetime, q_lps = simulated_series[[3]])) %>%
    mutate(site_code = 'w3', wy = target_wy)
#truth
run_out <- rbind(
    calculate_truth(raw_chem_list = simulated_series[[4]], q_df, period = period, flow_regime = 'base', cq = 'chemostatic'),
    run_out)
# apply
run_out <- rbind(
    apply_methods(chem_df, q_df, period = period,flow_regime = 'base', cq = 'chemostatic'),
    run_out)

### no pattern ####
chem_df <- coarsen_data(tibble(datetime = dn$datetime, con = simulated_series[[5]])) %>%
    filter(max(q_df$date) >= date,
           min(q_df$date) <= date) %>%
    mutate(site_code = 'w3', wy = target_wy)

#### under unaltered flow #####
q_df <- make_q_daily(tibble(datetime = dn$datetime, q_lps = simulated_series[[1]])) %>%
    mutate(site_code = 'w3', wy = target_wy)
#truth
run_out <- rbind(
    calculate_truth(raw_chem_list = simulated_series[[5]], q_df, period = period, flow_regime = 'unaltered', cq = 'none'),
    run_out)
# apply
run_out <- rbind(
    apply_methods(chem_df, q_df, period = period,flow_regime = 'unaltered', cq = 'none'),
    run_out)

#### under storm flow ####
q_df <- make_q_daily(tibble(datetime = dn$datetime, q_lps = simulated_series[[2]])) %>%
    mutate(site_code = 'w3', wy = target_wy)
#truth
run_out <- rbind(
    calculate_truth(raw_chem_list = simulated_series[[5]], q_df, period = period, flow_regime = 'storm', cq = 'none'),
    run_out)
# apply
run_out <- rbind(
    apply_methods(chem_df, q_df, period = period, flow_regime = 'storm', cq = 'none'),
    run_out)

#### under base flow ####
q_df <- make_q_daily(tibble(datetime = dn$datetime, q_lps = simulated_series[[3]])) %>%
    mutate(site_code = 'w3', wy = target_wy)
#truth
run_out <- rbind(
    calculate_truth(raw_chem_list = simulated_series[[4]], q_df, period = period, flow_regime = 'base', cq = 'none'),
    run_out)
# apply
run_out <- rbind(
    apply_methods(chem_df, q_df, period = period, flow_regime = 'base', cq = 'none'),
    run_out)

### strong enrich ####

#### under unaltered flow #####
chem_df <- coarsen_data(tibble(datetime = dn$datetime, con = en_unalt)) %>%
    filter(max(q_df$date) >= date,
           min(q_df$date) <= date) %>%
    mutate(site_code = 'w3', wy = target_wy)

q_df <- make_q_daily(tibble(datetime = dn$datetime, q_lps = simulated_series[[1]])) %>%
    mutate(site_code = 'w3', wy = target_wy)
#truth
run_out <- rbind(
    calculate_truth(raw_chem_list = en_unalt, q_df, period = period, flow_regime = 'unaltered', cq = 'enrich'),
    run_out)
# apply
run_out <- rbind(
    apply_methods(chem_df, q_df, period = period, flow_regime = 'unaltered', cq = 'enrich'),
    run_out)

#### under storm flow ####
chem_df <- coarsen_data(tibble(datetime = dn$datetime, con = en_storm)) %>%
    filter(max(q_df$date) >= date,
           min(q_df$date) <= date) %>%
    mutate(site_code = 'w3', wy = target_wy)

q_df <- make_q_daily(tibble(datetime = dn$datetime, q_lps = simulated_series[[2]])) %>%
    mutate(site_code = 'w3', wy = target_wy)
#truth
run_out <- rbind(
    calculate_truth(raw_chem_list = en_storm, period = period, q_df, flow_regime = 'storm', cq = 'enrich'),
    run_out)
# apply
run_out <- rbind(
    apply_methods(chem_df, q_df, period = period, flow_regime = 'storm', cq = 'enrich'),
    run_out)

#### under base flow ####
chem_df <- coarsen_data(tibble(datetime = dn$datetime, con = en_base)) %>%
    filter(max(q_df$date) >= date,
           min(q_df$date) <= date) %>%
    mutate(site_code = 'w3', wy = target_wy)

q_df <- make_q_daily(tibble(datetime = dn$datetime, q_lps = simulated_series[[3]])) %>%
    mutate(site_code = 'w3', wy = target_wy)
#truth
run_out <- rbind(
    calculate_truth(raw_chem_list = en_base, q_df, period = period, flow_regime = 'base', cq = 'enrich'),
    run_out)
# apply
run_out <- rbind(
    apply_methods(chem_df, q_df, period = period, flow_regime = 'base', cq = 'enrich'),
    run_out)

### simple dilution ####
#simple dilution

##### under unaltered flow #####
chem_df <- coarsen_data(tibble(datetime = dn$datetime, con = di_unalt)) %>%
    filter(max(q_df$date) >= date,
           min(q_df$date) <= date) %>%
    mutate(site_code = 'w3', wy = target_wy)

q_df <- make_q_daily(tibble(datetime = dn$datetime, q_lps = simulated_series[[1]])) %>%
    mutate(site_code = 'w3', wy = target_wy)
#truth
run_out <- rbind(
    calculate_truth(raw_chem_list = di_unalt, period = period, q_df, flow_regime = 'unaltered', cq = 'dilution'),
    run_out)
# apply
run_out <- rbind(
    apply_methods(chem_df, q_df, period = period, flow_regime = 'unaltered', cq = 'dilution'),
    run_out)

#### under storm flow ####
chem_df <- coarsen_data(tibble(datetime = dn$datetime, con = di_storm)) %>%
    filter(max(q_df$date) >= date,
           min(q_df$date) <= date) %>%
    mutate(site_code = 'w3', wy = target_wy)

q_df <- make_q_daily(tibble(datetime = dn$datetime, q_lps = simulated_series[[2]])) %>%
    mutate(site_code = 'w3', wy = target_wy)
#truth
run_out <- rbind(
    calculate_truth(raw_chem_list = di_storm, period = period, q_df, flow_regime = 'storm', cq = 'dilution'),
    run_out)
# apply
run_out <- rbind(
    apply_methods(chem_df, q_df, period = period, flow_regime = 'storm', cq = 'dilution'),
    run_out)

##### under base flow ####
chem_df <- coarsen_data(tibble(datetime = dn$datetime, con = di_base)) %>%
    filter(max(q_df$date) >= date,
           min(q_df$date) <= date) %>%
    mutate(site_code = 'w3', wy = target_wy)

q_df <- make_q_daily(tibble(datetime = dn$datetime, q_lps = simulated_series[[3]])) %>%
    mutate(site_code = 'w3', wy = target_wy)
#truth
run_out <- rbind(
    calculate_truth(raw_chem_list = di_base, q_df, period = period, flow_regime = 'base', cq = 'dilution'),
    run_out)
# apply
run_out <- rbind(
    apply_methods(chem_df, q_df, period = period, flow_regime = 'base', cq = 'dilution'),
    run_out)

### save out ####
loop_out <- run_out %>%
        mutate(runid = i) %>%
        rbind(., loop_out)
}
write_csv(loop_out, file = here('paper','ts_simulation', paste0(thin_freq,'Freq_',reps,'Reps20221221.csv')))
print(paste(thin_freq, ' done'))
}
save(simulated_series, file =  here('paper','ts_simulation', 'simulated_series.Rdata'))
