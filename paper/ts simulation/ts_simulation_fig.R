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

set.seed(53045)


source(here('source/flux_methods.R'))

thin_freq <- 'biweekly'
#thin_freq <- 'weekly'
#thin_freq <- 'monthly'
#period <- 'month'
period <- 'annual'
reps = 10


area <- 42.4
site_code = 'w3'
# loop start #####
if(period == 'annual'){
loop_out <- tibble(method = as.character(), estimate = as.numeric(),
                  flow = as.character(), cq = as.character(), runid = as.integer())
}
if(period == 'month'){
    loop_out <- tibble(method = as.character(), date = as.character(), estimate = as.numeric(),
                       flow = as.character(), cq = as.character(), runid = as.integer())
}

for(i in 1:reps){
#ts_len = 1000
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
           IS_spCond = na.approx(IS_spCond))

mean_q <- mean(dn$IS_discharge)
sum_q <- sum(dn$IS_discharge)

# max(dn$datetime)
# min(dn$datetime)
#
# max(dn$IS_NO3, na.rm = T)
# min(dn$IS_NO3, na.rm = T)
# sd(dn$IS_NO3, na.rm = T)
# mean(dn$IS_NO3, na.rm = T)
# hist(dn$IS_NO3)
#
# summary(lm(log10(dn$IS_NO3)~log10(dn$IS_discharge), data = dn))
#
# summary(lm(log10(dn$IS_spCond)~log10(dn$IS_discharge), data = dn))
#
# summary(lm(log10(dn$IS_FDOM)~log10(dn$IS_discharge), data = dn))


dn %>%
    select(datetime, IS_discharge, IS_NO3, IS_FDOM) %>%
    pivot_longer(cols = -datetime, values_to = 'val', names_to = 'var') %>%
    ggplot(., aes(x = datetime, y = val)) +
    #geom_point()+
    geom_line(lwd = 2)+
    facet_wrap(vars(var), ncol = 1, scales = 'free')+
    theme_classic()+
    theme(text = element_text(size = 20),
          axis.title = element_blank())

# plot(dn$datetime, dn$IS_discharge, type = 'l', lwd = 2)
# plot(dn$datetime, dn$IS_NO3, type = 'l', lwd = 2)
# plot(log10(dn$IS_discharge), log10(dn$IS_NO3), type = 'p', lwd = 2)


# DISCHARGE TS #####
#plot(dn$datetime, dn$IS_discharge, type = 'l', lwd = 2)

## fit ARIMA model to series; resample the residuals ####

fit <- auto.arima(xts(dn$IS_discharge, order.by = dn$datetime))
# methods(class='Arima')

#lines(dn$datetime, fit$fitted, col = 'blue', lwd = 2)

#ts_len = 1000
resampled_residuals = sample(fit$residuals,
                             #size = ts_len,
                             replace = TRUE)

## create random models
simulated_series = list()
# for(i in seq_len(10)){
#     Sys.sleep(1)
#     resampled_residuals = sample(fit$residuals, size = ts_len, replace = TRUE)
#     simulated_series[[i]] = d$val + resampled_residuals
#     lines(d$datetime, simulated_series[[i]], col = i)
# }

## make TS#####
###unaltered #####
reg <- dn$IS_discharge
reg[dn$IS_discharge < 1] = 1
simulated_series[[1]] = reg + resampled_residuals + 5
simulated_series[[1]][which(simulated_series[[1]] <= 0)] = 1
hold_factor <- (sum_q/sum(simulated_series[[1]]))
simulated_series[[1]] <- simulated_series[[1]]*hold_factor
#lines(dn$datetime, simulated_series[[1]], col = 'blue', type = 'l')

###stormflow dominated ####
storm <- dn$IS_discharge
storm[dn$IS_discharge < 1] = 1
simulated_series[[2]] = storm^1.5 + resampled_residuals + 5
simulated_series[[2]][which(simulated_series[[2]] <= 0)] = 1
hold_factor <- (sum_q/sum(simulated_series[[2]]))
simulated_series[[2]] <- simulated_series[[2]]*hold_factor
#lines(dn$datetime, simulated_series[[2]], col = 'red', type = 'l')

###baseflow dominated ####
base <- dn$IS_discharge
base[dn$IS_discharge < 1] = 1
simulated_series[[3]] = storm^0.8 + resampled_residuals + 5
simulated_series[[3]][which(simulated_series[[3]] <= 0)] = 1
hold_factor <- (sum_q/sum(simulated_series[[3]]))
simulated_series[[3]] <- simulated_series[[3]]*hold_factor
#lines(dn$datetime, simulated_series[[3]], col = 'green', type = 'l')


# plot(dn$datetime, dn$IS_discharge, type = 'l', lwd = 2)
# lines(dn$datetime, simulated_series[[1]], col = 'blue', type = 'l')
# lines(dn$datetime, simulated_series[[2]], col = 'red', type = 'l')
# lines(dn$datetime, simulated_series[[3]], col = 'green', type = 'l')

# CON TS ####

## make TS#####
# plot(dn$datetime, dn$IS_NO3, type = 'l', lwd = 2)

### chemostatic #####
# make random sampling function
rtnorm <- function(n, mean, sd, a = 0, b = 5){
    qnorm(runif(n, pnorm(a, mean, sd), pnorm(b, mean, sd)), mean, sd)
}
#### apply to make chemo #####
#simulated_series[[4]] <- rtnorm(n = nrow(dn), sd = (sd(dn$IS_NO3)), mean = mean(dn$IS_NO3))
simulated_series[[4]] <- rtnorm(n = nrow(dn), sd = 0.1, mean = 2)
plot(log10(simulated_series[[4]])~log10(simulated_series[[1]]))
#lines(dn$datetime, simulated_series[[4]], col = 'blue', type = 'l')

### no pattern ####
# apply to make no pattern ts
simulated_series[[5]] <- rtnorm(n = nrow(dn), sd = 0.5, mean = 2)
plot(log10(simulated_series[[5]])~log10(simulated_series[[1]]))
#lines(dn$datetime, simulated_series[[4]], col = 'blue', type = 'l')

# plot(simulated_series[[5]]~dn$IS_discharge, data = dn)

#### enriching ####
##### fit lm to sp ts ####
## commenting out fdom based values in favor of fully synthetic
# fit_fdom <- lm(log10(dn$IS_FDOM)~log10(dn$IS_discharge), data = dn)
# inter_range <- runif(1000, min = confint(fit_fdom)[1,1], max = confint(fit_fdom)[1,2])
# coef_range <- runif(1000, min = confint(fit_fdom)[2,1], max = confint(fit_fdom)[2,2])
# error_range <- rnorm(1000, mean = 0, sd = sd(dn$IS_FDOM)/2)
# simulated_series[[6]] <- as.numeric()
# ##### for all #####
# for(j in 1:length(simulated_series[[1]])){
#     inter <- sample(inter_range, size = 1)
#     slope <- sample(coef_range, size = 1)
#     error <- sample(error_range, size = 1)
#     q <- simulated_series[[1]][j]
#
#     pre_error_val <- 10^((log10(q)*slope)+inter)
#
#     eps <- pre_error_val+(error*(pre_error_val)/mean(dn$IS_FDOM))
#
#     simulated_series[[6]][j] <- pre_error_val + eps
#
# }
## fully synthetic effort
error_vec <- rnorm(length(simulated_series[[1]]), mean = 1, sd = 0.05)
simulated_series[[6]] <- (10^((log10(simulated_series[[1]])*1)+1))*error_vec
plot(log10(simulated_series[[6]])~log10(simulated_series[[1]]))

summary(lm(log10(simulated_series[[6]])~log10(simulated_series[[1]])))
#### two part dilution ####
## fully synthetic effort
error_vec <- rnorm(length(simulated_series[[1]]), mean = 1, sd = 0.05)
mild <- (10^((log10(simulated_series[[1]])*-0.1))+3)*error_vec
mean(mild)
simple <- (10^((log10(simulated_series[[1]])*-1)+2))*error_vec
mean(simple)
min(simple)
tibble(q = log10(simulated_series[[1]]),
       mild = log10(mild),
       heavy = log10(simple)) %>%
ggplot(.,aes(x = q))+
    geom_point(aes(y = mild))+
    geom_point(aes(y = heavy), color = 'red')

simulated_series[[7]]<-simulated_series[[1]]
low_thresh <- 1.45
high_thresh <- 1.47
for(j in 1:length(simulated_series[[1]])){

q <- log10(simulated_series[[1]][j])
if(q < low_thresh){simulated_series[[7]][j] <- mild[j]}

if(q > high_thresh){simulated_series[[7]][j] <- simple[j]}

if(q >= low_thresh &
   q <= high_thresh){
    choice <- sample(0:1, 1)
    if(choice == 0){simulated_series[[7]][j] <- mild[j]}
    if(choice == 1){simulated_series[[7]][j] <- simple[j]}
}

}

plot(log10(simulated_series[[7]])~log10(simulated_series[[1]]))
summary(lm(log10(simulated_series[[6]])~log10(simulated_series[[1]])))


# shelving this in favor of fully synthetic effort
# ##### fit lm to sp ts ####
# fit_cond <- lm(log10(dn$IS_spCond*0.06)~log10(dn$IS_discharge), data = dn)
# inter_range <- runif(1000, min = confint(fit_cond)[1,1], max = confint(fit_cond)[1,2])
# coef_range <- runif(1000, min = confint(fit_cond)[2,1], max = confint(fit_cond)[2,2])
# error_range <- rnorm(1000, mean = 0, sd = sd(dn$IS_spCond)/2)
# simulated_series[[7]] <- as.numeric()
# ##### for all #####
# for(j in 1:length(simulated_series[[1]])){
#     inter <- sample(inter_range, size = 1)
#     slope <- sample(coef_range, size = 1)
#     error <- sample(error_range, size = 1)
#     q <- simulated_series[[1]][j]
#
#     if(q <= 5){
#     pre_error_val <- 10^((log10(q)*slope)+inter)
#
#     eps <- pre_error_val+(error*(pre_error_val)/mean(dn$IS_spCond))
#
#     simulated_series[[7]][j] <- pre_error_val + eps
#     }
#     if(q > 5){
#         pre_error_val <- 10^((log10(q)*2*slope)+inter)
#
#         eps <- pre_error_val+(error*(pre_error_val)/mean(dn$IS_spCond))
#
#         simulated_series[[7]][j] <- pre_error_val + eps
#     }
#
# }


#plot(log10(dn$IS_FDOM)~log10(dn$IS_discharge), data = dn)

## check c:q #####
#plot(log10(simulated_series[[6]])~log10(simulated_series[[1]]))
#summary(lm(log10(simulated_series[[6]])~log10(simulated_series[[1]])))
#
# plot(log10(simulated_series[[6]])~log10(simulated_series[[2]]))
# summary(lm(log10(simulated_series[[6]])~log10(simulated_series[[2]])))
#
# plot(log10(simulated_series[[6]])~log10(simulated_series[[3]]))
# summary(lm(log10(simulated_series[[6]])~log10(simulated_series[[3]])))

# ESTIMATE FLUX #####
# coarsen function
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



make_q_daily <- function(q_df){
out <- q_df %>%
    group_by(lubridate::yday(datetime)) %>%
    summarize(date = date(datetime),
              q_lps = mean(q_lps)) %>%
    ungroup() %>%
    unique() %>%
    select(date, q_lps)
}

calculate_truth <- function(raw_chem_list, q_df, period = period, flow_regime = NULL, cq = NULL){
    if(period == 'annual'){
    chem_df <- tibble(datetime = dn$datetime, con = raw_chem_list) %>%
        group_by(lubridate::yday(datetime)) %>%
        summarize(date = date(datetime),
                  con = mean(con)) %>%
        ungroup() %>%
        unique() %>%
        select(date, con) %>%
        mutate(site_code = 'w3', wy = target_wy)

    q_df_add <- q_df %>%
        mutate(site_code = 'w3', wy = target_wy)

    out_val <- generate_residual_corrected_con(chem_df = chem_df, q_df = q_df_add, sitecol = 'site_code') %>%
        rename(datetime = date) %>%
        calculate_composite_from_rating_filled_df() %>%
        pull(flux)
    out <- tibble(method = 'truth', estimate = out_val,
                  flow = flow_regime, cq = cq)
    }
    if(period == 'month'){
    chem_df <- tibble(datetime = dn$datetime, con = raw_chem_list) %>%
        group_by(lubridate::yday(datetime)) %>%
        summarize(date = date(datetime),
                  con = mean(con)) %>%
        ungroup() %>%
        unique() %>%
        select(date, con) %>%
        mutate(site_code = 'w3', wy = target_wy)

    q_df_add <- q_df %>%
        mutate(site_code = 'w3', wy = target_wy)

    out_val <- generate_residual_corrected_con(chem_df = chem_df, q_df = q_df, sitecol = 'site_code') %>%
        select(datetime = date, q_lps, con, con_com, wy) %>%
        calculate_composite_from_rating_filled_df(., period = 'month') %>%
        mutate(method = 'truth',
               date_fixed = paste0(
                   str_split_fixed(as.character(date), '-', n = 3)[,1],
                   '-',
                   str_split_fixed(as.character(date), '-', n = 3)[,2]
               )
        ) %>%
        ungroup()%>%
        select(method, date = date_fixed,
               estimate = flux) %>%
        mutate(flow = flow_regime, cq = cq)

    out <- out_val
    }

    return(out)
}

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
if(period == 'annual'){
run_out <- tibble(method = as.character(), estimate = as.numeric(),
                  flow = as.character(), cq = as.character())
}
if(period == 'month'){
    run_out <- tibble(method = as.character(), date = as.character(), estimate = as.numeric(),
                      flow = as.character(), cq = as.character())
}
### chemostatic ####
q_df <- make_q_daily(tibble(datetime = dn$datetime, q_lps = simulated_series[[1]])) %>%
    mutate(site_code = 'w3', wy = target_wy)

chem_df <- coarsen_data(tibble(datetime = dn$datetime, con = simulated_series[[4]])) %>%
    #filter(max(q_df$date) >= date,
    #       min(q_df$date) <= date) %>%
    mutate(site_code = 'w3', wy = target_wy)

#truth
run_out <- rbind(
    calculate_truth(raw_chem_list = simulated_series[[4]], q_df, period = period, flow_regime = 'unaltered', cq = 'chemostatic'),
    run_out)
# apply
run_out <- rbind(
    apply_methods(chem_df, q_df, period = period, flow_regime = 'unaltered', cq = 'chemostatic'),
    run_out)

##### under storm domination ####
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

##### under base domination ####
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

##### under unaltered flow #####
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

##### under storm domination ####
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

##### under base domination ####
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
chem_df <- coarsen_data(tibble(datetime = dn$datetime, con = simulated_series[[6]])) %>%
    filter(max(q_df$date) >= date,
           min(q_df$date) <= date) %>%
    mutate(site_code = 'w3', wy = target_wy)

##### under unaltered flow #####
q_df <- make_q_daily(tibble(datetime = dn$datetime, q_lps = simulated_series[[1]])) %>%
    mutate(site_code = 'w3', wy = target_wy)
#truth
run_out <- rbind(
    calculate_truth(raw_chem_list = simulated_series[[6]], q_df, period = period, flow_regime = 'unaltered', cq = 'enrich'),
    run_out)
# apply
run_out <- rbind(
    apply_methods(chem_df, q_df, period = period, flow_regime = 'unaltered', cq = 'enrich'),
    run_out)

##### under storm domination ####
q_df <- make_q_daily(tibble(datetime = dn$datetime, q_lps = simulated_series[[2]])) %>%
    mutate(site_code = 'w3', wy = target_wy)
#truth
run_out <- rbind(
    calculate_truth(raw_chem_list = simulated_series[[6]], period = period, q_df, flow_regime = 'storm', cq = 'enrich'),
    run_out)
# apply
run_out <- rbind(
    apply_methods(chem_df, q_df, period = period, flow_regime = 'storm', cq = 'enrich'),
    run_out)

##### under base domination ####
q_df <- make_q_daily(tibble(datetime = dn$datetime, q_lps = simulated_series[[3]])) %>%
    mutate(site_code = 'w3', wy = target_wy)
#truth
run_out <- rbind(
    calculate_truth(raw_chem_list = simulated_series[[6]], q_df, period = period, flow_regime = 'base', cq = 'enrich'),
    run_out)
# apply
run_out <- rbind(
    apply_methods(chem_df, q_df, period = period, flow_regime = 'base', cq = 'enrich'),
    run_out)

### 2 part dilution ####
chem_df <- coarsen_data(tibble(datetime = dn$datetime, con = simulated_series[[7]])) %>%
    filter(max(q_df$date) >= date,
           min(q_df$date) <= date) %>%
    mutate(site_code = 'w3', wy = target_wy)

##### under unaltered flow #####
q_df <- make_q_daily(tibble(datetime = dn$datetime, q_lps = simulated_series[[1]])) %>%
    mutate(site_code = 'w3', wy = target_wy)
#truth
run_out <- rbind(
    calculate_truth(raw_chem_list = simulated_series[[7]], period = period, q_df, flow_regime = 'unaltered', cq = 'broken_dilution'),
    run_out)
# apply
run_out <- rbind(
    apply_methods(chem_df, q_df, period = period, flow_regime = 'unaltered', cq = 'broken_dilution'),
    run_out)

##### under storm domination ####
q_df <- make_q_daily(tibble(datetime = dn$datetime, q_lps = simulated_series[[2]])) %>%
    mutate(site_code = 'w3', wy = target_wy)
#truth
run_out <- rbind(
    calculate_truth(raw_chem_list = simulated_series[[7]], period = period, q_df, flow_regime = 'storm', cq = 'broken_dilution'),
    run_out)
# apply
run_out <- rbind(
    apply_methods(chem_df, q_df, period = period, flow_regime = 'storm', cq = 'broken_dilution'),
    run_out)

##### under base domination ####
q_df <- make_q_daily(tibble(datetime = dn$datetime, q_lps = simulated_series[[3]])) %>%
    mutate(site_code = 'w3', wy = target_wy)
#truth
run_out <- rbind(
    calculate_truth(raw_chem_list = simulated_series[[7]], q_df, period = period, flow_regime = 'base', cq = 'broken_dilution'),
    run_out)
# apply
run_out <- rbind(
    apply_methods(chem_df, q_df, period = period, flow_regime = 'base', cq = 'broken_dilution'),
    run_out)

### save out ####
loop_out <- run_out %>%
        mutate(runid = i) %>%
        rbind(., loop_out)
}
#save(loop_out, file = here('paper','ts simulation', 'biweekly.RData'))
#save(loop_out, file = here('paper','ts simulation', 'monthlyflux_biweeklythinning.RData'))
#load(file = here('paper','ts simulation', 'biweekly.RData'))
#load(file = here('paper','ts simulation', 'monthlyflux_biweeklythinning.RData'))

# Figure creation #####
# this is curently using the save/load objects to run
# eventually you will need to make this actually run with both
if(period == 'month'){
    loop_out_month <- loop_out %>%
        mutate(runid = paste0(date, '_', runid))
}
# # join monthly and annual data together
# loop_out <- loop_out_month %>%
#     select(-date)%>%
#     mutate(period = 'Month') %>%
#     rbind(., loop_out %>%
#               mutate(period = 'Year'))

### make header plots #####
side_ymin <- 0.01
side_ymax <- 1000
# unaltered q plot
p1 <- ggplot(dn, aes(x = date))+
        geom_line(aes(y = simulated_series[[1]])) +
        theme_classic()+
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.title.y=element_blank(),
          text = element_text(size = 20))+
    labs(title = 'Unaltered Flow')+
    scale_y_log10(limits = c(side_ymin,side_ymax))

p1

# storm q plot
p2 <- ggplot(dn, aes(x = date))+
    geom_line(aes(y = simulated_series[[2]])) +
    theme_classic()+
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          text = element_text(size = 20))+
    labs(title = 'Stormflow Dominated',
         y = 'Q (lps)')+
    scale_y_log10(limits = c(side_ymin,side_ymax))
p2

# base q plot
p3 <- ggplot(dn, aes(x = date))+
    geom_line(aes(y = simulated_series[[3]])) +
    theme_classic()+
    theme(axis.title.y=element_blank(),
          axis.title.x = element_blank(),
          text = element_text(size = 20),
          axis.text.x = element_text(angle = 45, hjust = 1)
          )+
    labs(title = 'Baseflow Dominated')+
    scale_y_log10(limits = c(side_ymin,side_ymax))+
    scale_x_continuous(breaks = c(as_date('2015-10-01'), as_date('2016-04-01'), as_date('2016-10-01')),
                       labels = c('10/2015', '4/2016', '10/2016'))


p3

# set common limits to top row graphs
top_row_breaks <- c(1e-2, 1, 1e2, 1e4)
top_row_labels <- c('0.01', '1', '100', '10,000')
top_row_ymax <- 1e4
top_row_ymin <- 1e-2
# chemo cq
p4 <- tibble(q = simulated_series[[1]], con = simulated_series[[4]]) %>%
    ggplot(aes(x = q, y = con)) +
    geom_point() +
    theme_classic()+
    scale_x_log10() +
    scale_y_log10(limits = c(top_row_ymin, top_row_ymax),
                  breaks = top_row_breaks,
                  labels = top_row_labels) +
    labs(title = 'Chemostatic',
         y = 'C (mg/L)')+
    theme(axis.title.x=element_blank(),
          text = element_text(size = 20))

p4

# no pattern cq
p5 <- tibble(q = simulated_series[[1]], con = simulated_series[[5]]) %>%
    ggplot(aes(x = q, y = con)) +
    geom_point() +
    theme_classic()+
    scale_x_log10() +
    scale_y_log10(limits = c(top_row_ymin, top_row_ymax),
                  breaks = top_row_breaks)+
    labs(title = 'No Pattern',
         x = 'Q (lps)')+
    theme(axis.title.y=element_blank(),
          text = element_text(size = 20))
p5

# enrich cq
p6 <- tibble(q = simulated_series[[1]], con = simulated_series[[6]]) %>%
    ggplot(aes(x = q, y = con)) +
    geom_point() +
    theme_classic()+
    scale_x_log10() +
    scale_y_log10(limits = c(top_row_ymin, top_row_ymax),
                  breaks = top_row_breaks)+
    labs(title = 'Enriching',
         x = 'Q')+
    theme(axis.title.y=element_blank(),
          axis.title.x=element_blank(),
          text = element_text(size = 20))
p6

# broekn dilute cq
p16 <- tibble(q = simulated_series[[1]], con = simulated_series[[7]]) %>%
    ggplot(aes(x = q, y = con)) +
    geom_point() +
    theme_classic()+
    scale_x_log10() +
    scale_y_log10(limits = c(top_row_ymin, top_row_ymax),
                  breaks = top_row_breaks)+
    labs(title = 'Two-Part Dilution',
         x = 'Q')+
    theme(axis.title.y=element_blank(),
          axis.title.x=element_blank(),
          text = element_text(size = 20))
p16

# broken dilution c:q

### make row 1 plots ####
ymin = -100
ymax = 100

plot_guts <- function(p){
    ggplot(p, aes(x = method, y = error))+
    geom_hline(yintercept = 0)+
    geom_boxplot(
        #aes(fill = period)
                 )+
    theme_classic()+
    theme(
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.title.y=element_blank(),
        text = element_text(size = 20),
        legend.position="none"
    )+
    ylim(ymin, ymax)
}

transform_loop_out <- function(loop_out){
    pivot_wider(loop_out, names_from = method, values_from = estimate,
                #id_cols = c(runid, period),
                id_cols = c(runid, thin_freq),
                values_fn = mean) %>%
    mutate(pw = ((pw-truth)/truth)*100,
           beale = ((beale - truth)/truth)*100,
           rating = ((rating-truth)/truth)*100,
           composite = ((composite - truth)/truth)*100) %>%
    select(-truth, -runid) %>%
    pivot_longer(cols = -thin_freq,
                 #cols = everything() ,
                 names_to = 'method', values_to = 'error')
}
#  spacer plot that is only legend
p0_data <- loop_out %>%
    transform_loop_out()


p0_data$method <- factor(p0_data$method, levels = c("pw", "beale", "rating", 'composite'))

p0 <-plot_guts(p0_data)+
    theme(legend.position = 'left',
          legend.title = element_text(size = 25),
          legend.text = element_text(size = 20)
          )+
    labs(fill = 'Frequency')

# now extract the legend
legend <- ggdraw(get_legend(p0))

# unaltered flow w/ chemo data
p7_data <- loop_out %>%
    filter(flow == 'unaltered',
           cq == 'chemostatic') %>%
    transform_loop_out()


p7_data$method <- factor(p7_data$method, levels = c("pw", "beale", "rating", 'composite'))

p7 <-plot_guts(p7_data)+
    theme(axis.title.y=element_text(size = 20))+
    labs(y = " ")

p7

# unaltered flow w/ no pattern data
p8_data <- loop_out %>%
    filter(flow == 'unaltered',
           cq == 'none') %>%
    transform_loop_out()

p8_data$method <- factor(p8_data$method, levels = c("pw", "beale", "rating", 'composite'))

p8 <- plot_guts(p8_data)
p8

# unaltered flow w/ enrich data
p9_data <- loop_out %>%
    filter(flow == 'unaltered',
           cq == 'enrich') %>%
    transform_loop_out()

p9_data$method <- factor(p9_data$method, levels = c("pw", "beale", "rating", 'composite'))

p9 <- plot_guts(p9_data)
p9

# unaltered flow w/ broken dilution data
p17_data <- loop_out %>%
    filter(flow == 'unaltered',
           cq == 'broken_dilution') %>%
    transform_loop_out()

p17_data$method <- factor(p17_data$method, levels = c("pw", "beale", "rating", 'composite'))

p17 <- plot_guts(p17_data)

p17

### make row 2 plots ####
# storm flow w/ chemo data
p10_data <- loop_out %>%
    filter(flow == 'storm',
           cq == 'chemostatic') %>%
    transform_loop_out()

p10_data$method <- factor(p10_data$method, levels = c("pw", "beale", "rating", 'composite'))

p10 <- plot_guts(p10_data)+
    theme(axis.title.y=element_text(angle = 90))+
    labs(y = 'Error (%)')
p10

# storm flow w/ no pattern data
p11_data <- loop_out %>%
    filter(flow == 'storm',
           cq == 'none') %>%
    transform_loop_out()

p11_data$method <- factor(p11_data$method, levels = c("pw", "beale", "rating", 'composite'))

p11 <- plot_guts(p11_data)

p11

# stormflow w/ enrich data
p12_data <- loop_out %>%
    filter(flow == 'storm',
           cq == 'enrich') %>%
    transform_loop_out()

p12_data$method <- factor(p12_data$method, levels = c("pw", "beale", "rating", 'composite'))

p12 <- plot_guts(p12_data)

p12

# stormflow w/ enrich data
p18_data <- loop_out %>%
    filter(flow == 'storm',
           cq == 'broken_dilution') %>%
    transform_loop_out()

p18_data$method <- factor(p18_data$method, levels = c("pw", "beale", "rating", 'composite'))

p18 <- plot_guts(p18_data)

p18

### make row 3 plots ####
# base flow w/ chemo data
p13_data <- loop_out %>%
    filter(flow == 'base',
           cq == 'chemostatic') %>%
    transform_loop_out()

p13_data$method <- factor(p13_data$method, levels = c("pw", "beale", "rating", 'composite'))

p13 <- plot_guts(p13_data)+
    theme(axis.title.y=element_text(size = 20))+
    labs(y = " ") +
    theme(axis.title.x=element_text(size = 20),
          axis.text.x=element_text(angle = 45, vjust = .5))+
    labs(x = " ")+
    scale_x_discrete(labels = c('LI', 'Beale', 'Rating', 'Composite'))
p13

# base flow w/ no pattern data
p14_data <- loop_out %>%
    filter(flow == 'base',
           cq == 'none') %>%
    transform_loop_out()

p14_data$method <- factor(p14_data$method, levels = c("pw", "beale", "rating", 'composite'))

p14 <- plot_guts(p14_data) +
    theme(axis.title.x=element_text(size = 20),
          axis.text.x=element_text(angle = 45, vjust = .5))+
    labs(x = 'Method') +
    scale_x_discrete(labels = c('LI', 'Beale', 'Rating', 'Composite'))
p14

# base flow w/ enrich data
p15_data <- loop_out %>%
    filter(flow == 'base',
           cq == 'enrich') %>%
    transform_loop_out()

p15_data$method <- factor(p15_data$method, levels = c("pw", "beale", "rating", 'composite'))

p15 <- plot_guts(p15_data) +
    theme(axis.title.x=element_text(size = 20),
          axis.text.x=element_text(angle = 45, vjust = .5))+
    labs(x = " ")+
    scale_x_discrete(labels = c('LI', 'Beale', 'Rating', 'Composite'))

p15

# base flow w/ broken diluiton data
p19_data <- loop_out %>%
    filter(flow == 'base',
           cq == 'broken_dilution') %>%
    transform_loop_out()

p19_data$method <- factor(p15_data$method, levels = c("pw", "beale", "rating", 'composite'))

p19 <- plot_guts(p19_data) +
    theme(axis.title.x=element_text(size = 20),
          axis.text.x=element_text(angle = 45, vjust = .5))+
    labs(x = " ")+
    scale_x_discrete(labels = c('LI', 'Beale', 'Rating', 'Composite'))

p19

## make mega fig ####
(legend+ theme(plot.margin = unit(c(0,30,0,0), "pt")) | p4 | p5 | p6 | p16)/
(p1+ theme(plot.margin = unit(c(0,30,0,0), "pt")) | p7 | p8 | p9 | p17)/
(p2+ theme(plot.margin = unit(c(0,30,0,0), "pt")) | p10 | p11 | p12 | p18)/
(p3+ theme(plot.margin = unit(c(0,30,0,0), "pt")) | p13 | p14 | p15 | p19)

#ggsave(filename = here('paper','ts simulation', 'pop.png'), width = 18, height = 9)
