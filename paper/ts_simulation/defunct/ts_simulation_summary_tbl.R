library(tidyverse)
library(kableExtra)

weekly <- read_csv(here('paper','ts simulation', 'weeklyFreq_100Reps.csv')) %>%
    mutate(freq = 'Weekly')
biweekly <- read_csv(here('paper','ts simulation', 'biweeklyFreq_100Reps.csv')) %>%
    mutate(freq = 'Biweekly')
monthly <- read_csv(here('paper','ts simulation', 'monthlyFreq_100Reps.csv')) %>%
    mutate(freq = 'Monthly')
loop_out <- rbind(weekly, biweekly, monthly) %>%
    mutate(freq = factor(freq, levels = c('Weekly', 'Biweekly', 'Monthly')))

transform_loop_out_tbl <- function(loop_out){
    pivot_wider(loop_out, names_from = method, values_from = estimate,
                #id_cols = c(runid, period),
                id_cols = c(runid, freq, flow, cq),
                values_fn = mean) %>%
        mutate(pw = ((pw-truth)/truth)*100,
               beale = ((beale - truth)/truth)*100,
               rating = ((rating-truth)/truth)*100,
               composite = ((composite - truth)/truth)*100) %>%
        select(-truth, -runid) %>%
        pivot_longer(cols = -c(freq, flow, cq),
                     #cols = everything() ,
                     names_to = 'method', values_to = 'error')
}


tbl_tbl <- loop_out %>%
    transform_loop_out_tbl() %>%
    group_by(method, cq, flow, freq) %>%
    summarize(mean = round(mean(error), digits = 2),
              sd = round(sd(error), digits = 2),
              quantile0 = round(quantile(error)[[1]], digits = 2),
              quantile25 = round(quantile(error)[[2]], digits = 2),
              quantile50 = round(quantile(error)[[3]], digits = 2),
              quantile75 = round(quantile(error)[[4]], digits = 2),
              quantile100 = round(quantile(error)[[5]], digits = 2))

tbl_tbl %>%
    filter(freq == 'Weekly') %>%
    select(-freq) %>%
    mutate(method = ifelse(method == 'beale', 'Beale', method),
           method = ifelse(method == 'pw', 'Linear Interpolation', method),
           method = ifelse(method == 'composite', 'Composite', method),
           method = ifelse(method == 'rating', 'Rating', method),
           cq = ifelse(cq == 'broken_dilution', 'Two-Part Dilution', cq),
           cq = ifelse(cq == 'enrich', 'Enriching', cq),
           cq = ifelse(cq == 'chemostatic', 'Chemostatic', cq),
           cq = ifelse(cq == 'none', 'No Pattern', cq),
           cq = ifelse(cq == 'broken_dilution', 'Two-Part Dilution', cq),
           flow = ifelse(flow == 'base', 'Baseflow', flow),
           flow = ifelse(flow == 'unaltered', 'Unaltered', flow),
           flow = ifelse(flow =='storm', 'Stormflow', flow)) %>%
    rename(Method = method,
           `C:Q` = cq,
           `Flow Regime` = flow,
           Mean = mean,
           SD = sd,
           Minimum = quantile0,
           `First Quartile` = quantile25 ,
           Median = quantile50 ,
           `Third Quartile` = quantile75,
           Maximum = quantile100) %>%
    mutate(Method = factor(Method, levels = c('Linear Interpolation', 'Beale', 'Rating', 'Composite'))) %>%
    arrange(`Flow Regime`, `C:Q`, Method) %>%
    kbl() %>%
    kable_classic_2(full_width = F, bootstrap_options = c("striped"))
