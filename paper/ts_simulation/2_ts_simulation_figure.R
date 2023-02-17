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
# read in data ######
d <- read_feather('C:/Users/gubbi/desktop/w3_sensor_wdisch.feather') %>%
    mutate(wy = water_year(datetime, origin = 'usgs'))

## Subset to 2016 wy
target_wy <- 2016
dn <- d %>%
    filter(wy == target_wy)
## read in output from 1_ts_simulation_analysis.R####
weekly <- read_csv(here('paper','ts_simulation', 'weeklyFreq_100Reps20221221.csv')) %>%
  mutate(freq = 'Weekly')
biweekly <- read_csv(here('paper','ts_simulation', 'biweeklyFreq_100Reps20221221.csv')) %>%
  mutate(freq = 'Biweekly')
monthly <- read_csv(here('paper','ts_simulation', 'monthlyFreq_100Reps20221221.csv')) %>%
  mutate(freq = 'Monthly')
loop_out <- rbind(weekly, biweekly, monthly) %>%
  mutate(freq = factor(freq, levels = c('Weekly', 'Biweekly', 'Monthly')))

## load simulated_series from 1_ts_simulation_analysis.R ####
load(here('paper','ts_simulation', 'simulated_series.Rdata'))

### make header plots #####
side_ymin <- 0.01
side_ymax <- 1000
side_breaks <- c(1e-1, 1e1, 1e3)
side_labels <- c('0.1', '10', '1,000')
x_axis_text_size = 25
y_axis_text_size = 25
# unaltered q plot
p1 <- ggplot(dn, aes(x = date))+
  geom_line(aes(y = simulated_series[[1]])) +
  theme_classic()+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.title.y=element_blank(),
        text = element_text(size = 20),
        axis.text.y = element_text(size = y_axis_text_size))+
  labs(title = 'Unaltered Flow')+
  scale_y_log10(limits = c(side_ymin,side_ymax),
                breaks = side_breaks,
                labels = side_labels)

p1

# storm q plot
p2 <- ggplot(dn, aes(x = date))+
  geom_line(aes(y = simulated_series[[2]])) +
  theme_classic()+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        text = element_text(size = 20),
        axis.text.y = element_text(size = y_axis_text_size),
        axis.title.y = element_text(size = y_axis_text_size))+
  labs(title = 'Stormflow Dominated',
       y = 'Q (lps)')+
  scale_y_log10(limits = c(side_ymin,side_ymax),
                breaks = side_breaks,
                labels = side_labels)
p2

# base q plot
p3 <- ggplot(dn, aes(x = date))+
  geom_line(aes(y = simulated_series[[3]])) +
  theme_classic()+
  theme(axis.title.y=element_blank(),
        axis.title.x = element_blank(),
        text = element_text(size = 20),
        axis.text.x = element_text(angle = 45, hjust = 1, size = x_axis_text_size),
        axis.text.y = element_text(size = y_axis_text_size)
  )+
  labs(title = 'Baseflow Dominated')+
  scale_y_log10(limits = c(side_ymin,side_ymax),
                breaks = side_breaks,
                labels = side_labels)+
  scale_x_continuous(breaks = c(as_date('2015-10-01'), as_date('2016-04-01'), as_date('2016-10-01')),
                     labels = c('10/2015', '4/2016', '10/2016'))


p3

# set common limits to top row graphs
top_row_breaks <- c(1e-2, 1, 1e2)
top_row_labels <- c('0.01', '1', '100')
top_row_ymax <- 1e2
top_row_ymin <- 1e-2
# chemo cq
p4 <- tibble(q = simulated_series[[1]], con = simulated_series[[4]]) %>%
  ggplot(aes(x = q, y = con)) +
  geom_point() +
  theme_classic()+
  scale_x_log10(breaks = side_breaks,
                labels = side_labels) +
  scale_y_log10(limits = c(top_row_ymin, top_row_ymax),
                breaks = top_row_breaks,
                labels = top_row_labels) +
  labs(title = 'Chemostatic',
       y = 'C (mg/L)')+
  theme(axis.title.x=element_blank(),
        text = element_text(size = 20),
        axis.text.y = element_text(size = y_axis_text_size),
        axis.text.x = element_text(size = x_axis_text_size),
        axis.title.y = element_text(size = y_axis_text_size))

p4

# no pattern cq
p5 <- tibble(q = simulated_series[[1]], con = simulated_series[[5]]) %>%
  ggplot(aes(x = q, y = con)) +
  geom_point() +
  theme_classic()+
  scale_x_log10(breaks = side_breaks,
                labels = side_labels) +
  scale_y_log10(limits = c(top_row_ymin, top_row_ymax),
                breaks = top_row_breaks,
                labels = top_row_labels)+
  labs(title = 'No Pattern',
       x = 'Q (lps)')+
  theme(axis.title.y=element_blank(),
        text = element_text(size = 20),
        axis.text.y = element_text(size = y_axis_text_size),
        axis.text.x = element_text(size = x_axis_text_size),
        axis.title.x = element_text(size = x_axis_text_size))
p5

# enrich cq
p6 <- tibble(q = simulated_series[[1]], con = simulated_series[[6]]) %>%
  ggplot(aes(x = q, y = con)) +
  geom_point() +
  theme_classic()+
  scale_x_log10(breaks = side_breaks,
                labels = side_labels) +
  scale_y_log10(limits = c(top_row_ymin, top_row_ymax),
                breaks = top_row_breaks,
                labels = top_row_labels)+
  labs(title = 'Enriching',
       x = 'Q')+
  theme(axis.title.y=element_blank(),
        axis.title.x=element_blank(),
        text = element_text(size = 20),
        axis.text.y = element_text(size = y_axis_text_size),
        axis.text.x = element_text(size = x_axis_text_size))
p6

# dilute cq
p16 <- tibble(q = simulated_series[[1]], con = simulated_series[[7]]) %>%
  ggplot(aes(x = q, y = con)) +
  geom_point() +
  theme_classic()+
  scale_x_log10(breaks = side_breaks,
                labels = side_labels) +
  scale_y_log10(limits = c(top_row_ymin, top_row_ymax),
                breaks = top_row_breaks,
                labels = top_row_labels)+
  labs(title = 'Dilution',
       x = 'Q')+
  theme(axis.title.y=element_blank(),
        axis.title.x=element_blank(),
        text = element_text(size = 20),
        axis.text.y = element_text(size = y_axis_text_size),
        axis.text.x = element_text(size = x_axis_text_size))
p16

# broken dilution c:q

### make row 1 plots ####
ymin = -75
ymax = 75

plot_guts <- function(p){
  ggplot(p, aes(x = method, y = error))+
    geom_hline(yintercept = 0)+
    geom_boxplot(
      aes(fill = freq)
    )+
    theme_classic()+
    theme(
      axis.title.x=element_blank(),
      axis.text.x=element_blank(),
      axis.title.y=element_blank(),
      text = element_text(size = 20),
      legend.position="none",
      axis.text.y = element_text(size = y_axis_text_size)
    )+
    scale_fill_manual(values = c('darkorange', 'gray', 'deepskyblue'))+
    ylim(ymin, ymax)
}

transform_loop_out <- function(loop_out){
  pivot_wider(loop_out, names_from = method, values_from = estimate,
              #id_cols = c(runid, period),
              id_cols = c(runid, freq),
              values_fn = mean) %>%
    mutate(pw = ((pw-truth)/truth)*100,
           beale = ((beale - truth)/truth)*100,
           rating = ((rating-truth)/truth)*100,
           composite = ((composite - truth)/truth)*100) %>%
    select(-truth, -runid) %>%
    pivot_longer(cols = -freq,
                 #cols = everything() ,
                 names_to = 'method', values_to = 'error')
}
#  spacer plot that is only legend
p0_data <- loop_out %>%
  transform_loop_out()


p0_data$method <- factor(p0_data$method, levels = c("pw", "beale", "rating", 'composite'))

p0 <-plot_guts(p0_data)+
  theme(legend.position = 'left',
        legend.title = element_text(size = 30),
        legend.text = element_text(size = 25),
        legend.key.size = unit(2.5, 'cm')
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
  theme(axis.title.y=element_text(size = 20),
        axis.text.y = element_text(size = y_axis_text_size))+
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
         cq == 'dilution') %>%
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
  theme(axis.title.y=element_text(angle = 90, size = y_axis_text_size))+
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

p12_data$error[p18_data$method == 'pw' | p12_data$method == 'beale'] <- NA
p12_data$method <- factor(p12_data$method, levels = c("pw", "beale", "rating", 'composite'))

p12 <- plot_guts(p12_data)

p12

# stormflow w/ enrich data
p18_data <- loop_out %>%
  filter(flow == 'storm',
         cq == 'dilution') %>%
  transform_loop_out()
p18_data$error[p18_data$method == 'pw' | p18_data$method == 'beale'] <- NA
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
  theme(axis.title.y=element_text(size = y_axis_text_size))+
  labs(y = " ") +
  theme(axis.title.x=element_text(size = x_axis_text_size),
        axis.text.x=element_text(angle = 45, vjust = .9, size = x_axis_text_size, hjust = 1))+
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
  theme(axis.title.x=element_text(size = 30),
        axis.text.x=element_text(angle = 45, vjust = .9, size = x_axis_text_size, hjust = 1))+
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
  theme(axis.title.x=element_text(size = x_axis_text_size),
        axis.text.x=element_text(angle = 45, vjust = .9, size = x_axis_text_size, hjust = 1))+
  labs(x = " ")+
  scale_x_discrete(labels = c('LI', 'Beale', 'Rating', 'Composite'))

p15

# base flow w/ broken diluiton data
p19_data <- loop_out %>%
  filter(flow == 'base',
         cq == 'dilution') %>%
  transform_loop_out()

p19_data$method <- factor(p15_data$method, levels = c("pw", "beale", "rating", 'composite'))

p19 <- plot_guts(p19_data) +
  theme(axis.title.x=element_text(size = x_axis_text_size),
        axis.text.x=element_text(angle = 45, vjust = .9, size = x_axis_text_size, hjust = 1))+
  labs(x = " ")+
  scale_x_discrete(labels = c('LI', 'Beale', 'Rating', 'Composite'))

p19

# make mega fig ####
(legend+ theme(plot.margin = unit(c(0,30,0,0), "pt")) | p4 | p5 | p6 | p16)/
  (p1+ theme(plot.margin = unit(c(0,30,0,0), "pt")) | p7 | p8 | p9 | p17)/
  (p2+ theme(plot.margin = unit(c(0,30,0,0), "pt")) | p10 | p11 | p12 | p18)/
  (p3+ theme(plot.margin = unit(c(0,30,0,0), "pt")) | p13 | p14 | p15 | p19)

ggsave(filename = here('paper','ts_simulation', 'pop_test.png'), width = 20, height = 16)

# make supp table 1 ####
pretty_flow <- tibble(flow = c('base', 'storm', 'unaltered'),
                      flow_pretty = c('Baseflow', 'Stormflow', 'Unaltered'))
pretty_cq <- tibble(cq = c('dilution', 'enrich', 'none', 'chemostatic'),
                    cq_pretty = c('Diluting', 'Enriching', 'No-Pattern', 'Chemostatic'))
pretty_method = tibble(method = c('pw', 'beale', 'rating', 'composite'),
                       method_pretty = c('Linear Interpolation', 'Beale', 'Rating', 'Composite'))



supp_tbl_pre <- loop_out %>%
  pivot_wider(names_from = method, values_from = estimate,
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
               names_to = 'method', values_to = 'error') %>%
  group_by(flow, cq, freq, method) %>%
  summarize(mean_error = round(mean(error), digits = 2),
            sd_error = round(sd(error), digits = 2),
            ci_95_low = round(mean_error-(sd_error*1.96), digits = 2),
            ci_95_hi = round(mean_error+(sd_error*1.96), digits = 2),
            min_error = round(min(error), digits = 2),
            max_error = round(max(error), digits = 2),
            iqr = IQR(error),
            median = round(median(error), digits = 2),
            n_outliers = sum((error > median+(iqr*1.5))+(error < median-(iqr*1.5))),
            ci = paste0(ci_95_low,', ', ci_95_hi)
  ) %>%
  left_join(.,pretty_cq, by = 'cq') %>%
  left_join(.,pretty_method, by = 'method') %>%
  left_join(.,pretty_flow, by = 'flow') %>%
  ungroup()%>%
  arrange(flow, cq, freq, method) %>%
  #filter(flow == target_flow) %>%
  select(`Flow Regime` = flow_pretty, `C:Q` = cq_pretty, Frequency = freq,
         Method = method_pretty, Mean = mean_error, SD = sd_error,
         `95% CI` = ci, Median = median, Minimum = min_error, Maximum = max_error, Outliers = n_outliers)

write_csv(supp_tbl_pre, file = here('paper','ts_simulation', 'supp_table.csv'))

