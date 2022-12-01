library(tidyverse)
library(feather)
library(here)



lo <- read_feather('data/ms/hbef/stream_chemistry/w3.feather') %>%
    rename(date = datetime)
hi <- read_feather('C:/Users/gubbi/desktop/w3_sensor_wdisch.feather') %>%
    select(date, IS_discharge, IS_spCond) %>%
    group_by(date) %>%
    summarize(q = mean(IS_discharge),
              sc = mean(IS_spCond),
              date = date) %>%
    unique()

dat <- left_join(hi, lo, by = 'date')

ggplot(dat, aes(x = sc, y = val)) +
           geom_point()+
           facet_wrap(~var, scales = 'free')

for(i in 1:length(unique(dat$var))){
solute = unique(dat$var)[i]

loop_dat <- dat %>%
        filter(var == solute)

fit <- summary(lm(val~sc, data = loop_dat))

ggplot(loop_dat, aes(x = sc, y = val))+
    geom_point()+
    geom_smooth(method = 'lm')+
    labs(x = 'SC',
         y = 'Concentration',
         title = solute,
         caption = paste('r-squared of ', round(fit$r.squared[[1]], digits = 3)))+
    theme(text = element_text(size = 20))

ggsave(here('paper/spcond investigation', paste(solute,'.png')))

}
