library(tidyverse)
library(feather)
library(here)
library(patchwork)

out_tbl <- tibble(domain = as.character(),
                  site = as.character(),
                  solute = as.character(),
                  year = as.integer(),
                  ms_rec = as.numeric(),
                  max_val = as.numeric(),
                  min_val = as.numeric(),
                  mean_val = as.numeric(),
                  rec_method = as.character())

domain_list <- list.files(here('paper', 'ms_application', 'ms_flux_12162022'))
domain_tbl <- tibble(domains = domain_list) %>%
    filter(!str_detect(domains, '.csv'),
           !str_detect(domains, '.feather'))
for(d in domain_tbl$domains){
#d <- domain_tbl$domains[10]
site_list <- list.files(here('paper', 'ms_application', 'ms_flux_12162022', d, 'stream_flux')) %>%
    tibble(sites = .) %>%
    filter(str_detect(sites, '.feather'))
for(s in unique(site_list$sites)){
#s <- site_list$sites[1]
site_dat <- read_feather(here('paper', 'ms_application', 'ms_flux_12162022', d, 'stream_flux', s)) %>%
    filter(method != 'wrtds',
           method != 'failed'#,
           #method != 'average'
           )
if(nrow(site_dat) == 0){next}else{

for(i in unique(site_dat$var)){
    #i = unique(site_dat$var)[1]

    solute_dat <- site_dat %>%
        filter(var == i)

for(j in unique(solute_dat$wy)){
        #j = unique(solute_dat$wy)[1]

        solute_year_dat <- solute_dat %>%
            filter(wy == j)

        ms_rec <- solute_year_dat$val[solute_year_dat$ms_recommended == 1]
        max_val <- max(solute_year_dat$val)
        min_val <- min(solute_year_dat$val)
        mean_val <- mean(solute_year_dat$val)
        rec_method <- solute_year_dat$method[solute_year_dat$ms_recommended == 1]


        loop_tbl <- tibble(domain = d,
               site = s,
               solute = i,
               year = j,
               ms_rec = ms_rec,
               max_val = max_val,
               min_val = min_val,
               mean_val = mean_val,
               rec_method = rec_method)

        out_tbl <- rbind(out_tbl, loop_tbl)

           } # end year
        } # end solute
}# end else
    } # end site
} # end domain

# make method range plot #####
n_frame <- out_tbl %>%
    mutate(wcs_max = ((max_val-ms_rec)/((ms_rec+max_val)*.5)*100),
            wcs_min = ((min_val-ms_rec)/((ms_rec+min_val)*.5)*100),
           pct_range = ((max_val-min_val)/((min_val+max_val)*.5)*100),
           wcs_range = wcs_max-wcs_min) %>%
    filter(solute == 'GN_NO3_N' |
           solute == 'GN_Ca')

mean(test$wcs_max)

mean(test$wcs_min)

ggplot(test, aes(x = rec_method)) +
    geom_bar()

p_box <- n_frame %>%
    select(solute, pct_range, wcs_range) %>%
    pivot_longer(cols = -solute, names_to = 'var', values_to = 'val') %>%
    filter(var == 'pct_range') %>%
    ggplot(aes(x = solute, y = val, color = solute))+
    geom_boxplot(size = 2)+
    scale_y_log10()+
    labs(x = 'Solute',
         y = 'Method Range (%)')+
    scale_x_discrete(labels = c('Calcium', 'Nitrate'))+
    scale_color_manual(labels = c('Calcium', 'Nitrate (as N)'),
                       values = c('red', 'blue'))+
    theme_classic()+
    theme(text = element_text(size = 20),
          legend.position = 0)


#ggsave(file = here('paper', 'ms_application', 'method_comp_fig.png'))

p_den <- out_tbl  %>%
     filter(solute == 'GN_NO3_N' |
                solute == 'GN_Ca') %>%
     select(ms_rec, solute) %>%
     ggplot(aes(x = ms_rec, color = solute))+
     geom_density(size = 2)+
    scale_x_log10(breaks = c(1e-2, 1, 1e2),
                  labels = c('0.01', '1', '100'))+
    scale_color_manual(labels = c('Calcium', 'Nitrate (as N)'),
                         values = c('red', 'blue'))+
    labs(x = 'Load (kg/ha/year)',
         y = 'Density',
         color = 'Solute')+
    theme_classic()+
    theme(text = element_text(size = 20))
ggsave(file = here('paper', 'ms_application', 'method_comp_den.png'))

p_den|p_box
ggsave(file = here('paper', 'ms_application', 'method_comp_combined.png'), width = 8, height = 6)
