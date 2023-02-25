library(tidyverse)
library(feather)
library(here)
library(patchwork)
library(ggthemes)

# initialize
out_tbl <- tibble(domain = as.character(),
                  site = as.character(),
                  solute = as.character(),
                  year = as.integer(),
                  ms_rec = as.numeric(),
                  max_val = as.numeric(),
                  min_val = as.numeric(),
                  mean_val = as.numeric(),
                  rec_method = as.character())

# gather domains
domain_list <- list.files(here('paper', 'macrosheds_application', 'ms_flux_12162022'))
domain_tbl <- tibble(domains = domain_list) %>%
    filter(!str_detect(domains, '.csv'),
           !str_detect(domains, '.feather'))

# loop through domains
for(d in domain_tbl$domains){
    site_list <- list.files(here('paper', 'macrosheds_application', 'ms_flux_12162022', d, 'stream_flux')) %>%
        tibble(sites = .) %>%
        filter(str_detect(sites, '.feather'))

    # loop through sites
    for(s in unique(site_list$sites)){
        site_dat <- read_feather(here('paper', 'macrosheds_application', 'ms_flux_12162022', d, 'stream_flux', s)) %>%
            filter(method != 'wrtds',
                   method != 'failed'#,
                   #method != 'average'
            )
        if(nrow(site_dat) == 0){next}else{
            # loop through variables
            for(i in unique(site_dat$var)){
                #i = unique(site_dat$var)[1]

                solute_dat <- site_dat %>%
                    filter(var == i)
                # loop through years
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

label_tbl <- tibble(solute = c('GN_Ca', 'GN_NO3_N'),
                    label = c('Calcium',
                              'Nitrate (as N)'))

out_tbl %>%
    filter(solute %in% c('GN_Ca', 'GN_NO3_N'),
           ms_rec > 0) %>% #need to fix this in rec pipeline
    left_join(., label_tbl, by = 'solute') %>%
    ggplot()+
        geom_histogram(aes(x = ms_rec, fill = solute), color = 'black', bins = 50)+
        scale_x_log10(breaks = c(0.01, 1, 1000),
                      labels = c('0.01', '1', '1000'))+
        facet_wrap(~label, ncol = 1)+
    labs(x = 'Load (kg/ha/year, log)',
         y = 'Count',
         title = 'Macrosheds Load Estimates')+
    scale_fill_manual(values = c('red', 'blue'))+
    theme_minimal()+
    theme(legend.position = 'none',
          text=element_text(size=20))
ggsave(file = here('paper', 'macrosheds_application', 'descriptive_hist.png'), width = 8, height = 6)

