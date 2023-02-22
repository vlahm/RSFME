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