# function to apply all methods #####
apply_methods_coarse <- function(chem_df, q_df){
    n = nrow(chem_df)

    out <- tibble(method = as.character(), estimate = as.numeric())#, se)
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
    return(out)
}

# create function to take every nth element from series. ####
nth_element <- function(vector, starting_position, n) {
    vector[seq(starting_position, length(vector), n)]
}