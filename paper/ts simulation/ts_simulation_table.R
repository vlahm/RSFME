weekly <- read_csv(here('paper','ts simulation', 'weeklyFreq_100Reps.csv')) %>%
    mutate(freq = 'Weekly')
biweekly <- read_csv(here('paper','ts simulation', 'biweeklyFreq_100Reps.csv')) %>%
    mutate(freq = 'Biweekly')
monthly <- read_csv(here('paper','ts simulation', 'monthlyFreq_100Reps.csv')) %>%
    mutate(freq = 'Monthly')
loop_out <- rbind(weekly, biweekly, monthly) %>%
    mutate(freq = factor(freq, levels = c('Weekly', 'Biweekly', 'Monthly')))

tbl_tbl <- loop_out %>%
    f