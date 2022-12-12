## this can be used to get clean versions of the data for the PHASE 2 MODEL
library(here)
library(tidyverse)
library(feather)
library(data.table)
library(lubridate)
library(gridExtra)
library(grid)

# Step 1: READING IN FULL AND SIMPLE CHEM AND Q
###############################################
simple_chem_and_Q <- read_feather(here("paper","hbef_corr_exploration", "simple_chem_and_Q_v2.feather"))
full_chem_and_Q <- read_feather(here("paper","hbef_corr_exploration", "full_chem_and_Q_v2.feather"))

simple_chem_and_Q <- read_csv(here("paper","hbef_corr_exploration", "HBEFdata_All_2022-11-17.csv"))

# read in custom functions to remove rows with NA
# -----------------------------------------------

all_na <- function(x) any(!is.na(x)) #removes any with all NA
any_na <- function(x) !any(is.na(x)) # removes columns with partial NAS

`%notin%` <- Negate(`%in%`)

# Step 2: for SIMPLE dataframe with 3 sensor columns
####################################################

# establish which variables are in the list
# -----------------------------------------
grab_vars <- grep("GN",colnames(simple_chem_and_Q),value = T)
sensor_vars <- grep("IS",colnames(simple_chem_and_Q),value = T)
static_vars <- colnames(simple_chem_and_Q)[which(colnames(simple_chem_and_Q) %notin% c(grab_vars, sensor_vars, c("site_name", "date")))]

static_vars_pivoted <- simple_chem_and_Q %>% select(site_name, !! static_vars) #make dataframe of only static variables


# choose a variable to use
# ------------------------
univariate <- "GN_Ca"

# this output (merged output) is what goes into the pre-processing function
merged_input <- simple_chem_and_Q %>% select(site_name, date, !!sensor_vars, !!univariate) %>% filter(complete.cases(.)) %>% left_join(., static_vars_pivoted, by = "site_name") %>% select_if(any_na)%>% distinct()

# Step 3: for FULL dataframe with 3 sensor columns
####################################################

grab_vars <- grep("GN",colnames(full_chem_and_Q),value = T)
sensor_vars <- grep("IS",colnames(full_chem_and_Q),value = T)
static_vars <- colnames(full_chem_and_Q)[which(colnames(full_chem_and_Q) %notin% c(grab_vars, sensor_vars, c("site_name", "date")))]

static_vars_pivoted <- full_chem_and_Q %>% select(site_name, !! static_vars) #make dataframe of only static variables


# Step 4: Make correlation plots -- not necessary
##################################################

# choose a variable to use
# ------------------------
univariate <- "GN_Ca"
# custom function for correlations
make_corr_plots <- function(univariate) {

  merged_input <- full_chem_and_Q %>% select(site_name, date, !!sensor_vars, !!univariate) %>%
    filter(complete.cases(.)) %>% left_join(., static_vars_pivoted, by = "site_name") %>%
    select_if(any_na) %>% distinct()

  plot_vars <- grep("IS", colnames(merged_input), value = T)



  if("IS_spCond" %in% plot_vars){
  a <- ggplot(data = merged_input) +
    geom_point(aes(x =IS_spCond, y = !!as.name(univariate))) + theme_bw()}
  else {a <- ggplot(data = data.frame()) +
    geom_point() + xlim(0, 10) + ylim(0, 100) + theme_bw()}

  if("IS_temp" %in% plot_vars){
  b <- ggplot(data = merged_input) +
    geom_point(aes(x =IS_temp, y = !!as.name(univariate))) + theme_bw()}
  else {b <- ggplot(data = data.frame()) +
    geom_point() + xlim(0, 10) + ylim(0, 100) + theme_bw()}

  if("IS_discharge" %in% plot_vars){
  c <- ggplot(data = merged_input) +
              geom_point(aes(x =IS_discharge, y = !!as.name(univariate))) + theme_bw()}
  else {c<- ggplot(data = data.frame()) +
    geom_point() + xlim(0, 10) + ylim(0, 100) + theme_bw()}

  if("IS_DO" %in% plot_vars){
  d <- ggplot(data = merged_input) +
    geom_point(aes(x =IS_DO, y = !!as.name(univariate))) + theme_bw()}
  else {d <- ggplot(data = data.frame()) +
    geom_point() + xlim(0, 10) + ylim(0, 100) + theme_bw()}

  if("IS_NO3_N" %in% plot_vars){
  e <- ggplot(data = merged_input) +
    geom_point(aes(x =IS_NO3_N, y = !!as.name(univariate))) + theme_bw()}
  else {e <- ggplot(data = data.frame()) +
    geom_point() + xlim(0, 10) + ylim(0, 100) + theme_bw()}

  if("IS_pH" %in% plot_vars){
  f <- ggplot(data = merged_input) +
    geom_point(aes(x =IS_pH, y = !!as.name(univariate))) + theme_bw()}
  else {f <- ggplot(data = data.frame()) +
    geom_point() + xlim(0, 10) + ylim(0, 100) + theme_bw()}

  return(grid.arrange(grobs = list(a,b, c, d, e, f), top= textGrob(univariate)))

}


make_corr_plots("GN_Ca", 'w3') # for one solute

#already run:
# grep("GN", colnames(full_chem_and_Q), value = T)
# pdf("correlation_plots_v1.pdf", width = 5, height = 6)
# lapply(grep("GN", colnames(full_chem_and_Q), value = T), make_corr_plots)
# dev.off()

