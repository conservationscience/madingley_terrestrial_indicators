
## REPOSITORY: https://github.com/conservationscience/madingley_terrestrial_indicators

rm(list = ls())

# Directory path to git repo

# Simone's deakin laptop

# cd "C:\\Users\\ssteven\\OneDrive - Deakin University\\Deakin\\Chapter_3_indicator_testing\\madingley_terrestrial_indicators"

# Data structure

#' Input data is structured hierarchically as follows:
#' # Location (1)
#' ## Scenarios (3)
#' ### Replicates (25)
#' 
#' Output indicator data is calculated at the same level, so each indicator should
#' have values for:
#' # Location (1)
#' ## Scenarios (3)
#' ### Replicates (25) 


# TODO LIST ----

# Libraries ----

## Data wrangling
library(tidyverse)
library(tidylog)

## Analysis
library(mgcv)
library(changepoint)

# Functions ----

#' Description

#' @param 
#' @param 
#' @param 
#' @return 

function_name <- function(param){

}

input <- data

smooth_gam <- function(input, smoothing_param) {
  
  input <- input %>% 
    mutate(abundance = abundance +
             (mean(abundance)*0.01))
  
  # mod <- gam(log10(abundance) ~ s(annual_time_step), sp = smoothing_param,  
  #            data = input, method = "REML")
  
  mod <- gam(log10(abundance) ~ s(annual_time_step, k = 4), sp = smoothing_param,  
             data = input, method = "REML")
  
  # plot(mod)
  # summary(mod)
  # gam.check(mod)
  
  modelled_abundance <- predict(mod, pred_annual_time_step = unique(input$annual_time_step))
  
  modelled_abundance <- as.data.frame(cbind(annual_time_step =unique(input$annual_time_step),
                                            transformed_abundance = exp(modelled_abundance)))
  
  newdata <- input %>% 
    merge(modelled_abundance, by = "annual_time_step")
  
  return(newdata)
  
  
}
# Set up paths ----

if (Sys.info()['nodename'] == "SIMONE-PC") {
  
   IndicatorsProject <- "N:/Quantitative-Ecology/Indicators-Project"
  
}  

if (Sys.info()['nodename'] == "ANALYTIX2") {
  
  SourceToModels <- "C:/Users/ssteven/Desktop/Serengeti-analysis"
  IndicatorsProject <- "N:/Quantitative-Ecology/Indicators-Project"
}

if (Sys.info()['nodename'] == "20FMPC0C6GH9") {
  
  SourceToModels <- "C:/Users/ssteven/Dropbox/Deakin/Serengeti-analysis"
  IndicatorsProject <- "N:/Quantitative-Ecology/Indicators-Project"
}

if (Sys.info()['nodename'] == "80VVPF0SSB53") {
  
  SourceToModels <- 'K:/BiodiversityIndicators/serengeti'
  IndicatorsProject <- "N:/Quantitative-Ecology/Indicators-Project"
}

# Inputs ----

# Get date to label outputs

today <- Sys.Date()

# Define mode (development == TRUE will only perform operations on a small subset
# of folders, not all outputs)

development_mode <- FALSE

# Specify universal input arguments

if (development_mode == FALSE) {
  
  impact_start <- 1100
  impact_end <- 1200
  burnin_months <- 1000*12 # in months
  n <- 12
  numboots <- 1000 # Rowland et al 2021 (uncertainty)
  start_time_step <- 1
  gen_timeframe <- 10 * 12
  interval <- 12
  
} else {
  
  impact_start <- 10
  impact_end <- 20
  burnin_months <- 1*12 # in months
  n <- 1 
  numboots <- 1000
  start_time_step <- 1
  gen_timeframe <- 10
  interval <- 2
  
}

indicators_project <- IndicatorsProject # File path for entire project directory

location <- 'Serengeti'

scenarios <- list("baseline", "land use", 
                  "carnivore harvesting", 
                  "herbivore harvesting")

# Set up output folders

analysis_inputs_folder <- file.path(indicators_project, 
                           "/Serengeti/Outputs_from_analysis_code/Analysis_inputs")

if( !dir.exists( file.path(analysis_inputs_folder) ) ) {
  dir.create( file.path(analysis_inputs_folder), recursive = TRUE )
  
}

analysis_outputs_folder <- file.path(indicators_project, 
                                     "/Serengeti/Outputs_from_analysis_code/Analysis_outputs")

if( !dir.exists( file.path(analysis_outputs_folder) ) ) {
  dir.create( file.path(analysis_outputs_folder), recursive = TRUE )
  
}

analysis_plots_folder <- file.path(indicators_project, 
                                     "/Serengeti/Outputs_from_analysis_code/Analysis_")

if( !dir.exists( file.path(analysis_plots_folder) ) ) {
  dir.create( file.path(analysis_plots_folder), recursive = TRUE )
  
}

# Load scenario level data ----

if (development_mode == TRUE) {
  
#indicators_all <- readRDS("N:/Quantitative-Ecology/Indicators-Project/Serengeti/Outputs_from_indicator_code/Indicator_outputs/2021-07-12_all_indicators_output_data.rds")
indicators_all <- readRDS("N:/Quantitative-Ecology/Indicators-Project/Serengeti/Outputs_from_indicator_code/Indicator_outputs/2021-07-13_all_indicators_output_data_list.rds")

} else if (development_mode == FALSE) {
  
indicators_all <- readRDS("N:/Quantitative-Ecology/Indicators-Project/Serengeti/Outputs_from_indicator_code/Indicator_outputs/2021-08-19_all_indicators_output_data_reps_averaged_list2.rds")

}

# Split the list by indicators

rli <- indicators_all[["RLI"]]

lpi <- indicators_all[["LPI"]]

harvested <- indicators_all[["total abundance harvested"]]


# CORRELATION ----

lpi_landuse <- lpi[[2]]
rli_landuse <- rli[[2]]
harvested_landuse <- harvested[[2]]

head(lpi_landuse)
dim(lpi_landuse)

cor(lpi_landuse$indicator_score, 
    harvested_landuse$indicator_score, method = "spearman")

cor(rli_landuse$indicator_score, harvested_landuse$indicator_score, 
    method = "spearman")

indicator_cor_scores <- list()

for (i in seq_along(indicators_all)) {
  
  # Get the indicator (all scenarios)
  
  indicator_scenarios <- indicators_all[[i]]
  
  print(indicator_scenarios[[1]]$indicator[1])
  
  harvested_scenarios <- indicators_all[["total abundance harvested"]]
  
  # Make a list to hold scenario correlation coefficients
  
  scenario_cor_scores <- list()
  
  for (j in seq_along(indicator_scenarios)) {
    
  print(scenarios[i])
    
  indicator <- indicator_scenarios[[j]] %>% 
               dplyr::select(annual_time_step, indicator_score) %>% 
               rename(indicator = indicator_score)
  
  harvested<- harvested_scenarios[[j]] %>% 
              dplyr::select(annual_time_step, indicator_score) %>% 
              rename(harvested = indicator_score)
  
  comparison <- indicator %>% 
                merge(harvested, by = "annual_time_step")
    
  cor <- cor(comparison$indicator, 
             comparison$harvested, method = "spearman")
  
  scenario_cor_scores[[j]] <- data.frame(indicator = indicator_scenarios[[j]]$indicator[1],
                                         correlation = cor,
                                         scenario = scenarios[[j]])
  
  }
  
  scenario_cor_df <- do.call(rbind, scenario_cor_scores)
  
  indicator_cor_scores[[i]] <- scenario_cor_df

}

correlation_dataframe <- do.call(rbind, indicator_cor_scores) %>% 
                         arrange(cor)

# BREAK POINT ANALYSIS ----


scenario_indicator_names <- names(indicators_all)

scenario_changepoint_summaries <- list()

for (i in seq_along(indicators_all)) {

# Get a single indicator data
  
single_indicator <- indicators_all[[i]]

indicator_changepoint_summaries <- list()

  for (j in seq_along(single_indicator)) {

  # Get a single scenario for the indicator
    
  data <- single_indicator[[j]]
    
  # Convert into a time series
  
  indicator_ts <- as.ts(data$indicator_score)
  
  # Calculate change points
  
  cpt <- cpt.mean(indicator_ts, method = "PELT", 
                   penalty = "AIC", 
                   pen.value = c(1,25))
  
  # Save break points
  
  indicator_changepoint_summaries[[j]] <- cpt
  
  }

scenario_changepoint_summaries[[i]] <- indicator_changepoint_summaries

}

plot(scenario_changepoint_summaries[[23]][[3]])
summary(scenario_changepoint_summaries[[23]][[3]])


# GAM ----

# Prepare a GAM for one indicator, for the land use scenario

harvested <- indicators_all[[23]][[2]] %>% 
              rename(autotrophs = indicator_score)

head(harvested)

input <- indicators_all[["LPI"]][[2]] %>%
         # Add a disturbance variable
         mutate(disturbance = ifelse(annual_time_step < 100, 0,
                              ifelse(annual_time_step > 99 & 
                                       annual_time_step < 200, 1, 0))) %>% 
         merge(harvested[c("annual_time_step", "autotrophs")], 
              by = "annual_time_step") %>% 
         mutate(auto_scaled = scale(autotrophs),
                lpi_scaled = scale(indicator_score))
head(input)

# mod <- gam(log10(abundance) ~ s(annual_time_step), sp = smoothing_param,  
#            data = input, method = "REML")

mod <- gam(lpi_scaled ~ s(auto_scaled, k = 20) + disturbance, 
           sp = 0.001,  
           data = input, method = "REML")

plot(mod)
summary(mod)
gam.check(mod)

modelled_abundance <- predict(mod, pred_annual_time_step = unique(input$annual_time_step))
head(modelled_abundance)

modelled_abundance <- as.data.frame(cbind(annual_time_step =unique(input$annual_time_step),
                                          transformed_abundance = modelled_abundance))

newdata <- input %>% 
  merge(modelled_abundance, by = "annual_time_step")

head(newdata)

ggplot(data = newdata) +
  geom_line(aes(x = annual_time_step, y = transformed_abundance)) +
  geom_line(aes(x = annual_time_step, y = lpi_scaled), col = "deep pink") 


# DERIVATIVES ----

# RED LIST INDEX ----

# Create the spline functions and calculate derivatives (1 per replicate)

rli_scenario_derivatives <- list()

for (i in seq_along(rli)) {
  
  # Get the outputs for a single scenario and create a spline function
    
    rli_spline <- splinefun(x = rli[[i]]$annual_time_step, 
                            y = rli[[i]]$indicator_score)
    
    # Use the spline function to calculate the first and second derivatives and
    # add them to our indicator score dataframes
    
    rli_scenario_derivatives[[i]] <- rli[[i]] %>% 
         mutate(first_derivative = rli_spline(seq(1, nrow(rli[[i]]), 1), 
                                              deriv = 1),
                second_derivative = rli_spline(seq(1, nrow(rli[[i]]), 1), 
                                               deriv = 2))
                                  
}

# Test plot RLI

rli_test <- rli_scenario_derivatives[[2]]
head(rli_test)

rli_annual <- ggplot(data = rli_test) +
  geom_line(aes(x = annual_time_step, y = first_derivative)) 
  #geom_line(aes(x = annual_time_step, y = second_derivative), col = "hot pink")

rli_annual

ggsave(file.path(analysis_outputs_folder, "rli_annual_2nd_dervs.png"),
       rli_annual,  device = "png")

# LIVING PLANET INDEX ----

# Create the spline functions and calculate derivatives (1 per replicate)

lpi_scenario_derivatives <- list()
lpi_replicate_derivatives <- list()

for (i in seq_along(lpi)) {
  
  # Get the outputs for a single scenario
  
  scenario <- lpi[[i]] 
  
  for (j in seq_along(scenario)) {
    
    # Get the outputs for a single replicate and create a spline function
    
    lpi_replicate_spline <- splinefun(x = scenario[[j]]$annual_time_step, 
                                      y = scenario[[j]]$indicator_score)
    
    # Use the spline function to calculate the first and second derivatives and
    # add them to our indicator score dataframes
    
    lpi_replicate_derivatives[[j]] <- scenario[[j]] %>% 
      mutate(first_derivative = lpi_replicate_spline(seq(1, 
                                                         nrow(scenario[[j]]), 1), 
                                                     deriv = 1),
             second_derivative = lpi_replicate_spline(seq(1, 
                                                          nrow(scenario[[j]]), 1), 
                                                      deriv = 2))
    
  }
  
  lpi_scenario_derivatives[[i]] <- lpi_replicate_derivatives
  
}

# Test plot LPI

lpi_test <- lpi_scenario_derivatives[[1]][[1]]
head(lpi_test)

lpi_annual <- ggplot(data = lpi_test) +
  #geom_line(aes(x = annual_time_step, y = first_derivative)) +
  geom_line(aes(x = annual_time_step, y = second_derivative), col = "hot pink")

ggsave(file.path(analysis_outputs_folder, "lpi_annual_2nd_dervs.png"),
       lpi_annual,  device = "png")

# Load  5 yr data ----

if (development_mode == TRUE) {
  
  #indicators_all <- readRDS("N:/Quantitative-Ecology/Indicators-Project/Serengeti/Outputs_from_indicator_code/Indicator_outputs/2021-07-12_all_indicators_output_data.rds")
  indicators_all <- readRDS("N:/Quantitative-Ecology/Indicators-Project/Serengeti/Outputs_from_indicator_code/Indicator_outputs/2021-07-13_all_indicators_output_data_list.rds")
  
} else if (development_mode == FALSE) {
  
  indicators_all <- readRDS("N:/Quantitative-Ecology/Indicators-Project/Serengeti/Outputs_from_indicator_code/Indicator_outputs/2021-08-04_all_indicators_output_data_list_5yr.rds")
  
}

# Split the list by indicators

rli <- indicators_all[["RLI"]]

lpi <- indicators_all[["LPI"]]

# RED LIST INDEX ----

# Create the spline functions and calculate derivatives (1 per replicate)

rli_scenario_derivatives <- list()
rli_replicate_derivatives <- list()

for (i in seq_along(rli)) {
  
  # Get the outputs for a single scenario
  
  scenario <- rli[[i]] 
  
  for (j in seq_along(scenario)) {
    
    # Get the outputs for a single replicate and create a spline function
    
    rli_replicate_spline <- splinefun(x = scenario[[j]]$annual_time_step, 
                                      y = scenario[[j]]$indicator_score)
    
    # Use the spline function to calculate the first and second derivatives and
    # add them to our indicator score dataframes
    
    rli_replicate_derivatives[[j]] <- scenario[[j]] %>% 
      mutate(first_derivative = rli_replicate_spline(seq(1, 
                                                         nrow(scenario[[j]]), 1), 
                                                     deriv = 1),
             second_derivative = rli_replicate_spline(seq(1, 
                                                          nrow(scenario[[j]]), 1), 
                                                      deriv = 2))
    
  }
  
  rli_scenario_derivatives[[i]] <- rli_replicate_derivatives
  
}

# Test plot RLI

rli_test <- rli_scenario_derivatives[[1]][[1]]
head(rli_test)

rli_5 <- ggplot(data = rli_test) +
  #geom_line(aes(x = annual_time_step, y = first_derivative)) +
  geom_line(aes(x = annual_time_step, y = second_derivative), col = "hot pink")

ggsave(file.path(analysis_outputs_folder, "rli_5yr_2nd_dervs.png"),
       rli_5,  device = "png")

# LIVING PLANET INDEX ----

# Create the spline functions and calculate derivatives (1 per replicate)

lpi_scenario_derivatives <- list()
lpi_replicate_derivatives <- list()

for (i in seq_along(lpi)) {
  
  # Get the outputs for a single scenario
  
  scenario <- lpi[[i]] 
  
  for (j in seq_along(scenario)) {
    
    # Get the outputs for a single replicate and create a spline function
    
    lpi_replicate_spline <- splinefun(x = scenario[[j]]$annual_time_step, 
                                      y = scenario[[j]]$indicator_score)
    
    # Use the spline function to calculate the first and second derivatives and
    # add them to our indicator score dataframes
    
    lpi_replicate_derivatives[[j]] <- scenario[[j]] %>% 
      mutate(first_derivative = lpi_replicate_spline(seq(1, 
                                                         nrow(scenario[[j]]), 1), 
                                                     deriv = 1),
             second_derivative = lpi_replicate_spline(seq(1, 
                                                          nrow(scenario[[j]]), 1), 
                                                      deriv = 2))
    
  }
  
  lpi_scenario_derivatives[[i]] <- lpi_replicate_derivatives
  
}

# Test plot LPI

lpi_test <- lpi_scenario_derivatives[[1]][[1]]
head(lpi_test)

lpi_5 <- ggplot(data = lpi_test) +
  #geom_line(aes(x = annual_time_step, y = first_derivative)) +
  geom_line(aes(x = annual_time_step, y = second_derivative), col = "hot pink")


ggsave(file.path(analysis_outputs_folder, "lpi_5yr_2nd_dervs.png"),
       lpi_5,  device = "png")
