
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

# Functions ----

#' Description

#' @param 
#' @param 
#' @param 
#' @return 

function_name <- function(param){

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

development_mode <- TRUE

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

# Load data ----

if (development_mode == TRUE) {
  
indicators_all <- readRDS("N:/Quantitative-Ecology/Indicators-Project/Serengeti/Outputs_from_indicator_code/Indicator_outputs/2021-07-12_all_indicators_output_data.rds")
indicators_all_list <- readRDS("N:/Quantitative-Ecology/Indicators-Project/Serengeti/Outputs_from_indicator_code/Indicator_outputs/2021-07-13_all_indicators_output_data_list.rds")

} else if (development_mode == FALSE) {
  
indicators_all <- readRDS("N:/Quantitative-Ecology/Indicators-Project/Serengeti/Outputs_from_indicator_code/Indicator_outputs/2021-07-12_all_indicators_output_data.rds")

}

# Split the list by indicators

rli <- indicators_all_list[["RLI"]]

lpi <- indicators_all_list[["LPI"]]

# RED LIST INDEX ----

# Create the spline functions and calculate derivatives (1 per replicate)

rli_scenario_derivatives <- list()
rli_replicate_derivatives <- list()

for (i in seq_along(rli)) {
  
  # Get the outputs for a single scenario
  
  scenario <- rli[[i]] 
  
  for (j in seq_along(scenario)) {
    
    # Get the outputs for a single replicate and create a spline function
    
    rli_replicate_spline <- splinefun(x = scenario[[j]]$time_step, 
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




