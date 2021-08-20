
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

#' TODO: Change LPI to geometric mean as per Mcrae 2017
#' TODO: Test sampling interval by trying different times of year
#' TODO: Check gen length equation
#' TODO: Check bootstrapping, looks weird for the LPI and maybe RLI (timesteps out of sync??)
#' TODO: Make LPI plots pretty
#' TODO: Make timesteps in plots annual (or whatever is relevant - when final sampling regime decided)

# Libraries ----

## Data wrangling

library(tidyverse)
library(rlist)
library(zoo)

# Functions ----

# Function to replace any non-extinct status following an extinction
# classification (resulting from juveniles growing into that massbin)
# with extinction
# https://stackoverflow.com/questions/45284226/how-to-replace-values-after-a-specific-event-in-a-data-frame

maintain_ex_status <- function(vec) { 
  vec_len = length(vec) 
  first_ex = match("EX", vec, nomatch = vec_len) 
  if(first_ex < vec_len) replace(vec, (first_ex+1):vec_len, "EX") 
  else vec 
}

# Similar function but replaces all non 0 values that occur after the first 
# 0 with 0

maintain_0_abundance <- function(vec) { 
  vec_len = length(vec) 
  first_0 = match(0, vec, nomatch = vec_len) 
  if(first_0 < vec_len) replace(vec, (first_0+1):vec_len, 0) 
  else vec 
}

# Scale values between 0 and 1

range01 <- function(x){(x-min(x, na.rm = TRUE))/(max(x, na.rm = TRUE)-min(x, na.rm = TRUE))}

# data <- rli_inputs
# numboots <- 5
# Function to calculate RLI

#' Calculate the RLI over time, from a dataframe of species and their 
#' IUCN Red List Categories 

#' @param data a data frame (should be named red list inputs) with columns:
#' group_id (representing functional group-massbin pseudo species), time_step
#' abundance, generation_length_yrs, "functional_group_index", "functional_group_name"
#' "mass_lower_g", "mass_upper_g", "massbin_g", timeframe, "diff" , "decline",
#' "rl_status",  "extinct", "rl_status"
#' @param numboots an integer representing the number of bootstraps to calculate
#' @param ci logical, do you want to include confidence intervals? default = FALSE
#' (note, selecting TRUE may increase processing time)
#' @return a dataframe of RLI index scores and confidence intervals over time

#data <- scenario_red_list_inputs_annual[[1]][[1]]

calculate_red_list_index <- function(data, numboots, ci = FALSE, replicate_num = NA){
  
  # Using equation from Butchart et al (2007) Improvements to the Red List Index
  
  require(tidyverse)
  
  # Remove data without RL status
  
  #data$redlist_assessment_year <- as.numeric(as.character(data$redlist_assessment_year))
  
  data <- data %>%
          filter(!is.na(rl_status)) %>%
          group_by(group_id) 
  
  head(data)
  
  # ecoregion <- as.factor(data$ecoregion_id[1])
  
  # Assign category weights
  
  weighted_data <- data %>%
    dplyr::mutate(rl_weight = ifelse(rl_status == "LC", 0,
                              ifelse(rl_status == "NT", 1,
                              ifelse(rl_status == "VU", 2,
                              ifelse(rl_status == "EN", 3,
                              ifelse(rl_status == "CR", 4,
                              ifelse(rl_status == "EX", 5, NA))))))) 
  head(weighted_data)
  dim(weighted_data)
 
  
  #weighted_data$RL_weight <- as.numeric(as.character(weighted_data$RL_weight))
  
  # Filter out rows with NE and DD
  weighted_data <- weighted_data %>%
                   filter(rl_status != "NE") %>%
                   filter(rl_status != "DD") %>%
                   filter(rl_status != "NA")
  
  dim(weighted_data)
  
  # Group data so the index is calculated for each functional group 
  # (would normally be taxa) for each year. If you run on a single group
  # it shouldn't matter, will just turn data into one big group
  
  grouped_data <- weighted_data %>% group_by(functional_group_name, annual_time_step)
  
  # Sum category weights for each group, in each timestep,
  # calculate number of species per group
  summed_weights <- summarise(grouped_data, 
                              total_weight = sum(rl_weight, na.rm = TRUE), # calc sum of all weights
                              total_count = n(),# calc number of species
                              .groups = "drop_last") %>%
                    mutate(total_count = max(total_count))  # Fix so it takes total number at beginning, otherwise n fluctuates between timesteps
  
  # Calculate RLI scores for each group, rounded to 3 decimal places
  
  index_scores <- summed_weights %>%
    mutate(RLI = 1 - (total_weight/(total_count * 5)), # actual RLI formula
           Criteria = "risk")

  if (ci == TRUE) {
  # Calculate confidence intervals via bootstrapping 
  # (see Rowland et al 2021 A guide to representing uncertainty)
  
  # Split by timestep - we want CI for each functional group, for each timestep 
  
  weighted_data_timestep_list <- split(weighted_data, weighted_data$annual_time_step)
  
  ## For each functional group (level 1)

  timestep_confidence_intervals <- list()
  
  for (i in seq_along(weighted_data_timestep_list)) {
    
    # Get single time-step then group by functional group
    
      grouped_timestep_data <- weighted_data_timestep_list[[i]] %>%
                             group_by(functional_group_name)
    
      time <- grouped_timestep_data$annual_time_step[1]
    
      boot <- list()
     # Calculate the bootstrap confidence intervals
      for (k in 1:numboots) {
        
        # Take k number of random samples from the weighted data
        replicate <- slice_sample(grouped_timestep_data, 
                                  prop = 1, replace = TRUE) %>%  # get random sample of rows and add to DF
                     mutate(replicate = k)
                   # label each replicate
        
        boot[[k]] <- replicate
        # Combine replicates into one dataframe
        
       # print(paste("Bootstrap", k, "of", numboots, "complete", sep =" "))
        
      }
    
      boot_reps <- do.call(rbind, boot)
      
      # Group by replicate
      #replicate_data <- group_by(boot_reps, replicate) # Group by replicate
        
      # Calculate the summary values needed to calc RLI for each replicate
      summed_weights_timestep_fg <- boot_reps %>%
                                    group_by(functional_group_name,replicate) %>% 
                                    summarise(total_weight = sum(rl_weight, na.rm = TRUE), # calc sum of all weights
                                              total_count = n(),# calc number of species
                                              .groups = "drop_last") %>%
                                    mutate(total_count = max(total_count))
        
      # Calculate the RLI score for each replicate
      rep_scores <- mutate(summed_weights_timestep_fg, 
                             RLI = 1 - (total_weight/(total_count * 5))) # actual RLI formula
        
      # Calculate the confidence intervals for each fg,
      ci_scores <- summarise(rep_scores, 
                               ci_lower = quantile(rep_scores$RLI, 
                                                   probs = 0.025),
                               ci_upper = quantile(rep_scores$RLI, 
                                                   probs = 0.975)) %>%
                   mutate(annual_time_step = time) 
      
      timestep_confidence_intervals[[i]] <- ci_scores
    
  }
  
  confidence_intervals <- do.call(rbind, timestep_confidence_intervals)
  
  red_list_scores <- index_scores %>%
                     merge(confidence_intervals, 
                           by = c("functional_group_name",
                                   "annual_time_step")) %>%
                     dplyr::select(functional_group_name, annual_time_step, ci_lower,
                                    RLI, ci_upper, everything()) %>% 
                     rename(indicator_score = RLI) %>% 
                     mutate(indicator = "RLI",
                            replicate = replicate_num)
                      
  
  return(red_list_scores)
  
  } else {
    
 red_list_scores <- index_scores  %>% 
                    rename(indicator_score = RLI) %>% 
                    mutate(indicator = "RLI",
                            replicate = replicate_num,
                            ci_lower = NA,
                            ci_upper = NA) %>% 
                    dplyr::select(functional_group_name, annual_time_step,
                                  ci_lower, indicator_score, ci_upper, total_weight,
                                  total_count, Criteria, indicator, replicate) 
                    
 
 return(red_list_scores)
 
    }
}

#' Return a line plot of the RLI scores over time faceted by functional group 

#' @param data a data frame (output from calculate_red_list_index) with columns:
#' functional_group, time_step, ci_lower, ci_upper, total.weight, total.count, RLI criteria
#' @param impact_start time_step impact began
#' @param impact_end time_step impact ended
#' @param ci logical, do you want to include confidence intervals? default = FALSE
#' (note, selecting TRUE may increase processing time)
#' @return a line plot of RLI over time for each functional group/taxa or whatever 
#' 
# data <- scenario_rli_outputs[[1]][[1]]

plot_red_list_index_by_group <- function(data, impact_start, impact_end, ci = FALSE) {
  
  require(ggplot2)
  require(viridis)
  
  if (ci == TRUE) {
  
  plot <- ggplot(data = data, aes(x = annual_time_step, y = indicator_score,
                           group = functional_group_name,
                           fill = functional_group_name)) +
    geom_line(aes(colour = functional_group_name)) +
    geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper),
                alpha = 0.4) +
    scale_fill_viridis_c() +
    scale_color_viridis_c() + 
    facet_wrap(~functional_group_name) +
    labs(x = "Time", 
         y = "Red List Index Score") +
    scale_y_continuous(limits = c(0,1)) +
    theme(panel.grid.major = element_blank(),
          axis.title = element_text(size = 18),
          axis.text = element_text(size = 18),
          panel.grid.minor = element_blank(),
          panel.background = element_rect(fill = "grey97"),
          axis.line = element_line(colour = "black"),
          legend.position = "none") +
    geom_vline(xintercept = impact_start, colour = "red") +
    geom_vline(xintercept = impact_end, colour = "blue")
  
  return(plot)
  
  } else {
    
    plot <- ggplot(data = data, aes(x = annual_time_step, y = indicator_score,
                                    group = functional_group_name,
                                    fill = functional_group_name)) +
      geom_line(aes(color = functional_group_name)) +
      scale_fill_viridis_d() +
      scale_color_viridis_d() + 
      facet_wrap(~functional_group_name) +
      labs(x = "Time", 
           y = "Red List Index Score") +
      scale_y_continuous(limits = c(0,1)) +
      theme(panel.grid.major = element_blank(),
            axis.title = element_text(size = 18),
            axis.text = element_text(size = 18),
            panel.grid.minor = element_blank(),
            panel.background = element_rect(fill = "grey97"),
            axis.line = element_line(colour = "black"),
            legend.position = "none") +
      geom_vline(xintercept = impact_start, colour = "red") +
      geom_vline(xintercept = impact_end, colour = "blue")
    
    return(plot)
    
  }

}

#' Return a line plot of the RLI scores over time for a single group 

#' @param data a data frame (output from calculate_red_list_index) with columns:
#' functional_group, time_step, ci_lower, ci_upper, total.weight, total.count, RLI criteria
#' @param impact_start time_step impact began
#' @param impact_end time_step impact ended
#' #' @param ci logical, do you want to include confidence intervals? default = FALSE
#' (note, selecting TRUE may increase processing time)
#' @return a line plot of RLI over time for each functional group/taxa or whatever 

plot_red_list_index <- function(data, impact_start, impact_end, ci = FALSE) {
  
  require(ggplot2)
  require(viridis)
  
  if (ci == TRUE) {
  
  plot <- ggplot(data = data, aes(x = annual_time_step, y = indicator_score)) +
    geom_line() +
    geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper),
                alpha = 0.4) +
    scale_fill_viridis_d() +
    scale_color_viridis_d() + 
    labs(x = "Time", 
         y = "Red List Index Score") + 
    scale_y_continuous(limits = c(0,1)) +
    theme(panel.grid.major = element_blank(),
          axis.title = element_text(size = 18),
          axis.text = element_text(size = 18),
          panel.grid.minor = element_blank(),
          panel.background = element_rect(fill = "grey97"),
          axis.line = element_line(colour = "black"),
          legend.position = "none") +
    geom_vline(xintercept = impact_start, colour = "red") +
    geom_vline(xintercept = impact_end, colour = "blue")
  
  return(plot)
  
  } else {
    
    plot <- ggplot(data = data, aes(x = annual_time_step, y = indicator_score)) +
      geom_line() +
      scale_fill_viridis_d() +
      scale_color_viridis_d() + 
      labs(x = "Time", 
           y = "Red List Index Score") + 
      scale_y_continuous(limits = c(0,1)) +
      theme(panel.grid.major = element_blank(),
            axis.title = element_text(size = 18),
            axis.text = element_text(size = 18),
            panel.grid.minor = element_blank(),
            panel.background = element_rect(fill = "grey97"),
            axis.line = element_line(colour = "black"),
            legend.position = "none") +
      geom_vline(xintercept = impact_start, colour = "red") +
      geom_vline(xintercept = impact_end, colour = "blue")
    
    return(plot)
  }
  
}

#' @param data a data frame with columns 'group id', 'time_step', 'abundance'
#' @param numboots an integer representing the number of bootstraps to calculate
#' @param ci logical, do you want to include confidence intervals? default = FALSE
#' (note, selecting TRUE may increase processing time)
#' @param replicate_num numeric/character identifier if calculating for more than one replicate
#' @return a dataframe of LPI index scores and confidence intervals over time

## Note: Have tested this code against the code Emily used on the LME ecopath
## data to see if they produce the same results, which they do.

#data <- scenario_lpi_inputs[[1]][[1]]

calculate_living_planet_index <- function(data, start_time_step = 1, ci = FALSE,
                                          numboots, replicate_num = NA){
  
  filtered_inputs <- data %>%
    # Remove timesteps if needed (if not make start_time_step = 1)
    filter(annual_time_step >= start_time_step) %>%
    # Group by virtual species (vs)
    group_by(group_id) %>% 
    # Remove vs that have no individuals at any timestep
    filter(!all(ave_abundance == 0)) %>% 
    # Save a copy of original abundance values
    # rename(abundance_original = abundance) %>%
    # Identify rows where abundance == 0, and change all subsequent years to 0 too
    # mutate(abundance = maintain_0_abundance(abundance_original)) 
    mutate(abundance = ave_abundance) %>% 
    dplyr::select(-ave_abundance)
  
  head(filtered_inputs)
  
  # Calculate LPI inputs and LPI
  ## Based on McRae, L., Loh. J., Bubb, P.J., Baillie, J.E.M., Kapos, V.,
  ## and Collen, B. 2008. The Living Planet Index - Guidance for
  ## National and Regional Use. UNEP-WCMC, Cambridge, UK.
  
  lpi_inputs <- filtered_inputs %>%
                # Calculate 1% of the mean population over time for each group
                group_by(group_id) %>%
                # Add 1% of the species mean abundance across all timesteps
                # so we can take the log in the next step
                mutate(abundance_adjusted = abundance +
                         (mean(abundance, na.rm = TRUE)*0.01)) %>%
                # Calculate the rate of change since the previous year (dt)
                # (equation 1 in Mcrae et al 2008)
                mutate(current_abundance = abundance_adjusted, #abundance at current timestep
                       previous_abundance = lag(abundance_adjusted, 1)) %>%  #abundance at previous timestep
                mutate(dt = log10(current_abundance/previous_abundance)) %>% # rate of change
                ungroup(.) %>%
                # Calculate mean dt across all groups, per timestep
                # (equation 4 in Mcrae et al 2008)
                group_by(annual_time_step) %>%
                summarise(mean_dt = mean(dt, na.rm = TRUE),
                          sd = sd(dt, na.rm = TRUE)) %>%
                # Add an empty LPI column to fill up in the next step
                mutate(LPI = NA)
  
  # Set the value of the LPI on the first time step to 1
  lpi_inputs$LPI[start_time_step] <- 1 

  # Annoyingly can't get a dplyr version of this to work, but whatever
  
       for (i in 1:(nrow(lpi_inputs) - 1)) {
         
         # T represents the current timestep we are calculating for, i represents
         # the previous timestep
         t <- i + 1
         
         # Multiple the LPI score from the previous timestep i by ten to the
         # power of dt in the current timestep t
         # Equation 5 in Mcrae et al 2008
         
         lpi_inputs$LPI[t] <- lpi_inputs$LPI[i] * (10 ^ lpi_inputs$mean_dt[t])
        
       
         }
  
  if(ci == TRUE) {
  # Now calculate the confidence intervals using bootstrapping as per 
  # Rowland et al 2021
  
  # Create n replicates of the original data (same number of rows but randomly
  # sampled, and with replacement)
  
  bootstrap_replicates <- list()
  
  for (j in 1:numboots) {
  
  # Make sure we get the same sample every time for reproducibility
    
  set.seed(j)
    
  # Produce a replicate dataset that is randomly sampled from the original 
  # dataset but with replacement
  
  rep <- filtered_inputs %>% 
         group_by(annual_time_step) %>% # so we take random samples stratified by timestep (otherwise end up with uneven sample sizes in each timestep)
         slice_sample(prop = 1, replace = TRUE) %>%  # get random sample of rows and add to DF
         mutate(replicate = j)
  
  # Calculate the mean rate of change across species
  
  rep_lpi_inputs <- rep %>%
    # Calculate 1% of the mean population over time for each group
    group_by(group_id) %>%
    # Add 1% of the species mean abundance across all timesteps 
    # so we can take the log in the next step
    mutate(abundance_adjusted = abundance + 
             (mean(abundance)*0.01)) %>%
    # Calculate the rate of change since the previous year (dt)
    # (equation 1 in Mcrae et al 2008) 
    mutate(current_abundance = abundance_adjusted, #abundance at current timestep
           previous_abundance = lag(abundance_adjusted, 1)) %>%  #abundance at previous timestep
    mutate(dt = log10(current_abundance/previous_abundance)) %>% # rate of change
    ungroup(.) %>% 
    # Calculate mean dt across all groups, per timestep
    # (equation 4 in Mcrae et al 2008)
    group_by(annual_time_step) %>% 
    summarise(mean_dt = mean(dt, na.rm = TRUE)) %>% 
    # Add an empty LPI column to fill up in the next step
    mutate(LPI = NA,
           replicate = j)
  
  # Set the value of the LPI on the first time step to 1
  rep_lpi_inputs$LPI[start_time_step] <- 1 
  
  # Calculate the LPI
  # Annoyingly can't get a dplyr version of this to work, but whatever
  
  for (i in 1:(nrow(rep_lpi_inputs) - 1)) {
    
    # T represents the current timestep we are calculating for, i represents
    # the previous timestep
    t <- i + 1
    
    # Multiple the LPI score from the previous timestep i by ten to the
    # power of dt in the current timestep t
    # Equation 5 in Mcrae et al 2008
    
    rep_lpi_inputs$LPI[t] <- rep_lpi_inputs$LPI[i] * (10 ^ rep_lpi_inputs$mean_dt[t])
    
    print(paste("bootstrap", j, "of", numboots, "complete", sep = " "))
    
  }
  
  bootstrap_replicates[[j]] <- rep_lpi_inputs
  
  }
  
  # Bind all the replicates into one DF
  
  bootstrap_replicates_df <- do.call(rbind, bootstrap_replicates)
  
  # Get the confidence intervals for each timestep
  
  ci_scores <- bootstrap_replicates_df %>% 
               group_by(annual_time_step) %>% # Get the LPI scores for all reps, for each timestep
               summarise(ci_lower = quantile(LPI,
                                             probs = 0.025), # Get the lower quantile of all LPI rep scores
                         ci_upper = quantile(LPI,
                                             probs = 0.975)) # Get the upper quantile
     
  # Merge the confidence intervals with the original LPI
  
  index_scores <- lpi_inputs %>% 
                  merge(ci_scores, by = "annual_time_step") %>% 
                  select(annual_time_step, LPI, ci_lower, ci_upper) %>% 
                  mutate(indicator = "LPI",
                         replicate = replicate_num) %>% 
                  rename(indicator_score = LPI)
  
  } else {
    
  index_scores <- lpi_inputs %>% 
    select(annual_time_step, LPI) %>%
    mutate(indicator = LPI,
           replicate = replicate_num,
           ci_lower = NA,
           ci_upper = NA) %>%
    rename(indicator_score = LPI)
  
  }
  
  return(index_scores)
 
}

# x <- calculate_living_planet_index(data, 1, FALSE, 1, "test")
# head(x)

#' @param data a data frame with columns 'time_step', 'mean_dt', 'LPI',
#' 'ci_lower' and 'ci_upper'
#' @param ci logical, do you want to plot confidence intervals? default = FALSE
#' (note, selecting TRUE may increase processing time)
#' @return a dataframe of LPI index scores and confidence intervals over time
#' 
plot_living_planet_index <- function(data, ci = FALSE) {
  
  if (ci == TRUE) {
  
    ggplot(data, aes(x = annual_time_step, y = indicator_score)) +
    geom_line()  +
    scale_y_continuous(limits = c(0,3)) +
    geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper),
                    alpha = 0.4)
  } else {
    
    ggplot(data, aes(x = annual_time_step, y = indicator_score)) +
      geom_line() +
      scale_y_continuous(limits = c(0, 3))
    
  }
}

# lpi <- calculate_living_planet_index(data, start_time_step = 1, ci = TRUE,
#                                      numboots = 1000)
# plot_living_planet_index(lpi, ci = TRUE)

# Set up paths ----

if (Sys.info()['nodename'] == "SIMONE-PC") {
  
  SourceToModels <- 'C:/Users/Simone/Dropbox/Deakin/Serengeti-analysis'
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

# disable scientific notation

options(scipen = 999)

# Get date to label outputs

today <- Sys.Date()

# Define mode (development == TRUE will only perform operations on a small subset
# of folders, not all outputs)

development_mode <- FALSE

# Specify universal input arguments

if (development_mode == FALSE) {
  

  burnin_months <- 1000*12 # in months
  n <- 12
  numboots <- 1000 # Rowland et al 2021 (uncertainty)
  start_time_step <- 1
  gen_timeframe <- 10 
  interval <- 12 
  # Don't adjust these
  max_timestep <- 300/(interval/12)
  impact_start <- max_timestep/3 * 1  #in years
  impact_end <- max_timestep/3 * 2  #in years
  
} else {
  
  # burnin_months <- 1 * 12 # in months
  # n <- 1 
  # numboots <- 100
  # start_time_step <- 1
  # gen_timeframe <- 10
  # interval <- 12
  # Don't adjust these
  # max_timestep <- 300/(interval/12)
  # impact_start <- max_timestep/3 * 1  #in years
  # impact_end <- max_timestep/3 * 2  #in years
  
  burnin_months <- 1000*12 # in months
  n <- 12
  numboots <- 2 # Rowland et al 2021 (uncertainty)
  start_time_step <- 1
  gen_timeframe <- 10 # in years
  interval <- 12 
  # Don't adjust these
  max_timestep <- 300/(interval/12)
  impact_start <- max_timestep/3 * 1  #in years
  impact_end <- max_timestep/3 * 2  #in years
  
}

indicators_project <- IndicatorsProject # File path for entire project directory

location <- 'Serengeti'

# Set up output folders

indicator_inputs_folder <- file.path(indicators_project, 
                           "/Serengeti/Outputs_from_indicator_code/Indicator_inputs")

if( !dir.exists( file.path(indicator_inputs_folder) ) ) {
  dir.create( file.path(indicator_inputs_folder), recursive = TRUE )
  
}

indicator_outputs_folder <- file.path(indicators_project, 
                                     "/Serengeti/Outputs_from_indicator_code/Indicator_outputs")

if( !dir.exists( file.path(indicator_outputs_folder) ) ) {
  dir.create( file.path(indicator_outputs_folder), recursive = TRUE )
  
}

indicator_plots_folder <- file.path(indicators_project, 
                                     "/Serengeti/Outputs_from_indicator_code/Indicator_plots")

if( !dir.exists( file.path(indicator_plots_folder) ) ) {
  dir.create( file.path(indicator_plots_folder), recursive = TRUE )
  
}

# Source functions ----

# note- would like to turn species_to_model repository into a package so you can eventually just run install_github( "conservationscience/species_to_models" )
# source( file.path( SourceToModels, "species_to_models", "madingley_get_groups.R" ) )
# source( file.path( SourceToModels, "species_to_models", "madingley_get_species_and_groups_key.R" ) )
# source( file.path( SourceToModels, "species_to_models", "madingley_process_trait_data.R" ) )
# source( file.path( SourceToModels, "species_to_models", "madingley_get_biomass_of_groups.R" ) )
# source( file.path( SourceToModels, "species_to_models", "madingley_get_abundance_of_groups.R" ) )
# source( file.path( SourceToModels, "species_to_models", "madingley_get_age_structure_data.R" ) )
# source( file.path( SourceToModels, "species_to_models", "madingley_get_autotroph_biomass.R" ) )
# 
# 
# source( file.path( SourceToModels, "model_outputs_to_indicator_inputs", 
#                    "1_process_outputs", "process_species_list.R") )
# source( file.path( SourceToModels, "model_outputs_to_indicator_inputs", 
#                    "1_process_outputs", "process_buildmodel_folder.R") )
# source( file.path( SourceToModels, "model_outputs_to_indicator_inputs", 
#                    "1_process_outputs", "process_output_functions.R") )
# 
# source( file.path( SourceToModels, "model_outputs_to_indicator_inputs", 
#                    "1_process_outputs", "plot_functional_groups.R") )

# source( file.path( SourceToModels, "model_outputs_to_indicator_inputs", 
#                    "2_prepare_inputs", "prepare_proportion_total_biomass_inputs.R") )
# 
# source( file.path( SourceToModels, "biodiversity_indicators", 
#                    "calculate_proportion_biomass.R") )

# Get file paths of processed model output data ---- 

processed_outputs_path <- file.path(IndicatorsProject, location,
                                    "Outputs_from_adaptor_code/map_of_life")

# Get lists of the directories at various levels needed

# * Scenario & simulation directories ----

processed_scenario_paths <- list.dirs(processed_outputs_path, recursive = FALSE)

if (development_mode == TRUE) {
  
  all_processed_scenario_paths <- processed_scenario_paths
  #processed_scenario_paths <- all_processed_scenario_paths[str_detect(all_processed_scenario_paths, "999_Test_runs")]
  processed_scenario_paths <- "N:\\Quantitative-Ecology\\Indicators-Project\\Serengeti\\Outputs_from_adaptor_code\\map_of_life\\666_Long_test_runs"
}


if (development_mode == FALSE) {
  
  processed_scenario_paths <- processed_scenario_paths[!str_detect(processed_scenario_paths, "999_Test_runs")] # Remove empty/unneeded scenarios
  
}

# Get scenario names so we can label outputs

scenarios <- lapply(processed_scenario_paths, basename)

## Get a list of the simulation directories within each scenario

processed_simulation_paths <- list()

for (i in seq_along(processed_scenario_paths)) {
  
  processed_simulation_paths[[i]] <- list.dirs(processed_scenario_paths[[i]], 
                                               recursive = FALSE)
  
}

processed_simulation_paths

# * Replicate paths ----

replicate_paths <- list()

scenario_replicate_paths <- list()

for (i in seq_along(processed_scenario_paths)) {
  
  scenario_path <- processed_scenario_paths[i] # Get the file path for one scenario
  
  scenario_simulation_paths <- list.dirs(scenario_path, recursive = FALSE) 
  
  # Get the file paths for each simulation within the scenario
  
  replicate_paths <- list()
  
  for (k in seq_along(scenario_simulation_paths)) { # For each simulation
    
    single_simulation_path <- scenario_simulation_paths[k]
    
    single_simulation_replicates <- list.files(single_simulation_path) #, pattern = "biomass.rds") # list the biomass files for that simulation and store
    
    replicate_paths[[k]] <- file.path(single_simulation_path, single_simulation_replicates)
    
  }
  
  scenario_replicate_paths[[i]] <- replicate_paths # Get the replicates for each simulation and store in correct scenario
  
}

scenario_replicate_paths

# Get groups ----
# Details for the virtual species or 'groups'
# Groups file is the same for all simulations, so can just pull it from whichever directory

groups <- readRDS(file.path(processed_simulation_paths[[1]][1], "groups.rds"))

# Get abundance ----

# read in the abundance and generation length data. Output will be a nested list 
# of scenarios (x3), each scenario containing a sub-list of abundance data for 
# each replicate, and a list of scenarios (x3), each containing a sub list of 
# generation length
# data for each replicate (x25)

scenario_generations_raw <- list()
scenario_abundance_raw <- list()
scenario_autotroph_raw <- list()

for (i in seq_along(scenario_replicate_paths)) {
  
scenario <- scenario_replicate_paths[[i]]

scenario <- flatten(scenario)

# Get generation length of all replicates

generation_files <- str_subset(scenario, "GenerationLengths")
generation_files <- generation_files[!str_detect(generation_files, ".png")]
generation_files <- generation_files[!str_detect(generation_files, ".csv")]

abundance_files <- str_subset(scenario, "abundance")
abundance_files <- abundance_files[!str_detect(abundance_files, ".png")]
abundance_files <- abundance_files[!str_detect(abundance_files, ".csv")]

autotroph_files <- str_subset(scenario, "autotroph")
autotroph_files <- autotroph_files[str_detect(autotroph_files, ".rds")]


# Check we've got the correct number of files for each (should match)

if (length(abundance_files) != length(generation_files)) {
  
  stop(paste("Number of abundance files and generation files in",
  scenarios[[i]], "do not match", sep = " "))

  }

# Read in the data
# Note, generation lengths are different sizes bc each replicate has slightly different number of groups present
scenario_generations_raw[[i]] <- lapply(generation_files, readRDS) 
scenario_abundance_raw[[i]] <- lapply(abundance_files, readRDS)
scenario_autotroph_raw[[i]] <- lapply(autotroph_files, readRDS)

}
scenario_autotroph_raw[[1]][[1]][1:5,1:5]
# scenario_generations_raw_all <- scenario_generations_raw
# scenario_abundance_raw_all <- scenario_abundance_raw
# 
# scenario_generations_raw <- list(scenario_generations_raw_all[[3]])
# scenario_abundance_raw<- list(scenario_abundance_raw_all[[3]])
# 
# scenarios <- list(scenarios[[3]])

# Remove burn-in ----

# Remove the burn in timesteps for abundance files (output should be same structure
# as input, a nested list of x3 scenarios, each with x25 replicate dataframes,
# but each dataframe will have fewer columns, no NA data and be absolute not log
# abundance densities)

scenario_abundance_formatted <- list()

for (i in seq_along(scenario_abundance_raw)) {
  
  # Get all replicates for one scenario
  abundance_reps <- scenario_abundance_raw[[i]]

  # Make a list to capture the outputs
  
  replicate_abundance_formatted <- list()

  # For each individual replicate
  
  for (j in seq_along(abundance_reps)) {
    
    rep <- abundance_reps[[j]]
      
    # Make correct numeric column names now
    
    col_names <- seq(1,ncol(rep),1)
    
    colnames(rep) <- col_names
    
    # Remove the burnin period 
    
    abundance_temp <- rep[,burnin_months:ncol(rep)]
    
    # Replace all -9999 (NA values) values with 0 because in this case we 
    # know that an NA is a true 0 (no abundance at that time step)
    
    abundance_temp <- na_if(abundance_temp, 0)
    
    # Exponentiate - model output abundance (from MassBinsOutput) is the 
    # log density of individuals as specified in the Madingley User Guide.
    
    abundance_temp <- 10 ^ abundance_temp
    
    replicate_abundance_formatted[[j]] <- abundance_temp
    
   }
  
  scenario_abundance_formatted[[i]] <- replicate_abundance_formatted
  
}

scenario_abundance_formatted[[1]][[1]][1:20,1:20]


scenario_auto_long <- list()

for (i in seq_along(scenario_autotroph_raw)) {
  
  # Get all replicates for one scenario
  auto_reps <- scenario_autotroph_raw[[i]]
  
  # Make a list to capture the outputs
  
  replicate_auto_formatted <- list()
  
  # For each individual replicate
  
  for (j in seq_along(auto_reps)) {
    
    rep <- auto_reps[[j]]
    
    # Make correct numeric column names now
    
    col_names <- seq(1,ncol(rep),1)
    
    colnames(rep) <- col_names
    
    # Remove the burnin period 
    
    auto_temp <- rep[,burnin_months:ncol(rep)]
    
    # Replace all -9999 (NA values) values with 0 because in this case we 
    # know that an NA is a true 0 (no abundance at that time step)
    
    auto_temp <- na_if(auto_temp, 0)
    
    # Convert to long format
    
    new_names <- colnames(auto_temp)
   
    # Convert from wide to long
    
    long <- auto_temp %>% 
      rownames_to_column(.) %>%
      pivot_longer(all_of(new_names)) %>%
      rename(monthly_time_step = name,
             abundance = value,
             group_id = rowname) %>%
      # mutate(abundance = 10 ^ abundance) %>% 
      mutate(monthly_time_step = as.numeric(monthly_time_step),
             abundance = as.numeric(abundance),
             replicate = j - 1,
             scenario = scenarios[[i]])  
      
    
     replicate_auto_formatted[[j]] <- long
    
  }
  
  scenario_auto_long[[i]] <- replicate_auto_formatted
  
}

head(scenario_auto_long[[1]][[1]])

# rm(abundance_reps, abundance_temp, replicate_abundance_formatted, rep)

# Check structure is still correct

length(scenario_abundance_formatted) == length(scenario_abundance_raw)

# Remove previous version to make space

rm(scenario_abundance_raw)

# Convert to long format ----

## NOTE: Now the burnin period has been removed, monthly time step now begins
## at 12000

## Pivot longer (group_id, timestep, abundance)

scenario_abundance_long <- list()

for (i in seq_along(scenario_abundance_formatted)) {
  
  # Get data for one scenario
  
  scenario_abundance <- scenario_abundance_formatted[[i]]
  
  # Make a list to catch replicate outputs
  
  replicate_abundance_long <- list()
  
  # For each individual replicate
  
  for (j in seq_along(scenario_abundance)) {
  # Rename columns as timesteps so we can use them as a numeric variable
  abundance_single <- as.data.frame(scenario_abundance[[j]])
  new_names <- colnames(abundance_single)
  # new_names <- str_remove(new_names, "V")
  # colnames(abundance_single) <- new_names
  
  # Convert from wide to long
  
  replicate_abundance_long[[j]] <- abundance_single %>% 
    rownames_to_column(.) %>%
    pivot_longer(all_of(new_names)) %>%
    rename(monthly_time_step = name,
           abundance = value,
           group_id = rowname) %>%
    mutate(monthly_time_step = as.numeric(monthly_time_step),
           abundance = as.numeric(abundance))
  
  }
 
   scenario_abundance_long[[i]] <- replicate_abundance_long

}

head(scenario_abundance_long[[1]][[1]])
test_group <- scenario_abundance_long[[1]][[1]] %>% filter(group_id == "10.42")
head(test_group)

# TEST CODE removing problematic carnivore reps ----

carnivore_scenario <- scenario_abundance_long[[3]]

carnivore_scenario_2 <- carnivore_scenario[1:30]

carnivore_scenario_3 <- carnivore_scenario_2[- c(21:26)]

scenario_abundance_long_og <- scenario_abundance_long

scenario_abundance_long[[3]] <- carnivore_scenario_3

# # Take a small representative sample
# 
# landuse_herbs <- scenario_abundance_long[[2]][[1]] %>% 
#                  filter(group_id == "10.38") %>% 
#                  mutate(abundance = ifelse(abundance == 0, NA, abundance))
# 
# write.csv(landuse_herbs, file.path(indicator_outputs_folder, "210817_example_time_series.csv"))
# 
# ggplot(data = landuse_herbs, aes(x = monthly_time_step, y = abundance)) +
#   geom_line()
# 
# ## ANNUAL SAMPLING ### ----
# 
# # Sample and get generation length ----
# 
# scenario_ab_gl_formatted_not_clean <- list()
# 
# for (i in seq_along(scenario_abundance_long)) {
#   
#   replicate_abundance <- scenario_abundance_long[[i]]
#   replicate_generations <- scenario_generations_raw[[i]]
#   
#   # Make a list to catch the outputs
#   
#   replicate_ab_gl_formatted <- list()
#  
#   # For each individual replicate
#   
#   for (j in seq_along(replicate_abundance)) {
#     
#     # Reduce size of the replicate generations dataframe or the merge won't work
#     gen_length <- replicate_generations[[j]] %>% 
#                   dplyr::select(group_id, generation_length_yrs, 
#                                 functional_group_name) %>% 
#                   distinct(.)
#     
#     # Add the generation length info to the abundance dataframe
#     replicate_ab_gl_formatted[[j]] <- replicate_abundance[[j]] %>%
#         merge(gen_length, by = "group_id") %>%
#         arrange(monthly_time_step, group_id) %>%
#     # Get the timeframe over which to assess decline (3 * gen length or 10 yrs,
#       # whichever is longer)
#     # Important - following lines assume an annual timeframe, will need to adjust if change interval
#     mutate(generation_by_three = generation_length_yrs * 3) %>% # Time over which to measure decline, 3 x gen length OR:
#     mutate(timeframe = ifelse(generation_by_three > gen_timeframe, # 10 years 
#                        round(generation_by_three), gen_timeframe)) %>%
#     dplyr::select(-generation_by_three) %>%
#     distinct(.) %>%
#     group_by(group_id) %>% 
#       # select rows that are multiples of the specified interval 
#       # (eg if interval is 12, it samples one month from every 12 (yearly))
#     slice(which(row_number() %% interval == 0)) %>% 
#     mutate(annual_time_step = seq(1,max_timestep,1)) # %>% 
#     
# 
#     print(paste("Replicate", j - 1, 
#                 "formatting complete", 
#                 sep = " "))
#     
#   }
#     
#   print(scenario[[i]])
#   print(length(replicate_ab_gl_formatted))
#   
#   scenario_ab_gl_formatted_not_clean[[i]] <- replicate_ab_gl_formatted
#     
# }
#   
# # Remove false extinctions ----
# 
# scenario_false_extinctions_removed <- list()
# 
# for (i in seq_along(scenario_ab_gl_formatted_not_clean)) {
#   
#   replicate_ab_gl <- scenario_ab_gl_formatted_not_clean[[i]]
#  
#   # Make a list to catch the outputs
#   
#   replicate_false_ex_removed <- list()
#   
#   # For each individual replicate
#   
#   for (j in seq_along(replicate_ab_gl)) {
#     
#     # Find the last time step where non-0 abundance occurred for each group
#     
#     temp2 <- replicate_ab_gl[[j]] %>% 
#       group_by(group_id) %>% 
#       filter(abundance > 0) %>% 
#       dplyr::select(group_id, annual_time_step, abundance) %>% 
#       filter(annual_time_step == max(annual_time_step)) %>% 
#       dplyr::select(group_id, annual_time_step) %>% 
#       rename(last_abundance = annual_time_step)
#     
#     # Add the year of last positive abundance number as a column to the data    
#     temp3 <- replicate_ab_gl[[j]] %>% 
#       merge(temp2, by = c("group_id"), all = TRUE)
#     
#     # Use the last positive abundance year and current abundance value to determine
#     # if a zero abundance is a true extinction or just a missing value (false extinction)
#     temp4 <- temp3 %>%
#       group_by(group_id) %>%
#       mutate(true_extinction = ifelse(abundance == 0 &
#                                         annual_time_step < last_abundance,
#                                       "false extinction",
#                                       ifelse(abundance > 0 &
#                                                annual_time_step < last_abundance,
#                                              "not extinct",
#                                              ifelse(abundance == 0 &
#                                                       annual_time_step >= last_abundance,
#                                                     "true extinction", "not extinct")))) %>%
#       filter(true_extinction != "false extinction") %>%
#       group_by(group_id) %>%
#       arrange(annual_time_step)
#     
#     # Add massbin index
#     data <- temp4 %>% 
#             merge(groups[c("group_id", "bodymass_index", "mass_lower")],
#                   by = "group_id") %>% 
#             arrange(functional_group_name,
#                     bodymass_index, annual_time_step)
#             
#     # 
#     # data <- temp3 %>% 
#     #   group_by(group_id) %>% 
#     #   # Identify false extinctions (where abundance = 0 but it's just missing 
#     #   # data/cohorts moving massbins)
#     #   mutate(true_extinction = ifelse(abundance == 0 & 
#     #                            annual_time_step < last_abundance,
#     #                            "false extinction",
#     #                            ifelse(abundance > 0 & 
#     #                            annual_time_step < last_abundance,
#     #                            "not extinct",
#     #                            ifelse(abundance == 0 & 
#     #                            annual_time_step >= last_abundance,
#     #                            "true extinction", "not extinct")))) %>% 
#     #   #filter(true_extinction != "false extinction") %>% 
#     #   # Convert the false zeroes to NA
#     #   mutate(abundance = ifelse(true_extinction == "false extinction",
#     #                             NA, abundance)) %>%
#     #   group_by(group_id) %>% 
#     #   arrange(annual_time_step)
#     
#     
#     # Check if there are any carnivorous endotherms
#     
#     check <- data %>% 
#       group_by(functional_group_name) %>% 
#       summarise(present = sum(abundance)) %>% 
#       filter(functional_group_name == "carnivore endotherm") %>% 
#       dplyr::select(present) %>% 
#       pull(.)
#     
#     
#     print(paste("Replicate", j - 1, 
#                 "formatting complete", 
#                 sep = " "))
#     
#     # Replace data with 0 if no carnivores
#     
#     if(length(check) == 0) {
#       
#       data <- NULL
#       
#       print(paste("Replicate", j - 1, 
#                   "removed because no carnivorous endotherms are present", 
#                   sep = " "))
#       
#     }
#     
#     replicate_false_ex_removed[[j]] <- data
#     
#   }
#   
#   print(scenario[[i]])
#   print(length(replicate_ab_gl))
#   
#   scenario_false_extinctions_removed[[i]] <- replicate_false_ex_removed
#   
# }
# 
# # Remove replicates with no carnivorous endotherms ----
# 
# scenario_ab_gl_formatted <- list()
# 
# for (i in seq_along(scenario_false_extinctions_removed)) {
#   
#   replicate_not_clean <- scenario_false_extinctions_removed[[i]]
#   
#   scenario_ab_gl_formatted[[i]] <- list.clean(replicate_not_clean)
#   
# }
# 
# 
# # TEST CODE ----
# # test_group <- "14.17.35"
# # 
# # # Completely formatted long data
# # length(scenario_false_extinctions_removed[[1]])
# # 
# # long_formatted <- scenario_false_extinctions_removed[[1]][[2]]%>% 
# #   filter(group_id == test_group)
# # 
# # # Not formatted long data
# # length(scenario_abundance_long[[1]])
# # 
# # long_unformatted <- scenario_abundance_long[[1]][[2]] %>% 
# #                     filter(group_id == test_group)
# # 
# # # Not formatted wide data
# # x <- scenario_abundance_formatted[[1]][[2]]
# # wide_formatted <- as.matrix(scenario_abundance_formatted[[1]][[2]][test_group,])
# # 
# # compare <- long_unformatted %>% 
# #   rename(og_abundance = abundance) %>% 
# #   merge(long_formatted[c("monthly_time_step",
# #                          "annual_time_step", 
# #                          "abundance")],
# #         by = "monthly_time_step",
# #         all = TRUE) %>%
# #   rename(new_abundance = abundance) %>% 
# #   mutate(mean_abundance = rollmean(og_abundance, 10, fill = NA)) %>% 
# #   mutate(og_abundance_NA = ifelse(og_abundance == 0, NA, og_abundance)) %>% 
# #   mutate(mean_abundance_NA = rollapplyr(og_abundance_NA, 120, mean, by = 120, 
# #                                         partial = TRUE, na.rm = TRUE, 
# #                                         align = "left", fill = "extend"))
# # 
# # compare <- cbind(wide_formatted, compare)
# # 
# # compare_long <- compare %>% 
# #                 pivot_longer(c(wide_formatted, og_abundance,
# #                              new_abundance, mean_abundance,
# #                              og_abundance_NA, mean_abundance_NA)) 
# # 
# # # test_nums_na <- c(1, 2, 3 , NA, 1, 2, 3, NA,1, 2, 3, NA,1, 2, 3, NA,1, 2, 3, NA)
# # # groups <- c(1,1,1,1,2,2,2,2,3,3,3,3,4,4,4,4,5,5,5,5)
# # # 
# # # test_nums <- as.data.frame(cbind(test_nums_na, groups)) 
# # # 
# # # test_nums_sum <- test_nums %>% 
# # #              group_by(groups) %>% 
# # #              summarise(mean = mean(test_nums_na, na.rm = TRUE))
# # # 
# # # mean_extend <- rollmean(test_nums_na, 4, partial = TRUE, fill = "extend")
# # # 
# # # mean_fill <- rollmean(test_nums_na, 4, fill = NA)
# # # 
# # # mean_narm <- rollapplyr(test_nums_na, 4, mean, partial = TRUE, na.rm = TRUE, align = "center")
# # # 
# # # test <- cbind(test_nums, mean_narm)
# # 
# # 
# # ggplot(data = compare) +
# #   geom_line(aes(x = monthly_time_step, y = wide_formatted), 
# #             col = "yellow") +
# #   geom_line(aes(x = monthly_time_step, y = og_abundance), 
# #             alpha = 0.3, col = "blue") +
# #   geom_point(aes(x = monthly_time_step, y = new_abundance), 
# #             col = "red") 
# # 
# # ggplot(data = compare) +
# #   geom_line(aes(x = monthly_time_step, y = mean_abundance), 
# #             col = "deep pink") 
# # 
# # ggplot(data = compare_long) +
# #   geom_line(aes(x = monthly_time_step, y = value), 
# #             col = "red") +
# #   facet_wrap(~ name)
# # 
# # 
# # compare_annual <- compare %>% 
# #                   slice(which(row_number() %% 60 == 0))
# #   
# # compare_long_annual <- compare_annual %>% 
# #   pivot_longer(c(wide_formatted, og_abundance,
# #                  new_abundance, mean_abundance,
# #                  og_abundance_NA, mean_abundance_NA)) 
# # 
# # 
# # ggplot(data = compare_annual) +
# #   geom_line(aes(x = monthly_time_step, y = wide_formatted), 
# #             col = "yellow") +
# #   geom_line(aes(x = monthly_time_step, y = og_abundance), 
# #             alpha = 0.3, col = "blue") +
# #   geom_point(aes(x = monthly_time_step, y = new_abundance), 
# #              col = "red") 
# # 
# # ggplot(data = compare_long_annual) +
# #   geom_line(aes(x = monthly_time_step, y = value), 
# #             col = "purple") +
# #   facet_wrap(~ name)
# 
# # Identify and deal with weird mass bins that blink in and out
# 
# # scenario_abundance_clean <- list()
# # 
# # for (i in seq_along(scenario_ab_gl_formatted)) {
# # 
# #   replicate_allgroups <- scenario_ab_gl_formatted[[i]]
# #   replicate_gens <- scenario_generations_raw[[i]]
# # 
# #   reps_out <- list()
# # 
# #   for (j in seq_along(replicate_allgroups)) {
# # 
# #     data <- replicate_allgroups[[j]]
# # 
# #     gen <- replicate_gens[[j]] %>%
# #            dplyr::select(group_id, mass_lower_g) %>%
# #            distinct(.)
# # 
# #    # Determine which groups were there at beginning (post burnin)
# #     temp <- data %>%
# #       group_by(group_id) %>%
# #       filter(monthly_time_step == min(monthly_time_step)) %>%
# #       mutate(first_appearance = annual_time_step,
# #              beginning = ifelse(first_appearance == 1,
# #                                 TRUE, FALSE)) %>%
# #       dplyr::select(group_id, first_appearance, beginning)
# # 
# #     reps_out[[j]] <- data %>%
# #         merge(temp, by = "group_id") %>%
# #         merge(gen, by = "group_id") %>%
# #         arrange(mass_lower_g) %>%
# #         tidylog::filter(beginning == TRUE) # Note 14% is highest percentage of data removed by this line
# # 
# #     rm(temp)
# #   }
# # 
# #   scenario_abundance_clean[[i]] <- reps_out
# # 
# # }
# 
# # * Smooth abundance ----
# 
# scenario_abundance_clean <- scenario_ab_gl_formatted
# 
# 
# ave_window <- 10
# 
# scenario_smoothed_abundance <- list()
# 
# for (i in seq_along(scenario_abundance_clean)) {
# 
#   # Get replicate data for a single scenario
# 
#   replicate_ab_gl <- scenario_abundance_clean[[i]]
# 
#   replicate_smoothed_abundance <- list()
# 
#   # For each individual replicate
# 
#   for (j in seq_along(replicate_ab_gl)) {
# 
#     group_ab_gl <- replicate_ab_gl[[j]]
# 
#     group_list <- split(group_ab_gl, group_ab_gl$group_id)
# 
#     group_smoothed_abundance <- list()
# 
#     for (k in seq_along(group_list)) {
# 
#       group_df <- group_list[[k]]
# 
#       if (is.na(sum(group_df$abundance))) {
# 
#       group_smoothed_abundance[[k]] <- NULL
# 
#       } else {
# 
#       group_smoothed_abundance[[k]] <- group_df %>%
#                                        arrange(annual_time_step) %>%
#                                        mutate(ave_abundance = rollapply(abundance,
#                                                               ave_window,
#                                                               mean,
#                                                               na.rm = TRUE,
#                                                               partial = TRUE),
#                                               ave_abundance = ifelse(ave_abundance < 1,
#                                                                      0, ave_abundance))
# 
#       print(k)
# 
#       }
#     }
# 
#     all_groups_smooth <- do.call(rbind,group_smoothed_abundance)
# 
#     replicate_smoothed_abundance[[j]] <- all_groups_smooth
# 
#     print(j)
#   }
# 
#   scenario_smoothed_abundance[[i]] <- replicate_smoothed_abundance
# 
#   print(i)
# 
# }
# 
# check <- scenario_smoothed_abundance[[1]][[1]]
# head(check)
# # RED LIST INDEX ----
# 
# # * Create folders ----
# 
# 
# rli_inputs_folder <- file.path(indicator_inputs_folder, "RLI_inputs", today)
# 
# if( !dir.exists( file.path(rli_inputs_folder) ) ) {
#   dir.create( file.path(rli_inputs_folder), recursive = TRUE )
#   
# }
# 
# rli_outputs_folder <- file.path(indicator_outputs_folder, "RLI_outputs", today)
# 
# if( !dir.exists( file.path(rli_outputs_folder) ) ) {
#   dir.create( file.path(rli_outputs_folder), recursive = TRUE )
#   
# }
# 
# rli_plots_folder <- file.path(indicator_plots_folder, "RLI_plots", today)
# 
# if( !dir.exists( file.path(rli_plots_folder) ) ) {
#   dir.create( file.path(rli_plots_folder), recursive = TRUE )
#   
# }
# 
# ## Referring to the thresholds quote under Criterion A, Reason 1 (declines
# ## are the result of reversible pressures) according to:
# ## https://portals.iucn.org/library/sites/library/files/documents/RL-2001-001-2nd.pdf
# 
# 
# # * Assign Red List Categories ----
# 
# scenario_red_list_data <- list()
# 
# #for (i in seq_along(scenario_ab_gl_formatted)) {
# for (i in seq_along(scenario_smoothed_abundance)) {
#   
#   # Get replicate data for a single scenario
#   
#   replicate_ab_gl <- scenario_smoothed_abundance[[i]]
#   
#   print(paste("Processing scenario", scenarios[[i]], sep = " "))
#   
#   replicate_red_list_data <- list()
#   
#   # For each individual replicate
#   
#   for (j in seq_along(replicate_ab_gl)) {
#     
#   print(paste("Processing replicate", j, sep = " "))
#   
#   # Split by functional group, because we calculate RLI for different
#   # functional groups then aggregate later (as per Butchart etal 2010),
#   # except we are using functional groups as proxies for taxa (eg mammals, birds, 
#   # reptiles) used in real world RLI calcs
#   
#   status_inputs <- split(replicate_ab_gl[[j]], 
#                          replicate_ab_gl[[j]]$group_id)
#   
#   # Make a list to hold output for each individual massbin-func-group (ie virtual spp)
#   
#   group_red_list_data <- list()
#   
#   for (k in seq_along(status_inputs)) {
#     
#     print(paste("Processing group", names(status_inputs)[[k]], sep = " "))
#     
#     group_red_list_data[[k]] <- status_inputs[[k]] %>%
#       group_by(group_id) %>%
#       arrange(monthly_time_step) %>%
#       # calculate the difference in abundance over 10 yrs or 3 generation lengths
#       # (specified by 'timeframe' column). Its okay to take the first value of 
#       # timeframe bc the dataframe is grouped by group_id, and timeframe only changes
#       # between and not within group_ids
#       # mutate(diff = (abundance - dplyr::lag(abundance, timeframe[1]))) %>%
#       mutate(diff = (ave_abundance - dplyr::lag(ave_abundance, timeframe[1]))) %>%
#       # calculate the rate of change
#       # mutate(decline = diff/dplyr::lag(abundance, timeframe[1])) %>% 
#       mutate(decline = diff/dplyr::lag(ave_abundance, timeframe[1])) %>% 
#       # assign red list risk status based on decline 
#       mutate(rl_status = ifelse(decline > -0.40, "LC",
#                          ifelse(decline <= -0.40 & decline > -0.50, "NT", # Where did this and LC thresholds come from?
#                          ifelse(decline <= -0.50 & decline > -0.70, "VU",
#                          ifelse(decline <= -0.70 & decline > -0.90, "EN",
#                          ifelse(decline <= -0.90 & decline > -1, "CR",
#                          ifelse(decline <= -1, "EX", "NA"))))))) %>%
#       arrange(group_id, monthly_time_step) %>%
#       # Replace all non-ex status with ex after first occurrence 
#       # mutate(extinct = match("EX", rl_status)) %>%
#       mutate(extinct = ifelse(rl_status == "EX", 1, 0)) %>% 
#       # mutate(rl_status = with(., ave(rl_status, 
#       #                                         FUN=maintain_ex_status)))
#       #mutate(rl_status = rl_status) %>% 
#       group_by(group_id)
# 
#   }
#   
#   print(paste("replicate", j, "from", scenarios[[i]], "complete", sep = " "))
#   
#   replicate_red_list_df <- do.call(rbind, group_red_list_data)
#    
#   replicate_red_list_data[[j]] <- replicate_red_list_df
#    
#    # Save the inputs
#    
#    saveRDS(replicate_red_list_df,
#            file.path(rli_inputs_folder,
#                      paste(today, scenarios[[i]], "replicate", j,
#                            "RLI_input_data.rds", sep = "_")))
# 
#    write.csv(replicate_red_list_df,
#            file.path(rli_inputs_folder,
#                      paste(today, scenarios[[i]], "replicate", j,
#                            "RLI_input_data.csv", sep = "_")))
#   
#   
#   }
# 
#   scenario_red_list_data[[i]] <- replicate_red_list_data
#   
# }
# 
# # Check we have correct structure still
# length(scenario_red_list_data) == length(scenario_ab_gl_formatted)
# length(scenario_red_list_data[[1]]) == length(scenario_ab_gl_formatted[[1]])
# 
# # Have a quick look at the outputs
# 
# rli_inputs <- scenario_red_list_data[[1]][[1]]
# tail(rli_inputs)
# 
# write.csv(rli_inputs, file.path(indicator_outputs_folder, "rli_input_example_annual.csv"))
# 
# # Plot some results to check they're not completely whack
# 
# ## Get one group to check how their status changes over time relative to how
# ## their abundance changes
# 
# # group_id_select <- "13.16.17" # Shows example of 'resurrected' virtual spp
# # # group_id_select <- "10.40"
# # 
# # data <- rli_inputs %>% dplyr::filter(group_id == group_id_select)
# # 
# # ggplot(data, aes(x = time_step, y = abundance)) +
# #   geom_line() +
# #   geom_text(aes(label= rl_status,
# #                 col = rl_status),hjust=0, vjust=0)
# 
# 
# 
# 
# # * Calculate RLI ----
# 
# # RLI by individual functional groups
# 
# scenario_fg_rli_outputs <- list()
# 
# for (i in seq_along(scenario_red_list_data)) {
#   
#   replicate_red_list_inputs <- scenario_red_list_data[[i]]
#   
#   replicate_fg_rli_outputs <- list()
#   
#   for (j in seq_along(replicate_red_list_inputs)) {
#     
#   replicate_rli <- calculate_red_list_index(
#     replicate_red_list_inputs[[j]], numboots, ci = FALSE) %>%
#     mutate(replicate = j)
#   
#   replicate_fg_rli_outputs[[j]] <- replicate_rli 
#   
#   # saveRDS(replicate_fg_rli_outputs[[j]],
#   #         file.path(rli_outputs_folder,
#   #                   paste(today, scenarios, "replicate", j,
#   #                         "RLI_func_group_output_data.rds",
#   #                         sep = "_")))
# 
#   write.csv(replicate_fg_rli_outputs[[j]],
#             file.path(rli_outputs_folder,
#                       paste(today, scenarios[[i]], "RLI_func_group_output_data.rds",
#                             sep = "_")))
#   
#   print(paste("RLI for replicate", j, "complete", sep = " "))
#   
#   }
# 
#   scenario_fg_rli_outputs[[i]] <- replicate_fg_rli_outputs
# 
# }
#   
#   
# x <- scenario_fg_rli_outputs[[1]][[3]]
# head(x)
# 
# # Mean RLI aggregated across groups
# 
# scenario_rli_outputs <- list()
# 
# for (i in seq_along(scenario_fg_rli_outputs)) {
#   
#   # Get replicates for a single scenario
#   replicate_rli_fg <- scenario_fg_rli_outputs[[i]]
#   
#   replicate_rli_outputs <- list()
#   
#   # Aggregate RLI across functional groups for each replicate
#   for (j in seq_along(replicate_rli_fg)) {
#   
#    if ("ci_lower" %in% names(replicate_rli_fg[[j]])) {
#     
#    replicate_rli_outputs[[j]] <- replicate_rli_fg[[j]] %>%
#                                  group_by(annual_time_step) %>%
#                                  summarise(indicator_score = mean(indicator_score),
#                                            ci_lower = mean(ci_lower),
#                                            ci_upper = mean(ci_upper)) %>%
#                                  mutate(indicator = "RLI",
#                                         replicate = j)
#    } else {
#      
#    replicate_rli_outputs[[j]] <- replicate_rli_fg[[j]] %>%
#                                  group_by(annual_time_step) %>%
#                                  summarise(indicator_score = mean(indicator_score)) %>%
#                                  mutate(indicator = "RLI",
#                                         replicate = j)
#    }
# 
#   # saveRDS(replicate_rli_outputs[[j]],
#   #       file.path(rli_outputs_folder,
#   #                 paste(today, scenarios[[i]], "replicate", j,
#   #                       "RLI_aggregate_output_data.rds",
#   #                       sep = "_")))
#   # 
#   # write.csv(replicate_rli_outputs[[j]],
#   #           file.path(rli_outputs_folder,
#   #                     paste(today, scenarios[[i]], "replicate", j,
#   #                           "RLI_aggregate_output_data.rds",
#   #                           sep = "_")))
# 
#   }
#   
#   scenario_rli_outputs[[i]] <- replicate_rli_outputs
# 
# }
# 
# head(scenario_rli_outputs)[[1]][[1]]
# 
# 
# 
# # * Plot RLI ----
# 
# ## By functional group
# 
# scenario_fg_rli_plots <- list()
# 
# for (i in seq_along(scenario_fg_rli_outputs)) {
#   
#   replicate_fg_rli <- scenario_fg_rli_outputs[[i]]
#   
#   replicate_fg_rli_plots <- list()
#   
#   for (j in seq_along(replicate_fg_rli)) {
# 
#   replicate_fg_rli_plots[[j]] <-  plot_red_list_index_by_group(
#                                       replicate_fg_rli[[j]],
#                                       impact_start,
#                                       impact_end,
#                                       ci = FALSE)
# 
#   ggsave(file.path(rli_plots_folder, paste(today, scenarios[[i]], "replicate", j,
#                 "RLI_by_functional_group_annual.png",
#                 sep = "_")),
#        replicate_fg_rli_plots[[j]],  device = "png")
# 
#   }
#   
# scenario_fg_rli_plots[[i]] <- replicate_fg_rli_plots
# 
# }
# 
# scenario_fg_rli_plots[[1]][[8]]
# 
# 
# # Small test to see if averaging indicator scores after works better (it doesn't)
# x <- scenario_rli_outputs[[3]][[5]]
# x <- x[-1,]
# 
# x <- x %>% 
#      mutate(x = rollmean(indicator_score, 10, na.pad = TRUE))
# 
# ggplot(x, aes(x = annual_time_step, y = x))+
#   geom_line()
# 
# # RLI with all functional groups aggregated
# # i.e. mean of each 'taxa' RLI as per Butchart et al (2010) 'Indicators of
# # recent declines'
# 
# scenario_rli_plots <- list()
# 
# for (i in seq_along(scenario_rli_outputs)) {
#   
#   replicate_rli <- scenario_rli_outputs[[i]]
#   
#   replicate_rli_plots <- list()
#   
#   for (j in seq_along(replicate_rli)) {
# 
#     replicate_rli_plots[[j]] <- plot_red_list_index(replicate_rli[[j]],
#                                                    impact_start, 
#                                                    impact_end,
#                                                    ci = TRUE)
# 
# 
#     ggsave(file.path(rli_plots_folder, paste(today, scenarios[[i]], 
#                                              "replicate", j, 
#                                              "RLI_aggregated_annual.png",
#                                              sep = "_")),
#            replicate_rli_plots[[j]],  device = "png")                                   
# 
#   }
#   
#   scenario_rli_plots[[i]] <- replicate_rli_plots
# 
# }
# 
# i <- 1
# i <- i+1
# scenario_rli_plots[[1]][[i]]
# 
# # * Plot all replicates together ----
# 
# ## Collapse input data so RLI scores for all replicates in one scenario exist in 
# ## a single data frame
# 
# scenario_rli_outputs_aggregated <- list()
# 
# for (i in seq_along(scenario_rli_outputs)) {
#   
#   scenario_rli_outputs_aggregated[[i]] <- do.call(rbind, 
#                                                   scenario_rli_outputs[[i]]) %>%
#                                           mutate(scenario = scenarios[[i]]) 
#   
#   
#   scenario_mean_rli <- scenario_rli_outputs_aggregated[[i]] %>%
#                        group_by(annual_time_step) %>%
#                        summarise(indicator_score = mean(indicator_score),
#                                  ci_lower = mean(ci_lower),
#                                  ci_upper = mean(ci_upper)) %>%
#                        mutate(indicator = "RLI",
#                               replicate = 0,
#                               scenario = scenarios[[i]]) # Replicate 0 will always be the mean
#   
#   scenario_rli_outputs_aggregated[[i]] <-  rbind(scenario_rli_outputs_aggregated[[i]],
#                                                  scenario_mean_rli) %>%
#                                            mutate(replicate = as.factor(replicate)) %>%
#                                            mutate(level = ifelse(replicate == 0,
#                                                                  "Mean RLI", 
#                                                                  "Replicate RLI"),
#                                                   scenario = scenarios[[i]])
#   
# }
# 
# head(scenario_rli_outputs_aggregated[[1]])
# tail(scenario_rli_outputs_aggregated[[1]])
# 
# # Plot all together
# 
# scenario_rli_plots_aggregated <- list()
# 
# for (i in seq_along(scenario_rli_outputs_aggregated)) {
# 
# scenario_rli_plots_aggregated[[i]] <- ggplot(data = scenario_rli_outputs_aggregated[[i]], 
#        aes(x = annual_time_step, y = indicator_score, group = replicate,
#            color = level)) +
#   geom_line() +
#   scale_color_manual(values = c("black", "gray62")) + 
#   labs(x = "Time", 
#        y = "Red List Index Score") +
#   theme(panel.grid.major = element_blank(),
#         axis.title = element_text(size = 18),
#         axis.text = element_text(size = 18),
#         panel.grid.minor = element_blank(),
#         panel.background = element_rect(fill = "grey97"),
#         axis.line = element_line(colour = "black")) +
#   geom_vline(xintercept = impact_start, colour = "red") +
#   geom_vline(xintercept = impact_end, colour = "blue")
# 
# }
# 
# scenario_rli_plots_aggregated[[1]]
# 
# # LIVING PLANET INDEX ----
# 
# # * Create folders ----
# 
# lpi_inputs_folder <- file.path(indicator_inputs_folder, "LPI_inputs", today)
# 
# if( !dir.exists( file.path(lpi_inputs_folder) ) ) {
#   dir.create( file.path(lpi_inputs_folder), recursive = TRUE )
#   
# }
# 
# lpi_outputs_folder <- file.path(indicator_outputs_folder, "LPI_outputs", today)
# 
# if( !dir.exists( file.path(lpi_outputs_folder) ) ) {
#   dir.create( file.path(lpi_outputs_folder), recursive = TRUE )
#   
# }
# 
# lpi_plots_folder <- file.path(indicator_plots_folder, "LPI_plots", today)
# 
# if( !dir.exists( file.path(lpi_plots_folder) ) ) {
#   dir.create( file.path(lpi_plots_folder), recursive = TRUE )
#   
# }
# 
# # TEMP CODE ---
# ## Look at the data we are dealing with
# 
# # data <- scenario_abundance_long[[1]][[1]]
# # 
# # head(data)
# # 
# # ggplot(data, aes(x = time_step, y = abundance,
# #                  col = group_id)) +
# #           geom_line()  + 
# #           geom_text(aes(label= group_id),hjust=0, vjust=0) +
# #           theme(legend.position = "none")
# 
# # * Sample data ----
# 
# scenario_lpi_inputs <- list()
# 
# 
# for (i in seq_along(scenario_smoothed_abundance)) {
#   
#   # Get replicates for a single scenario
#  # replicate_abundance_long <- scenario_smoothed_abundance[[i]]
#   
#   replicate_abundance_long <- scenario_smoothed_abundance[[i]]
#   
#   replicate_lpi_inputs <- list()
#   # Calculate the LPI for each replicate within the scenario
#   for (j in seq_along(replicate_abundance_long)) {
# 
#   replicate_lpi_inputs[[j]] <- replicate_abundance_long[[j]] %>% 
#                                dplyr::select(group_id, annual_time_step, 
#                                              ave_abundance)
#     
#     
#   }
#   
#   scenario_lpi_inputs[[i]] <- replicate_lpi_inputs
# 
# }
# 
# # lpi_input <- scenario_lpi_inputs[[1]][[2]]
# # head(lpi_input)
# # write.csv(lpi_input, file.path(indicator_outputs_folder, "lpi_input_example_annual.csv"))
# 
# # * Calculate LPI ----
# 
# # Retain naming convention, the LPI just takes the abundance dataframes we
# # already formatted while making the RLI inputs
# 
# # scenario_lpi_inputs <- scenario_abundance_long
# 
# # Loop through each scenario and replicate and calculate the LPI per rep
# 
# scenario_lpi_outputs <- list()
# 
# for (i in seq_along(scenario_lpi_inputs)) {
#   
#   # Get replicates for a single scenario
#   replicate_lpi_inputs <- scenario_lpi_inputs[[i]]
#   
#   replicate_lpi_outputs <- list()
#   
#  # Calculate the LPI for each replicate within the scenario
#   for (j in seq_along(replicate_lpi_inputs)) {
#     
#   replicate_lpi_outputs[[j]] <- calculate_living_planet_index(
#     
#     replicate_lpi_inputs[[j]], start_time_step, ci = FALSE, numboots, j
#   ) 
#     
#   # Save the output LPI data as a csv and rds
#   
#     saveRDS(replicate_lpi_outputs[[j]],
#           file.path(lpi_outputs_folder,
#                     paste(today, scenarios[[i]], "replicate", j,
#                           "LPI_output_data_annual.rds",
#                           sep = "_")))
# 
#     write.csv(replicate_lpi_outputs[[j]],
#               file.path(lpi_outputs_folder,
#                         paste(today, scenarios[[i]], "replicate", j,
#                               "LPI_output_data_annual.rds",
#                               sep = "_")))
#     
#   }
#   
#   scenario_lpi_outputs[[i]] <- replicate_lpi_outputs
#   
# }
# 
# head(scenario_lpi_outputs)[[1]][[1]]
# 
# # * Aggregate all LPI scores ----
# 
# ## Collapse input data so LPI scores for all replicates in one scenario exist in 
# ## a single data frame
# 
# scenario_lpi_outputs_aggregated <- list()
# 
# for (i in seq_along(scenario_lpi_outputs)) {
#   
#   scenario_lpi_outputs_aggregated[[i]] <- do.call(rbind, 
#                                                   scenario_lpi_outputs[[i]]) %>%
#                                           mutate(scenario = scenarios[[i]]) 
#   
#   scenario_mean_lpi <- scenario_lpi_outputs_aggregated[[i]] %>%
#     group_by(annual_time_step) %>%
#     summarise(indicator_score = mean(indicator_score),
#               ci_lower = mean(ci_lower),
#               ci_upper = mean(ci_upper)) %>%
#     mutate(replicate = 0,# Replicate 0 will always be the mean
#            indicator = "LPI",
#            scenario = scenarios[[i]]) 
#   
#   scenario_lpi_outputs_aggregated[[i]] <-  rbind(scenario_lpi_outputs_aggregated[[i]],
#                                                  scenario_mean_lpi) %>%
#     mutate(replicate = as.factor(replicate)) %>%
#     mutate(level = ifelse(replicate == 0,
#                           "Mean LPI", 
#                           "Replicate LPI"))
#   
# }
# 
# head(scenario_lpi_outputs_aggregated[[1]])
# tail(scenario_lpi_outputs_aggregated[[1]])
# 
# # * Plot LPI replicates individually ----
# 
# scenario_lpi_plots <- list()
# 
# for (i in seq_along(scenario_lpi_outputs)) {
#   
#   replicate_lpi <- scenario_lpi_outputs[[i]]
#   replicate_lpi_plots <- list()
#   
#   for (j in seq_along(replicate_lpi)) {
#     
#     replicate_lpi_plots[[j]] <- plot_living_planet_index(replicate_lpi[[j]],
#                                                          ci = FALSE)
#     
#     
#     ggsave(file.path(lpi_plots_folder, paste(today, scenarios[[i]], 
#                                              "replicate", j, 
#                                              "LPI_aggregated_annual.png",
#                                              sep = "_")),
#            replicate_lpi_plots[[j]],  device = "png")                                   
#     
#   }
#   
#   scenario_lpi_plots[[i]] <- replicate_lpi_plots
#   
# }
# 
# i <- 1
# scenario_lpi_plots[[1]][[i]]
# 
# i <- i + 1
# scenario_lpi_plots[[1]][[i]]
# 
# # * Plot all replicates together ----
# 
# scenario_lpi_plots_aggregated <- list()
# 
# for (i in seq_along(scenario_lpi_outputs_aggregated)){
#   
#   scenario_lpi_plots_aggregated[[i]] <- ggplot(data = scenario_lpi_outputs_aggregated[[i]], 
#                                                aes(x = annual_time_step, 
#                                                    y = indicator_score, 
#                                                    group = replicate,
#                                                    color = level)) +
#     geom_line() +
#     scale_color_manual(values = c("black", "gray62")) + 
#     labs(x = "Time", 
#          y = "Living Planet Index Score") +
#     theme(panel.grid.major = element_blank(),
#           axis.title = element_text(size = 18),
#           axis.text = element_text(size = 18),
#           panel.grid.minor = element_blank(),
#           panel.background = element_rect(fill = "grey97"),
#           axis.line = element_line(colour = "black")) +
#     geom_vline(xintercept = impact_start, colour = "red") +
#     geom_vline(xintercept = impact_end, colour = "blue")
#   
# }
# 
# scenario_lpi_plots_aggregated[[4]]
# 
# # Combine indicators ----
# 
# all_indicators_list <- list(scenario_rli_outputs,
#                             scenario_lpi_outputs)
# 
# names(all_indicators_list) <- c("RLI", "LPI")
# 
# saveRDS(all_indicators_list,
#         file.path(indicator_outputs_folder,
#                   paste(today, "all_indicators_output_data_list_annual.rds",
#                         sep = "_")))
# 
# all_lpi <- do.call(rbind, scenario_lpi_outputs_aggregated) %>% 
#            filter(replicate != 0) # Remove the mean so we just have replicates 
# 
# all_rli <- do.call(rbind, scenario_rli_outputs_aggregated) %>% 
#            filter(replicate != 0) # Remove the mean so we just have replicates
# 
# all_indicators <- rbind(all_lpi, all_rli)
# 
# saveRDS(all_indicators,
#         file.path(indicator_outputs_folder,
#                   paste(today, "all_indicators_output_data_annual.rds",
#                         sep = "_")))
# 
# write.csv(all_indicators,
#           file.path(indicator_outputs_folder,
#                     paste(today, "all_indicators_output_data_annual.csv",
#                           sep = "_")))
# 
# ## 5 YEAR SAMPLING ### ----
# # * Merge abundance and generation length data ----
# 
# interval <- 12 * 5
# max_timestep <- 300/(interval/12)
# impact_start <- max_timestep/3 * 1  #in years
# impact_end <- max_timestep/3 * 2  #in years
# 
# 
# scenario_ab_gl_formatted_not_clean <- list()
# #scenario_ab_gl_removed <- list() # For replicates we removed from analysis
# 
# for (i in seq_along(scenario_abundance_long)) {
#   
#   replicate_abundance <- scenario_abundance_long[[i]]
#   replicate_generations <- scenario_generations_raw[[i]]
#   
#   # Make a list to catch the outputs
#   
#   replicate_ab_gl_formatted <- list()
#   
#   # For each individual replicate
#   
#   for (j in seq_along(replicate_abundance)) {
#     
#     # Reduce size of the replicate generations dataframe or the merge won't work
#     gen_length <- replicate_generations[[j]] %>% 
#       dplyr::select(group_id, generation_length_yrs, 
#                     functional_group_name) %>% 
#       distinct(.)
#     
#     # Add the generation length info to the abundance dataframe
#     temp1 <- replicate_abundance[[j]] %>%
#       merge(gen_length, by = "group_id") %>%
#       arrange(monthly_time_step, group_id) %>%
#       # Important - following lines assume an annual timeframe, will need to adjust if change interval
#       mutate(generation_by_three = generation_length_yrs * 3) %>% # Time over which to measure decline, 3 x gen length OR:
#       mutate(timeframe = ifelse(generation_by_three > gen_timeframe, # 10 years 
#                                 round(generation_by_three), gen_timeframe)) %>%
#       dplyr::select(-generation_by_three) %>%
#       distinct(.) %>%
#       group_by(group_id) %>% 
#       # select rows that are multiples of the specified interval 
#       # (eg if interval is 12, it samples one month from every 12 (yearly))
#       slice(which(row_number() %% interval == 0)) %>% 
#       mutate(annual_time_step = seq(1,max_timestep,1)) # %>% 
#     
#     # Find the last time step where non-0 abundance occurred for each group
#     
#     temp2 <- temp1 %>% 
#       group_by(group_id) %>% 
#       filter(abundance > 0) %>% 
#       dplyr::select(group_id, annual_time_step, abundance) %>% 
#       filter(annual_time_step == max(annual_time_step)) %>% 
#       dplyr::select(group_id, annual_time_step) %>% 
#       rename(last_abundance = annual_time_step)
#     
#     # Add the year of last positive abundance number as a column to the data    
#     temp3 <- temp1 %>% 
#       merge(temp2, by = c("group_id"), all = TRUE)
#     
#     # Use the last positive abundance year and current abundance value to determine
#     # if a zero abundance is a true extinction or just a missing value (false extinction)
#     data <- temp3 %>%
#       group_by(group_id) %>%
#       mutate(true_extinction = ifelse(abundance == 0 &
#                                         annual_time_step < last_abundance,
#                                       "false extinction",
#                                       ifelse(abundance > 0 &
#                                                annual_time_step < last_abundance,
#                                              "not extinct",
#                                              ifelse(abundance == 0 &
#                                                       annual_time_step >= last_abundance,
#                                                     "true extinction", "not extinct")))) %>%
#       filter(true_extinction != "false extinction") %>%
#       group_by(group_id) %>%
#       arrange(annual_time_step)
#     # 
#     # data <- temp3 %>% 
#     #   group_by(group_id) %>% 
#     #   # Identify false extinctions (where abundance = 0 but it's just missing 
#     #   # data/cohorts moving massbins)
#     #   mutate(true_extinction = ifelse(abundance == 0 & 
#     #                            annual_time_step < last_abundance,
#     #                            "false extinction",
#     #                            ifelse(abundance > 0 & 
#     #                            annual_time_step < last_abundance,
#     #                            "not extinct",
#     #                            ifelse(abundance == 0 & 
#     #                            annual_time_step >= last_abundance,
#     #                            "true extinction", "not extinct")))) %>% 
#     #   #filter(true_extinction != "false extinction") %>% 
#     #   # Convert the false zeroes to NA
#     #   mutate(abundance = ifelse(true_extinction == "false extinction",
#     #                             NA, abundance)) %>%
#     #   group_by(group_id) %>% 
#     #   arrange(annual_time_step)
#     
#     
#     # Check if there are any carnivorous endotherms
#     
#     check <- data %>% 
#       group_by(functional_group_name) %>% 
#       summarise(present = sum(abundance)) %>% 
#       filter(functional_group_name == "carnivore endotherm") %>% 
#       dplyr::select(present) %>% 
#       pull(.)
#     
#     
#     print(paste("Replicate", j - 1, 
#                 "formatting complete", 
#                 sep = " "))
#     
#     # Replace data with 0 if no carnivores
#     
#     if(length(check) == 0) {
#       
#       data <- NULL
#       
#       print(paste("Replicate", j - 1, 
#                   "removed because no carnivorous endotherms are present", 
#                   sep = " "))
#       
#     }
#     
#     replicate_ab_gl_formatted[[j]] <- data
#     
#   }
#   
#   print(scenario[[i]])
#   print(length(replicate_ab_gl_formatted))
#   
#   scenario_ab_gl_formatted_not_clean[[i]] <- replicate_ab_gl_formatted
#   
# }
# 
# 
# # Remove empty replicates (couldn't get this to work in previous loop)
# 
# scenario_ab_gl_formatted <- list()
# 
# for (i in seq_along(scenario_ab_gl_formatted_not_clean)) {
#   
#   replicate_not_clean <- scenario_ab_gl_formatted_not_clean[[i]]
#   
#   scenario_ab_gl_formatted[[i]] <- list.clean(replicate_not_clean)
#   
# }
# 
# test <- scenario_ab_gl_formatted[[1]][[1]]
# 
# # Identify and deal with weird mass bins that blink in and out
# 
# scenario_abundance_clean <- list()
# 
# for (i in seq_along(scenario_ab_gl_formatted)) {
#   
#   replicate_allgroups <- scenario_ab_gl_formatted[[i]]
#   replicate_gens <- scenario_generations_raw[[i]]
#   
#   reps_out <- list()
#   
#   for (j in seq_along(replicate_allgroups)) {
#     
#     data <- replicate_allgroups[[j]]
#     
#     gen <- replicate_gens[[j]] %>%
#       dplyr::select(group_id, mass_lower_g) %>%
#       distinct(.)
#     
#     # Determine which groups were there at beginning (post burnin)
#     temp <- data %>%
#       group_by(group_id) %>%
#       filter(monthly_time_step == min(monthly_time_step)) %>%
#       mutate(first_appearance = annual_time_step,
#              beginning = ifelse(first_appearance == 1,
#                                 TRUE, FALSE)) %>%
#       dplyr::select(group_id, first_appearance, beginning)
#     
#     reps_out[[j]] <- data %>%
#       merge(temp, by = "group_id") %>%
#       merge(gen, by = "group_id") %>%
#       arrange(mass_lower_g) %>%
#       tidylog::filter(beginning == TRUE) # Note 14% is highest percentage of data removed by this line
#     
#     rm(temp)
#   }
#   
#   scenario_abundance_clean[[i]] <- reps_out
#   
# }
# 
# # * Smooth abundance ----
# 
# #scenario_abundance_clean <- scenario_ab_gl_formatted
# 
# # Test smoothing function parameters on group being harvested
# 
# ave_window <- 10
# 
# scenario_smoothed_abundance <- list()
# 
# for (i in seq_along(scenario_abundance_clean)) {
#   
#   # Get replicate data for a single scenario
#   
#   replicate_ab_gl <- scenario_abundance_clean[[i]]
#   
#   replicate_smoothed_abundance <- list()
#   
#   # For each individual replicate
#   
#   for (j in seq_along(replicate_ab_gl)) {
#     
#     group_ab_gl <- replicate_ab_gl[[j]]
#     
#     group_list <- split(group_ab_gl, group_ab_gl$group_id)
#     
#     group_smoothed_abundance <- list()
#     
#     for (k in seq_along(group_list)) {
#       
#       group_df <- group_list[[k]]
#       
#       if (is.na(sum(group_df$abundance))) {
#         
#         group_smoothed_abundance[[k]] <- NULL
#         
#       } else {
#         
#         group_smoothed_abundance[[k]] <- group_df %>%
#           arrange(annual_time_step) %>%
#           mutate(ave_abundance = rollmean(abundance,
#                                           ave_window,
#                                           fill = NA),
#                  ave_abundance = ifelse(ave_abundance < 1,
#                                         0, ave_abundance))
#         
#         print(k)
#         
#       }
#     }
#     
#     all_groups_smooth <- do.call(rbind,group_smoothed_abundance)
#     
#     replicate_smoothed_abundance[[j]] <- all_groups_smooth
#     
#     print(j)
#   }
#   
#   scenario_smoothed_abundance[[i]] <- replicate_smoothed_abundance
#   
#   print(i)
#   
# }
# 
# check <- scenario_smoothed_abundance[[1]][[1]]
# head(check)
# 
# # * Assign Red List Categories ----
# 
# scenario_red_list_data <- list()
# 
# #for (i in seq_along(scenario_ab_gl_formatted)) {
# for (i in seq_along(scenario_smoothed_abundance)) {
#   
#   # Get replicate data for a single scenario
#   
#   replicate_ab_gl <- scenario_smoothed_abundance[[i]]
#   
#   print(paste("Processing scenario", scenarios[[i]], sep = " "))
#   
#   replicate_red_list_data <- list()
#   
#   # For each individual replicate
#   
#   for (j in seq_along(replicate_ab_gl)) {
#     
#     print(paste("Processing replicate", j, sep = " "))
#     
#     # Split by functional group, because we calculate RLI for different
#     # functional groups then aggregate later (as per Butchart etal 2010),
#     # except we are using functional groups as proxies for taxa (eg mammals, birds, 
#     # reptiles) used in real world RLI calcs
#     
#     status_inputs <- split(replicate_ab_gl[[j]], 
#                            replicate_ab_gl[[j]]$group_id)
#     
#     # Make a list to hold output for each individual massbin-func-group (ie virtual spp)
#     
#     group_red_list_data <- list()
#     
#     for (k in seq_along(status_inputs)) {
#       
#       print(paste("Processing group", names(status_inputs)[[k]], sep = " "))
#       
#       group_red_list_data[[k]] <- status_inputs[[k]] %>%
#         group_by(group_id) %>%
#         arrange(monthly_time_step) %>%
#         # calculate the difference in abundance over 10 yrs or 3 generation lengths
#         # (specified by 'timeframe' column). Its okay to take the first value of 
#         # timeframe bc the dataframe is grouped by group_id, and timeframe only changes
#         # between and not within group_ids
#         # mutate(diff = (abundance - dplyr::lag(abundance, timeframe[1]))) %>%
#         mutate(diff = (ave_abundance - dplyr::lag(ave_abundance, timeframe[1]))) %>%
#         # calculate the rate of change
#         # mutate(decline = diff/dplyr::lag(abundance, timeframe[1])) %>% 
#         mutate(decline = diff/dplyr::lag(ave_abundance, timeframe[1])) %>% 
#         # assign red list risk status based on decline 
#         mutate(rl_status = ifelse(decline > -0.40, "LC",
#                                   ifelse(decline <= -0.40 & decline > -0.50, "NT", # Where did this and LC thresholds come from?
#                                          ifelse(decline <= -0.50 & decline > -0.70, "VU",
#                                                 ifelse(decline <= -0.70 & decline > -0.90, "EN",
#                                                        ifelse(decline <= -0.90 & decline > -1, "CR",
#                                                               ifelse(decline <= -1, "EX", "NA"))))))) %>%
#         arrange(group_id, monthly_time_step) %>%
#         # Replace all non-ex status with ex after first occurrence 
#         # mutate(extinct = match("EX", rl_status)) %>%
#         mutate(extinct = ifelse(rl_status == "EX", 1, 0)) %>% 
#         # mutate(rl_status = with(., ave(rl_status, 
#         #                                         FUN=maintain_ex_status)))
#         #mutate(rl_status = rl_status) %>% 
#         group_by(group_id)
#       
#     }
#     
#     print(paste("replicate", j, "from", scenarios[[i]], "complete", sep = " "))
#     
#     replicate_red_list_df <- do.call(rbind, group_red_list_data)
#     
#     replicate_red_list_data[[j]] <- replicate_red_list_df
#     
#     # Save the inputs
#     
#     saveRDS(replicate_red_list_df,
#             file.path(rli_inputs_folder,
#                       paste(today, scenarios[[i]], "replicate", j,
#                             "RLI_input_data.rds", sep = "_")))
#     
#     write.csv(replicate_red_list_df,
#               file.path(rli_inputs_folder,
#                         paste(today, scenarios[[i]], "replicate", j,
#                               "RLI_input_data.csv", sep = "_")))
#     
#     
#   }
#   
#   scenario_red_list_data[[i]] <- replicate_red_list_data
#   
# }
# 
# # Check we have correct structure still
# length(scenario_red_list_data) == length(scenario_ab_gl_formatted)
# length(scenario_red_list_data[[1]]) == length(scenario_ab_gl_formatted[[1]])
# 
# # Have a quick look at the outputs
# 
# rli_inputs <- scenario_red_list_data[[1]][[1]]
# tail(rli_inputs)
# 
# write.csv(rli_inputs, file.path(indicator_outputs_folder, "rli_input_example_annual.csv"))
# 
# # Plot some results to check they're not completely whack
# 
# ## Get one group to check how their status changes over time relative to how
# ## their abundance changes
# 
# # group_id_select <- "13.16.17" # Shows example of 'resurrected' virtual spp
# # # group_id_select <- "10.40"
# # 
# # data <- rli_inputs %>% dplyr::filter(group_id == group_id_select)
# # 
# # ggplot(data, aes(x = time_step, y = abundance)) +
# #   geom_line() +
# #   geom_text(aes(label= rl_status,
# #                 col = rl_status),hjust=0, vjust=0)
# 
# 
# 
# 
# # * Calculate RLI ----
# 
# # RLI by individual functional groups
# 
# scenario_fg_rli_outputs <- list()
# 
# for (i in seq_along(scenario_red_list_data)) {
#   
#   replicate_red_list_inputs <- scenario_red_list_data[[i]]
#   
#   replicate_fg_rli_outputs <- list()
#   
#   for (j in seq_along(replicate_red_list_inputs)) {
#     
#     replicate_rli <- calculate_red_list_index(
#       replicate_red_list_inputs[[j]], numboots, ci = FALSE) %>%
#       mutate(replicate = j)
#     
#     replicate_fg_rli_outputs[[j]] <- replicate_rli 
#     
#     # saveRDS(replicate_fg_rli_outputs[[j]],
#     #         file.path(rli_outputs_folder,
#     #                   paste(today, scenarios, "replicate", j,
#     #                         "RLI_func_group_output_data.rds",
#     #                         sep = "_")))
#     
#     write.csv(replicate_fg_rli_outputs[[j]],
#               file.path(rli_outputs_folder,
#                         paste(today, scenarios[[i]], "RLI_func_group_output_data.rds",
#                               sep = "_")))
#     
#     print(paste("RLI for replicate", j, "complete", sep = " "))
#     
#   }
#   
#   scenario_fg_rli_outputs[[i]] <- replicate_fg_rli_outputs
#   
# }
# 
# 
# x <- scenario_fg_rli_outputs[[1]][[3]]
# head(x)
# 
# # Mean RLI aggregated across groups
# 
# scenario_rli_outputs <- list()
# 
# for (i in seq_along(scenario_fg_rli_outputs)) {
#   
#   # Get replicates for a single scenario
#   replicate_rli_fg <- scenario_fg_rli_outputs[[i]]
#   
#   replicate_rli_outputs <- list()
#   
#   # Aggregate RLI across functional groups for each replicate
#   for (j in seq_along(replicate_rli_fg)) {
#     
#     if ("ci_lower" %in% names(replicate_rli_fg[[j]])) {
#       
#       replicate_rli_outputs[[j]] <- replicate_rli_fg[[j]] %>%
#         group_by(annual_time_step) %>%
#         summarise(indicator_score = mean(indicator_score),
#                   ci_lower = mean(ci_lower),
#                   ci_upper = mean(ci_upper)) %>%
#         mutate(indicator = "RLI",
#                replicate = j)
#     } else {
#       
#       replicate_rli_outputs[[j]] <- replicate_rli_fg[[j]] %>%
#         group_by(annual_time_step) %>%
#         summarise(indicator_score = mean(indicator_score)) %>%
#         mutate(indicator = "RLI",
#                replicate = j)
#     }
#     
#     # saveRDS(replicate_rli_outputs[[j]],
#     #       file.path(rli_outputs_folder,
#     #                 paste(today, scenarios[[i]], "replicate", j,
#     #                       "RLI_aggregate_output_data.rds",
#     #                       sep = "_")))
#     # 
#     # write.csv(replicate_rli_outputs[[j]],
#     #           file.path(rli_outputs_folder,
#     #                     paste(today, scenarios[[i]], "replicate", j,
#     #                           "RLI_aggregate_output_data.rds",
#     #                           sep = "_")))
#     
#   }
#   
#   scenario_rli_outputs[[i]] <- replicate_rli_outputs
#   
# }
# 
# head(scenario_rli_outputs)[[1]][[1]]
# 
# 
# 
# # * Plot RLI ----
# 
# ## By functional group
# 
# scenario_fg_rli_plots <- list()
# 
# for (i in seq_along(scenario_fg_rli_outputs)) {
#   
#   replicate_fg_rli <- scenario_fg_rli_outputs[[i]]
#   
#   replicate_fg_rli_plots <- list()
#   
#   for (j in seq_along(replicate_fg_rli)) {
#     
#     replicate_fg_rli_plots[[j]] <-  plot_red_list_index_by_group(
#       replicate_fg_rli[[j]],
#       impact_start,
#       impact_end,
#       ci = FALSE)
#     
#     ggsave(file.path(rli_plots_folder, paste(today, scenarios[[i]], "replicate", j,
#                                              "RLI_by_functional_group_5yr.png",
#                                              sep = "_")),
#            replicate_fg_rli_plots[[j]],  device = "png")
#     
#   }
#   
#   scenario_fg_rli_plots[[i]] <- replicate_fg_rli_plots
#   
# }
# 
# scenario_fg_rli_plots[[1]][[8]]
# 
# 
# # Small test to see if averaging indicator scores after works better (it doesn't)
# x <- scenario_rli_outputs[[3]][[5]]
# x <- x[-1,]
# 
# x <- x %>% 
#   mutate(x = rollmean(indicator_score, 10, na.pad = TRUE))
# 
# ggplot(x, aes(x = annual_time_step, y = x))+
#   geom_line()
# 
# # RLI with all functional groups aggregated
# # i.e. mean of each 'taxa' RLI as per Butchart et al (2010) 'Indicators of
# # recent declines'
# 
# scenario_rli_plots <- list()
# 
# for (i in seq_along(scenario_rli_outputs)) {
#   
#   replicate_rli <- scenario_rli_outputs[[i]]
#   
#   replicate_rli_plots <- list()
#   
#   for (j in seq_along(replicate_rli)) {
#     
#     replicate_rli_plots[[j]] <- plot_red_list_index(replicate_rli[[j]],
#                                                     impact_start, 
#                                                     impact_end,
#                                                     ci = TRUE)
#     
#     
#     ggsave(file.path(rli_plots_folder, paste(today, scenarios[[i]], 
#                                              "replicate", j, 
#                                              "RLI_aggregated_5yr.png",
#                                              sep = "_")),
#            replicate_rli_plots[[j]],  device = "png")                                   
#     
#   }
#   
#   scenario_rli_plots[[i]] <- replicate_rli_plots
#   
# }
# 
# i <- 1
# i <- i+1
# scenario_rli_plots[[1]][[i]]
# 
# # * Plot all replicates together ----
# 
# ## Collapse input data so RLI scores for all replicates in one scenario exist in 
# ## a single data frame
# 
# scenario_rli_outputs_aggregated <- list()
# 
# for (i in seq_along(scenario_rli_outputs)) {
#   
#   scenario_rli_outputs_aggregated[[i]] <- do.call(rbind, 
#                                                   scenario_rli_outputs[[i]]) %>%
#     mutate(scenario = scenarios[[i]]) 
#   
#   
#   scenario_mean_rli <- scenario_rli_outputs_aggregated[[i]] %>%
#     group_by(annual_time_step) %>%
#     summarise(indicator_score = mean(indicator_score),
#               ci_lower = mean(ci_lower),
#               ci_upper = mean(ci_upper)) %>%
#     mutate(indicator = "RLI",
#            replicate = 0,
#            scenario = scenarios[[i]]) # Replicate 0 will always be the mean
#   
#   scenario_rli_outputs_aggregated[[i]] <-  rbind(scenario_rli_outputs_aggregated[[i]],
#                                                  scenario_mean_rli) %>%
#     mutate(replicate = as.factor(replicate)) %>%
#     mutate(level = ifelse(replicate == 0,
#                           "Mean RLI", 
#                           "Replicate RLI"),
#            scenario = scenarios[[i]])
#   
# }
# 
# head(scenario_rli_outputs_aggregated[[1]])
# tail(scenario_rli_outputs_aggregated[[1]])
# 
# # Plot all together
# 
# scenario_rli_plots_aggregated <- list()
# 
# for (i in seq_along(scenario_rli_outputs_aggregated)) {
#   
#   scenario_rli_plots_aggregated[[i]] <- ggplot(data = scenario_rli_outputs_aggregated[[i]], 
#                                                aes(x = annual_time_step, 
#                                                    y = indicator_score, 
#                                                    group = replicate,
#                                                    color = replicate)) +
#     geom_line(aes(linetype = level)) +
#     #scale_color_manual(values = c("black", "gray62")) + 
#     labs(x = "Time", 
#          y = "Red List Index Score") +
#     theme(panel.grid.major = element_blank(),
#           axis.title = element_text(size = 18),
#           axis.text = element_text(size = 18),
#           panel.grid.minor = element_blank(),
#           panel.background = element_rect(fill = "grey97"),
#           axis.line = element_line(colour = "black")) +
#     geom_vline(xintercept = impact_start, colour = "red") +
#     geom_vline(xintercept = impact_end, colour = "blue") +
#     scale_y_continuous(limits = c(0,1))
#   
# }
# 
# scenario_rli_plots_aggregated[[1]]
# 
# # LIVING PLANET INDEX ----
# 
# # * Create folders ----
# 
# lpi_inputs_folder <- file.path(indicator_inputs_folder, "LPI_inputs", today)
# 
# if( !dir.exists( file.path(lpi_inputs_folder) ) ) {
#   dir.create( file.path(lpi_inputs_folder), recursive = TRUE )
#   
# }
# 
# lpi_outputs_folder <- file.path(indicator_outputs_folder, "LPI_outputs", today)
# 
# if( !dir.exists( file.path(lpi_outputs_folder) ) ) {
#   dir.create( file.path(lpi_outputs_folder), recursive = TRUE )
#   
# }
# 
# lpi_plots_folder <- file.path(indicator_plots_folder, "LPI_plots", today)
# 
# if( !dir.exists( file.path(lpi_plots_folder) ) ) {
#   dir.create( file.path(lpi_plots_folder), recursive = TRUE )
#   
# }
# 
# # TEMP CODE ---
# ## Look at the data we are dealing with
# 
# # data <- scenario_abundance_long[[1]][[1]]
# # 
# # head(data)
# # 
# # ggplot(data, aes(x = time_step, y = abundance,
# #                  col = group_id)) +
# #           geom_line()  + 
# #           geom_text(aes(label= group_id),hjust=0, vjust=0) +
# #           theme(legend.position = "none")
# 
# # * Sample data ----
# 
# scenario_lpi_inputs <- list()
# 
# 
# for (i in seq_along(scenario_smoothed_abundance)) {
#   
#   # Get replicates for a single scenario
#   # replicate_abundance_long <- scenario_smoothed_abundance[[i]]
#   
#   replicate_abundance_long <- scenario_smoothed_abundance[[i]]
#   
#   replicate_lpi_inputs <- list()
#   # Calculate the LPI for each replicate within the scenario
#   for (j in seq_along(replicate_abundance_long)) {
#     
#     replicate_lpi_inputs[[j]] <- replicate_abundance_long[[j]] %>%
#       dplyr::select(group_id, annual_time_step,
#                     ave_abundance)
#     
#   #  lpi_incomplete  <- replicate_abundance_long[[j]] %>%
#   #     dplyr::select(group_id, annual_time_step,
#   #                   ave_abundance)
#   #  
#   # replicate_lpi_inputs[[j]] <- lpi_incomplete[complete.cases(lpi_incomplete),]
#   
#   }
#   
#   scenario_lpi_inputs[[i]] <- replicate_lpi_inputs
#   
# }
# 
# # lpi_input <- scenario_lpi_inputs[[1]][[2]]
# # head(lpi_input)
# # write.csv(lpi_input, file.path(indicator_outputs_folder, "lpi_input_example_5yr.csv"))
# 
# # * Calculate LPI ----
# 
# # Retain naming convention, the LPI just takes the abundance dataframes we
# # already formatted while making the RLI inputs
# 
# # scenario_lpi_inputs <- scenario_abundance_long
# 
# # Loop through each scenario and replicate and calculate the LPI per rep
# 
# scenario_lpi_outputs <- list()
# 
# for (i in seq_along(scenario_lpi_inputs)) {
#   
#   # Get replicates for a single scenario
#   replicate_lpi_inputs <- scenario_lpi_inputs[[i]]
#   
#   replicate_lpi_outputs <- list()
#   
#   # Calculate the LPI for each replicate within the scenario
#   for (j in seq_along(replicate_lpi_inputs)) {
#     
#     replicate_lpi_outputs[[j]] <- calculate_living_planet_index(
#       
#       replicate_lpi_inputs[[j]], start_time_step, ci = FALSE, numboots, j
#     ) 
#     
#     # Save the output LPI data as a csv and rds
#     
#     saveRDS(replicate_lpi_outputs[[j]],
#             file.path(lpi_outputs_folder,
#                       paste(today, scenarios[[i]], "replicate", j,
#                             "LPI_output_data_5yr.rds",
#                             sep = "_")))
#     
#     write.csv(replicate_lpi_outputs[[j]],
#               file.path(lpi_outputs_folder,
#                         paste(today, scenarios[[i]], "replicate", j,
#                               "LPI_output_data_5yr.rds",
#                               sep = "_")))
#     
#   }
#   
#   scenario_lpi_outputs[[i]] <- replicate_lpi_outputs
#   
# }
# 
# head(scenario_lpi_outputs)[[1]][[1]]
# x <- scenario_lpi_outputs[[1]][[1]]
# # * Aggregate all LPI scores ----
# 
# ## Collapse input data so LPI scores for all replicates in one scenario exist in 
# ## a single data frame
# 
# scenario_lpi_outputs_aggregated <- list()
# 
# for (i in seq_along(scenario_lpi_outputs)) {
#   
#   scenario_lpi_outputs_aggregated[[i]] <- do.call(rbind, 
#                                                   scenario_lpi_outputs[[i]]) %>%
#     mutate(scenario = scenarios[[i]]) 
#   
#   scenario_mean_lpi <- scenario_lpi_outputs_aggregated[[i]] %>%
#     group_by(annual_time_step) %>%
#     summarise(indicator_score = mean(indicator_score),
#               ci_lower = mean(ci_lower),
#               ci_upper = mean(ci_upper)) %>%
#     mutate(replicate = 0,# Replicate 0 will always be the mean
#            indicator = "LPI",
#            scenario = scenarios[[i]]) 
#   
#   scenario_lpi_outputs_aggregated[[i]] <-  rbind(scenario_lpi_outputs_aggregated[[i]],
#                                                  scenario_mean_lpi) %>%
#     mutate(replicate = as.factor(replicate)) %>%
#     mutate(level = ifelse(replicate == 0,
#                           "Mean LPI", 
#                           "Replicate LPI"))
#   
# }
# 
# head(scenario_lpi_outputs_aggregated[[1]])
# tail(scenario_lpi_outputs_aggregated[[1]])
# 
# # * Plot LPI replicates individually ----
# 
# scenario_lpi_plots <- list()
# 
# for (i in seq_along(scenario_lpi_outputs)) {
#   
#   replicate_lpi <- scenario_lpi_outputs[[i]]
#   replicate_lpi_plots <- list()
#   
#   for (j in seq_along(replicate_lpi)) {
#     
#     replicate_lpi_plots[[j]] <- plot_living_planet_index(replicate_lpi[[j]],
#                                                          ci = FALSE)
#     
#     
#     ggsave(file.path(lpi_plots_folder, paste(today, scenarios[[i]], 
#                                              "replicate", j, 
#                                              "LPI_aggregated_5yr.png",
#                                              sep = "_")),
#            replicate_lpi_plots[[j]],  device = "png")                                   
#     
#   }
#   
#   scenario_lpi_plots[[i]] <- replicate_lpi_plots
#   
# }
# 
# i <- 1
# scenario_lpi_plots[[1]][[i]]
# 
# i <- i + 1
# scenario_lpi_plots[[1]][[i]]
# 
# # * Plot all replicates together ----
# 
# scenario_lpi_plots_aggregated <- list()
# 
# for (i in seq_along(scenario_lpi_outputs_aggregated)){
#   
#   scenario_lpi_plots_aggregated[[i]] <- ggplot(data = scenario_lpi_outputs_aggregated[[i]], 
#                                                aes(x = annual_time_step, 
#                                                    y = indicator_score, 
#                                                    color = replicate)) +
#     geom_line(aes(linetype = level)) +
#     #scale_color_manual(values = c("black", "gray62")) + 
#     labs(x = "Time", 
#          y = "Living Planet Index Score") +
#     theme(panel.grid.major = element_blank(),
#           axis.title = element_text(size = 18),
#           axis.text = element_text(size = 18),
#           panel.grid.minor = element_blank(),
#           panel.background = element_rect(fill = "grey97"),
#           axis.line = element_line(colour = "black")) +
#     geom_vline(xintercept = impact_start, colour = "red") +
#     geom_vline(xintercept = impact_end, colour = "blue")
#   
# }
# 
# scenario_lpi_plots_aggregated[[1]]
# 
# # Combine indicators ----
# 
# all_indicators_list <- list(scenario_rli_outputs,
#                             scenario_lpi_outputs)
# 
# names(all_indicators_list) <- c("RLI", "LPI")
# 
# saveRDS(all_indicators_list,
#         file.path(indicator_outputs_folder,
#                   paste(today, "all_indicators_output_data_list_5yr.rds",
#                         sep = "_")))
# 
# all_lpi <- do.call(rbind, scenario_lpi_outputs_aggregated) %>% 
#   filter(replicate != 0) # Remove the mean so we just have replicates 
# 
# all_rli <- do.call(rbind, scenario_rli_outputs_aggregated) %>% 
#   filter(replicate != 0) # Remove the mean so we just have replicates
# 
# all_indicators <- rbind(all_lpi, all_rli)
# 
# saveRDS(all_indicators,
#         file.path(indicator_outputs_folder,
#                   paste(today, "all_indicators_output_data_5yr.rds",
#                         sep = "_")))
# 
# write.csv(all_indicators,
#           file.path(indicator_outputs_folder,
#                     paste(today, "all_indicators_output_data_5yr.csv",
#                           sep = "_")))
# 
# ## MERGE FUNCTIONAL GROUPS ----
# 
# calculate_red_list_index2 <- function(data, numboots, ci = FALSE, replicate_num = NA){
#   
#   # Using equation from Butchart et al (2007) Improvements to the Red List Index
#   
#   require(tidyverse)
#   
#   # Remove data without RL status
#   
#   #data$redlist_assessment_year <- as.numeric(as.character(data$redlist_assessment_year))
#   
#   data <- data %>%
#     filter(!is.na(rl_status)) %>%
#     group_by(group_id) 
#   
#   head(data)
#   
#   # ecoregion <- as.factor(data$ecoregion_id[1])
#   
#   # Assign category weights
#   
#   weighted_data <- data %>%
#     dplyr::mutate(rl_weight = ifelse(rl_status == "LC", 0,
#                               ifelse(rl_status == "NT", 1,
#                               ifelse(rl_status == "VU", 2,
#                               ifelse(rl_status == "EN", 3,
#                               ifelse(rl_status == "CR", 4,
#                               ifelse(rl_status == "EX", 5, NA))))))) 
#   head(weighted_data)
#   dim(weighted_data)
#   
#   
#   #weighted_data$RL_weight <- as.numeric(as.character(weighted_data$RL_weight))
#   
#   # Filter out rows with NE and DD
#   weighted_data <- weighted_data %>%
#     filter(rl_status != "NE") %>%
#     filter(rl_status != "DD") %>%
#     filter(rl_status != "NA")
#   
#   dim(weighted_data)
#   
#   # Group data so the index is calculated for each functional group 
#   # (would normally be taxa) for each year. If you run on a single group
#   # it shouldn't matter, will just turn data into one big group
#   
#   grouped_data <- weighted_data %>% group_by(group_id, annual_time_step)
#   
#   # Sum category weights for each group, in each timestep,
#   # calculate number of species per group
#   summed_weights <- summarise(grouped_data, 
#                               total_weight = sum(rl_weight, na.rm = TRUE), # calc sum of all weights
#                               total_count = n(),# calc number of species
#                               .groups = "drop_last") %>%
#     mutate(total_count = max(total_count))  # Fix so it takes total number at beginning, otherwise n fluctuates between timesteps
#   
#   # Calculate RLI scores for each group, rounded to 3 decimal places
#   
#   index_scores <- summed_weights %>%
#     mutate(RLI = 1 - (total_weight/(total_count * 5)), # actual RLI formula
#            Criteria = "risk")
#   
#   if (ci == TRUE) {
#     # Calculate confidence intervals via bootstrapping 
#     # (see Rowland et al 2021 A guide to representing uncertainty)
#     
#     # Split by timestep - we want CI for each functional group, for each timestep 
#     
#     weighted_data_timestep_list <- split(weighted_data, weighted_data$annual_time_step)
#     
#     ## For each functional group (level 1)
#     
#     timestep_confidence_intervals <- list()
#     
#     for (i in seq_along(weighted_data_timestep_list)) {
#       
#       # Get single time-step then group by functional group
#       
#       grouped_timestep_data <- weighted_data_timestep_list[[i]] %>%
#         group_by(group_id)
#       
#       time <- grouped_timestep_data$annual_time_step[1]
#       
#       boot <- list()
#       # Calculate the bootstrap confidence intervals
#       for (k in 1:numboots) {
#         
#         # Take k number of random samples from the weighted data
#         replicate <- slice_sample(grouped_timestep_data, 
#                                   prop = 1, replace = TRUE) %>%  # get random sample of rows and add to DF
#           mutate(replicate = k)
#         # label each replicate
#         
#         boot[[k]] <- replicate
#         # Combine replicates into one dataframe
#         
#         # print(paste("Bootstrap", k, "of", numboots, "complete", sep =" "))
#         
#       }
#       
#       boot_reps <- do.call(rbind, boot)
#       
#       # Group by replicate
#       #replicate_data <- group_by(boot_reps, replicate) # Group by replicate
#       
#       # Calculate the summary values needed to calc RLI for each replicate
#       summed_weights_timestep_fg <- boot_reps %>%
#         group_by(group_id,replicate) %>% 
#         summarise(total_weight = sum(rl_weight, na.rm = TRUE), # calc sum of all weights
#                   total_count = n(),# calc number of species
#                   .groups = "drop_last") %>%
#         mutate(total_count = max(total_count))
#       
#       # Calculate the RLI score for each replicate
#       rep_scores <- mutate(summed_weights_timestep_fg, 
#                            RLI = 1 - (total_weight/(total_count * 5))) # actual RLI formula
#       
#       # Calculate the confidence intervals for each fg,
#       ci_scores <- summarise(rep_scores, 
#                              ci_lower = quantile(rep_scores$RLI, 
#                                                  probs = 0.025),
#                              ci_upper = quantile(rep_scores$RLI, 
#                                                  probs = 0.975)) %>%
#         mutate(annual_time_step = time) 
#       
#       timestep_confidence_intervals[[i]] <- ci_scores
#       
#     }
#     
#     confidence_intervals <- do.call(rbind, timestep_confidence_intervals)
#     
#     red_list_scores <- index_scores %>%
#       merge(confidence_intervals, 
#             by = c("group_id",
#                    "annual_time_step")) %>%
#       dplyr::select(group_id, annual_time_step, ci_lower,
#                     RLI, ci_upper, everything()) %>% 
#       rename(indicator_score = RLI) %>% 
#       mutate(indicator = "RLI",
#              replicate = replicate_num)
#     
#     
#     return(red_list_scores)
#     
#   } else {
#     
#     red_list_scores <- index_scores  %>% 
#       rename(indicator_score = RLI) %>% 
#       mutate(indicator = "RLI",
#              replicate = replicate_num,
#              ci_lower = NA,
#              ci_upper = NA) %>% 
#       dplyr::select(group_id, annual_time_step,
#                     ci_lower, indicator_score, ci_upper, total_weight,
#                     total_count, Criteria, indicator, replicate) 
#     
#     
#     return(red_list_scores)
#     
#   }
# }
# 
# burnin_months <- 1000*12 # in months
# n <- 12
# numboots <- 1000 # Rowland et al 2021 (uncertainty)
# start_time_step <- 1
# gen_timeframe <- 10 
# interval <- 12 
# # Don't adjust these
# max_timestep <- 300/(interval/12)
# impact_start <- max_timestep/3 * 1  #in years
# impact_end <- max_timestep/3 * 2  #in years
# 
# # * Merge abundance and generation length data ----
# 
# scenario_ab_gl_formatted_not_clean <- list()
# 
# for (i in seq_along(scenario_abundance_long)) {
#   
#   replicate_abundance <- scenario_abundance_long[[i]]
#   replicate_generations <- scenario_generations_raw[[i]]
#   
#   # Make a list to catch the outputs
#   
#   replicate_ab_gl_formatted <- list()
#   
#   # For each individual replicate
#   
#   for (j in seq_along(replicate_abundance)) {
#     
#     # Reduce size of the replicate generations dataframe or the merge won't work
#     gen_length <- replicate_generations[[j]] %>% 
#       dplyr::select(group_id, generation_length_yrs, 
#                     functional_group_name) %>% 
#       distinct(.)
#     
#     # Add the generation length info to the abundance dataframe
#     temp1 <- replicate_abundance[[j]] %>%
#       merge(gen_length, by = "group_id") %>%
#       arrange(monthly_time_step, group_id) %>%
#       # Important - following lines assume an annual timeframe, will need to adjust if change interval
#       mutate(generation_by_three = generation_length_yrs * 3) %>% # Time over which to measure decline, 3 x gen length OR:
#       mutate(timeframe = ifelse(generation_by_three > gen_timeframe, # 10 years 
#                                 round(generation_by_three), gen_timeframe)) %>%
#       dplyr::select(-generation_by_three) %>%
#       distinct(.) %>%
#       group_by(group_id) %>% 
#       # select rows that are multiples of the specified interval 
#       # (eg if interval is 12, it samples one month from every 12 (yearly))
#       slice(which(row_number() %% interval == 0)) %>% 
#       mutate(annual_time_step = seq(1,max_timestep,1)) %>% 
#       dplyr::select(- generation_length_yrs) %>% 
#       distinct(.) %>% 
#       group_by(functional_group_name, annual_time_step) %>% 
#       mutate(fg_abundance = sum(abundance, na.rm = TRUE)) %>% 
#       dplyr::select(-group_id, -abundance) %>% 
#       rename(abundance = fg_abundance,
#              group_id = functional_group_name)
#     
#     head(temp1)
#     
#     # Find the last time step where non-0 abundance occurred for each group
#     
#     temp2 <- temp1 %>% 
#       group_by(group_id) %>% 
#       filter(abundance > 0) %>% 
#       dplyr::select(group_id, annual_time_step, abundance) %>% 
#       filter(annual_time_step == max(annual_time_step)) %>% 
#       dplyr::select(group_id, annual_time_step) %>% 
#       rename(last_abundance = annual_time_step)
#     
#     # Add the year of last positive abundance number as a column to the data    
#     temp3 <- temp1 %>% 
#       merge(temp2, by = c("group_id"), all = TRUE)
#     
#     # Use the last positive abundance year and current abundance value to determine
#     # if a zero abundance is a true extinction or just a missing value (false extinction)
#     data <- temp3 %>%
#       group_by(group_id) %>%
#       mutate(true_extinction = ifelse(abundance == 0 &
#                                    annual_time_step < last_abundance,
#                                       "false extinction",
#                                       ifelse(abundance > 0 &
#                                                annual_time_step < last_abundance,
#                                              "not extinct",
#                                              ifelse(abundance == 0 &
#                                                       annual_time_step >= last_abundance,
#                                                     "true extinction", "not extinct")))) %>%
#       filter(true_extinction != "false extinction") %>%
#       group_by(group_id) %>%
#       arrange(annual_time_step) %>% 
#       distinct(.)
#     
#     # data <- temp3 %>% 
#     #   group_by(group_id) %>% 
#     #   # Identify false extinctions (where abundance = 0 but it's just missing 
#     #   # data/cohorts moving massbins)
#     #   mutate(true_extinction = ifelse(abundance == 0 & 
#     #                            annual_time_step < last_abundance,
#     #                            "false extinction",
#     #                            ifelse(abundance > 0 & 
#     #                            annual_time_step < last_abundance,
#     #                            "not extinct",
#     #                            ifelse(abundance == 0 & 
#     #                            annual_time_step >= last_abundance,
#     #                            "true extinction", "not extinct")))) %>% 
#     #   #filter(true_extinction != "false extinction") %>% 
#     #   # Convert the false zeroes to NA
#     #   mutate(abundance = ifelse(true_extinction == "false extinction",
#     #                             NA, abundance)) %>%
#     #   group_by(group_id) %>% 
#     #   arrange(annual_time_step)
#     
#     
#     # Check if there are any carnivorous endotherms
#     
#     check <- data %>% 
#       group_by(group_id) %>% 
#       summarise(present = sum(abundance, na.rm = TRUE)) %>% 
#       filter(group_id == "carnivore endotherm") %>% 
#       dplyr::select(present) %>% 
#       pull(.)
#     
#     
#     print(paste("Replicate", j - 1, 
#                 "formatting complete", 
#                 sep = " "))
#     
#     # Replace data with 0 if no carnivores
#     
#     if(length(check == 0)) {
#       
#       data <- NULL
#       
#       print(paste("Replicate", j - 1, 
#                   "removed because no carnivorous endotherms are present", 
#                   sep = " "))
#       
#     }
#     
#     replicate_ab_gl_formatted[[j]] <- data
#     
#   }
#   
#   print(scenario[[i]])
#   print(length(replicate_ab_gl_formatted))
#   
#   scenario_ab_gl_formatted_not_clean[[i]] <- replicate_ab_gl_formatted
#   
# }
# 
# 
# # Remove empty replicates (couldn't get this to work in previous loop)
# 
# scenario_ab_gl_formatted <- list()
# 
# for (i in seq_along(scenario_ab_gl_formatted_not_clean)) {
#   
#   replicate_not_clean <- scenario_ab_gl_formatted_not_clean[[i]]
#   
#   scenario_ab_gl_formatted[[i]] <- list.clean(replicate_not_clean)
#   
# }
# 
# test <- scenario_ab_gl_formatted[[1]][[1]]
# head(test)
# 
# # * Smooth abundance ----
# 
# scenario_abundance_clean <- scenario_ab_gl_formatted
# 
# # Test smoothing function parameters on group being harvested
# 
# ave_window <- 10
# 
# scenario_smoothed_abundance <- list()
# 
# for (i in seq_along(scenario_abundance_clean)) {
#   
#   # Get replicate data for a single scenario
#   
#   replicate_ab_gl <- scenario_abundance_clean[[i]]
#   
#   replicate_smoothed_abundance <- list()
#   
#   # For each individual replicate
#   
#   for (j in seq_along(replicate_ab_gl)) {
#     
#     group_ab_gl <- replicate_ab_gl[[j]]
#     
#     group_list <- split(group_ab_gl, group_ab_gl$group_id)
#     
#     group_smoothed_abundance <- list()
#     
#     for (k in seq_along(group_list)) {
#       
#       group_df <- group_list[[k]]
#       
#       if (is.na(sum(group_df$abundance))) {
#         
#         group_smoothed_abundance[[k]] <- NULL
#         
#       } else {
#         
#         group_smoothed_abundance[[k]] <- group_df %>%
#           arrange(annual_time_step) %>%
#           mutate(ave_abundance = rollmean(abundance,
#                                           ave_window,
#                                           fill = "extend"),
#                  ave_abundance = ifelse(ave_abundance < 1,
#                                         0, ave_abundance))
#         
#         print(k)
#         
#       }
#     }
#     
#     all_groups_smooth <- do.call(rbind,group_smoothed_abundance)
#     
#     replicate_smoothed_abundance[[j]] <- all_groups_smooth
#     
#     print(j)
#   }
#   
#   scenario_smoothed_abundance[[i]] <- replicate_smoothed_abundance
#   
#   print(i)
#   
# }
# 
# check <- scenario_smoothed_abundance[[1]][[1]]
# 
# # * Assign Red List Categories ----
# 
# scenario_red_list_data <- list()
# 
# #for (i in seq_along(scenario_ab_gl_formatted)) {
# for (i in seq_along(scenario_smoothed_abundance)) {
#   
#   # Get replicate data for a single scenario
#   
#   replicate_ab_gl <- scenario_smoothed_abundance[[i]]
#   
#   print(paste("Processing scenario", scenarios[[i]], sep = " "))
#   
#   replicate_red_list_data <- list()
#   
#   # For each individual replicate
#   
#   for (j in seq_along(replicate_ab_gl)) {
#     
#     print(paste("Processing replicate", j, sep = " "))
#     
#     # Split by functional group, because we calculate RLI for different
#     # functional groups then aggregate later (as per Butchart etal 2010),
#     # except we are using functional groups as proxies for taxa (eg mammals, birds, 
#     # reptiles) used in real world RLI calcs
#     
#     status_inputs <- split(replicate_ab_gl[[j]], 
#                            replicate_ab_gl[[j]]$group_id)
#     
#     # Make a list to hold output for each individual massbin-func-group (ie virtual spp)
#     
#     group_red_list_data <- list()
#     
#     for (k in seq_along(status_inputs)) {
#       
#       print(paste("Processing group", names(status_inputs)[[k]], sep = " "))
#       
#       group_red_list_data[[k]] <- status_inputs[[k]] %>%
#         group_by(group_id) %>%
#         arrange(monthly_time_step) %>%
#         # calculate the difference in abundance over 10 yrs or 3 generation lengths
#         # (specified by 'timeframe' column). Its okay to take the first value of 
#         # timeframe bc the dataframe is grouped by group_id, and timeframe only changes
#         # between and not within group_ids
#         # mutate(diff = (abundance - dplyr::lag(abundance, timeframe[1]))) %>%
#         mutate(diff = (ave_abundance - dplyr::lag(ave_abundance, timeframe[1]))) %>%
#         # calculate the rate of change
#         # mutate(decline = diff/dplyr::lag(abundance, timeframe[1])) %>% 
#         mutate(decline = diff/dplyr::lag(ave_abundance, timeframe[1])) %>% 
#         # assign red list risk status based on decline 
#         mutate(rl_status = ifelse(decline > -0.40, "LC",
#                            ifelse(decline <= -0.40 & decline > -0.50, "NT", # Where did this and LC thresholds come from?
#                            ifelse(decline <= -0.50 & decline > -0.70, "VU",
#                            ifelse(decline <= -0.70 & decline > -0.90, "EN",
#                            ifelse(decline <= -0.90 & decline > -1, "CR",
#                            ifelse(decline <= -1, "EX", "NA"))))))) %>%
#         arrange(group_id, monthly_time_step) %>%
#         # Replace all non-ex status with ex after first occurrence 
#         # mutate(extinct = match("EX", rl_status)) %>%
#         mutate(extinct = ifelse(rl_status == "EX", 1, 0)) %>% 
#         # mutate(rl_status = with(., ave(rl_status, 
#         #                                         FUN=maintain_ex_status)))
#         #mutate(rl_status = rl_status) %>% 
#         group_by(group_id)
#       
#     }
#     
#     print(paste("replicate", j, "from", scenarios[[i]], "complete", sep = " "))
#     
#     replicate_red_list_df <- do.call(rbind, group_red_list_data)
#     
#     replicate_red_list_data[[j]] <- replicate_red_list_df
#     
#     # Save the inputs
#     
#     saveRDS(replicate_red_list_df,
#             file.path(rli_inputs_folder,
#                       paste(today, scenarios[[i]], "replicate", j,
#                             "RLI_input_data_func_groups.rds", sep = "_")))
#     
#     write.csv(replicate_red_list_df,
#               file.path(rli_inputs_folder,
#                         paste(today, scenarios[[i]], "replicate", j,
#                               "RLI_input_data_func_groups.csv", sep = "_")))
#     
#     
#   }
#   
#   scenario_red_list_data[[i]] <- replicate_red_list_data
#   
# }
# 
# # Check we have correct structure still
# length(scenario_red_list_data) == length(scenario_ab_gl_formatted)
# length(scenario_red_list_data[[1]]) == length(scenario_ab_gl_formatted[[1]])
# 
# # Have a quick look at the outputs
# 
# rli_inputs <- scenario_red_list_data[[1]][[1]]
# tail(rli_inputs)
# 
# write.csv(rli_inputs, file.path(indicator_outputs_folder, 
#                                 "rli_input_example_func_groups.csv"))
# 
# # Plot some results to check they're not completely whack
# 
# ## Get one group to check how their status changes over time relative to how
# ## their abundance changes
# 
# # group_id_select <- "13.16.17" # Shows example of 'resurrected' virtual spp
# # # group_id_select <- "10.40"
# # 
# # data <- rli_inputs %>% dplyr::filter(group_id == group_id_select)
# # 
# # ggplot(data, aes(x = time_step, y = abundance)) +
# #   geom_line() +
# #   geom_text(aes(label= rl_status,
# #                 col = rl_status),hjust=0, vjust=0)
# 
# 
# 
# 
# # * Calculate RLI ----
# 
# # RLI by individual functional groups
# 
# scenario_fg_rli_outputs <- list()
# 
# for (i in seq_along(scenario_red_list_data)) {
#   
#   replicate_red_list_inputs <- scenario_red_list_data[[i]]
#   
#   replicate_fg_rli_outputs <- list()
#   
#   for (j in seq_along(replicate_red_list_inputs)) {
#     
#     replicate_rli <- calculate_red_list_index2(
#       replicate_red_list_inputs[[j]], numboots, ci = FALSE) %>%
#       mutate(replicate = j)
#     
#     replicate_fg_rli_outputs[[j]] <- replicate_rli 
#     
#     # saveRDS(replicate_fg_rli_outputs[[j]],
#     #         file.path(rli_outputs_folder,
#     #                   paste(today, scenarios, "replicate", j,
#     #                         "RLI_func_group_output_data.rds",
#     #                         sep = "_")))
#     
#     write.csv(replicate_fg_rli_outputs[[j]],
#               file.path(rli_outputs_folder,
#                         paste(today, scenarios[[i]], 
#                               "RLI_func_group_output_data_func_groups.rds",
#                               sep = "_")))
#     
#     print(paste("RLI for replicate", j, "complete", sep = " "))
#     
#   }
#   
#   scenario_fg_rli_outputs[[i]] <- replicate_fg_rli_outputs
#   
# }
# 
# 
# x <- scenario_fg_rli_outputs[[1]][[3]]
# head(x)
# 
# # Mean RLI aggregated across groups
# 
# scenario_rli_outputs <- list()
# 
# for (i in seq_along(scenario_fg_rli_outputs)) {
#   
#   # Get replicates for a single scenario
#   replicate_rli_fg <- scenario_fg_rli_outputs[[i]]
#   
#   replicate_rli_outputs <- list()
#   
#   # Aggregate RLI across functional groups for each replicate
#   for (j in seq_along(replicate_rli_fg)) {
#     
#     if ("ci_lower" %in% names(replicate_rli_fg[[j]])) {
#       
#       replicate_rli_outputs[[j]] <- replicate_rli_fg[[j]] %>%
#         group_by(annual_time_step) %>%
#         summarise(indicator_score = mean(indicator_score),
#                   ci_lower = mean(ci_lower),
#                   ci_upper = mean(ci_upper)) %>%
#         mutate(indicator = "RLI",
#                replicate = j)
#     } else {
#       
#       replicate_rli_outputs[[j]] <- replicate_rli_fg[[j]] %>%
#         group_by(annual_time_step) %>%
#         summarise(indicator_score = mean(indicator_score)) %>%
#         mutate(indicator = "RLI",
#                replicate = j)
#     }
#     
#     # saveRDS(replicate_rli_outputs[[j]],
#     #       file.path(rli_outputs_folder,
#     #                 paste(today, scenarios[[i]], "replicate", j,
#     #                       "RLI_aggregate_output_data.rds",
#     #                       sep = "_")))
#     # 
#     # write.csv(replicate_rli_outputs[[j]],
#     #           file.path(rli_outputs_folder,
#     #                     paste(today, scenarios[[i]], "replicate", j,
#     #                           "RLI_aggregate_output_data.rds",
#     #                           sep = "_")))
#     
#   }
#   
#   scenario_rli_outputs[[i]] <- replicate_rli_outputs
#   
# }
# 
# head(scenario_rli_outputs)[[1]][[1]]
# 
# 
# 
# # * Plot RLI ----
# 
# ## By functional group
# 
# scenario_fg_rli_plots <- list()
# 
# for (i in seq_along(scenario_fg_rli_outputs)) {
#   
#   replicate_fg_rli <- scenario_fg_rli_outputs[[i]]
#   
#   replicate_fg_rli_plots <- list()
#   
#   for (j in seq_along(replicate_fg_rli)) {
#     
#     replicate_fg_rli_plots[[j]] <-  plot_red_list_index_by_group(
#       replicate_fg_rli[[j]],
#       impact_start,
#       impact_end,
#       ci = FALSE)
#     
#     ggsave(file.path(rli_plots_folder, paste(today, scenarios[[i]], "replicate", j,
#                                              "RLI_by_functional_group_func_groups.png",
#                                              sep = "_")),
#            replicate_fg_rli_plots[[j]],  device = "png")
#     
#   }
#   
#   scenario_fg_rli_plots[[i]] <- replicate_fg_rli_plots
#   
# }
# 
# scenario_fg_rli_plots[[1]][[8]]
# 
# # RLI with all functional groups aggregated
# # i.e. mean of each 'taxa' RLI as per Butchart et al (2010) 'Indicators of
# # recent declines'
# 
# scenario_rli_plots <- list()
# 
# for (i in seq_along(scenario_rli_outputs)) {
#   
#   replicate_rli <- scenario_rli_outputs[[i]]
#   
#   replicate_rli_plots <- list()
#   
#   for (j in seq_along(replicate_rli)) {
#     
#     replicate_rli_plots[[j]] <- plot_red_list_index(replicate_rli[[j]],
#                                                     impact_start, 
#                                                     impact_end,
#                                                     ci = TRUE)
#     
#     
#     ggsave(file.path(rli_plots_folder, paste(today, scenarios[[i]], 
#                                              "replicate", j, 
#                                              "RLI_aggregated_func_groups.png",
#                                              sep = "_")),
#            replicate_rli_plots[[j]],  device = "png")                                   
#     
#   }
#   
#   scenario_rli_plots[[i]] <- replicate_rli_plots
#   
# }
# 
# i <- 1
# i <- i+1
# scenario_rli_plots[[1]][[i]]
# 
# # * Plot all replicates together ----
# 
# ## Collapse input data so RLI scores for all replicates in one scenario exist in 
# ## a single data frame
# 
# scenario_rli_outputs_aggregated <- list()
# 
# for (i in seq_along(scenario_rli_outputs)) {
#   
#   scenario_rli_outputs_aggregated[[i]] <- do.call(rbind, 
#                                                   scenario_rli_outputs[[i]]) %>%
#     mutate(scenario = scenarios[[i]]) 
#   
#   
#   scenario_mean_rli <- scenario_rli_outputs_aggregated[[i]] %>%
#     group_by(annual_time_step) %>%
#     summarise(indicator_score = mean(indicator_score),
#               ci_lower = mean(ci_lower),
#               ci_upper = mean(ci_upper)) %>%
#     mutate(indicator = "RLI",
#            replicate = 0,
#            scenario = scenarios[[i]]) # Replicate 0 will always be the mean
#   
#   scenario_rli_outputs_aggregated[[i]] <-  rbind(scenario_rli_outputs_aggregated[[i]],
#                                                  scenario_mean_rli) %>%
#     mutate(replicate = as.factor(replicate)) %>%
#     mutate(level = ifelse(replicate == 0,
#                           "Mean RLI", 
#                           "Replicate RLI"),
#            scenario = scenarios[[i]])
#   
# }
# 
# head(scenario_rli_outputs_aggregated[[1]])
# tail(scenario_rli_outputs_aggregated[[1]])
# 
# # Plot all together
# 
# scenario_rli_plots_aggregated <- list()
# 
# for (i in seq_along(scenario_rli_outputs_aggregated)) {
#   
#   scenario_rli_plots_aggregated[[i]] <- ggplot(data = scenario_rli_outputs_aggregated[[i]], 
#                                                aes(x = annual_time_step, y = indicator_score, group = replicate,
#                                                    color = level)) +
#     geom_line() +
#     scale_color_manual(values = c("black", "gray62")) + 
#     labs(x = "Time", 
#          y = "Red List Index Score") +
#     theme(panel.grid.major = element_blank(),
#           axis.title = element_text(size = 18),
#           axis.text = element_text(size = 18),
#           panel.grid.minor = element_blank(),
#           panel.background = element_rect(fill = "grey97"),
#           axis.line = element_line(colour = "black")) +
#     geom_vline(xintercept = impact_start, colour = "red") +
#     geom_vline(xintercept = impact_end, colour = "blue")
#   
# }
# 
# scenario_rli_plots_aggregated[[1]]
# 
# # LIVING PLANET INDEX ----
# 
# # * Create folders ----
# 
# lpi_inputs_folder <- file.path(indicator_inputs_folder, "LPI_inputs", today)
# 
# if( !dir.exists( file.path(lpi_inputs_folder) ) ) {
#   dir.create( file.path(lpi_inputs_folder), recursive = TRUE )
#   
# }
# 
# lpi_outputs_folder <- file.path(indicator_outputs_folder, "LPI_outputs", today)
# 
# if( !dir.exists( file.path(lpi_outputs_folder) ) ) {
#   dir.create( file.path(lpi_outputs_folder), recursive = TRUE )
#   
# }
# 
# lpi_plots_folder <- file.path(indicator_plots_folder, "LPI_plots", today)
# 
# if( !dir.exists( file.path(lpi_plots_folder) ) ) {
#   dir.create( file.path(lpi_plots_folder), recursive = TRUE )
#   
# }
# 
# # TEMP CODE ---
# ## Look at the data we are dealing with
# 
# # data <- scenario_abundance_long[[1]][[1]]
# # 
# # head(data)
# # 
# # ggplot(data, aes(x = time_step, y = abundance,
# #                  col = group_id)) +
# #           geom_line()  + 
# #           geom_text(aes(label= group_id),hjust=0, vjust=0) +
# #           theme(legend.position = "none")
# 
# # * Sample data ----
# 
# scenario_lpi_inputs <- list()
# 
# for (i in seq_along(scenario_smoothed_abundance)) {
#   
#   # Get replicates for a single scenario
#   # replicate_abundance_long <- scenario_smoothed_abundance[[i]]
#   
#   replicate_abundance_long <- scenario_smoothed_abundance[[i]]
#   
#   replicate_lpi_inputs <- list()
#   # Calculate the LPI for each replicate within the scenario
#   for (j in seq_along(replicate_abundance_long)) {
#     
#     replicate_lpi_inputs[[j]] <- replicate_abundance_long[[j]] %>% 
#       dplyr::select(group_id, annual_time_step, 
#                     ave_abundance)
#     
#     
#   }
#   
#   scenario_lpi_inputs[[i]] <- replicate_lpi_inputs
#   
# }
# 
# # lpi_input <- scenario_lpi_inputs[[1]][[2]]
# # head(lpi_input)
# # write.csv(lpi_input, file.path(indicator_outputs_folder, "lpi_input_example_annual.csv"))
# 
# # * Calculate LPI ----
# 
# # Retain naming convention, the LPI just takes the abundance dataframes we
# # already formatted while making the RLI inputs
# 
# # scenario_lpi_inputs <- scenario_abundance_long
# 
# # Loop through each scenario and replicate and calculate the LPI per rep
# 
# scenario_lpi_outputs <- list()
# 
# for (i in seq_along(scenario_lpi_inputs)) {
#   
#   # Get replicates for a single scenario
#   replicate_lpi_inputs <- scenario_lpi_inputs[[i]]
#   
#   replicate_lpi_outputs <- list()
#   
#   # Calculate the LPI for each replicate within the scenario
#   for (j in seq_along(replicate_lpi_inputs)) {
#     
#     replicate_lpi_outputs[[j]] <- calculate_living_planet_index(
#       
#       replicate_lpi_inputs[[j]], start_time_step, ci = FALSE, numboots, j
#     ) 
#     
#     # Save the output LPI data as a csv and rds
#     
#     saveRDS(replicate_lpi_outputs[[j]],
#             file.path(lpi_outputs_folder,
#                       paste(today, scenarios[[i]], "replicate", j,
#                             "LPI_output_data_func_groups.rds",
#                             sep = "_")))
#     
#     write.csv(replicate_lpi_outputs[[j]],
#               file.path(lpi_outputs_folder,
#                         paste(today, scenarios[[i]], "replicate", j,
#                               "LPI_output_data_func_groups.rds",
#                               sep = "_")))
#     
#   }
#   
#   scenario_lpi_outputs[[i]] <- replicate_lpi_outputs
#   
# }
# 
# head(scenario_lpi_outputs)[[1]][[1]]
# 
# # * Aggregate all LPI scores ----
# 
# ## Collapse input data so LPI scores for all replicates in one scenario exist in 
# ## a single data frame
# 
# scenario_lpi_outputs_aggregated <- list()
# 
# for (i in seq_along(scenario_lpi_outputs)) {
#   
#   scenario_lpi_outputs_aggregated[[i]] <- do.call(rbind, 
#                                                   scenario_lpi_outputs[[i]]) %>%
#     mutate(scenario = scenarios[[i]]) 
#   
#   scenario_mean_lpi <- scenario_lpi_outputs_aggregated[[i]] %>%
#     group_by(annual_time_step) %>%
#     summarise(indicator_score = mean(indicator_score),
#               ci_lower = mean(ci_lower),
#               ci_upper = mean(ci_upper)) %>%
#     mutate(replicate = 0,# Replicate 0 will always be the mean
#            indicator = "LPI",
#            scenario = scenarios[[i]]) 
#   
#   scenario_lpi_outputs_aggregated[[i]] <-  rbind(scenario_lpi_outputs_aggregated[[i]],
#                                                  scenario_mean_lpi) %>%
#     mutate(replicate = as.factor(replicate)) %>%
#     mutate(level = ifelse(replicate == 0,
#                           "Mean LPI", 
#                           "Replicate LPI"))
#   
# }
# 
# head(scenario_lpi_outputs_aggregated[[1]])
# tail(scenario_lpi_outputs_aggregated[[1]])
# 
# # * Plot LPI replicates individually ----
# 
# scenario_lpi_plots <- list()
# 
# for (i in seq_along(scenario_lpi_outputs)) {
#   
#   replicate_lpi <- scenario_lpi_outputs[[i]]
#   replicate_lpi_plots <- list()
#   
#   for (j in seq_along(replicate_lpi)) {
#     
#     replicate_lpi_plots[[j]] <- plot_living_planet_index(replicate_lpi[[j]],
#                                                          ci = FALSE)
#     
#     
#     ggsave(file.path(lpi_plots_folder, paste(today, scenarios[[i]], 
#                                              "replicate", j, 
#                                              "LPI_aggregated_func_groups.png",
#                                              sep = "_")),
#            replicate_lpi_plots[[j]],  device = "png")                                   
#     
#   }
#   
#   scenario_lpi_plots[[i]] <- replicate_lpi_plots
#   
# }
# 
# i <- 1
# scenario_lpi_plots[[1]][[i]]
# 
# i <- i + 1
# scenario_lpi_plots[[1]][[i]]
# 
# # * Plot all replicates together ----
# 
# scenario_lpi_plots_aggregated <- list()
# 
# for (i in seq_along(scenario_lpi_outputs_aggregated)){
#   
#   scenario_lpi_plots_aggregated[[i]] <- ggplot(data = scenario_lpi_outputs_aggregated[[i]], 
#                                                aes(x = annual_time_step, 
#                                                    y = indicator_score, 
#                                                    group = replicate,
#                                                    color = level)) +
#     geom_line() +
#     scale_color_manual(values = c("black", "gray62")) + 
#     labs(x = "Time", 
#          y = "Living Planet Index Score") +
#     theme(panel.grid.major = element_blank(),
#           axis.title = element_text(size = 18),
#           axis.text = element_text(size = 18),
#           panel.grid.minor = element_blank(),
#           panel.background = element_rect(fill = "grey97"),
#           axis.line = element_line(colour = "black")) +
#     geom_vline(xintercept = impact_start, colour = "red") +
#     geom_vline(xintercept = impact_end, colour = "blue")
#   
# }
# 
# scenario_lpi_plots_aggregated[[4]]
# 
# # Combine indicators ----
# 
# all_indicators_list <- list(scenario_rli_outputs,
#                             scenario_lpi_outputs)
# 
# names(all_indicators_list) <- c("RLI", "LPI")
# 
# saveRDS(all_indicators_list,
#         file.path(indicator_outputs_folder,
#                   paste(today, "all_indicators_output_data_list_func_groups.rds",
#                         sep = "_")))
# 
# all_lpi <- do.call(rbind, scenario_lpi_outputs_aggregated) %>% 
#   filter(replicate != 0) # Remove the mean so we just have replicates 
# 
# all_rli <- do.call(rbind, scenario_rli_outputs_aggregated) %>% 
#   filter(replicate != 0) # Remove the mean so we just have replicates
# 
# all_indicators <- rbind(all_lpi, all_rli)
# 
# saveRDS(all_indicators,
#         file.path(indicator_outputs_folder,
#                   paste(today, "all_indicators_output_data_func_groups.rds",
#                         sep = "_")))
# 
# write.csv(all_indicators,
#           file.path(indicator_outputs_folder,
#                     paste(today, "all_indicators_output_data_func_groups.csv",
#                           sep = "_")))
# 
# ## UPDATED TO KEEP NA TIMESTEPS AND CHANGED RL CATEGORY THRESHOLDS ----
# 
# # Get generation length ----
# 
# scenario_ab_gl_formatted_not_clean <- list()
# 
# for (i in seq_along(scenario_abundance_long)) {
#   
#   replicate_abundance <- scenario_abundance_long[[i]]
#   replicate_generations <- scenario_generations_raw[[i]]
#   
#   # Make a list to catch the outputs
#   
#   replicate_ab_gl_formatted <- list()
#   
#   # For each individual replicate
#   
#   for (j in seq_along(replicate_abundance)) {
#     
#     # Reduce size of the replicate generations dataframe or the merge won't work
#     gen_length <- replicate_generations[[j]] %>% 
#       dplyr::select(group_id, generation_length_yrs, 
#                     functional_group_name) %>% 
#       distinct(.)
#     
#     # Add the generation length info to the abundance dataframe
#     replicate_ab_gl_formatted[[j]] <- replicate_abundance[[j]] %>%
#       merge(gen_length, by = "group_id") %>%
#       arrange(monthly_time_step, group_id) %>%
#       # Get the timeframe over which to assess decline (3 * gen length or 10 yrs,
#       # whichever is longer)
#       # Important - following lines assume an annual timeframe, will need to adjust if change interval
#       mutate(generation_by_three = generation_length_yrs * 3) %>% # Time over which to measure decline, 3 x gen length OR:
#       mutate(timeframe = generation_by_three) %>%
#       #round the time frame to whole years
#       mutate(timeframe = round(timeframe)) %>% 
#       dplyr::select(-generation_by_three) %>%
#       distinct(.) %>%
#       group_by(group_id) %>% 
#       # select rows that are multiples of the specified interval 
#       # (eg if interval is 12, it samples one month from every 12 (yearly))
#       slice(which(row_number() %% interval == 0)) %>% 
#       mutate(annual_time_step = seq(1,max_timestep,1)) # %>% 
#     
#     
#     print(paste("Replicate", j - 1, 
#                 "formatting complete", 
#                 sep = " "))
#     
#   }
#   
#   print(scenarios[[i]])
#   print(length(replicate_ab_gl_formatted))
#   
#   scenario_ab_gl_formatted_not_clean[[i]] <- replicate_ab_gl_formatted
#   
# }
# 
# # Remove false extinctions ----
# 
# scenario_false_extinctions_removed <- list()
# 
# for (i in seq_along(scenario_ab_gl_formatted_not_clean)) {
#   
#   replicate_ab_gl <- scenario_ab_gl_formatted_not_clean[[i]]
#   
#   # Make a list to catch the outputs
#   
#   replicate_false_ex_removed <- list()
#   
#   # For each individual replicate
#   
#   for (j in seq_along(replicate_ab_gl)) {
#     
#     # Find the last time step where non-0 abundance occurred for each group
#     
#     temp2 <- replicate_ab_gl[[j]] %>% 
#       group_by(group_id) %>% 
#       filter(abundance > 0) %>% 
#       dplyr::select(group_id, annual_time_step, abundance) %>% 
#       filter(annual_time_step == max(annual_time_step)) %>% 
#       dplyr::select(group_id, annual_time_step) %>% 
#       rename(last_abundance = annual_time_step)
#     
#     # Add the year of last positive abundance number as a column to the data    
#     temp3 <- replicate_ab_gl[[j]] %>% 
#       merge(temp2, by = c("group_id"), all = TRUE)
#     
#     # Use the last positive abundance year and current abundance value to determine
#     # if a zero abundance is a true extinction or just a missing value (false extinction)
#     temp4 <- temp3 %>%
#       group_by(group_id) %>%
#       mutate(true_extinction = ifelse(abundance == 0 &
#                                       annual_time_step < last_abundance,
#                                       "false extinction",
#                                ifelse(abundance > 0 &
#                                       annual_time_step < last_abundance,
#                                       "not extinct",
#                                ifelse(abundance == 0 &
#                                       annual_time_step >= last_abundance,
#                                      "true extinction", "not extinct")))) %>%
#       #filter(true_extinction != "false extinction") %>%
#       mutate(abundance = ifelse(true_extinction == "false extinction",
#                                   NA, abundance)) %>%
#       group_by(group_id) %>%
#       arrange(annual_time_step)
#     
#     # Add massbin index
#     data <- temp4 %>% 
#       merge(groups[c("group_id", "bodymass_index", "mass_lower")],
#             by = "group_id") %>% 
#       arrange(functional_group_name,
#               bodymass_index, annual_time_step)
#     
# 
#     # Check if there are any carnivorous endotherms
#     
#     check <- data %>% 
#       group_by(functional_group_name) %>% 
#       summarise(present = sum(abundance, na.rm = TRUE)) %>% 
#       filter(functional_group_name == "carnivore endotherm") %>% 
#       dplyr::select(present) %>% 
#       pull(.)
#     
#     
#     print(paste("Replicate", j - 1, 
#                 "formatting complete", 
#                 sep = " "))
#     
#     # Replace data with 0 if no carnivores
#     
#     if(check == 0) {
#       
#       data <- NULL
#       
#       print(paste("Replicate", j - 1, 
#                   "removed because no carnivorous endotherms are present", 
#                   sep = " "))
#       
#     }
#     
#     replicate_false_ex_removed[[j]] <- data
#     
#   }
#   
#   print(scenario[[i]])
#   print(length(replicate_ab_gl))
#   
#   scenario_false_extinctions_removed[[i]] <- replicate_false_ex_removed
#   
# }
# 
# test <- scenario_false_extinctions_removed[[1]][[1]]
# head(test)
# 
# # Remove replicates with no carnivorous endotherms ----
# 
# scenario_ab_gl_formatted <- list()
# 
# for (i in seq_along(scenario_false_extinctions_removed)) {
#   
#   replicate_not_clean <- scenario_false_extinctions_removed[[i]]
#   
#   scenario_ab_gl_formatted[[i]] <- list.clean(replicate_not_clean)
#   
# }
# 
# test <- scenario_ab_gl_formatted[[1]][[1]]
# head(test)
# 
# # * Smooth abundance ----
# 
# scenario_abundance_clean <- scenario_ab_gl_formatted
# 
# ave_window <- 10
# 
# scenario_smoothed_abundance <- list()
# 
# for (i in seq_along(scenario_abundance_clean)) {
#   
#   # Get replicate data for a single scenario
#   
#   replicate_ab_gl <- scenario_abundance_clean[[i]]
#   
#   replicate_smoothed_abundance <- list()
#   
#   # For each individual replicate
#   
#   for (j in seq_along(replicate_ab_gl)) {
#     
#     group_ab_gl <- replicate_ab_gl[[j]]
#     
#     group_list <- split(group_ab_gl, group_ab_gl$group_id)
#     
#     group_smoothed_abundance <- list()
#     
#     for (k in seq_along(group_list)) {
#       
#       group_df <- group_list[[k]]
#       
#       # check if the group has any abundance values, make it null if not
#       
#       if (sum(group_df$abundance, na.rm = TRUE) == 0) {
# 
#         group_smoothed_abundance[[k]] <- NULL
# 
#       } else {
#         
#         group_smoothed_abundance[[k]] <- group_df %>%
#           arrange(annual_time_step) %>%
#           mutate(ave_abundance = rollapply(abundance,
#                                            ave_window,
#                                            mean,
#                                            na.rm = TRUE,
#                                            partial = TRUE),
#                  ave_abundance = ifelse(ave_abundance < 1,
#                                         0, ave_abundance))
#         
#         print(k)
#         
#       }
#     }
#     
#     all_groups_smooth <- do.call(rbind,group_smoothed_abundance)
#     
#     replicate_smoothed_abundance[[j]] <- all_groups_smooth
#     
#     print(j)
#   }
#   
#   scenario_smoothed_abundance[[i]] <- replicate_smoothed_abundance
#   
#   print(i)
#   
# }
# 
# check <- scenario_smoothed_abundance[[1]][[1]]
# head(check)
# 
# checkgroup <- check %>% filter(group_id == "13.16.27")
# 
# # RED LIST INDEX ----
# 
# # * Create folders ----
# 
# 
# rli_inputs_folder <- file.path(indicator_inputs_folder, "RLI_inputs", today)
# 
# if( !dir.exists( file.path(rli_inputs_folder) ) ) {
#   dir.create( file.path(rli_inputs_folder), recursive = TRUE )
#   
# }
# 
# rli_outputs_folder <- file.path(indicator_outputs_folder, "RLI_outputs", today)
# 
# if( !dir.exists( file.path(rli_outputs_folder) ) ) {
#   dir.create( file.path(rli_outputs_folder), recursive = TRUE )
#   
# }
# 
# rli_plots_folder <- file.path(indicator_plots_folder, "RLI_plots", today)
# 
# if( !dir.exists( file.path(rli_plots_folder) ) ) {
#   dir.create( file.path(rli_plots_folder), recursive = TRUE )
#   
# }
# 
# ## Referring to the thresholds quote under Criterion A, Reason 1 (declines
# ## are the result of reversible pressures) according to:
# ## https://portals.iucn.org/library/sites/library/files/documents/RL-2001-001-2nd.pdf
# 
# 
# # * Assign Red List Categories ----
# 
# scenario_red_list_data <- list()
# 
# #for (i in seq_along(scenario_ab_gl_formatted)) {
# for (i in seq_along(scenario_smoothed_abundance)) {
#   
#   # Get replicate data for a single scenario
#   
#   replicate_ab_gl <- scenario_smoothed_abundance[[i]]
#   
#   print(paste("Processing scenario", scenarios[[i]], sep = " "))
#   
#   replicate_red_list_data <- list()
#   
#   # For each individual replicate
#   
#   for (j in seq_along(replicate_ab_gl)) {
#     
#     print(paste("Processing replicate", j, sep = " "))
#     
#     # Split by functional group, because we calculate RLI for different
#     # functional groups then aggregate later (as per Butchart etal 2010),
#     # except we are using functional groups as proxies for taxa (eg mammals, birds, 
#     # reptiles) used in real world RLI calcs
#     
#     status_inputs <- split(replicate_ab_gl[[j]], 
#                            replicate_ab_gl[[j]]$group_id)
#     
#     # Make a list to hold output for each individual massbin-func-group (ie virtual spp)
#     
#     group_red_list_data <- list()
#     
#     for (k in seq_along(status_inputs)) {
#       
#       print(paste("Processing group", names(status_inputs)[[k]], sep = " "))
#       
#       group_red_list_data[[k]] <- status_inputs[[k]] %>%
#         group_by(group_id) %>%
#         arrange(monthly_time_step) %>%
#         # calculate the difference in abundance over 10 yrs or 3 generation lengths
#         # (specified by 'timeframe' column). Its okay to take the first value of 
#         # timeframe bc the dataframe is grouped by group_id, and timeframe only changes
#         # between and not within group_ids
#         # mutate(diff = (abundance - dplyr::lag(abundance, timeframe[1]))) %>%
#         mutate(diff = (ave_abundance - dplyr::lag(ave_abundance, 10))) %>%
#         # Using the formula from p 35 (Complex patterns of decline) Guidelines 
#         # for Using the IUCN Red List Categories and Criteria v14 August 2019 
#         mutate(decline = 1 - ave_abundance/dplyr::lag(ave_abundance, 10)) %>%
#         mutate(decline = ifelse(ave_abundance == 0, NA, decline)) %>% 
#         # calculate the rate of change
#         # mutate(decline = diff/dplyr::lag(abundance, timeframe[1])) %>% 
#         # mutate(decline = diff/dplyr::lag(ave_abundance, 10)) %>%
#         # mutate(prev = dplyr::lag(ave_abundance, 10)) %>% 
#         # assign red list risk status based on decline 
#         # Using the thresholds from p 16 Categories A2 - A4 Guidelines 
#         # for Using the IUCN Red List Categories and Criteria v14 August 2019
#         mutate(rl_status = ifelse(decline < 0.20, "LC",
#                                   ifelse(decline >= 0.20 & decline < 0.30, "NT", # Where did this and LC thresholds come from?
#                                   ifelse(decline >= 0.30 & decline < 0.50, "VU",
#                                   ifelse(decline >= 0.50 & decline < 0.80, "EN",
#                                   ifelse(decline >= 0.80, "CR",
#                                   ifelse(decline == NA, "EX", "TBD"))))))) %>%
#         arrange(group_id, monthly_time_step) %>%
#         # Replace all non-ex status with ex after first occurrence 
#         # mutate(extinct = match("EX", rl_status)) %>%
#         mutate(extinct = ifelse(rl_status == "EX", 1, 0)) %>% 
#         # mutate(rl_status = with(., ave(rl_status, 
#         #                                         FUN=maintain_ex_status)))
#         #mutate(rl_status = rl_status) %>% 
#         group_by(group_id)
#       
#     }
#     
#     print(paste("replicate", j, "from", scenarios[[i]], "complete", sep = " "))
#     
#     replicate_red_list_df <- do.call(rbind, group_red_list_data)
#     
#     replicate_red_list_data[[j]] <- replicate_red_list_df
#     
#     # Save the inputs
#     
#     saveRDS(replicate_red_list_df,
#             file.path(rli_inputs_folder,
#                       paste(today, scenarios[[i]], "replicate", j,
#                             "RLI_input_data.rds", sep = "_")))
#     
#     write.csv(replicate_red_list_df,
#               file.path(rli_inputs_folder,
#                         paste(today, scenarios[[i]], "replicate", j,
#                               "RLI_input_data.csv", sep = "_")))
#     
#     
#   }
#   
#   scenario_red_list_data[[i]] <- replicate_red_list_data
#   
# }
# 
# # Check we have correct structure still
# length(scenario_red_list_data) == length(scenario_ab_gl_formatted)
# length(scenario_red_list_data[[1]]) == length(scenario_ab_gl_formatted[[1]])
# 
# # Have a quick look at the outputs
# 
# rli_inputs <- scenario_red_list_data[[1]][[1]]
# tail(rli_inputs)
# 
# write.csv(rli_inputs, file.path(indicator_outputs_folder, "rli_input_example_annual.csv"))
# 
# rli_inputs_group <- rli_inputs %>% filter(group_id == "13.16.27")
# 
# #rli_inputs_group <- x %>% filter(group_id == "13.16.27")
# 
# ggplot(data = rli_inputs_group) +
#   geom_path(aes(x = annual_time_step, y = ave_abundance)) +
#   theme(legend.position = "none") +
#   geom_text(aes(x = annual_time_step, y = ave_abundance, label = rl_status))
# 
# # * Get harvested group only ----
# 
# scenario_harvested_groups <- list()
# 
# for (i in seq_along(scenario_red_list_data)) {
#   
#   # Get replicate data for a single scenario
#   
#   replicate_rl_data <- scenario_red_list_data[[i]]
#   
#   harvested <- list()
#   
#   # For each individual replicate
#   
#   for (j in seq_along(replicate_rl_data)) {
#   
#   if  (scenarios[[i]] == "000_Baseline") {
#       
#       harvested[[j]]  <- replicate_rl_data[[j]] %>% 
#         filter(functional_group_name == "carnivore endotherm"|
#                  functional_group_name == "carnivore ectotherm") %>% 
#         filter(mass_lower == 100000)
#   
#       } else if (scenarios[[i]] == "100_Land_use") {
#     
#   harvested[[j]]  <- replicate_rl_data[[j]] %>% 
#          filter(functional_group_name == "herbivore endotherm"|
#                   functional_group_name == "herbivore ectotherm") 
#     
#   } else if (scenarios[[i]] == "200_Harvesting_carnivores") {
#     
#   harvested[[j]] <- replicate_rl_data[[j]] %>% 
#       filter(functional_group_name == "carnivore endotherm"|
#                functional_group_name == "carnivore ectotherm") %>% 
#       filter(mass_lower == 100000)
#     
#   } else if (scenarios[[i]] == "300_Harvesting_herbivores") {
#     
#   harvested[[j]]  <- replicate_rl_data[[j]] %>% 
#       filter(functional_group_name == "herbivore endotherm"|
#              functional_group_name == "herbivore ectotherm") %>% 
#       filter(mass_lower == 100000)
#     
#   } 
#  
# }
#     
#   scenario_harvested_groups[[i]] <- harvested
#     
# }
# 
# harvested_group <- scenario_harvested_groups[[1]][[1]]
# 
# #rli_inputs_group <- x %>% filter(group_id == "13.16.27")
# 
# ggplot(data = harvested_group) +
#   geom_line(aes(x = annual_time_step, y = ave_abundance, col = group_id)) +
#   theme(legend.position = "none") +
#   geom_text(aes(x = annual_time_step, y = ave_abundance, label = rl_status))
# 
# # Plot some results to check they're not completely whack
# 
# ## Get one group to check how their status changes over time relative to how
# ## their abundance changes
# 
# # group_id_select <- "13.16.17" # Shows example of 'resurrected' virtual spp
# # # group_id_select <- "10.40"
# # 
# # data <- rli_inputs %>% dplyr::filter(group_id == group_id_select)
# # 
# # ggplot(data, aes(x = time_step, y = abundance)) +
# #   geom_line() +
# #   geom_text(aes(label= rl_status,
# #                 col = rl_status),hjust=0, vjust=0)
# 
# 
# 
# 
# # * Calculate RLI ----
# 
# # RLI by individual functional groups
# 
# scenario_fg_rli_outputs <- list()
# 
# for (i in seq_along(scenario_red_list_data)) {
#   
#   replicate_red_list_inputs <- scenario_red_list_data[[i]]
#   
#   replicate_fg_rli_outputs <- list()
#   
#   for (j in seq_along(replicate_red_list_inputs)) {
#     
#     replicate_rli <- calculate_red_list_index(
#       replicate_red_list_inputs[[j]], numboots, ci = FALSE) %>%
#       mutate(replicate = j)
#     
#     replicate_fg_rli_outputs[[j]] <- replicate_rli 
#     
#     # saveRDS(replicate_fg_rli_outputs[[j]],
#     #         file.path(rli_outputs_folder,
#     #                   paste(today, scenarios, "replicate", j,
#     #                         "RLI_func_group_output_data.rds",
#     #                         sep = "_")))
#     
#     # write.csv(replicate_fg_rli_outputs[[j]],
#     #           file.path(rli_outputs_folder,
#     #                     paste(today, scenarios[[i]], "RLI_func_group_output_data.rds",
#     #                           sep = "_")))
#     
#     print(paste("RLI for replicate", j, "complete", sep = " "))
#     
#   }
#   
#   scenario_fg_rli_outputs[[i]] <- replicate_fg_rli_outputs
#   
# }
# 
# 
# x <- scenario_fg_rli_outputs[[1]][[3]]
# head(x)
# 
# # Mean RLI aggregated across groups
# 
# scenario_rli_outputs <- list()
# 
# for (i in seq_along(scenario_fg_rli_outputs)) {
#   
#   # Get replicates for a single scenario
#   replicate_rli_fg <- scenario_fg_rli_outputs[[i]]
#   
#   replicate_rli_outputs <- list()
#   
#   # Aggregate RLI across functional groups for each replicate
#   for (j in seq_along(replicate_rli_fg)) {
#     
#     if ("ci_lower" %in% names(replicate_rli_fg[[j]])) {
#       
#       replicate_rli_outputs[[j]] <- replicate_rli_fg[[j]] %>%
#         group_by(annual_time_step) %>%
#         summarise(indicator_score = mean(indicator_score),
#                   ci_lower = mean(ci_lower),
#                   ci_upper = mean(ci_upper)) %>%
#         mutate(indicator = "RLI",
#                replicate = j)
#     } else {
#       
#       replicate_rli_outputs[[j]] <- replicate_rli_fg[[j]] %>%
#         group_by(annual_time_step) %>%
#         summarise(indicator_score = mean(indicator_score)) %>%
#         mutate(indicator = "RLI",
#                replicate = j)
#     }
#     
#     # saveRDS(replicate_rli_outputs[[j]],
#     #       file.path(rli_outputs_folder,
#     #                 paste(today, scenarios[[i]], "replicate", j,
#     #                       "RLI_aggregate_output_data.rds",
#     #                       sep = "_")))
#     # 
#     # write.csv(replicate_rli_outputs[[j]],
#     #           file.path(rli_outputs_folder,
#     #                     paste(today, scenarios[[i]], "replicate", j,
#     #                           "RLI_aggregate_output_data.rds",
#     #                           sep = "_")))
#     
#   }
#   
#   scenario_rli_outputs[[i]] <- replicate_rli_outputs
#   
# }
# 
# head(scenario_rli_outputs)[[1]][[1]]
# 
# # * Plot RLI ----
# 
# ## By functional group
# 
# scenario_fg_rli_plots <- list()
# 
# for (i in seq_along(scenario_fg_rli_outputs)) {
#   
#   replicate_fg_rli <- scenario_fg_rli_outputs[[i]]
#   
#   replicate_fg_rli_plots <- list()
#   
#   for (j in seq_along(replicate_fg_rli)) {
#     
#     replicate_fg_rli_plots[[j]] <-  plot_red_list_index_by_group(
#       replicate_fg_rli[[j]],
#       impact_start,
#       impact_end,
#       ci = FALSE)
#     
#     ggsave(file.path(rli_plots_folder, paste(today, scenarios[[i]], "replicate", j,
#                                              "RLI_by_functional_group_annual_updated.png",
#                                              sep = "_")),
#            replicate_fg_rli_plots[[j]],  device = "png")
#     
#   }
#   
#   scenario_fg_rli_plots[[i]] <- replicate_fg_rli_plots
#   
# }
# 
# scenario_fg_rli_plots[[1]][[2]]
# 
# 
# # Small test to see if averaging indicator scores after works better (it doesn't)
# x <- scenario_rli_outputs[[1]][[5]]
# x <- x[-1,]
# 
# x <- x %>% 
#   mutate(x = rollmean(indicator_score, 10, na.pad = TRUE))
# 
# ggplot(x, aes(x = annual_time_step, y = x))+
#   geom_line()
# 
# # RLI with all functional groups aggregated
# # i.e. mean of each 'taxa' RLI as per Butchart et al (2010) 'Indicators of
# # recent declines'
# 
# scenario_rli_plots <- list()
# 
# for (i in seq_along(scenario_rli_outputs)) {
#   
#   replicate_rli <- scenario_rli_outputs[[i]]
#   
#   replicate_rli_plots <- list()
#   
#   for (j in seq_along(replicate_rli)) {
#     
#     replicate_rli_plots[[j]] <- plot_red_list_index(replicate_rli[[j]],
#                                                     impact_start, 
#                                                     impact_end,
#                                                     ci = TRUE)
#     
#     
#     ggsave(file.path(rli_plots_folder, paste(today, scenarios[[i]], 
#                                              "replicate", j, 
#                                              "RLI_aggregated_annual_updated.png",
#                                              sep = "_")),
#            replicate_rli_plots[[j]],  device = "png")                                   
#     
#   }
#   
#   scenario_rli_plots[[i]] <- replicate_rli_plots
#   
# }
# 
# i <- 1
# i <- i+1
# scenario_rli_plots[[1]][[i]]
# 
# # * Plot all replicates together ----
# 
# ## Collapse input data so RLI scores for all replicates in one scenario exist in 
# ## a single data frame
# 
# scenario_rli_outputs_aggregated <- list()
# 
# for (i in seq_along(scenario_rli_outputs)) {
#   
#   scenario_rli_outputs_aggregated[[i]] <- do.call(rbind, 
#                                                   scenario_rli_outputs[[i]]) %>%
#     mutate(scenario = scenarios[[i]]) 
#   
#   
#   scenario_mean_rli <- scenario_rli_outputs_aggregated[[i]] %>%
#     group_by(annual_time_step) %>%
#     summarise(indicator_score = mean(indicator_score),
#               ci_lower = mean(ci_lower),
#               ci_upper = mean(ci_upper)) %>%
#     mutate(indicator = "RLI",
#            replicate = 0,
#            scenario = scenarios[[i]]) # Replicate 0 will always be the mean
#   
#   scenario_rli_outputs_aggregated[[i]] <-  rbind(scenario_rli_outputs_aggregated[[i]],
#                                                  scenario_mean_rli) %>%
#     mutate(replicate = as.factor(replicate)) %>%
#     mutate(level = ifelse(replicate == 0,
#                           "Mean RLI", 
#                           "Replicate RLI"),
#            scenario = scenarios[[i]])
#   
# }
# 
# head(scenario_rli_outputs_aggregated[[1]])
# tail(scenario_rli_outputs_aggregated[[1]])
# 
# # Plot all together
# 
# scenario_rli_plots_aggregated <- list()
# 
# for (i in seq_along(scenario_rli_outputs_aggregated)) {
#   
#   scenario_rli_plots_aggregated[[i]] <- ggplot(data = scenario_rli_outputs_aggregated[[i]], 
#                                                aes(x = annual_time_step, y = indicator_score, group = replicate,
#                                                    color = level)) +
#     geom_line() +
#     scale_color_manual(values = c("black", "gray62")) + 
#     labs(x = "Time", 
#          y = "Red List Index Score") +
#     theme(panel.grid.major = element_blank(),
#           axis.title = element_text(size = 18),
#           axis.text = element_text(size = 18),
#           panel.grid.minor = element_blank(),
#           panel.background = element_rect(fill = "grey97"),
#           axis.line = element_line(colour = "black")) +
#     geom_vline(xintercept = impact_start, colour = "red") +
#     geom_vline(xintercept = impact_end, colour = "blue")
#   
# }
# 
# scenario_rli_plots_aggregated[[1]]
# 
# # LIVING PLANET INDEX ----
# 
# # * Create folders ----
# 
# lpi_inputs_folder <- file.path(indicator_inputs_folder, "LPI_inputs", today)
# 
# if( !dir.exists( file.path(lpi_inputs_folder) ) ) {
#   dir.create( file.path(lpi_inputs_folder), recursive = TRUE )
#   
# }
# 
# lpi_outputs_folder <- file.path(indicator_outputs_folder, "LPI_outputs", today)
# 
# if( !dir.exists( file.path(lpi_outputs_folder) ) ) {
#   dir.create( file.path(lpi_outputs_folder), recursive = TRUE )
#   
# }
# 
# lpi_plots_folder <- file.path(indicator_plots_folder, "LPI_plots", today)
# 
# if( !dir.exists( file.path(lpi_plots_folder) ) ) {
#   dir.create( file.path(lpi_plots_folder), recursive = TRUE )
#   
# }
# 
# # TEMP CODE ---
# ## Look at the data we are dealing with
# 
# # data <- scenario_abundance_long[[1]][[1]]
# # 
# # head(data)
# # 
# # ggplot(data, aes(x = time_step, y = abundance,
# #                  col = group_id)) +
# #           geom_line()  + 
# #           geom_text(aes(label= group_id),hjust=0, vjust=0) +
# #           theme(legend.position = "none")
# 
# # * Sample data ----
# 
# scenario_lpi_inputs <- list()
# 
# 
# for (i in seq_along(scenario_smoothed_abundance)) {
#   
#   # Get replicates for a single scenario
#   # replicate_abundance_long <- scenario_smoothed_abundance[[i]]
#   
#   replicate_abundance_long <- scenario_smoothed_abundance[[i]]
#   
#   replicate_lpi_inputs <- list()
#   # Calculate the LPI for each replicate within the scenario
#   for (j in seq_along(replicate_abundance_long)) {
#     
#     replicate_lpi_inputs[[j]] <- replicate_abundance_long[[j]] %>% 
#       dplyr::select(group_id, annual_time_step, 
#                     ave_abundance)
#     
#     
#   }
#   
#   scenario_lpi_inputs[[i]] <- replicate_lpi_inputs
#   
# }
# 
# # lpi_input <- scenario_lpi_inputs[[1]][[2]]
# # head(lpi_input)
# # write.csv(lpi_input, file.path(indicator_outputs_folder, "lpi_input_example_annual.csv"))
# 
# # * Calculate LPI ----
# 
# # Retain naming convention, the LPI just takes the abundance dataframes we
# # already formatted while making the RLI inputs
# 
# # scenario_lpi_inputs <- scenario_abundance_long
# 
# # Loop through each scenario and replicate and calculate the LPI per rep
# 
# scenario_lpi_outputs <- list()
# 
# for (i in seq_along(scenario_lpi_inputs)) {
#   
#   # Get replicates for a single scenario
#   replicate_lpi_inputs <- scenario_lpi_inputs[[i]]
#   
#   replicate_lpi_outputs <- list()
#   
#   # Calculate the LPI for each replicate within the scenario
#   for (j in seq_along(replicate_lpi_inputs)) {
#     
#     replicate_lpi_outputs[[j]] <- calculate_living_planet_index(
#       
#       replicate_lpi_inputs[[j]], start_time_step, ci = FALSE, numboots, j
#     ) 
#     
#     # Save the output LPI data as a csv and rds
#     
#     saveRDS(replicate_lpi_outputs[[j]],
#             file.path(lpi_outputs_folder,
#                       paste(today, scenarios[[i]], "replicate", j,
#                             "LPI_output_data_annual_updated.rds",
#                             sep = "_")))
#     
#     write.csv(replicate_lpi_outputs[[j]],
#               file.path(lpi_outputs_folder,
#                         paste(today, scenarios[[i]], "replicate", j,
#                               "LPI_output_data_annual_updated.rds",
#                               sep = "_")))
#     
#   }
#   
#   scenario_lpi_outputs[[i]] <- replicate_lpi_outputs
#   
# }
# 
# head(scenario_lpi_outputs)[[1]][[1]]
# 
# # * Aggregate all LPI scores ----
# 
# ## Collapse input data so LPI scores for all replicates in one scenario exist in 
# ## a single data frame
# 
# scenario_lpi_outputs_aggregated <- list()
# 
# for (i in seq_along(scenario_lpi_outputs)) {
#   
#   scenario_lpi_outputs_aggregated[[i]] <- do.call(rbind, 
#                                                   scenario_lpi_outputs[[i]]) %>%
#     mutate(scenario = scenarios[[i]]) 
#   
#   scenario_mean_lpi <- scenario_lpi_outputs_aggregated[[i]] %>%
#     group_by(annual_time_step) %>%
#     summarise(indicator_score = mean(indicator_score),
#               ci_lower = mean(ci_lower),
#               ci_upper = mean(ci_upper)) %>%
#     mutate(replicate = 0,# Replicate 0 will always be the mean
#            indicator = "LPI",
#            scenario = scenarios[[i]]) 
#   
#   scenario_lpi_outputs_aggregated[[i]] <-  rbind(scenario_lpi_outputs_aggregated[[i]],
#                                                  scenario_mean_lpi) %>%
#     mutate(replicate = as.factor(replicate)) %>%
#     mutate(level = ifelse(replicate == 0,
#                           "Mean LPI", 
#                           "Replicate LPI"))
#   
# }
# 
# head(scenario_lpi_outputs_aggregated[[1]])
# tail(scenario_lpi_outputs_aggregated[[1]])
# 
# # * Plot LPI replicates individually ----
# 
# scenario_lpi_plots <- list()
# 
# for (i in seq_along(scenario_lpi_outputs)) {
#   
#   replicate_lpi <- scenario_lpi_outputs[[i]]
#   replicate_lpi_plots <- list()
#   
#   for (j in seq_along(replicate_lpi)) {
#     
#     replicate_lpi_plots[[j]] <- plot_living_planet_index(replicate_lpi[[j]],
#                                                          ci = FALSE)
#     
#     
#     ggsave(file.path(lpi_plots_folder, paste(today, scenarios[[i]], 
#                                              "replicate", j, 
#                                              "LPI_aggregated_annual_updated.png",
#                                              sep = "_")),
#            replicate_lpi_plots[[j]],  device = "png")                                   
#     
#   }
#   
#   scenario_lpi_plots[[i]] <- replicate_lpi_plots
#   
# }
# 
# i <- 1
# scenario_lpi_plots[[1]][[i]]
# 
# i <- i + 1
# scenario_lpi_plots[[1]][[i]]
# 
# # * Plot all replicates together ----
# 
# scenario_lpi_plots_aggregated <- list()
# 
# for (i in seq_along(scenario_lpi_outputs_aggregated)){
#   
#   scenario_lpi_plots_aggregated[[i]] <- ggplot(data = scenario_lpi_outputs_aggregated[[i]], 
#                                                aes(x = annual_time_step, 
#                                                    y = indicator_score, 
#                                                    group = replicate,
#                                                    color = level)) +
#     geom_line() +
#     scale_color_manual(values = c("black", "gray62")) + 
#     labs(x = "Time", 
#          y = "Living Planet Index Score") +
#     theme(panel.grid.major = element_blank(),
#           axis.title = element_text(size = 18),
#           axis.text = element_text(size = 18),
#           panel.grid.minor = element_blank(),
#           panel.background = element_rect(fill = "grey97"),
#           axis.line = element_line(colour = "black")) +
#     geom_vline(xintercept = impact_start, colour = "red") +
#     geom_vline(xintercept = impact_end, colour = "blue")
#   
# }
# 
# scenario_lpi_plots_aggregated[[1]]
# 
# # Combine indicators ----
# 
# all_indicators_list <- list(scenario_rli_outputs,
#                             scenario_lpi_outputs)
# 
# names(all_indicators_list) <- c("RLI", "LPI")
# 
# saveRDS(all_indicators_list,
#         file.path(indicator_outputs_folder,
#                   paste(today, "all_indicators_output_data_list_annual_updated.rds",
#                         sep = "_")))
# 
# all_lpi <- do.call(rbind, scenario_lpi_outputs_aggregated) %>% 
#   filter(replicate != 0) # Remove the mean so we just have replicates 
# 
# all_rli <- do.call(rbind, scenario_rli_outputs_aggregated) %>% 
#   filter(replicate != 0) # Remove the mean so we just have replicates
# 
# all_indicators <- rbind(all_lpi, all_rli)
# 
# saveRDS(all_indicators,
#         file.path(indicator_outputs_folder,
#                   paste(today, "all_indicators_output_data_annual_updated.rds",
#                         sep = "_")))
# 
# write.csv(all_indicators,
#           file.path(indicator_outputs_folder,
#                     paste(today, "all_indicators_output_data_annual_updated.csv",
#                           sep = "_")))
# 
# ## TAKE MONTHLY NOT ANNUAL MEAN ----
# 
# # Get generation length ----
# 
# scenario_ab_gl_formatted_not_clean <- list()
# 
# for (i in seq_along(scenario_abundance_long)) {
#   
#   replicate_abundance <- scenario_abundance_long[[i]]
#   replicate_generations <- scenario_generations_raw[[i]]
#   
#   # Make a list to catch the outputs
#   
#   replicate_ab_gl_formatted <- list()
#   
#   # For each individual replicate
#   
#   for (j in seq_along(replicate_abundance)) {
#     
#     # Reduce size of the replicate generations dataframe or the merge won't work
#     gen_length <- replicate_generations[[j]] %>% 
#       dplyr::select(group_id, generation_length_yrs, 
#                     functional_group_name) %>% 
#       distinct(.)
#     
#     # Add the generation length info to the abundance dataframe
#     replicate_ab_gl_formatted[[j]] <- replicate_abundance[[j]] %>%
#       merge(gen_length, by = "group_id") %>%
#       arrange(monthly_time_step, group_id) %>%
#       # Get the timeframe over which to assess decline (3 * gen length or 10 yrs,
#       # whichever is longer)
#       # Important - following lines assume an annual timeframe, will need to adjust if change interval
#       mutate(generation_by_three = generation_length_yrs * 3) %>% # Time over which to measure decline, 3 x gen length OR:
#       mutate(timeframe = generation_by_three) %>%
#       #round the time frame to whole years
#       mutate(timeframe = round(timeframe)) %>% 
#       dplyr::select(-generation_by_three) %>%
#       distinct(.) %>%
#       group_by(group_id) %>% 
#       # Add info about bodymass so can easily sort
#       merge(groups[c("group_id", "bodymass_index", "mass_lower")], by = "group_id") %>% 
#       arrange(monthly_time_step)
#       # select rows that are multiples of the specified interval 
#       # (eg if interval is 12, it samples one month from every 12 (yearly))
#       # slice(which(row_number() %% interval == 0)) %>% 
#       # mutate(annual_time_step = seq(1,max_timestep,1)) # %>% 
#     
#     
#     print(paste("Replicate", j - 1, 
#                 "formatting complete", 
#                 sep = " "))
#     
#   }
#   
#   print(scenario[[i]])
#   print(length(replicate_ab_gl_formatted))
#   
#   scenario_ab_gl_formatted_not_clean[[i]] <- replicate_ab_gl_formatted
#   
# }
# 
# test <- scenario_ab_gl_formatted_not_clean[[1]][[1]]
# head(test)
# # Remove false extinctions ----
# 
# scenario_false_extinctions_removed <- list()
# 
# for (i in seq_along(scenario_ab_gl_formatted_not_clean)) {
#   
#   replicate_ab_gl <- scenario_ab_gl_formatted_not_clean[[i]]
#   
#   # Make a list to catch the outputs
#   
#   replicate_false_ex_removed <- list()
#   
#   # For each individual replicate
#   
#   for (j in seq_along(replicate_ab_gl)) {
#     
#     # Find the last time step where non-0 abundance occurred for each group
#     
#     temp2 <- replicate_ab_gl[[j]] %>% 
#       group_by(group_id) %>% 
#       filter(abundance > 0) %>% 
#       dplyr::select(group_id, monthly_time_step, abundance) %>% 
#       filter(monthly_time_step == max(monthly_time_step)) %>% 
#       dplyr::select(group_id, monthly_time_step) %>% 
#       rename(last_abundance = monthly_time_step)
#     
#     # Add the year of last positive abundance number as a column to the data    
#     temp3 <- replicate_ab_gl[[j]] %>% 
#       merge(temp2, by = c("group_id"), all = TRUE)
#     
#     # Use the last positive abundance year and current abundance value to determine
#     # if a zero abundance is a true extinction or just a missing value (false extinction)
#     data <- temp3 %>%
#       group_by(group_id) %>%
#       mutate(true_extinction = ifelse(abundance == 0 &
#                                         monthly_time_step < last_abundance,
#                                       "false extinction",
#                                       ifelse(abundance > 0 &
#                                                monthly_time_step < last_abundance,
#                                              "not extinct",
#                                              ifelse(abundance == 0 &
#                                                       monthly_time_step >= last_abundance,
#                                                     "true extinction", "not extinct")))) %>%
#       #filter(true_extinction != "false extinction") %>%
#       mutate(abundance = ifelse(true_extinction == "false extinction",
#                                 NA, abundance)) %>%
#       group_by(group_id) %>%
#       arrange(group_id, monthly_time_step)
#     
#     # Check if there are any carnivorous endotherms
#     
#     check <- data %>% 
#       group_by(functional_group_name) %>% 
#       summarise(present = sum(abundance, na.rm = TRUE)) %>% 
#       filter(functional_group_name == "carnivore endotherm") %>% 
#       dplyr::select(present) %>% 
#       pull(.)
#     
#     
#     print(paste("Replicate", j - 1, 
#                 "formatting complete", 
#                 sep = " "))
#     
#     # Replace data with 0 if no carnivores
#     
#     if(check == 0) {
#       
#       data <- NULL
#       
#       print(paste("Replicate", j - 1, 
#                   "removed because no carnivorous endotherms are present", 
#                   sep = " "))
#       
#     }
#     
#     replicate_false_ex_removed[[j]] <- data
#     
#   }
#   
#   print(scenario[[i]])
#   print(length(replicate_ab_gl))
#   
#   scenario_false_extinctions_removed[[i]] <- replicate_false_ex_removed
#   
# }
# 
# test <- scenario_false_extinctions_removed[[1]][[1]]
# head(test)
# 
# # Remove replicates with no carnivorous endotherms ----
# 
# scenario_ab_gl_formatted <- list()
# 
# for (i in seq_along(scenario_false_extinctions_removed)) {
#   
#   replicate_not_clean <- scenario_false_extinctions_removed[[i]]
#   
#   scenario_ab_gl_formatted[[i]] <- list.clean(replicate_not_clean)
#   
# }
# 
# test <- scenario_ab_gl_formatted[[1]][[1]]
# head(test)
# 
# # Smooth abundance ----
# 
# scenario_abundance_clean <- scenario_ab_gl_formatted
# 
# ave_window <- 120
# 
# scenario_smoothed_abundance <- list()
# 
# for (i in seq_along(scenario_abundance_clean)) {
#   
#   # Get replicate data for a single scenario
#   
#   replicate_ab_gl <- scenario_abundance_clean[[i]]
#   
#   replicate_smoothed_abundance <- list()
#   
#   # For each individual replicate
#   
#   for (j in seq_along(replicate_ab_gl)) {
#     
#     group_ab_gl <- replicate_ab_gl[[j]]
#     
#     group_list <- split(group_ab_gl, group_ab_gl$group_id)
#     
#     group_smoothed_abundance <- list()
#     
#     for (k in seq_along(group_list)) {
#       
#       group_df <- group_list[[k]]
#       
#       # check if the group has any abundance values, make it null if not
#       
#       if (sum(group_df$abundance, na.rm = TRUE) == 0) {
#         
#         group_smoothed_abundance[[k]] <- NULL
#         
#       } else {
#         
#         group_smoothed_abundance[[k]] <- group_df %>%
#           arrange(monthly_time_step) %>%
#           mutate(ave_abundance = rollapply(abundance,
#                                            ave_window,
#                                            mean,
#                                            na.rm = TRUE,
#                                            partial = TRUE),
#                  ave_abundance = ifelse(ave_abundance < 1,
#                                         0, ave_abundance))
#         
#         print(k)
#         
#       }
#     }
#     
#     all_groups_smooth <- do.call(rbind,group_smoothed_abundance)
#     
#     replicate_smoothed_abundance[[j]] <- all_groups_smooth
#     
#     print(j)
#   }
#   
#   scenario_smoothed_abundance[[i]] <- replicate_smoothed_abundance
#   
#   print(i)
#   
# }
# 
# check <- scenario_smoothed_abundance[[1]][[1]]
# head(check)
# 
# checkgroup <- check %>% filter(group_id == "13.16.27")
# 
# # RED LIST INDEX ----
# 
# # * Create folders ----
# 
# 
# rli_inputs_folder <- file.path(indicator_inputs_folder, "RLI_inputs", today)
# 
# if( !dir.exists( file.path(rli_inputs_folder) ) ) {
#   dir.create( file.path(rli_inputs_folder), recursive = TRUE )
#   
# }
# 
# rli_outputs_folder <- file.path(indicator_outputs_folder, "RLI_outputs", today)
# 
# if( !dir.exists( file.path(rli_outputs_folder) ) ) {
#   dir.create( file.path(rli_outputs_folder), recursive = TRUE )
#   
# }
# 
# rli_plots_folder <- file.path(indicator_plots_folder, "RLI_plots", today)
# 
# if( !dir.exists( file.path(rli_plots_folder) ) ) {
#   dir.create( file.path(rli_plots_folder), recursive = TRUE )
#   
# }
# 
# # * Take an annual sample ----
# 
# scenario_smoothed_abundance_annual <- list()
# 
# for (i in seq_along(scenario_smoothed_abundance)) {
#   
#   # Get replicate data for a single scenario
#   
#   replicate <- scenario_smoothed_abundance[[i]]
#   
#   replicate_smoothed_annual <- list()
#   
#   # For each individual replicate
#   
#   for (j in seq_along(replicate)) {
#     
#     replicate_smoothed_annual[[j]] <- replicate[[j]] %>% 
#       slice(which(row_number() %% interval == 0)) %>% 
#       mutate(annual_time_step = seq(1,max_timestep,1))
#     
# 
#   }
#   
#   scenario_smoothed_abundance_annual[[i]] <- replicate_smoothed_annual
#   
# }
# 
# check <- scenario_smoothed_abundance_annual[[1]][[1]]
# dim(check)
# 
# ## Referring to the thresholds quote under Criterion A, Reason 1 (declines
# ## are the result of reversible pressures) according to:
# ## https://portals.iucn.org/library/sites/library/files/documents/RL-2001-001-2nd.pdf
# 
# 
# # * Assign Red List Categories ----
# 
# scenario_red_list_data <- list()
# 
# #for (i in seq_along(scenario_ab_gl_formatted)) {
# for (i in seq_along(scenario_smoothed_abundance_annual)) {
#   
#   # Get replicate data for a single scenario
#   
#   replicate_ab_gl <- scenario_smoothed_abundance_annual[[i]]
#   
#   print(paste("Processing scenario", scenarios[[i]], sep = " "))
#   
#   replicate_red_list_data <- list()
#   
#   # For each individual replicate
#   
#   for (j in seq_along(replicate_ab_gl)) {
#     
#     print(paste("Processing replicate", j, sep = " "))
#     
#     # Split by functional group, because we calculate RLI for different
#     # functional groups then aggregate later (as per Butchart etal 2010),
#     # except we are using functional groups as proxies for taxa (eg mammals, birds, 
#     # reptiles) used in real world RLI calcs
#     
#     status_inputs <- split(replicate_ab_gl[[j]], 
#                            replicate_ab_gl[[j]]$group_id)
#     
#     # Make a list to hold output for each individual massbin-func-group (ie virtual spp)
#     
#     group_red_list_data <- list()
#     
#     for (k in seq_along(status_inputs)) {
#       
#       print(paste("Processing group", names(status_inputs)[[k]], sep = " "))
#       
#       group_red_list_data[[k]] <- status_inputs[[k]] %>%
#         group_by(group_id) %>%
#         arrange(monthly_time_step) %>%
#         # calculate the difference in abundance over 10 yrs or 3 generation lengths
#         # (specified by 'timeframe' column). Its okay to take the first value of 
#         # timeframe bc the dataframe is grouped by group_id, and timeframe only changes
#         # between and not within group_ids
#         # mutate(diff = (abundance - dplyr::lag(abundance, timeframe[1]))) %>%
#         mutate(diff = (ave_abundance - dplyr::lag(ave_abundance, 10))) %>%
#         # Using the formula from p 35 (Complex patterns of decline) Guidelines 
#         # for Using the IUCN Red List Categories and Criteria v14 August 2019 
#         mutate(decline = 1 - ave_abundance/dplyr::lag(ave_abundance, 10)) %>%
#         mutate(decline = ifelse(ave_abundance == 0, NA, decline)) %>% 
#         # calculate the rate of change
#         # mutate(decline = diff/dplyr::lag(abundance, timeframe[1])) %>% 
#         # mutate(decline = diff/dplyr::lag(ave_abundance, 10)) %>%
#         # mutate(prev = dplyr::lag(ave_abundance, 10)) %>% 
#         # assign red list risk status based on decline 
#         # Using the thresholds from p 16 Categories A2 - A4 Guidelines 
#         # for Using the IUCN Red List Categories and Criteria v14 August 2019
#         mutate(rl_status = ifelse(decline < 0.20, "LC",
#                                   ifelse(decline >= 0.20 & decline < 0.30, "NT", # Where did this and LC thresholds come from?
#                                   ifelse(decline >= 0.30 & decline < 0.50, "VU",
#                                   ifelse(decline >= 0.50 & decline < 0.80, "EN",
#                                   ifelse(decline >= 0.80, "CR",
#                                   ifelse(decline == NA, "EX", "TBD"))))))) %>%
#         arrange(group_id, monthly_time_step) %>%
#         # Replace all non-ex status with ex after first occurrence 
#         # mutate(extinct = match("EX", rl_status)) %>%
#         mutate(extinct = ifelse(rl_status == "EX", 1, 0)) %>% 
#         # mutate(rl_status = with(., ave(rl_status, 
#         #                                         FUN=maintain_ex_status)))
#         #mutate(rl_status = rl_status) %>% 
#         group_by(group_id)
#       
#     }
#     
#     print(paste("replicate", j, "from", scenarios[[i]], "complete", sep = " "))
#     
#     replicate_red_list_df <- do.call(rbind, group_red_list_data)
#     
#     replicate_red_list_data[[j]] <- replicate_red_list_df
#     
#     # Save the inputs
#     
#     saveRDS(replicate_red_list_df,
#             file.path(rli_inputs_folder,
#                       paste(today, scenarios[[i]], "replicate", j,
#                             "RLI_input_data_monthly_smoothing.rds", sep = "_")))
#     
#     write.csv(replicate_red_list_df,
#               file.path(rli_inputs_folder,
#                         paste(today, scenarios[[i]], "replicate", j,
#                               "RLI_input_data_monthly_smoothing.csv", sep = "_")))
#     
#     
#   }
#   
#   scenario_red_list_data[[i]] <- replicate_red_list_data
#   
# }
# 
# # Check we have correct structure still
# length(scenario_red_list_data) == length(scenario_ab_gl_formatted)
# length(scenario_red_list_data[[1]]) == length(scenario_ab_gl_formatted[[1]])
# 
# # Have a quick look at the outputs
# 
# rli_inputs <- scenario_red_list_data[[1]][[1]]
# tail(rli_inputs)
# 
# write.csv(rli_inputs, file.path(indicator_outputs_folder, "rli_input_example_monthly_smoothing.csv"))
# 
# rli_inputs_group <- rli_inputs %>% filter(group_id == "13.16.27")
# 
# #rli_inputs_group <- x %>% filter(group_id == "13.16.27")
# 
# ggplot(data = rli_inputs_group) +
#   geom_path(aes(x = annual_time_step, y = ave_abundance)) +
#   theme(legend.position = "none") +
#   geom_text(aes(x = annual_time_step, y = ave_abundance, label = rl_status))
# 
# 
# # * Take coarser sample ----
# 
# sample_interval <- 1 # make a different number than one to actually sample
# sample_max_timestep <- 300/sample_interval
# 
# scenario_redlist_data_sampled <- list()
# 
# for (i in seq_along(scenario_red_list_data)) {
#   
#   # Get replicate data for a single scenario
#   
#   replicate <- scenario_red_list_data[[i]]
#   
#   replicate_sampled <- list()
#   
#   # For each individual replicate
#   
#   for (j in seq_along(replicate)) {
#     
#     replicate_sampled[[j]] <- replicate[[j]] %>% 
#       slice(which(row_number() %% sample_interval == 0)) %>% 
#       mutate(annual_time_step = seq(1,sample_max_timestep,1))
#     
#     
#   }
#   
#   scenario_redlist_data_sampled[[i]] <- replicate_sampled
#   
# }
# 
# test <- scenario_redlist_data_sampled[[1]][[1]]
# test_group <- test %>% filter(group_id == "13.16.27")
# dim(test_group)
# 
# # * Get harvested group only ----
# 
# scenario_harvested_groups <- list()
# 
# for (i in seq_along(scenario_redlist_data_sampled)) {
#   
#   # Get replicate data for a single scenario
#   
#   replicate_rl_data <- scenario_redlist_data_sampled[[i]]
#   
#   harvested <- list()
#   
#   # For each individual replicate
#   
#   for (j in seq_along(replicate_rl_data)) {
#     
#     
#     if (scenarios[[i]] == "100_Land_Use") {
#       
#       harvested[[j]]  <- replicate_rl_data[[j]] %>% 
#         filter(functional_group_name == "herbivore endotherm"|
#                  functional_group_name == "herbivore ectotherm") 
#       
#     } else if (scenarios[[i]] == "200_Harvesting_carnivores") {
#       
#       harvested[[j]] <- replicate_rl_data[[j]] %>% 
#         filter(functional_group_name == "carnivore endotherm" & mass_lower == 10000|
#                  functional_group_name == "carnivore ectotherm" & mass_lower == 10000) 
#       
#     } else if (scenarios[[i]] == "300_Harvesting_herbivores") {
#       
#       harvested[[j]]  <- replicate_rl_data[[j]] %>% 
#         filter(functional_group_name == "herbivore endotherm" & mass_lower == 10000|
#                  functional_group_name == "herbivore ectotherm" & mass_lower == 10000) 
#       
#     } else if (scenarios[[i]] == "000_Baseline") {
#       
#       harvested[[j]]  <- replicate_rl_data[[j]] %>% 
#         filter(functional_group_name == "carnivore endotherm"|
#                  functional_group_name == "carnivore ectotherm") %>% 
#         filter(mass_lower == 1000)
#       
#     }
#     
#   }
#   
#   scenario_harvested_groups[[i]] <- harvested
#   
# }
# 
# harvested_rep <- scenario_harvested_groups[[1]][[1]]
# 
# ggplot(data = harvested_rep) +
#   geom_line(aes(x = annual_time_step, y = ave_abundance, col = group_id)) +
#   theme(legend.position = "none") +
#   geom_text(aes(x = annual_time_step, y = ave_abundance, label = rl_status))
# 
# # * Plot harvested groups ----
# 
# harvested_plots_folder <- file.path(indicator_plots_folder, "harvested_plots", today)
# 
# if( !dir.exists( file.path(harvested_plots_folder) ) ) {
#   dir.create( file.path(harvested_plots_folder), recursive = TRUE )
#   
# }
# 
# 
# scenario_harvested_plots <- list()
# 
# for ( i in seq_along(scenario_harvested_groups)) {
#   
#   replicate_harvest <- scenario_harvested_groups[[i]]
#   
#   replicate_harvested_plots <- list()
#   
#   # For each individual replicate
#   
#   for (j in seq_along(replicate_harvest)) {
#     
#     replicate_harvested_plots[[j]] <- ggplot(data = replicate_harvest[[j]]) +
#       geom_smooth(aes(x = annual_time_step, 
#                       y = abundance, 
#                       col = functional_group_name)) +
#       theme(legend.position = "bottom") +
#       labs(title = paste(scenarios[[i]], "harvested groups", sep = " "))
#     
#     ggsave(file.path(harvested_plots_folder, paste(today, scenarios[[i]], 
#                                                    "replicate", j - 1,
#                                                    "harvested_reps_averaged.png",
#                                                    sep = "_")),
#            replicate_harvested_plots[[j]],  device = "png")
#   }
#   
#   scenario_harvested_plots[[i]] <- replicate_harvested_plots 
#   
# }
# 
# # Plot some results to check they're not completely whack
# 
# ## Get one group to check how their status changes over time relative to how
# ## their abundance changes
# 
# # group_id_select <- "13.16.17" # Shows example of 'resurrected' virtual spp
# # # group_id_select <- "10.40"
# # 
# # data <- rli_inputs %>% dplyr::filter(group_id == group_id_select)
# # 
# # ggplot(data, aes(x = time_step, y = abundance)) +
# #   geom_line() +
# #   geom_text(aes(label= rl_status,
# #                 col = rl_status),hjust=0, vjust=0)
# 
# 
# 
# 
# # * Calculate RLI ----
# 
# # RLI by individual functional groups
# 
# scenario_fg_rli_outputs <- list()
# 
# for (i in seq_along(scenario_redlist_data_sampled)) {
#   
#   replicate_red_list_inputs <- scenario_redlist_data_sampled[[i]]
#   
#   replicate_fg_rli_outputs <- list()
#   
#   for (j in seq_along(replicate_red_list_inputs)) {
#     
#     replicate_rli <- calculate_red_list_index(
#       replicate_red_list_inputs[[j]], numboots, ci = FALSE) %>%
#       mutate(replicate = j)
#     
#     replicate_fg_rli_outputs[[j]] <- replicate_rli 
#     
#     # saveRDS(replicate_fg_rli_outputs[[j]],
#     #         file.path(rli_outputs_folder,
#     #                   paste(today, scenarios, "replicate", j,
#     #                         "RLI_func_group_output_data.rds",
#     #                         sep = "_")))
#     
#     write.csv(replicate_fg_rli_outputs[[j]],
#               file.path(rli_outputs_folder,
#                         paste(today, scenarios[[i]], "RLI_func_group_output_data_monthly_smoothing.rds",
#                               sep = "_")))
#     
#     print(paste("RLI for replicate", j, "complete", sep = " "))
#     
#   }
#   
#   scenario_fg_rli_outputs[[i]] <- replicate_fg_rli_outputs
#   
# }
# 
# 
# x <- scenario_fg_rli_outputs[[1]][[3]]
# head(x)
# 
# # Mean RLI aggregated across groups
# 
# scenario_rli_outputs <- list()
# 
# for (i in seq_along(scenario_fg_rli_outputs)) {
#   
#   # Get replicates for a single scenario
#   replicate_rli_fg <- scenario_fg_rli_outputs[[i]]
#   
#   replicate_rli_outputs <- list()
#   
#   # Aggregate RLI across functional groups for each replicate
#   for (j in seq_along(replicate_rli_fg)) {
#     
#     if ("ci_lower" %in% names(replicate_rli_fg[[j]])) {
#       
#       replicate_rli_outputs[[j]] <- replicate_rli_fg[[j]] %>%
#         group_by(annual_time_step) %>%
#         summarise(indicator_score = mean(indicator_score),
#                   ci_lower = mean(ci_lower),
#                   ci_upper = mean(ci_upper)) %>%
#         mutate(indicator = "RLI",
#                replicate = j)
#     } else {
#       
#       replicate_rli_outputs[[j]] <- replicate_rli_fg[[j]] %>%
#         group_by(annual_time_step) %>%
#         summarise(indicator_score = mean(indicator_score)) %>%
#         mutate(indicator = "RLI",
#                replicate = j)
#     }
#     
#     # saveRDS(replicate_rli_outputs[[j]],
#     #       file.path(rli_outputs_folder,
#     #                 paste(today, scenarios[[i]], "replicate", j,
#     #                       "RLI_aggregate_output_data.rds",
#     #                       sep = "_")))
#     # 
#     # write.csv(replicate_rli_outputs[[j]],
#     #           file.path(rli_outputs_folder,
#     #                     paste(today, scenarios[[i]], "replicate", j,
#     #                           "RLI_aggregate_output_data.rds",
#     #                           sep = "_")))
#     
#   }
#   
#   scenario_rli_outputs[[i]] <- replicate_rli_outputs
#   
# }
# 
# head(scenario_rli_outputs)[[1]][[1]]
# 
# # * Plot RLI ----
# 
# ## By functional group
# 
# scenario_fg_rli_plots <- list()
# 
# for (i in seq_along(scenario_fg_rli_outputs)) {
#   
#   replicate_fg_rli <- scenario_fg_rli_outputs[[i]]
#   
#   replicate_fg_rli_plots <- list()
#   
#   for (j in seq_along(replicate_fg_rli)) {
#     
#     replicate_fg_rli_plots[[j]] <-  plot_red_list_index_by_group(
#       replicate_fg_rli[[j]],
#       impact_start,
#       impact_end,
#       ci = FALSE)
#     
#     ggsave(file.path(rli_plots_folder, paste(today, scenarios[[i]], "replicate", j,
#                                              "RLI_by_functional_group_monthly_smoothing.png",
#                                              sep = "_")),
#            replicate_fg_rli_plots[[j]],  device = "png")
#     
#   }
#   
#   scenario_fg_rli_plots[[i]] <- replicate_fg_rli_plots
#   
# }
# 
# scenario_fg_rli_plots[[1]][[5]]
# 
# 
# # Small test to see if averaging indicator scores after works better (it doesn't)
# x <- scenario_rli_outputs[[1]][[5]]
# x <- x[-1,]
# 
# x <- x %>% 
#   mutate(x = rollmean(indicator_score, 10, na.pad = TRUE))
# 
# ggplot(x, aes(x = annual_time_step, y = x))+
#   geom_line()
# 
# # RLI with all functional groups aggregated
# # i.e. mean of each 'taxa' RLI as per Butchart et al (2010) 'Indicators of
# # recent declines'
# 
# scenario_rli_plots <- list()
# 
# for (i in seq_along(scenario_rli_outputs)) {
#   
#   replicate_rli <- scenario_rli_outputs[[i]]
#   
#   replicate_rli_plots <- list()
#   
#   for (j in seq_along(replicate_rli)) {
#     
#     replicate_rli_plots[[j]] <- plot_red_list_index(replicate_rli[[j]],
#                                                     impact_start, 
#                                                     impact_end,
#                                                     ci = TRUE)
#     
#     
#     ggsave(file.path(rli_plots_folder, paste(today, scenarios[[i]], 
#                                              "replicate", j, 
#                                              "RLI_aggregated_monthly_smoothing.png",
#                                              sep = "_")),
#            replicate_rli_plots[[j]],  device = "png")                                   
#     
#   }
#   
#   scenario_rli_plots[[i]] <- replicate_rli_plots
#   
# }
# 
# i <- 1
# i <- i+1
# scenario_rli_plots[[1]][[i]]
# 
# # * Plot all replicates together ----
# 
# ## Collapse input data so RLI scores for all replicates in one scenario exist in 
# ## a single data frame
# 
# scenario_rli_outputs_aggregated <- list()
# 
# for (i in seq_along(scenario_rli_outputs)) {
#   
#   scenario_rli_outputs_aggregated[[i]] <- do.call(rbind, 
#                                                   scenario_rli_outputs[[i]]) %>%
#     mutate(scenario = scenarios[[i]]) 
#   
#   
#   scenario_mean_rli <- scenario_rli_outputs_aggregated[[i]] %>%
#     group_by(annual_time_step) %>%
#     summarise(indicator_score = mean(indicator_score),
#               ci_lower = mean(ci_lower),
#               ci_upper = mean(ci_upper)) %>%
#     mutate(indicator = "RLI",
#            replicate = 0,
#            scenario = scenarios[[i]]) # Replicate 0 will always be the mean
#   
#   scenario_rli_outputs_aggregated[[i]] <-  rbind(scenario_rli_outputs_aggregated[[i]],
#                                                  scenario_mean_rli) %>%
#     mutate(replicate = as.factor(replicate)) %>%
#     mutate(level = ifelse(replicate == 0,
#                           "Mean RLI", 
#                           "Replicate RLI"),
#            scenario = scenarios[[i]])
#   
# }
# 
# head(scenario_rli_outputs_aggregated[[1]])
# tail(scenario_rli_outputs_aggregated[[1]])
# 
# # Plot all together
# 
# scenario_rli_plots_aggregated <- list()
# 
# for (i in seq_along(scenario_rli_outputs_aggregated)) {
#   
#   scenario_rli_plots_aggregated[[i]] <- ggplot(data = scenario_rli_outputs_aggregated[[i]], 
#                                                aes(x = annual_time_step, 
#                                                    y = indicator_score, 
#                                                    group = replicate,
#                                                    color = level)) +
#     geom_line() +
#     scale_color_manual(values = c("black", "gray62")) + 
#     labs(x = "Time", 
#          y = "Red List Index Score") +
#     theme(panel.grid.major = element_blank(),
#           axis.title = element_text(size = 18),
#           axis.text = element_text(size = 18),
#           panel.grid.minor = element_blank(),
#           panel.background = element_rect(fill = "grey97"),
#           axis.line = element_line(colour = "black")) +
#     geom_vline(xintercept = impact_start, colour = "red") +
#     geom_vline(xintercept = impact_end, colour = "blue")
#   
# }
# 
# scenario_rli_plots_aggregated[[1]]
# 
# # LIVING PLANET INDEX ----
# 
# # * Create folders ----
# 
# lpi_inputs_folder <- file.path(indicator_inputs_folder, "LPI_inputs", today)
# 
# if( !dir.exists( file.path(lpi_inputs_folder) ) ) {
#   dir.create( file.path(lpi_inputs_folder), recursive = TRUE )
#   
# }
# 
# lpi_outputs_folder <- file.path(indicator_outputs_folder, "LPI_outputs", today)
# 
# if( !dir.exists( file.path(lpi_outputs_folder) ) ) {
#   dir.create( file.path(lpi_outputs_folder), recursive = TRUE )
#   
# }
# 
# lpi_plots_folder <- file.path(indicator_plots_folder, "LPI_plots", today)
# 
# if( !dir.exists( file.path(lpi_plots_folder) ) ) {
#   dir.create( file.path(lpi_plots_folder), recursive = TRUE )
#   
# }
# 
# # TEMP CODE ---
# ## Look at the data we are dealing with
# 
# # data <- scenario_abundance_long[[1]][[1]]
# # 
# # head(data)
# # 
# # ggplot(data, aes(x = time_step, y = abundance,
# #                  col = group_id)) +
# #           geom_line()  + 
# #           geom_text(aes(label= group_id),hjust=0, vjust=0) +
# #           theme(legend.position = "none")
# 
# # * Sample data ----
# 
# scenario_lpi_inputs <- list()
# 
# for (i in seq_along(scenario_redlist_data_sampled)) {
#   
#   # Get replicates for a single scenario
#   # replicate_abundance_long <- scenario_smoothed_abundance[[i]]
#   
#   replicate_abundance_long <- scenario_redlist_data_sampled[[i]]
#   
#   replicate_lpi_inputs <- list()
#   # Calculate the LPI for each replicate within the scenario
#   for (j in seq_along(replicate_abundance_long)) {
#     
#     replicate_lpi_inputs[[j]] <- replicate_abundance_long[[j]] %>% 
#       dplyr::select(group_id, annual_time_step, 
#                     ave_abundance)
#     
#     
#   }
#   
#   scenario_lpi_inputs[[i]] <- replicate_lpi_inputs
#   
# }
# 
# # lpi_input <- scenario_lpi_inputs[[1]][[2]]
# # head(lpi_input)
# # write.csv(lpi_input, file.path(indicator_outputs_folder, "lpi_input_example_annual.csv"))
# 
# # * Calculate LPI ----
# 
# # Retain naming convention, the LPI just takes the abundance dataframes we
# # already formatted while making the RLI inputs
# 
# # scenario_lpi_inputs <- scenario_abundance_long
# 
# # Loop through each scenario and replicate and calculate the LPI per rep
# 
# scenario_lpi_outputs <- list()
# 
# for (i in seq_along(scenario_lpi_inputs)) {
#   
#   # Get replicates for a single scenario
#   replicate_lpi_inputs <- scenario_lpi_inputs[[i]]
#   
#   replicate_lpi_outputs <- list()
#   
#   # Calculate the LPI for each replicate within the scenario
#   for (j in seq_along(replicate_lpi_inputs)) {
#     
#     replicate_lpi_outputs[[j]] <- calculate_living_planet_index(
#       
#       replicate_lpi_inputs[[j]], start_time_step, ci = FALSE, numboots, j
#     ) 
#     
#     # Save the output LPI data as a csv and rds
#     
#     saveRDS(replicate_lpi_outputs[[j]],
#             file.path(lpi_outputs_folder,
#                       paste(today, scenarios[[i]], "replicate", j,
#                             "LPI_output_data_monthly_smoothing.rds",
#                             sep = "_")))
#     
#     write.csv(replicate_lpi_outputs[[j]],
#               file.path(lpi_outputs_folder,
#                         paste(today, scenarios[[i]], "replicate", j,
#                               "LPI_output_data_monthly_smoothing.rds",
#                               sep = "_")))
#     
#   }
#   
#   scenario_lpi_outputs[[i]] <- replicate_lpi_outputs
#   
# }
# 
# head(scenario_lpi_outputs)[[1]][[1]]
# 
# # * Aggregate all LPI scores ----
# 
# ## Collapse input data so LPI scores for all replicates in one scenario exist in 
# ## a single data frame
# 
# scenario_lpi_outputs_aggregated <- list()
# 
# for (i in seq_along(scenario_lpi_outputs)) {
#   
#   scenario_lpi_outputs_aggregated[[i]] <- do.call(rbind, 
#                                                   scenario_lpi_outputs[[i]]) %>%
#     mutate(scenario = scenarios[[i]]) 
#   
#   scenario_mean_lpi <- scenario_lpi_outputs_aggregated[[i]] %>%
#     group_by(annual_time_step) %>%
#     summarise(indicator_score = mean(indicator_score),
#               ci_lower = mean(ci_lower),
#               ci_upper = mean(ci_upper)) %>%
#     mutate(replicate = 0,# Replicate 0 will always be the mean
#            indicator = "LPI",
#            scenario = scenarios[[i]]) 
#   
#   scenario_lpi_outputs_aggregated[[i]] <-  rbind(scenario_lpi_outputs_aggregated[[i]],
#                                                  scenario_mean_lpi) %>%
#     mutate(replicate = as.factor(replicate)) %>%
#     mutate(level = ifelse(replicate == 0,
#                           "Mean LPI", 
#                           "Replicate LPI"))
#   
# }
# 
# head(scenario_lpi_outputs_aggregated[[1]])
# tail(scenario_lpi_outputs_aggregated[[1]])
# 
# # * Plot LPI replicates individually ----
# 
# scenario_lpi_plots <- list()
# 
# for (i in seq_along(scenario_lpi_outputs)) {
#   
#   replicate_lpi <- scenario_lpi_outputs[[i]]
#   replicate_lpi_plots <- list()
#   
#   for (j in seq_along(replicate_lpi)) {
#     
#     replicate_lpi_plots[[j]] <- plot_living_planet_index(replicate_lpi[[j]],
#                                                          ci = FALSE)
#     
#     
#     ggsave(file.path(lpi_plots_folder, paste(today, scenarios[[i]], 
#                                              "replicate", j, 
#                                              "LPI_aggregated_monthly_smoothing.png",
#                                              sep = "_")),
#            replicate_lpi_plots[[j]],  device = "png")                                   
#     
#   }
#   
#   scenario_lpi_plots[[i]] <- replicate_lpi_plots
#   
# }
# 
# i <- 1
# scenario_lpi_plots[[1]][[i]]
# 
# i <- i + 1
# scenario_lpi_plots[[1]][[i]]
# 
# # * Plot all replicates together ----
# 
# scenario_lpi_plots_aggregated <- list()
# 
# for (i in seq_along(scenario_lpi_outputs_aggregated)){
#   
#   scenario_lpi_plots_aggregated[[i]] <- ggplot(data = scenario_lpi_outputs_aggregated[[i]], 
#                                                aes(x = annual_time_step, 
#                                                    y = indicator_score, 
#                                                    group = replicate,
#                                                    color = level)) +
#     geom_line() +
#     scale_color_manual(values = c("black", "gray62")) + 
#     labs(x = "Time", 
#          y = "Living Planet Index Score") +
#     theme(panel.grid.major = element_blank(),
#           axis.title = element_text(size = 18),
#           axis.text = element_text(size = 18),
#           panel.grid.minor = element_blank(),
#           panel.background = element_rect(fill = "grey97"),
#           axis.line = element_line(colour = "black")) +
#     geom_vline(xintercept = impact_start, colour = "red") +
#     geom_vline(xintercept = impact_end, colour = "blue")
#   
# }
# 
# scenario_lpi_plots_aggregated[[1]]
# 
# # Combine indicators ----
# 
# all_indicators_list <- list(scenario_rli_outputs,
#                             scenario_lpi_outputs)
# 
# names(all_indicators_list) <- c("RLI", "LPI", "harvested")
# 
# saveRDS(all_indicators_list,
#         file.path(indicator_outputs_folder,
#                   paste(today, "all_indicators_output_data_list_monthly_smoothing.rds",
#                         sep = "_")))
# 
# all_lpi <- do.call(rbind, scenario_lpi_outputs_aggregated) %>% 
#   filter(replicate != 0) # Remove the mean so we just have replicates 
# 
# all_rli <- do.call(rbind, scenario_rli_outputs_aggregated) %>% 
#   filter(replicate != 0) # Remove the mean so we just have replicates
# 
# all_indicators <- rbind(all_lpi, all_rli)
# 
# saveRDS(all_indicators,
#         file.path(indicator_outputs_folder,
#                   paste(today, "all_indicators_output_data_monthly_smoothing.rds",
#                         sep = "_")))
# 
# write.csv(all_indicators,
#           file.path(indicator_outputs_folder,
#                     paste(today, "all_indicators_output_data__monthly_smoothing.csv",
#                           sep = "_")))
# 
# # TAKE A 5 YR SAMPLE ----
# 
# # * Take coarser sample ----
# 
# sample_interval <- 5 # make a different number than one to actually sample
# sample_max_timestep <- 300/sample_interval
# 
# scenario_redlist_data_sampled <- list()
# 
# for (i in seq_along(scenario_red_list_data)) {
#   
#   # Get replicate data for a single scenario
#   
#   replicate <- scenario_red_list_data[[i]]
#   
#   replicate_sampled <- list()
#   
#   # For each individual replicate
#   
#   for (j in seq_along(replicate)) {
#     
#     replicate_sampled[[j]] <- replicate[[j]] %>% 
#       slice(which(row_number() %% sample_interval == 0)) %>% 
#       mutate(annual_time_step = seq(1,sample_max_timestep,1))
#     
#     
#   }
#   
#   scenario_redlist_data_sampled[[i]] <- replicate_sampled
#   
# }
# 
# test <- scenario_redlist_data_sampled[[1]][[1]]
# test_group <- test %>% filter(group_id == "13.16.27")
# dim(test_group)
# 
# # * Get harvested group only ----
# 
# scenario_harvested_groups <- list()
# 
# for (i in seq_along(scenario_redlist_data_sampled)) {
#   
#   # Get replicate data for a single scenario
#   
#   replicate_rl_data <- scenario_redlist_data_sampled[[i]]
#   
#   harvested <- list()
#   
#   # For each individual replicate
#   
#   for (j in seq_along(replicate_rl_data)) {
#     
#     
#     if (scenarios[[i]] == "100_Land_Use") {
#       
#       harvested[[j]]  <- replicate_rl_data[[j]] %>% 
#         filter(functional_group_name == "herbivore endotherm"|
#                  functional_group_name == "herbivore ectotherm") 
#       
#     } else if (scenarios[[i]] == "200_Harvesting_carnivores") {
#       
#       harvested[[j]] <- replicate_rl_data[[j]] %>% 
#         filter(functional_group_name == "carnivore endotherm" & mass_lower == 10000|
#                  functional_group_name == "carnivore ectotherm" & mass_lower == 10000) 
#       
#     } else if (scenarios[[i]] == "300_Harvesting_herbivores") {
#       
#       harvested[[j]]  <- replicate_rl_data[[j]] %>% 
#         filter(functional_group_name == "herbivore endotherm" & mass_lower == 10000|
#                  functional_group_name == "herbivore ectotherm" & mass_lower == 10000) 
#       
#     } else if (scenarios[[i]] == "000_Baseline") {
#       
#       harvested[[j]]  <- replicate_rl_data[[j]] %>% 
#         filter(functional_group_name == "carnivore endotherm"|
#                  functional_group_name == "carnivore ectotherm") %>% 
#         filter(mass_lower == 1000)
#       
#     }
#     
#   }
#   
#   scenario_harvested_groups[[i]] <- harvested
#   
# }
# 
# harvested_rep <- scenario_harvested_groups[[1]][[1]]
# 
# ggplot(data = harvested_rep) +
#   geom_line(aes(x = annual_time_step, y = ave_abundance, col = group_id)) +
#   theme(legend.position = "none") +
#   geom_text(aes(x = annual_time_step, y = ave_abundance, label = rl_status))
# 
# # * Plot harvested groups ----
# 
# harvested_plots_folder <- file.path(indicator_plots_folder, "harvested_plots", today)
# 
# if( !dir.exists( file.path(harvested_plots_folder) ) ) {
#   dir.create( file.path(harvested_plots_folder), recursive = TRUE )
#   
# }
# 
# 
# scenario_harvested_plots <- list()
# 
# for ( i in seq_along(scenario_harvested_groups)) {
#   
#   replicate_harvest <- scenario_harvested_groups[[i]]
#   
#   replicate_harvested_plots <- list()
#   
#   # For each individual replicate
#   
#   for (j in seq_along(replicate_harvest)) {
#   
#     replicate_harvested_plots[[j]] <- ggplot(data = replicate_harvest[[j]]) +
#     geom_smooth(aes(x = annual_time_step, 
#                     y = abundance, 
#                     col = functional_group_name)) +
#     theme(legend.position = "bottom") +
#     labs(title = paste(scenarios[[i]], "harvested groups", sep = " "))
#   
#   ggsave(file.path(harvested_plots_folder, paste(today, scenarios[[i]], 
#                                                  "replicate", j - 1,
#                                                  "harvested_reps_averaged.png",
#                                                  sep = "_")),
#          replicate_harvested_plots[[j]],  device = "png")
#   }
# 
#   scenario_harvested_plots[[i]] <- replicate_harvested_plots 
#   
# }
# 
# # Plot some results to check they're not completely whack
# 
# ## Get one group to check how their status changes over time relative to how
# ## their abundance changes
# 
# # group_id_select <- "13.16.17" # Shows example of 'resurrected' virtual spp
# # # group_id_select <- "10.40"
# # 
# # data <- rli_inputs %>% dplyr::filter(group_id == group_id_select)
# # 
# # ggplot(data, aes(x = time_step, y = abundance)) +
# #   geom_line() +
# #   geom_text(aes(label= rl_status,
# #                 col = rl_status),hjust=0, vjust=0)
# 
# 
# 
# 
# # * Calculate RLI ----
# 
# # RLI by individual functional groups
# 
# scenario_fg_rli_outputs <- list()
# 
# for (i in seq_along(scenario_redlist_data_sampled)) {
#   
#   replicate_red_list_inputs <- scenario_redlist_data_sampled[[i]]
#   
#   replicate_fg_rli_outputs <- list()
#   
#   for (j in seq_along(replicate_red_list_inputs)) {
#     
#     replicate_rli <- calculate_red_list_index(
#       replicate_red_list_inputs[[j]], numboots, ci = FALSE) %>%
#       mutate(replicate = j)
#     
#     replicate_fg_rli_outputs[[j]] <- replicate_rli 
#     
#     # saveRDS(replicate_fg_rli_outputs[[j]],
#     #         file.path(rli_outputs_folder,
#     #                   paste(today, scenarios, "replicate", j,
#     #                         "RLI_func_group_output_data.rds",
#     #                         sep = "_")))
#     
#     write.csv(replicate_fg_rli_outputs[[j]],
#               file.path(rli_outputs_folder,
#                         paste(today, scenarios[[i]], "RLI_func_group_output_data_5yrs.rds",
#                               sep = "_")))
#     
#     print(paste("RLI for replicate", j, "complete", sep = " "))
#     
#   }
#   
#   scenario_fg_rli_outputs[[i]] <- replicate_fg_rli_outputs
#   
# }
# 
# 
# x <- scenario_fg_rli_outputs[[1]][[3]]
# head(x)
# 
# # Mean RLI aggregated across groups
# 
# scenario_rli_outputs <- list()
# 
# for (i in seq_along(scenario_fg_rli_outputs)) {
#   
#   # Get replicates for a single scenario
#   replicate_rli_fg <- scenario_fg_rli_outputs[[i]]
#   
#   replicate_rli_outputs <- list()
#   
#   # Aggregate RLI across functional groups for each replicate
#   for (j in seq_along(replicate_rli_fg)) {
#     
#     if ("ci_lower" %in% names(replicate_rli_fg[[j]])) {
#       
#       replicate_rli_outputs[[j]] <- replicate_rli_fg[[j]] %>%
#         group_by(annual_time_step) %>%
#         summarise(indicator_score = mean(indicator_score),
#                   ci_lower = mean(ci_lower),
#                   ci_upper = mean(ci_upper)) %>%
#         mutate(indicator = "RLI",
#                replicate = j)
#     } else {
#       
#       replicate_rli_outputs[[j]] <- replicate_rli_fg[[j]] %>%
#         group_by(annual_time_step) %>%
#         summarise(indicator_score = mean(indicator_score)) %>%
#         mutate(indicator = "RLI",
#                replicate = j)
#     }
#     
#     # saveRDS(replicate_rli_outputs[[j]],
#     #       file.path(rli_outputs_folder,
#     #                 paste(today, scenarios[[i]], "replicate", j,
#     #                       "RLI_aggregate_output_data.rds",
#     #                       sep = "_")))
#     # 
#     # write.csv(replicate_rli_outputs[[j]],
#     #           file.path(rli_outputs_folder,
#     #                     paste(today, scenarios[[i]], "replicate", j,
#     #                           "RLI_aggregate_output_data.rds",
#     #                           sep = "_")))
#     
#   }
#   
#   scenario_rli_outputs[[i]] <- replicate_rli_outputs
#   
# }
# 
# head(scenario_rli_outputs)[[1]][[1]]
# 
# # * Plot RLI ----
# 
# ## By functional group
# 
# scenario_fg_rli_plots <- list()
# 
# for (i in seq_along(scenario_fg_rli_outputs)) {
#   
#   replicate_fg_rli <- scenario_fg_rli_outputs[[i]]
#   
#   replicate_fg_rli_plots <- list()
#   
#   for (j in seq_along(replicate_fg_rli)) {
#     
#     replicate_fg_rli_plots[[j]] <-  plot_red_list_index_by_group(
#       replicate_fg_rli[[j]],
#       20,
#       40,
#       ci = FALSE)
#     
#     ggsave(file.path(rli_plots_folder, paste(today, scenarios[[i]], "replicate", j,
#                                              "RLI_by_functional_group_5yrs.png",
#                                              sep = "_")),
#            replicate_fg_rli_plots[[j]],  device = "png")
#     
#   }
#   
#   scenario_fg_rli_plots[[i]] <- replicate_fg_rli_plots
#   
# }
# 
# scenario_fg_rli_plots[[1]][[5]]
# 
# 
# # Small test to see if averaging indicator scores after works better (it doesn't)
# x <- scenario_rli_outputs[[1]][[5]]
# x <- x[-1,]
# 
# x <- x %>% 
#   mutate(x = rollmean(indicator_score, 10, na.pad = TRUE))
# 
# ggplot(x, aes(x = annual_time_step, y = x))+
#   geom_line()
# 
# # RLI with all functional groups aggregated
# # i.e. mean of each 'taxa' RLI as per Butchart et al (2010) 'Indicators of
# # recent declines'
# 
# scenario_rli_plots <- list()
# 
# for (i in seq_along(scenario_rli_outputs)) {
#   
#   replicate_rli <- scenario_rli_outputs[[i]]
#   
#   replicate_rli_plots <- list()
#   
#   for (j in seq_along(replicate_rli)) {
#     
#     replicate_rli_plots[[j]] <- plot_red_list_index(replicate_rli[[j]],
#                                                     20, 
#                                                     40,
#                                                     ci = TRUE)
#     
#     
#     ggsave(file.path(rli_plots_folder, paste(today, scenarios[[i]], 
#                                              "replicate", j, 
#                                              "RLI_aggregated_5yrs.png",
#                                              sep = "_")),
#            replicate_rli_plots[[j]],  device = "png")                                   
#     
#   }
#   
#   scenario_rli_plots[[i]] <- replicate_rli_plots
#   
# }
# 
# i <- 1
# i <- i+1
# scenario_rli_plots[[1]][[i]]
# 
# # * Plot all replicates together ----
# 
# ## Collapse input data so RLI scores for all replicates in one scenario exist in 
# ## a single data frame
# 
# scenario_rli_outputs_aggregated <- list()
# 
# for (i in seq_along(scenario_rli_outputs)) {
#   
#   scenario_rli_outputs_aggregated[[i]] <- do.call(rbind, 
#                                                   scenario_rli_outputs[[i]]) %>%
#     mutate(scenario = scenarios[[i]]) 
#   
#   
#   scenario_mean_rli <- scenario_rli_outputs_aggregated[[i]] %>%
#     group_by(annual_time_step) %>%
#     summarise(indicator_score = mean(indicator_score),
#               ci_lower = mean(ci_lower),
#               ci_upper = mean(ci_upper)) %>%
#     mutate(indicator = "RLI",
#            replicate = 0,
#            scenario = scenarios[[i]]) # Replicate 0 will always be the mean
#   
#   scenario_rli_outputs_aggregated[[i]] <-  rbind(scenario_rli_outputs_aggregated[[i]],
#                                                  scenario_mean_rli) %>%
#     mutate(replicate = as.factor(replicate)) %>%
#     mutate(level = ifelse(replicate == 0,
#                           "Mean RLI", 
#                           "Replicate RLI"),
#            scenario = scenarios[[i]])
#   
# }
# 
# head(scenario_rli_outputs_aggregated[[1]])
# tail(scenario_rli_outputs_aggregated[[1]])
# 
# # Plot all together
# 
# scenario_rli_plots_aggregated <- list()
# 
# for (i in seq_along(scenario_rli_outputs_aggregated)) {
#   
#   scenario_rli_plots_aggregated[[i]] <- ggplot(data = scenario_rli_outputs_aggregated[[i]], 
#                                                aes(x = annual_time_step, 
#                                                    y = indicator_score, 
#                                                    group = replicate,
#                                                    color = level)) +
#     geom_line() +
#     scale_color_manual(values = c("black", "gray62")) + 
#     labs(x = "Time", 
#          y = "Red List Index Score") +
#     theme(panel.grid.major = element_blank(),
#           axis.title = element_text(size = 18),
#           axis.text = element_text(size = 18),
#           panel.grid.minor = element_blank(),
#           panel.background = element_rect(fill = "grey97"),
#           axis.line = element_line(colour = "black")) +
#     geom_vline(xintercept = 20, colour = "red") +
#     geom_vline(xintercept = 40, colour = "blue")
#   
# }
# 
# scenario_rli_plots_aggregated[[1]]
# 
# # LIVING PLANET INDEX ----
# 
# # * Create folders ----
# 
# lpi_inputs_folder <- file.path(indicator_inputs_folder, "LPI_inputs", today)
# 
# if( !dir.exists( file.path(lpi_inputs_folder) ) ) {
#   dir.create( file.path(lpi_inputs_folder), recursive = TRUE )
#   
# }
# 
# lpi_outputs_folder <- file.path(indicator_outputs_folder, "LPI_outputs", today)
# 
# if( !dir.exists( file.path(lpi_outputs_folder) ) ) {
#   dir.create( file.path(lpi_outputs_folder), recursive = TRUE )
#   
# }
# 
# lpi_plots_folder <- file.path(indicator_plots_folder, "LPI_plots", today)
# 
# if( !dir.exists( file.path(lpi_plots_folder) ) ) {
#   dir.create( file.path(lpi_plots_folder), recursive = TRUE )
#   
# }
# 
# # TEMP CODE ---
# ## Look at the data we are dealing with
# 
# # data <- scenario_abundance_long[[1]][[1]]
# # 
# # head(data)
# # 
# # ggplot(data, aes(x = time_step, y = abundance,
# #                  col = group_id)) +
# #           geom_line()  + 
# #           geom_text(aes(label= group_id),hjust=0, vjust=0) +
# #           theme(legend.position = "none")
# 
# # * Sample data ----
# 
# scenario_lpi_inputs <- list()
# 
# for (i in seq_along(scenario_redlist_data_sampled)) {
#   
#   # Get replicates for a single scenario
#   # replicate_abundance_long <- scenario_smoothed_abundance[[i]]
#   
#   replicate_abundance_long <- scenario_redlist_data_sampled[[i]]
#   
#   replicate_lpi_inputs <- list()
#   # Calculate the LPI for each replicate within the scenario
#   for (j in seq_along(replicate_abundance_long)) {
#     
#     replicate_lpi_inputs[[j]] <- replicate_abundance_long[[j]] %>% 
#       dplyr::select(group_id, annual_time_step, 
#                     ave_abundance)
#     
#     
#   }
#   
#   scenario_lpi_inputs[[i]] <- replicate_lpi_inputs
#   
# }
# 
# # lpi_input <- scenario_lpi_inputs[[1]][[2]]
# # head(lpi_input)
# # write.csv(lpi_input, file.path(indicator_outputs_folder, "lpi_input_example_annual.csv"))
# 
# # * Calculate LPI ----
# 
# # Retain naming convention, the LPI just takes the abundance dataframes we
# # already formatted while making the RLI inputs
# 
# # scenario_lpi_inputs <- scenario_abundance_long
# 
# # Loop through each scenario and replicate and calculate the LPI per rep
# 
# scenario_lpi_outputs <- list()
# 
# for (i in seq_along(scenario_lpi_inputs)) {
#   
#   # Get replicates for a single scenario
#   replicate_lpi_inputs <- scenario_lpi_inputs[[i]]
#   
#   replicate_lpi_outputs <- list()
#   
#   # Calculate the LPI for each replicate within the scenario
#   for (j in seq_along(replicate_lpi_inputs)) {
#     
#     replicate_lpi_outputs[[j]] <- calculate_living_planet_index(
#       
#       replicate_lpi_inputs[[j]], start_time_step, ci = FALSE, numboots, j
#     ) 
#     
#     # Save the output LPI data as a csv and rds
#     
#     saveRDS(replicate_lpi_outputs[[j]],
#             file.path(lpi_outputs_folder,
#                       paste(today, scenarios[[i]], "replicate", j,
#                             "LPI_output_data_5yrs.rds",
#                             sep = "_")))
#     
#     write.csv(replicate_lpi_outputs[[j]],
#               file.path(lpi_outputs_folder,
#                         paste(today, scenarios[[i]], "replicate", j,
#                               "LPI_output_data_5yrs.rds",
#                               sep = "_")))
#     
#   }
#   
#   scenario_lpi_outputs[[i]] <- replicate_lpi_outputs
#   
# }
# 
# head(scenario_lpi_outputs)[[1]][[1]]
# 
# # * Aggregate all LPI scores ----
# 
# ## Collapse input data so LPI scores for all replicates in one scenario exist in 
# ## a single data frame
# 
# scenario_lpi_outputs_aggregated <- list()
# 
# for (i in seq_along(scenario_lpi_outputs)) {
#   
#   scenario_lpi_outputs_aggregated[[i]] <- do.call(rbind, 
#                                                   scenario_lpi_outputs[[i]]) %>%
#     mutate(scenario = scenarios[[i]]) 
#   
#   scenario_mean_lpi <- scenario_lpi_outputs_aggregated[[i]] %>%
#     group_by(annual_time_step) %>%
#     summarise(indicator_score = mean(indicator_score),
#               ci_lower = mean(ci_lower),
#               ci_upper = mean(ci_upper)) %>%
#     mutate(replicate = 0,# Replicate 0 will always be the mean
#            indicator = "LPI",
#            scenario = scenarios[[i]]) 
#   
#   scenario_lpi_outputs_aggregated[[i]] <-  rbind(scenario_lpi_outputs_aggregated[[i]],
#                                                  scenario_mean_lpi) %>%
#     mutate(replicate = as.factor(replicate)) %>%
#     mutate(level = ifelse(replicate == 0,
#                           "Mean LPI", 
#                           "Replicate LPI"))
#   
# }
# 
# head(scenario_lpi_outputs_aggregated[[1]])
# tail(scenario_lpi_outputs_aggregated[[1]])
# 
# # * Plot LPI replicates individually ----
# 
# scenario_lpi_plots <- list()
# 
# for (i in seq_along(scenario_lpi_outputs)) {
#   
#   replicate_lpi <- scenario_lpi_outputs[[i]]
#   replicate_lpi_plots <- list()
#   
#   for (j in seq_along(replicate_lpi)) {
#     
#     replicate_lpi_plots[[j]] <- plot_living_planet_index(replicate_lpi[[j]],
#                                                          ci = FALSE)
#     
#     
#     ggsave(file.path(lpi_plots_folder, paste(today, scenarios[[i]], 
#                                              "replicate", j, 
#                                              "LPI_aggregated_5yrs.png",
#                                              sep = "_")),
#            replicate_lpi_plots[[j]],  device = "png")                                   
#     
#   }
#   
#   scenario_lpi_plots[[i]] <- replicate_lpi_plots
#   
# }
# 
# i <- 1
# scenario_lpi_plots[[1]][[i]]
# 
# i <- i + 1
# scenario_lpi_plots[[1]][[i]]
# 
# # * Plot all replicates together ----
# 
# scenario_lpi_plots_aggregated <- list()
# 
# for (i in seq_along(scenario_lpi_outputs_aggregated)){
#   
#   scenario_lpi_plots_aggregated[[i]] <- ggplot(data = scenario_lpi_outputs_aggregated[[i]], 
#                                                aes(x = annual_time_step, 
#                                                    y = indicator_score, 
#                                                    group = replicate,
#                                                    color = level)) +
#     geom_line() +
#     scale_color_manual(values = c("black", "gray62")) + 
#     labs(x = "Time", 
#          y = "Living Planet Index Score") +
#     theme(panel.grid.major = element_blank(),
#           axis.title = element_text(size = 18),
#           axis.text = element_text(size = 18),
#           panel.grid.minor = element_blank(),
#           panel.background = element_rect(fill = "grey97"),
#           axis.line = element_line(colour = "black")) +
#     geom_vline(xintercept = 20, colour = "red") +
#     geom_vline(xintercept = 40, colour = "blue")
#   
# }
# 
# scenario_lpi_plots_aggregated[[1]]
# 
# # Combine indicators ----
# 
# all_indicators_list <- list(scenario_rli_outputs,
#                             scenario_lpi_outputs)
# 
# names(all_indicators_list) <- c("RLI", "LPI")
# 
# saveRDS(all_indicators_list,
#         file.path(indicator_outputs_folder,
#                   paste(today, "all_indicators_output_data_list_5yrs.rds",
#                         sep = "_")))
# 
# all_lpi <- do.call(rbind, scenario_lpi_outputs_aggregated) %>% 
#   filter(replicate != 0) # Remove the mean so we just have replicates 
# 
# all_rli <- do.call(rbind, scenario_rli_outputs_aggregated) %>% 
#   filter(replicate != 0) # Remove the mean so we just have replicates
# 
# all_indicators <- rbind(all_lpi, all_rli)
# 
# saveRDS(all_indicators,
#         file.path(indicator_outputs_folder,
#                   paste(today, "all_indicators_output_data_5yrs.rds",
#                         sep = "_")))
# 
# write.csv(all_indicators,
#           file.path(indicator_outputs_folder,
#                     paste(today, "all_indicators_output_data__5yrs.csv",
#                           sep = "_")))

# ## AVERAGE REPLICATES ----

# Get generation length ----

scenario_ab_gl_formatted_not_clean <- list()

for (i in seq_along(scenario_abundance_long)) {
  
  replicate_abundance <- scenario_abundance_long[[i]]
  replicate_generations <- scenario_generations_raw[[i]]
  
  # Make a list to catch the outputs
  
  replicate_ab_gl_formatted <- list()
  
  # For each individual replicate
  
  for (j in seq_along(replicate_abundance)) {
    
    # Reduce size of the replicate generations dataframe or the merge won't work
    gen_length <- replicate_generations[[j]] %>% 
      dplyr::select(group_id, generation_length_yrs, 
                    functional_group_name) %>% 
      distinct(.)
    
    # Add the generation length info to the abundance dataframe
    replicate_ab_gl_formatted[[j]] <- replicate_abundance[[j]] %>%
      merge(gen_length, by = "group_id") %>%
      arrange(monthly_time_step, group_id) %>%
      # Get the timeframe over which to assess decline (3 * gen length or 10 yrs,
      # whichever is longer)
      # Important - following lines assume an annual timeframe, will need to adjust if change interval
      mutate(generation_by_three = generation_length_yrs * 3) %>% # Time over which to measure decline, 3 x gen length OR:
      mutate(timeframe = generation_by_three) %>%
      #round the time frame to whole years
      mutate(timeframe = round(timeframe)) %>% 
      dplyr::select(-generation_by_three) %>%
      distinct(.) %>%
      group_by(group_id) %>% 
      # Add info about bodymass so can easily sort
      merge(groups[c("group_id", "bodymass_index", "mass_lower")], 
            by = "group_id") %>% 
      arrange(monthly_time_step)
    # select rows that are multiples of the specified interval 
    # (eg if interval is 12, it samples one month from every 12 (yearly))
    # slice(which(row_number() %% interval == 0)) %>% 
    # mutate(annual_time_step = seq(1,max_timestep,1)) # %>% 
    
    
    print(paste("Replicate", j - 1, 
                "formatting complete", 
                sep = " "))
    
  }
  
  print(scenario[[i]])
  print(length(replicate_ab_gl_formatted))
  
  scenario_ab_gl_formatted_not_clean[[i]] <- replicate_ab_gl_formatted
  
}

#today <- "2021-08-10"
test <- scenario_ab_gl_formatted_not_clean[[1]][[1]]
head(test)
# Remove false extinctions ----

scenario_false_extinctions_removed <- list()

for (i in seq_along(scenario_ab_gl_formatted_not_clean)) {
  
  replicate_ab_gl <- scenario_ab_gl_formatted_not_clean[[i]]
  
  # Make a list to catch the outputs
  
  replicate_false_ex_removed <- list()
  
  # For each individual replicate
  
  for (j in seq_along(replicate_ab_gl)) {
    
    # Find the last time step where non-0 abundance occurred for each group
    
    temp2 <- replicate_ab_gl[[j]] %>% 
      group_by(group_id) %>% 
      filter(abundance > 0) %>% 
      dplyr::select(group_id, monthly_time_step, abundance) %>% 
      filter(monthly_time_step == max(monthly_time_step)) %>% 
      dplyr::select(group_id, monthly_time_step) %>% 
      rename(last_abundance = monthly_time_step)
    
    # Add the year of last positive abundance number as a column to the data    
    temp3 <- replicate_ab_gl[[j]] %>% 
      merge(temp2, by = c("group_id"), all = TRUE)
    
    # Use the last positive abundance year and current abundance value to determine
    # if a zero abundance is a true extinction or just a missing value (false extinction)
    data <- temp3 %>%
      group_by(group_id) %>%
      mutate(true_extinction = ifelse(abundance == 0 &
                                        monthly_time_step < last_abundance,
                                      "false extinction",
                                      ifelse(abundance > 0 &
                                               monthly_time_step < last_abundance,
                                             "not extinct",
                                             ifelse(abundance == 0 &
                                                      monthly_time_step >= last_abundance,
                                                    "true extinction", "not extinct")))) %>%
      #filter(true_extinction != "false extinction") %>%
      mutate(abundance = ifelse(true_extinction == "false extinction",
                                NA, abundance)) %>%
      group_by(group_id) %>%
      arrange(group_id, monthly_time_step)
    
    # Check if there are any carnivorous endotherms
    
    check <- data %>% 
      group_by(functional_group_name) %>% 
      filter(monthly_time_step > 13200) %>% 
      summarise(present = sum(abundance, na.rm = TRUE)) %>% 
      filter(functional_group_name == "carnivore endotherm") %>% 
      dplyr::select(present) %>% 
      pull(.)
    
    
    print(paste("Replicate", j - 1, 
                "formatting complete", 
                sep = " "))
    
    # Replace data with 0 if no carnivores
    
    if(check == 0) {
      
      data <- NULL
      
      print(paste("Replicate", j - 1, 
                  "removed because no carnivorous endotherms are present", 
                  sep = " "))
      
    }
    
    replicate_false_ex_removed[[j]] <- data
    
  }
  
  print(scenario[[i]])
  print(length(replicate_ab_gl))
  
  scenario_false_extinctions_removed[[i]] <- replicate_false_ex_removed
  
}

test <- scenario_false_extinctions_removed[[1]][[2]]
head(test)
any(is.nan(test$abundance))

# Remove replicates with no carnivorous endotherms ----

scenario_ab_gl_formatted <- list()

for (i in seq_along(scenario_false_extinctions_removed)) {
  
  replicate_not_clean <- scenario_false_extinctions_removed[[i]]
  
  scenario_ab_gl_formatted[[i]] <- list.clean(replicate_not_clean)
  
}

test <- scenario_ab_gl_formatted[[1]][[1]]
head(test)
any(is.nan(test$abundance))


# Take a random sample of 25 replicates ----

# scenario_ab_gl_formatted_25 <- list()
# scenario_auto_long_25 <- list()
# 
# for (i in seq_along(scenario_ab_gl_formatted)) {
# 
#   if(scenarios[[i]] == "000_Baseline") {
# 
#   scenario_ab_gl_formatted_25[[i]] <- scenario_ab_gl_formatted[[i]]
#   scenario_auto_long_25[[i]] <- scenario_auto[[i]]
# 
#   } else {
# 
#   replicates <- scenario_ab_gl_formatted[[i]]
#   replicates_auto <- scenario_auto[[i]]
# 
#   set.seed(159)
# 
#   scenario_ab_gl_formatted_25[[i]] <- sample(scenario_ab_gl_formatted[[i]], 25)
#   scenario_auto_long_25[[i]] <- sample(scenario_auto[[i]], 25) # Will this sample the same reps??
# 
#   }
# 
# }

scenario_ab_gl_formatted_25 <- scenario_ab_gl_formatted
scenario_auto_long_25 <- scenario_auto_long

# Average across replicates ----

scenario_averaged <- list()
scenario_auto_averaged <- list()

for (i in seq_along(scenario_ab_gl_formatted_25)) {
  
  # Heterotrophs
  
  replicates <- scenario_ab_gl_formatted_25[[i]]
  
  n <- length(replicates)
  
  replicate_df <- do.call(rbind, replicates)
  
  scenario_averaged[[i]] <- replicate_df %>% 
                            mutate(replicate = j - 1) %>% 
                            group_by(group_id, monthly_time_step) %>% 
                            mutate(mean_abundance = mean(abundance, na.rm = TRUE),
                                   # 95% CIs using this info: https://www.cyclismo.org/tutorial/R/confidence.html
                                   sd = sd(abundance, na.rm = TRUE),
                                   error = qt(0.975, df = n - 1) * 
                                              sd/sqrt(n),
                                   lower_ci = mean_abundance - error, 
                                   upper_ci = mean_abundance + error,  
                                   mean_gen_length = mean(generation_length_yrs, na.rm = TRUE),
                                   mean_timeframe = round(mean(timeframe, na.rm = TRUE))) %>%
                            dplyr::select(-abundance, -generation_length_yrs, -timeframe,
                                          -true_extinction, - replicate, - last_abundance) %>% 
                            rename(abundance = mean_abundance,
                                   generation_length_yrs = mean_gen_length,
                                   timeframe = mean_timeframe) %>% 
                            distinct(.) %>% 
                            ungroup(.)
  
  #Autotrophs
  
  auto_replicates <- scenario_auto_long_25[[i]]
  
  n2 <- length(auto_replicates)
  
  auto_replicate_df <- do.call(rbind, auto_replicates)
  
  scenario_auto_averaged[[i]] <- auto_replicate_df %>% 
    group_by(group_id, monthly_time_step) %>% 
    mutate(mean_abundance = mean(abundance, na.rm = TRUE),
           # 95% CIs using this info: https://www.cyclismo.org/tutorial/R/confidence.html
           sd = sd(abundance, na.rm = TRUE),
           error = qt(0.975, df = n - 1) * 
             sd/sqrt(n),
           lower_ci = mean_abundance - error, 
           upper_ci = mean_abundance + error) %>%
    dplyr::select(-abundance, - replicate) %>% 
    rename(abundance = mean_abundance) %>% 
    distinct(.) %>% 
    ungroup(.)
  
  
  print(paste(scenarios[[i]], "replicates averaged"))
  
}

test2 <- scenario_averaged[[4]]
head(test2)
any(is.nan(test2$abundance))
group <- test2 %>%  filter(group_id == "15.18.10")

test3 <- scenario_auto_averaged[[4]]
head(test3)
# Remove blinking groups ----

scenario_abundance_clean <- list()

for (i in seq_along(scenario_averaged)) {
  
  
    data <- scenario_averaged[[i]]
    
    # gen <- replicate_gens[[j]] %>%
    #   dplyr::select(group_id, mass_lower_g) %>%
    #   distinct(.)
    
    # data <- data %>% filter(group_id == "15.18.21") # For testing
    
    # Determine which groups were there at beginning (post burnin)
    temp <- data %>%
      group_by(group_id) %>%
      filter(abundance != is.nan(abundance)) %>% 
      filter(monthly_time_step == min(monthly_time_step)) %>%
      mutate(first_appearance = monthly_time_step,
             beginning = ifelse(first_appearance == 12000,
                                TRUE, FALSE)) %>%
      dplyr::select(group_id, first_appearance, beginning)
    
    scenario_abundance_clean[[i]] <- data %>%
      merge(temp, by = "group_id") %>%
            arrange(mass_lower) %>%
      tidylog::filter(beginning == TRUE) # Note 14% is highest percentage of data removed by this line
    
    rm(temp)
}

# Plot average abundance ----

scenario_abundance_line_plots <- list()
scenario_abundance_smooth_plots <- list()

for (i in seq_along(scenario_abundance_clean)) {
  
  plot_data <- scenario_abundance_clean[[i]] %>% 
               group_by(monthly_time_step, functional_group_name) %>% 
               mutate(total_abundance = mean(abundance, na.rm = TRUE),
                      total_ci_lower = mean(lower_ci, na.rm = TRUE),
                      total_ci_upper = mean(upper_ci, na.rm = TRUE)) %>% 
               dplyr::select(monthly_time_step, total_abundance,
                             total_ci_lower, total_ci_upper) %>% 
               distinct(.) %>% 
               ungroup(.) %>% 
               slice(which(row_number() %% interval == 0))
  
  scenario_abundance_line_plots[[i]] <- ggplot(plot_data) +
    geom_line(aes(x = monthly_time_step, y = total_abundance)) +
    # geom_ribbon(aes(x = monthly_time_step, 
    #                 ymin = total_ci_lower, 
    #                 ymax = total_ci_upper,
    #                 fill = "band"),
    #             alpha = 0.8) +
    theme(axis.text.y = element_blank()) +
    labs(title = scenarios[[i]]) 
  
  scenario_abundance_smooth_plots[[i]] <- ggplot(plot_data) +
    geom_smooth(aes(x = monthly_time_step, y = total_abundance)) +
    theme(axis.text.y = element_blank()) +
    labs(title = scenarios[[i]]) 
    
}

scenario_abundance_line_plots[[1]]
scenario_abundance_smooth_plots[[4]]

# Plot functional group average abundance - doesn't work well yet

scenario_fg_line_plots <- list()
scenario_fg_smooth_plots <- list()

for (i in seq_along(scenario_abundance_clean)) {
  
 plot_data <- scenario_abundance_clean[[i]] %>% 
    group_by(functional_group_name, monthly_time_step) %>% 
    mutate(total_abundance = sum(abundance, na.rm = TRUE)) %>%
    dplyr::select(monthly_time_step, total_abundance) %>% 
    distinct(.) %>% 
    ungroup(.) %>% 
    slice(which(row_number() %% interval == 0))
 
  scenario_fg_line_plots[[i]] <- ggplot(plot_data) +
    geom_line(aes(x = monthly_time_step, y = total_abundance)) +
    labs(title = scenarios[[i]]) +
    facet_wrap(~ functional_group_name)
  
  scenario_fg_smooth_plots[[i]] <- ggplot(plot_data) +
    geom_smooth(aes(x = monthly_time_step, y = total_abundance)) +
    theme(axis.text.y = element_blank()) +
    labs(title = scenarios[[i]]) +
    facet_wrap(~ functional_group_name)
  
}

scenario_fg_line_plots[[1]]
scenario_fg_smooth_plots[[4]]

# # Smooth abundance ----

ave_window <- 12

scenario_smoothed_abundance <- list()

for (i in seq_along(scenario_abundance_clean)) {

    group_list <- split(scenario_abundance_clean[[i]], 
                        scenario_abundance_clean[[i]]$group_id)

    group_smoothed_abundance <- list()

    for (j in seq_along(group_list)) {

      group_df <- group_list[[j]]

      # check if the group has any abundance values, make it null if not

      if (sum(group_df$abundance, na.rm = TRUE) == 0) {

        group_smoothed_abundance[[j]] <- NULL

      } else {

        group_smoothed_abundance[[j]] <- group_df %>%
          arrange(monthly_time_step) %>%
          mutate(ave_abundance = rollapply(abundance,
                                           ave_window,
                                           mean,
                                           na.rm = TRUE,
                                           partial = TRUE,
                                           align = "left"),
                 ave_abundance = ifelse(ave_abundance < 1,
                                        0, ave_abundance))

        print(j)

      }
    }

    scenario_smoothed_abundance[[i]]  <- do.call(rbind,
                                                 group_smoothed_abundance)
}

saveRDS(indicator_outputs_folder, "scenario_smoothed_abundance.rds")

check <- scenario_smoothed_abundance[[1]]
head(check)

checkgroup <- check %>% filter(group_id == "13.16.27")

# Smooth auto abundance

scenario_smoothed_auto_abundance <- list()

for (i in seq_along(scenario_auto_averaged)) {
  
  group_list <- split(scenario_auto_averaged[[i]], 
                      scenario_auto_averaged[[i]]$group_id)
  
  group_smoothed_auto <- list()
  
  for (j in seq_along(group_list)) {
    
    group_df <- group_list[[j]]
    
    # check if the group has any abundance values, make it null if not
    
    if (sum(group_df$abundance, na.rm = TRUE) == 0) {
      
      group_smoothed_auto[[j]] <- NULL
      
    } else {
      
      group_smoothed_auto[[j]] <- group_df %>%
        arrange(monthly_time_step) %>%
        mutate(ave_abundance = rollapply(abundance,
                                         ave_window,
                                         mean,
                                         na.rm = TRUE,
                                         partial = TRUE,
                                         align = "left"),
               ave_abundance = ifelse(ave_abundance < 1,
                                      0, ave_abundance))
      
      print(j)
      
    }
  }
  
  scenario_smoothed_auto_abundance[[i]] <- do.call(rbind, group_smoothed_auto)
}

check <- scenario_smoothed_auto_abundance[[1]]
head(check)

# RED LIST INDEX ----

# * Create folders ----

rli_inputs_folder <- file.path(indicator_inputs_folder, "RLI_inputs", today)

if( !dir.exists( file.path(rli_inputs_folder) ) ) {
  dir.create( file.path(rli_inputs_folder), recursive = TRUE )
  
}

rli_outputs_folder <- file.path(indicator_outputs_folder, "RLI_outputs", today)

if( !dir.exists( file.path(rli_outputs_folder) ) ) {
  dir.create( file.path(rli_outputs_folder), recursive = TRUE )
  
}

rli_plots_folder <- file.path(indicator_plots_folder, "RLI_plots", today)

if( !dir.exists( file.path(rli_plots_folder) ) ) {
  dir.create( file.path(rli_plots_folder), recursive = TRUE )
  
}

# * Take an annual sample ----

interval <- 12

scenario_annual <- list()

for (i in seq_along(scenario_smoothed_abundance)) {
  
  scenario_groups <- split(scenario_smoothed_abundance[[i]], 
                           scenario_smoothed_abundance[[i]]$group_id)
  
  groups_annual <- list()
  
 
  for (j in seq_along(scenario_groups)) {
    
    groups_annual[[j]] <- scenario_groups[[j]] %>% 
      slice(which(row_number() %% interval == 0)) %>%
      mutate(annual_time_step = seq(1,max_timestep,1))
    
    }
  
  groups_annual_df <- do.call(rbind, groups_annual)
  
  scenario_annual[[i]] <- groups_annual_df
  
}

check <- scenario_annual[[1]]
head(check)
dim(check)

# Sample autotrophs

scenario_auto_annual <- list()

for (i in seq_along(scenario_smoothed_abundance)) {
  
  scenario_auto_groups <- split(scenario_smoothed_auto_abundance[[i]], 
                                scenario_smoothed_auto_abundance[[i]]$group_id)
  
  groups_auto_annual <- list()
  
  for (j in seq_along(scenario_auto_groups)) {
    
 
    groups_auto_annual[[j]] <- scenario_auto_groups[[j]] %>% 
      slice(which(row_number() %% interval == 0)) %>%
      mutate(annual_time_step = seq(1,max_timestep,1))
    
  }
  
  groups_annual_auto_df <- do.call(rbind, groups_auto_annual)
  
  scenario_auto_annual[[i]] <- groups_annual_auto_df

}

check <- scenario_auto_annual[[1]]
head(check)
length(scenario_auto_annual)

## Referring to the thresholds quote under Criterion A, Reason 1 (declines
## are the result of reversible pressures) according to:
## https://portals.iucn.org/library/sites/library/files/documents/RL-2001-001-2nd.pdf


# * Assign Red List Categories ----

scenario_red_list_data <- list()

#for (i in seq_along(scenario_ab_gl_formatted)) {
for (i in seq_along(scenario_annual)) {
  
  # Get replicate data for a single scenario
  
  print(paste("Processing scenario", scenarios[[i]], sep = " "))
  
  # Split by functional group, because we calculate RLI for different
  # functional groups then aggregate later (as per Butchart etal 2010),
  # except we are using functional groups as proxies for taxa (eg mammals, birds, 
  # reptiles) used in real world RLI calcs
  
  status_inputs <- split(scenario_annual[[i]], 
                         scenario_annual[[i]]$group_id)
  
  # Make a list to hold output for each individual massbin-func-group (ie virtual spp)
  
  group_red_list_data <- list()
  
  for (j in seq_along(status_inputs)) {
    
    print(paste("Processing group", names(status_inputs)[[j]], sep = " "))
    
    group_red_list_data[[j]] <- status_inputs[[j]] %>%
      group_by(group_id) %>%
      arrange(monthly_time_step) %>%
      # calculate the difference in abundance over 10 yrs or 3 generation lengths
      # (specified by 'timeframe' column). Its okay to take the first value of 
      # timeframe bc the dataframe is grouped by group_id, and timeframe only changes
      # between and not within group_ids
      # mutate(diff = (abundance - dplyr::lag(abundance, timeframe[1]))) %>%
      mutate(diff = (ave_abundance - dplyr::lag(ave_abundance, timeframe[1]))) %>%
      # Using the formula from p 35 (Complex patterns of decline) Guidelines 
      # for Using the IUCN Red List Categories and Criteria v14 August 2019 
      mutate(decline = 1 - ave_abundance/dplyr::lag(ave_abundance, timeframe[1])) %>%
      mutate(decline = ifelse(ave_abundance == 0, NA, decline)) %>% 
      # calculate the rate of change
      # mutate(decline = diff/dplyr::lag(abundance, timeframe[1])) %>% 
      # mutate(decline = diff/dplyr::lag(ave_abundance, 10)) %>%
      # mutate(prev = dplyr::lag(ave_abundance, 10)) %>% 
      # assign red list risk status based on decline 
      # Using the thresholds from p 16 Categories A2 - A4 Guidelines 
      # for Using the IUCN Red List Categories and Criteria v14 August 2019
      mutate(rl_status = ifelse(decline < 0.20, "LC",
                         ifelse(decline >= 0.20 & decline < 0.30, "NT", # Where did this and LC thresholds come from?
                         ifelse(decline >= 0.30 & decline < 0.50, "VU",
                         ifelse(decline >= 0.50 & decline < 0.80, "EN",
                         ifelse(decline >= 0.80, "CR",
                         ifelse(is.na(decline), "EX", "TBD"))))))) %>%
      arrange(group_id, monthly_time_step) %>%
      # Replace all non-ex status with ex after first occurrence 
      # mutate(extinct = match("EX", rl_status)) %>%
      mutate(rl_status = ifelse(ave_abundance == 0, "EX", rl_status)) %>% 
      mutate(extinct = ifelse(rl_status == "EX", 1, 0)) %>% 
      # mutate(rl_status = with(., ave(rl_status, 
      #                                         FUN=maintain_ex_status)))
      #mutate(rl_status = rl_status) %>% 
      group_by(group_id)
    
  }
  
  scenario_red_list_df <- do.call(rbind, group_red_list_data)
  
  scenario_red_list_data[[i]] <- scenario_red_list_df
  
  # Save the inputs
  
  # saveRDS(scenario_red_list_df,
  #         file.path(rli_inputs_folder,
  #                   paste(today, scenarios[[i]], "replicate", j,
  #                         "RLI_input_data_monthly_smoothing.rds", sep = "_")))
  # 
  # write.csv(replicate_red_list_df,
  #           file.path(rli_inputs_folder,
  #                     paste(today, scenarios[[i]], "replicate", j,
  #                           "RLI_input_data_monthly_smoothing.csv", sep = "_")))
  
}


# Have a quick look at the outputs

rli_inputs <- scenario_red_list_data[[1]]
tail(rli_inputs)
dim(rli_inputs)

write.csv(rli_inputs, file.path(indicator_outputs_folder, 
                                "rli_input_example_averaged_reps.csv"))

rli_inputs_group <- rli_inputs %>% filter(group_id == "13.16.27")

ggplot(data = rli_inputs_group) +
  geom_path(aes(x = annual_time_step, y = ave_abundance)) +
  theme(legend.position = "none") +
  geom_text(aes(x = annual_time_step, y = ave_abundance, label = rl_status))


# * Take coarser sample ----
## Heterotrophs

sample_interval <- 1 # interval between samples in years
sample_max_timestep <- 300/sample_interval

scenario_redlist_data_sampled <- list()

for (i in seq_along(scenario_red_list_data)) {
    
  # Sample
   sampled <- scenario_red_list_data[[i]] %>% 
     group_by(group_id) %>%
    slice(which(row_number() %% sample_interval == 0)) %>% 
    mutate(annual_time_step = seq(1,sample_max_timestep,1))
  
  # Remove groups that have no biomass in sampled years or you get weird outputs
   
   scenario_redlist_data_sampled[[i]] <- sampled %>% 
           group_by(group_id) %>% 
           mutate(total_abundance = sum(abundance, na.rm = TRUE),
                  exists = ifelse(total_abundance == 0 , FALSE, TRUE)) %>% 
           filter(exists == TRUE)
  
}

test <- scenario_redlist_data_sampled[[1]]
test_group <- test %>% filter(group_id == "13.16.27")
dim(test_group) # should have 300 rows

## Autotrophs

scenario_auto_sampled <- list()

for (i in seq_along(scenario_auto_annual)) {
  
  # Sample
  scenario_auto_sampled[[i]] <- scenario_auto_annual[[i]] %>% 
    group_by(group_id) %>% 
    slice(which(row_number() %% sample_interval == 0)) %>% 
    mutate(annual_time_step = seq(1,sample_max_timestep,1))
  
}

test <- scenario_auto_sampled[[1]]
test_group <- test %>% filter(group_id == "autotrophs")
dim(test_group) # should have 300 rows

# * Get harvested groups only ----

scenario_harvested_groups <- list()

for (i in seq_along(scenario_redlist_data_sampled)) {
  
  scenario <- scenarios[[i]]
  
  if (scenario == "000_Baseline") {
    
    # No disturbance so just get whatever groups
    
    harvested <- scenario_redlist_data_sampled[[i]] %>% 
                                      mutate(scenario = scenarios[[i]]) %>% 
                                      group_by(monthly_time_step) %>% 
                                      mutate(indicator_score = sum(abundance, na.rm = TRUE),
                                             indicator = "total abundance harvested",
                                             ci_lower = NA,
                                             ci_upper = NA,
                                             replicate = NA) %>% 
                                      ungroup(.) %>%                                 
                                      dplyr::select(annual_time_step, 
                                                    indicator_score,
                                                    ci_lower,
                                                    ci_upper,
                                                    indicator,
                                                    replicate,
                                                    scenario)%>% 
      distinct(.)
    
  } else if (scenario == "100_Land_use") {
    
    # Get the autotroph abundance
    
    harvested <- scenario_auto_sampled[[i]] %>% 
        filter(group_id == "autotrophs") %>% 
      mutate(scenario = scenarios[[i]]) %>%
      group_by(monthly_time_step) %>% 
      #mutate(ab_scaled = range01(abundance),
       mutate(indicator_score = sum(abundance, na.rm = TRUE),
             indicator = "total abundance harvested",
             ci_lower = NA,
             ci_upper = NA,
             replicate = NA) %>% 
      ungroup(.) %>% 
      dplyr::select(annual_time_step, 
                    indicator_score,
                    ci_lower,
                    ci_upper,
                    indicator,
                    replicate,
                    scenario) %>% 
      distinct(.)
    
  } else if (scenario == "200_Harvesting_carnivores") {
    
    # Get the carnivorous endotherms 100 - 200kg
    # Harvest doesn't seem to affect the ectotherms???
    
    harvested <- scenario_redlist_data_sampled[[i]] %>% 
      filter(functional_group_name == "carnivore endotherm" &
             group_id == "11.67"|
             functional_group_name == "carnivore endotherm" &
             group_id == "11.68")  %>% 
      mutate(scenario = scenarios[[i]]) %>% 
      ungroup(.) %>% 
      dplyr::select(annual_time_step, abundance, scenario) %>% 
      distinct(.) %>% 
      mutate(ab_scaled = range01(abundance)) %>% 
      group_by(annual_time_step) %>% 
      #mutate(indicator_score = sum(ab_scaled, na.rm = TRUE),
       mutate(indicator_score = sum(abundance, na.rm = TRUE),
              indicator = "total abundance harvested",
             ci_lower = NA,
             ci_upper = NA,
             replicate = NA) %>% 
      ungroup(.) %>% 
      dplyr::select(annual_time_step, 
                    indicator_score,
                    ci_lower,
                    ci_upper,
                    indicator,
                    replicate,
                    scenario) %>% 
      distinct(.)
     
    
  } else if (scenario == "300_Harvesting_herbivores") {
    
    # Get the herbivorous endotherms 100 - 200kg
    # Harvest doesn't seem to affect the ectotherms???
    
    harvested <- scenario_redlist_data_sampled[[i]] %>% 
      filter(functional_group_name == "herbivore endotherm" &
             group_id == "10.67"|
             functional_group_name == "herbivore endotherm" &
             group_id == "10.68")   %>% 
      mutate(scenario = scenarios[[i]]) %>% 
      ungroup(.) %>% 
      dplyr::select(annual_time_step, abundance, scenario) %>% 
      distinct(.) %>% 
      #mutate(ab_scaled = range01(abundance)) %>% 
      group_by(annual_time_step) %>% 
      mutate(indicator_score = sum(abundance, na.rm = TRUE),
             indicator = "total abundance harvested",
             ci_lower = NA,
             ci_upper = NA,
             replicate = NA) %>% 
      ungroup(.) %>% 
      dplyr::select(annual_time_step, 
                    indicator_score,
                    ci_lower,
                    ci_upper,
                    indicator,
                    replicate,
                    scenario) %>% 
      distinct(.)

  } 
  
  scenario_harvested_groups[[i]] <- harvested
}

data <- scenario_harvested_groups[[2]]
head(data)
dim(data)

# * Plot harvested groups ----

harvested_plots_folder <- file.path(indicator_plots_folder, 
                                    "harvested_plots", 
                                    today)

if( !dir.exists( file.path(harvested_plots_folder) ) ) {
  dir.create( file.path(harvested_plots_folder), recursive = TRUE )
  
}

harvested_plots <- list()

for ( i in seq_along(scenario_harvested_groups)) {
  
harvested_plots[[i]] <- ggplot(data = scenario_harvested_groups[[i]]) +
  geom_line(aes(x = annual_time_step, 
                y = log(indicator_score))) +
  theme(legend.position = "bottom") +
  labs(title = paste(scenarios[[i]], "harvested groups", sep = " "))
  
ggsave(file.path(harvested_plots_folder, paste(today, scenarios[[i]],
                                         "harvested_reps_averaged.png",
                                         sep = "_")),
       harvested_plots[[i]],  device = "png")

}

harvested_plots[[1]]
harvested_plots[[2]]
harvested_plots[[3]]
harvested_plots[[4]]

# * Get abundance of functional groups ----

scenario_fg_abundance <- list()

for (i in seq_along(scenario_redlist_data_sampled)) {
  
  scenario_fg_abundance[[i]] <- scenario_redlist_data_sampled[[i]] %>% 
          mutate(scenario = scenarios[[i]]) %>% 
          ungroup(.) %>% 
          dplyr::select(annual_time_step, functional_group_name, 
                        abundance, scenario) %>% 
          distinct(.) %>% 
          mutate(ab_scaled = range01(abundance)) %>% 
          group_by(annual_time_step, functional_group_name) %>% 
          mutate(indicator_score = mean(ab_scaled, na.rm = TRUE),
                 indicator = paste(functional_group_name, "mean abundance", sep = " "),
                 ci_lower = NA,
                 ci_upper = NA,
                 replicate = NA) %>% 
          ungroup(.) %>% 
          dplyr::select(annual_time_step, 
                        indicator_score,
                        ci_lower,
                        ci_upper,
                        indicator,
                        replicate,
                        scenario) %>% 
          distinct(.)
    
}

fg <- scenario_fg_abundance[[2]]
head(fg)

data <- scenario_fg_abundance[[4]]
head(data)

ggplot(data = data) +
  geom_line(aes(x = annual_time_step, 
                y = log(indicator_score),
                col = indicator)) 

# Plot abundance of functional groups

scenario_fg_abundance_plots <- list()

for ( i in seq_along(scenario_fg_abundance)) {
  
  groups <- split(scenario_fg_abundance[[i]], 
                  scenario_fg_abundance[[i]]$indicator)
  
  group_plots <- list()
  
  for ( j in seq_along(groups)) {
    
    group_plots[[j]] <- ggplot(data = groups[[j]]) +
    geom_line(aes(x = annual_time_step, 
                  y = indicator_score)) +
    theme(legend.position = "bottom") +
    labs(title = paste(scenarios[[i]], "functional groups", sep = " ")) +
    facet_wrap(~ indicator)
  
  ggsave(file.path(harvested_plots_folder, paste(today, scenarios[[i]],
                                                 groups[[j]]$indicator,
                                                 "functional_group_abundance.png",
                                                 sep = "_")),
         harvested_plots[[i]],  device = "png")
  }
  
  scenario_fg_abundance_plots[[i]] <- group_plots
  
}

scenario_fg_abundance_plots[[1]]

i <- 1
scenario_fg_abundance_plots[[2]][[i]]
i <- i + 1
scenario_fg_abundance_plots[[2]][[i]]
i <- 1
scenario_fg_abundance_plots[[3]][[i]]
i <- i + 1
scenario_fg_abundance_plots[[3]][[i]]
i <- 1
scenario_fg_abundance_plots[[4]][[i]]
i <- i + 1
scenario_fg_abundance_plots[[4]][[i]]


head(scenario_redlist_data_sampled)

# Plot the individual groups in each functional group

scenario_biomass_plots <- list()

for ( i in seq_along(scenario_redlist_data_sampled)) {

data <- scenario_redlist_data_sampled[[i]] %>% 
        filter(bodymass_index > 60)

fg_data <- split(data, data$functional_group_name)

fg_biomass_plots <- list()

for (j in seq_along (fg_data)) {
  
  fg <- fg_data[[j]]$functional_group_name[1]
  
fg_biomass_plots[[j]] <- ggplot(data = fg_data[[j]]) +
  geom_smooth(aes(x = annual_time_step, 
                y = log(abundance),
                col = group_id)) +
  geom_text(data = subset(fg_data[[j]], annual_time_step == 300),
              aes(x = annual_time_step,
                y = log(abundance),
                label = group_id, colour = group_id)) +
  theme(legend.position = "none") +
  #facet_wrap( ~ functional_group_name) +
  labs(title = paste(fg, scenarios[[i]]), sep = " ")

ggsave(file.path(harvested_plots_folder, paste(today, scenarios[[i]], fg,
                                         "log_abundance_all_species.png",
                                         sep = "_")),
       fg_biomass_plots[[j]],  device = "png")

  }

scenario_biomass_plots[[i]] <- fg_biomass_plots

}

scenario_biomass_plots[[1]]
scenario_biomass_plots[[2]]
scenario_biomass_plots[[3]]
scenario_biomass_plots[[4]]
            
# * Get proportion extinct ----

scenario_extinctions <- list()
scenario_rl_status_plots <- list()

for (i in seq_along(scenario_redlist_data_sampled)) {
  
  scenario_extinctions[[i]] <- scenario_redlist_data_sampled[[i]] %>% 
          group_by(annual_time_step, rl_status) %>% 
          filter(rl_status != is.na(rl_status)) %>% 
          summarise(test = n())
  
  scenario_rl_status_plots[[i]] <- ggplot(scenario_extinctions[[i]], 
                                          aes(x = annual_time_step, 
                                              y = test, 
                                              fill = rl_status)) +
    geom_bar(position = "stack", stat = "identity") +
    scale_fill_viridis_d()
  
  ggsave(file.path(rli_plots_folder, paste(today, scenarios[[i]], 
                                           "red_list_status.png",
                                                 sep = "_")),
         scenario_rl_status_plots[[i]],  device = "png")
  
}

scenario_rl_status_plots[[1]]
scenario_rl_status_plots[[2]]
scenario_rl_status_plots[[3]]
scenario_rl_status_plots[[4]]

# * Calculate RLI ----

# RLI by individual functional groups

scenario_fg_rli_outputs <- list()

for (i in seq_along(scenario_redlist_data_sampled)) {
  
  scenario_fg_rli_outputs[[i]] <- calculate_red_list_index(
    scenario_redlist_data_sampled[[i]], numboots, ci = FALSE) %>%
    mutate(scenario = scenarios[[i]]) 
  
}

x <- scenario_fg_rli_outputs[[1]]
head(x)

# Mean RLI aggregated across groups

scenario_rli_outputs <- list()

for (i in seq_along(scenario_fg_rli_outputs)) {
  
  if ("ci_lower" %in% names(scenario_fg_rli_outputs[[i]])) {
    
    scenario_rli_outputs[[i]] <- scenario_fg_rli_outputs[[i]] %>%
      group_by(annual_time_step) %>%
      summarise(indicator_score = mean(indicator_score),
                ci_lower = mean(ci_lower),
                ci_upper = mean(ci_upper)) %>%
      mutate(indicator = "RLI",
             replicate = j)
  } else {
    
    scenario_rli_outputs[[i]] <- scenario_fg_rli_outputs[[i]] %>%
      group_by(annual_time_step) %>%
      summarise(indicator_score = mean(indicator_score)) %>%
      mutate(indicator = "RLI",
             replicate = j)
  }
  
}

head(scenario_rli_outputs)[[1]]

# * Plot RLI ----

## By functional group

scenario_fg_rli_plots <- list()

for (i in seq_along(scenario_fg_rli_outputs)) {
  
  scenario_fg_rli_plots[[i]] <-  plot_red_list_index_by_group(
    scenario_fg_rli_outputs[[i]],
    impact_start,
    impact_end,
    ci = FALSE)
  
  ggsave(file.path(rli_plots_folder, paste(today, scenarios[[i]], 
                                           "RLI_by_functional_group_reps_averaged.png",
                                           sep = "_")),
         scenario_fg_rli_plots[[i]],  device = "png")
  
}

scenario_fg_rli_plots[[4]]


# RLI with all functional groups aggregated
# i.e. mean of each 'taxa' RLI as per Butchart et al (2010) 'Indicators of
# recent declines'

scenario_rli_plots <- list()

for (i in seq_along(scenario_rli_outputs)) {
  
  
  scenario_rli_plots[[i]] <- plot_red_list_index(scenario_rli_outputs[[i]],
                                                 impact_start, 
                                                 impact_end,
                                                 ci = TRUE)
  
  ggsave(file.path(rli_plots_folder, paste(today, scenarios[[i]], 
                                           "RLI_aggregated_averaged.png",
                                           sep = "_")),
         scenario_rli_plots[[i]],  device = "png")
  
}

scenario_rli_plots[[1]]
scenario_rli_plots[[2]]
scenario_rli_plots[[3]]
scenario_rli_plots[[4]]

# LIVING PLANET INDEX ----

# * Create folders ----

lpi_inputs_folder <- file.path(indicator_inputs_folder, "LPI_inputs", today)

if( !dir.exists( file.path(lpi_inputs_folder) ) ) {
  dir.create( file.path(lpi_inputs_folder), recursive = TRUE )
  
}

lpi_outputs_folder <- file.path(indicator_outputs_folder, "LPI_outputs", today)

if( !dir.exists( file.path(lpi_outputs_folder) ) ) {
  dir.create( file.path(lpi_outputs_folder), recursive = TRUE )
  
}

lpi_plots_folder <- file.path(indicator_plots_folder, "LPI_plots", today)

if( !dir.exists( file.path(lpi_plots_folder) ) ) {
  dir.create( file.path(lpi_plots_folder), recursive = TRUE )
  
}

# * Subset data ----

scenario_lpi_inputs <- list()

for (i in seq_along(scenario_redlist_data_sampled)) {
  
  scenario_lpi_inputs[[i]] <- scenario_redlist_data_sampled[[i]] %>% 
    dplyr::select(group_id, annual_time_step, 
                  ave_abundance, functional_group_name) # %>% 
    # note - the carnivorous endotherms cause the time series to look like they're declining even before impact
    # filter(functional_group_name != "carnivore endotherm")
  
  
}

# lpi_input <- scenario_lpi_inputs[[1]]
# head(lpi_input)
# write.csv(lpi_input, file.path(indicator_outputs_folder, "lpi_input_example_annual.csv"))

# * Calculate LPI ----

# Retain naming convention, the LPI just takes the abundance dataframes we
# already formatted while making the RLI inputs

# scenario_lpi_inputs <- scenario_abundance_long

# Loop through each scenario and replicate and calculate the LPI per rep

# ** All groups ----

scenario_lpi_outputs <- list()

for (i in seq_along(scenario_lpi_inputs)) {
  
  scenario_lpi_outputs[[i]] <- calculate_living_planet_index(
    
    scenario_lpi_inputs[[i]], start_time_step, ci = FALSE, numboots, j
    
  ) 
  
}

head(scenario_lpi_outputs)[[1]]

# * Plot LPI replicates individually ----

scenario_lpi_plots <- list()

for (i in seq_along(scenario_lpi_outputs)) {
  
  scenario_lpi_plots[[i]] <- plot_living_planet_index(scenario_lpi_outputs[[i]],
                                                      ci = FALSE)
  
  
  ggsave(file.path(lpi_plots_folder, paste(today, scenarios[[i]], 
                                           "LPI_aggregated_reps_averaged.png",
                                           sep = "_")),
         scenario_lpi_plots[[i]],  device = "png")                                   
  
}

scenario_lpi_plots[[1]]
scenario_lpi_plots[[2]]
scenario_lpi_plots[[3]]
scenario_lpi_plots[[4]]

# ** Functional groups ----

scenario_lpi_fg_outputs <- list()

for (i in seq_along(scenario_lpi_inputs)) {
  
  fg_data <- split(scenario_lpi_inputs[[i]], 
                   scenario_lpi_inputs[[i]]$functional_group_name)

  fg_lpi_outputs <- list()
  
  for ( j in seq_along(fg_data)) {
  
    fg_lpi_outputs[[j]] <- calculate_living_planet_index(
      
      fg_data[[j]], start_time_step, ci = FALSE, numboots, 
      fg_data[[j]]$functional_group_name[1])
  }
  
  scenario_lpi_fg_outputs[[i]] <- do.call(rbind,fg_lpi_outputs)
  
}

head(scenario_lpi_fg_outputs[[1]])

# Plot fg

scenario_fg_lpi_plots <- list()

for (i in seq_along(scenario_lpi_fg_outputs)) {
  
  scenario_fg_lpi_plots[[i]] <- ggplot(data = scenario_lpi_fg_outputs[[i]],
                                       aes(x = annual_time_step, 
                                           y = indicator_score,
                                           col = replicate)) +
                                geom_line() +
                                facet_wrap( ~ replicate)
  
  ggsave(file.path(lpi_plots_folder, paste(today, scenarios[[i]], 
                                           "LPI_func_groups_reps_averaged.png",
                                           sep = "_")),
         scenario_fg_lpi_plots[[i]],  device = "png")   
  
}

scenario_fg_lpi_plots[[1]]
scenario_fg_lpi_plots[[2]]
scenario_fg_lpi_plots[[3]]
scenario_fg_lpi_plots[[4]]


# ** Trophic groups ----

scenario_lpi_tg_outputs <- list()

for (i in seq_along(scenario_lpi_inputs)) {
  
  data <- scenario_lpi_inputs[[i]] %>% 
          mutate(trophic_group = word(functional_group_name, 1))
  
  head(data)
  
  tg_data <- split(data, data$trophic_group)
  
  tg_lpi_outputs <- list()
  
  for ( j in seq_along(tg_data)) {
    
    tg_lpi_outputs[[j]] <- calculate_living_planet_index(
      
      tg_data[[j]], start_time_step, ci = FALSE, numboots, 
      tg_data[[j]]$trophic_group[1])
  }
  
  scenario_lpi_tg_outputs[[i]] <- do.call(rbind,tg_lpi_outputs)
  
}

head(scenario_lpi_tg_outputs[[1]])

# Plot trophic groups

scenario_tg_lpi_plots <- list()

for (i in seq_along(scenario_lpi_tg_outputs)) {
  
  scenario_tg_lpi_plots[[i]] <- ggplot(data = scenario_lpi_tg_outputs[[i]],
                                       aes(x = annual_time_step, 
                                           y = indicator_score,
                                           col = replicate)) +
    geom_line() +
    facet_wrap( ~ replicate)
  
  ggsave(file.path(lpi_plots_folder, paste(today, scenarios[[i]], 
                                           "LPI_troph_groups_reps_averaged.png",
                                           sep = "_")),
         scenario_tg_lpi_plots[[i]],  device = "png")   
  
}

scenario_tg_lpi_plots[[1]]
scenario_tg_lpi_plots[[2]]
scenario_tg_lpi_plots[[3]]
scenario_tg_lpi_plots[[4]]

# Combine indicators ----

## Standardise outputs

for ( i in seq_along(scenario_rli_outputs)) {
  
  scenario_rli_outputs[[i]] <- scenario_rli_outputs[[i]] %>% 
                               mutate(scenario = scenarios[[i]]) %>%
    ungroup(.) %>% 
    dplyr::select(annual_time_step, indicator_score, ci_lower, ci_upper,
                  indicator, scenario)
}

head(scenario_rli_outputs[[2]])

for ( i in seq_along(scenario_fg_rli_outputs)) {
  
  scenario_fg_rli_outputs[[i]] <- scenario_fg_rli_outputs[[i]] %>% 
    mutate(scenario = scenarios[[i]],
           indicator = paste(functional_group_name, "RLI", sep = " ")) %>%
    ungroup(.) %>% 
    dplyr::select(annual_time_step, indicator_score, ci_lower, ci_upper,
                  indicator, scenario)
}

head(scenario_fg_rli_outputs[[2]])

for ( i in seq_along(scenario_lpi_outputs)) {
  
  scenario_lpi_outputs[[i]] <- scenario_lpi_outputs[[i]] %>% 
    mutate(scenario = scenarios[[i]],
           indicator = "LPI") %>% 
    dplyr::select(annual_time_step, indicator_score, ci_lower, ci_upper,
                  indicator, scenario)
}

head(scenario_lpi_outputs[[2]])

for ( i in seq_along(scenario_lpi_fg_outputs)) {
  
  scenario_lpi_fg_outputs[[i]] <- scenario_lpi_fg_outputs[[i]] %>% 
    mutate(scenario = scenarios[[i]],
           indicator = paste(replicate, "LPI", sep = " ")) %>% 
    dplyr::select(annual_time_step, indicator_score, ci_lower, ci_upper,
                  indicator, scenario)
}

head(scenario_lpi_fg_outputs[[2]])

for ( i in seq_along(scenario_lpi_tg_outputs)) {
  
  scenario_lpi_tg_outputs[[i]] <- scenario_lpi_tg_outputs[[i]] %>% 
    mutate(scenario = scenarios[[i]],
           indicator = paste(replicate, "LPI", sep = " ")) %>% 
    dplyr::select(annual_time_step, indicator_score, ci_lower, ci_upper,
                  indicator, scenario)
}

head(scenario_lpi_tg_outputs[[2]])

head(scenario_fg_abundance[[2]])

for ( i in seq_along(scenario_fg_abundance)) {
  
  scenario_fg_abundance[[i]] <- scenario_fg_abundance[[i]] %>% 
        dplyr::select(annual_time_step, indicator_score, ci_lower, ci_upper,
                  indicator, scenario)
}

head(scenario_fg_abundance[[2]])

head(scenario_harvested_groups[[2]])

for ( i in seq_along(scenario_harvested_groups)) {
  
  scenario_harvested_groups[[i]] <- scenario_harvested_groups[[i]] %>% 
    dplyr::select(annual_time_step, indicator_score, ci_lower, ci_upper,
                  indicator, scenario)
}

head(scenario_harvested_groups[[2]])


names(scenario_rli_outputs[[1]]) == names(scenario_fg_rli_outputs[[1]])

names(scenario_fg_rli_outputs[[1]]) == names(scenario_lpi_outputs[[1]])

names(scenario_lpi_outputs[[1]]) == names(scenario_lpi_fg_outputs[[1]])

names(scenario_lpi_fg_outputs[[1]]) == names(scenario_lpi_tg_outputs[[1]])

names(scenario_lpi_tg_outputs[[1]]) == names(scenario_fg_abundance[[1]])

names(scenario_fg_abundance[[1]]) == names(scenario_harvested_groups[[1]])

test <- split(scenario_fg_rli_outputs, scenario_fg_rli_outputs$indicator)

# Save outputs

all_indicators_list <- list(scenario_rli_outputs,
                            scenario_fg_rli_outputs,
                            scenario_lpi_outputs,
                            scenario_lpi_fg_outputs,
                            scenario_lpi_tg_outputs,
                            scenario_fg_abundance,
                            scenario_harvested_groups)

names(all_indicators_list) <- c("RLI",
                                "RLI functional groups",
                                "LPI",
                                "LPI functional groups",
                                "LPI trophic groups",
                                "abundance functional groups",
                                "abundance harvested groups")

saveRDS(all_indicators_list,
        file.path(indicator_outputs_folder,
                  paste(today, "all_indicators_output_data_list_reps_averaged.rds",
                        sep = "_")))

all_indicators_all_scenarios <- flatten(all_indicators_list)
all_indicators <- do.call(rbind, all_indicators_all_scenarios)

saveRDS(all_indicators,
        file.path(indicator_outputs_folder,
                  paste(today, "all_indicators_output_data_reps_averaged.rds",
                        sep = "_")))

write.csv(all_indicators,
          file.path(indicator_outputs_folder,
                    paste(today, "all_indicators_output_data_reps_averaged.csv",
                          sep = "_")))

# Reshuffle the list so the top level of the list is indicator, then scenario

new_split <- split(all_indicators, all_indicators$indicator)

new_indicator_list <- list()
indicator_names <- list()

for (i in seq_along(new_split)) {
  
  indicator_data <- new_split[[i]]
  
  new_indicator_list[[i]] <- split(indicator_data, indicator$scenario)
  
  indicator_names[i] <- indicator_data$indicator[1]

}

names(new_indicator_list) <- indicator_names

saveRDS(new_indicator_list,
        file.path(indicator_outputs_folder,
                  paste(today, "all_indicators_output_data_reps_averaged_list2.rds",
                        sep = "_")))
