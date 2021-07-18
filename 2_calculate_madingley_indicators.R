
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
library(tidylog)

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
    scale_fill_viridis_d() +
    scale_color_viridis_d() + 
    facet_wrap(~functional_group_name) +
    labs(x = "Time", 
         y = "Red List Index Score") +
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
      geom_line(aes(colour = functional_group_name)) +
      scale_fill_viridis_d() +
      scale_color_viridis_d() + 
      facet_wrap(~functional_group_name) +
      labs(x = "Time", 
           y = "Red List Index Score") +
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



calculate_living_planet_index <- function(data, start_time_step = 1, ci = FALSE,
                                          numboots, replicate_num = NA){
  
  filtered_inputs <- data %>%
    # Remove timesteps if needed (if not make start_time_step = 1)
    filter(annual_time_step >= start_time_step) %>%
    # Group by virtual species (vs)
    group_by(group_id) %>% 
    # Remove vs that have no individuals at any timestep
    filter(!all(abundance == 0)) %>% 
    # Save a copy of original abundance values
    rename(abundance_original = abundance) %>%
    # Identify rows where abundance == 0, and change all subsequent years to 0 too
    # mutate(abundance = maintain_0_abundance(abundance_original)) 
    mutate(abundance = abundance_original)
  
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
    geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper),
                    alpha = 0.4)
  } else {
    
    ggplot(data, aes(x = annual_time_step, y = indicator_score)) +
      geom_line()
    
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

# Get date to label outputs

today <- Sys.Date()

# Define mode (development == TRUE will only perform operations on a small subset
# of folders, not all outputs)

development_mode <- TRUE

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

# Check we've got the correct number of files for each (should match)

if (length(abundance_files) != length(generation_files)) {
  
  stop("Number of abundance files and generation files do not match")

  }

# Read in the data
# Note, generation lengths are different sizes bc each replicate has slightly different number of groups present
scenario_generations_raw[[i]] <- lapply(generation_files, readRDS) 
scenario_abundance_raw[[i]] <- lapply(abundance_files, readRDS)

}

# x <- scenario_abundance_raw[[1]][[1]]

# Remove burn-in ----

# Remove the burn in timesteps for abundance files (output should be same structure
# as input, a nested list of x3 scenarios, each with x25 replicate dataframes,
# but each dataframe will have fewer columns, no NA data and be absolute not log
# abundance densities)

scenario_abundance_formatted <- list()
replicate_abundance_formatted <- list()

for (i in seq_along(scenario_abundance_raw)) {
  
  # Get all replicates for one scenario
  abundance_reps <- scenario_abundance_raw[[i]]
  
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

rm(abundance_reps, abundance_temp, replicate_abundance_formatted, rep)

# Check structure is still correct

length(scenario_abundance_formatted) == length(scenario_abundance_raw)

# Remove previous version to make space

rm(scenario_abundance_raw)

# Convert to long format ----

## NOTE: Now the burnin period has been removed, what was previously monthly timestep
## 1200 is now monthly time step 1.

## Pivot longer (group_id, timestep, abundance)

scenario_abundance_long <- list()
replicate_abundance_long <- list()

for (i in seq_along(scenario_abundance_formatted)) {
  
  scenario_abundance <- scenario_abundance_formatted[[i]]
  
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

#' # Mean generation length ----
#' 
#' #' TODO: Check this function doesn't exclude groups not present in all 
#' #' replicates
#' 
#' mean_gen_lengths_scenarios <- list()
#' 
#' for (i in seq_along(generation)) {
#'   
#' 
#'   # Select only group id and gen length, reduce rows (input has a
#'   # new row for each cohort) for one scenario
#'   
#'    scenario_gen_lengths <- lapply(generation[[i]], function(x){
#'        
#'     x %>%
#'     dplyr::select(group_id, generation_length_yrs) %>%
#'       unique(.)}
#'     )
#'   
#'   # Convert individual replicates into a single dataframe where each rep is 
#'   # a column (suppress warning about duplicate column names bc we fix it 
#'   # in next step)
#'    
#'   mean_gen_lengths <- suppressWarnings(Reduce(function(x, y) merge(x, y, by = "group_id"), 
#'                              scenario_gen_lengths))
#'   
#'   # Fix up column names for each replicate so they're unique
#'   all_reps <- colnames(mean_gen_lengths)
#'   all_reps <- all_reps[!grepl('group_id', all_reps)]
#'   all_reps_unique <- make.unique(all_reps)
#'   names(mean_gen_lengths) <- c("group_id", all_reps_unique)
#'   head(mean_gen_lengths)
#'   
#'   # Calculate the mean generation length for each group across replicates
#'   
#'   mean_gen_lengths_scenarios[[i]] <- mean_gen_lengths %>%
#'                       mutate(mean_gen_length = 
#'                                rowMeans(mean_gen_lengths[all_reps_unique])) %>%
#'                       dplyr::select(group_id, mean_gen_length)
#' }
#' 
#' rm(mean_gen_lengths, scenario_gen_lengths)

# RED LIST INDEX ----

# * Create folders ----

rli_inputs_folder <- file.path(indicator_inputs_folder, "RLI_inputs")

if( !dir.exists( file.path(rli_inputs_folder) ) ) {
  dir.create( file.path(rli_inputs_folder), recursive = TRUE )
  
}

rli_outputs_folder <- file.path(indicator_outputs_folder, "RLI_outputs")

if( !dir.exists( file.path(rli_outputs_folder) ) ) {
  dir.create( file.path(rli_outputs_folder), recursive = TRUE )
  
}

rli_plots_folder <- file.path(indicator_plots_folder, "RLI_plots")

if( !dir.exists( file.path(rli_plots_folder) ) ) {
  dir.create( file.path(rli_plots_folder), recursive = TRUE )
  
}

## Referring to the thresholds quote under Criterion A, Reason 1 (declines
## are the result of reversible pressures) according to:
## https://portals.iucn.org/library/sites/library/files/documents/RL-2001-001-2nd.pdf

# * Merge abundance and generation length data ----

rm(formatted)

scenario_ab_gl_formatted <- list()
replicate_ab_gl_formatted <- list()

for (i in seq_along(scenario_abundance_long)) {
  
  replicate_abundance <- scenario_abundance_long[[i]]
  replicate_generations <- scenario_generations_raw[[i]]
  
  # For each individual replicate
  
  for (j in seq_along(replicate_abundance)) {
    
    # Reduce size of the replicate generations dataframe or the merge won't work
    gen_length <- replicate_generations[[j]] %>% 
                  dplyr::select(group_id, generation_length_yrs, functional_group_name) %>% 
                  distinct(.)
    
    # Add the generation length info to the abundance dataframe
    temp1 <- replicate_abundance[[j]] %>%
        merge(gen_length, by = "group_id") %>%
        arrange(monthly_time_step, group_id) %>%
    # Important - following lines assume an annual timeframe, will need to adjust if change interval
    mutate(generation_by_three = generation_length_yrs * 3) %>% # Time over which to measure decline, 3 x gen length OR:
    mutate(timeframe = ifelse(generation_by_three > gen_timeframe, # 10 years 
                       round(generation_by_three), gen_timeframe)) %>%
    dplyr::select(-generation_by_three) %>%
    distinct(.) %>%
    group_by(group_id) %>% 
      # select rows that are multiples of the specified interval 
      # (eg if interval is 12, it samples one month from every year)
    slice(which(row_number() %% interval == 0)) %>% 
    mutate(annual_time_step = seq(1,max_timestep,1)) # %>% 
    
    # Find the last time step where non-0 abundance occurred for each group
    
    temp2 <- temp1 %>% 
            group_by(group_id) %>% 
            filter(abundance > 0) %>% 
            dplyr::select(group_id, annual_time_step, abundance) %>% 
            filter(annual_time_step == max(annual_time_step)) %>% 
            dplyr::select(group_id, annual_time_step) %>% 
            rename(last_abundance = annual_time_step)

    # Add the year of last positive abundance number as a column to the data    
    temp3 <- temp1 %>% 
           merge(temp2, by = c("group_id"), all = TRUE)
    
    
    replicate_ab_gl_formatted[[j]] <- temp3 %>% 
               group_by(group_id) %>% 
               mutate(true_extinction = ifelse(abundance == 0 & 
                                               annual_time_step < last_abundance,
                      "false extinction",
                      ifelse(abundance > 0 & 
                               annual_time_step < last_abundance,
                             "not extinct",
                             ifelse(abundance == 0 & 
                                    annual_time_step >= last_abundance,
                                    "true extinction", "not extinct")))) %>% 
               filter(true_extinction != "false extinction") %>% 
               group_by(group_id) %>% 
               arrange(annual_time_step)
    

    # if (length(formatted_replicate) == 0) {
    #   
    #   rm(formatted_replicate)
    #   
    # } else {
    
    # add formatted replicate here
      
    # Check if there are any carnivorous endotherms
    
    check <- replicate_ab_gl_formatted[[j]] %>% 
             group_by(functional_group_name) %>% 
             summarise(present = sum(abundance)) %>% 
             filter(functional_group_name == "carnivore endotherm") %>% 
             dplyr::select(present) %>% 
             pull(.)
    
    # if they are present, keep the replicate
    
    if (length(check != 0)) {
      
      scenario_ab_gl_formatted[[i]] <- replicate_ab_gl_formatted
      
    } else {
      
    # if they are not, remove that replicate
      
      rm(replicate_ab_gl_formatted)
      
      print(paste("Replicate", i, "removed because no carnivorous endotherms are present", sep = " "))
      
      
    }
  }
}  

head(scenario_ab_gl_formatted[[1]][[1]])
dim(scenario_ab_gl_formatted[[1]][[1]])
length(unique(scenario_ab_gl_formatted[[1]][[1]]$group_id))

formatted <- scenario_ab_gl_formatted[[1]][[3]]

# * Sample data ----

# scenario_red_list_inputs_annual <- list()
# replicate_red_list_inputs_annual <- list()
# 
# for (i in seq_along(scenario_red_list_data)) {
#   
#   # Get replicate data for a single scenario
#   
#   replicate_red_list_data <- scenario_red_list_data[[i]]
#   
#   # For each individual replicate
#   
#   for (j in seq_along(replicate_red_list_data)) {
#     
#     replicate_red_list_inputs_annual[[j]] <-  replicate_red_list_data[[j]] %>%
#       group_by(group_id) %>% 
#       # select rows that are multiples of the specified interval 
#       # (eg if interval is 12, it samples one month from every year)
#       slice(which(row_number() %% interval == 0)) %>% 
#       mutate(annual_time_step = seq(1,max_timestep,1)) #see if it works just hard coding in number of time steps
#     
#   }
#   
#   scenario_red_list_inputs_annual[[i]] <- replicate_red_list_inputs_annual
#   
# }
# 
# y <- scenario_red_list_data[[1]][[1]]
# x <- scenario_red_list_inputs_annual[[1]][[1]]

# * Assign Red List Categories ----

scenario_red_list_data <- list()
replicate_red_list_data <- list()

for (i in seq_along(scenario_ab_gl_formatted)) {
  
  # Get replicate data for a single scenario
  
  replicate_ab_gl <- scenario_ab_gl_formatted[[i]]
  
  # For each individual replicate
  
  for (j in seq_along(replicate_ab_gl)) {
  
  # Split by functional group, because we calculate RLI for different
  # functional groups then aggregate later (as per Butchart etal 2010),
  # except we are using functional groups as proxies for taxa (eg mammals, birds, 
  # reptiles) used in real world RLI calcs
  
  status_inputs <- split(replicate_ab_gl[[j]], 
                         replicate_ab_gl[[j]]$group_id)
  
  # Make a list to hold output for each individual massbin-func-group (ie virtual spp)
  
  group_red_list_data <- list()
  
  
  for (k in seq_along(status_inputs)) {
    
    group_red_list_data[[k]] <- status_inputs[[k]] %>%
      group_by(group_id) %>%
      arrange(monthly_time_step) %>%
      # calculate the difference in abundance over 10 yrs or 3 generation lengths
      # (specified by 'timeframe' column). Its okay to take the first value of 
      # timeframe bc the dataframe is grouped by group_id, and timeframe only changes
      # between and not within group_ids
      mutate(diff = (abundance - dplyr::lag(abundance, timeframe[1]))) %>%
      # calculate the rate of change
      mutate(decline = diff/dplyr::lag(abundance, timeframe[1])) %>% 
      # assign red list risk status based on decline 
      mutate(rl_status = ifelse(decline > -0.40, "LC",
                                ifelse(decline <= -0.40 & decline > -0.50, "NT", # Where did this and LC thresholds come from?
                                       ifelse(decline <= -0.50 & decline > -0.70, "VU",
                                              ifelse(decline <= -0.70 & decline > -0.90, "EN",
                                                     ifelse(decline <= -0.90 & decline > -1, "CR",
                                                            ifelse(decline <= -1, "EX", "NA"))))))) %>%
      arrange(group_id, monthly_time_step) %>%
      # Replace all non-ex status with ex after first occurrence 
      # mutate(extinct = match("EX", rl_status)) %>%
      mutate(extinct = ifelse(rl_status == "EX", 1, 0)) %>% 
      # mutate(rl_status = with(., ave(rl_status, 
      #                                         FUN=maintain_ex_status)))
      #mutate(rl_status = rl_status) %>% 
      group_by(group_id)
    
    # group_red_list_data[[k]] <- status_inputs[[k]] %>%
    # group_by(group_id) %>%
    # arrange(monthly_time_step) %>%
    # # calculate the difference in abundance over 10 yrs or 3 generation lengths
    # # (specified by 'timeframe' column). Its okay to take the first value of 
    # # timeframe bc the dataframe is grouped by group_id, and timeframe only changes
    # # between and not within group_ids
    # mutate(diff = (abundance_adjusted - dplyr::lag(abundance_adjusted, timeframe[1]))) %>%
    # # calculate the rate of change
    # mutate(decline = diff/dplyr::lag(abundance_adjusted, timeframe[1])) %>% 
    # # assign red list risk status based on decline 
    # mutate(rl_status = ifelse(decline > -0.40, "LC",
    #                           ifelse(decline <= -0.40 & decline > -0.50, "NT", # Where did this and LC thresholds come from?
    #                           ifelse(decline <= -0.50 & decline > -0.70, "VU",
    #                           ifelse(decline <= -0.70 & decline > -0.90, "EN",
    #                           ifelse(decline <= -0.90 & decline > -1, "CR",
    #                           ifelse(decline <= -1, "EX", "NA"))))))) %>%
    #              arrange(group_id, monthly_time_step) %>%
    #  # Replace all non-ex status with ex after first occurrence 
    #  # mutate(extinct = match("EX", rl_status)) %>%
    #  mutate(extinct = ifelse(rl_status == "EX", 1, 0)) %>% 
    #  # mutate(rl_status = with(., ave(rl_status, 
    #  #                                         FUN=maintain_ex_status)))
    #  #mutate(rl_status = rl_status) %>% 
    #  group_by(group_id) #%>% 
     # # Check whether red list status is temporary or if it repeats
     # mutate(goes_extinct = ifelse(annual_time_step == max(annual_time_step) &
     #                              rl_status == "EX", TRUE, FALSE),
     #        previous_status = lag(rl_status, 1),
     #        next_status = lead(rl_status, 1),
     #        true_extinction = ifelse(goes_extinct == TRUE &
     #                                 rl_status == "EX" &
     #                                 previous_status == "EX",
     #                                 1,
     #                                 ifelse(goes_extinct == FALSE &
     #                                        rl_status == "EX",
     #                                 0, 2)),
     #        rl_status = ifelse(true_extinction == 0, NA , rl_status))
  
    # Print a message if extinctions have been adjusted
    
    if(!any(is.na(group_red_list_data[[k]]$extinct))) {
      
      print(paste("Note, non-extinct status replaced with extinct status for group", 
                  group_red_list_data[[k]]$group_id[1], sep = " "))
    }
   
  }
  
   replicate_red_list_df <- do.call(rbind, group_red_list_data)
   
   replicate_red_list_data[[j]] <- replicate_red_list_df
   
   # Save the inputs
   
   saveRDS(replicate_red_list_df,
           file.path(rli_inputs_folder,
                     paste(today, scenarios, "replicate", j,
                           "RLI_input_data.rds", sep = "_")))

   write.csv(replicate_red_list_df,
           file.path(rli_inputs_folder,
                     paste(today, scenarios, "replicate", j,
                           "RLI_input_data.csv", sep = "_")))
  
  
  }

  scenario_red_list_data[[i]] <- replicate_red_list_data
  
}

# Check we have correct structure still
length(scenario_red_list_data) == length(scenario_ab_gl_formatted)
length(scenario_red_list_data[[1]]) == length(scenario_ab_gl_formatted[[1]])

# Have a quick look at the outputs

rli_inputs <- scenario_red_list_data[[1]][[1]]
tail(rli_inputs)

write.csv(rli_inputs, file.path(indicator_outputs_folder, "rli_input_example_annual.csv"))

# Plot some results to check they're not completely whack

## Get one group to check how their status changes over time relative to how
## their abundance changes

# group_id_select <- "13.16.17" # Shows example of 'resurrected' virtual spp
# # group_id_select <- "10.40"
# 
# data <- rli_inputs %>% dplyr::filter(group_id == group_id_select)
# 
# ggplot(data, aes(x = time_step, y = abundance)) +
#   geom_line() +
#   geom_text(aes(label= rl_status,
#                 col = rl_status),hjust=0, vjust=0)




# * Calculate RLI ----

# RLI by individual functional groups

scenario_fg_rli_outputs <- list()
replicate_fg_rli_outputs <- list()

for (i in seq_along(scenario_red_list_data)) {
  
  replicate_red_list_inputs <- scenario_red_list_data[[i]]
  
  for (j in seq_along(replicate_red_list_inputs)) {
    
  replicate_rli <- calculate_red_list_index(
    replicate_red_list_inputs[[j]], numboots, ci = TRUE) %>%
    mutate(replicate = j)
  
  replicate_fg_rli_outputs[[j]] <- replicate_rli 
  
  # saveRDS(replicate_fg_rli_outputs[[j]],
  #         file.path(rli_outputs_folder,
  #                   paste(today, scenarios, "replicate", j,
  #                         "RLI_func_group_output_data.rds",
  #                         sep = "_")))

  write.csv(replicate_fg_rli_outputs[[j]],
            file.path(rli_outputs_folder,
                      paste(today, scenarios, "RLI_func_group_output_data.rds",
                            sep = "_")))
  
  print(paste("RLI for replicate", j, "complete", sep = " "))
  
  }

  scenario_fg_rli_outputs[[i]] <- replicate_fg_rli_outputs

}
  
  
x <- scenario_fg_rli_outputs[[1]][[3]]
head(x)

# Mean RLI aggregated across groups

scenario_rli_outputs <- list()
replicate_rli_outputs <- list()

for (i in seq_along(scenario_fg_rli_outputs)) {
  
  # Get replicates for a single scenario
  replicate_rli_fg <- scenario_fg_rli_outputs[[i]]
  
  # Aggregate RLI across functional groups for each replicate
  for (j in seq_along(replicate_rli_fg)) {
  
   if ("ci_lower" %in% names(replicate_rli_fg[[j]])) {
    
   replicate_rli_outputs[[j]] <- replicate_rli_fg[[j]] %>%
                                 group_by(annual_time_step) %>%
                                 summarise(indicator_score = mean(indicator_score),
                                           ci_lower = mean(ci_lower),
                                           ci_upper = mean(ci_upper)) %>%
                                 mutate(indicator = "RLI",
                                        replicate = j)
   } else {
     
   replicate_rli_outputs[[j]] <- replicate_rli_fg[[j]] %>%
                                 group_by(annual_time_step) %>%
                                 summarise(indicator_score = mean(indicator_score)) %>%
                                 mutate(indicator = "RLI",
                                        replicate = j)
   }

  # saveRDS(replicate_rli_outputs[[j]],
  #       file.path(rli_outputs_folder,
  #                 paste(today, scenarios[[i]], "replicate", j,
  #                       "RLI_aggregate_output_data.rds",
  #                       sep = "_")))
  # 
  # write.csv(replicate_rli_outputs[[j]],
  #           file.path(rli_outputs_folder,
  #                     paste(today, scenarios[[i]], "replicate", j,
  #                           "RLI_aggregate_output_data.rds",
  #                           sep = "_")))

  }
  
  scenario_rli_outputs[[i]] <- replicate_rli_outputs

}

head(scenario_rli_outputs)[[1]][[1]]


# * Plot RLI ----

## By functional group

scenario_fg_rli_plots <- list()
replicate_fg_rli_plots <- list()

for (i in seq_along(scenario_fg_rli_outputs)) {
  
  replicate_fg_rli <- scenario_fg_rli_outputs[[i]]
  
  for (j in seq_along(replicate_fg_rli)) {

  replicate_fg_rli_plots[[j]] <-  plot_red_list_index_by_group(
                                      replicate_fg_rli[[j]],
                                      impact_start,
                                      impact_end,
                                      ci = TRUE)

  ggsave(file.path(rli_plots_folder, paste(today, scenarios[[i]], "replicate", j,
                "RLI_by_functional_group.png",
                sep = "_")),
       replicate_fg_rli_plots[[i]],  device = "png")

  }
  
scenario_fg_rli_plots[[i]] <- replicate_fg_rli_plots

}

scenario_fg_rli_plots[[1]][[3]]

# RLI with all functional groups aggregated
# i.e. mean of each 'taxa' RLI as per Butchart et al (2010) 'Indicators of
# recent declines'

scenario_rli_plots <- list()
replicate_rli_plots <- list()

for (i in seq_along(scenario_rli_outputs)) {
  
  replicate_rli <- scenario_rli_outputs[[i]]
  
  for (j in seq_along(replicate_rli)) {

    replicate_rli_plots[[j]] <- plot_red_list_index(replicate_rli[[j]],
                                                   impact_start, 
                                                   impact_end,
                                                   ci = TRUE)


    ggsave(file.path(rli_plots_folder, paste(today, scenarios[[i]], 
                                             "replicate", j, 
                                             "RLI_aggregated.png",
                                             sep = "_")),
           replicate_rli_plots[[i]],  device = "png")                                   

  }
  
  scenario_rli_plots[[i]] <- replicate_rli_plots

}

scenario_rli_plots[[1]][[1]]

# * Plot all replicates together ----

## Collapse input data so RLI scores for all replicates in one scenario exist in 
## a single data frame

scenario_rli_outputs_aggregated <- list()

for (i in seq_along(scenario_rli_outputs)) {
  
  scenario_rli_outputs_aggregated[[i]] <- do.call(rbind, 
                                                  scenario_rli_outputs[[i]]) %>%
                                          mutate(scenario = scenarios[[i]]) 
  
  
  scenario_mean_rli <- scenario_rli_outputs_aggregated[[i]] %>%
                       group_by(annual_time_step) %>%
                       summarise(indicator_score = mean(indicator_score),
                                 ci_lower = mean(ci_lower),
                                 ci_upper = mean(ci_upper)) %>%
                       mutate(indicator = "RLI",
                              replicate = 0,
                              scenario = scenarios[[i]]) # Replicate 0 will always be the mean
  
  scenario_rli_outputs_aggregated[[i]] <-  rbind(scenario_rli_outputs_aggregated[[i]],
                                                 scenario_mean_rli) %>%
                                           mutate(replicate = as.factor(replicate)) %>%
                                           mutate(level = ifelse(replicate == 0,
                                                                 "Mean RLI", 
                                                                 "Replicate RLI"),
                                                  scenario = scenarios[[i]])
  
}

head(scenario_rli_outputs_aggregated[[1]])
tail(scenario_rli_outputs_aggregated[[1]])

# Plot all together

scenario_rli_plots_aggregated <- list()

for (i in seq_along(scenario_rli_outputs_aggregated)){

scenario_rli_plots_aggregated[[i]] <- ggplot(data = scenario_rli_outputs_aggregated[[i]], 
       aes(x = annual_time_step, y = indicator_score, group = replicate,
           color = level)) +
  geom_line() +
  scale_color_manual(values = c("black", "gray62")) + 
  labs(x = "Time", 
       y = "Red List Index Score") +
  theme(panel.grid.major = element_blank(),
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 18),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "grey97"),
        axis.line = element_line(colour = "black")) +
  geom_vline(xintercept = impact_start, colour = "red") +
  geom_vline(xintercept = impact_end, colour = "blue")

}

scenario_rli_plots_aggregated[[1]]

# LIVING PLANET INDEX ----

# * Create folders ----

lpi_inputs_folder <- file.path(indicator_inputs_folder, "LPI_inputs")

if( !dir.exists( file.path(lpi_inputs_folder) ) ) {
  dir.create( file.path(lpi_inputs_folder), recursive = TRUE )
  
}

lpi_outputs_folder <- file.path(indicator_outputs_folder, "LPI_outputs")

if( !dir.exists( file.path(lpi_outputs_folder) ) ) {
  dir.create( file.path(lpi_outputs_folder), recursive = TRUE )
  
}

lpi_plots_folder <- file.path(indicator_plots_folder, "LPI_plots")

if( !dir.exists( file.path(lpi_plots_folder) ) ) {
  dir.create( file.path(lpi_plots_folder), recursive = TRUE )
  
}

# TEMP CODE ---
## Look at the data we are dealing with

# data <- scenario_abundance_long[[1]][[1]]
# 
# head(data)
# 
# ggplot(data, aes(x = time_step, y = abundance,
#                  col = group_id)) +
#           geom_line()  + 
#           geom_text(aes(label= group_id),hjust=0, vjust=0) +
#           theme(legend.position = "none")

# * Sample data ----

scenario_lpi_inputs <- list()
replicate_lpi_inputs <- list()

for (i in seq_along(scenario_abundance_long)) {
  
  # Get replicates for a single scenario
  replicate_abundance_long <- scenario_abundance_long[[i]]
  
  # Calculate the LPI for each replicate within the scenario
  for (j in seq_along(replicate_abundance_long)) {

  replicate_lpi_inputs[[j]] <- replicate_abundance_long[[j]] %>%
    group_by(group_id) %>% 
    # select rows that are multiples of the specified interval 
    # (eg if interval is 12, it samples one month from every year)
    slice(which(row_number() %% interval == 0)) %>% 
    mutate(annual_time_step = seq(1,max_timestep,1)) #see if it works just hard coding in number of time steps
  
  }
  
  scenario_lpi_inputs[[i]] <- replicate_lpi_inputs

}

# lpi_input <- scenario_lpi_inputs[[1]][[2]]
# head(lpi_input)
# write.csv(lpi_input, file.path(indicator_outputs_folder, "lpi_input_example_annual.csv"))

# * Calculate LPI ----

# Retain naming convention, the LPI just takes the abundance dataframes we
# already formatted while making the RLI inputs

# scenario_lpi_inputs <- scenario_abundance_long

# Loop through each scenario and replicate and calculate the LPI per rep

scenario_lpi_outputs <- list()
replicate_lpi_outputs <- list()

for (i in seq_along(scenario_lpi_inputs)) {
  
  # Get replicates for a single scenario
  replicate_lpi_inputs <- scenario_lpi_inputs[[i]]
  
 # Calculate the LPI for each replicate within the scenario
  for (j in seq_along(replicate_lpi_inputs)) {
    
  replicate_lpi_outputs[[j]] <- calculate_living_planet_index(
    
    replicate_lpi_inputs[[j]], start_time_step, ci = TRUE, numboots, j
  ) 
    
  # Save the output LPI data as a csv and rds
  
    saveRDS(replicate_lpi_outputs[[j]],
          file.path(lpi_outputs_folder,
                    paste(today, scenarios[[i]], "replicate", j,
                          "LPI_output_data.rds",
                          sep = "_")))

    write.csv(replicate_lpi_outputs[[j]],
              file.path(lpi_outputs_folder,
                        paste(today, scenarios[[i]], "replicate", j,
                              "LPI_output_data.rds",
                              sep = "_")))
    
  }
  
  scenario_lpi_outputs[[i]] <- replicate_lpi_outputs
  
}

head(scenario_lpi_outputs)[[1]][[1]]

# * Aggregate all LPI scores ----

## Collapse input data so LPI scores for all replicates in one scenario exist in 
## a single data frame

scenario_lpi_outputs_aggregated <- list()

for (i in seq_along(scenario_lpi_outputs)) {
  
  scenario_lpi_outputs_aggregated[[i]] <- do.call(rbind, 
                                                  scenario_lpi_outputs[[i]]) %>%
                                          mutate(scenario = scenarios[[i]]) 
  
  scenario_mean_lpi <- scenario_lpi_outputs_aggregated[[i]] %>%
    group_by(annual_time_step) %>%
    summarise(indicator_score = mean(indicator_score),
              ci_lower = mean(ci_lower),
              ci_upper = mean(ci_upper)) %>%
    mutate(replicate = 0,# Replicate 0 will always be the mean
           indicator = "LPI",
           scenario = scenarios[[i]]) 
  
  scenario_lpi_outputs_aggregated[[i]] <-  rbind(scenario_lpi_outputs_aggregated[[i]],
                                                 scenario_mean_lpi) %>%
    mutate(replicate = as.factor(replicate)) %>%
    mutate(level = ifelse(replicate == 0,
                          "Mean LPI", 
                          "Replicate LPI"))
  
}

head(scenario_lpi_outputs_aggregated[[1]])
tail(scenario_lpi_outputs_aggregated[[1]])

# * Plot LPI replicates individually ----

scenario_lpi_plots <- list()
replicate_lpi_plots <- list()

for (i in seq_along(scenario_lpi_outputs)) {
  
  replicate_lpi <- scenario_lpi_outputs[[i]]
  
  for (j in seq_along(replicate_lpi)) {
    
    replicate_lpi_plots[[j]] <- plot_living_planet_index(replicate_lpi[[j]],
                                                         ci = TRUE)
    
    
    ggsave(file.path(lpi_plots_folder, paste(today, scenarios[[i]], 
                                             "replicate", j, 
                                             "LPI_aggregated.png",
                                             sep = "_")),
           replicate_lpi_plots[[i]],  device = "png")                                   
    
  }
  
  scenario_lpi_plots[[i]] <- replicate_lpi_plots
  
}

i <- 1
scenario_lpi_plots[[1]][[i]]

i <- i + 1
scenario_lpi_plots[[1]][[i]]

# * Plot all replicates together ----

scenario_lpi_plots_aggregated <- list()

for (i in seq_along(scenario_lpi_outputs_aggregated)){
  
  scenario_lpi_plots_aggregated[[i]] <- ggplot(data = scenario_lpi_outputs_aggregated[[i]], 
                                               aes(x = annual_time_step, 
                                                   y = indicator_score, 
                                                   group = replicate,
                                                   color = level)) +
    geom_line() +
    scale_color_manual(values = c("black", "gray62")) + 
    labs(x = "Time", 
         y = "Living Planet Index Score") +
    theme(panel.grid.major = element_blank(),
          axis.title = element_text(size = 18),
          axis.text = element_text(size = 18),
          panel.grid.minor = element_blank(),
          panel.background = element_rect(fill = "grey97"),
          axis.line = element_line(colour = "black")) +
    geom_vline(xintercept = impact_start, colour = "red") +
    geom_vline(xintercept = impact_end, colour = "blue")
  
}

scenario_lpi_plots_aggregated[[1]]

# Combine indicators ----

all_indicators_list <- list(scenario_rli_outputs,
                            scenario_lpi_outputs)

names(all_indicators_list) <- c("RLI", "LPI")

saveRDS(all_indicators_list,
        file.path(indicator_outputs_folder,
                  paste(today, "all_indicators_output_data_list.rds",
                        sep = "_")))

all_lpi <- do.call(rbind, scenario_lpi_outputs_aggregated) %>% 
           filter(replicate != 0) # Remove the mean so we just have replicates 

all_rli <- do.call(rbind, scenario_rli_outputs_aggregated) %>% 
           filter(replicate != 0) # Remove the mean so we just have replicates

all_indicators <- rbind(all_lpi, all_rli)

saveRDS(all_indicators,
        file.path(indicator_outputs_folder,
                  paste(today, "all_indicators_output_data.rds",
                        sep = "_")))

write.csv(all_indicators,
          file.path(indicator_outputs_folder,
                    paste(today, "all_indicators_output_data.csv",
                          sep = "_")))

