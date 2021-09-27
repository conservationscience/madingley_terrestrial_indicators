
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

## Data visualisation

library(ggridges)
library(cowplot)
library(ggThemeAssist)
library(visdat)

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
      ci_scores <- rep_scores %>% 
        arrange(rep_scores) %>% 
      summarise(ci_lower = quantile(rep_scores$RLI, 
                                                   probs = 0.025),
                               ci_upper = quantile(rep_scores$RLI, 
                                                   probs = 0.975)) %>%
                   mutate(annual_time_step = time) %>% 
        arrange(annual_time_step)
      
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

# For all taxa at once

# data <- scenario_redlist_data_sampled[[2]]

calculate_red_list_index2 <- function(data, numboots, ci = FALSE, replicate_num = NA){
  
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
  
  grouped_data <- weighted_data %>% group_by(annual_time_step)
  
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
      
      ci_scores <- rep_scores %>% 
        arrange(rep_scores) %>% 
        summarise(ci_lower = quantile(rep_scores$RLI, 
                                      probs = 0.025),
                  ci_upper = quantile(rep_scores$RLI, 
                                      probs = 0.975)) %>%
        mutate(annual_time_step = time) %>% 
        arrange(annual_time_step)
      
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
      dplyr::select(annual_time_step,
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
  
  plot <- ggplot() +
    geom_line(data = data, aes(x = annual_time_step, y = indicator_score)) +
    # geom_point(data = raw, aes(x = annual_time_step, y = standardised_abundance), 
    #            alpha = 0.2) +
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
    
    plot <- ggplot() +
      geom_line(data = data, aes(x = annual_time_step, y = indicator_score)) +
      # geom_point(data = raw, aes(x = annual_time_step, y = standardised_abundance), 
      #            alpha = 0.2) +
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
    rename(ave_abundance = abundance) %>% 
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
  

  #burnin_months <- 1000*12 # in months
  burnin_months <- 900*12 # in months - burnin complete after 1000 but need the extra years to get RL status
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
                           "/Serengeti/Outputs_from_indicator_code/Indicator_inputs",
                           today)

if( !dir.exists( file.path(indicator_inputs_folder) ) ) {
  dir.create( file.path(indicator_inputs_folder), recursive = TRUE )
  
}

indicator_outputs_folder <- file.path(indicators_project, 
                                     "/Serengeti/Outputs_from_indicator_code/Indicator_outputs",
                                     today)

if( !dir.exists( file.path(indicator_outputs_folder) ) ) {
  dir.create( file.path(indicator_outputs_folder), recursive = TRUE )
  
}

general_indicator_outputs_folder <- file.path(indicators_project, 
                                      "/Serengeti/Outputs_from_indicator_code/Indicator_outputs/general",
                                      today)

if( !dir.exists( file.path(general_indicator_outputs_folder) ) ) {
  dir.create( file.path(general_indicator_outputs_folder), recursive = TRUE )
  
}

indicator_plots_folder <- file.path(indicators_project, 
                                     "/Serengeti/Outputs_from_indicator_code/Indicator_plots",
                                    today)

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

# Get groups & species ----
# Details for the virtual species or 'groups'
# Groups file is the same for all simulations, so can just pull it from whichever directory

groups <- readRDS(file.path(processed_simulation_paths[[1]][1], "groups.rds"))

key <- readRDS(file.path(processed_simulation_paths[[1]][1], "species_and_groups_key.rds"))

species <- readRDS("N:/Quantitative-Ecology/Indicators-Project/Serengeti/Outputs_from_adaptor_code/map_of_life/comparable_taxa.rds")

# Add group_id to species

species <- species %>% 
           merge(key, by = "species_id")

# Just take the first species matched with each group id as an example
example_species <- species %>% 
                 group_by(group_id) %>% 
                 #slice(1) %>% 
                 dplyr::select(group_id, common_name, accepted_name, nutrition_source,
                                 thermoregulation, bodymass)

all_species <- split(species, species$group_id)

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

# Get modelled generation lengths ----

# Get the maximum generation length for each group ID across all replicates and 
# scenarios

scenario_gen_length_dfs <- list()

for (i in seq_along(scenario_generations_raw)) {
  
  scenario_gen_length_dfs[[i]] <- do.call(rbind, scenario_generations_raw[[i]]) %>% 
        group_by(group_id) %>% 
        summarise(mean_gen_length = mean(generation_length_yrs, na.rm = TRUE),
                  median_gen_length = median(generation_length_yrs, na.rm = TRUE),
                  max_gen_length = max(generation_length_yrs, na.rm = TRUE))
    
}

gen_lengths <- do.call(rbind, scenario_gen_length_dfs) %>% 
               group_by(group_id) %>% 
               summarise(mean_gen_length = mean(mean_gen_length, na.rm = TRUE),
               median_gen_length = median(median_gen_length, na.rm = TRUE),
               max_gen_length = max(max_gen_length, na.rm = TRUE))


model_species <- gen_lengths %>% 
  merge(example_species, by = "group_id", all = TRUE) %>% 
  merge(groups[c("group_id", "nutrition_source","endo_ectotherm",
                 "mass_lower", "mass_upper")], 
        by = "group_id") %>% 
  dplyr::select(-nutrition_source.x, -thermoregulation) %>% 
  rename(thermoregulation = endo_ectotherm,
         nutrition_source = nutrition_source.y) %>% 
  mutate(matched = ifelse(is.na(mean_gen_length),
                          "real spp not represented in model",
                          ifelse(is.na(bodymass),
                                 "mod spp not represented in real spp list",
                                 "model spp matches real spp"))) %>% 
  mutate(bodymass_grams = bodymass) %>% 
  dplyr::select(-bodymass) %>% 
  dplyr::select(group_id, mean_gen_length, median_gen_length, max_gen_length,
                common_name, accepted_name, nutrition_source, thermoregulation, 
                matched, bodymass_grams, mass_lower, mass_upper)

write.csv(model_species, file.path(general_indicator_outputs_folder,
                                   "model_species.csv"))

# Get calculated generation lengths ----

generation_lengths <- read.csv(file.path(indicators_project, location,
                                         "Inputs_to_adaptor_code",
                                         "Generation_lengths",
                                         "Generation_lengths.csv"),
                               colClasses = c(GL_FIN = "numeric"))
head(generation_lengths)

generation_lengths <- generation_lengths %>% 
                      rename(group_id = FG,
                             generation_length_yrs = GL_FIN) %>% 
                      select(group_id, generation_length_yrs) 

class(generation_lengths$generation_length_yrs)

# SUBSET FOR TESTING ----
## (optional, leave commented out to continue applying to full dataset)

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
rm(scenario_abundance_raw)

# Do the same for the autotrophs

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
x <- scenario_auto_long[[1]][[1]]
rm(scenario_autotroph_raw)

# rm(abundance_reps, abundance_temp, replicate_abundance_formatted, rep)

# Check structure is still correct

length(scenario_abundance_formatted) == length(scenario_abundance_raw)


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

rm(scenario_abundance_formatted)

# ## AVERAGE REPLICATES - DECOMPOSED ----

# Merge abundance and generation length ----

scenario_ab_gl_formatted_not_clean <- list()

for (i in seq_along(scenario_abundance_long)) {
  
  replicate_abundance <- scenario_abundance_long[[i]]
  #replicate_generations <- scenario_generations_raw[[i]]
  
  # Make a list to catch the outputs
  
  replicate_ab_gl_formatted <- list()
  
  # For each individual replicate
  
  for (j in seq_along(replicate_abundance)) {
    
    # Reduce size of the replicate generations dataframe or the merge won't work
    # gen_length <- replicate_generations[[j]] %>% 
    #   dplyr::select(group_id, generation_length_yrs, 
    #                 functional_group_name) %>% 
    #   distinct(.)
    
    # Add the generation length info to the abundance dataframe
    replicate_ab_gl_formatted[[j]] <- replicate_abundance[[j]] %>%
      merge(generation_lengths, by = "group_id") %>%
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
      merge(groups[c("group_id", "functional_group_name","bodymass_index", "mass_lower")], 
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
  
  rm(replicate_ab_gl_formatted)
  
}

test <- scenario_ab_gl_formatted_not_clean[[1]][[1]]
head(test)

rm(scenario_abundance_long)

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

rm(scenario_ab_gl_formatted_not_clean)

# Remove replicates with no carnivorous endotherms ----

scenario_ab_gl_formatted <- list()

for (i in seq_along(scenario_false_extinctions_removed)) {
  
  replicate_not_clean <- scenario_false_extinctions_removed[[i]]
  
  scenario_ab_gl_formatted[[i]] <- list.clean(replicate_not_clean)
  
}

test <- scenario_ab_gl_formatted[[1]][[1]]
head(test)
any(is.nan(test$abundance))
rm(scenario_false_extinctions_removed)


saveRDS(scenario_ab_gl_formatted,
        file.path(general_indicator_outputs_folder,
                  "formatted_abundance_1.rds"))

# scenario_ab_gl_formatted <- readRDS(file.path(indicator_outputs_folder,
#                   paste(today, "formatted_abundance_1.rds", sep = "_")))


# Truncate and decompose trends ----

scenario_decomposed <- list()
scenario_auto_decomposed <- list()

for (i in seq_along(scenario_ab_gl_formatted)) {

  # Heterotrophs

  # Get reps for one scenario

  replicates <- scenario_ab_gl_formatted[[i]]

  replicates_decomposed <- list()
 
  for (j in seq_along(replicates)) {

    replicate <- replicates[[j]]
    
    # Remove groups with no abundance
    
    replicate <- replicate %>% 
            group_by(group_id) %>% 
            filter(!all(is.na(abundance)))
    
    # Split into individual groups
    
    rep_groups <- split(replicate, replicate$group_id)

    groups_decomposed <- list()

    for (k in seq_along(rep_groups)) {
      
      # Check the time series is long enough
      
      ts_length <- nrow(rep_groups[[k]] %>% filter(abundance > 0))
      
      # If not, add a null element to the list
      
      if (ts_length < 24) {
      
      groups_decomposed[[k]] <- NULL
        
      } else { #otherwise proceed
      
      # Get bottom and top 10% percentile limits
      abundance_lims <- quantile(rep_groups[[k]]$abundance,
                              probs=c(.10, .90), na.rm = TRUE)

      # ggplot(data = rep_groups[[i]], aes(x = monthly_time_step,
      #                                   y = abundance)) +
      #   geom_point()
      
      print(rep_groups[[k]]$group_id[1])
      
      # Decompose the seasonal and random noise from the trend using a moving
      # average https://www.rdocumentation.org/packages/stats/versions/3.6.2/topics/decompose
      # Chops off the first and last 5 observations
      
      with_vals <- rep_groups[[k]] %>%
        mutate(abundance_0 = abundance) %>% #Save un-altered abundance
        mutate(abundance = ifelse(abundance_0 > abundance_lims[2], NA, #Create truncated abundance
                                  ifelse(abundance_0 < abundance_lims[1],NA,
                                         abundance_0))) %>%
        # If extinction has occurred replace any +ve or NA with 0
        mutate(abundance = ifelse(true_extinction == "true extinction",
                                  0, abundance)) %>% 
        filter(!is.na(abundance)) %>% 
        mutate(abundance_1 = abundance, #save truncated abundance
               timeseries = ts(abundance_1, frequency = 12),# convert it to timeseries obj
               abundance = decompose(timeseries)$trend) # make decomposed trend the abundance variable
      
      without_vals <- rep_groups[[k]] %>%
        mutate(abundance_0 = abundance) %>% #Save un-altered abundance
        mutate(abundance = ifelse(abundance_0 > abundance_lims[2], NA,
                                  ifelse(abundance_0 < abundance_lims[1],NA,
                                         abundance_0))) %>%
        filter(is.na(abundance)) %>% 
        mutate(abundance_1 = abundance,
               timeseries = NA,
               abundance = abundance_1)
      
      
      groups_decomposed[[k]] <- rbind(with_vals, without_vals) %>% 
                                arrange(monthly_time_step) %>%                         
                                slice(-1) %>% # Make sure we have exactly 3600 timesteps
                                slice(1:3600) %>% 
                                mutate(annual_time_step = rep(1:300, each = 12),
                                       replicate = j - 1) %>% # Add annual timestep
                                arrange(monthly_time_step) 
      
      # ggplot(data = x, aes(x = monthly_time_step,
      #                                    y = abundance)) +
      #   geom_point()

      }
    }
   
   #clear out null elements
      
   groups_decomposed <- list.clean(groups_decomposed)
   
   # Bind and add to list
   
   df <- do.call(rbind, groups_decomposed)
   
   replicates_decomposed[[j]] <- df
   
  }
  
  scenario_decomposed[[i]] <- replicates_decomposed

  print(paste(scenarios[[i]], "replicates decomposed"))
  
  rm(replicate, replicates_decomposed, groups_decomposed)

}

saveRDS(scenario_decomposed,
        file.path(general_indicator_outputs_folder,
                  "decomposed_abundance_2.rds"))

# head(scenario_decomposed[[1]][[1]])

# Repeat for autotrophs

scenario_auto_decomposed <- list()

for (i in seq_along(scenario_auto_long)) {
  
  # Get reps for one scenario
  
  replicates <- scenario_auto_long[[i]]
  
  replicates_decomposed <- list()
  
  for (j in seq_along(replicates)) {
    
    replicate <- replicates[[j]]
    
    # Remove groups with no abundance
    
    replicate <- replicate %>% 
      group_by(group_id) %>% 
      filter(!all(is.na(abundance)))
    
    # Split into individual groups
    
    rep_groups <- split(replicate, replicate$group_id)
    
    groups_decomposed <- list()
    
    for (k in seq_along(rep_groups)) {
      
      # Check the time series is long enough
      
      ts_length <- nrow(rep_groups[[k]] %>% filter(abundance > 0))
      
      # If not, add a null element to the list
      
      if (ts_length < 24) {
        
        groups_decomposed[[k]] <- NULL
        
      } else { #otherwise proceed
        
        # Get bottom and top 10% percentile limits
        abundance_lims <- quantile(rep_groups[[k]]$abundance,
                                   probs=c(.10, .90), na.rm = TRUE)
        
        # ggplot(data = rep_groups[[i]], aes(x = monthly_time_step,
        #                                   y = abundance)) +
        #   geom_point()
        
        print(rep_groups[[k]]$group_id[1])
        
        # Decompose the seasonal and random noise from the trend using a moving
        # average https://www.rdocumentation.org/packages/stats/versions/3.6.2/topics/decompose
        # Chops off the first and last 5 observations
        
        with_vals <- rep_groups[[k]] %>%
          mutate(abundance_0 = abundance) %>% #Save un-altered abundance
          mutate(abundance = ifelse(abundance_0 > abundance_lims[2], NA, #Create truncated abundance
                                    ifelse(abundance_0 < abundance_lims[1],NA,
                                           abundance_0))) %>%
          filter(!is.na(abundance)) %>% 
          mutate(abundance_1 = abundance, #save truncated abundance
                 timeseries = ts(abundance_1, frequency = 12),# convert it to timeseries obj
                 abundance = decompose(timeseries)$trend) # make decomposed trend the abundance variable
        
        without_vals <- rep_groups[[k]] %>%
          mutate(abundance_0 = abundance) %>% #Save un-altered abundance
          mutate(abundance = ifelse(abundance_0 > abundance_lims[2], NA,
                                    ifelse(abundance_0 < abundance_lims[1],NA,
                                           abundance_0))) %>%
          filter(is.na(abundance)) %>% 
          mutate(abundance_1 = abundance,
                 timeseries = NA,
                 abundance = abundance_1)
        
        
        groups_decomposed[[k]] <- rbind(with_vals, without_vals) %>% 
          arrange(monthly_time_step) %>%                         
          #slice(-1) %>% # Make sure we have exactly 3600 timesteps
          slice(1:3600) %>% 
          mutate(annual_time_step = rep(1:300, each = 12),
                 replicate = j - 1) %>% # Add annual timestep
          arrange(monthly_time_step) 
        
        # ggplot(data = x, aes(x = monthly_time_step,
        #                                    y = abundance)) +
        #   geom_point()
        
      }
    }
    
    #clear out null elements
    
    groups_decomposed <- list.clean(groups_decomposed)
    
    # Bind and add to list
    
    df <- do.call(rbind, groups_decomposed)
    
    replicates_decomposed[[j]] <- df
    
  }
  
  scenario_auto_decomposed[[i]] <- replicates_decomposed
  
  print(paste(scenarios[[i]], "replicates decomposed"))
  
  rm(replicate, replicates_decomposed, groups_decomposed)
  
}

saveRDS(scenario_auto_decomposed,
        file.path(general_indicator_outputs_folder,
                  "decomposed_auto_abundance_2.rds"))


# Take a random sample of 25 replicates ----

scenario_decomposed_25 <- list()
scenario_auto_decomposed_25 <- list()

for (i in seq_along(scenario_decomposed)) {

  if(scenarios[[i]] == "000_Baseline") {

  scenario_decomposed_25[[i]] <- scenario_decomposed[[i]]
  scenario_auto_decomposed_25[[i]] <- scenario_auto_decomposed[[i]]

  } else {

  replicates <- scenario_decomposed[[i]]
  replicates_auto <- scenario_auto_decomposed[[i]]

  set.seed(159)

  scenario_decomposed_25[[i]] <- sample(scenario_decomposed[[i]], 25)
  scenario_auto_decomposed_25[[i]] <- sample(scenario_auto_decomposed[[i]], 25) # Will this sample the same reps??

  }

}

# Average across replicates ----

scenario_averaged <- list()
scenario_auto_averaged <- list()

# for (i in seq_along(scenario_ab_gl_formatted_25)) {

  for (i in seq_along(scenario_decomposed_25)) {
  
  # Heterotrophs
  
  # replicates <- scenario_ab_gl_formatted_25[[i]]
  
  replicates <- scenario_decomposed_25[[i]]
  
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
                                          -true_extinction, - replicate, - last_abundance,
                                          -abundance_0,-abundance_1, - timeseries) %>% 
                            rename(abundance = mean_abundance,
                                   generation_length_yrs = mean_gen_length,
                                   timeframe = mean_timeframe) %>% 
                            distinct(.) %>% 
                            ungroup(.)
  
  rm(replicate_df, replicates, n)
  
  #Autotrophs
  
  auto_replicates <- scenario_auto_decomposed_25[[i]]
  
  n2 <- length(auto_replicates)
  
  auto_replicate_df <- do.call(rbind, auto_replicates)
  
  scenario_auto_averaged[[i]] <- auto_replicate_df %>% 
    group_by(group_id, monthly_time_step) %>% 
    mutate(mean_abundance = mean(abundance, na.rm = TRUE),
           # 95% CIs using this info: https://www.cyclismo.org/tutorial/R/confidence.html
           sd = sd(abundance, na.rm = TRUE),
           error = qt(0.975, df = n2 - 1) * 
             sd/sqrt(n2),
           lower_ci = mean_abundance - error, 
           upper_ci = mean_abundance + error) %>%
    dplyr::select(-abundance, - replicate) %>% 
    rename(abundance = mean_abundance) %>% 
    distinct(.) %>% 
    ungroup(.)
  
  rm(auto_replicate_df, auto_replicates, n2)
  
  print(paste(scenarios[[i]], "replicates averaged"))
  
}

baseline <- scenario_averaged[[1]]
head(baseline)
any(is.nan(baseline$abundance))

saveRDS(scenario_averaged,
        file.path(general_indicator_outputs_folder,
                  paste("averaged_abundance_3.rds", sep = "_")))

baseline <- scenario_auto_averaged[[1]]
head(baseline)
any(is.nan(baseline$abundance))

saveRDS(scenario_auto_averaged,
        file.path(general_indicator_outputs_folder,
                  "averaged_auto_abundance_3.rds"))

# Check for false extinctions

## Get groups that have experienced 0 abundance

# zeroes <- baseline %>% 
#           group_by(group_id) %>% 
#           filter(any(abundance == 0))
# 
# zeroes_list <- split(zeroes, zeroes$group_id)
# 
# b <- b+1
# x <- zeroes_list[[b]]

group_test <- "13.16.9"
group <- baseline %>%  filter(group_id == group_test)
head(group)

# 13 has vals
h <- h + 1
og_abundance <- scenario_decomposed[[1]][[13]] %>% filter(group_id == group_test)
dim(og_abundance)
print(h)

ggplot() +
  geom_line(data = og_abundance, 
             aes(x = monthly_time_step, y = abundance_0), col = "yellow") +
  geom_line(data = og_abundance,
             aes(x = monthly_time_step, y = abundance_1), alpha = 0.9, col = "light blue") +
  geom_line(data = group, aes(x = monthly_time_step, y = abundance), col = "red") +
  labs(title = "averaged") +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank())

test3 <- scenario_auto_averaged[[1]]
head(test3)

# saveRDS(scenario_averaged, file.path(indicator_outputs_folder, 
#                                      paste(today, "scenarios_averaged.RDS")))

# Remove blinking groups ----

scenario_abundance_clean <- list()

for (i in seq_along(scenario_averaged)) {

    # data <- data %>% filter(group_id == "15.18.21") # For testing

    # Determine which groups were there at beginning (post burnin)
    temp<- scenario_averaged[[i]] %>%
      group_by(group_id) %>%
      filter(abundance != is.nan(abundance)) %>%
      filter(monthly_time_step == min(monthly_time_step)) %>%
      mutate(first_appearance = monthly_time_step,
             beginning = ifelse(first_appearance <= 12006,
                                TRUE, FALSE)) %>% 
      dplyr::select(group_id, first_appearance, beginning)
    
    scenario_abundance_clean[[i]] <- temp %>% 
                                     merge(scenario_averaged[[i]]) %>% 
                                     tidylog::filter(beginning == TRUE)
    
    rm(temp)
}

sc_av <- scenario_averaged[[1]]
sc_ab <- scenario_abundance_clean[[1]]

# Check we have the same number of species at beginning and end

x1 <- scenario_abundance_clean[[1]] %>% filter(monthly_time_step == 12006)
x300 <- scenario_abundance_clean[[1]] %>% filter(monthly_time_step == max(monthly_time_step))

unique(x1$group_id) == unique(x300$group_id)
x1x300 <- cbind(x1$group_id,x300$group_id)

# Remove false extinctions from averaged vals ----

scenario_remove_false_extinctions_2 <- list()

# for (i in seq_along(scenario_averaged)) {
#   
#   # Identify the last time-step a species has abundance values > 0
#   temp2 <- scenario_averaged[[i]] %>% 

for (i in seq_along(scenario_abundance_clean)) {
  
  # Identify the last time-step a species has abundance values > 0
  temp2 <- scenario_abundance_clean[[i]] %>%
           group_by(group_id) %>% 
           slice(1:3600) %>% 
           filter(abundance > 0) %>% 
           dplyr::select(group_id, monthly_time_step, abundance) %>% 
           filter(monthly_time_step == max(monthly_time_step)) %>% 
           dplyr::select(group_id, monthly_time_step) %>% 
           rename(last_abundance = monthly_time_step) %>% 
           ungroup()
  
  dim(temp2)
  
  # Add the year of last positive abundance number as a column to the data    
  temp3 <- scenario_abundance_clean[[i]] %>% 
    #remove existing last abundance val bc it's leftover from before averaging
    # dplyr::select(-last_abundance) %>% 
    merge(temp2, by = "group_id", all = TRUE)
  
  dim(temp3)
  
  rm(temp2)
  
  # Use the last positive abundance year and current abundance value to determine
  # if a zero abundance is a true extinction or just a missing value (false extinction)
  scenario_remove_false_extinctions_2[[i]] <- temp3 %>%
    group_by(group_id) %>%
    slice(1:3600) %>% 
    mutate(true_extinction = ifelse(abundance == 0 &
                                    monthly_time_step < last_abundance,
                                    "false extinction",
                                    ifelse(abundance > 0 &
                                           monthly_time_step < last_abundance,
                                    "not extinct",
                                    ifelse(abundance == 0 &
                                          monthly_time_step >= last_abundance,
                                    "true extinction", "not extinct")))) %>%
    mutate(abundance = ifelse(true_extinction == "false extinction",
                              NA, abundance)) %>%
    group_by(group_id) %>%
    arrange(group_id, monthly_time_step)  %>% 
    ungroup()
  
  rm(temp3)
}

x1 <- scenario_remove_false_extinctions_2[[1]] %>% 
      filter(monthly_time_step == 12006)

x300 <- scenario_remove_false_extinctions_2[[1]] %>% 
        filter(monthly_time_step == max(monthly_time_step))

dim(x1)[1] == dim(x300)[1]

saveRDS(scenario_remove_false_extinctions_2,
        file.path(general_indicator_outputs_folder,
                  "weird_groups_removed_abundance_4.rds"))

# # Smooth abundance ----

ave_window <- 12 * 5

scenario_smoothed_abundance <- list()

for (i in seq_along(scenario_remove_false_extinctions_2)) {

    group_list <- split(scenario_remove_false_extinctions_2[[i]], 
                        scenario_remove_false_extinctions_2[[i]]$group_id)

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
    
    rm(group_smoothed_abundance, group_df)
}

x1 <- scenario_smoothed_abundance[[1]] %>% filter(monthly_time_step == 12006)
x300 <- scenario_smoothed_abundance[[1]] %>% filter(monthly_time_step == max(monthly_time_step))
xx1 <- x1$group_id
xx300 <- x300$group_id

unique(x1$group_id) == unique(x300$group_id)

saveRDS(scenario_smoothed_abundance,
        file.path(general_indicator_outputs_folder,
                  "smoothed_abundance_5.rds"))

## Autotrophs

ave_window <- 12 * 5

scenario_auto_smoothed_abundance <- list()

for (i in seq_along(scenario_auto_averaged)) {
  
  group_list <- split(scenario_auto_averaged[[i]], 
                      scenario_auto_averaged[[i]]$group_id)
  
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
  
  scenario_auto_smoothed_abundance[[i]]  <- do.call(rbind,
                                               group_smoothed_abundance)
  
  rm(group_smoothed_abundance, group_df)
}

saveRDS(scenario_auto_smoothed_abundance,
        file.path(general_indicator_outputs_folder,
                  "smoothed_auto_abundance_5.rds"))

# Plot data transformation ----

## Get an example group for each scenario (random for baseline, mid size herb for land use,
## big carn for carnivores, big herb for herbivores)
example_group <- c("10.59", "10.40", "11.68", "10.68")

# Raw abundance
r <- 1 # Pick a replicate

## * Decomposition ----

### Plot the change from raw data to truncated and decomposed

decomposition_plots <- list()

for ( i in seq_along(scenario_ab_gl_formatted)) {

group_raw <- scenario_ab_gl_formatted[[i]][[r]] %>% 
             filter(group_id == example_group[[i]]) %>% 
             mutate(data_processing_level = "1 - raw for single replicate")

# Decomposed abundance

group_decomposed <- scenario_decomposed[[i]][[r]] %>% 
                    filter(group_id == example_group[[i]]) %>% 
                    mutate(data_processing_level = "2 - 80th percentile & decomposed for single replicate")

plot_data <- rbind(group_raw, group_decomposed)

decomposition_plots[[i]] <- ggplot(data = plot_data, aes(x = monthly_time_step, 
                                                         y = abundance,
                                                         col = data_processing_level)) +
  geom_point() +
  scale_color_manual(values = c("light blue", "pink")) +
  theme(panel.grid.major = element_line(linetype = "blank"), 
        panel.grid.minor = element_line(linetype = "blank"), 
        panel.background = element_rect(fill = "gray49"),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        legend.position = "bottom") +
  labs(title = scenarios[[i]])

rm(group_raw, group_decomposed, plot_data)

}

decomposition_plots[[1]]
decomposition_plots[[2]]
decomposition_plots[[3]]
decomposition_plots[[4]]

## * Smoothing ----

smoothing_plots <- list()

for ( i in seq_along(scenario_averaged)) {
  
  group_decomposed <- scenario_decomposed[[i]][[r]] %>% 
    filter(group_id == example_group[[i]]) %>% 
    mutate(data_processing_level = "2 - 80th percentile & decomposed for one replicate")
  
  # Compare to averaged (not smoothed)
  group_averaged <- scenario_averaged[[i]] %>% 
    filter(group_id == example_group[[i]]) %>% 
    mutate(data_processing_level = "3 - Averaged across replicates")
  
  # Compare to averaged and false extinctions gone
  group_smoothed <- scenario_smoothed_abundance[[i]] %>% 
    filter(group_id == example_group[[i]]) %>% 
    mutate(data_processing_level = "4 - Smoothed over five year rollin window") %>% 
    dplyr::select(-abundance) %>% 
    rename(abundance = ave_abundance)
  
  plot_data <- rbind(group_decomposed, group_averaged, group_smoothed)

  smoothing_plots[[i]] <-  ggplot() +
    geom_point(data = plot_data, aes(x = monthly_time_step, 
                                     y = abundance,
                                     col = data_processing_level),
               alpha = 0.2) +
    scale_color_manual(values = c("light blue", "purple", "yellow")) +
    theme(panel.grid.major = element_line(linetype = "blank"), 
          panel.grid.minor = element_line(linetype = "blank"), 
          panel.background = element_rect(fill = "gray49"),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          legend.position = "bottom") +
    labs(title = scenarios[[i]])
  
  
  rm(group_smoothed, group_averaged, plot_data, group_decomposed)
  
}

smoothing_plots[[1]]
smoothing_plots[[2]]
smoothing_plots[[3]]
smoothing_plots[[4]]


# # REPEAT PROCESSING WITH ONLY 5 REPS ----
# 
# # Take a random sample of 5 replicates ----
# 
# scenario_ab_gl_formatted_5 <- list()
# #scenario_auto_long_5 <- list()
# 
# for (i in seq_along(scenario_ab_gl_formatted)) {
# 
#   replicates <- scenario_ab_gl_formatted[[i]]
#   #replicates_auto <- scenario_auto[[i]]
# 
#   set.seed(159)
# 
#   scenario_ab_gl_formatted_5[[i]] <- sample(scenario_ab_gl_formatted[[i]], 5)
#   #scenario_auto_long_5[[i]] <- sample(scenario_auto[[i]], 5) # Will this sample the same reps??
# 
# }
# 
# 
# # Truncate and decompose trends ----
# 
# scenario_decomposed_5 <- list()
# #scenario_auto_decomposed <- list()
# 
# for (i in seq_along(scenario_ab_gl_formatted_5)) {
#   
#   # Heterotrophs
#   
#   # Get reps for one scenario
#   
#   replicates <- scenario_ab_gl_formatted_5[[i]]
#   
#   replicates_decomposed <- list()
#   
#   for (j in seq_along(replicates)) {
#     
#     replicate <- replicates[[j]]
#     
#     # Remove groups with no abundance
#     
#     replicate <- replicate %>% 
#       group_by(group_id) %>% 
#       filter(!all(is.na(abundance)))
#     
#     # Split into individual groups
#     
#     rep_groups <- split(replicate, replicate$group_id)
#     
#     groups_decomposed <- list()
#     
#     for (k in seq_along(rep_groups)) {
#       
#       # Check the time series is long enough
#       
#       ts_length <- nrow(rep_groups[[k]] %>% filter(abundance > 0))
#       
#       # If not, add a null element to the list
#       
#       if (ts_length < 24) {
#         
#         groups_decomposed[[k]] <- NULL
#         
#       } else { #otherwise proceed
#         
#         # Get bottom and top 10% percentile limits
#         abundance_lims <- quantile(rep_groups[[k]]$abundance,
#                                    probs=c(.10, .90), na.rm = TRUE)
#         
#         # ggplot(data = rep_groups[[i]], aes(x = monthly_time_step,
#         #                                   y = abundance)) +
#         #   geom_point()
#         
#         print(rep_groups[[k]]$group_id[1])
#         
#         # Decompose the seasonal and random noise from the trend using a moving
#         # average https://www.rdocumentation.org/packages/stats/versions/3.6.2/topics/decompose
#         # Chops off the first and last 5 observations
#         
#         with_vals <- rep_groups[[k]] %>%
#           mutate(abundance_0 = abundance) %>% #Save un-altered abundance
#           mutate(abundance = ifelse(abundance_0 > abundance_lims[2], NA, #Create truncated abundance
#                                     ifelse(abundance_0 < abundance_lims[1],NA,
#                                            abundance_0))) %>%
#           # If extinction has occurred replace any +ve or NA with 0
#           mutate(abundance = ifelse(true_extinction == "true extinction",
#                                     0, abundance)) %>% 
#           filter(!is.na(abundance)) %>% 
#           mutate(abundance_1 = abundance, #save truncated abundance
#                  timeseries = ts(abundance_1, frequency = 12),# convert it to timeseries obj
#                  abundance = decompose(timeseries)$trend) # make decomposed trend the abundance variable
#         
#         without_vals <- rep_groups[[k]] %>%
#           mutate(abundance_0 = abundance) %>% #Save un-altered abundance
#           mutate(abundance = ifelse(abundance_0 > abundance_lims[2], NA,
#                                     ifelse(abundance_0 < abundance_lims[1],NA,
#                                            abundance_0))) %>%
#           filter(is.na(abundance)) %>% 
#           mutate(abundance_1 = abundance,
#                  timeseries = NA,
#                  abundance = abundance_1)
#         
#         
#         groups_decomposed[[k]] <- rbind(with_vals, without_vals) %>% 
#           arrange(monthly_time_step) %>%                         
#           slice(-1) %>% # Make sure we have exactly 3600 timesteps
#           slice(1:3600) %>% 
#           mutate(annual_time_step = rep(1:300, each = 12)) %>% # Add annual timestep
#           arrange(monthly_time_step) 
#         
#         # ggplot(data = x, aes(x = monthly_time_step,
#         #                                    y = abundance)) +
#         #   geom_point()
#         
#       }
#     }
#     
#     #clear out null elements
#     
#     groups_decomposed <- list.clean(groups_decomposed)
#     
#     # Bind and add to list
#     
#     df <- do.call(rbind, groups_decomposed)
#     
#     replicates_decomposed[[j]] <- df
#     
#   }
#   
#   scenario_decomposed_5[[i]] <- replicates_decomposed
#   
#   print(paste(scenarios[[i]], "replicates decomposed"))
#   
#   rm(replicate, replicates_decomposed, groups_decomposed)
#   
# }
# 
# vis_dat(scenario_decomposed_5[[1]][[1]][1:50,])
# vis_dat(scenario_decomposed_5[[2]][[1]][1:50,])
# vis_dat(scenario_decomposed_5[[3]][[1]][1:50,])
# vis_dat(scenario_decomposed_5[[4]][[1]][1:50,])
# 
# # Average across replicates ----
# 
# scenario_averaged_5 <- list()
# #scenario_auto_averaged <- list()
# 
# # for (i in seq_along(scenario_ab_gl_formatted_25)) {
# 
# for (i in seq_along(scenario_decomposed_5)) {
#   
#   # Heterotrophs
#   
#   # replicates <- scenario_ab_gl_formatted_25[[i]]
#   
#   replicates <- scenario_decomposed_5[[i]]
#   
#   n <- length(replicates)
#   
#   replicate_df <- do.call(rbind, replicates)
#   
#   scenario_averaged_5[[i]] <- replicate_df %>% 
#     mutate(replicate = j - 1) %>% 
#     group_by(group_id, monthly_time_step) %>% 
#     mutate(mean_abundance = mean(abundance, na.rm = TRUE),
#            # 95% CIs using this info: https://www.cyclismo.org/tutorial/R/confidence.html
#            sd = sd(abundance, na.rm = TRUE),
#            error = qt(0.975, df = n - 1) * 
#              sd/sqrt(n),
#            lower_ci = mean_abundance - error, 
#            upper_ci = mean_abundance + error,  
#            mean_gen_length = mean(generation_length_yrs, na.rm = TRUE),
#            mean_timeframe = round(mean(timeframe, na.rm = TRUE))) %>%
#     dplyr::select(-abundance, -generation_length_yrs, -timeframe,
#                   -true_extinction, - replicate, - last_abundance,
#                   -abundance_0,-abundance_1, - timeseries) %>% 
#     rename(abundance = mean_abundance,
#            generation_length_yrs = mean_gen_length,
#            timeframe = mean_timeframe) %>% 
#     distinct(.) %>% 
#     ungroup(.)
#   
#   rm(replicate_df, replicates, n)
#   
#   # #Autotrophs
#   # 
#   # auto_replicates <- scenario_auto_long[[i]]
#   # 
#   # n2 <- length(auto_replicates)
#   # 
#   # auto_replicate_df <- do.call(rbind, auto_replicates)
#   # 
#   # scenario_auto_averaged[[i]] <- auto_replicate_df %>% 
#   #   group_by(group_id, monthly_time_step) %>% 
#   #   mutate(mean_abundance = mean(abundance, na.rm = TRUE),
#   #          # 95% CIs using this info: https://www.cyclismo.org/tutorial/R/confidence.html
#   #          sd = sd(abundance, na.rm = TRUE),
#   #          error = qt(0.975, df = n2 - 1) * 
#   #            sd/sqrt(n2),
#   #          lower_ci = mean_abundance - error, 
#   #          upper_ci = mean_abundance + error) %>%
#   #   dplyr::select(-abundance, - replicate) %>% 
#   #   rename(abundance = mean_abundance) %>% 
#   #   distinct(.) %>% 
#   #   ungroup(.)
#   # 
#   # rm(auto_replicate_df, auto_replicates, n2)
#   
#   print(paste(scenarios[[i]], "replicates averaged"))
#   
# }
# 
# baseline <- scenario_averaged_5[[1]]
# 
# head(baseline)
# tail(baseline)
# 
# vis_dat(scenario_averaged_5[[1]][1:50,])
# vis_dat(scenario_averaged_5[[2]][1:50,])
# vis_dat(scenario_averaged_5[[3]][1:50,])
# vis_dat(scenario_averaged_5[[4]][1:50,])
# 
# 
# # saveRDS(scenario_averaged, file.path(indicator_outputs_folder, 
# #                                      paste(today, "scenarios_averaged.RDS")))
# 
# # Remove blinking groups ----
# 
# scenario_abundance_clean_5 <- list()
# 
# for (i in seq_along(scenario_averaged_5)) {
#   
#   # data <- data %>% filter(group_id == "15.18.21") # For testing
#   
#   # Determine which groups were there at beginning (post burnin)
#   temp<- scenario_averaged_5[[i]] %>%
#     group_by(group_id) %>%
#     filter(abundance != is.nan(abundance)) %>%
#     filter(monthly_time_step == min(monthly_time_step)) %>%
#     mutate(first_appearance = monthly_time_step,
#            beginning = ifelse(first_appearance <= 12006,
#                               TRUE, FALSE)) %>% 
#     dplyr::select(group_id, first_appearance, beginning)
#   
#   scenario_abundance_clean_5[[i]] <- temp %>% 
#     merge(scenario_averaged_5[[i]]) %>% 
#     tidylog::filter(beginning == TRUE)
#   
#   rm(temp)
#   
# }
# 
# vis_dat(scenario_abundance_clean_5[[1]][1:50,])
# vis_dat(scenario_abundance_clean_5[[1]][1:50,])
# vis_dat(scenario_abundance_clean_5[[1]][1:50,])
# vis_dat(scenario_abundance_clean_5[[1]][1:50,])
# 
# # Check we have the same number of species at beginning and end
# 
# x1 <- scenario_abundance_clean_5[[1]] %>% filter(monthly_time_step == 12006)
# x300 <- scenario_abundance_clean_5[[1]] %>% filter(monthly_time_step == max(monthly_time_step))
# 
# unique(x1$group_id) == unique(x300$group_id)
# x1x300 <- cbind(x1$group_id,x300$group_id)
# 
# # Remove false extinctions from averaged vals ----
# 
# scenario_remove_false_extinctions_2_5 <- list()
# 
# # for (i in seq_along(scenario_averaged)) {
# #   
# #   # Identify the last time-step a species has abundance values > 0
# #   temp2 <- scenario_averaged[[i]] %>% 
# 
# for (i in seq_along(scenario_abundance_clean_5)) {
#   
#   # Identify the last time-step a species has abundance values > 0
#   temp2 <- scenario_abundance_clean_5[[i]] %>%
#     group_by(group_id) %>% 
#     slice(1:3600) %>% 
#     filter(abundance > 0) %>% 
#     dplyr::select(group_id, monthly_time_step, abundance) %>% 
#     filter(monthly_time_step == max(monthly_time_step)) %>% 
#     dplyr::select(group_id, monthly_time_step) %>% 
#     rename(last_abundance = monthly_time_step) %>% 
#     ungroup()
#   
#   dim(temp2)
#   
#   # Add the year of last positive abundance number as a column to the data    
#   temp3 <- scenario_abundance_clean_5[[i]] %>% 
#     #remove existing last abundance val bc it's leftover from before averaging
#     # dplyr::select(-last_abundance) %>% 
#     merge(temp2, by = "group_id", all = TRUE)
#   
#   dim(temp3)
#   
#   rm(temp2)
#   
#   # Use the last positive abundance year and current abundance value to determine
#   # if a zero abundance is a true extinction or just a missing value (false extinction)
#   scenario_remove_false_extinctions_2_5[[i]] <- temp3 %>%
#     group_by(group_id) %>%
#     slice(1:3600) %>% 
#     mutate(true_extinction = ifelse(abundance == 0 &
#                                       monthly_time_step < last_abundance,
#                                     "false extinction",
#                                     ifelse(abundance > 0 &
#                                              monthly_time_step < last_abundance,
#                                            "not extinct",
#                                            ifelse(abundance == 0 &
#                                                     monthly_time_step >= last_abundance,
#                                                   "true extinction", "not extinct")))) %>%
#     mutate(abundance = ifelse(true_extinction == "false extinction",
#                               NA, abundance)) %>%
#     group_by(group_id) %>%
#     arrange(group_id, monthly_time_step)  %>% 
#     ungroup()
#   
#   rm(temp3)
# }
# 
# vis_dat(scenario_remove_false_extinctions_2_5[[1]][1:50,])
# vis_dat(scenario_remove_false_extinctions_2_5[[2]][1:50,])
# vis_dat(scenario_remove_false_extinctions_2_5[[3]][1:50,])
# vis_dat(scenario_remove_false_extinctions_2_5[[4]][1:50,])
# 
# x1 <- scenario_remove_false_extinctions_2_5[[1]] %>% 
#   filter(monthly_time_step == 12006)
# 
# x300 <- scenario_remove_false_extinctions_2_5[[1]] %>% 
#   filter(monthly_time_step == max(monthly_time_step))
# 
# dim(x1)[1] == dim(x300)[1]
# 
# baseline3 <- scenario_remove_false_extinctions_2_5[[1]]
# vis_dat(baseline3)
# 
# # # Smooth abundance ----
# 
# ave_window <- 12 * 5
# 
# scenario_smoothed_abundance_5 <- list()
# 
# for (i in seq_along(scenario_remove_false_extinctions_2_5)) {
#   
#   group_list <- split(scenario_remove_false_extinctions_2_5[[i]], 
#                       scenario_remove_false_extinctions_2_5[[i]]$group_id)
#   
#   group_smoothed_abundance <- list()
#   
#   for (j in seq_along(group_list)) {
#     
#     group_df <- group_list[[j]]
#     
#     # check if the group has any abundance values, make it null if not
#     
#     if (sum(group_df$abundance, na.rm = TRUE) == 0) {
#       
#       group_smoothed_abundance[[j]] <- NULL
#       
#     } else {
#       
#       group_smoothed_abundance[[j]] <- group_df %>%
#         arrange(monthly_time_step) %>%
#         mutate(ave_abundance = rollapply(abundance,
#                                          ave_window,
#                                          mean,
#                                          na.rm = TRUE,
#                                          partial = TRUE,
#                                          align = "left"),
#                ave_abundance = ifelse(ave_abundance < 1,
#                                       0, ave_abundance))
#       
#       print(j)
#       
#     }
#   }
#   
#   scenario_smoothed_abundance_5[[i]]  <- do.call(rbind,
#                                                group_smoothed_abundance)
#   
#   rm(group_smoothed_abundance, group_df)
# }
# 
# x1 <- scenario_smoothed_abundance_5[[1]] %>% filter(monthly_time_step == 12006)
# x300 <- scenario_smoothed_abundance_5[[1]] %>% filter(monthly_time_step == max(monthly_time_step))
# xx1 <- x1$group_id
# xx300 <- x300$group_id
# 
# unique(x1$group_id) == unique(x300$group_id)
# 
# vis_dat( scenario_smoothed_abundance_5[[1]][1:50,])
# vis_dat( scenario_smoothed_abundance_5[[2]][1:50,])
# vis_dat( scenario_smoothed_abundance_5[[3]][1:50,])
# vis_dat(scenario_remove_false_extinctions_2_5[[3]][1:50,])
# vis_dat( scenario_smoothed_abundance_5[[4]][1:50,])
# vis_dat(scenario_remove_false_extinctions_2_5[[4]][1:50,])
# 
# 
# saveRDS(scenario_smoothed_abundance_5,
#         file.path(indicator_outputs_folder,
#                   paste(today, "scenario_smoothed_abundance_5.rds", sep = "_")))


  
# CHECKPOINT - CAN LOAD PROCESSED DATA FROM HERE ----

input_date <- "2021-09-22"

gen_folder <- "N:/Quantitative-Ecology/Indicators-Project//Serengeti/Outputs_from_indicator_code/Indicator_outputs/general"

# Load data

# scenario_ab_gl_formatted <- readRDS(file.path(gen_folder,input_date,
#                                               "formatted_abundance_1.rds"))
# 
# scenario_decomposed <- readRDS(file.path(gen_folder,input_date,
#                                          "decomposed_abundance_2.rds"))
# 
# scenario_averaged <- readRDS(file.path(gen_folder,input_date,
#                                           "averaged_abundance_3.rds"))
# 
# scenario_filtered <- readRDS(file.path(gen_folder,input_date,
#                                        paste(input_date,
#                                          "weird_groups_removed_abundance_4.rds"))
# 
# scenario_smoothed_abundance <- readRDS(file.path(gen_folder,input_date,
#                                           "smoothed_abundance_5.rds"))
# 
# scenario_auto_smoothed_abundance <- readRDS(file.path(gen_folder,input_date,
#                                                  "smoothed_auto_abundance_5.rds"))
# 
# scenario_harvested <- readRDS(file.path(gen_folder,input_date,
#                                              "scenario_density_harvested.rds"))
# 
# scenario_smoothed_5_reps <- readRDS(file.path(gen_folder,input_date,
#                                         "scenario_smoothed_abundance_5.rds"))

# RED LIST INDEX ANNUAL ----

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
    
    # groups_annual[[j]] <- scenario_groups[[j]] %>% 
    #   slice(-n()) %>% 
    #   slice(which(row_number() %% interval == 0))  
    
     groups_annual[[j]] <- scenario_groups[[j]] %>% 
                           group_by(annual_time_step) %>% 
       filter(ave_abundance == first(na.omit(ave_abundance))) %>% 
       dplyr::select(-monthly_time_step, - true_extinction,
                     -abundance, -sd, - error, -lower_ci, -upper_ci) %>% 
       distinct(.) %>% 
       ungroup(.) %>% 
       group_by(group_id, annual_time_step) %>% 
       slice(1) %>% 
       ungroup(.)
    
    }
  
  groups_annual_df <- do.call(rbind, groups_annual)
  
  scenario_annual[[i]] <- groups_annual_df
  
  rm(groups_annual, groups_annual_df)
  
}

check <- scenario_annual[[1]] %>% filter(group_id == "11.63")
head(check)
dim(check)
any(duplicated(check$annual_time_step))

ggplot() +
  geom_point(data = check, aes(x = annual_time_step, y = ave_abundance)) 

x <- scenario_annual[[1]]
x1 <- x %>% filter(annual_time_step == 1)
x300 <- x %>% filter(annual_time_step == 300)

# Sample autotrophs

scenario_auto_annual <- list()

for (i in seq_along(scenario_auto_smoothed_abundance)) {
  
  scenario_auto_groups <- split(scenario_auto_smoothed_abundance[[i]], 
                                scenario_auto_smoothed_abundance[[i]]$group_id)
  
  groups_auto_annual <- list()
  
  for (j in seq_along(scenario_auto_groups)) {
    
    groups_auto_annual[[j]] <- scenario_auto_groups[[j]] %>% 
      group_by(annual_time_step) %>%
      arrange(monthly_time_step) %>% 
      slice(1)
    
    # groups_auto_annual[[j]] <- scenario_auto_groups[[j]] %>% 
    #   slice(which(row_number() %% interval == 0)) %>%
    #   mutate(annual_time_step = seq(1,max_timestep,1))
    
  }
  
  groups_annual_auto_df <- do.call(rbind, groups_auto_annual)
  
  scenario_auto_annual[[i]] <- groups_annual_auto_df
  
  rm(groups_auto_annual, groups_annual_auto_df)

}

check_auto <- scenario_auto_annual[[1]]
head(check_auto)
length(scenario_auto_annual)

saveRDS(scenario_annual,
        file.path(general_indicator_outputs_folder,
                   "scenario_annual_6.rds"))

saveRDS(scenario_auto_annual,
        file.path(general_indicator_outputs_folder,
                  "scenario_auto_annual_6.rds"))

## Referring to the thresholds quote under Criterion A, Reason 1 (declines
## are the result of reversible pressures) according to:
## https://portals.iucn.org/library/sites/library/files/documents/RL-2001-001-2nd.pdf


# * Assign Red List Categories ----

scenario_red_list_data <- list()

for (i in seq_along(scenario_annual)) {
  
  # Get replicate data for a single scenario
  
  print(paste("Processing scenario", scenarios[[i]], sep = " "))
  
  # Get one scenario dataframe
  
  status_inputs <- split(scenario_annual[[i]], 
                         scenario_annual[[i]]$group_id)
  
  # Make a list to hold output for each individual massbin-func-group (ie virtual spp)
  
  group_red_list_data <- list()
  
  for (j in seq_along(status_inputs)) {
    
    print(paste("Processing group", names(status_inputs)[[j]], sep = " "))
    
    group_red_list_data[[j]] <- status_inputs[[j]] %>%
      # calculate the difference in abundance over 10 yrs or 3 generation lengths
      # (specified by 'timeframe' column). Its okay to take the first value of 
      # timeframe bc the dataframe is grouped by group_id, and timeframe only changes
      # between and not within group_ids
      mutate(gen_length_3 = round(generation_length_yrs * 3)) %>% 
      mutate(timeframe = ifelse(gen_length_3 > 10,
                                gen_length_3, 10)) %>%
      group_by(group_id) %>%
      arrange(annual_time_step) %>%
      mutate(diff = (ave_abundance - dplyr::lag(ave_abundance, timeframe[1]))) %>%
      # Using the formula from p 35 (Complex patterns of decline) Guidelines 
      # for Using the IUCN Red List Categories and Criteria v14 August 2019 
      mutate(decline = 1 - ave_abundance/dplyr::lag(ave_abundance, timeframe[1])) %>%
      mutate(decline = ifelse(ave_abundance == 0, NA, decline)) %>% 
      # calculate the rate of change
      # assign red list risk status based on decline 
      # Using the thresholds from p 16 Categories A2 - A4 Guidelines 
      # for Using the IUCN Red List Categories and Criteria v14 August 2019
      mutate(rl_status = ifelse(decline < 0.20, "LC",
                         ifelse(decline >= 0.20 & decline < 0.30, "NT", 
                         ifelse(decline >= 0.30 & decline < 0.50, "VU",
                         ifelse(decline >= 0.50 & decline < 0.80, "EN",
                         ifelse(decline >= 0.80, "CR",
                         ifelse(is.na(decline), "EX", "TBD"))))))) %>%
      arrange(group_id, annual_time_step) %>%
      # Replace all non-ex status with ex after first occurrence 
      mutate(rl_status = ifelse(ave_abundance == 0, "EX", rl_status),
             rl_status = ifelse(is.na(rl_status), lag(rl_status, 1),
                                rl_status)) %>% 
      mutate(extinct = ifelse(rl_status == "EX", 1, 0)) #%>% 
      #group_by(group_id) #%>% 
      #filter(annual_time_step != c(296, 297, 298, 299, 300)) # remove last 5 timesteps 
    
  }
  
  scenario_red_list_df <- do.call(rbind, group_red_list_data)
  
  scenario_red_list_data[[i]] <- scenario_red_list_df
  
}


# Have a quick look at the outputs

rli_inputs <- scenario_red_list_data[[1]]
tail(rli_inputs)
dim(rli_inputs)

# write.csv(rli_inputs, file.path(indicator_outputs_folder, 
#                                 "rli_input_example_averaged_reps.csv"))

rli_inputs_group <- rli_inputs %>% filter(group_id == "10.68")

ggplot(data = rli_inputs_group) +
  geom_path(aes(x = annual_time_step, y = ave_abundance)) +
  theme(legend.position = "none") +
  geom_text(aes(x = annual_time_step, y = ave_abundance, label = rl_status))


# * Take coarser sample ----
## Heterotrophs

# We actually odn't need to take a coarser sample here for annual calcs so 
# just rename the object

scenario_redlist_data_sampled <- scenario_red_list_data

saveRDS(scenario_redlist_data_sampled,
        file.path(general_indicator_outputs_folder, 
                  "scenario_redlist_data_annual_7.rds"))

## Autotrophs

scenario_auto_sampled <- scenario_auto_annual

#####

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
# 
# x <- scenario_redlist_data_sampled[[1]]
# x1 <- x %>% filter(annual_time_step == 1)
# x300 <- x %>% filter(annual_time_step == 300)
# 
# # * Calculate RLI ----
# 
# # RLI by individual functional groups
# 
# scenario_fg_rli_outputs <- list()
# 
# for (i in seq_along(scenario_redlist_data_sampled)) {
#   
#   scenario_fg_rli_outputs[[i]] <- calculate_red_list_index(
#     scenario_redlist_data_sampled[[i]], numboots, ci = FALSE) %>%
#     mutate(scenario = scenarios[[i]]) 
#   
# }
# 
# x <- scenario_fg_rli_outputs[[1]]
# head(x)
# 
# # Mean RLI aggregated across groups
# 
# scenario_rli_outputs <- list()
# 
# for (i in seq_along(scenario_fg_rli_outputs)) {
#   
#   if ("ci_lower" %in% names(scenario_fg_rli_outputs[[i]])) {
#     
#     scenario_rli_outputs[[i]] <- scenario_fg_rli_outputs[[i]] %>%
#       group_by(annual_time_step) %>%
#       summarise(indicator_score = mean(indicator_score),
#                 ci_lower = mean(ci_lower),
#                 ci_upper = mean(ci_upper)) %>%
#       mutate(indicator = "RLI",
#              replicate = j)
#   } else {
#     
#     scenario_rli_outputs[[i]] <- scenario_fg_rli_outputs[[i]] %>%
#       group_by(annual_time_step) %>%
#       summarise(indicator_score = mean(indicator_score)) %>%
#       mutate(indicator = "RLI",
#              replicate = j)
#   }
#   
# }
# 
# head(scenario_rli_outputs)[[1]]
# 
# # * Get the mean abundance of all groups in each time step
# 
# scenario_mean_abundance <- list()
# 
# for (i in seq_along(scenario_redlist_data_sampled)) {
# 
#   scenario_mean_abundance[[i]] <- scenario_redlist_data_sampled[[i]] %>% 
#   ungroup() %>% 
#   group_by(annual_time_step) %>% 
#   mutate(mean_group_abundance = mean(ave_abundance, na.rm = TRUE)) %>% 
#   ungroup() %>% 
#   dplyr::select(annual_time_step, mean_group_abundance) %>% 
#   distinct(.)
# 
#   scenario_mean_abundance[[i]]$standardised_abundance <- range01(scenario_mean_abundance[[i]]$mean_group_abundance)  
# 
# }
# 
# x <- scenario_mean_abundance[[1]]
# x <- scenario_redlist_data_sampled[[1]]
# 
# # * Plot RLI ----
# 
# ## By functional group
# 
# scenario_fg_rli_plots <- list()
# 
# for (i in seq_along(scenario_fg_rli_outputs)) {
#   
#   scenario_fg_rli_plots[[i]] <-  plot_red_list_index_by_group(
#     scenario_fg_rli_outputs[[i]],
#     impact_start,
#     impact_end,
#     ci = FALSE) 
#   
#   ggsave(file.path(rli_plots_folder, paste(today, scenarios[[i]], 
#                                            "RLI_by_functional_group_reps_averaged.png",
#                                            sep = "_")),
#          scenario_fg_rli_plots[[i]],  device = "png")
#   
# }
# 
# scenario_fg_rli_plots[[4]]
# 
# 
# # RLI with all functional groups aggregated
# # i.e. mean of each 'taxa' RLI as per Butchart et al (2010) 'Indicators of
# # recent declines'
# 
# scenario_rli_plots <- list()
# 
# for (i in seq_along(scenario_rli_outputs)) {
#   
#   
#   scenario_rli_plots[[i]] <- plot_red_list_index(scenario_rli_outputs[[i]],
#                                                  impact_start, 
#                                                  impact_end,
#                                                  ci = FALSE) +
#                              labs(title = "Annual sampling, RLI by func. groups then averaged")
#   
#   ggsave(file.path(rli_plots_folder, paste(today, scenarios[[i]],
#                                            "RLI_aggregated_averaged.png",
#                                            sep = "_")),
#          scenario_rli_plots[[i]],  device = "png")
#   
# }
# 
# scenario_rli_plots[[1]]
# scenario_rli_plots[[2]]
# scenario_rli_plots[[3]]
# scenario_rli_plots[[4]]
# 
# RED LIST INDEX ANNUAL V2 ----

## Lumping all species together to calculate RLI first up, rather than calculating
## by functional group then taking the mean

# * Calculate RLI ----

# RLI across all groups

scenario_all_rli_outputs <- list()

for (i in seq_along(scenario_redlist_data_sampled)) {

  scenario_all_rli_outputs[[i]] <- calculate_red_list_index2(
    scenario_redlist_data_sampled[[i]], numboots, ci = FALSE) %>%
    mutate(scenario = scenarios[[i]])

}

x <- scenario_all_rli_outputs[[1]]
head(x)

# * Plot RLI ----

# Using V2

scenario_all_rli_plots <- list()

for (i in seq_along(scenario_all_rli_outputs)) {

  scenario_all_rli_plots[[i]] <- plot_red_list_index(scenario_all_rli_outputs[[i]],
                                                     impact_start,
                                                     impact_end,
                                                     ci = FALSE) +
    labs(title = "Annual sampling, RLI all groups")

  ggsave(file.path(rli_plots_folder, paste(today, scenarios[[i]],
                                           "RLI_aggregated_averaged_annual_v2.png",
                                           sep = "_")),
         scenario_all_rli_plots[[i]],  device = "png")

}

scenario_all_rli_plots[[1]]
scenario_all_rli_plots[[2]]
scenario_all_rli_plots[[3]]
scenario_all_rli_plots[[4]]
# 
# # RED LIST INDEX ANNUAL V3 ----
# 
# # Combining all species as above, but only including species larger than 10kg
# 
# # * Calculate RLI ----
# 
# # RLI with only large species
# 
# scenario_large_spp_rli_outputs <- list()
# 
# for (i in seq_along(scenario_redlist_data_sampled)) {
#   
#   large_spp <- scenario_redlist_data_sampled[[i]] %>% 
#               filter(mass_lower > 10000)
#   
#   scenario_large_spp_rli_outputs[[i]] <- calculate_red_list_index2(
#     large_spp, numboots, ci = FALSE) %>%
#     mutate(scenario = scenarios[[i]]) 
#   
#   rm(large_spp)
#   
# }
# 
# x <- scenario_large_spp_rli_outputs[[1]]
# head(x)
# 
# # * Plot RLI ----
# 
# # Using V2
# 
# scenario_large_spp_rli_plots <- list()
# 
# for (i in seq_along(scenario_large_spp_rli_outputs)) {
#   
#   scenario_large_spp_rli_plots[[i]] <- plot_red_list_index(scenario_large_spp_rli_outputs[[i]],
#                                                      impact_start, 
#                                                      impact_end,
#                                                      ci = FALSE) +
#     labs(title = "Annual sampling, RLI of all spp > 10kg")
#   
#   ggsave(file.path(rli_plots_folder, paste(today, scenarios[[i]],
#                                            "RLI_aggregated_averaged_annual_lgspp_v3.png",
#                                            sep = "_")),
#          scenario_large_spp_rli_plots[[i]],  device = "png")
#   
# }
# 
# scenario_large_spp_rli_plots[[1]]
# scenario_large_spp_rli_plots[[2]]
# scenario_large_spp_rli_plots[[3]]
# scenario_large_spp_rli_plots[[4]]


# RED LIST INDEX 5 YRS ----

# * Take coarser sample ----
## Heterotrophs

sample_interval <- 5 # interval between samples in years
sample_max_timestep <- 300/sample_interval

scenario_redlist_data_sampled_5yr <- list()

for (i in seq_along(scenario_red_list_data)) {
  
  # Sample
  
  set <- seq(5,295,5)
 
  scenario_redlist_data_sampled_5yr[[i]] <- scenario_red_list_data[[i]] %>% 
    group_by(group_id) %>%
    #Calculate how many timesteps of data available for that species
    mutate(n_timesteps = n()) %>% 
    #Add a 5 yr interval column for sampling
    mutate(annual_time_step = rep(set, each = sample_interval,
                          length.out = n_timesteps[1])) %>% 
    # Group by species and interval
    group_by(group_id, annual_time_step) %>%
    # Take the first observation
    slice(1)
  
}

test <- scenario_redlist_data_sampled_5yr[[1]]
test_group <- test %>% filter(group_id == "10.32")
dim(test_group) 


## Autotrophs

scenario_auto_sampled_5yr <- list()

for (i in seq_along(scenario_auto_annual)) {
  
  # Sample
  scenario_auto_sampled_5yr[[i]] <- scenario_auto_annual[[i]] %>% 
    group_by(group_id) %>% 
    slice(which(row_number() %% sample_interval == 0)) 
  
}

test <- scenario_auto_sampled_5yr[[1]]
test_group <- test %>% filter(group_id == "autotrophs")
dim(test_group) # should have 300 rows

saveRDS(scenario_auto_sampled_5yr,
        file.path(general_indicator_outputs_folder, 
                  "scenario_redlist_data_auto_annual_8.rds"))

# # * Get harvested groups only ----
# 
# scenario_harvested_groups_5yr <- list()
# 
# for (i in seq_along(scenario_redlist_data_sampled_5yr)) {
#   
#   scenario <- scenarios[[i]]
#   
#   if (scenario == "000_Baseline") {
#     
#     # No disturbance so just get whatever groups
#     
#     harvested <- scenario_redlist_data_sampled_5yr[[i]] %>% 
#       mutate(scenario = scenarios[[i]]) %>% 
#       group_by(annual_time_step) %>% 
#       mutate(indicator_score = sum(ave_abundance, na.rm = TRUE),
#              indicator = "total abundance harvested",
#              ci_lower = NA,
#              ci_upper = NA,
#              replicate = NA) %>% 
#       ungroup(.) %>%                                 
#       dplyr::select(annual_time_step, 
#                     indicator_score,
#                     ci_lower,
#                     ci_upper,
#                     indicator,
#                     replicate,
#                     scenario)%>% 
#       distinct(.)
#     
#   } else if (scenario == "100_Land_use") {
#     
#     # Get the autotroph abundance
#     
#     harvested <- scenario_auto_sampled_5yr[[i]] %>% 
#       filter(group_id == "autotrophs") %>% 
#       mutate(scenario = scenarios[[i]]) %>%
#       group_by(annual_time_step) %>% 
#       #mutate(ab_scaled = range01(abundance),
#       mutate(indicator_score = sum(ave_abundance, na.rm = TRUE),
#              indicator = "total abundance harvested",
#              ci_lower = NA,
#              ci_upper = NA,
#              replicate = NA) %>% 
#       ungroup(.) %>% 
#       dplyr::select(annual_time_step, 
#                     indicator_score,
#                     ci_lower,
#                     ci_upper,
#                     indicator,
#                     replicate,
#                     scenario) %>% 
#       distinct(.)
#     
#   } else if (scenario == "200_Harvesting_carnivores") {
#     
#     # Get the carnivorous endotherms 100 - 200kg
#     # Harvest doesn't seem to affect the ectotherms???
#     
#     harvested <- scenario_redlist_data_sampled_5yr[[i]] %>% 
#       filter(functional_group_name == "carnivore endotherm" &
#                group_id == "11.67"|
#                functional_group_name == "carnivore endotherm" &
#                group_id == "11.68")  %>% 
#       mutate(scenario = scenarios[[i]]) %>% 
#       ungroup(.) %>% 
#       dplyr::select(annual_time_step, ave_abundance, scenario) %>% 
#       distinct(.) %>% 
#       mutate(ab_scaled = range01(ave_abundance)) %>% 
#       group_by(annual_time_step) %>% 
#       #mutate(indicator_score = sum(ab_scaled, na.rm = TRUE),
#       mutate(indicator_score = sum(ave_abundance, na.rm = TRUE),
#              indicator = "total abundance harvested",
#              ci_lower = NA,
#              ci_upper = NA,
#              replicate = NA) %>% 
#       ungroup(.) %>% 
#       dplyr::select(annual_time_step, 
#                     indicator_score,
#                     ci_lower,
#                     ci_upper,
#                     indicator,
#                     replicate,
#                     scenario) %>% 
#       distinct(.)
#     
#     
#   } else if (scenario == "300_Harvesting_herbivores") {
#     
#     # Get the herbivorous endotherms 100 - 200kg
#     # Harvest doesn't seem to affect the ectotherms???
#     
#     harvested <- scenario_redlist_data_sampled_5yr[[i]] %>% 
#       filter(functional_group_name == "herbivore endotherm" &
#                group_id == "10.67"|
#                functional_group_name == "herbivore endotherm" &
#                group_id == "10.68")   %>% 
#       mutate(scenario = scenarios[[i]]) %>% 
#       ungroup(.) %>% 
#       dplyr::select(annual_time_step, ave_abundance, scenario) %>% 
#       distinct(.) %>% 
#       #mutate(ab_scaled = range01(abundance)) %>% 
#       group_by(annual_time_step) %>% 
#       mutate(indicator_score = sum(ave_abundance, na.rm = TRUE),
#              indicator = "total abundance harvested",
#              ci_lower = NA,
#              ci_upper = NA,
#              replicate = NA) %>% 
#       ungroup(.) %>% 
#       dplyr::select(annual_time_step, 
#                     indicator_score,
#                     ci_lower,
#                     ci_upper,
#                     indicator,
#                     replicate,
#                     scenario) %>% 
#       distinct(.)
#     
#   } 
#   
#   scenario_harvested_groups_5yr[[i]] <- harvested
# }

# any(is.null(scenario_harvested_groups_5yr))
# 
# harv5yr <- scenario_harvested_groups_5yr[[1]]
# head(harv5yr)
# dim(harv5yr)

# # * Get proportion extinct ----
# 
# test2 <- scenario_redlist_data_sampled[[1]] %>% 
#          group_by(annual_time_step) %>% 
#          summarise(x = n())
# 
# test <- scenario_redlist_data_sampled_5yr[[1]] %>% 
#         group_by(annual_time_step) %>% 
#         summarise(x = n())
# 
# 
# scenario_extinctions_5yr <- list()
# scenario_rl_status_plots_5yr <- list()
# 
# for (i in seq_along(scenario_redlist_data_sampled_5yr)) {
#   
#   scenario_extinctions_5yr[[i]] <- scenario_redlist_data_sampled_5yr[[i]] %>% 
#     group_by(annual_time_step, rl_status) %>% 
#     filter(rl_status != is.na(rl_status)) %>% 
#     summarise(test = n())
#   
#   scenario_rl_status_plots_5yr[[i]] <- ggplot(scenario_extinctions_5yr[[i]], 
#                                           aes(x = annual_time_step, 
#                                               y = test, 
#                                               fill = rl_status)) +
#     geom_bar(position = "stack", stat = "identity") +
#     scale_fill_viridis_d()
#   
#   ggsave(file.path(rli_plots_folder, paste(today, scenarios[[i]], 
#                                            "red_list_status.png",
#                                            sep = "_")),
#          scenario_rl_status_plots_5yr[[i]],  device = "png")
#   
# }
# 
# scenario_rl_status_plots_5yr[[1]]
# scenario_rl_status_plots_5yr[[2]]
# scenario_rl_status_plots_5yr[[3]]
# scenario_rl_status_plots_5yr[[4]]
# 
# # * Calculate RLI ----
# 
# # RLI by individual functional groups
# 
# scenario_fg_rli_outputs_5yr <- list()
# 
# for (i in seq_along(scenario_redlist_data_sampled_5yr)) {
#   
#   scenario_fg_rli_outputs_5yr[[i]] <- calculate_red_list_index(
#     scenario_redlist_data_sampled_5yr[[i]], numboots, ci = FALSE) %>%
#     mutate(scenario = scenarios[[i]]) 
#   
# }
# 
# x <- scenario_fg_rli_outputs_5yr[[4]]
# head(x)
# 
# # Mean RLI aggregated across groups
# 
# scenario_rli_outputs_5yr <- list()
# 
# for (i in seq_along(scenario_fg_rli_outputs_5yr)) {
#   
#   if ("ci_lower" %in% names(scenario_fg_rli_outputs_5yr[[i]])) {
#     
#     scenario_rli_outputs_5yr[[i]] <- scenario_fg_rli_outputs_5yr[[i]] %>%
#       group_by(annual_time_step) %>%
#       summarise(indicator_score = mean(indicator_score),
#                 ci_lower = mean(ci_lower),
#                 ci_upper = mean(ci_upper)) %>%
#       mutate(indicator = "RLI",
#              replicate = j)
#   } else {
#     
#     scenario_rli_outputs_5yr[[i]] <- scenario_fg_rli_outputs_5yr[[i]] %>%
#       group_by(annual_time_step) %>%
#       summarise(indicator_score = mean(indicator_score)) %>%
#       mutate(indicator = "RLI",
#              replicate = j)
#   }
#   
# }
# 
# head(scenario_rli_outputs_5yr)[[3]]
# 
# # * Plot RLI ----
# 
# ## By functional group
# 
# scenario_rli_plots_fg_5yr <- list()
# 
# for (i in seq_along(scenario_fg_rli_outputs_5yr)) {
#   
#   scenario_rli_plots_fg_5yr[[i]] <-  plot_red_list_index_by_group(
#     scenario_fg_rli_outputs_5yr[[i]],
#     impact_start,
#     impact_end,
#     ci = FALSE)
#   
#   ggsave(file.path(rli_plots_folder, paste(today, scenarios[[i]], 
#                                            "RLI_by_functional_group_reps_averaged_5yr.png",
#                                            sep = "_")),
#          scenario_rli_outputs_5yr[[i]],  device = "png")
#   
# }
# 
# scenario_rli_plots_fg_5yr[[1]]
# scenario_rli_plots_fg_5yr[[2]]
# scenario_rli_plots_fg_5yr[[3]]
# scenario_rli_plots_fg_5yr[[4]]
# 
# # RLI with all functional groups aggregated
# # i.e. mean of each 'taxa' RLI as per Butchart et al (2010) 'Indicators of
# # recent declines'
# 
# scenario_rli_plots_5yr <- list()
# 
# for (i in seq_along(scenario_rli_outputs_5yr)) {
#   
#   
#   scenario_rli_plots_5yr[[i]] <- plot_red_list_index(scenario_rli_outputs_5yr[[i]],
#                                                  impact_start, 
#                                                  impact_end,
#                                                  ci = FALSE) +
#     labs(title = "5yr sampling, RLI by func. groups then averaged")
#   
#   ggsave(file.path(rli_plots_folder, paste(today, scenarios[[i]], 
#                                            "RLI_aggregated_averaged_5yr.png",
#                                            sep = "_")),
#          scenario_rli_plots_5yr[[i]],  device = "png")
#   
# }
# 
# scenario_rli_plots_5yr[[1]]
# scenario_rli_plots_5yr[[2]]
# scenario_rli_plots_5yr[[3]]
# scenario_rli_plots_5yr[[4]]


# RED LIST INDEX 5 YRS V2 ----

## Assign the red list categories after sampling instead of before

## Referring to the thresholds quote under Criterion A, Reason 1 (declines
## are the result of reversible pressures) according to:
## https://portals.iucn.org/library/sites/library/files/documents/RL-2001-001-2nd.pdf


# * Assign Red List Categories ----

scenario_red_list_data_v2 <- list()

#for (i in seq_along(scenario_ab_gl_formatted)) {
for (i in seq_along(scenario_redlist_data_sampled_5yr)) {
  
  # Get replicate data for a single scenario
  
  print(paste("Processing scenario", scenarios[[i]], sep = " "))
  
  # Split by functional group, because we calculate RLI for different
  # functional groups then aggregate later (as per Butchart etal 2010),
  # except we are using functional groups as proxies for taxa (eg mammals, birds, 
  # reptiles) used in real world RLI calcs
  
  status_inputs <- split(scenario_redlist_data_sampled_5yr[[i]], 
                         scenario_redlist_data_sampled_5yr[[i]]$group_id)
  
  # Make a list to hold output for each individual massbin-func-group (ie virtual spp)
  
  group_red_list_data <- list()
  
  for (j in seq_along(status_inputs)) {
    
    print(paste("Processing group", names(status_inputs)[[j]], sep = " "))
    
    group_red_list_data[[j]] <- status_inputs[[j]] %>%
      # calculate the difference in abundance over 10 yrs or 3 generation lengths
      # (specified by 'timeframe' column). Its okay to take the first value of 
      # timeframe bc the dataframe is grouped by group_id, and timeframe only changes
      # between and not within group_ids
      mutate(gen_length_3 = round(generation_length_yrs * 3)) %>% 
      mutate(timeframe = ifelse(gen_length_3 > 10,
                                gen_length_3, 10)) %>%
      # Convert the yearly timeframe to deal with 5 yr intervals
      mutate(timeframe = round(timeframe/5)) %>% 
      group_by(group_id) %>%
      arrange(annual_time_step) %>%
      mutate(diff = (ave_abundance - dplyr::lag(ave_abundance, timeframe[1]))) %>%
      # Using the formula from p 35 (Complex patterns of decline) Guidelines 
      # for Using the IUCN Red List Categories and Criteria v14 August 2019 
      mutate(decline = 1 - ave_abundance/dplyr::lag(ave_abundance, timeframe[1])) %>%
      mutate(decline = ifelse(ave_abundance == 0, NA, decline)) %>% 
      # calculate the rate of change
      # assign red list risk status based on decline 
      # Using the thresholds from p 16 Categories A2 - A4 Guidelines 
      # for Using the IUCN Red List Categories and Criteria v14 August 2019
      mutate(rl_status = ifelse(decline < 0.20, "LC",
                                ifelse(decline >= 0.20 & decline < 0.30, "NT", 
                                       ifelse(decline >= 0.30 & decline < 0.50, "VU",
                                              ifelse(decline >= 0.50 & decline < 0.80, "EN",
                                                     ifelse(decline >= 0.80, "CR",
                                                            ifelse(is.na(decline), "EX", "TBD"))))))) %>%
      arrange(group_id, annual_time_step) %>%
      # Replace all non-ex status with ex after first occurrence 
      mutate(rl_status = ifelse(ave_abundance == 0, "EX", rl_status),
             rl_status = ifelse(is.na(rl_status), lag(rl_status, 1),
                                rl_status)) %>% 
      mutate(extinct = ifelse(rl_status == "EX", 1, 0)) #%>% 
    #group_by(group_id) #%>% 
    #filter(annual_time_step != c(296, 297, 298, 299, 300)) # remove last 5 timesteps 
    
  }
  
  scenario_red_list_df <- do.call(rbind, group_red_list_data)
  
  scenario_red_list_data_v2[[i]] <- scenario_red_list_df
  
  rm(scenario_red_list_df, group_red_list_data, status_inputs)
  
  # Save the inputs
  
}

saveRDS(scenario_red_list_data_v2,
        file.path(general_indicator_outputs_folder, 
                  "scenario_redlist_data_5yr_8.rds"))

# Have a quick look at the outputs

rli_inputs <- scenario_red_list_data_v2[[4]]
tail(rli_inputs)
dim(rli_inputs)

rli_inputs_group <- rli_inputs %>% filter(group_id == "10.59")

ggplot(data = rli_inputs_group) +
  geom_path(aes(x = annual_time_step, y = ave_abundance)) +
  theme(legend.position = "none") +
  geom_text(aes(x = annual_time_step, y = ave_abundance, label = rl_status))

# * Get proportion extinct ----

scenario_extinctions <- list()
scenario_rl_status_plots <- list()

for (i in seq_along(scenario_red_list_data_v2)) {
  
  scenario_extinctions[[i]] <- scenario_red_list_data_v2[[i]] %>% 
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

# RLI across all groups

scenario_all_rli_outputs_v2 <- list()

for (i in seq_along(scenario_red_list_data_v2)) {
  
  scenario_all_rli_outputs_v2[[i]] <- calculate_red_list_index2(
    scenario_red_list_data_v2[[i]], numboots, ci = FALSE) %>%
    mutate(scenario = scenarios[[i]],
           indicator = "RLI 5y all spp") 
  
}

x <- scenario_all_rli_outputs_v2[[3]]
head(x)

# * Plot RLI ----

# Using V2

scenario_all_rli_plots_v2 <- list()

for (i in seq_along(scenario_all_rli_outputs_v2)) {
  
  scenario_all_rli_plots_v2[[i]] <- plot_red_list_index(scenario_all_rli_outputs_v2[[i]],
                                                 impact_start, 
                                                 impact_end,
                                                 ci = FALSE)
  
  ggsave(file.path(rli_plots_folder, paste(today, scenarios[[i]],
                                           "RLI_aggregated_averaged_5yr2.png",
                                           sep = "_")),
         scenario_all_rli_plots_v2[[i]],  device = "png")
  
}

scenario_all_rli_plots_v2[[1]]
scenario_all_rli_plots_v2[[2]]
scenario_all_rli_plots_v2[[3]]
scenario_all_rli_plots_v2[[4]]

# RED LIST INDEX 5YRS V3 ----

# * Calculate RLI ----

# RLI with only large species

scenario_large_spp_rli_outputs_5yr <- list()

for (i in seq_along(scenario_redlist_data_sampled_5yr)) {

  large_spp <- scenario_redlist_data_sampled_5yr[[i]] %>%
    tidylog::filter(mass_lower > 10000)

  scenario_large_spp_rli_outputs_5yr[[i]] <- calculate_red_list_index2(
    large_spp, numboots, ci = FALSE) %>%
    mutate(scenario = scenarios[[i]],
           indicator = "RLI 5y large spp")

  rm(large_spp)

}

x <- scenario_large_spp_rli_outputs_5yr[[2]]
head(x)

# * Plot RLI ----

# Using V2

scenario_large_spp_rli_plots_5yr <- list()

for (i in seq_along(scenario_large_spp_rli_outputs_5yr)) {

  scenario_large_spp_rli_plots_5yr[[i]] <- plot_red_list_index(scenario_large_spp_rli_outputs_5yr[[i]],
                                                           impact_start,
                                                           impact_end,
                                                           ci = FALSE)

  ggsave(file.path(rli_plots_folder, paste(today, scenarios[[i]],
                                           "RLI_aggregated_averaged_5yr_lgspp_v3.png",
                                           sep = "_")),
         scenario_large_spp_rli_plots_5yr[[i]],  device = "png")

}

scenario_large_spp_rli_plots_5yr[[1]]
scenario_all_rli_plots_v2[[2]]
scenario_large_spp_rli_plots_5yr[[2]]
scenario_all_rli_plots_v2[[3]]
scenario_large_spp_rli_plots_5yr[[3]]
scenario_all_rli_plots_v2[[4]]
scenario_large_spp_rli_plots_5yr[[4]]


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
                  ave_abundance, functional_group_name) %>% 
    rename(abundance = ave_abundance)
  
  
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
    
  ) %>% 
    mutate(indicator = "LPI")
  
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

dev.off()

# # ** Functional groups ----
# 
# scenario_lpi_fg_outputs <- list()
# 
# for (i in seq_along(scenario_lpi_inputs)) {
#   
#   fg_data <- split(scenario_lpi_inputs[[i]], 
#                    scenario_lpi_inputs[[i]]$functional_group_name)
# 
#   fg_lpi_outputs <- list()
#   
#   for ( j in seq_along(fg_data)) {
#   
#     fg_lpi_outputs[[j]] <- calculate_living_planet_index(
#       
#       fg_data[[j]], start_time_step, ci = FALSE, numboots, 
#       fg_data[[j]]$functional_group_name[1])
#   }
#   
#   scenario_lpi_fg_outputs[[i]] <- do.call(rbind,fg_lpi_outputs)
#   
# }
# 
# head(scenario_lpi_fg_outputs[[1]])
# 
# # Plot fg
# 
# scenario_fg_lpi_plots <- list()
# 
# for (i in seq_along(scenario_lpi_fg_outputs)) {
#   
#   scenario_fg_lpi_plots[[i]] <- ggplot(data = scenario_lpi_fg_outputs[[i]],
#                                        aes(x = annual_time_step, 
#                                            y = indicator_score,
#                                            col = replicate)) +
#                                 geom_line() +
#                                 facet_wrap( ~ replicate)
#   
#   ggsave(file.path(lpi_plots_folder, paste(today, scenarios[[i]], 
#                                            "LPI_func_groups_reps_averaged.png",
#                                            sep = "_")),
#          scenario_fg_lpi_plots[[i]],  device = "png")   
#   
# }
# 
# scenario_fg_lpi_plots[[1]]
# scenario_fg_lpi_plots[[2]]
# scenario_fg_lpi_plots[[3]]
# scenario_fg_lpi_plots[[4]]
# 
# dev.off()
# 
# # ** Trophic groups ----
# 
# scenario_lpi_tg_outputs <- list()
# 
# for (i in seq_along(scenario_lpi_inputs)) {
#   
#   data <- scenario_lpi_inputs[[i]] %>% 
#           mutate(trophic_group = word(functional_group_name, 1))
#   
#   head(data)
#   
#   tg_data <- split(data, data$trophic_group)
#   
#   tg_lpi_outputs <- list()
#   
#   for ( j in seq_along(tg_data)) {
#     
#     tg_lpi_outputs[[j]] <- calculate_living_planet_index(
#       
#       tg_data[[j]], start_time_step, ci = FALSE, numboots, 
#       tg_data[[j]]$trophic_group[1])
#   }
#   
#   scenario_lpi_tg_outputs[[i]] <- do.call(rbind,tg_lpi_outputs)
#   
# }
# 
# head(scenario_lpi_tg_outputs[[1]])
# 
# # Plot trophic groups lpi
# 
# scenario_tg_lpi_plots <- list()
# 
# for (i in seq_along(scenario_lpi_tg_outputs)) {
#   
#   scenario_tg_lpi_plots[[i]] <- ggplot(data = scenario_lpi_tg_outputs[[i]],
#                                        aes(x = annual_time_step, 
#                                            y = indicator_score,
#                                            col = replicate)) +
#     geom_line() +
#     facet_wrap( ~ replicate)
#   
#   ggsave(file.path(lpi_plots_folder, paste(today, scenarios[[i]], 
#                                            "LPI_troph_groups_reps_averaged.png",
#                                            sep = "_")),
#          scenario_tg_lpi_plots[[i]],  device = "png")   
#   
# }
# 
# scenario_tg_lpi_plots[[1]]
# scenario_tg_lpi_plots[[2]]
# scenario_tg_lpi_plots[[3]]
# scenario_tg_lpi_plots[[4]]
# 
# head(scenario_lpi_tg_outputs[[1]])


# MODEL OUTPUTS ----

# * Get harvested groups only ----

scenario_harvested_groups <- list()

for (i in seq_along(scenario_redlist_data_sampled)) {
  
  scenario <- scenarios[[i]]
  
  if (scenario == "000_Baseline") {
    
    # No disturbance so just get whatever groups
    
    harvested <- scenario_redlist_data_sampled[[i]] %>% 
      mutate(scenario = scenarios[[i]]) %>% 
      group_by(annual_time_step) %>% 
      mutate(indicator_score = sum(ave_abundance, na.rm = TRUE),
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
      group_by(annual_time_step) %>% 
      #mutate(ab_scaled = range01(abundance),
      mutate(indicator_score = sum(ave_abundance, na.rm = TRUE),
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
                    scenario)
    
    #TEMP FIX ----
    
    harvested[1,2] <- harvested[2,2]
    
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
      dplyr::select(annual_time_step, ave_abundance, scenario) %>% 
      distinct(.) %>% 
      #mutate(ab_scaled = range01(abundance)) %>% 
      group_by(annual_time_step) %>% 
      #mutate(indicator_score = sum(ab_scaled, na.rm = TRUE),
      mutate(indicator_score = sum(ave_abundance, na.rm = TRUE),
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
      dplyr::select(annual_time_step, ave_abundance, scenario) %>% 
      distinct(.) %>% 
      group_by(annual_time_step) %>% 
      mutate(indicator_score = sum(ave_abundance, na.rm = TRUE),
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

rm(x)
x<- scenario_harvested_groups[[4]]
head(x)
dim(x)

# * Calculate density of harvested groups ----

scenario_density_harvested <- list()

for (i in seq_along(scenario_harvested_groups)) {
  
  scenario_density_harvested[[i]] <-  scenario_harvested_groups[[i]] %>% 
    distinct(.) %>%
    mutate(indicator_score = indicator_score/first(indicator_score),
           density = indicator_score/sum(indicator_score, na.rm = TRUE)) 
}

head(scenario_density_harvested[[2]])

ggplot(scenario_density_harvested[[2]], aes(annual_time_step, indicator_score)) +
  geom_point()

saveRDS(scenario_density_harvested,
        file.path(indicator_outputs_folder,
                  paste(today, "scenario_density_harvested.rds", sep = "_")))


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
    group_by(functional_group_name, annual_time_step) %>% 
    mutate(indicator_score = mean(ave_abundance, na.rm = TRUE),
           indicator = "total abundance",
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
                  functional_group_name,
                  scenario) %>% 
    distinct(.) %>%
    group_by(functional_group_name)  %>% 
    mutate(indicator_score = indicator_score/first(indicator_score),
           density = indicator_score/sum(indicator_score, na.rm = TRUE)) %>% 
    ungroup(.)
  
}

data <- scenario_fg_abundance[[3]]
head(data)

ggplot(data = data) +
  geom_line(aes(x = annual_time_step, 
                y = density,
                col = functional_group_name)) 

# Plot abundance of functional groups

func_group_abundance_plots <- list()

for(i in seq_along(scenario_fg_abundance)) {
  
  library(ggridges)
  
  scenario_name <- scenarios[[i]]
  
  func_group_abundance_plots[[i]] <- ggplot(scenario_fg_abundance[[i]], 
                                            aes(x = annual_time_step, 
                                                y = functional_group_name,
                                                height = density,
                                                group = functional_group_name,
                                                color = functional_group_name,
                                                fill = functional_group_name)) + 
    geom_density_ridges(stat = "identity", scale = 3, alpha = 0.6) +
    theme_classic() +
    theme(panel.background = element_rect(fill = "snow2"),
          legend.position = "none",
          strip.text = element_text(face="bold", size=9)) +
    geom_vline(xintercept = 100, linetype = "dashed") +
    annotate(x=100,y=+Inf,label="Impact start",vjust=2,geom="label",
             size = 2) +
    geom_vline(xintercept = 200, linetype = "dashed") +
    annotate(x=200,y=+Inf,label="Impact end",vjust=2,geom="label",
             size = 2) +
    labs(x = "Annual time step",
         y = "Model output - relative abundance density",
         title = scenario_name) + 
    theme(plot.title = element_text(size = 7)) 
  
  ggsave(file.path(indicator_outputs_folder, paste(today, scenarios[[i]],
                                                   "functional_group_abundance.png",
                                                   sep = "_")),
         func_group_abundance_plots[[i]],  device = "png")
  
}

func_group_abundance_plots[[1]]
func_group_abundance_plots[[2]]
func_group_abundance_plots[[3]]
func_group_abundance_plots[[4]]

# * Plot density of species ----

scenario_species_abundance <- list()

for ( i in seq_along(scenario_redlist_data_sampled)) {
  
  scenario_species_abundance[[i]] <- scenario_redlist_data_sampled[[i]] %>% 
    mutate(scenario = scenarios[[i]]) %>% 
    ungroup(.) %>% 
    group_by(group_id, annual_time_step) %>% 
    mutate(indicator_score = mean(ave_abundance, na.rm = TRUE),
           indicator = "total abundance",
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
                  group_id,
                  functional_group_name,
                  scenario) %>% 
    distinct(.) %>%
    group_by(group_id)  %>% 
    mutate(indicator_score = indicator_score/first(indicator_score),
           density = indicator_score/sum(indicator_score, na.rm = TRUE)) %>% 
    ungroup(.)
  
}

scenario_species_abundance_plots <- list()

for(i in seq_along(scenario_species_abundance)) {
  
  split_data <- split(scenario_species_abundance[[i]], 
                      scenario_species_abundance[[i]]$functional_group_name)
  fg_plots <- list()
  for (j in seq_along(split_data)){
    
    fg <- split_data[[j]]$functional_group_name[1]
    
    fg_plots[[j]] <- ggplot(split_data[[j]], 
                            aes(x = annual_time_step, 
                                y = group_id,
                                height = density,
                                group = group_id,
                                color = group_id,
                                fill = group_id)) + 
      geom_density_ridges(stat = "identity", scale = 3, alpha = 0.6) +
      theme_classic() +
      theme(panel.background = element_rect(fill = "snow2"),
            legend.position = "none",
            strip.text = element_text(face="bold", size=9)) +
      geom_vline(xintercept = 100, linetype = "dashed") +
      annotate(x=100,y=+Inf,label="Impact start",vjust=2,geom="label",
               size = 3) +
      geom_vline(xintercept = 200, linetype = "dashed") +
      annotate(x=200,y=+Inf,label="Impact end",vjust=2,geom="label",
               size = 3) +
      labs(x = "Annual time step",
           y = "Model output - relative abundance density",
           title = paste(scenarios[[i]], fg, sep = " ")) + 
      theme(plot.title = element_text(size = 15))  + 
      theme(axis.text.x = element_text(size = 1))
    
    ggsave(file.path(indicator_outputs_folder, paste(today, scenarios[[i]], fg,
                                                     "species_group_abundance.png",
                                                     sep = "_")),
           fg_plots[[j]],  device = "png", width = 30, height = 25, units= "in")
    
  }
  
  scenario_species_abundance_plots[[i]] <- fg_plots
  
}

saveRDS(scenario_species_abundance_plots,
        file.path(indicator_plots_folder, paste(
          today,"species_abundance_plots.rds", sep = "_")))

bl_spp_density <- plot_grid(scenario_species_abundance_plots[[1]][[1]],
                            scenario_species_abundance_plots[[1]][[2]],
                            scenario_species_abundance_plots[[1]][[3]],
                            scenario_species_abundance_plots[[1]][[4]],
                            scenario_species_abundance_plots[[1]][[5]],
                            scenario_species_abundance_plots[[1]][[6]],
                            nrow = 2, align = "v")

bl_spp_density
ggsave(file.path(indicator_outputs_folder, paste(today, scenarios[[i]],
                                                 "species_group_abundance.png",
                                                 sep = "_")),
       bl_spp_density,  device = "png", width = 40, height = 20, units= "in")

lu_spp_density <- plot_grid(scenario_species_abundance_plots[[2]][[1]],
                            scenario_species_abundance_plots[[2]][[2]],
                            scenario_species_abundance_plots[[2]][[3]],
                            scenario_species_abundance_plots[[2]][[4]],
                            scenario_species_abundance_plots[[2]][[5]],
                            scenario_species_abundance_plots[[2]][[6]],
                            nrow = 2, align = "v")

lu_spp_density
ggsave(file.path(indicator_outputs_folder, paste(today, scenarios[[i]],
                                                 "species_group_abundance.png",
                                                 sep = "_")),
       lu_spp_density,  device = "png", width = 40, height = 20, units= "in")

ca_spp_density <- plot_grid(scenario_species_abundance_plots[[3]][[1]],
                            scenario_species_abundance_plots[[3]][[2]],
                            scenario_species_abundance_plots[[3]][[3]],
                            scenario_species_abundance_plots[[3]][[4]],
                            scenario_species_abundance_plots[[3]][[5]],
                            scenario_species_abundance_plots[[3]][[6]],
                            nrow = 2, align = "v")

ca_spp_density

ggsave(file.path(indicator_outputs_folder, paste(today, scenarios[[i]],
                                                 "species_group_abundance.png",
                                                 sep = "_")),
       ca_spp_density,  device = "png", width = 40, height = 20, units= "in")


he_spp_density <- plot_grid(scenario_species_abundance_plots[[4]][[1]],
                            scenario_species_abundance_plots[[4]][[2]],
                            scenario_species_abundance_plots[[4]][[3]],
                            scenario_species_abundance_plots[[4]][[4]],
                            scenario_species_abundance_plots[[4]][[5]],
                            scenario_species_abundance_plots[[4]][[6]],
                            nrow = 2, align = "v")

ggsave(file.path(indicator_outputs_folder, paste(today, scenarios[[i]],
                                                 "species_group_abundance.png",
                                                 sep = "_")),
       he_spp_density,  device = "png", width = 40, height = 20, units= "in",
       dpi = 400)

# * Get abundance of trophic groups ----

scenario_tg_abundance <- list()

for (i in seq_along(scenario_redlist_data_sampled)) {
  
  heterotrophs <- scenario_redlist_data_sampled[[i]] %>% 
    ungroup(.) %>% 
    mutate(trophic_group = word(functional_group_name, 1),
           scenario = scenarios[[i]]) %>% 
    group_by(trophic_group, annual_time_step) %>% 
    mutate(indicator_score = mean(abundance, na.rm = TRUE),
           indicator = "total abundance",
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
                  scenario, 
                  trophic_group) %>% 
    distinct(.) %>%
    group_by(trophic_group)  %>% 
    mutate(indicator_score = indicator_score/first(indicator_score),
           density = indicator_score/sum(indicator_score, na.rm = TRUE)) %>% 
    ungroup(.)
  
  autotrophs <- scenario_auto_sampled[[i]] %>% 
    ungroup(.) %>% 
    filter(group_id == "autotrophs") %>% 
    mutate(trophic_group = "autotroph",
           scenario = scenarios[[i]]) %>% 
    group_by(trophic_group, annual_time_step) %>% 
    mutate(indicator_score = mean(abundance, na.rm = TRUE),
           indicator = "total abundance",
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
                  scenario, 
                  trophic_group) %>% 
    distinct(.) %>%
    group_by(trophic_group) %>% 
    mutate(indicator_score = indicator_score/first(indicator_score),
           density = indicator_score/sum(indicator_score, na.rm = TRUE)) %>% 
    ungroup(.)
  
  x <- rbind(heterotrophs, autotrophs)
  
  scenario_tg_abundance[[i]] <- x
}

# * Plot abundance of trophic groups ----

trophic_group_abundance_plots <- list()

for(i in seq_along(scenario_tg_abundance)) {
  
  
  if(scenarios[[i]] == "000_Baseline") {
    
    scenario_palette <- c("grey20", "grey50", "grey70", "grey80")
    
  }
  
  if(scenarios[[i]] == "100_Land_use") {
    
    scenario_palette <- c("#440154FF", "#481568FF", "#482677FF", "#453781FF")
    
    scenario_name <- "1 - Land use scenario"
    
  }
  
  if(scenarios[[i]] == "200_Harvesting_carnivores") {
    
    scenario_palette <- c("#39558CFF", "#32648EFF", "#2D718EFF", "#287D8EFF")
    
    scenario_name <- "2 - Carnivore harvesting scenario"
    
  }
  
  if(scenarios[[i]] == "300_Harvesting_herbivores") {
    
    scenario_palette <- c("#1F968BFF", "#20A386FF", "#29AF7FFF", "#3CBC75FF")
    
    scenario_name <- "3 - Herbivore harvesting scenario"
    
  }
  
  
  trophic_group_abundance_plots[[i]] <- ggplot(scenario_tg_abundance[[i]], 
                                               aes(x = annual_time_step, 
                                                   y = trophic_group,
                                                   height = density,
                                                   group = trophic_group,
                                                   color = trophic_group,
                                                   fill = trophic_group)) + 
    geom_density_ridges(stat = "identity", scale = 3, alpha = 0.6) +
    scale_color_manual(values = scenario_palette) +
    scale_fill_manual(values = scenario_palette) +
    theme_classic() +
    theme(panel.background = element_rect(fill = "snow2"),
          legend.position = "none",
          strip.text = element_text(face="bold", size=9)) +
    geom_vline(xintercept = 100, linetype = "dashed") +
    annotate(x=100,y=+Inf,label="Impact start",vjust=2,geom="label",
             size = 2) +
    geom_vline(xintercept = 200, linetype = "dashed") +
    annotate(x=200,y=+Inf,label="Impact end",vjust=2,geom="label",
             size = 2) +
    labs(x = "Annual time step",
         y = "Model output - relative abundance density",
         title = scenario_name) + 
    theme(plot.title = element_text(size = 7)) 
  
  ggsave(file.path(indicator_outputs_folder, paste(today, scenarios[[i]], 
                                                   "trophic_group_abundance.png",
                                                   sep = "_")),
         trophic_group_abundance_plots[[i]],  device = "png")
  
}

trophic_group_abundance_plots[[1]]
trophic_group_abundance_plots[[2]]
trophic_group_abundance_plots[[3]]
trophic_group_abundance_plots[[4]]

saveRDS(trophic_group_abundance_plots,
        file.path(indicator_plots_folder, paste(today, "trophic_group_density.rds",
                                                sep = "_")))

fig <- plot_grid(trophic_group_abundance_plots[[2]], 
                 trophic_group_abundance_plots[[3]],
                 trophic_group_abundance_plots[[4]],
                 align = "H", 
                 ncol = 3, rel_heights = c(1/8, 1/8, 1/8))

fig

ggsave(file.path(indicator_plots_folder, paste(today, "trophic_group_density.png",
                                               sep = "_")),
       fig,  device = "png")


# COMBINE INDICATORS ----

## Standardise outputs

for ( i in seq_along(scenario_large_spp_rli_outputs_5yr)) {
  
  scenario_large_spp_rli_outputs_5yr[[i]] <- scenario_large_spp_rli_outputs_5yr[[i]] %>%
    ungroup(.) %>%
    dplyr::select(annual_time_step, indicator_score, ci_lower, ci_upper,
                  indicator, scenario)
}

head(scenario_large_spp_rli_outputs_5yr[[1]])

for ( i in seq_along(scenario_all_rli_outputs_v2)) {

  scenario_all_rli_outputs_v2[[i]] <- scenario_all_rli_outputs_v2[[i]] %>%
    ungroup(.) %>%
    dplyr::select(annual_time_step, indicator_score, ci_lower, ci_upper,
                  indicator, scenario)
}

head(scenario_all_rli_outputs_v2[[1]])


##### ----
# for ( i in seq_along(scenario_rli_outputs_5yr)) {
#   
#   scenario_rli_outputs_5yr[[i]] <- scenario_rli_outputs_5yr[[i]] %>% 
#     mutate(scenario = scenarios[[i]]) %>%
#     ungroup(.) %>% 
#     dplyr::select(annual_time_step, indicator_score, ci_lower, ci_upper,
#                   indicator, scenario) %>% 
#     mutate(indicator = "RLI 5yr")
# }
# 
for ( i in seq_along(scenario_all_rli_outputs)) {

  scenario_all_rli_outputs[[i]] <-scenario_all_rli_outputs[[i]] %>%
                               mutate(scenario = scenarios[[i]]) %>%
    ungroup(.) %>%
    dplyr::select(annual_time_step, indicator_score, ci_lower, ci_upper,
                  indicator, scenario) %>% 
    mutate(indicator = "RLI annual")
}
# 
# head(scenario_rli_outputs[[2]])
# 
# for ( i in seq_along(scenario_fg_rli_outputs)) {
#   
#   scenario_fg_rli_outputs[[i]] <- scenario_fg_rli_outputs[[i]] %>% 
#     mutate(scenario = scenarios[[i]],
#            indicator = paste(functional_group_name, "RLI", sep = " ")) %>%
#     ungroup(.) %>% 
#     dplyr::select(annual_time_step, indicator_score, ci_lower, ci_upper,
#                   indicator, scenario)
# }
# 
# head(scenario_fg_rli_outputs[[2]])
#####

for ( i in seq_along(scenario_lpi_outputs)) {
  
  scenario_lpi_outputs[[i]] <- scenario_lpi_outputs[[i]] %>% 
    mutate(scenario = scenarios[[i]],
           indicator = "LPI") %>% 
    dplyr::select(annual_time_step, indicator_score, ci_lower, ci_upper,
                  indicator, scenario)
}

head(scenario_lpi_outputs[[2]])

#####
# for ( i in seq_along(scenario_lpi_fg_outputs)) {
#   
#   scenario_lpi_fg_outputs[[i]] <- scenario_lpi_fg_outputs[[i]] %>% 
#     mutate(scenario = scenarios[[i]],
#            indicator = paste(replicate, "LPI", sep = " ")) %>% 
#     dplyr::select(annual_time_step, indicator_score, ci_lower, ci_upper,
#                   indicator, scenario)
# }
# 
# head(scenario_lpi_fg_outputs[[2]])
# 
# for ( i in seq_along(scenario_lpi_tg_outputs)) {
#   
#   scenario_lpi_tg_outputs[[i]] <- scenario_lpi_tg_outputs[[i]] %>% 
#     mutate(scenario = scenarios[[i]],
#            indicator = paste(replicate, "LPI", sep = " ")) %>% 
#     dplyr::select(annual_time_step, indicator_score, ci_lower, ci_upper,
#                   indicator, scenario)
# }
# 
# head(scenario_lpi_tg_outputs[[2]])
# 
# head(scenario_fg_abundance[[2]])
# 
# for ( i in seq_along(scenario_fg_abundance)) {
#   
#   scenario_fg_abundance[[i]] <- scenario_fg_abundance[[i]] %>% 
#         dplyr::select(annual_time_step, indicator_score, ci_lower, ci_upper,
#                   indicator, scenario)
# }
# 
# head(scenario_fg_abundance[[2]])
# 
# head(scenario_harvested_groups[[2]])
#####

for ( i in seq_along(scenario_harvested_groups)) {
  
  scenario_harvested_groups[[i]] <- scenario_harvested_groups[[i]] %>% 
    dplyr::select(annual_time_step, indicator_score, ci_lower, ci_upper,
                  indicator, scenario)
}

head(scenario_harvested_groups[[2]])

names(scenario_rli_outputs_5yr[[1]]) == names(scenario_large_spp_rli_outputs_5yr[[1]])

identical(names(scenario_all_rli_outputs_v2[[1]]), 
          names(scenario_large_spp_rli_outputs_5yr[[1]]))

# Save outputs

# all_indicators_list <- list(scenario_rli_outputs_5yr,
#                             scenario_rli_outputs,
#                             scenario_fg_rli_outputs,
#                             scenario_lpi_outputs,
#                             scenario_lpi_fg_outputs,
#                             scenario_lpi_tg_outputs,
#                             scenario_fg_abundance,
#                             scenario_harvested_groups)
# 
# names(all_indicators_list) <- c("RLI 5yr",
#                                 "RLI",
#                                 "RLI functional groups",
#                                 "LPI",
#                                 "LPI functional groups",
#                                 "LPI trophic groups",
#                                 "abundance functional groups",
#                                 "abundance harvested groups")

all_indicators_list <- list(scenario_all_rli_outputs_v2, #5yrs, all spp
                            scenario_large_spp_rli_outputs_5yr, #5yrs, large spp
                            scenario_all_rli_outputs, #annual
                            scenario_lpi_outputs,
                            scenario_harvested_groups)

names(all_indicators_list) <- c("RLI all",
                                "RLI large spp",
                                "RLI annual",
                                "LPI",
                                "abundance harvested groups")

saveRDS(all_indicators_list,
        file.path(general_indicator_outputs_folder,
                 "all_indicators_output_list.rds"))

all_indicators_all_scenarios <- flatten(all_indicators_list)
all_indicators <- do.call(rbind, all_indicators_all_scenarios)
unique(all_indicators$indicator) # Check we have them all

saveRDS(all_indicators,
        file.path(general_indicator_outputs_folder,
                   "all_indicators_output_dataframe.rds"))

write.csv(all_indicators,
          file.path(general_indicator_outputs_folder,
                    "all_indicators_output_dataframe.csv"))

# Reshuffle the list so the top level of the list is indicator, then scenario

new_split <- split(all_indicators, all_indicators$indicator)

new_indicator_list <- list()
indicator_names <- list()

for (i in seq_along(new_split)) {
  
  indicator_data <- new_split[[i]]
  
  new_indicator_list[[i]] <- split(indicator_data, indicator_data$scenario)
  
  indicator_names[i] <- indicator_data$indicator[1]

}

names(new_indicator_list) <- indicator_names

saveRDS(new_indicator_list,
        file.path(general_indicator_outputs_folder,
                  "all_indicators_output_list_v2.rds"))
