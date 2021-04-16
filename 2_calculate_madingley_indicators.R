
## REPOSITORY: https://github.com/Simonestevo/marine_to_terrestrial

rm(list = ls())

# Directory path

# cd "C:\\Users\\ssteven\\Dropbox\\Deakin\\Serengeti-analysis\\marine_to_terrestrial"
# cd "C:\\Users\\ssteven\\Desktop\\Serengeti-analysis\\marine_to_terrestrial"
# cd "C:\\Users\\ssteven\Desktop\\git_repos\\marine_to_terrestrial

# TODO LIST ----

#' TODO: CHECK WHY MEAN ABUNDANCE OF EVERY GROUP IN YR 2 = 0 IN TEST DATA
#' TODO: Change LPI to geometric mean as per Mcrae 2017
#' TODO: Test sampling interval by trying different times of year
#' TODO: Check gen length equation
#' TODO: Create final output that combines indicator values for each replicate
#' TODO: Check bootstrapping, looks weird (timesteps out of sync??)

# Libraries ----

## Data wrangling

library(tidyverse)
library(rlpi)

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

maintain_0_abundance <- function(vec) { 
  vec_len = length(vec) 
  first_0 = match(0, vec, nomatch = vec_len) 
  if(first_0 < vec_len) replace(vec, (first_0+1):vec_len, 0) 
  else vec 
}

# data <- rli_inputs
# numboots <- 5
# Function to calculate RLI
calculate_red_list_index <- function(data, numboots, ci = FALSE){
  
  # Using equation from Butchart et al (2007) Improvements to the Red List Index
  
  require(tidyverse)
  
  # Remove data without RL status
  
  #data$redlist_assessment_year <- as.numeric(as.character(data$redlist_assessment_year))
  
  data <- data %>%
          filter(!is.na(rl_status_adjusted)) %>%
          group_by(group_id) 
  
  head(data)
  
  # ecoregion <- as.factor(data$ecoregion_id[1])
  
  # Assign category weights
  
  weighted_data <- data %>%
    dplyr::mutate(rl_weight = ifelse(rl_status_adjusted == "LC", 0,
                              ifelse(rl_status_adjusted == "NT", 1,
                              ifelse(rl_status_adjusted == "VU", 2,
                              ifelse(rl_status_adjusted == "EN", 3,
                              ifelse(rl_status_adjusted == "CR", 4,
                              ifelse(rl_status_adjusted == "EX", 5, NA))))))) 
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
  
  grouped_data <- weighted_data %>% group_by(functional_group_name, time_step)
  
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
  
  # Split by timestep - we want CI for each functional group, for each timestep (I think)
  
  weighted_data_timestep_list <- split(weighted_data, weighted_data$time_step)
  
  ## For each functional group (level 1)

  timestep_confidence_intervals <- list()
  
  for (i in seq_along(weighted_data_timestep_list)) {
    
    # Get single time-step then group by functional group
    
      grouped_timestep_data <- weighted_data_timestep_list[[i]] %>%
                             group_by(functional_group_name)
    
      time <- grouped_timestep_data$time_step[1]
    
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
                                                   probs = 0.05),
                               ci_upper = quantile(rep_scores$RLI, 
                                                   probs = 0.95)) %>%
                   mutate(time_step = time) 
      
      timestep_confidence_intervals[[i]] <- ci_scores
    
  }
  
  confidence_intervals <- do.call(rbind, timestep_confidence_intervals)
  
  red_list_scores <- index_scores %>%
                     merge(confidence_intervals, 
                           by = c("functional_group_name",
                                   "time_step")) %>%
                     dplyr::select(functional_group_name, time_step, ci_lower,
                                    RLI, ci_upper, everything()) 
  
  return(red_list_scores)
  
  } else {
    
 red_list_scores <- index_scores 
 return(red_list_scores)
 
    }
}

  

#' Return a line plot of the RLI scores over time faceted by functional group 

#' @param data a data frame (output from calculate_red_list_index) with columns:
#' functional_group, time_step, ci_lower, ci_upper, total.weight, total.count, RLI criteria
#' @param impact_start time_step impact began
#' @param impact_end time_step impact ended
#' @return a line plot of RLI over time for each functional group/taxa or whatever 

plot_red_list_index_by_group <- function(data, impact_start, impact_end, ci = FALSE) {
  
  require(ggplot2)
  require(viridis)
  
  if (ci == TRUE) {
  
  plot <- ggplot(data = data, aes(x = time_step, y = RLI,
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
    
    plot <- ggplot(data = data, aes(x = time_step, y = RLI,
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
#' @return a line plot of RLI over time for each functional group/taxa or whatever 

plot_red_list_index <- function(data, impact_start, impact_end, ci = FALSE) {
  
  require(ggplot2)
  require(viridis)
  
  if (ci == TRUE) {
  
  plot <- ggplot(data = data, aes(x = time_step, y = RLI)) +
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
    
    plot <- ggplot(data = data, aes(x = time_step, y = RLI)) +
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


data <- scenario_abundance_long[[1]][[1]]

calculate_living_planet_index <- function(data, start_time_step, ci = FALSE){
  
  filtered_inputs <- data %>%
    # Remove timesteps if needed (if not make start_time_step = 1)
    filter(time_step >= start_time_step) %>%
    # Group by virtual species (vs)
    group_by(group_id) %>% 
    # Remove vs that have no individuals at any timestep
    filter(!all(abundance == 0)) %>% 
    # Save a copy of original abundance values
    rename(abundance_original = abundance) %>%
    # Identify rows where abundance == 0, and change all subsequent years to 0 too
    mutate(abundance = maintain_0_abundance(abundance_original)) 
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
                mutate(dt = log10(abundance_adjusted)/ #abundance at current timestep
                            lag(log10(abundance_adjusted), 1)) %>% #abundance at previous timestep
                # Set the initial dt as 1
                mutate(dt = ifelse(time_step == start_time_step, 1, dt)) 
 
  
       head(lpi_inputs)
       group <- 11.43
       test_data <- lpi_inputs %>% filter(group_id == group)
       ggplot(test_data, aes(x = time_step, y = abundance_original)) +
         geom_line() 
  
  # Take the mean rate of change across all species, in each timestep
  # (equation 4 in Mcrae et al 2008)
  
  lpi_inputs_aggregated <- lpi_inputs %>% 
                           group_by(time_step) %>%
                           summarise(mean_dt = mean(dt, na.rm = TRUE),
                                     mean_abundance = mean(abundance_adjusted,
                                                           na.rm = TRUE)) %>%
                           ungroup(.)
               
  head(lpi_inputs_aggregated)
  
  # ggplot(lpi_inputs_aggregated, aes(x = time_step, y = mean_dt)) +
  #   geom_line()
  # 
  # data <- lpi_inputs_aggregated[3:nrow(lpi_inputs_aggregated),]
  # ggplot(data, aes(x = time_step, y = mean_abundance)) +
  #   geom_line()
  
  # Converted averaged rates of change into the LPI                 
  index_scores <- lpi_inputs_aggregated %>%
                  # Set the first time step in LPI values to = 1
                  mutate(LPI = ifelse(time_step == start_time_step, 1,
                                1.5)) 
                  # Next set the next step as LPI in the previous timestep
                  # multiplied by ten to the power of mean_dt
                  
  
  head(index_scores)

  ggplot(index_scores, aes(x = time_step, y = LPI)) +
    geom_line()
  
  return(index_scores)
 
}

data <- index_scores
head(data)

plot_living_planet_index <- function(data) {
  
  ggplot(data, aes(x = time_step, y = LPI)) +
    geom_line()  
  
}

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
# of folders not all outputs)

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
  numboots <- 5
  start_time_step <- 10
  gen_timeframe <- 10
  interval <- 2
  
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

# Get file paths ---- 

processed_outputs_path <- file.path(IndicatorsProject, location,
                                    "Outputs_from_adaptor_code/map_of_life")

# Get lists of the directories at various levels needed

# * Scenario & simulation directories ----

processed_scenario_paths <- list.dirs(processed_outputs_path, recursive = FALSE)

if (development_mode == TRUE) {
  
  all_processed_scenario_paths <- processed_scenario_paths
  processed_scenario_paths <- all_processed_scenario_paths[str_detect(all_processed_scenario_paths, "999_Test_runs")]
  
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

groups <- readRDS(file.path(processed_simulation_paths[[1]][2], "groups.rds"))

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

rm(abundance_reps, abundance_temp, replicate_abundance_formatted, rep)

# Check structure is still correct

length(scenario_abundance_formatted) == length(scenario_abundance_raw)

# Remove previous version to make space

rm(scenario_abundance_raw)

# Convert to long format ----

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
  new_names <- str_remove(new_names, "V")
  colnames(abundance_single) <- new_names
  
  # Convert from wide to long
  replicate_abundance_long[[j]] <- abundance_single %>% 
    rownames_to_column(.) %>%
    pivot_longer(all_of(new_names)) %>%
    rename(time_step = name,
           abundance = value,
           group_id = rowname) %>%
    mutate(time_step = as.numeric(time_step),
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

scenario_ab_gl_formatted <- list()
replicate_ab_gl_formatted <- list()

for (i in seq_along(scenario_abundance_long)) {
  
  replicate_abundance <- scenario_abundance_long[[i]]
  replicate_generations <- scenario_generations_raw[[i]]
  
  # For each individual replicate
  
  for (j in seq_along(replicate_abundance)) {

    replicate_ab_gl_formatted[[j]] <- replicate_abundance[[j]] %>%
        merge(replicate_generations[[j]], by = "group_id") %>%
    arrange(time_step, group_id) %>%
    mutate(generation_by_three = generation_length_yrs * 3) %>% # Time over which to measure decline, 3 x gen length OR:
    mutate(timeframe = ifelse(generation_by_three > gen_timeframe, # 10 years 
                              round(generation_by_three), gen_timeframe)) %>%
    dplyr::select(-generation_by_three, -parent_cohort_id, -adult_mass_g) %>%
    distinct(.)
  
  }
  
  scenario_ab_gl_formatted[[i]] <- replicate_ab_gl_formatted

}

head(scenario_ab_gl_formatted[[1]][[1]])
length(unique(scenario_ab_gl_formatted[[1]][[1]]$group_id))

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
    arrange(time_step) %>%
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
                 arrange(group_id, time_step) %>%
     # Replace all non-ex status with ex after first occurrence 
     mutate(extinct = match("EX", rl_status)) %>%
     mutate(rl_status_adjusted = with(., ave(rl_status, 
                                             FUN=maintain_ex_status)))
  
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

# Plot some results to check they're not completely whack

## Get one group to check how their status changes over time relative to how
## their abundance changes

group_id_select <- "13.16.17" # Shows example of 'resurrected' virtual spp
# group_id_select <- "10.40"

data <- rli_inputs %>% dplyr::filter(group_id == group_id_select)

ggplot(data, aes(x = time_step, y = abundance)) +
  geom_line() +
  geom_text(aes(label= rl_status_adjusted,
                col = rl_status_adjusted),hjust=0, vjust=0)

# * Sample data ----

scenario_red_list_inputs_annual <- list()
replicate_red_list_inputs_annual <- list()

for (i in seq_along(scenario_red_list_data)) {
  
  # Get replicate data for a single scenario
  
  replicate_red_list_data <- scenario_red_list_data[[i]]
  
  # For each individual replicate
  
  for (j in seq_along(replicate_red_list_data)) {
    
  replicate_red_list_inputs_annual[[j]] <-  replicate_red_list_data[[j]] %>%
                                     slice(which(row_number() %% interval == 1))
  
  }
  
  scenario_red_list_inputs_annual[[i]] <- replicate_red_list_inputs_annual
  
}


# * Calculate RLI ----

# RLI by individual functional groups

scenario_fg_rli_outputs <- list()
replicate_fg_rli_outputs <- list()

for (i in seq_along(scenario_red_list_inputs_annual)) {
  
  replicate_red_list_inputs <- scenario_red_list_inputs_annual[[i]]
  
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
  
  
x <- scenario_fg_rli_outputs[[1]][[1]]
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
                                 group_by(time_step) %>%
                                 summarise(RLI = mean(RLI),
                                           ci_lower = mean(ci_lower),
                                           ci_upper = mean(ci_upper)) %>%
                                 mutate(replicate = j)
   } else {
     
   replicate_rli_outputs[[j]] <- replicate_rli_fg[[j]] %>%
                                 group_by(time_step) %>%
                                 summarise(RLI = mean(RLI)) %>%
                                 mutate(replicate = j)
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
                                      impact_end)

  ggsave(file.path(rli_plots_folder, paste(today, scenarios[[i]], "replicate", j,
                "RLI_by_functional_group.png",
                sep = "_")),
       replicate_fg_rli_plots[[i]],  device = "png")

  }
  
scenario_fg_rli_plots[[i]] <- replicate_fg_rli_plots

}

scenario_fg_rli_plots[[1]][[1]]

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
                                                   impact_end)


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
                                                  scenario_rli_outputs[[i]]) 
  
  scenario_mean_rli <- scenario_rli_outputs_aggregated[[i]] %>%
                       group_by(time_step) %>%
                       summarise(RLI = mean(RLI)) %>%
                       mutate(replicate = 0) # Replicate 0 will always be the mean
  
  scenario_rli_outputs_aggregated[[i]] <-  rbind(scenario_rli_outputs_aggregated[[i]],
                                                 scenario_mean_rli) %>%
                                           mutate(replicate = as.factor(replicate)) %>%
                                           mutate(level = ifelse(replicate == 0,
                                                                 "Mean RLI", 
                                                                 "Replicate RLI"))
  
}

head(scenario_rli_outputs_aggregated[[1]])
tail(scenario_rli_outputs_aggregated[[1]])

scenario_rli_plots_aggregated <- list()

for (i in seq_along(scenario_rli_outputs_aggregated)){

scenario_rli_plots_aggregated[[i]] <- ggplot(data = scenario_rli_outputs_aggregated[[i]], 
       aes(x = time_step, y = RLI, group = replicate,
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

data <- abundance_long[[1]]

head(data)

ggplot(filtered_inputs, aes(x = time_step, y = abundance,
                 col = group_id)) +
          geom_line()  + 
          geom_text(aes(label= group_id),hjust=0, vjust=0) +
          theme(legend.position = "none")

# * Calculate LPI ----

scenario_lpi_outputs <- list()

for (i in seq_along(abundance_long)){
  
scenario_lpi_outputs[[i]] 
  
                             

 }
  

}

