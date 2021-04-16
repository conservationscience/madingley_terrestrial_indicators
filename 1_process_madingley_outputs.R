
# https://github.com/Simonestevo/marine_to_terrestrial

# cd "C:\\Users\\ssteven\\Dropbox\\Deakin\\Serengeti-analysis\\marine_to_terrestrial"
# cd "C:\\Users\\ssteven\\Desktop\\Serengeti-analysis\\marine_to_terrestrial"

# TODO: Add an if statement to species list code so it only runs if species list
# hasn't been built (currently just commented out)
# TODO: Add some details about the data structure up here, then make sure everything
# is named consistently according to it's hierarchy and step in the process

# Libraries ----

library(tidyverse)
#library(raster)
library(data.table)
#install_github("conservationscience/functionaltraits")
library(functionaltraits)
library(ncdf4)
# devtools::install_github("ropensci/taxizedb")
library(taxizedb)


# Define mode (development == TRUE will only perform operations on a small subset
# of folders not all outputs)

development_mode <- TRUE

# Define whether you need to generate a new species list

make_species_list <- FALSE

# Specify universal input arguments

if (development_mode == FALSE) {
  
  impact_start <- 1100
  impact_end <- 1200
  burnin <- 1000*12 # in months
  n <- 12
  
} else {
  
  impact_start <- 10
  impact_end <- 20
  burnin <- 1*12 # in months
  n <- 1 
  
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

# Source functions ----

# note- would like to turn species_to_model repository into a package so you can eventually just run install_github( "conservationscience/species_to_models" )
source( file.path( SourceToModels, "species_to_models", "madingley_get_groups.R" ) )
source( file.path( SourceToModels, "species_to_models", "madingley_get_species_and_groups_key.R" ) )
source( file.path( SourceToModels, "species_to_models", "madingley_process_trait_data.R" ) )
source( file.path( SourceToModels, "species_to_models", "madingley_get_biomass_of_groups.R" ) )
source( file.path( SourceToModels, "species_to_models", "madingley_get_abundance_of_groups.R" ) )
source( file.path( SourceToModels, "species_to_models", "madingley_get_age_structure_data.R" ) )
source( file.path( SourceToModels, "species_to_models", "madingley_get_autotroph_biomass.R" ) )


source( file.path( SourceToModels, "model_outputs_to_indicator_inputs", 
                   "1_process_outputs", "process_species_list.R") )
source( file.path( SourceToModels, "model_outputs_to_indicator_inputs", 
                   "1_process_outputs", "process_buildmodel_folder.R") )
source( file.path( SourceToModels, "model_outputs_to_indicator_inputs", 
                   "1_process_outputs", "process_output_functions.R") )

source( file.path( SourceToModels, "model_outputs_to_indicator_inputs", 
                   "1_process_outputs", "plot_functional_groups.R") )

source( file.path( SourceToModels, "model_outputs_to_indicator_inputs", 
                   "2_prepare_inputs", "prepare_proportion_total_biomass_inputs.R") )

source( file.path( SourceToModels, "biodiversity_indicators", 
                   "calculate_proportion_biomass.R") )

get_numbers <- function(big, small) {
  
  new <- str_remove(big, paste(small,"/", sep = ""))
  new <- str_remove(new, "_BuildModel")

}

# PART 1 - BUILD SPECIES LIST ----

databases <- functionaltraits::Databases$new( file.path( IndicatorsProject, 
                                                         "functionaltraits_data" ) )

if( !databases$ready() ) {
  print( "Downloading databases to ")
  databases$initialise()
} else {
  print( "Databases already downloaded, in ")
}
print( databases$dir )

species_list <- read.table( 
  file.path( IndicatorsProject, "Serengeti", "Inputs_to_adaptor_code", 
             "Species_lists", "map_of_life.csv" ), 
  sep = ",", header = TRUE, stringsAsFactors = FALSE, quote = ""
)

species_list <- species_list$Scientific.Name
species_list <- unique( species_list )

# run the code below for every new species list you add
# or every time you update the list of species
# need to change the example directory folders before you use it

if(make_species_list == TRUE) {

process_species_list(
   species_list,
   databases,
   file.path( IndicatorsProject, "Serengeti", "Outputs_from_adaptor_code",
              "map_of_life" )
)

} else {
  
  print("Using existing species list")
  
}

# PART 2 - PROCESS MODEL OUTPUTS ----

# Get model output file locations ----

indicators_project <- IndicatorsProject # File path for entire project directory

location <- 'Serengeti' # Modelled location you want to process

# scenarios <- c("Baseline", "Climate_change", "Harvesting_carnivores", 
#                "Harvesting_herbivores", "Land_use", "Test_runs")

model_outputs_path <- file.path(IndicatorsProject, location,
                      "Inputs_to_adaptor_code/Madingley_simulation_outputs")

scenario_paths <- list.dirs(model_outputs_path, recursive = FALSE)

scenarios <- lapply(scenario_paths, basename)


# Get lists of the directories at various levels needed

## Scenario directories


if (development_mode == FALSE) {

#scenario_paths <- scenario_paths[!str_detect(scenario_paths, "Climate_change")] # Remove empty/unneeded scenarios
scenario_paths <- scenario_paths[!str_detect(scenario_paths, "999_Test_runs")] # Remove empty/unneeded scenarios

}

if (development_mode == TRUE) {
  
  #scenario_paths <- scenario_paths[!str_detect(scenario_paths, "Climate_change")] # Remove empty/unneeded scenarios
  scenario_paths <- scenario_paths[str_detect(scenario_paths, "999_Test_runs")] # Remove empty/unneeded scenarios
  
}

scenarios <- lapply(scenario_paths, basename)

scenarios

## Get a list of the simulation directories within each scenario

simulation_paths <- list()

for (i in seq_along(scenario_paths)) {
  
  simulation_paths[[i]] <- list.dirs(scenario_paths[[i]], recursive = FALSE)
  
}

simulation_numbers <- list()

for (i in seq_along(simulation_paths)) {
  
  x <- simulation_paths[[i]]
  
  y <- sapply(x, get_numbers, scenario_paths[i])
  
  y <- unname(y)
  
  simulation_numbers[[i]] <- y
  
}

simulation_folder_names <- list()

for (i in seq_along(simulation_paths)) {
  
  x <- simulation_paths[[i]]
  
  y <- sapply(x, str_remove, paste(scenario_paths[i],"/", sep = ""))
  
  y <- unname(y)
  
  simulation_folder_names[[i]] <- y
  
}

## Get a list of the directories for each simulation results containing model output

simulation_results_paths <- list()

for (i in seq_along(simulation_paths)) {
  
  sub_list <- simulation_paths[i]
  
  for (j in seq_along(sub_list)) {
    
    sub_folders <- list.dirs(sub_list[[j]], recursive = FALSE)
    
    rep_folder <- sub_folders[!str_detect(sub_folders, "input")]
    
    rep_folder <- paste(rep_folder, "/", sep = "")
    
  }
  
  simulation_results_paths[[i]] <- rep_folder
  
  rm(sub_list)
  
}

# Check all our folders match

simulation_paths
scenarios
simulation_results_paths
simulation_numbers
simulation_folder_names

# Process outputs ----

# Loop through each scenario and pull out the folders containing simulation runs

for (i in seq_along(simulation_paths)) {
  
scenario_simulation_paths <- simulation_paths[[i]]

for (j in seq_along(scenario_simulation_paths)) {

# For each simulation run, convert netcdf files into .rds and .csv files & save 
  
process_buildmodel_folder(
  
  scenario_simulation_paths[j], file.path(IndicatorsProject, location, 
                                          "Outputs_from_adaptor_code/map_of_life", 
                                          scenarios[[i]])
  )
 }
}


# Plot biomass and location ----

for (i in seq_along(simulation_results_paths)) {
  
  for (j in seq_along(simulation_results_paths[[i]])) {
  
save_plot(
  
  simulation_results_paths[[i]][j],
  file.path(IndicatorsProject, location, "Outputs_from_adaptor_code/map_of_life",
            scenarios[i], simulation_folder_names[[i]][j]),
            simulation_numbers[[i]][j],"Africa"
  
 )
  }
}

# Remove juveniles, get generation lengths and save plot of abundance in each 
# in each bodymass bin for each functional group

##### WARNING - SLOW CODE, USES A LOT OF MEMORY, RUN ON HPC #########
#' TODO: Save into different folders than the other processed outputs
#' 



for (i in seq_along(simulation_paths)) {
  
  for (j in seq_along(simulation_paths[[i]])) {
    
     get_generation_lengths(simulation_paths[[i]][j], 
                            file.path(IndicatorsProject, location, 
                                      "Outputs_from_adaptor_code/map_of_life",
                                      scenarios[i], 
                                      simulation_folder_names[[i]][j]))

   }
}



# get_age_structure_data("N:/Quantitative-Ecology/Indicators-Project/Serengeti/Inputs_to_adaptor_code/Madingley_simulation_outputs/Baseline/001_BuildModel", 
#                        "N:/Quantitative-Ecology/Indicators-Project/Serengeti/Outputs_from_adaptor_code/map_of_life/Madingley_simulation_outputs/Baseline/001_BuildModel",
#                        scenarios[1],
#                        simulation_folder_names[[1]][[1]],
#                        burnin)

# PART 3 - VISUALISE MODEL OUTPUTS ----

#' TODO: Will need to add somewhere code to extract adults

# indicators_project <- IndicatorsProject # File path for entire project directory
# location <- 'Serengeti' # Modelled location you want to process
# 
# scenarios <- c("Baseline", "Climate_change", "Harvesting_carnivores", 
#                "Harvesting_herbivores", "Land_use", "Test_runs")

# Get processed output locations 

processed_outputs_path <- file.path(IndicatorsProject, location,
                                "Outputs_from_adaptor_code/map_of_life")

# Get lists of the directories at various levels needed

## Scenario directories

processed_scenario_paths <- list.dirs(processed_outputs_path, recursive = FALSE)

if (development_mode == TRUE) {
  
  all_processed_scenario_paths <- processed_scenario_paths
  processed_scenario_paths <- all_processed_scenario_paths[str_detect(all_processed_scenario_paths, "999_Test_runs")]

}


if (development_mode == FALSE) {
  
  # processed_scenario_paths <- processed_scenario_paths[!str_detect(processed_scenario_paths, "Climate_change")] # Remove empty/unneeded scenarios
  processed_scenario_paths <- processed_scenario_paths[!str_detect(processed_scenario_paths, "999_Test_runs")] # Remove empty/unneeded scenarios
  
}

## Get a list of the simulation directories within each scenario

processed_simulation_paths <- list()

for (i in seq_along(processed_scenario_paths)) {
  
  processed_simulation_paths[[i]] <- list.dirs(processed_scenario_paths[[i]], recursive = FALSE)
  
}

# Plot simulation functional_groups ----

for (i in seq_along(processed_simulation_paths)) {
  
  processed_simulation_paths_single_scenario <- processed_simulation_paths[[i]]
  simulation_numbers_single_scenario <- simulation_numbers[[i]]
  
  for (j in seq_along(processed_simulation_paths_single_scenario)) {
  
  plot_simulation_functional_groups(processed_simulation_paths_single_scenario[[j]], 
                         simulation_numbers_single_scenario[[j]],  12)
  
  }
}

# Plot scenario functional groups ----

n <- 12*9

scenario_plots <- list()

for (i in seq_along(processed_scenario_paths)) {

scenario_plots[[i]] <- plot_single_scenario_functional_groups(processed_scenario_paths[[i]], 
                       file.path(IndicatorsProject, location, 
                                 "\\Outputs_from_indicator_code\\Indicator_plots\\mean_biomass"),
                                scenarios[[i]], 
                                burnin, n) 
}

# Plot all scenarios functional groups ----


all_scenarios_plot <- plot_all_scenarios_functional_groups(file.path(
                      IndicatorsProject, location, 
                      "\\Outputs_from_indicator_code\\Indicator_plots\\mean_biomass")) 


# PART 4 - PREPARE INDICATOR INPUTS ----

# Proportion of total biomass ----

# Get the file paths for the biomass replicates

replicate_paths <- list()

biomass_replicate_paths <- list()

for (i in seq_along(processed_scenario_paths)) {
  
  scenario_path <- processed_scenario_paths[i] # Get the file path for one scenario
  
  scenario_simulation_paths <- list.dirs(scenario_path, recursive = FALSE) 
  
  # Get the file paths for each simulation within the scenario
  
  replicate_paths <- list()

  for (k in seq_along(scenario_simulation_paths)) { # For each simulation
    
    single_simulation_path <- scenario_simulation_paths[k]
    
    single_simulation_replicates <- list.files(single_simulation_path, pattern = "biomass.rds") # list the biomass files for that simulation and store
    
    replicate_paths[[k]] <- file.path(single_simulation_path, single_simulation_replicates)
  
  }
  
  biomass_replicate_paths[[i]] <- replicate_paths # Get the replicates for each simulation and store in correct scenario
  
}


# Prepare inputs 

## Proportion of total biomass inputs
#' TODO: Check if the outputs are in months or years
#' 

indicator <- "proportion_total_biomass"
interval <- 12
func <- mean

proportion_biomass_inputs <- list()
proportion_biomass_simulations <- list()


for (i in seq_along(biomass_replicate_paths)) { 

  biomass_scenarios <- biomass_replicate_paths[[i]] # Get the input files for a single scenario

  proportion_biomass_simulations <- list()
  
    for (j in seq_along(biomass_scenarios)) { 

    biomass_simulations <- biomass_scenarios[[j]] # Get the input files within a single simulation of the single scenario
   
     proportion_biomass_replicates <- list()
    
     for (k in seq_along(biomass_simulations)) { 

      replicate_numbers <- 0:(length(biomass_simulations) - 1)

      proportion_biomass_replicates[[k]] <- prepare_proportion_total_biomass_inputs( # Save the single replicate output

        biomass_simulations[k], # Get the single replicate input file within a single simulation within a single scenario

        file.path(IndicatorsProject, location,
                  "Outputs_from_indicator_code/Indicator_inputs",
                  indicator, scenarios[i]),
        simulation_numbers[[i]][j], replicate_numbers[k],
        burnin, interval, func )
      
    }

    proportion_biomass_simulations[[j]] <- proportion_biomass_replicates # Add all the replicates into one simulation folder

  }

  proportion_biomass_inputs[[i]] <- proportion_biomass_simulations # Add all the simulation folders to their scenario folder

}

# PART 5 - CALCULATE INDICATOR VALUES FOR EACH REPLICATE ----

# Proportion of total biomass ----

#' TODO: Should this output just be one dataframe at the end now? Is there any
#' benefit to maintaining the file structure at this point?

proportion_biomass_outputs <- list()
proportion_biomass_replicate_outputs <- list()
proportion_biomass_simulation_outputs <- list()


for (i in seq_along(proportion_biomass_inputs)) { # Scenario
  
  scenario_proportion_biomass_input <- proportion_biomass_inputs[[i]]
  
  proportion_biomass_simulation_outputs <- list()
  
  for (j in seq_along(scenario_proportion_biomass_input)) { # Simulation folder
    
    simulation_proportion_biomass_input <- scenario_proportion_biomass_input[[j]]
    
    replicate_numbers <- 0:(length(simulation_proportion_biomass_input) - 1)
    
    proportion_biomass_replicate_outputs <- list()
    
    for (k in seq_along(simulation_proportion_biomass_input)) { # Replicates within simulation folder
      
      replicate_proportion_biomass_input <- simulation_proportion_biomass_input[[k]]
      
      proportion_biomass_replicate_outputs[[k]] <- calculate_proportion_biomass(
        
      replicate_proportion_biomass_input, 
      
      file.path(IndicatorsProject, location, 
                "Outputs_from_indicator_code/Indicator_outputs", 
                indicator, scenarios[i]),
      simulation_numbers[[i]][[j]], replicate_numbers[k],  
      10 )
      
      
    }
    
    proportion_biomass_simulation_outputs[[j]] <- proportion_biomass_replicate_outputs
    
  }
  
  proportion_biomass_outputs[[i]] <- proportion_biomass_simulation_outputs
  
}

# PART 6 - AGGREGATE REPLICATE INDICATOR VALUES ----

# Proportion of total biomass ----

# Flatten the list structure so it contains one list for each scenario, containing
# all replicate indicator results for that scenario (removing the 'Simulation' 
# level of organisation)

proportion_of_biomass_scenario_replicates <- list()

for (i in seq_along(proportion_biomass_outputs)) { # for each scenario

for (j in seq_along(proportion_biomass_inputs[[i]])) { # for each simulation

pbm_simulations <- proportion_biomass_outputs[[i]] # get the replicates from each simulation

pbm_scenarios <- flatten(pbm_simulations) # pool replicates so they're no longer separated by simulation number

  }

proportion_of_biomass_scenario_replicates[[i]] <- pbm_scenarios

}

# Loop through indicator replicates for each scenario, and calculate the mean 
# indicator values over time for that scenario

scenario_proportion_of_biomass_final <- list()

for (i in seq_along(proportion_of_biomass_scenario_replicates)) {
  
  data <- do.call(rbind, proportion_of_biomass_scenario_replicates[[i]])
  
  summary_data <- data %>%
                  group_by(year) %>%
                  dplyr::summarise(mean_relative_proportion_of_biomass = 
                                   mean(relative_proportion_of_biomass),
                                   sd_pbm = sd(relative_proportion_of_biomass),
                                  n = n(),
                                  se_pbm = sd_pbm/sqrt(n),
                                  lower_bound = 
                                  mean_relative_proportion_of_biomass - 1.96 * se_pbm,
                                  upper_bound = 
                                  mean_relative_proportion_of_biomass + 1.96 * se_pbm) %>%
                  mutate(scenario = scenarios[[i]])
 
  
  scenario_proportion_of_biomass_final[[i]] <- summary_data
  
  rm(data, summary_data)
  
}

all_proportion_of_biomass_final <- do.call(rbind, scenario_proportion_of_biomass_final)

# PART 7 - PLOT AGGREGATED INDICATOR VALUES ----

# Proportion of total biomass ----

#' TODO: Save plots
#' TODO: Check the correct way to calculate CI (this is just a guess)
#' TODO: Check correct way to determine pre exploitation period (mean over time period
#' or just one time point)
#' 

# Individual scenario plots

pre_exploitation_period <- 10

scenario_proportion_of_biomass_plots <- list()

for (i in seq_along(scenario_proportion_of_biomass_final)) {
  
   scenario_proportion_of_biomass_plots[[i]] <- plot_total_proportion_biomass(
    scenario_proportion_of_biomass_final[[i]],
    file.path(IndicatorsProject, location, 
              "Outputs_from_indicator_code/Indicator_plots", 
              indicator), scenarios[i],
    pre_exploitation_period, 100, 200
  )
}

# All scenarios in one plot

dev.off()

plot_all_scenarios_total_proportion_biomass(all_proportion_of_biomass_final,
                                      file.path(IndicatorsProject, location, 
                                                "Outputs_from_indicator_code/Indicator_plots", 
                                                indicator), scenarios,
                                      pre_exploitation_period, 100, 200)
