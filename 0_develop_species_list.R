rm(list = ls())
# Set up ----

library(sf)
library(tidyverse)
#install.packages("rnaturalearth")
library(rnaturalearth)
#install.packages("rnaturalearthdata")
library(rnaturalearthdata)
library(rgeos)
library(taxizedb)

# Inputs ----

species_inputs <- "N:/Quantitative-Ecology/Simone/extinction_test/inputs"

outputs <- "N:\\Quantitative-Ecology\\Indicators-Project\\Serengeti\\Outputs_from_analysis_code\\temp"

get_species_in_domain <- function(range_map, domain, class_name = NA) {
  
   ranges_domain <- range_map[st_intersects(range_map, domain) %>% 
                                                 lengths > 0,]
  
  # Standardise the data so we can add it to the species data easily
  
  if (class_name == "Birds") {
    
    species_in_domain <- ranges_domain %>%
      dplyr::select(-PRESENCE) %>%
      rename(Scientific.Name = SCINAME,
             id_no = SISID) %>%
      mutate(source = "birdlife_international",
             Taxonomic.Group = class_name,
             Family = NA,
             Common.Name = NA) %>% 
      select(Scientific.Name, Common.Name, Family, Taxonomic.Group, source)
    
  } else {
    
    species_in_domain <- as.data.frame(ranges_domain %>%
                                            dplyr::select(binomial, 
                                                          class, family)) %>%
      mutate(source = "iucn_redlist",
             family = str_to_title(family),
             Common.Name = NA,
             Taxonomic.Group = class_name) %>% 
      rename(Scientific.Name = binomial,
             Family = family) %>% 
      select(Scientific.Name, Common.Name, Family, Taxonomic.Group, source)
    
  }
  
  return(species_in_domain)
  
}

# Load data ----

serengeti <- st_read("N:\\Quantitative-Ecology\\Indicators-Project\\Serengeti\\Outputs_from_analysis_code\\Analysis_inputs\\Serengeti_Ecosystem\\v3_serengeti_ecosystem.shp")

plot(st_geometry(serengeti))

## Get species crs info

amphibian_rangemap_dir <- file.path(species_inputs, "redlist_amphibian_range_maps")

# Read in the range map

amphibian_ranges <- st_read(amphibian_rangemap_dir)

## read in the coords

coords_df <- read.csv("N:\\Quantitative-Ecology\\Indicators-Project\\Serengeti\\Inputs_to_adaptor_code\\Madingley_simulation_outputs\\100_Land_use\\101_BuildModel\\landuse_1012019-12-12_8.47.41\\SpecificLocations.csv")

# Make domain polygons ----

coords_sf <- st_as_sf(x = coords_df, #update here
                    coords = c("Longitude", "Latitude"))

## Projected to match serengeti map

lon = c(722208.2, 722208.2 + 110000)
lat = c(9668249.5,9668249.5 + 110000 )

Poly_Coord_df = data.frame(lon, lat)

poly <- Poly_Coord_df %>% 
  st_as_sf(coords = c("lon", "lat"), 
           crs = 21036) %>% 
  st_bbox() %>% 
  st_as_sfc()

st_crs(serengeti) == st_crs(poly)

plot(poly)
st_bbox(poly)
st_write(poly, file.path(outputs, "study_site_utm.shp"))

## Geographic to match rangemaps

lon2 = c(35, 36)
lat2 = c(-3,-2)

Poly_Coord_df2 = data.frame(lon2, lat2)

poly_wgs <- Poly_Coord_df2 %>% 
  st_as_sf(coords = c("lon2", "lat2"), 
           crs = 4326) %>% 
  st_bbox() %>% 
  st_as_sfc()

st_crs(amphibian_ranges) == st_crs(poly_wgs)

st_bbox(poly_wgs)
st_write(poly_wgs, file.path(outputs, "study_site_wgs.shp"))

# Get species ---- 

# * Amphibians ----

  # Remove unnecessary columns as well
  
  amphibians <- get_species_in_domain(amphibian_ranges, 
                                      poly_wgs,
                                      "Amphibians")
  
  amphibians$Taxonomic.Group <- tolower(amphibians$Taxonomic.Group)
  
  write.csv(amphibians, file.path(outputs, "serengeti_amphibians_df.csv"))
  
  # * Mammals ----

  mammal_rangemap_dir <- file.path(species_inputs, "redlist_mammal_range_maps")
  
  # Read in the rangemap
  
  mammal_ranges <- st_read(mammal_rangemap_dir)
  
  mammals <- get_species_in_domain(mammal_ranges, 
                        poly_wgs,
                        "Mammals")
  
  rm(mammal_ranges)
  
  mammals$Taxonomic.Group <- tolower(mammals$Taxonomic.Group)
  
  #head(mammals)
  
  write.csv(mammals, file.path(outputs, "serengeti_mammals_df.csv"))
  
# * Reptiles ----
  
  ## Don't need to worry about processing the points bc none are in the serengeti

  reptile_rangemap_dir <- file.path(species_inputs, "redlist_reptile_range_maps")
  
  # Read in the rangemap
  
  reptile_ranges <- st_read(reptile_rangemap_dir)
  

  reptiles <- get_species_in_domain(reptile_ranges, 
                                    poly_wgs,
                                    "Reptiles")
  
  rm(reptile_ranges)
  
  reptiles_df <- reptiles %>% 
    dplyr::select(-geometry)
  
  write.csv(reptiles_df, file.path(outputs, "serengeti_reptiles_df.csv"))
  
  
  # * Birds ----
  
  ## Note that bird maps come from birdlife international, not iucn, so
  ## they are stored in a geodatabase with slightly different geometry types, so
  ## require a couple of extra steps to process
  
  bird_rangemap_dir <- file.path(species_inputs, "birdlife_avian_range_maps","BOTW.gdb")

## WARNING SLOW CODE - bird ranges are huge, take ages to load
  
# Read in the rangemap
  
bird_ranges <- st_read(bird_rangemap_dir, layer = "All_Species")
  
names(bird_ranges)

# Get only the columns we need
  
bird_ranges <- bird_ranges %>%
               select(SISID,
                     SCINAME, 
                     PRESENCE, 
                     Shape)
            
# bird_ranges contains both multipolygon and multisurface
# geometry types, which alot of sf functions don't like
  
# Get only multipolygon geoms bc don't know how to intersect the multisurface geoms
  
bird_ranges <- bird_ranges %>%
               filter(st_geometry_type(Shape) == "MULTIPOLYGON")
  
## WARNING  - SLOW CODE - takes like an hour, so annoying
  
birds <- get_species_in_domain(bird_ranges, 
                               poly_wgs, "Birds")
  
rm(bird_ranges)
  
birds_df <- birds %>% 
            st_drop_geometry() %>% 
            distinct(.)
  
  
head(birds_df)
  
write.csv(birds_df, file.path(outputs, "serengeti_birds_df.csv"))
  

# Consolidate species list ----
  
serengeti_iucn_species_list <- rbind(amphibians, 
                                     mammals,
                                     birds_df,
                                     reptiles) 

# Get MoL species list

serengeti_mol_species_list <- read.csv("N:/Quantitative-Ecology/Indicators-Project/Serengeti/Inputs_to_adaptor_code/Species_lists/map_of_life.csv")

serengeti_mol_species_list <- serengeti_mol_species_list %>% 
                              mutate(source = "map of life")

unique(serengeti_mol_species_list$Taxonomic.Group)

# Combine sources

serengeti_species_list <- rbind(serengeti_iucn_species_list,
                                serengeti_mol_species_list)

## How many species?

length(unique(serengeti_species_list$Scientific.Name))

serengeti_species_list <- serengeti_species_list %>% 
                          distinct(.)

##Save list

write.csv(serengeti_species_list, 
          file.path("N:/Quantitative-Ecology/Indicators-Project/Serengeti/Inputs_to_adaptor_code/Species_lists", 
                    "serengeti_species_list_mol_iucn.csv"))

# Get functional traits ----

SourceToModels <- "C:/Users/ssteven/Desktop/git_repos"
IndicatorsProject <- "N:/Quantitative-Ecology/Indicators-Project"

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
             "Species_lists", "serengeti_species_list_mol_iucn.csv" ), 
  sep = ",", header = TRUE, stringsAsFactors = FALSE, quote = ""
)

species_list <- species_list$Scientific.Name
species_list <- unique( species_list )

# run the code below for every new species list you add
# or every time you update the list of species
# need to change the example directory folders before you use it

  process_species_list(
    species_list,
    databases,
    "N:\\Quantitative-Ecology\\Indicators-Project\\Serengeti\\Outputs_from_analysis_code\\supporting_info_plots_folder"
  )
  
        