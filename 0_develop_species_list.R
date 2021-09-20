
# Set up ----

library(sf)
library(tidyverse)
#install.packages("rnaturalearth")
library(rnaturalearth)
#install.packages("rnaturalearthdata")
library(rnaturalearthdata)
library(rgeos)

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
             Common.Name = NA) %>% 
      rename(Scientific.Name = binomial,
             Family = family,
             Taxonomic.Group = class_name) %>% 
      select(Scientific.Name, Common.Name, Family, Taxonomic.Group, source)
    
  }
  
  return(species_in_domain)
  
}

# Make domain polygon ----

## read in the coords

coords_df <- read.csv("N:\\Quantitative-Ecology\\Indicators-Project\\Serengeti\\Inputs_to_adaptor_code\\Madingley_simulation_outputs\\100_Land_use\\101_BuildModel\\landuse_1012019-12-12_8.47.41\\SpecificLocations.csv")

coords_sf <- st_as_sf(x = coords_df, #update here
                    coords = c("Longitude", "Latitude"))

## Get species crs info

amphibian_rangemap_dir <- file.path(species_inputs, "redlist_amphibian_range_maps")

# Read in the range map

amphibian_ranges <- st_read(amphibian_rangemap_dir)

# Check if the CRS is projected (will be FALSE or NA)
st_is_longlat(coords_sf)

# Set the CRS and check again (should be TRUE now)
coords_sf <- st_set_crs(coords_sf, st_crs(amphibian_ranges)) # Set to native crs

st_is_longlat(coords_sf)
st_crs(coords_sf)
st_write(coords_sf, file.path(outputs, "coords.shp"))

# Transform to a projected CRS for Tanzania UTM 35S

coords_sf_utm <- st_transform(coords_sf, 21035)
st_crs(coords_sf)

# Conver to an sp object
coords_sp <- as(coords_sf, Class = "Spatial")
summary(coords_sp)


study_site_sp <- gBuffer(coords_sp, width = 110000/2, capStyle = "SQUARE")

study_site_sf_utm <- st_as_sf(study_site_sp)

study_site_sf <- st_transform(study_site_sf_utm, st_crs(amphibian_ranges))

st_write(study_site_sf, file.path(outputs, "study_site.shp"))

world <- ne_countries(scale = "medium", returnclass = "sf")
st_write(world, file.path(outputs,"world.shp"))

countries <- world %>% 
             filter(geounit == "Tanzania"|
                      geounit == "Kenya")

st_write(countries ,file.path(outputs, "kenya_tanzania.shp"))


# Get species ---- 

# * Amphibians ----

  # Remove unnecessary columns as well
  
  amphibians <- get_species_in_domain(amphibian_ranges, 
                                      study_site_sf,
                                      "Amphibians")
  
  amphibians$Taxonomic.Group <- tolower(amphibians$Taxonomic.Group)
  
  write.csv(amphibians, file.path(outputs, "serengeti_amphibians_df.csv"))
  
  # * Mammals ----

  mammal_rangemap_dir <- file.path(species_inputs, "redlist_mammal_range_maps")
  
  # Read in the rangemap
  
  mammal_ranges <- st_read(mammal_rangemap_dir)
  
  mammals <- get_species_in_domain(mammal_ranges, 
                        study_site_sf,
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
                                    study_site_sf,
                                    "Reptiles")
  
  rm(reptile_ranges)
  
   st_write(reptiles, file.path(outputs, "serengeti_reptiles_spatial.shp"))
  
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
                                 study_site_sf, "Birds")
  
rm(bird_ranges)
  
birds_df <- birds %>% 
            st_drop_geometry() %>% 
            distinct(.)
  
  
head(birds_df)
  
write.csv(birds_df, file.path(outputs, "serengeti_birds_df.csv"))
  

# Consolidate species list ----
  
serengeti_iucn_species_list <- rbind(amphibians_df, 
                                     mammals_df,
                                     birds_df,
                                     reptiles_df) 

# Get MoL species list

serengeti_mol_species_list <- read.csv("N:/Quantitative-Ecology/Indicators-Project/Serengeti/Inputs_to_adaptor_code/Species_lists/map_of_life.csv")

serengeti_mol_species_list <- serengeti_mol_species_list %>% 
                              mutate(source = "map of life")

unique(serengeti_mol_species_list$Taxonomic.Group)

# Combine sources

serengeti_species_list <- rbind(serengeti_iucn_species_list,
                                serengeti_mol_species_list)

## How many species?

length(unique(serengeti_species_list$binomial))

##Save list

write.csv(serengeti_species_list, file.path("N:\Quantitative-Ecology\Indicators-Project\Serengeti\Inputs_to_adaptor_code\Species_lists", 
                                            paste(today, "serengeti_species_list.csv")))