
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
  
  if (class_name == "bird") {
    
    species_in_domain <- ranges_domain %>%
      dplyr::select(-Shape, PRESENCE) %>%
      rename(Scientific.Name = SCINAME,
             id_no = SISID) %>%
      mutate(source = "birdlife_international") 
    
  } else {
    
    species_in_domain <- as.data.frame(ranges_domain %>%
                                            dplyr::select(binomial, 
                                                          class, family)) %>%
      mutate(source = "iucn_redlist",
             family = tolower(family),
             Common.Name = NA) %>% 
      rename(Scientific.Name = binomial,
             Family = family,
             Taxonomic.Group = class) %>% 
      select(Scientific.Name, Common.Name, Family, Taxonomic.Group, source)
    
  }
  
  return(species_in_domain)
  
}

# Make domain polygon ----

## Get africa for reference

coords_sf <- st_as_sf(x = coords_df, #update here
                    coords = c("Longitude", "Latitude"))

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

  amphibian_rangemap_dir <- file.path(inputs, "redlist_amphibian_range_maps")
  
  # Read in the range map
  
  amphibian_ranges <- st_read(amphibian_rangemap_dir)
  
  # Remove unnecessary columns as well
  
  amphibians <- get_species_in_domain(amphibian_ranges, 
                                      study_site_sf,
                                      "amphibian")
  
  amphibians$Taxonomic.Group <- tolower(amphibians$Taxonomic.Group)
  
  write.csv(amphibians, file.path(outputs, "serengeti_amphibians_df.csv"))
  
  # * Mammals ----

  mammal_rangemap_dir <- file.path(inputs, "redlist_mammal_range_maps")
  
  # Read in the rangemap
  
  mammal_ranges <- st_read(mammal_rangemap_dir)
  
  mammals <- get_species_in_domain(mammal_ranges, 
                        study_site_sf,
                        "mammals")
  
  rm(mammal_ranges)
  
  mammals$Taxonomic.Group <- tolower(mammals$Taxonomic.Group)
  
  #head(mammals)
  
  write.csv(mammals, file.path(outputs, "serengeti_mammals_df.csv"))
  
# * Birds ----

  class <- "birds"

## Note that bird maps come from birdlife international, not iucn, so
## they are stored in a geodatabase with slightly different geometry types, so
## require a couple of extra steps to process

  bird_rangemap_dir <- file.path(inputs, "birdlife_avian_range_maps","BOTW.gdb")
  
  # Read in the rangemap
  
  bird_ranges <- st_read(bird_rangemap_dir, layer = "All_Species")
  
  # Remove unneccessary columns 
  
  bird_ranges <- bird_ranges %>%
    select(SISID,
           SCINAME, 
           PRESENCE, 
           Shape)
  
  # bird_ranges contains both multipolygon and multisurface
  # geometry types, which means other st functions won't work.
  
  # Get only multisurface geoms
  #' TODO: Figure out how to recast these to multipolygons
  
  # bird_ranges_ms <- bird_ranges %>%
  #   filter(st_geometry_type(Shape) == "MULTISURFACE")
  
  # Get only multipolygon geoms
  
  bird_ranges <- bird_ranges %>%
                 filter(st_geometry_type(Shape) == "MULTIPOLYGON")
  
  
  birds <- get_species_in_domain(bird_ranges, 
                        study_site_sf, "bird")
  
  rm(bird_ranges)
  
  st_write(birds, file.path(outputs, "serengeti_birds_spatial.shp"))
  
  birds_df <- birds %>% 
              dplyr::select(-Shape, -PRESENCE) %>% 
              mutate(class = "AVES") %>% 
              dplyr::select(id_no, binomial, class, source) %>% 
              st_drop_geometry()
  
  
  head(birds_df)
  
  write.csv(birds_df, file.path(outputs, "serengeti_birds_df.csv"))
 

# * Reptiles ----

  reptile_rangemap_dir <- file.path(inputs, "redlist_reptile_range_maps")
  
  # Read in the rangemap
  
  reptile_ranges <- st_read(reptile_rangemap_dir)
  

  reptiles <- get_species_in_domain(reptile_ranges, 
                                    study_site_sf,
                                    "reptile")
  
  rm(reptile_ranges)
  
  # Add in the point data
  
  reptile_points <- read.csv(file.path(inputs, "redlist_reptile_range_maps",
                                       "REPTILES_points.csv"))
  # Select necessary columns
  
  reptile_points <- reptile_points %>%
    dplyr::select(binomial, latitude, longitude, category) 
  
  # Convert to sf object and set crs
  
  reptile_points_sf <- st_as_sf(reptile_points, coords = c('longitude', 
                                                           'latitude'), 
                                crs = st_crs(study_site_sf))
  
  # Get ecoregions the points fall within
  
  reptile_points <- st_intersection(reptile_points_sf, 
                                           study_site_sf,
                                           "reptile")
  nrow(reptile_points)
  
  
  st_write(reptiles, file.path(outputs, "serengeti_reptiles_spatial.shp"))
  
  reptiles_df <- reptiles %>% 
    dplyr::select(-geometry)
  
  write.csv(reptiles_df, file.path(outputs, "serengeti_reptiles_df.csv"))

# Consolidate species list ----
  
serengeti_iucn_species_list <- rbind(amphibians_df, 
                                     mammals_df,
                                     birds_df,
                                     reptiles_df) 

# Get MoL species list

serengeti_mol_species_list <- read.csv("N:/Quantitative-Ecology/Indicators-Project/Serengeti/Inputs_to_adaptor_code/Species_lists/map_of_life.csv")

serengeti_mol_species_list <- serengeti_mol_species_list %>% 
                              mutate(source = "map of life")

# Combine sources

serengeti_species_list <- rbind(serengeti_iucn_species_list,
                                serengeti_mol_species_list)

## How many species?

length(unique(serengeti_species_list$binomial))

##Save list

write.csv(serengeti_species_list, file.path("N:\Quantitative-Ecology\Indicators-Project\Serengeti\Inputs_to_adaptor_code\Species_lists", 
                                            paste(today, "serengeti_species_list.csv")))