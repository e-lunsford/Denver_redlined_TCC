#-------------------------------------------------------------------------------
# title: Setting area of interest for analysis
# author: E Lunsford  
# date: 2025-08-14 
#  
# This code is to set the area of interest for further spatial coding and 
# map making.                                                              
#
# Last Run: 08/14/2025 and code was in working order 
# using R 4.5.1 and RStudio 2025.05.1+513 
#-------------------------------------------------------------------------------


#################################################################################
# Load Libraries
#################################################################################
library(broom)
library(colorspace)
library(ggspatial)
library(ggthemes)
library(ggplot2)
library(keyring)
library(knitr)
library(lubridate)
library(purrr)
library(RAQSAPI)
library(raster)
library(sf)
library(stringr)
library(tidycensus)
library(tidyverse)
library(tigris)
library(units)

#################################################################################
# Setup background geometry for analysis.
## Get Colorado shp boundaries and set bounding boxes.
#################################################################################

# Colorado Boundary
co_nad83 <- tigris::states(year = 2021, cb = T) %>%
  filter(GEOID == "08")

# Change CRS to WGS84
co_wgs84 <- tigris::states(year = 2021, cb = T) %>%
  filter(GEOID == "08") %>%
  sf::st_transform(crs = 4326)

# Colorado counties
co_counties_nad83 <- tigris::counties(state = "08",
                                      year = 2021,
                                      cb = T)

# Change CRS to WGS84
co_counties_wgs84 <- co_counties_nad83 %>%
  sf::st_transform(crs = 4326)

# Filter Denver county
den_county_wgs <- co_counties_wgs84 %>%
  filter(NAME == "Denver")

# Obtain ACS 2021 census block group boundaries
den_cbg_nad83 <- tigris::block_groups(state = "08",
                                      county = "Denver",
                                      year = 2021,
                                      cb = T)

den_cbg_wgs84 <- tigris::block_groups(state = "08",
                                        county = "Denver",
                                        year = 2021,
                                        cb = T) %>%
  # Set CRS to wgs84
  sf::st_transform(crs = 4326)

# Save as RData for future use
save(file = "Data/den_cbg_nad83.RData",
     x = den_cbg_nad83)

save(file = "Data/den_cbg_wgs84.RData",
     x = den_cbg_wgs84)





