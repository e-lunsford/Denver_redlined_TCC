#-------------------------------------------------------------------------------
# title: Obtain ACS, tree canopy, and redlining data
# author: E Lunsford  
# date: 2025-08-15
#  
# This code is to obtain ACS 5 year estimates for 2017-2021 data from ACS, 
# NLCD 2021 USFS Tree Canopy Cover (CONUS) from mrlc, and redlining data from
# mapping inequality.
#
#
# Last Run: 08/15/2025 and code was in working order 
# using R 4.5.1 and RStudio 2025.05.1+513 
#
#-------------------------------------------------------------------------------


#################################################################################
# Load Libraries.
#################################################################################
library(broom)
library(exactextractr)
library(FedData)
library(ggmap)
library(ggplot2)
library(ggspatial)
library(gstat)
library(ggthemes)
library(keyring)
library(knitr)
library(purrr)
library(raster)
library(sf)
library(stringr)
library(terra)
library(tidycensus)
library(tidyverse)
library(tigris)
library(units)

#################################################################################
#
# Obtain ACS Data
#
# Methods
# https://www2.census.gov/programs-surveys/acs/tech_docs/statistical_testing/2023_Instructions_for_Stat_Testing_ACS.pdf
# 
# SE for each population age group is derived from MOE with Z = 1.645.
# Standard Error = Margin of Error / 1.645
# MOE(A+B) = SQRT(MOE_A + MOE_B)
#
#################################################################################

# Load keyring to obtain saved ACS key
keyring::key_list()


# Load 2021 ACS5 variable list
acs5_var_list <- tidycensus::load_variables(year = 2021,
                                            dataset = "acs5",
                                            cache = T)

# Create list of wanted codes
acs_vars <- c(total_pop = "B02001_001",
              nonhisp_black_pop = "B03002_004",
              nonhisp_white_pop = "B03002_003",
              hisp_pop = "B03002_012",
              income_below = "B17010_002",
              income_at = "B17010_022")

# Get Denver ACS
den_cbg_acs5_2021 <- tidycensus::get_acs(geography = "cbg", #census block group level
                                         variables = acs_vars, # variable list
                                         geometry = TRUE,
                                         key = keyring::key_get(service = "get_acs"),
                                         state = "Colorado",
                                         county = "Denver",
                                         year = 2021, #2017-2021
                                         survey = "acs5") #default of get acs is 5

# Save ACS data to outputs folder
readr::write_csv(file = "Data/den_acs5_2021_nad83.csv",
                 x = den_cbg_acs5_2021)

#################################################################################
# Pivot wider for summaries
#################################################################################

den_cbg_acs5_2021_wide <- den_cbg_acs5_2021 %>%
  pivot_wider(id_cols = c(GEOID, NAME, geometry),
              names_from = variable,
              values_from = c(estimate, moe))

#################################################################################
# Calculate totals and standard errors
#################################################################################

# Set z value
Z = 1.645

pop_table <- den_cbg_acs5_2021_wide %>%
  mutate(
    # population estimates and SEs
    total_pop = estimate_total_pop,
    total_pop_se = moe_total_pop / Z,
    
    # NHW
    nhw_pop = estimate_nonhisp_white_pop,
    nhw_pop_se = moe_nonhisp_white_pop / Z,
    nhw_per = (nhw_pop / total_pop) * 100,
    
    # NHB
    nhb_pop = estimate_nonhisp_black_pop,
    nhb_pop_se = moe_nonhisp_black_pop / Z,
    nhb_per = (nhb_pop / total_pop) * 100,
    
    # Hispanic/Latino
    hisp_pop = estimate_hisp_pop,
    hisp_pop_se = moe_hisp_pop / Z,
    hisp_per = (hisp_pop/total_pop) * 100,
    
    # Poverty (RSME for SE)
    pov_pop = estimate_income_below + estimate_income_at,
    pov_pop_se = sqrt((moe_income_below / Z)^2 + (moe_income_at / Z)^2),
    pov_per = (pov_pop/total_pop) * 100) 

#################################################################################
# Calculate CIs for population data
##################################################################################

pop_table2 <- pop_table %>%
  mutate(
    # Total pop
    total_pop_lower = total_pop - (1.96 * total_pop_se),
    total_pop_upper = total_pop + (1.96 * total_pop_se),
    
    # NHW
    nhw_pop_lower = nhw_pop - (1.96 * nhw_pop_se),
    nhw_pop_upper = nhw_pop + (1.96 * nhw_pop_se),
    
    # NHB
    nhb_pop_lower = nhb_pop - (1.96 * nhb_pop_se),
    nhw_pop_upper = nhb_pop + (1.96 * nhb_pop_se),
    
    # Hisp
    hisp_pop_lower = hisp_pop - (1.96 * hisp_pop_se),
    hisp_pop_upper = hisp_pop + (1.96 * hisp_pop_se),
    
    # Pov
    pov_pop_lower = pov_pop - (1.96 * pov_pop_se),
    pov_pop_upper = pov_pop + (1.96 * pov_pop_se)
  )

pop_table3 <- pop_table2 %>%
  mutate(across(
    .cols = ends_with("_lower"),
    .fns = ~ ifelse(.x < 0, 0, .x)
  ))

# Save data for future use
save(file = "Data/population_tables.RData",
     x = pop_table3)

#################################################################################
#
# Obtain tree Data
#
# Use FedData package to directly import tree canopy for area of interest.
# get_nlcd function returns a SpatRaster NCLD data cropped to template.
# 
#################################################################################

# Load area of interest from "01" code.
load(file = "Data/den_cbg_nad83.RData")
load(file = "Data/den_cbg_wgs84.RData")

# Change CRS to match nlcd CRS (NAD 83 / Conus Albers (EPSG: 5070))
den_5070 <- sf::st_transform(x = den_cbg_nad83,
                             crs = 5070)

# Change sf object to a Spat Vector
den_5070_sv <- terra::vect(x = den_5070)

# Get canopy data
nlcd_tree_2021 <- FedData::get_nlcd(template = den_5070_sv, # Polygon template
                                    label = "tree", # Name
                                    year = 2021,
                                    dataset = "canopy")

# Basic plot to view data
terra::plot(nlcd_tree_2021)

# Export raster data
writeRaster(x = nlcd_tree_2021,
            filename = "Data/nlcd_tree_2021.tif",
            filetype = "GTiff",
            overwrite = TRUE)

#################################################################################
# Explore nlcd characteristics
#
# For reference:
# Geographic Data Science with R (Chapter 6.1).
# https://bookdown.org/mcwimberly/gdswr-book/raster-geospatial-data---continuous.html
#
#################################################################################

# Load in FedData data
nlcd_tree_2021 <- rast("Data/nlcd_tree_2021.tif")

class(nlcd_tree_2021)

# Number of grid cells
ncell(nlcd_tree_2021)

# Number of grid rows
nrow(nlcd_tree_2021)

# Number of grid columns
ncol(nlcd_tree_2021)

# Number of raster layers
nlyr(nlcd_tree_2021)

# Grid cell size
res(nlcd_tree_2021)

# Abbreviated summary of CRS
crs(nlcd_tree_2021, describe = TRUE)


################################################################################
# Extract tree canopy statistics by census blocks
################################################################################

# Get class of both objects
class(nlcd_tree_2021)
class(den_5070_sv)

# Create a SpatExtent
denver_extent <- terra::ext(den_5070_sv)

# Crop tree data to Denver
cropped_raster <- crop(nlcd_tree_2021, denver_extent)


# Extract nlcd raster values to CBG using exactextractr package
sum_stat_tree_values <- as.data.frame(
  exact_extract(x = cropped_raster, 
                y = den_5070, 
                fun = c("mean", "sum",
                        "count", "median",
                        "min", "max",
                        "variance", "stdev")))

# Add "nlcd" to column names for future merging
colnames(sum_stat_tree_values) <- paste0("nlcd_", colnames(sum_stat_tree_values))

# Assign GEOIDs from CBG SF
sum_stat_tree_values$GEOID <- den_5070$GEOID

# Merge extracted NLCD to CBG SF
den_5070_tree <- den_5070 %>%
  left_join(sum_stat_tree_values)

# Save Tree CBG SF data
save(file = "Data/den_5070_tree.RData",
     x = den_5070_tree)

# Transform to NAD83
den_nad83_tree <- st_transform(x = den_5070_tree,
                               crs = 4269)

# Save data
save(file = "Data/den_nad83_tree.RData",
     x = den_nad83_tree)

# Transform to WGS84
den_wgs84_tree <- den_nad83_tree %>%
  st_transform(crs = 4326)

# Save data
save(file = "Data/den_wgs84_tree.RData",
     x = den_wgs84_tree)


#################################################################################
#
# Load redlining data
#
# Download redlining (HOLC) shapefile for City of Denver.
# Nelson, Robert K., LaDale Winling, et al. "Mapping Inequality: Redlining in New Deal America." 
# Edited by Robert K. Nelson and Edward L. Ayers. American Panorama: An Atlas of United States History, 2023.
# https://dsl.richmond.edu/panorama/redlining.  
#
# Note that the redlining data is WGS84
#################################################################################

# Load downloaded data
redline_wgs84 <- sf::st_read("Data/CODenver1938/cartodb-query.shp") %>%
  st_make_valid()

# Save data
save(file = "Data/redline_wgs84.RData",
     x = redline_wgs84)

# Transform CRS to match ACS and tree data
redline_nad83 <- redline_wgs84 %>%
  st_transform(crs = 4269)

# Save data
save(file = "Data/redline_nad83.RData",
     x = redline_nad83)

# Map original HOLC neighborhoods with ggplot and add Denver boundaries
ggplot() +
  geom_sf(data = pop_table3,
          fill = NA,
          colour = "black") +
  geom_sf(data = redline_nad83,
          aes(fill = holc_grade)) +
  scale_fill_manual(values = 
                      c("olivedrab",
                        "lightblue",
                        "lightgoldenrodyellow",
                        "firebrick",
                        "grey")) +
  theme_void()

# Get summary of redlining data
redline_nad83 %>%
  count(holc_grade)





