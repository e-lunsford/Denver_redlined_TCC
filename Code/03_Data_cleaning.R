#-------------------------------------------------------------------------------
# title: Data cleaning
# author: E Lunsford  
# date: 2025-08-18
#  
# This code is to create one dataset for analysis by merging together
# population, tree canopy, and redlining data
#
#
# Last Run: 08/18/2025 and code was in working order 
# using R 4.5.1 and RStudio 2025.05.1+513 
#
#-------------------------------------------------------------------------------

#################################################################################
# Load Libraries.
#################################################################################
library(broom)
library(ggplot2)
library(ggspatial)
library(gstat)
library(ggthemes)
library(knitr)
library(raster)
library(sf)
library(stringr)
library(terra)
library(tidycensus)
library(tidyverse)
library(tigris)

#################################################################################
# Load data
#################################################################################

# ACS data
load(file = "Data/population_tables.RData")

# Tree canopy data
load(file = "Data/den_nad83_tree.RData")

# Redlining data
load(file = "Data/redline_nad83.RData")


#################################################################################
# Combine ACS & TCC data
#################################################################################

# Assign CRS to modern formatting
acs_sf <- sf::st_set_crs(pop_table3, st_crs(pop_table3))
tree_sf <- sf::st_set_crs(den_nad83_tree, st_crs(den_nad83_tree))

# Make sure ACS and tree CRS match
if (st_crs(acs_sf) != st_crs(tree_sf)) {
  tree_sf <- st_transform(tree_sf, st_crs(acs_sf))
}

# Merge tree data into acs data
merged_sf <- acs_sf %>%
  left_join(st_drop_geometry(tree_sf), by = "GEOID")

#################################################################################
# Prepare redlining data for merge
#################################################################################

# Check CRS of new merged sf and redlining data
redlining_sf <- sf::st_set_crs(redline_nad83, st_crs(redline_nad83))

if (st_crs(merged_sf) != st_crs(redlining_sf)) {
  redlining_sf <- st_transform(redlining_sf, st_crs(merged_sf))
}

# Intersect census block polygons with holc data
holc_col <- "holc_grade"

intersect_sf <- st_intersection(
  merged_sf %>% select(GEOID, geometry),
  redlining_sf %>% select(all_of(holc_col))
)

# Calculate intersection area
intersect_sf <- intersect_sf %>%
  mutate(in_area = st_area(geometry))

# Calculate total area per block
block_area <- merged_sf %>%
  mutate(total_area = st_area(geometry)) %>%
  st_drop_geometry() %>%
  select(GEOID, total_area)

# Calculate area by HOLC grade
by_grade <- intersect_sf %>%
  st_drop_geometry() %>%
  group_by(GEOID, !!sym(holc_col)) %>%
  summarise(area_sum = sum(as.numeric(in_area)), .groups = "drop") 

# Add NA category
covered_sf <- by_grade %>%
  group_by(GEOID) %>%
  summarize(covered_area = sum(area_sum), .groups = "drop")


na_rows <- block_area %>%
  left_join(covered_sf, by = "GEOID") %>%
  mutate(covered_area = coalesce(covered_area, 0),
         area_sum = pmax(as.numeric(total_area) - covered_area, 0),
         !!holc_col := "NA") %>%
  select(GEOID, !!sym(holc_col), area_sum)
    
# Combine covered and na holc 
grouped_holc <- bind_rows(by_grade, na_rows) %>%
  left_join(block_area, by = "GEOID") %>%
  mutate(holc_pct_cover = 100 * area_sum / as.numeric(total_area))


# Pivot so each HOLC grade is a column
holc_pct_wide <- grouped_holc %>%
  select(GEOID, holc_grade, holc_pct_cover) %>%
    pivot_wider(names_from = !!sym(holc_col),
              values_from = holc_pct_cover,
              values_fill = 0,
              names_prefix = "pct_holc_")

# Calculate majority and worst holc grades by geoid
holc_pct_major <- holc_pct_wide %>%
  mutate(
    # Majority HOLC: pick column with max pct
    majority_holc = case_when(
      pct_holc_A == pmax(pct_holc_A, pct_holc_B, pct_holc_C, pct_holc_D, pct_holc_NA, na.rm = TRUE) ~ "A",
      pct_holc_B == pmax(pct_holc_A, pct_holc_B, pct_holc_C, pct_holc_D, pct_holc_NA, na.rm = TRUE) ~ "B",
      pct_holc_C == pmax(pct_holc_A, pct_holc_B, pct_holc_C, pct_holc_D, pct_holc_NA, na.rm = TRUE) ~ "C",
      pct_holc_D == pmax(pct_holc_A, pct_holc_B, pct_holc_C, pct_holc_D, pct_holc_NA, na.rm = TRUE) ~ "D",
      pct_holc_NA == pmax(pct_holc_A, pct_holc_B, pct_holc_C, pct_holc_D, pct_holc_NA, na.rm = TRUE) ~ "NA",
      TRUE ~ NA_character_
    ),
    # Worst HOLC: follow severity order (D > C > B > A > NA)
    worst_holc = case_when(
      pct_holc_D > 0 ~ "D",
      pct_holc_C > 0 ~ "C",
      pct_holc_B > 0 ~ "B",
      pct_holc_A > 0 ~ "A",
      pct_holc_NA > 0 ~ "NA",
      TRUE ~ NA_character_
    )
  )

# Set HOLC category levels
category_levels <- c("A", "B", "C", "D", "NA")

# Finalize data by ordering the levels of HOLC data
holc_pct_final <- holc_pct_major %>%
  mutate(majority_holc = factor(majority_holc, levels = category_levels),
         worst_holc = factor(worst_holc, levels = category_levels))


#################################################################################
# Merge final redlining data to acs & tcc data
#################################################################################

combined_sf <- merged_sf %>%
  left_join(holc_pct_final,
            by = "GEOID")

# Check final data with a map
ggplot() +
  geom_sf(data = combined_sf,
          fill = NA,
          colour = "black") +
  geom_sf(data = combined_sf,
          aes(fill = worst_holc)) +
  scale_fill_manual(values =
                      c("olivedrab",
                        "lightblue",
                        "lightgoldenrodyellow",
                        "firebrick",
                        "grey")) +
  theme_void()

# Check final data with a map
ggplot() +
  geom_sf(data = combined_sf,
          fill = NA,
          colour = "black") +
  geom_sf(data = combined_sf,
          aes(fill = majority_holc)) +
  scale_fill_manual(values =
                      c("olivedrab",
                        "lightblue",
                        "lightgoldenrodyellow",
                        "firebrick",
                        "grey")) +
  theme_void()


# Save data for analysis
acs_tcc_holc_final <- combined_sf

save(file = "Data/acs_tcc_holc_final.RData",
     x = acs_tcc_holc_final)

readr::write_csv(file = "Data/acs_tcc_holc_final.csv",
                 x = acs_tcc_holc_final)


# Table 1: Sum Stats
den_sum_stat <- acs_tree_final_op1 %>%
  st_drop_geometry() %>%
  filter(total_pop > 0,
         !is.na(nlcd_mean)) %>%
  group_by(holc_grade) %>%
  reframe("% Hispanic/Latino" = mean(hisp_per),
          "% Non-Hispanic Black" = mean(nhb_per),
          "% Non-Hispanic White" = mean(nhw_per),
          "% Poverty" = mean(pov_per),
          "Total Population" = sum(total_pop),
          "Average Tree Canopy Coverage" = mean(nlcd_mean),
          #"Min TCC" = min(nlcd_min),
          "Sample Size (n)" = n())


# Table: Sum Stats for trees
den_tree_sum_stat <- acs_tree_final_op1 %>%
  st_drop_geometry() %>%
  filter(!is.na(nlcd_mean)) %>%
  dplyr::group_by(holc_grade) %>%
  reframe("mean" = mean(nlcd_mean),
          "med of mean" = median(nlcd_mean),
          "min of mean" = min(nlcd_mean),
          "max of mean" = max(nlcd_mean),
          "avg med" = mean(nlcd_median),
          "avg min" = mean(nlcd_min),
          "avg max" = mean(nlcd_max),
          "sample size (n)" = n()
  )
