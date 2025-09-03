#-------------------------------------------------------------------------------
# title: Data cleaning
# author: E Lunsford  
# date: 2025-09-01
#  
# This code is to create one dataset for analysis by merging together
# population, tree canopy, and redlining data
#
#
# Last Run: 09/01/2025 and code was in working order 
# using R 4.5.1 and RStudio 2025.05.1+513 
#
#-------------------------------------------------------------------------------

# Last update: updated code to match new finalized data from 02

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
library(dplyr)

#################################################################################
# Load data
#################################################################################

# ACS data
load(file = "Data/population_tables.RData")
pop_table2 <- pop_table3

# Tree canopy data
load(file = "Data/den_nad83_tree.RData")

# Redlining data
load(file = "Data/redline_nad83.RData")


#################################################################################
# Combine ACS & TCC data
#################################################################################

# Assign CRS to modern formatting
acs_sf <- sf::st_set_crs(pop_table2, st_crs(pop_table2))
tree_sf <- sf::st_set_crs(den_nad83_tree, st_crs(den_nad83_tree))

# Make sure ACS and tree CRS match
if (st_crs(acs_sf) != st_crs(tree_sf)) {
  tree_sf <- st_transform(tree_sf, st_crs(acs_sf))
}

# Merge tree data into acs data
merged_sf <- acs_sf %>%
  left_join(st_drop_geometry(tree_sf), by = "GEOID") 

nrow(merged_sf)
# 571
length(unique(merged_sf$GEOID))
# 571

#################################################################################
# Prepare redlining data for merge
#################################################################################

# Check CRS of new merged sf and redlining data
redlining_sf <- sf::st_set_crs(redline_nad83, st_crs(redline_nad83))

if (st_crs(merged_sf) != st_crs(redlining_sf)) {
  redlining_sf <- st_transform(redlining_sf, st_crs(merged_sf))
}

# Intersect cbg with holc data
holc_col <- "holc_grade"

intersect_sf <- st_intersection(
  merged_sf %>% dplyr::select(GEOID, geometry),
  redlining_sf %>% dplyr::select(all_of(holc_col))
)

# Calculate intersection area
intersect_sf <- intersect_sf %>%
  dplyr::mutate(in_area = st_area(geometry))

# Calculate total area per block
block_area <- merged_sf %>%
  dplyr::mutate(total_area = st_area(geometry)) %>%
  st_drop_geometry() %>%
  dplyr::select(GEOID, total_area)

# Calculate area by HOLC grade
by_grade <- intersect_sf %>%
  st_drop_geometry() %>%
  dplyr::group_by(GEOID, !!sym(holc_col)) %>%
  summarise(area_sum = sum(as.numeric(in_area)), .groups = "drop") 

# Add NA category
covered_sf <- by_grade %>%
  dplyr::group_by(GEOID) %>%
  summarize(covered_area = sum(area_sum), .groups = "drop")

na_rows <- block_area %>%
  left_join(covered_sf, by = "GEOID") %>%
  dplyr::mutate(covered_area = coalesce(covered_area, 0),
         area_sum = pmax(as.numeric(total_area) - covered_area, 0),
         !!holc_col := "no_holc") %>%
  dplyr::select(GEOID, !!sym(holc_col), area_sum)
    
# Combine covered and na holc 
grouped_holc <- bind_rows(by_grade, na_rows) %>%
  left_join(block_area, by = "GEOID") %>%
  dplyr::mutate(holc_pct_cover = 100 * area_sum / as.numeric(total_area))

# Pivot so each HOLC grade is a column
holc_pct_wide <- grouped_holc %>%
  dplyr::select(GEOID, holc_grade, holc_pct_cover) %>%
  tidyr::pivot_wider(names_from = !!sym(holc_col),
              values_from = holc_pct_cover,
              values_fill = 0,
              names_prefix = "pct_holc_")

# Calculate majority and worst holc grades by geoid
holc_pct_major <- holc_pct_wide %>%
  dplyr::mutate(
    # Majority HOLC: pick column with max pct
    majority_holc = case_when(
      pct_holc_A == pmax(pct_holc_A, pct_holc_B, pct_holc_C, pct_holc_D, pct_holc_no_holc, na.rm = TRUE) ~ "A",
      pct_holc_B == pmax(pct_holc_A, pct_holc_B, pct_holc_C, pct_holc_D, pct_holc_no_holc, na.rm = TRUE) ~ "B",
      pct_holc_C == pmax(pct_holc_A, pct_holc_B, pct_holc_C, pct_holc_D, pct_holc_no_holc, na.rm = TRUE) ~ "C",
      pct_holc_D == pmax(pct_holc_A, pct_holc_B, pct_holc_C, pct_holc_D, pct_holc_no_holc, na.rm = TRUE) ~ "D",
      pct_holc_no_holc == pmax(pct_holc_A, pct_holc_B, pct_holc_C, pct_holc_D, pct_holc_no_holc, na.rm = TRUE) ~ "no_holc",
      TRUE ~ NA_character_
    ),
    # Worst HOLC: follow severity order (D > C > B > A > NA)
    worst_holc = case_when(
      pct_holc_D > 0 ~ "D",
      pct_holc_C > 0 ~ "C",
      pct_holc_B > 0 ~ "B",
      pct_holc_A > 0 ~ "A",
      pct_holc_no_holc > 0 ~ "no_holc",
      TRUE ~ NA_character_
    )
  )

# Set HOLC category levels
category_levels <- c("A", "B", "C", "D", "no_holc")

# Finalize data by ordering the levels of HOLC data
holc_pct_final <- holc_pct_major %>%
  dplyr::mutate(majority_holc = factor(majority_holc, levels = category_levels),
         worst_holc = factor(worst_holc, levels = category_levels))
  
nrow(holc_pct_final)
length(unique(holc_pct_final$GEOID))

#################################################################################
# Merge final redlining data to acs & tcc data
#################################################################################

combined_sf <- merged_sf %>%
  left_join(holc_pct_final,
            by = "GEOID")

nrow(combined_sf)
length(unique(combined_sf$GEOID))

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
save(file = "Data/combined_sf.RData",
     x = combined_sf)

readr::write_csv(file = "Data/combined_sf.csv",
                 x = combined_sf)

#------------------------------------------------------------------------------
sample_test <- combined_sf %>%
  select(GEOID,nlcd_mean,nlcd_sum,nlcd_count,nlcd_median, nlcd_min,nlcd_max,
         test_nlcd_mean,test_nlcd_sum,test_nlcd_count,test_nlcd_median, test_nlcd_min,test_nlcd_max)

summary(sample_test)

# Table 1: Sum Stats
den_sum_stat <- combined_sf %>%
  st_drop_geometry() %>%
  group_by(worst_holc) %>%
  reframe("% Hispanic/Latino" = mean(hisp_per),
          "% Non-Hispanic Black" = mean(nhb_per),
          "% Non-Hispanic White" = mean(nhw_per),
          "% Poverty" = mean(pov_per),
          "Total Population" = sum(total_pop),
          "Zero pop (n)" = (sum(total_pop == 0)),
          "Zero Pop (%)" = (sum(total_pop == 0) / n()) * 100,
          "mean" = mean(nlcd_mean),
          "med of mean" = median(nlcd_mean),
          "min of mean" = min(nlcd_mean),
          "max of mean" = max(nlcd_mean),
          "avg med" = mean(nlcd_median),
          "avg min" = mean(nlcd_min),
          "avg max" = mean(nlcd_max),
          "sample size (n)" = n(),
          "NA count" = sum(is.na(nlcd_mean)),
          "NA %" = (sum(is.na(nlcd_mean)) / n()) * 100)

write.csv(x = den_sum_stat, 
          file = "Outputs/Sum_stats.csv")
