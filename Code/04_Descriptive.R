#-------------------------------------------------------------------------------
# title: Descriptive exploration
# author: E Lunsford  
# date: 2025-09-03
#  
# This code is to explore that merged final dataset
#
#
# Last Run: 09/03/2025 and code was in working order 
# using R 4.5.1 and RStudio 2025.05.1+513 
#
#-------------------------------------------------------------------------------

# Last update: updated code to match new merged data from 03

#################################################################################
# Load Libraries.
#################################################################################
library(broom)
library(ggplot2)
library(colorspace)
library(ggspatial)
library(gstat)
library(ggthemes)
library(knitr)
library(raster)
library(sf)
library(stringr)
library(terra)
library(tidyverse)

#################################################################################
# Load data
#################################################################################

load(file = "Data/combined_sf.RData")

#################################################################################
# Finalize data
# Options: 
#   1) without DIA
#   2) without DIA and no population
#   3) without DIA, no population, and no HOLC grade
#################################################################################

# Option 1) remove DIA

dia_id <- c("080319800011")

df_op1 <- combined_sf %>%
  filter(!GEOID %in% dia_id)

# Option 2) remove DIA and cbg with 0 population

df_op2 <- combined_sf %>%
  filter(!GEOID %in% dia_id,
         total_pop > 0)

# Option 3) remove DIA, cbg with 0 population, and no HOLC grade

holc_grade_id <- c("NA")

df_op3 <- combined_sf %>%
  filter(!GEOID %in% dia_id,
         total_pop > 0,
         !worst_holc == "no_holc")

#################################################################################
# Set map theme.
#################################################################################
map_theme <- theme(
  #aspect.ratio = 1,
  text  = element_text(size = 12, color = 'black'),
  #panel.spacing.y = unit(0,"cm"),
  #panel.spacing.x = unit(0.25, "lines"),
  panel.grid.minor = element_line(color = "transparent"),
  panel.grid.major = element_line(color = "transparent"),
  panel.border = element_blank(),
  panel.background=element_blank(),
  axis.ticks = element_blank(),
  axis.text = element_blank(),
  # legend.position = c(0.1,0.1),
  # plot.margin = grid::unit(c(0,0,0,0), "mm"),
  legend.key = element_blank(),
  #legend.background = element_rect(fill='transparent'),
  plot.margin = unit(c(0, 0, 0, 0), "cm"),
  panel.spacing = unit(c(0, 0, 0, 0), "cm"),
  legend.position = "inside",
  legend.position.inside = c(0.84, 0.4),
  legend.title = element_text(size = 12),
  legend.background = element_rect(fill = "white"),
  legend.spacing.y = unit(0.2, 'cm'),
  legend.margin = margin(0.1,0,0,0, unit="cm")
)

den_col <- "darkgrey"
ges_col <- "#e41a1c"
oth_col <- "#4daf4a"
hwy_col <- "darkblue"

blue <- "#377eb8"
red <- "#e41a1c"
green <- "#4daf4a"
purple <- "#984ea3"
orange <- "#ff7f00"  
yellow <- "#ffff33"

#################################################################################
# Plot 1: Tree Canopy
#################################################################################
ggplot() +
  geom_sf(data = df_op1,
          inherit.aes = F,
          mapping = aes(fill = nlcd_mean)) +
  labs(title = "2021 Tree Canopy",
       x = "",
       y = "",
       fill = "Percent") +
  scale_fill_continuous_sequential(palette = "Viridis") +
  map_theme +
  annotation_scale(location = "br") +
  annotation_north_arrow(aes(location = "tl"),
                         height = unit(1, "cm"), 
                         width = unit(1, "cm"),
                         pad_x = unit(0.5, "cm"),
                         pad_y = unit(1, "cm")) +
  theme(plot.margin = margin(t = 10,  # Top margin
                             r = 20,  # Right margin
                             b = 10,  # Bottom margin
                             l = 10)) # Left margin


#################################################################################
# Plot 2: Present day CBG with HOLC
#################################################################################

ggplot() +
  geom_sf(data = df_op1,
          inherit.aes = F,
          fill = NA,
          linewidth = 1) +
  geom_sf(data = df_op1,
          inherit.aes = F,
          aes(fill = worst_holc)) +
  scale_fill_manual(values =
                      c("olivedrab",
                        "lightblue",
                        "lightgoldenrodyellow",
                        "firebrick",
                        "grey")) +
  coord_sf(expand = F) +
  map_theme +
  annotation_scale(location = "br") +
  annotation_north_arrow(aes(location = "tl"),
                         height = unit(1, "cm"), 
                         width = unit(1, "cm"),
                         pad_x = unit(0.5, "cm"),
                         pad_y = unit(1, "cm")) +
  theme(plot.margin = margin(t = 10,  # Top margin
                             r = 20,  # Right margin
                             b = 10,  # Bottom margin
                             l = 10)) + # Left margin
  labs(title = "Former Redlining Status (HOLC Grades) by Census Block Groups",
       fill = "HOLC \n Grade") +
  theme(plot.margin = margin(t = 20,  # Top margin
                             r = 50,  # Right margin
                             b = 40,  # Bottom margin
                             l = 10)) # Left margin


#################################################################################
# Explore the data
#################################################################################

# Normality plots
set.seed(123)
sample_n(df_op2, 10)

ggplot() +
  geom_histogram(data = df_op2,
                 mapping = aes(x = nlcd_mean,
                               y = after_stat(density)),
                 color = 1,
                 fill = "grey") +
  geom_density() +
  theme_clean() +
  labs(x = "Mean TCC",
       y = "Density")


ggplot(data = df_op2,
       mapping = aes(sample = nlcd_mean)) +
  geom_qq() +
  geom_qq_line() +
  theme_bw() +
  labs(title = "Tree Canopy Coverage (%)",
       x = "Theoretical",
       y = "Percent") +
  theme(plot.margin = margin(t = 20,  # Top margin
                             r = 10,  # Right margin
                             b = 20,  # Bottom margin
                             l = 10)) # Left margin

ggplot(data = df_op2,
       mapping = aes(sample = nlcd_mean,
                     colour = factor(worst_holc))) +
  geom_qq() +
  geom_qq_line() +
  theme_bw() +
  labs(title = "Tree Canopy Coverage (%) by HOLC Grade",
       x = "Theoretical",
       y = "Percent",
       color = "HOLC Grade") +
  theme(plot.margin = margin(t = 20,  # Top margin
                             r = 10,  # Right margin
                             b = 20,  # Bottom margin
                             l = 10)) # Left margin 

# Shapiro-Wilks for normality
shapiro.test(df_op2$nlcd_mean)

# Levene Test
library(car)

leveneTest(nlcd_mean ~ worst_holc,
           data = df_op2)
           
# Boxplot
ggplot(data = df_op2,
       mapping = aes(x = worst_holc,
                     y = nlcd_mean)) +
  geom_boxplot(aes(fill = worst_holc),
               show.legend = F) +
  scale_fill_manual(values = c("olivedrab",
                               "lightblue",
                               "lightgoldenrod",
                               "firebrick",
                               "grey")) +
  labs(x = "HOLC Grade",
       y = "Percent Tree Canopy",
       title = "Distribution of tree canopy by HOLC grades",
       caption = "Sources: Mapping Inequality, NLCD") +
  theme(text = element_text(size = 10)) +
  theme_few() +
  theme(plot.margin = margin(t = 20,  # Top margin
                             r = 10,  # Right margin
                             b = 20,  # Bottom margin
                             l = 10)) # Left margin

# Hispanic/Latino
# Levene Test for percent Hispanic/Latino ~ holc grade
leveneTest(hisp_per ~ worst_holc,
           data = df_op2)

# Boxplot of distribution
boxplot(hisp_per ~ worst_holc,
        data = df_op2,
        main = "Percent hisp by HOLC Grade",
        xlab = "HOLC",
        ylab = "hisp",
        col = "steelblue",
        border = "black")

# NHB
# Levene Test for percent non-Hispanic Black ~ holc grade
leveneTest(nhb_per ~ worst_holc,
           data = df_op2)

# Boxplot of distribution
boxplot(nhb_per ~ worst_holc,
        data = df_op2,
        main = "Percent non-hisp black by HOLC Grade",
        xlab = "HOLC",
        ylab = "non-hisp black",
        col = "steelblue",
        border = "black")

# NHW
# Levene Test for percent poverty ~ holc grade
leveneTest(nhw_per ~ worst_holc,
           data = df_op2)

# Box plot of distribution
boxplot(nhw_per ~ worst_holc,
        data = df_op2,
        main = "Percent NHW by HOLC Grade",
        xlab = "HOLC",
        ylab = "NHW",
        col = "steelblue",
        border = "black")

# Poverty
# Levene Test for percent non-Hispanic white ~ holc grade
leveneTest(pov_per ~ worst_holc,
           data = df_op2)

# Box plot of distribution
boxplot(pov_per ~ worst_holc,
        data = df_op2,
        main = "Percent poverty by HOLC Grade",
        xlab = "HOLC",
        ylab = "poverty",
        col = "steelblue",
        border = "black")


# Explore some more
ggplot(data = df_op2,
       mapping = aes(x = nlcd_mean,
                     y = hisp_per)) +
  geom_jitter(aes(color = worst_holc),
              size = 2) +
  scale_color_manual(values = c("olivedrab",
                                "lightblue",
                                "lightgoldenrod",
                                "darkred",
                                "grey"))+
  labs(x = "Percent Tree Canopy",
       y = "Percent Hispanic/Latino",
       color = "HOLC Grade",
       title = "Tree Canopy \n by Percent Hispanic/Latino and HOLC Grade",
       caption = "Source: Mapping Inequality, 5-year American Community Survey") +
  theme(plot.margin = margin(t = 20,  # Top margin
                             r = 10,  # Right margin
                             b = 20,  # Bottom margin
                             l = 10)) # Left margin 

