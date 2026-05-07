rm(list=ls())

library(c("sf", "spsurvey"))

source("lab_paths.R")
setwd("Dropbox (University of Oregon)/ca_sensors")
source("src/grts_prep.R")

## load csv here, colnames should be location, cell, vegetation

# Simulate fake data to test function
## South Fork grid:
## rows A:R, columns 1:9

south_fork_cells <- expand.grid(
  grid_number = 1:9,
  grid_letter = LETTERS[1:18]
)
south_fork_cells$location <- "South Fork"
south_fork_cells$cell <- paste0(
  south_fork_cells$grid_number,
  south_fork_cells$grid_letter
)

## Emerald Queen grid:
## rows A:F, columns 1:9

emerald_queen_cells <- expand.grid(
  grid_number = 1:9,
  grid_letter = LETTERS[1:6]
)
emerald_queen_cells$location <- "Emerald Queen"
emerald_queen_cells$cell <- paste0(
  emerald_queen_cells$grid_number,
  emerald_queen_cells$grid_letter
)

## Combine locations into one sampling frame
cells <- rbind(
  south_fork_cells,
  emerald_queen_cells
)

## Add fake vegetation categories.

cells$vegetation <- sample(
  x = c("forest", "shrub", "serpentine", "riparian"),
  size = nrow(cells),
  replace = TRUE,
  prob = c(0.45, 0.25, 0.20, 0.10)
)

## Keep only the columns needed for sampling.
cells <- cells[, c("location", "cell", "vegetation")]


##  Make the grid selections for both locations
draw_today <- choose_grts_cells_by_location(
  cell_table = cells,
  n_by_location = c(
    "South Fork" = 10,
    "Emerald Queen" = 6
  ),
  cell_col = "cell",
  location_col = "location",
  veg_col = "vegetation",
  stratify_by_vegetation = TRUE,
  vegetation_allocation = "proportional",
  ## seed = 123, ## leave commented out or change each survey
  n_over = 3  ## extra draw in case some cannot be accessed due to
  ## farm activities etc.
)

draw_today$selected[, c(
  "location",
  "cell",
  "vegetation",
  "survey_order_within_location"
)]

table(
  draw_today$selected$location,
  draw_today$selected$vegetation
)
