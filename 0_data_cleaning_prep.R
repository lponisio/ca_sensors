# Objective:
# The following script reads in raw data from CA Sensors 2025 field season, cleans the data,
#   saves/exports a copy into the dropbox folder, and generates summary stats and visualizations
#   of these data. The data included are:
#   - 1) bee marking data
#   - 2) floral survey data
#   - 3) climate data **NOT READY/CODED**
#   - 4) camera trap observation data **NOT READY/CODED**

# Set working directory
setwd("~/")
## run the lab_paths script to appropriately set path for your computer (must be set up for each lab member)
## once amended with your computer's information, lab_paths.R can be saved to your home directory
source("lab_paths.R")
local.path
## generate path from your folder to the CA sensors folder
dir.CASensors <- file.path(local.path,"ca_sensors_saved")
## change working directory
setwd(dir.CASensors)

# ---------------- #
#  Load libraries  #
# ---------------- #

library(dplyr)
library(tidyr)
library(readxl)
library(lubridate)
library(tidyverse)
library(MetBrewer)
#install.packages(c("googlesheets4", "googledrive"))
library(googlesheets4)
library(googledrive)
library(vegan)

# ------------------------------------ #
#  Read in and clean bee marking data  #
# ------------------------------------ #

## --- Read in data from drive --- ##
# This will be changed once data is finalized - data will be moved into dropbox ca_sensors_saved
# Authenticate google drive - only once per session
# this will prompt you to login to google
#gs4_auth()

# read in marking data
marks.raw <- read_xlsx("./data/raw/CASensors_2026_beeMarking.xlsx", sheet = "Marking data") %>%
  rename("col.date" = "date",
         "recapture_YN" = "recapture? (y/n)",
         "grid_cell" = "grid_capture_location") %>%
  filter(is.na(col.date) == F, # filter out empty rows attached to bottom of spreadsheet
         is.na(Aruco_num) == F, # filter out missing codes
         site == "EQN") %>%
  mutate(capture_time = as.numeric(capture_time))

# clean up the data (this will be updated as we find more oddities)
## DECISION POINT: all ambiguous vos/calig are put into vosnesenskii until further notice
marks.clean <- marks.raw %>%
  mutate(col.date = ymd(gsub("\\.", "-", col.date)), # reformat date w/lubridate
         recapture_YN = toupper(recapture_YN), # change recaptures to consistent case
         recapture_YN = case_when(recapture_YN == "N *" ~ "N",
                                  .default = recapture_YN),
         floral_host = case_when(floral_host == "N/A" ~ NA,
                                 grepl(pattern = "none", floral_host) == T ~ NA,
                                       .default = floral_host),
         bee_sp_id = case_when(bee_sp_id == "insularis *" ~ "insularis",
                               bee_sp_id == "flavifrons*" ~ "flavifrons",
                               bee_sp_id == "vosnesenskii_caliginosus" ~ "vosnesenskii",
                               bee_sp_id == "melanopygus *" ~ "melanopygus",
                               bee_sp_id == "insularis *" ~ "insularis",
                               bee_sp_id == "vosnesenskii_caliginosus_insularis" ~ "vosnesenskii",
                               bee_sp_id == "No ID" ~ NA,
                               .default = bee_sp_id
                               ),
         Aruco_num = as.numeric(Aruco_num)
  ) %>%
  filter(is.na(Aruco_num) == F,
         is.na(bee_sp_id) == F)

# Write the cleaned data to a csv file. Commented out to avoid overwriting.
write.csv(marks.clean,
          file = "../ca_sensors_saved/data/cleaned/CASensors_BeeMarking2026.csv",
          row.names = F)

# generate a summary of captures by species
bee.sp.summary <- marks.clean %>%
  group_by(bee_sp_id, recapture_YN) %>%
  summarise(n.rows = n()) %>%
  pivot_wider(names_from = recapture_YN, values_from = n.rows, values_fill = 0) %>%
  ungroup() %>%
  rename("Marked" = "N",
         "Recaptured" = "Y") %>%
  arrange(desc(Marked))

# commented out, used for summary stats/basic updates
# write.csv(bee.sp.summary, file = "~/Desktop/CASensors_2026_beeCaptureSummary.csv",
#           row.names = F)

# visualize the bee captures by species across the season
all.bees.agg <- marks.clean %>%
  mutate(case_when(bee_sp_id == "vos/calig" ~ "vosnesenskii",
                   .default = bee_sp_id)) %>%
  filter(caste_sex != "q") %>% # remove queens
  group_by(bee_sp_id, col.date) %>%
  summarize(n.captures = n()) %>% ungroup %>% # get number of captures by species and date
  group_by(bee_sp_id) %>% filter(is.na(bee_sp_id) == F) %>%
  ungroup()

bee.obs.lines <- ggplot(data = all.bees.agg, aes(x = col.date, y = n.captures,
                                                 color = bee_sp_id)) +
  geom_point() +
  geom_line() +
  scale_color_manual(values = MetBrewer::met.brewer("Homer1", length(unique(all.bees.agg$bee_sp_id)))) +
  facet_wrap(vars(bee_sp_id)) +
  theme_bw() +
  labs(x = "Date", y = "Number of Bees Captured", color = "Bombus species") +
  theme(legend.position = "none")

bee.obs.lines

ggsave(plot = bee.obs.lines, units = "in", width = 7, height = 5, device = "png",
       file = "../ca_sensors_saved/figures/CASensors_2026_bombus_obs_summary.png")
# -------------------------------------- #
#  Read in and clean floral survey data  #
# -------------------------------------- #

# Read in floral data
flowers.raw <- read_xlsx("./data/raw/CASensors_2026_FloweringPlantSurveys_Entry.xlsx",
                         sheet = "Floral data") %>%
  rename("col.date" = "date") # rename the date column

## Generate a species list of plants for fixing plant names
# flw.sp <- data.frame(plant_species = unique(flowers.raw$plant_species)) %>%
#   mutate(capital_count = str_count(plant_species, "[A-Z]"),
#          four_plus_caps = capital_count >= 4)
# write.csv(flw.sp, file = "~/Desktop/CASensors_2025_unique_plant_sp.csv", row.names = F)

# Generate a cleaned floral survey dataframe
# not much cleaning to do outside of name changes which were completed outside of R
flowers.clean <- flowers.raw %>%
  mutate(floral_abundance_logBins = case_when(num_flowers == "1" ~ "1-10",
                                              num_flowers == "2" ~ "11-100",
                                              num_flowers == "3" ~ "101-1000",
                                              num_flowers == "4" ~ "1001-10000",
                                              num_flowers == "5" ~ ">10000"),
         num_flowers = as.numeric(num_flowers)
  ) %>%
  filter(is.na(num_flowers) == F,
         site == "EQN")

# Write the cleaned data to a csv file. Commented out to avoid overwriting.
write.csv(flowers.clean,
          file = "./data/cleaned/CASensors_Flowers_clean2026.csv",
          row.names = F)

# --------------------------------- #
#  Read in and clean climatic data  #
# --------------------------------- #

# Load data
climate.raw <- read_xlsx("CASensors2025_BeeMarking_raw.xlsx",
                         sheet = "survey conditions") %>%
  rename("col.date" = "date") # rename the date column
str(climate.raw) 

climate.clean <- climate.raw %>%
  mutate(col.date = ymd(gsub("\\.", "-", col.date)), # reformat date w/lubridate
         across(
           where(is.character),
           ~ na_if(str_to_lower(.x), "na")
           ),
         # Below code is pulled from ChatGPT
         # It just fixes the mixed unit temperature recordings, converting them to celsius
         across(c(start_temp, end_temp),
           ~ {val  <- as.numeric(str_extract(.x, "-?\\d+\\.*\\d*")) # extract the numeric portion 
             unit <- str_to_lower(str_extract(.x, "[a-z]")) # extract the unit character (c or f)
             round(if_else(unit == "f", (val - 32) * 5/9, val), 1) # convert F → C, leave C as-is, then round to 1 decimal
           },
           .names = "{.col}" # overwrite column names
           ),
         across(c(survey_start, survey_end),
                ~{as.numeric(.x)},
                .names = "{.col}"),
         across(
           c(survey_start, survey_end),
           ~ {
             # convert Excel fractional day → POSIXct time (date arbitrary)
             t <- as.POSIXct(.x * 86400,
                             origin = "1970-01-01", tz = "UTC")
             # add 12h to values earlier than 08:00 (field times that were PM)
             t <- ifelse(format(t, "%H:%M") < "08:00",
                         t + 12 * 3600,
                         t)
             # return HH:MM (24-hour) string
             format(as.POSIXct(t, origin = "1970-01-01", tz = "UTC"),
                    "%H:%M")
           },
           .names = "{.col}"),
         site = toupper(site)
         ) %>%
  filter(grepl("missing", x = notes) == F)

# Write the cleaned data to a csv file. Commented out to avoid overwriting.
write.csv(climate.clean,
          file = "../cleaned/CASensors_climate2025_clean.csv",
          row.names = F)

# ------------------------------- #
#  Read in and clean effort data  #
# ------------------------------- #

effort.raw <- read_sheet("https://docs.google.com/spreadsheets/d/1orl3jILZ-L5GrLEYgVyq2-gm7KAX7-KdhIdBmg8HGjk/edit?gid=0#gid=0",
                         col_types = "dDcccTTc") %>%
  select(date, site, surveyor, grid_cell)

effort.clean <- effort.raw %>%
  filter(is.na(date) == F)

# Write the cleaned data to a csv file. Commented out to avoid overwriting.
write.csv(effort.clean,
          file = "./data/cleaned/CASensors_Effort_clean.csv",
          row.names = F)

# --------------------------------------- #
#  Read in and clean camera station data  #
# --------------------------------------- #

cam.station.raw <- read_xlsx("./data/raw/CASensors_2026_cameraStationMetadata.xlsx") %>%
  rename("grid_cell" = "Grid cell",
         "canopy_cover" = "canopy cover")

cam.station.clean <- cam.station.raw %>%
  select(grid_cell, lat, long, canopy_cover)

# Write the cleaned data to a csv file. Commented out to avoid overwriting.
write.csv(cam.station.clean,
          file = "./data/cleaned/CASensors_canopyCover_cleaned.csv",
          row.names = F)
