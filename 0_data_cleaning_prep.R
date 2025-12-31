# load libraries
library(tidyverse)
library(readxl)
library(lubridate)

# Set working directory
setwd("~/")
## run the lab_paths script to appropriately set path for your computer (must be set up for each lab member)
## once amended with your computer's information, lab_paths.R can be saved to your home directory
source("lab_paths.R")

## generate path from your folder to the BMMA folder
dir.CASensors <- file.path("~/University\ of\ Oregon\ Dropbox/Nevin\ Cullen/CA Sensors")
## change working directory
setwd(dir.CASensors)

# for now we'll just setwd to a local directory until Lauren Makes a CA sensors folder on DB
setwd("~/Dropbox (Personal)/ca_sensors_saved/data/raw/")



# ------------------------------------ #
#  Read in and clean bee marking data  #
# ------------------------------------ #

# Read in the data
sensors.raw <- read_xlsx("CASensors2025_BeeMarking_raw.xlsx", sheet = "data sheet") %>%
  rename("col.date" = "date",
         "multiple_capture_YN" = "multiple_capture?")

# inspect the bee marking data
str(sensors.raw)

# clean up the data a little bit (this will be updated as we find more oddities)
sensors.clean <- sensors.raw %>%
  mutate(col.date = ymd(gsub("\\.", "-", col.date)), # reformat date w/lubridate
         recapture_YN = toupper(recapture_YN), # change recaptures to consistent case
         multiple_capture_YN = toupper(
           case_when(is.na(multiple_capture_YN) ~ "N",
                     .default = multiple_capture_YN) # add N values to the recapture column
         ),
         grid_cell_final = case_when(grid_cell_final == "n/a" ~ NA,
                                     grid_cell_final == "na" ~ NA,
                                     .default = grid_cell_final),
         flower_sp_visited = case_when(flower_sp_visited == "NA*" ~ NA,
                                       flower_sp_visited == "na" ~ NA,
                                       .default = flower_sp_visited),
         capture_time = as.numeric(case_when(capture_time == "na" ~ NA,
                                             .default = capture_time)),
         # This next bit of mutate code is from ChatGPT to fix capture times
         # It converts the decimal time to hh:mm and corrects PM captures
         capture_time_posix = as.POSIXct(capture_time * 86400,
                                         origin = "1970-01-01", tz = "UTC"),
         # adjust any time before 08:00 (really PM field times) by adding 12 hours
         capture_time_posix = if_else(
           format(capture_time_posix, "%H:%M") < "08:00",
           capture_time_posix + 12 * 3600,
           capture_time_posix),
         capture_time = format(capture_time_posix, "%H:%M")
  ) %>%
  select(-capture_time_posix, -start_time, -bee_sp_id,
         -caste_sex, -notes, -notes_update, -og_grid_location) %>% # remove columns not used for analysis
  filter(is.na(bee_sp_id_Final) == F, # remove all unidentifiable bee species
         grepl("\\?", x = qr_num1) == F, # remove any oddball qr codes
  )

# Write the cleaned data to a csv file. Commented out to avoid overwriting.
write.csv(sensors.clean,
          file = "~/Desktop/Ponisio Lab/CA Sensors/data/CASensors_BeeMarking2025_clean.csv",
          row.names = F)

# generate a summary of captures by species
bee.sp.summary <- sensors.clean %>%
  group_by(bee_sp_id_Final, recapture_YN, site) %>%
  summarise(n.rows = n())

# -------------------------------------- #
#  Read in and clean floral survey data  #
# -------------------------------------- #

# Load data
flowers.raw <- read_xlsx("CASensors2025_FlowerSurveys_raw.xlsx",
                         sheet = "data sheet")%>%
  rename("col.date" = "date") # rename the date column
str(flowers.raw) 

flowers.clean <- flowers.raw %>%
  rename("abund_code" = "num_flowers") %>%
  mutate(col.date = ymd(gsub("\\.", "-", col.date)), # reformat date w/lubridate
         abund_code = case_when(tolower(abund_code) == "present" ~ "1",
                                .default = abund_code),
         floral_abundance_logBins = case_when(abund_code == "1" ~ "1-10",
                                              abund_code == "2" ~ "11-100",
                                              abund_code == "3" ~ "101-1000",
                                              abund_code == "4" ~ "1001-10000",
                                              abund_code == "5" ~ ">10000")
  ) %>%
  filter(plant_species != "*see photos") %>%
  select(!notes, notes, -c(provisional_grid_location, ORIG_SPP_CODE))

# Write the cleaned data to a csv file. Commented out to avoid overwriting.
write.csv(flowers.clean,
          file = "../cleaned/CASensors_Flowers_clean.csv",
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