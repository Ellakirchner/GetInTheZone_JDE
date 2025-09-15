#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Author: Ella Kirchner 
# Project name: RS Get in the Zone
# Date Last Updated: # Thu Sep  4 11:51:43 2025
# Purpose: Generate all input datasets required for Zoning paper analysis
# Input: Various datasets that are not published here to ensure data privacy. 
#         The input datasets are all saved in one folder and jointly loaded into R (see line 88 onwards).
# Output Files: 
          # - 6 .RData files with spatial data based on oaf crop cut locations (ndvi, evi, temp, chirp, soil, all with suffix "_oaf")
          # - 2 .RData files with temperature data based on oaf crop cut locations, meantemp and gdd separately
          # - 6 .RData files with spatial data based on 50k geopoints on maize crop land (ndvi, evi, temp, chirp, soil, all with suffix "_50k")
          # - 6 .RData files with spatial data based on both (oaf + 50k) (ndvi, evi, temp, chirp, soil, all with suffix "_2")
          # - Raw oaf data (.RData)
          # - (preprocessed) OAF data, with admin units merged to it and obs filtered for analysis (.RData)
          # - 4 .RData files on admin data (level 1 and level 3, each for full country and clipped to ROI)
          # - 2 .RData files on temperature to do descriptive statistics on temp
# ReadMe: 1) Set up a project folder and adapt the filepath (line 32, placeholder) to that location
#         2) In the project folder, make a subfolder called "raw_data" and save all input files from GDrive in that data folder 
#         3) In the project folder, make a subfolder called "datasets" (which will be used to save the processed datasets)
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# %%%%%%%%%%%%%%%%%%%%%%%%%%%
# Set Up 
# %%%%%%%%%%%%%%%%%%%%%%%%%%%

# 1. Clean out workspace
rm(list = ls(all = TRUE))

# 2. Set Up File path & directory: 

# This should be the location of the project folder
file_path<- "put_your_file_path_here"
  
filepath <- function(x){paste0(pathfile, x)}

setwd(file_path)

# Save all the datasets from the GDrive folder in the raw_data folder in your project folder
raw_data <- paste0(file_path, "/raw_data")

# Create a folder called "datasets" in your project folder for the processed datasets
output_data <- paste0(file_path, "/datasets")

# 3. Set Seed for replication 
set.seed(123456789)

# %%%%%%%%%%%%%%%%%%%%%%%%%%%
# Load libraries 
# %%%%%%%%%%%%%%%%%%%%%%%%%%%
library(dplyr)
library(terra)
library(sp)
library(sf)
library(proxy)
library(geodata) 
library(janitor)
library(ggplot2)
library(tidyverse)
library(lubridate)

# %%%%%%%%%%%%%%%%%%%%%%%%%%%
# 1 Load Data ---------
# %%%%%%%%%%%%%%%%%%%%%%%%%%%

# 1.1 Kenya admin data ----
# (only admin3 needed as it also contains all information from admin1 and admin2):

admin1_kenya <- st_as_sf(gadm("GADM", level=1, country="KEN"))
admin1_kenya <- clean_names(admin1_kenya)# Clean column names
admin1_kenya$name_concat <- paste(admin1_kenya$name_1, admin1_kenya$name_2, admin1_kenya$name_3, sep = "_")# Create a new variable by concatenating columns

admin3_kenya <- st_as_sf(gadm("GADM", level=3, country="KEN"))
admin3_kenya <- clean_names(admin3_kenya)# Clean column names
admin3_kenya$name_concat <- paste(admin3_kenya$name_1, admin3_kenya$name_2, admin3_kenya$name_3, sep = "_")# Create a new variable by concatenating columns

# Define ROI
# Focus on the following counties (admin1)
districts_to_keep <- c("Bungoma", "Busia", "Homa Bay", "Kakamega", 
                       "Kericho", "Kisumu", "Kisii", "Migori", "Nandi", "Nyamira", 
                       "Siaya", "Trans Nzoia", "Uasin Gishu", "Vihiga")

# Clip Admin datasets to ROI
class(admin1_kenya)
admin1_kenya <- st_as_sf(admin1_kenya)
admin1_clipped <- admin1_kenya %>% filter(name_1 %in% districts_to_keep)
admin3_kenya <- st_as_sf(admin3_kenya)
admin3_clipped <- admin3_kenya %>% filter(name_1 %in% districts_to_keep)

# 1.2 Read in all raw datafiles ---- 
# List all csv files in the folder
files <- list.files(path = raw_data, pattern = "\\.csv$", full.names = TRUE)

# Loop through each file and read it into the global environment
for (file in files) {
  name <- tools::file_path_sans_ext(basename(file))  # Get filename without extension
  data <- read.csv(file)                             # Read the CSV file
  assign(name, data, envir = .GlobalEnv)             # Create a variable in the environment
}

# 1.3 OAF Data ----
if (exists("oaf_new_data_kenya") && ncol(oaf_new_data_kenya) == 118) {
  oaf_raw <- oaf_new_data_kenya
} else {
  oaf_raw <- read.csv(paste0(raw_data, "/oaf_new_data_kenya.csv"), header = TRUE, sep = ";")
}

# reduce dataset to required variables and use it as working dataset (named: oaf): 
oaf <- oaf_raw[c("unique_id", "district", "yield_kg_pa", "yield_kg_ph", "year", "field_latitude", "field_longitude", "intercrop", "plot_acres", "plot_hectares")]

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# A: Prepare data ---- 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# 0) Clean OAF data ----
# Clean OAF data:
# - Filter for regions 
# - Filter for Yield>0 

# Remove rows with missing values in columns "long" and "lat" from working dataset (named: oaf)
oaf <- oaf[!is.na(oaf$field_longitude) & !is.na(oaf$field_latitude), ]

# Drop all observations that are duplicates in terms of lat, lon and year: 
oaf_dup <- oaf %>%
  group_by(field_latitude, field_longitude, year) %>%
  mutate (n_samelatlongyear = n(),
          has_dupes = ifelse(n_samelatlongyear>=2, 1,0),
          id_latlonyear=paste0(field_latitude, "_", field_longitude, "_",year)) %>%
  ungroup() %>%
  arrange (-n_samelatlongyear, id_latlonyear)

oaf <- left_join(oaf, oaf_dup %>% select(unique_id, n_samelatlongyear, id_latlonyear), by="unique_id")
tabyl(oaf$n_samelatlongyear)

oaf <- oaf %>% filter(n_samelatlongyear==1)
oaf <- oaf %>% select(-n_samelatlongyear)

# Add Admin data to oaf data 
# Assign projection system to the dataset (same crs as for the geodata that are to be merged): 
desired_crs <- "+proj=longlat +datum=WGS84 +no_defs"

oaf_geo <- st_as_sf(oaf, coords = c("field_longitude", "field_latitude"), crs = desired_crs)
st_crs(oaf_geo) 
st_crs(admin3_kenya)

oaf <- st_join(oaf_geo, admin3_kenya, join=st_within)

# From the admin dataset, drop excess information:
oaf <- oaf %>%
  select(-nl_name_1, -nl_name_2, -nl_name_3, -varname_3,  
         -engtype_3, -hasc_3, -gid_0, -gid_1, -gid_2, -gid_3, -cc_3, -type_3)

# Filter the working dataset to keep only the desired districts 
oaf <- oaf %>% filter(name_1 %in% districts_to_keep)
table(oaf$name_1)

# Check how many of the obs with zero yield are doing intercropping:
oaf_zero <- oaf %>% 
  mutate(yield = yield_kg_ph/1000) %>%
  filter(yield==0)
table(oaf_zero$intercrop)


# Drop observations with crop cut values of zero 
oaf <- oaf %>%
  mutate(yield = yield_kg_ph/1000)%>%
  filter(yield>0)

# 1) Merge datasets by indicator and point database ----
# point database refers to whether the spatial data were extracted at the 50k random geopoints or at the exact locations of the OAF crop cuts 
# we aim for separate datasets per point database and a joined dataset

merge_by_keywords <- function(keywords) {
  # Find matching object names
  all_objects <- ls(envir = .GlobalEnv)
  matching_names <- all_objects
  
  # Filter for objects containing all keywords
  for (keyword in keywords) {
    matching_names <- grep(keyword, matching_names, value = TRUE)
  }
  
  print("Datasets matching all keywords:")
  print(matching_names)
  
  # Get the first dataset as base and remove duplicates
  base_name <- matching_names[1]
  print(paste("Using", base_name, "as base dataset"))
  
  base_df <- get(base_name)
  
  # Check for and report duplicates in base dataset
  n_total <- nrow(base_df)
  n_unique <- length(unique(base_df$.geo))
  n_duplicates <- n_total - n_unique
  
  print(paste("Base dataset has", n_duplicates, "duplicate .geo values out of", n_total, "rows"))
  
  # Keep only the first occurrence of each unique .geo value in base dataset
  base_df <- base_df %>% 
    group_by(.geo) %>%
    slice(1) %>%
    ungroup()
  
  print(paste("Reduced base dataset from", n_total, "to", nrow(base_df), "rows"))
  
  # Join all other datasets to the base using the .geo column
  for (i in 2:length(matching_names)) {
    current_name <- matching_names[i]
    print(paste("Joining:", current_name))
    
    current_df <- get(current_name)
    
    # Keep only the first occurrence of each .geo value in the current dataset too
    # This prevents many-to-one joins
    current_df <- current_df %>%
      group_by(.geo) %>%
      slice(1) %>%
      ungroup() %>%
      select(-matches("system.index|classification|constant|unique"))
    
    # Join to base dataframe
    base_df <- base_df %>%
      left_join(current_df, by = ".geo")
    print(paste("Completed joining:", current_name))
  }
  
  # Clean up the final dataframe
  base_df <- base_df %>%
    select(-matches("system.index|classification|constant|unique"))
  
  # Format column names as dates
  date_variable_mask <- grepl("^X\\d{4}\\.\\d{2}\\.\\d{2}$", names(base_df))
  names(base_df) <- ifelse(
    date_variable_mask, 
    as.character(as.Date(substring(names(base_df), 2, 11), format='%Y.%m.%d')),
    names(base_df)
  )
  
  # Filter for growing season
  base_df <- base_df %>%
    select(-matches("^\\d{4}-(01|02|03|10|11|12)-\\d{2}$")) 
    # we do not filter directly for the season of interest because we need the dates before and after for potential imputations of NAs
    # we aggregate temp and chirps data anyway, soil data is time-invariant and for evi and ndvi we drop April and September after replacing NAs
        
  # Convert to spatial format
  base_df$geometry <- geojsonsf::geojson_sf(base_df$.geo)$geometry
  base_df <- base_df %>% select(-.geo)
  base_sf <- st_as_sf(base_df)
  print("Converted to spatial format")
  
  return(base_sf)
}

merge_by_keywords_oaf <- function(keywords) {
  # Find matching object names
  all_objects <- ls(envir = .GlobalEnv)
  matching_names <- all_objects
  
  # Filter for objects containing all keywords
  for (keyword in keywords) {
    matching_names <- grep(keyword, matching_names, value = TRUE)
  }
  
  print("Datasets matching all keywords:")
  print(matching_names)
  
  # Get the first dataset as base and remove duplicates
  base_name <- matching_names[1]
  print(paste("Using", base_name, "as base dataset"))
  
  base_df <- get(base_name)
  
  # Check for and report duplicates in base dataset
  n_total <- nrow(base_df)
  n_unique <- length(unique(base_df$unique_id))
  n_duplicates <- n_total - n_unique
  
  print(paste("Base dataset has", n_duplicates, "duplicate unique_id values out of", n_total, "rows"))
  
  # Keep only the first occurrence of each unique .geo value in base dataset
  base_df <- base_df %>% 
    group_by(unique_id) %>%
    slice(1) %>%
    ungroup()
  
  print(paste("Reduced base dataset from", n_total, "to", nrow(base_df), "rows"))
  
  # Join all other datasets to the base using the .geo column
  for (i in 2:length(matching_names)) {
    current_name <- matching_names[i]
    print(paste("Joining:", current_name))
    
    current_df <- get(current_name)
    
    # Keep only the first occurrence of each .geo value in the current dataset too
    # This prevents many-to-one joins
    current_df <- current_df %>%
      group_by(unique_id) %>%
      slice(1) %>%
      ungroup() %>%
      select(-matches("system.index|classification|constant|geo"))
    
    # Join to base dataframe
    base_df <- base_df %>%
      left_join(current_df, by = "unique_id")
    print(paste("Completed joining:", current_name))
  }
  
  # Clean up the final dataframe
  base_df <- base_df %>%
    select(-matches("system.index|classification|constant"))
  
  # Format column names as dates
  date_variable_mask <- grepl("^X\\d{4}\\.\\d{2}\\.\\d{2}$", names(base_df))
  names(base_df) <- ifelse(
    date_variable_mask, 
    as.character(as.Date(substring(names(base_df), 2, 11), format='%Y.%m.%d')),
    names(base_df)
  )
  
  # Filter for growing season
  base_df <- base_df %>%
    select(-matches("^\\d{4}-(01|02|03|10|11|12)-\\d{2}$")) 
  # we do not filter directly for the season of interest because we need the dates before and after for potential imputations of NAs
  # we aggregate temp and chirps data anyway, soil data is time-invariant and for evi and ndvi we drop April and September after replacing NAs
  
  # Convert to spatial format
  base_df$geometry <- geojsonsf::geojson_sf(base_df$.geo)$geometry
  base_df <- base_df %>% select(-.geo)
  base_sf <- st_as_sf(base_df)
  print("Converted to spatial format")
  
  return(base_sf)
}


# Apply the function to each keyword-pointdatabase combination: 
chirp_oaf <- merge_by_keywords_oaf("chirpstack_oaf")
evi_oaf <- merge_by_keywords_oaf("evistack_oaf") 
ndvi_oaf <- merge_by_keywords_oaf("ndvistack_oaf") 
gdd_oaf <- merge_by_keywords_oaf("gdd_oaf_2")
chirp_50k <- merge_by_keywords("chirpstack_50k") 
evi_50k <- merge_by_keywords("evistack_50k") 
ndvi_50k <- merge_by_keywords("ndvistack_50k") 
gdd_50k <- merge_by_keywords("gdd_50k_2")
temp_50k <- merge_by_keywords("meantemp_2") 


# 2) Prep Temp data ----
# Variable names of temp and gdd differ, hence, these datasets were not yet filtered for the season 
# For Temp and GDD, aggregation to higher level required as daily values result in too many variables
# As datasets also vary between point database (oaf vs 50k), we implement two similar (but different) approaches
# We start with the 50k datapoints (as we will need the oaf datasets for the table calculated in #11)

# 2.1) 50k geopoints ----
names(gdd_50k) <- names(gdd_50k) %>%
  gsub("_GDD$|^X", "", .) %>%
  # For date variables (8 digits), convert to YYYY-MM-DD format
  ifelse(grepl("^\\d{8}$", .), 
         format(as.Date(., format = "%Y%m%d"), "%Y-%m-%d"),
         .)
names(temp_50k) <- names(temp_50k) %>%
  gsub("_temperature_2m$|^X", "", .) %>%
  # For date variables (8 digits), convert to YYYY-MM-DD format
  ifelse(grepl("^\\d{8}$", .), 
         format(as.Date(., format = "%Y%m%d"), "%Y-%m-%d"),
         .)

# Check for missings before aggregating:
anyNA(gdd_50k) # no, there are no NAs
# Check for missings before aggregating:
anyNA(temp_50k) # no, there are no NAs

# Aggregate GDD variables over 10 day intervals: 
gdd_short <- gdd_50k
# Drop geometry function of the geometry column to be able to do the aggregation 
# (The geometry column will be added again after the aggregation)
gdd_short_df <- gdd_short %>% 
  st_drop_geometry() %>%
  mutate(geometry = st_geometry(gdd_short))  
# Aggregate over intervals: 
for (year in 2000:2020) {
  month05_a <- paste0(year, "-05-", sprintf("%02d", 1:15))
  month05_b <- paste0(year, "-05-", sprintf("%02d", 16:31))
  month06_a <- paste0(year, "-06-", sprintf("%02d", 1:15))
  month06_b <- paste0(year, "-06-", sprintf("%02d", 16:30))
  month07_a <- paste0(year, "-07-", sprintf("%02d", 1:15))
  month07_b <- paste0(year, "-07-", sprintf("%02d", 16:31))
  month08_a <- paste0(year, "-08-", sprintf("%02d", 1:15))
  month08_b <- paste0(year, "-08-", sprintf("%02d", 16:31))
  
  gdd_short_df <- gdd_short_df %>% 
    mutate(
      !!paste0("gdd_", year, "_05_a") := rowSums(select(., month05_a)),
      !!paste0("gdd_", year, "_05_b") := rowSums(select(., month05_b)),
      !!paste0("gdd_", year, "_06_a") := rowSums(select(., month06_a)),
      !!paste0("gdd_", year, "_06_b") := rowSums(select(., month06_b)),
      !!paste0("gdd_", year, "_07_a") := rowSums(select(., month07_a)),
      !!paste0("gdd_", year, "_07_b") := rowSums(select(., month07_b)),
      !!paste0("gdd_", year, "_08_a") := rowSums(select(., month08_a)),
      !!paste0("gdd_", year, "_08_b") := rowSums(select(., month08_b)),
    )
}

# Keep only the aggregated columns 
gdd_short_df <- gdd_short_df %>%
  select(geometry, ends_with("_a"), ends_with("_b")) #, ends_with("_c"))

# Add geometry column again
gdd_short <- st_as_sf(gdd_short_df, sf_column_name = "geometry")

# Aggregate MeanTemp Variables over 10 day intervals: 
meantemp_short <- temp_50k # generate new dataset for the aggregated variables 
# Drop geometry for aggregation (geometry will later be added again)
meantemp_short_df <- meantemp_short %>% 
  st_drop_geometry() %>%
  mutate(geometry = st_geometry(meantemp_short))
# Aggregate columns
for (year in 2000:2020) {
  month05_a <- paste0(year, "-05-", sprintf("%02d", 1:15))
  month05_b <- paste0(year, "-05-", sprintf("%02d", 16:31))
  month06_a <- paste0(year, "-06-", sprintf("%02d", 1:15))
  month06_b <- paste0(year, "-06-", sprintf("%02d", 16:30))
  month07_a <- paste0(year, "-07-", sprintf("%02d", 1:15))
  month07_b <- paste0(year, "-07-", sprintf("%02d", 16:31))
  month08_a <- paste0(year, "-08-", sprintf("%02d", 1:15))
  month08_b <- paste0(year, "-08-", sprintf("%02d", 16:31))
  
  meantemp_short_df <- meantemp_short_df %>% 
    mutate(
      !!paste0("mean_", year, "_05_a") := rowMeans(select(., month05_a)),
      !!paste0("mean_", year, "_05_b") := rowMeans(select(., month05_b)),
      !!paste0("mean_", year, "_06_a") := rowMeans(select(., month06_a)),
      !!paste0("mean_", year, "_06_b") := rowMeans(select(., month06_b)),
      !!paste0("mean_", year, "_07_a") := rowMeans(select(., month07_a)),
      !!paste0("mean_", year, "_07_b") := rowMeans(select(., month07_b)),
      !!paste0("mean_", year, "_08_a") := rowMeans(select(., month08_a)),
      !!paste0("mean_", year, "_08_b") := rowMeans(select(., month08_b)),
    )
}

# Keep only the aggregated columns 
meantemp_short_df <- meantemp_short_df %>%
  select(geometry, ends_with("_a"), ends_with("_b"), ends_with("_c"))

# Add geometry column again
meantemp_short <- st_as_sf(meantemp_short_df, sf_column_name = "geometry")

# Merge aggregated datasets of GDD, Meantemp: 
temp_50k <- st_join(meantemp_short, gdd_short, left = TRUE)

# 2.2) oaf crop cut locations ---- 
# Check for missings before aggregating:
    anyNA(gdd_oaf) # yes, there are NAs 
    # Check magnitude of missings
    colSums(is.na(gdd_oaf)) # only 1NA 
    # Check which row (is it always the same?)
    which(!complete.cases(st_drop_geometry(gdd_oaf)))# yes, all NAs from one row
    # Drop rows for which all variables, except for the geometry, are NA:
    gdd_oaf <- gdd_oaf %>%
      filter(!if_all(-c(geometry, unique_id), is.na))
    # Check if that was the case for the observation with the missing values: 
    anyNA(gdd_oaf) #yes, no NA left. 

# Separate out variables by indicator (GDD vs mean temp) before aggregating: 
gdd <- gdd_oaf %>% 
  select(geometry, unique_id, ends_with("_GDD"))
meantemp <- gdd_oaf %>% 
  select(geometry, unique_id, ends_with("_temperature_2m"))

# Rename variables into date format to be able to filter for them later
names(gdd) <- names(gdd) %>%
  gsub("_GDD$|^X", "", .) %>%
  # For date variables (8 digits), convert to YYYY-MM-DD format
  ifelse(grepl("^\\d{8}$", .), 
         format(as.Date(., format = "%Y%m%d"), "%Y-%m-%d"),
         .)
names(meantemp) <- names(meantemp) %>%
  gsub("_temperature_2m$|^X", "", .) %>%
  # For date variables (8 digits), convert to YYYY-MM-DD format
  ifelse(grepl("^\\d{8}$", .), 
         format(as.Date(., format = "%Y%m%d"), "%Y-%m-%d"),
         .)

# Aggregate GDD variables over 2 intervals per month: 
gdd_short <- gdd
# Drop geometry function of the geometry column to be able to do the aggregation 
# (The geometry column will be added again after the aggregation)
gdd_short_df <- gdd_short %>% 
  st_drop_geometry() %>%
  mutate(geometry = st_geometry(gdd_short))  
# Aggregate over intervals: 
for (year in 2000:2020) {
  #year_str <- sprintf("%02d", year - 2000)
  month05_a <- paste0(year, "-05-", sprintf("%02d", 1:15))
  month05_b <- paste0(year, "-05-", sprintf("%02d", 16:31))
  month06_a <- paste0(year, "-06-", sprintf("%02d", 1:15))
  month06_b <- paste0(year, "-06-", sprintf("%02d", 16:30))
  month07_a <- paste0(year, "-07-", sprintf("%02d", 1:15))
  month07_b <- paste0(year, "-07-", sprintf("%02d", 16:31))
  month08_a <- paste0(year, "-08-", sprintf("%02d", 1:15))
  month08_b <- paste0(year, "-08-", sprintf("%02d", 16:31))
  
  gdd_short_df <- gdd_short_df %>% 
    mutate(
      !!paste0("gdd_", year, "_05_a") := rowSums(select(., month05_a)),
      !!paste0("gdd_", year, "_05_b") := rowSums(select(., month05_b)),
      !!paste0("gdd_", year, "_06_a") := rowSums(select(., month06_a)),
      !!paste0("gdd_", year, "_06_b") := rowSums(select(., month06_b)),
      !!paste0("gdd_", year, "_07_a") := rowSums(select(., month07_a)),
      !!paste0("gdd_", year, "_07_b") := rowSums(select(., month07_b)),
      !!paste0("gdd_", year, "_08_a") := rowSums(select(., month08_a)),
      !!paste0("gdd_", year, "_08_b") := rowSums(select(., month08_b)),
    )
}

# Keep only the aggregated columns 
gdd_short_df <- gdd_short_df %>%
  select(geometry, unique_id, ends_with("_a"), ends_with("_b"))#, ends_with("_c")) #we don't keep unique_id here to avoid duplication later on

# Aggregate MeanTemp Variables over 2 intervals per month: 
meantemp_short <- meantemp # generate new dataset for the aggregated variables 
# Drop geometry for aggregation (geometry will later be added again)
meantemp_short_df <- meantemp_short %>% 
  st_drop_geometry() %>%
  mutate(geometry = st_geometry(meantemp_short))
# Aggregate columns
for (year in 2000:2020) {
  year_str <- sprintf("%02d", year - 2000)
  month05_a <- paste0(year, "-05-", sprintf("%02d", 1:15))
  month05_b <- paste0(year, "-05-", sprintf("%02d", 16:31))
  month06_a <- paste0(year, "-06-", sprintf("%02d", 1:15))
  month06_b <- paste0(year, "-06-", sprintf("%02d", 16:30))
  month07_a <- paste0(year, "-07-", sprintf("%02d", 1:15))
  month07_b <- paste0(year, "-07-", sprintf("%02d", 16:31))
  month08_a <- paste0(year, "-08-", sprintf("%02d", 1:15))
  month08_b <- paste0(year, "-08-", sprintf("%02d", 16:31))
  
  meantemp_short_df <- meantemp_short_df %>% 
    mutate(
      !!paste0("mean_", year, "_05_a") := rowMeans(select(., month05_a)),
      !!paste0("mean_", year, "_05_b") := rowMeans(select(., month05_b)),
      !!paste0("mean_", year, "_06_a") := rowMeans(select(., month06_a)),
      !!paste0("mean_", year, "_06_b") := rowMeans(select(., month06_b)),
      !!paste0("mean_", year, "_07_a") := rowMeans(select(., month07_a)),
      !!paste0("mean_", year, "_07_b") := rowMeans(select(., month07_b)),
      !!paste0("mean_", year, "_08_a") := rowMeans(select(., month08_a)),
      !!paste0("mean_", year, "_08_b") := rowMeans(select(., month08_b)),
    )
}

# Keep only the aggregated columns 
meantemp_short_df <- meantemp_short_df %>%
  select(unique_id, ends_with("_a"), ends_with("_b"))#, ends_with("_c"))

# Merge both datasets based on unique_id:
temp_oaf <- gdd_short_df %>% left_join(meantemp_short_df, by="unique_id")

# Add geometry column again
temp_oaf <- st_as_sf(temp_oaf, sf_column_name = "geometry")



# 3) Prep Soil data ----
# Previous data prep section doesn't apply to soil as the datasets do not require any merging

soil_oaf <- soil_oaf %>% 
  group_by(unique_id) %>%
  slice(1) %>%
  ungroup() %>%
  select(-matches("system.index|classification|constant")) %>%
  mutate(geometry = geojsonsf::geojson_sf(.geo)$geometry) %>%
  select(-.geo) %>%
  sf::st_as_sf()

soil_50k <- soil_50k %>% 
  group_by(.geo) %>%
  slice(1) %>%
  ungroup() %>%
  select(-matches("system.index|classification|constant|unique")) %>%
  mutate(geometry = geojsonsf::geojson_sf(.geo)$geometry) %>%
  select(-.geo) %>%
  sf::st_as_sf()

# 4) Prep CHIRPS data ----
# Dataset every 5 days, so much higher temporal resolution than NDVI/EVI/aggregated temp data 
# Aggregation to two variables per month 

# Check for missings in both datasets
anyNA(chirp_oaf) # yes 
# Drop observations without data except for location and unique_id: 
chirp_oaf <- chirp_oaf %>%
  filter(!if_all(-c(geometry, unique_id), is.na))

# Check again 
anyNA(chirp_oaf) # no NAs left
anyNA(chirp_50k)


# 4.1) oaf crop cut locations ---- 
# Drop geometry for merging first: #TODO: EB notes "Columns `2016-05-01`, `2016-05-06`, and `2016-05-11` don't exist."
chirp_oaf_sum <- chirp_oaf %>% 
  st_drop_geometry() %>%
  mutate(geometry = st_geometry(chirp_oaf))  
# Aggregate over intervals: 
for (year in 2000:2020) {
  month05_a <- c(paste0(year, "-05-01"), paste0(year, "-05-06"), paste0(year, "-05-11"))
  month05_b <- c(paste0(year, "-05-16"), paste0(year, "-05-21"), paste0(year, "-05-26"))
  month06_a <- c(paste0(year, "-06-01"), paste0(year, "-06-06"), paste0(year, "-06-11"))
  month06_b <- c(paste0(year, "-06-16"), paste0(year, "-06-21"), paste0(year, "-06-26"))
  month07_a <- c(paste0(year, "-07-01"), paste0(year, "-07-06"), paste0(year, "-07-11"))
  month07_b <- c(paste0(year, "-07-16"), paste0(year, "-07-21"), paste0(year, "-07-26"))
  month08_a <- c(paste0(year, "-08-01"), paste0(year, "-08-06"), paste0(year, "-08-11"))
  month08_b <- c(paste0(year, "-08-16"), paste0(year, "-08-21"), paste0(year, "-08-26"))
  
  chirp_oaf_sum <- chirp_oaf_sum %>% 
    mutate(
      !!paste0(year, "_05_a") := rowSums(select(., month05_a)),
      !!paste0(year, "_05_b") := rowSums(select(., month05_b)),
      !!paste0(year, "_06_a") := rowSums(select(., month06_a)),
      !!paste0(year, "_06_b") := rowSums(select(., month06_b)),
      !!paste0(year, "_07_a") := rowSums(select(., month07_a)),
      !!paste0(year, "_07_b") := rowSums(select(., month07_b)),
      !!paste0(year, "_08_a") := rowSums(select(., month08_a)),
      !!paste0(year, "_08_b") := rowSums(select(., month08_b)),
    )
}
# Keep only the aggregated columns 
chirp_oaf_sum <- chirp_oaf_sum %>%
  select(geometry, unique_id, ends_with("_a"), ends_with("_b"))#, ends_with("_c"))

# Add geometry column again
chirp_oaf <- st_as_sf(chirp_oaf_sum, sf_column_name = "geometry")

# 4.2) 50k geopoints ----
# Drop geometry for merging first: 
chirp_50k_sum <- chirp_50k %>% 
  st_drop_geometry() %>%
  mutate(geometry = st_geometry(chirp_50k))  
# Aggregate over intervals: 
for (year in 2000:2020) {# for (year in 2000:2020) {
  month05_a <- c(paste0(year, "-05-01"), paste0(year, "-05-06"), paste0(year, "-05-11"))
  month05_b <- c(paste0(year, "-05-16"), paste0(year, "-05-21"), paste0(year, "-05-26"))
  month06_a <- c(paste0(year, "-06-01"), paste0(year, "-06-06"), paste0(year, "-06-11"))
  month06_b <- c(paste0(year, "-06-16"), paste0(year, "-06-21"), paste0(year, "-06-26"))
  month07_a <- c(paste0(year, "-07-01"), paste0(year, "-07-06"), paste0(year, "-07-11"))
  month07_b <- c(paste0(year, "-07-16"), paste0(year, "-07-21"), paste0(year, "-07-26"))
  month08_a <- c(paste0(year, "-08-01"), paste0(year, "-08-06"), paste0(year, "-08-11"))
  month08_b <- c(paste0(year, "-08-16"), paste0(year, "-08-21"), paste0(year, "-08-26"))
  
  chirp_50k_sum <- chirp_50k_sum %>% 
    mutate(
      !!paste0(year, "_05_a") := rowSums(select(., month05_a)),
      !!paste0(year, "_05_b") := rowSums(select(., month05_b)),
      !!paste0(year, "_06_a") := rowSums(select(., month06_a)),
      !!paste0(year, "_06_b") := rowSums(select(., month06_b)),
      !!paste0(year, "_07_a") := rowSums(select(., month07_a)),
      !!paste0(year, "_07_b") := rowSums(select(., month07_b)),
      !!paste0(year, "_08_a") := rowSums(select(., month08_a)),
      !!paste0(year, "_08_b") := rowSums(select(., month08_b)),
    )
}
# Keep only the aggregated columns 
chirp_50k_sum <- chirp_50k_sum %>%
  select(geometry, ends_with("_a"), ends_with("_b"))

# Add geometry column again
chirp_50k <- st_as_sf(chirp_50k_sum, sf_column_name = "geometry")


# 5) Filter oaf point based datasets for ROI ----
# Some crop cuts were conducted outside our ROI and are still included in dataset
# Filtering done based on clipped admin areas

# Function to filter oaf point datasets based on admin1_clipped boundary
filter_by_boundary <- function(sf_object, boundary) {
  if (!inherits(sf_object, "sf") || !inherits(boundary, "sf")) {
    stop("Both inputs must be sf objects")
  }
  # Apply spatial filter
  filtered <- sf::st_filter(sf_object, boundary)
  
  # Remove observation with unique_id == "c6dc4134" (justification is in script 04_calcDesc, as this geolocation is on water)
  if ("unique_id" %in% names(filtered)) {
    filtered <- filtered[filtered$unique_id != "c6dc4134", ]
  }
  
  return(filtered)
}

# Save unfiltered data (in case there is need to revert back to them)
temp_oaf_unfiltered <- temp_oaf
ndvi_oaf_unfiltered <- ndvi_oaf
evi_oaf_unfiltered <- evi_oaf
soil_oaf_unfiltered <- soil_oaf
chirp_oaf_unfiltered <- chirp_oaf

# Apply function to all oaf point datasets
temp_oaf <- filter_by_boundary(temp_oaf, admin1_clipped)
ndvi_oaf <- filter_by_boundary(ndvi_oaf, admin1_clipped)
evi_oaf <- filter_by_boundary(evi_oaf, admin1_clipped)
soil_oaf <- filter_by_boundary(soil_oaf, admin1_clipped)
chirp_oaf <- filter_by_boundary(chirp_oaf, admin1_clipped)

# Visual check
# All points should be within the clipped area 
ggplot() +
  geom_sf(data = temp_oaf,color = "green", size = 0.2) +
  geom_sf(data = ndvi_oaf,color = "blue", size = 0.2) +
  geom_sf(data = evi_oaf,color = "red", size = 0.2) +
  geom_sf(data = soil_oaf,color = "yellow", size = 0.2) +
  geom_sf(data = chirp_oaf,color = "orange", size = 0.2) +
  geom_sf(data = admin1_clipped, fill = NA, color = "black") +
  theme_minimal() +
  labs(title = "Points within ROI?") # EB looks good 6/18/25

# 6) Deal with NAs----
# 6.1) Check first ----
# For Temp and CHIRPS we tested for missingness prior to aggregating the variables 
# After the aggregation, we also do not expect any missing data
anyNA(temp_oaf)
anyNA(chirp_oaf)
anyNA(temp_50k)
anyNA(chirp_50k)

# For all other indicators there are NAs in the data 
anyNA(ndvi_oaf)  #yes
anyNA(evi_oaf)   #yes
anyNA(soil_oaf)  #yes
anyNA(ndvi_50k)           #yes
anyNA(evi_50k)            #yes
anyNA(soil_50k)           #yes

# Copy of dataset with missings: 
ndvi_oaf_withNA <- ndvi_oaf # keep copy in case we need to revert back to it, but work with ndvi_c dataset
ndvi_50k_withNA <- ndvi_50k
evi_oaf_withNA <- evi_oaf
evi_50k_withNA <- evi_50k
soil_oaf_withNA <- soil_oaf
soil_50k_withNA <- soil_50k

# For all of them, drop obs for which all variables (except geometry and unique_id, if included) are NA
ndvi_oaf <- ndvi_oaf %>% filter(!if_all(-c(geometry, unique_id), is.na))
evi_oaf <- evi_oaf %>% filter(!if_all(-c(geometry, unique_id), is.na))
soil_oaf <- soil_oaf %>% filter(!if_all(-c(geometry, unique_id), is.na))
ndvi_50k <- ndvi_50k[!apply(st_drop_geometry(ndvi_50k), 1, function(row) all(is.na(row))), ]
evi_50k <- evi_50k[!apply(st_drop_geometry(evi_50k), 1, function(row) all(is.na(row))), ]
soil_50k <- soil_50k[!apply(st_drop_geometry(soil_50k), 1, function(row) all(is.na(row))), ]

# Check again for NAs: 
anyNA(ndvi_oaf)  #yes
anyNA(evi_oaf)   #yes
anyNA(soil_oaf)  #no
anyNA(ndvi_50k)           #yes
anyNA(evi_50k)            #yes
anyNA(soil_50k)           #yes

# Check where the missings are
colSums(is.na(ndvi_oaf))
colSums(is.na(evi_oaf))
colSums(is.na(ndvi_50k)) 
colSums(is.na(evi_50k)) 
colSums(is.na(soil_50k)) # only one obs in one variable
# For NDVI and EVI: Peak of NAs end of August 2020!

# 6.2) Replace timeseries NAs ----
# Replace missing values with average on date before and after:
replace_na_withavg <- function(df) {
  
  # Identify columns that are dates
  date_cols <- names(df)[!is.na(as.Date(names(df), format = "%Y-%m-%d"))]
  # Make sure they are sorted chronologically
  date_cols <- date_cols[order(as.Date(date_cols))]
  
  # We replace first and last columns with column averages 
  # First and last column are in April and September and will be dropped later
  # This is just done to ensure that all variables have a previous and later value to take the average from and hence, that no NAs remain in final output
  
  # Identify first and last date columns
  first_col <- date_cols[1]
  last_col <- date_cols[length(date_cols)]
  # Replace NAs in first column with column mean
  df[[first_col]][is.na(df[[first_col]])] <- mean(df[[first_col]], na.rm = TRUE)
  # Replace NAs in last column with column mean
  df[[last_col]][is.na(df[[last_col]])] <- mean(df[[last_col]], na.rm = TRUE)
  
  # Replace missing values with average on date before and after
  for (i in seq_along(date_cols)) {
    
    # Get the current column (ndvi value)
    col <- df[[date_cols[i]]]
    
    # Get previous and next columns 
    # prev_2 <- if (i > 2) df[[date_cols[i - 2]]] else NA
    prev_1 <- if (i > 1) df[[date_cols[i - 1]]] else NA
    next_1 <- if (i < length(date_cols)) df[[date_cols[i + 1]]] else NA
    # next_2 <- if (i < length(date_cols) - 1) df[[date_cols[i + 2]]] else NA

    # Impute NA values using the average of previous and next columns, row by row
    for (j in 1:nrow(df)) {
      if (is.na(col[j])) {
        # Use the average of available values
        # avg <- mean(c(prev_2[j], prev_1[j], next_1[j], next_2[j]), na.rm = TRUE)
        avg <- mean(c(prev_1[j], next_1[j]), na.rm = TRUE)
        col[j] <- avg
      }
    }
    # Assign the updated column back to the dataframe
    df[[date_cols[i]]] <- col
  }
  # Filter out the months of April and September as they are no longer needed after the imputation
  df <- df %>%
    select(-matches("^\\d{4}-(04|09)-\\d{2}$")) 
  
  return(df)
}

# Apply function
ndvi_oaf <- replace_na_withavg(ndvi_oaf)
evi_oaf <- replace_na_withavg(evi_oaf)
ndvi_50k <- replace_na_withavg(ndvi_50k)
evi_50k <- replace_na_withavg(evi_50k)

anyNA(ndvi_oaf)  
anyNA(evi_oaf)   
anyNA(ndvi_50k)          
anyNA(evi_50k)            


# 6.3) Replace NAs in soil data ---- 
# replace missings with column average
soil_50k[] <- lapply(soil_50k, function(col) {
  if (is.numeric(col)) {
    col[is.na(col)] <- mean(col, na.rm = TRUE)
  }
  col
})
anyNA(soil_50k)


# 7) Generate "all" indicators dataframe ----
# 7.1) oaf crop cut database ----
all_oaf <- temp_oaf %>%
  as.data.frame() %>%
  rename_with(~paste0("temp_", .), -c(geometry, unique_id)) %>%
  
  left_join(evi_oaf %>%
      as.data.frame() %>%
      rename_with(~paste0("evi_", .), -unique_id),
      by = "unique_id") %>%
  
  left_join(ndvi_oaf %>%
      as.data.frame() %>%
      rename_with(~paste0("ndvi_", .), -unique_id),
      by = "unique_id") %>%
  
  left_join(chirp_oaf %>%
      as.data.frame() %>%
      rename_with(~paste0("chirp_", .), -unique_id),
      by = "unique_id") %>%
  
  left_join(soil_oaf %>%
      as.data.frame() %>%
      rename_with(~paste0("soil_", .), -unique_id),
      by = "unique_id") %>%
  select(-matches("_geometry")) # This drops all the additional "geometry"-variables but keeps the one without prefix

# Convert back to sf
all_oaf <- st_as_sf(all_oaf, sf_column_name = "geometry")

# 7.2) 50k geopoints ----
all_50k <- temp_50k %>%
  as.data.frame() %>%
  rename_with(~paste0("temp_", .), -geometry) %>%
  
  left_join(evi_50k %>%
              as.data.frame() %>%
              rename_with(~paste0("evi_", .), -geometry),
            by = "geometry") %>%
  
  left_join(ndvi_50k %>%
              as.data.frame() %>%
              rename_with(~paste0("ndvi_", .), -geometry),
            by = "geometry") %>%
  
  left_join(chirp_50k %>%
              as.data.frame() %>%
              rename_with(~paste0("chirp_", .), -geometry),
            by = "geometry") %>%
  
  left_join(soil_50k %>%
              as.data.frame() %>%
              rename_with(~paste0("soil_", .), -geometry),
            by = "geometry")

# Convert back to sf
all_50k <- st_as_sf(all_50k, sf_column_name = "geometry")


# 9) Append datasets based on different point bases  ----
# The appended datasets are used for clustering. We do not need to have observations on the crop cut level nor do we need the unique_id.
# Clean oaf data first:
# - drop the unique_id column 
# - remove observations with the same geolocation 
datasets <- c('temp_oaf', 'chirp_oaf', 'ndvi_oaf', 'evi_oaf','soil_oaf', 'all_oaf')
for (df in datasets) {
df_clean <- get(df) %>%
  filter(unique_id!= "c6dc4134")%>% # this is the result of the outlier testing in lines 445... of script 04_calcDescriptives (where we test for outliers)
  select(-unique_id) %>%
  group_by(geometry)%>%
  slice(1) %>%
  ungroup()
new_name <- paste0(df, "_clean")
assign(new_name, df_clean)
}

# Ensure the datasets are on the same CRS
identical(st_crs(soil_oaf_clean), st_crs(soil_50k))
identical(st_crs(temp_oaf_clean), st_crs(temp_50k))
identical(st_crs(chirp_oaf_clean), st_crs(chirp_50k))
identical(st_crs(ndvi_oaf_clean), st_crs(ndvi_50k))
identical(st_crs(evi_oaf_clean), st_crs(evi_50k))
identical(st_crs(all_oaf_clean), st_crs(all_50k))

# Append the "_oaf" dataframes to the other dataframes
temp_2 <- rbind(temp_50k, temp_oaf_clean)
ndvi_2 <- rbind(ndvi_50k, ndvi_oaf_clean)
evi_2 <- rbind(evi_50k, evi_oaf_clean)
soil_2 <- rbind(soil_50k, soil_oaf_clean)
chirp_2 <- rbind(chirp_50k, chirp_oaf_clean)
all_2 <- rbind(all_50k, all_oaf_clean)

#10) Adapt datasets for sumstats ----
# we need meantemp and gdd as dataframes to create tab2, but as of now they're 
# too large to be pushed to git. Hence, we split the years. 
meantemp_raw <- meantemp %>% 
  select(unique_id, matches("^(2000|2001|2002|2003|2004|2005|2006|2007|2008|2009|2010|2011|2012|2013|2014|2015)-(05|06|07|08)-\\d{2}"))
gdd_raw <- gdd %>% 
  select(unique_id, matches("^(2000|2001|2002|2003|2004|2005|2006|2007|2008|2009|2010|2011|2012|2013|2014|2015)-(05|06|07|08)-\\d{2}"))

# Filter for OAF crop cut locations: 

# Therefore, define new dataset just to check locations
geo_oaf <- oaf %>% 
  select(unique_id) %>%
  st_drop_geometry() %>%
  filter(unique_id!= "c6dc4134")%>% # this is the result of the outlier testing in lines 445... of script 04_calcDescriptives (where we test for outliers)
  mutate(exist = 1)

# And then clean the datasets:
meantemp_clean <- meantemp_raw %>%
  group_by(geometry) %>%
  slice(1) %>%
  ungroup() %>%
  st_drop_geometry() %>%
  left_join(geo_oaf, by = "unique_id") %>%
  filter(!is.na(exist)) # to keep only those obs that were also in geo_oaf

gdd_clean <- gdd_raw %>%
  group_by(geometry) %>%
  slice(1) %>%
  ungroup() %>%
  st_drop_geometry() %>%
  left_join(geo_oaf, by = "unique_id") %>%
  filter(!is.na(exist)) # to keep only those obs that were also in geo_oaf

# 11) Save datasets ----

# Resulting datasets:

# Spatial datasets based on oaf crop cut locations:
# Years: 2000-2020
# Suited for 'oracle' scenario
save(ndvi_oaf, file = paste0(output_data, "/ndvi_oaf.RData"))
save(evi_oaf, file = paste0(output_data, "/evi_oaf.RData"))
save(temp_oaf, file = paste0(output_data, "/temp_oaf.RData"))
save(soil_oaf, file = paste0(output_data, "/soil_oaf.RData"))
save(chirp_oaf, file = paste0(output_data, "/chirp_oaf.RData"))
save(all_oaf, file = paste0(output_data, "/all_oaf.RData"))
# Needed for Sumstats on spatial data:
save(meantemp_clean, file = paste0(output_data, "/meantemp.RData"))
save(gdd_clean, file = paste0(output_data, "/gdd.RData"))

# Spatial datasets based on 50k geopoints:
# Years: 2000-2020
# Primarily used as an input for the joined datasets below
save(ndvi_50k, file = paste0(output_data, "/ndvi_50k.RData"))
save(evi_50k, file = paste0(output_data, "/evi_50k.RData"))
save(temp_50k, file = paste0(output_data, "/temp_50k.RData"))
save(soil_50k, file = paste0(output_data, "/soil_50k.RData"))
save(chirp_50k, file = paste0(output_data, "/chirp_50k.RData"))
save(all_50k, file = paste0(output_data, "/all_50k.RData"))

# Spatial datasets based on both (oaf (minus location duplicates) + 50k):
# Years: 2000-2020
# Needed for Clustering Procedure
save(ndvi_2, file = paste0(output_data, "/ndvi_2.RData"))
save(evi_2, file = paste0(output_data, "/evi_2.RData"))
save(temp_2, file = paste0(output_data, "/temp_2.RData"))
save(soil_2, file = paste0(output_data, "/soil_2.RData"))
save(chirp_2, file = paste0(output_data, "/chirp_2.RData"))
save(all_2, file = paste0(output_data, "/all_2.RData"))

# OAF datasets:
# OAF data (merged with admin units, filtered for ROI, duplicates dropped, zero yields filtered out, not winsorized, no clusters yet)
# 1) raw imported OAF data:
save(oaf_raw, file = paste0(output_data, "/oaf_raw.RData"))
# 2) OAF data (merged with admin units, filtered for ROI, duplicates dropped, zero yields filtered out, not winsorized, no clusters yet)
save(oaf, file = paste0(output_data, "/oaf_noclusters.RData"))

# Admin datasets:
save(admin3_clipped, file = paste0(output_data, "/admin3_clipped.RData"))
save(admin1_clipped, file = paste0(output_data, "/admin1_clipped.RData"))
save(admin3_kenya, file = paste0(output_data, "/admin3_kenya.RData"))
save(admin1_kenya, file = paste0(output_data, "/admin1_kenya.RData"))
