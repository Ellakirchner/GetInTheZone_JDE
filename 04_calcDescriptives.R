#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Author: Ella Kirchner, Elinor Benami, Andrew Hobbs 
# Project name: RS Get in the Zone
# Date Last Updated: # Fri Sep 12 11:03:04 2025
# Purpose:  Generate descriptive statistics on input data (crop cut and spatial data) as well as various specific numbers mentioned throughout the manuscript and figures for the appendix
# Input:  all to be loaded from the RProj
          # - 6 .RData files with spatial data based on oaf crop cut locations (ndvi, evi, temp, chirp, soil, all with suffix "_oaf")
          # - 2 .RData files with temperature data based on oaf crop cut locations, meantemp and gdd separately
          # - 4 .RData files on admin data (level 1 and level 3, each for full country and clipped to ROI)
          # - (processed) OAF dataset with cluster assignment of k=14
          # - (processed) OAF dataset with only chirps cluster assignment based only on OAF points for clustering
          # - 6 .RData files on insurance zone assignments on ward level (one .RData for every agri-environmental indicator, prefix "dynamic_")
# Output Files: figures and tables; no datasets 
          # - Table 1 (Descriptives on Crop cut data)
          # - Table 2 (Descriptives on spatial data)
          # - Table A.2 (Exemplary distribution of observations by cluster)
          # - Fig. A.1 (Yield Distributions)
          # - Fig. A.5 (Illustration of different zoning approaches)
          # - Fig. A.6 (Number of Insurance Zones considered vs. designed)
# ReadMe: 1) Adapt file path to work in R-project, then loading of datasets should work
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# %%%%%%%%%%%%%%%%%%%%%%%%%%%
# Set Up 
# %%%%%%%%%%%%%%%%%%%%%%%%%%%

# 1. Clean out workspace
rm(list = ls(all = TRUE))

# 2. Set Up File path & directory: 
file_path <- "put_your_file_path_here" #File path to project folder
setwd(file_path)

saving <- paste0(file_path,"/figs")

# 3. Set Seed for replication 
set.seed(123456789)

# 4. Set date 
date_today <- as.character(Sys.Date())


# %%%%%%%%%%%%%%%%%%%%%%%%%%%
# Load libraries 
# %%%%%%%%%%%%%%%%%%%%%%%%%%%
library(dplyr)
library(tidyr)
library(xtable)
library(sf)  
library(tidyverse)
library(RColorBrewer)
library(cowplot)

# %%%%%%%%%%%%%%%%%%%%%%%%%%%
# Load datasets ---- 
# %%%%%%%%%%%%%%%%%%%%%%%%%%%
# Admin data
load("datasets/admin3_clipped.RData")
load("datasets/admin3_kenya.RData")
load("datasets/admin1_clipped.RData")
load("datasets/admin1_kenya.RData")

# Anonymized OAF data for replication:
load("datasets/oaf_anonym.RData")
oaf <- oaf_anonym

# # Spatial data
# load("datasets/ndvi_oaf.RData")
# load("datasets/evi_oaf.RData")
# load("datasets/chirp_oaf.RData")
# load("datasets/temp_oaf.RData")
# load("datasets/soil_oaf.RData")
# load("datasets/all_oaf.RData")
# load("datasets/meantemp.RData")
# load("datasets/gdd.RData")

# Datasets required for section 6 of this script 
load("datasets/dynamic_chirps.RData")
load("datasets/dynamic_temp.RData")
load("datasets/dynamic_ndvi.RData")
load("datasets/dynamic_evi.RData")
load("datasets/dynamic_soil.RData")
load("datasets/dynamic_all.RData")
load("datasets/oaf_woafdataclusters_anonym.RData") # Needed for panel c in Fig. A.5
oaf_woafdataclusters <- oaf_woafdataclusters_anonym

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Calc Descriptives
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%
oaf_allobs <- oaf # unfiltered dataset needed for table 2 and for appendix (fig a.5)

# Key figures on admin units ----
# see script 01_prepData 

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 1 Descriptives on Crop cut data ----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Script to generate Tab.1 of the manuscript. 

# Winsorization 
# Winsorization is also done in analysis, hence also for this desc stats table
# (roughly 8t is the maximum realistic yield)
quantiles <- quantile(oaf$yield, c(0.01, 0.05, 0.1, 0.25, 0.5,0.75,0.9,0.95, 0.99))
quantiles # 99% equals roughly 8t, hence, we use 99% 
# Winsorize symmetrically on both ends (1%, 99%):
oaf$yield <- pmin(oaf$yield, quantile(oaf$yield, 0.99))
oaf$yield <- pmax(oaf$yield, quantile(oaf$yield, 0.01))
summary(oaf$yield)

# Drop those wards that do not have data in every year: 
# Filter to keep only name_concat values that appear in all years
all_years <- unique(oaf$year)

oaf <- oaf %>%
  group_by(name_concat) %>%
  filter(length(unique(year)) == length(all_years)) %>% # drops 7k observations 
  ungroup()

# Generate table
tab1<- oaf %>%
  st_drop_geometry()%>%
  filter(year %in% 2016:2020) %>%
  group_by(year) %>%
  summarize(
    `Total on-field crop cuts` = n(),
    `Average yield (t/ha)` = mean(yield, na.rm = TRUE),
    `SD` = sd(yield, na.rm = TRUE),
    .groups = 'drop'
  ) %>%
  pivot_longer(-year, names_to = "Variable", values_to = "value") %>%
  pivot_wider(names_from = year, values_from = value, names_prefix = "") %>%
  select(Variable, `2016`, `2017`, `2018`, `2019`, `2020`)

#Table output in latex
tab1_output <- xtable(tab1)
print(tab1_output, include.rownames = FALSE)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 2 Yield distributions ----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Script to generate figure A.1 presented in the appendix to the manuscript.

clean_chart_clutter_explore <- 
  theme(    
    panel.grid.major = element_blank(),      # Remove panel grid lines
    panel.grid.minor = element_blank(),      # Remove panel grid lines
    panel.background = element_blank(),      # Remove panel background
    axis.line = element_line(colour = "grey"),       # Add axis line
    axis.title.y = element_text(angle = 0, vjust = 0.5),      # Rotate y axis so don't have to crank head
    legend.position="bottom"
  ) 

# Calculate year means to add them to graph
year_means <- oaf %>%
  group_by(year) %>%
  summarise(mean_yield = mean(yield, na.rm = TRUE), .groups = 'drop')

yield_lines <- ggplot(data = oaf, aes(x = yield, color = factor(year))) +
  geom_density(size = 1.2) +
  geom_vline(data = year_means, 
             aes(xintercept = mean_yield, color = factor(year)),
             linetype = "dashed", size = 1) +
  scale_color_manual(values = c("2016" = "#2B3B90", "2017" = "gold2", "2018" = "#82A67D", "2019" = "#8E14AD", "2020" = "#c47500")) +
  labs(title = "Distribution of Yields by Year",
       subtitle = "Yields winsorized to 1st and 99th percentile",
       x = "Yield (t/ha)",
       y = "Density",
       color = "Year") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  clean_chart_clutter_explore+
  theme(plot.title = element_text(size = 13, face = "bold"),
        plot.subtitle = element_text(size = 12, face = "plain"))
yield_lines 

# Save plot 
ggsave(paste0(saving,"/yield_distribution.png"), plot= yield_lines, width=8, height= 4)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 3 Descriptives on plots and income ----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# This section calculates specific descriptive figures mentioned throughout the manuscript.

summary(oaf$plot_acres)
summary(oaf$plot_hectares)
oaf$plot_hectares <- oaf$plot_acres * 0.405 
summary(oaf$plot_hectares)

# Income Parameters
price_per_kg = 40
kes_per_usd =  105
price_kg_usd = price_per_kg/kes_per_usd
price_t_usd=price_kg_usd*1000

oaf$income <- oaf$yield * price_t_usd
oaf$income_plot <- oaf$plot_hectares * oaf$income

summary(oaf$income) # Average income for instance mentioned in footnote 12. 
summary(oaf$income_plot)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 4 Descriptives on area size ----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# This section calculates the figures mentioned in footnote 20 (appendix A.4.1). 

# Check current projection system 
st_crs(admin1_clipped)

# Reproject to UTM Zone 37S for Kenya
admin1_area <- st_transform(admin1_clipped, crs = "EPSG:32637")

# Calculate areas in square kilometers
admin1_area$area_sqkm <- as.numeric(st_area(admin1_area)) / 1e6

# Results (total and avg area)
total_area_sqkm <- sum(admin1_area$area_sqkm)
average_area_sqkm <- mean(admin1_area$area_sqkm)

cat("Total area:", total_area_sqkm, "sq km\n")
cat("Average area:", average_area_sqkm, "sq km\n")

# Area of wards with data: 

# Filter oaf data and reduce to wards of interest 
all_years <- unique(oaf$year)
oaf_filtered <- oaf %>%
  group_by(name_concat) %>%
  filter(length(unique(year)) == length(all_years)) %>% # drops 7k observations
  ungroup()
wardstokeep <- oaf_filtered%>% 
  select(name_concat) %>% 
  st_drop_geometry() %>% 
  group_by(name_concat) %>% 
  slice(1)%>%
  ungroup()

admin3_area <- wardstokeep %>% 
  left_join(admin3_clipped, by = "name_concat")%>% 
  st_as_sf()
st_crs(admin3_clipped)

# Reproject to UTM Zone 37S for Kenya
admin3_area <- st_transform(admin3_area, crs = "EPSG:32637")

# Calculate areas in square kilometers
admin3_area$area_sqkm <- as.numeric(st_area(admin3_area)) / 1e6

# Results (total and avg area)
total_area_sqkm_wardswithdata <- sum(admin3_area$area_sqkm)
avg_k14_wards <- print(total_area_sqkm_wardswithdata/14)

cat("Total area:", total_area_sqkm_wardswithdata, "sq km\n")
cat("Average area when grouping only wards:", avg_k14_wards, "sq km\n")

# #%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# # 5 Desc Stats on Spatial data ---- 
# #%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# ### NOT possible on anonymized data. 
# 
# # This section computes the statistics presented in Table 2 of the manuscript. 
# 
# # Make sure that only those observation also considered in the analysis are included in these descriptive data: 
# geo_oaf <- oaf_allobs %>% 
#   select(unique_id) %>%
#   st_drop_geometry() %>%
#   filter(unique_id!= "c6dc4134")%>% # based on outlier testing below (see lines 413-431)
#   mutate(exist = 1)
# 
# ### Temperature: ----
# # The temperatures we are using for the clustering are cumulative across 10 day intervals, hence, need to use base data. 
# # The datasets are called meantemp and gdd
# 
# # Calculate table values for mean temp
# # First calculate seasonal temps per year
# seasonal_temps <- meantemp_clean %>%
#   select(unique_id, matches("^(2000|2001|2002|2003|2004|2005|2006|2007|2008|2009|2010|2011|2012|2013|2014|2015)-(05|06|07|08)-\\d{2}")) %>%
#   pivot_longer(-unique_id, names_to = "date", values_to = "temp") %>%
#   mutate(year = as.numeric(str_extract(date, "^\\d{4}"))) %>%
#   group_by(unique_id, year) %>%
#   summarise(seasonal_temp = mean(temp, na.rm = TRUE), .groups = 'drop')
# 
# # Create wide format with year columns + summary stats
# meantemp_tab <- seasonal_temps %>%
#   pivot_wider(names_from = year, values_from = seasonal_temp, names_prefix = "temp_") %>%
#   rowwise() %>%
#   mutate(
#     meantemp_mean = mean(c_across(starts_with("temp_")), na.rm = TRUE),
#     meantemp_min = min(c_across(starts_with("temp_")), na.rm = TRUE),
#     meantemp_max = max(c_across(starts_with("temp_")), na.rm = TRUE),
#     meantemp_sd = sd(c_across(starts_with("temp_")), na.rm = TRUE)
#   ) %>%
#   ungroup()
# 
# seasonal_gdd <- gdd_clean %>%
#   select(unique_id, matches("^(2000|2001|2002|2003|2004|2005|2006|2007|2008|2009|2010|2011|2012|2013|2014|2015)-(05|06|07|08)-\\d{2}")) %>%
#   pivot_longer(-unique_id, names_to = "date", values_to = "gdd") %>%
#   mutate(year = as.numeric(str_extract(date, "^\\d{4}"))) %>%
#   group_by(unique_id, year) %>%
#   summarise(seasonal_gdd = sum(gdd, na.rm = TRUE), .groups = 'drop')
# 
# # Create wide format with year columns + summary stats
# gdd_tab <- seasonal_gdd %>%
#   pivot_wider(names_from = year, values_from = seasonal_gdd, names_prefix = "gdd_") %>%
#   rowwise() %>%
#   mutate(
#     gdd_mean = mean(c_across(starts_with("gdd_")), na.rm = TRUE),
#     gdd_min = min(c_across(starts_with("gdd_")), na.rm = TRUE),
#     gdd_max = max(c_across(starts_with("gdd_")), na.rm = TRUE),
#     gdd_sd = sd(c_across(starts_with("gdd_")), na.rm = TRUE)
#   ) %>%
#   ungroup()
# 
# ### CHIRPS ----
# # Here we can take the chirps_oaf data as it is a cumulative measure 
# # We already cleaned this dataset for duplicates before (in data prep script), but dropped unique_id, hence clean again and keep unique_id
# chirp_oaf_clean <- chirp_oaf %>%
#   group_by(geometry) %>%
#   slice(1) %>%
#   ungroup() %>%
#   st_drop_geometry() %>%
#   left_join(geo_oaf, by = "unique_id") %>%
#   filter(!is.na(exist)) # to keep only those obs that were also in geo_oaf
# 
# # Same approach as for temp data: 
# seasonal_chirp <- chirp_oaf_clean %>%
#   select(unique_id, matches("^200[0-9]|^201[0-5]")) %>%
#   #  select(unique_id, matches("^(2000|2001|2002|2003|2004|2005|2006|2007|2008|2009|2010|2011|2012|2013|2014|2015)-")) %>%
#   pivot_longer(-unique_id, names_to = "date", values_to = "chirp") %>%
#   mutate(year = as.numeric(str_extract(date, "^\\d{4}"))) %>%
#   group_by(unique_id, year) %>%
#   summarise(seasonal_chirp = sum(chirp, na.rm = TRUE), .groups = 'drop')
# 
# # Create wide format with year columns + summary stats
# chirp_tab <- seasonal_chirp %>%
#   pivot_wider(names_from = year, values_from = seasonal_chirp, names_prefix = "chirp_") %>%
#   rowwise() %>%
#   mutate(
#     chirp_mean = mean(c_across(starts_with("chirp_")), na.rm = TRUE),
#     chirp_min = min(c_across(starts_with("chirp_")), na.rm = TRUE),
#     chirp_max = max(c_across(starts_with("chirp_")), na.rm = TRUE),
#     chirp_sd = sd(c_across(starts_with("chirp_")), na.rm = TRUE)
#   ) %>%
#   ungroup()
# 
# 
# ### NDVI: ----
# # We did not aggregate anything for this variable. 
# # The average greenness is not of primary interest here, but rather the peak NDVI
# 
# # We drop duplicates first (but keep unique_id, not like before)
# ndvi_oaf_clean <- ndvi_oaf %>%
#   group_by(geometry) %>%
#   slice(1) %>%
#   ungroup() %>%
#   st_drop_geometry() %>%
#   left_join(geo_oaf, by = "unique_id") %>%
#   filter(!is.na(exist)) # to keep only those obs that were also in geo_oaf
# 
# # Same approach as before
# seasonal_ndvi <- ndvi_oaf_clean %>%
#   select(unique_id, matches("^(2000|2001|2002|2003|2004|2005|2006|2007|2008|2009|2010|2011|2012|2013|2014|2015)-(05|06|07|08)-\\d{2}")) %>%
#   pivot_longer(-unique_id, names_to = "date", values_to = "ndvi") %>%
#   mutate(year = as.numeric(str_extract(date, "^\\d{4}"))) %>%
#   group_by(unique_id, year) %>%
#   summarise(seasonal_ndvi = max(ndvi, na.rm = TRUE), .groups = 'drop')
# 
# ndvi_tab <- seasonal_ndvi %>%
#   pivot_wider(names_from = year, values_from = seasonal_ndvi, names_prefix = "ndvi_") %>%
#   mutate(across(starts_with("ndvi_"), ~ .x / 10000)) %>%  # Divide all ndvi_ columns by 10000
#   rowwise() %>%
#   mutate(
#     ndvi_mean = mean(c_across(starts_with("ndvi_")), na.rm = TRUE),
#     ndvi_min = min(c_across(starts_with("ndvi_")), na.rm = TRUE),
#     ndvi_max = max(c_across(starts_with("ndvi_")), na.rm = TRUE),
#     ndvi_sd = sd(c_across(starts_with("ndvi_")), na.rm = TRUE)
#   ) %>%
#   ungroup()
# 
# 
# ### EVI: ----
# # We did not aggregate anything for this variable. 
# # The average greeness is not of primary interest here, but rather the peak evi
# 
# # We drop duplicates first (but keep unique_id, not like before)
# evi_oaf_clean <- evi_oaf %>%
#   group_by(geometry) %>%
#   slice(1) %>%
#   ungroup() %>%
#   st_drop_geometry() %>%
#   left_join(geo_oaf, by = "unique_id") %>%
#   filter(!is.na(exist)) # to keep only those obs that were also in geo_oaf
# 
# # Same approach as before
# seasonal_evi <- evi_oaf_clean %>%
#   select(unique_id, matches("^(2000|2001|2002|2003|2004|2005|2006|2007|2008|2009|2010|2011|2012|2013|2014|2015)-(05|06|07|08)-\\d{2}")) %>%
#   pivot_longer(-unique_id, names_to = "date", values_to = "evi") %>%
#   mutate(year = as.numeric(str_extract(date, "^\\d{4}"))) %>%
#   group_by(unique_id, year) %>%
#   summarise(seasonal_evi = max(evi, na.rm = TRUE), .groups = 'drop')
# 
# evi_tab <- seasonal_evi %>%
#   pivot_wider(names_from = year, values_from = seasonal_evi, names_prefix = "evi_") %>%
#   mutate(across(starts_with("evi_"), ~ .x / 10000)) %>%  # Divide all evi_ columns by 10000
#   rowwise() %>%
#   mutate(
#     evi_mean = mean(c_across(starts_with("evi_")), na.rm = TRUE),
#     evi_min = min(c_across(starts_with("evi_")), na.rm = TRUE), 
#     evi_max = max(c_across(starts_with("evi_")), na.rm = TRUE),
#     evi_sd = sd(c_across(starts_with("evi_")), na.rm = TRUE)
#   ) %>%
#   ungroup()
# 
# ### Soil:  ----
# soil_tab <- soil_oaf %>%
#   group_by(geometry) %>%
#   slice(1) %>% 
#   ungroup() %>%
#   st_drop_geometry() %>%
#   rename_with(~ paste0("soil_", .x), .cols = -unique_id) %>%
#   left_join(geo_oaf, by = "unique_id") %>%
#   filter(!is.na(exist)) # to keep only those obs that were also in geo_oaf
# 
# ### Assess negative VI outliers ----
#     # # Sort values according to EVI/NDVI descending: 
#     seasonal_ndvi <- seasonal_ndvi %>%
#       arrange(seasonal_ndvi)
#     head(seasonal_ndvi$seasonal_ndvi, 25)
#     print(head(seasonal_ndvi[,c("unique_id")],25),n=25)
# 
#     seasonal_evi <-seasonal_evi %>%
#       arrange(seasonal_evi)
#     head(seasonal_evi$seasonal_evi, 25)
#     print(head(seasonal_evi[, c("unique_id")], 25), n = 25)
#  
#     # # c6dc4134 and af143dce appear quite often in both dataframes --> double check locations: 
#     # oaf_allobs[oaf_allobs$unique_id %in% c("c6dc4134", "af143dce"), "geometry"]
#     # # --> c6dc4134 is on water! Drop observation from entire analysis. 
#     # # Check if that unique id is part of the filtered oaf data used for the analysis: 
#     # oaf[oaf$unique_id %in% c("c6dc4134"),] # No value does not exist. 
#     # 
#     # # Hence, we only need to exclude that value from the clustering and the descriptives, not from the welfare analysis. 
# 
# ### Create table ----
# # Merge datasets for table: 
# desc <- meantemp_tab %>% 
#   left_join(gdd_tab, by="unique_id")%>%
#   left_join(chirp_tab, by="unique_id")%>%
#   left_join(ndvi_tab, by="unique_id")%>%
#   left_join(evi_tab, by="unique_id")%>%
#   left_join(soil_tab, by="unique_id")%>%
#   select(ends_with("_mean"),
#          ends_with("_min"),
#          ends_with("_max"),
#          ends_with("_sd"),
#          starts_with("soil")
#   )
# 
# # Name all variable prefixes
# var_prefixes <- c("meantemp", "gdd", "chirp", "ndvi", "evi") 
# 
# # Create summary for each variable type
# summary_rows <- map_dfr(var_prefixes, ~ {
#   prefix <- .x
#   tibble(
#     variable = prefix,
#     mean = mean(desc[[paste0(prefix, "_mean")]], na.rm = TRUE),
#     min = min(desc[[paste0(prefix, "_min")]], na.rm = TRUE),
#     max = max(desc[[paste0(prefix, "_max")]], na.rm = TRUE),
#     sd = sd(desc[[paste0(prefix, "_mean")]], na.rm = TRUE)
#   )
# })
# 
# # Add soil variables 
# soil_summary <- desc %>%
#   select(starts_with("soil")) %>%
#   summarise(across(everything(), list(
#     mean = ~ mean(.x, na.rm = TRUE),
#     min = ~ min(.x, na.rm = TRUE),
#     max = ~ max(.x, na.rm = TRUE),
#     sd = ~ sd(.x, na.rm = TRUE)
#   ))) %>%
#   pivot_longer(everything()) %>%
#   separate(name, into = c("variable", "statistic"), sep = "_(?=[^_]*$)") %>%
#   pivot_wider(names_from = statistic, values_from = value)
# 
# # Combine both
# summary_table <- bind_rows(summary_rows, soil_summary)  
# 
# # Adapt variable names: 
# summary_table <- summary_table %>%
#   mutate(
#     variable = case_when(
#       variable == "gdd" ~ "Average sum of growing degree days per season",
#       variable == "meantemp" ~ "Average daily mean temperature (°C) during season",
#       variable == "chirp" ~ "Average cumulative precipitation (mm) per season",
#       variable == "ndvi" ~ "Average peak NDVI value per season",
#       variable == "evi" ~ "Average peak EVI value per season",
#       variable == "soil_al" ~ "Aluminum",
#       variable == "soil_awcp30" ~ "Available water capacity (cm water/cm soil)",
#       variable == "soil_erzd" ~ "Effective root zone depth (cm)",
#       variable == "soil_soc" ~ "Soil organic carbon (dg/kg)",
#       variable == "soil_totN" ~ "Total nitrogen content (cg/kg)", 
#       TRUE ~ variable  # Keep original name if not matched
#     )
#   )
# 
# desc_spatialtable <- xtable(summary_table)
# print(desc_spatialtable, include.rownames = FALSE)
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 6 Appendix on Zone size/number trade-off ----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# This section creates the figures and table in Appendix A.4.1. 

### Gen Fig. A.6 (linegraph) ----
# Define list of wards
all_years <- unique(oaf$year)

wardstokeep <- oaf%>% 
  select(name_concat) %>% 
  st_drop_geometry() %>% 
  group_by(name_concat) %>% 
  slice(1)%>%
  ungroup()
  
count_distinct_values <- function(df, dataset_name) {
  # Filter for wards of interest: 
  data <- wardstokeep %>%
    left_join(df, by="name_concat")
  
  # Find all cluster columns
  cluster_cols <- names(data)[grepl("^cluster_\\d+$", names(data))]
  
  # Sort columns by cluster number
  cluster_nums <- as.numeric(gsub("cluster_", "", cluster_cols))
  cluster_cols <- cluster_cols[order(cluster_nums)]
  
  # Count distinct values for each cluster column
  distinct_counts <- data %>%
    select(all_of(cluster_cols)) %>%
    summarise(across(everything(), ~ n_distinct(.x, na.rm = TRUE))) %>%
    mutate(dataset = dataset_name) %>%
    select(dataset, everything())
  
  return(distinct_counts)
}

# Process each dataset
datasets <- list(
  "dynamic_evi" = dynamic_evi,
  "dynamic_ndvi" = dynamic_ndvi,
  "dynamic_temp" = dynamic_temp, 
  "dynamic_chirps" = dynamic_chirps, 
  "dynamic_soil" = dynamic_soil, 
  "dynamic_all" = dynamic_all
)

# Apply function to each dataset and combine results
zone_counts <- map_dfr(names(datasets), function(name) {
  count_distinct_values(datasets[[name]], name)
})


# Prepare data for visualization
zoneplot_data <- zone_counts %>%
  pivot_longer(cols = starts_with("cluster_"), 
               names_to = "cluster_var", 
               values_to = "distinct_values") %>%
  mutate(cluster_size = as.numeric(gsub("cluster_", "", cluster_var))) %>%
  filter(!is.na(distinct_values)) %>%
  mutate(dataset = case_when(
    dataset == "dynamic_ndvi" ~ "NDVI",
    dataset == "dynamic_evi" ~ "EVI", 
    dataset == "dynamic_temp" ~ "Temperature",
    dataset == "dynamic_chirps" ~ "Precipitation", 
    dataset == "dynamic_soil" ~ "Soil",
    dataset == "dynamic_all" ~ "Combined",
    TRUE ~ dataset  # keeps any other names unchanged
  ))

# Create the main line graph visualization
# Add administrative boundary reference points
admin_points <- data.frame(
  cluster_size = c(14, 100),
  distinct_values = c(14, 67),
  label = c("Admin 1", "Admin 2")
)

zone_plot <- ggplot(zoneplot_data, aes(x = cluster_size, y = distinct_values, color = dataset)) +
  geom_line(linewidth = 1.2) +
  geom_point(data = admin_points, 
             aes(x = cluster_size, y = distinct_values), 
             color = "black", size = 3, shape = 17, inherit.aes = FALSE) +
  annotate("text", x = 14, y = 18, label = "Admin 1", color = "black") +
  annotate("text", x = 100, y = 64, label = "Admin 2", color = "black") +
  scale_color_manual(name = "Agri-environmental\nInput", 
  values = c(
    "Precipitation" = "#2E86C1",
    "NDVI" = "#1E8449",       
    "EVI" = "#82C09A",        
    "Temperature" = "#C3332B",
    "Soil" = "#8B4843", 
    "Combined" = "gold"
  )) +
  labs(
    title = "Development of Number of Insurance Zones incl. in Analysis",
    x = "Data-driven Insurance Zones in Study Region",
    y = "Insurance Zones\nwith Crop Cut Data\nfor Our Analysis",
    color = "Dataset"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust=0.5),
    axis.title.y = element_text(angle = 0, vjust = 0.5, hjust = 0.5),
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    panel.grid.major.x = element_line(color = "grey90", linewidth = 0.5)  # Control major grid lines
  ) +
  scale_x_continuous(
    breaks = c(5, seq(10, 100, by = 10)), 
    limits = c(5, 105),
    expand = c(0, 0))
print(zone_plot)

# Save Fig.A.6
ggsave(paste0(saving, "/appendix_numberofzones_",date_today, ".png"), plot = zone_plot, width = 9, height = 4)

### Gen maps for Fig. A.5 ----
# Define the color palette
cluster_colors <- c(
  "#2B3B90",  # Deep Navy
  "#27AE60",  # Forest Green
  "#8E14AD",  # Burgundy
  "#9d5D7E",  # Slate Blue
  "gold2",  # Teal
  "#c47500",  # Rust Orange
  "#82A67D",  # Sage Green
  "#9b150B",  # Dusty Rose
  "#7D6608",  # Warm Gray
  "#16A085",  # Muted Purple
  "#935116",  # Olive Brown
  "#E8743B",  # Soft Coral
  "#5499C7",  # Steel Blue
  "#6E797F"   # Taupe
)

# # 1 Map illustrating simply the points
# 
# ### NOT possible with anonymized data. 
# 
# # Create dataset with binary indicator whether it is in the final dataset or not
# oaf_filteredmap1 <- oaf%>% 
#   mutate(exist = 1) %>% 
#   select(unique_id, exist) %>% 
#   st_drop_geometry()
# oaf_map1 <- oaf_allobs %>% 
#   left_join(oaf_filteredmap1, by = "unique_id")
# map1 <- ggplot() +
#   geom_sf(data = oaf_map1, aes(color = factor(ifelse(exist == 1, "1", "NA"))), size = 1) +
#   # geom_sf(data= oaf_map1, aes(), color = "black", size = 1) + 
#   geom_sf(data = admin3_clipped, fill = NA, color = "grey35", linewidth = 0.6) +
#   geom_sf(data = admin1_clipped, fill = NA, color = "black", linewidth = 1) +
#   scale_color_manual(
#     values = c("1" = "black", "NA" = "lightblue"),
#     guide = "none"
#   ) +
#   scale_fill_manual(guide = "none") +
#   theme_minimal()+
#   labs(title = "a) Crop Cut Locations in Study Region\n")+
#   theme(
#     plot.title.position = "plot",   # align title to full plot area
#     plot.title = element_text(hjust = 0)  # left-align title itself
#   )
# map1 

# 2 Map & Tab of CHIRPS zones across study region 

# Add geometry of wards to dynamic zone dataset:
dynamic_chirps_wward <- dynamic_chirps %>% 
  left_join(admin3_clipped, by= "name_concat")%>%
  st_sf()

# Make map2: 
map2 <- ggplot()+
  geom_sf(data = dynamic_chirps_wward, aes(fill = as.factor(cluster_14)), color = NA, size = 0) +
  geom_sf(data = admin3_clipped, fill = NA, color = "grey35", linewidth = 0.6) +
  geom_sf(data = admin1_clipped, fill = NA, color = "black", linewidth = 0.9) +
  scale_fill_manual(values = cluster_colors, na.value = "lightgrey", guide = "none")  +
  theme_minimal()+
  labs(title = "b) 14 Insurance Zones (based on Precipitation)\nSpanning Entire Study Region")+
  theme(
    plot.title.position = "plot",   # align title to full plot area
    plot.title = element_text(hjust = 0)  # left-align title itself
  )

map2

# 3 Map & Tab of CHIRPS zones only on oaf data locations 

# Create dataset used for this graphical illustration:
# Extract cluster assignment from saved oaf dataset (ensuring it is only the wards of interest and that the clustering approach was only based on those datapoints)
wards_reduced <- oaf_woafdataclusters %>% 
  group_by(name_concat) %>%
  slice(1) %>%
  ungroup()%>%
  st_drop_geometry()

# Add ward geometry features:
map3_data <-admin3_clipped %>%
  left_join(wards_reduced, by="name_concat")

#Check if correct number of wards has/has no cluster assignment
table(map3_data$c.chirp)
sum(is.na(map3_data$c.chirp))

# Create map
map3 <- ggplot() +
  geom_sf(data = map3_data, aes(fill = as.factor(c.chirp)), color = NA, size = 0) +
  geom_sf(data = admin3_clipped, fill = NA, color = "grey35", linewidth = 0.6) +
  geom_sf(data = admin1_clipped, fill = NA, color = "black", linewidth = 0.9) +
  scale_fill_manual(values = cluster_colors, na.value = "lightgrey", guide = "none")  +
  theme_minimal()+
  labs(title = "c) 14 Insurance Zones (based on Precipitation)\nSpanning only Wards with Crop Cut Locations")+
  theme(
    plot.title.position = "plot",   # align title to full plot area
    plot.title = element_text(hjust = 0)  # left-align title itself
  )
map3

# Legend for the graphs: 
map3_wlegend <- ggplot() +
  geom_sf(data = map3_data, aes(fill = as.factor(c.chirp)), color = NA, size = 0) +
  geom_sf(data = admin3_clipped, fill = NA, color = "grey35", linewidth = 0.6) +
  geom_sf(data = admin1_clipped, fill = NA, color = "black", linewidth = 0.9) +
  scale_fill_manual(name = "Cluster", values = cluster_colors, na.value = "lightgrey")  +
  theme_minimal()
legend_inszone_maps <- get_legend(map3_wlegend)

# # Combine the plots
# 
# ### NOT possible as first graph is missing due to anonymization. 
# 
# spacer <- rectGrob(gp = gpar(col = "white", fill = "white"))
# overview <- grid.arrange(
#   map1, spacer, map2, spacer,map3, spacer, legend_inszone_maps,
#   ncol = 7,
#   nrow = 1,
#   widths = c(1, 0.1, 1, 0.1, 1,0.1,0.2)
#   #heights = c(0.1,1, 0.1) # If space around plot needed
# )
# overview


### Gen Table on obs/zone ----
variables <- c("c.chirp", "c.evi", "c.ndvi", "c.temp", "c.soil", "c.all", "name_1_num")
oaf_freqtable <- oaf %>% 
  mutate(name_1_num = as.numeric(as.factor(name_1)))

freq_table <- data.frame(row.names = 1:14)

for (var in variables) {
  freq_table[[var]] <- table(factor(oaf_freqtable[[var]], levels = 5:18))
}

freq_table <- freq_table %>% 
  mutate(Cluster = 1:14,
         Precipitation = as.integer(c.chirp),
         EVI = as.integer(c.evi), 
         NDVI = as.integer(c.ndvi), 
         Temperature = as.integer(c.temp), 
         Soil = as.integer(c.soil), 
         Combined = as.integer(c.all), 
         Admin1 = as.integer(name_1_num))%>% 
  select(-starts_with("c."), -name_1_num)%>% 
  select(Cluster, everything()) %>% 
  mutate(Cluster = as.character(Cluster))
  

total_row <- freq_table %>%
  select(-Cluster) %>%  # exclude Cluster for summing
  summarise(across(everything(), sum)) %>%
  mutate(Cluster = "Total") %>%
  select(Cluster, everything())

# Bind total row to the original table
freq_table <- bind_rows(freq_table, total_row)

tabapp_output <- xtable(freq_table)
print(tabapp_output, include.rownames = FALSE)

