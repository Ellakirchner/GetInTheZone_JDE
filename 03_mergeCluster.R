#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Author: Ella Kirchner, Elinor Benami, Andrew Hobbs 
# Project name: RS Get in the Zone
# Date Last Updated: # Fri Sep 12 10:12:45 2025 
# Purpose: Perform clustering, transform clustering into insurance zones, visualize process, assess cluster stability and add/save insurance zone assignment
# Input:  all to be loaded from the RProj
        # - Six .RData files, one for each agri-environmental input, based on 50k points and OAF crop cuts (suffix "_2")
        # - Three .RData files for specific visualizations:
        #     - Two .RData files on soil data for separate point databases (suffix "_oaf", suffix "_50k")
        #     - One .RData file on precip data for oaf points only (suffix "_oaf")
        # - Four .RData files on admin data
        # - One .RData file on OAF data (prior to/without insurance zone assignments)
# Output Files: 
        # - (processed) OAF dataset with cluster assignment of k=14
        # - (processed) OAF dataset with only chirps cluster assignment based only on OAF points for clustering
        # - Insurance zone assignments on ward level (one .RData for every agri-environmental indicator, prefix "dynamic_")
        # - Insurance zone assignments on ward level using lowertail approach (three .RData for every agri-environmental indicator, suffixes "lowertail", "lowerlowertail", "log")
        # - 6 .Rdata files that contain cluster assignment at gridcell and at ward level, for k=14 only
        # - Fig. A.3: NMI heatmap on pairwise comparison at k=14
        # - Tab. A.1:Table for cluster stability assessment over time
# ReadMe: 1) Adapt file path to work in R-project, then loading of datasets should work
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# %%%%%%%%%%%%%%%%%%%%%%%%%%%
# Set Up 
# %%%%%%%%%%%%%%%%%%%%%%%%%%%

# 1. Clean out workspace
rm(list = ls(all = TRUE))

# 2. Set Up File path & directory: 
file_path <-  "put_your_file_path_here" # File path to project folder

#setwd(file_path)
saving <- paste0(file_path,"/figs")

# 3. Set Seed for replication 
set.seed(123456789)

# %%%%%%%%%%%%%%%%%%%%%%%%%%%
# Load libraries 
# %%%%%%%%%%%%%%%%%%%%%%%%%%%
library(dplyr)
library(ClustGeo)
library(spdep)
library(proxy)
library(tictoc)
library(sf)  
library(tidyverse)
library(RColorBrewer)
library(gtable)
library(cowplot)
library(grid)
library(gridExtra)

# For cluster stability assessment
library(mclust)
library(aricode)
library(writexl)

# Parallelize
library(parallel)

# %%%%%%%%%%%%%%%%%%%%%%%%%%%
# Load datasets ---- 
# %%%%%%%%%%%%%%%%%%%%%%%%%%%

# Admin data
load("datasets/admin3_clipped.RData")
load("datasets/admin3_kenya.RData")
load("datasets/admin1_clipped.RData")
load("datasets/admin1_kenya.RData")

# OAF 
load("datasets/oaf_noclusters.RData")

# Spatial data
load("datasets/ndvi_2.RData")
load("datasets/evi_2.RData")
load("datasets/chirp_2.RData")
load("datasets/temp_2.RData")
load("datasets/soil_2.RData")
load("datasets/all_2.RData")

# Dataset needed to create Figs for Appendix:
load("datasets/chirp_oaf.RData")
load("datasets/soil_oaf.RData")
load("datasets/soil_50k.RData")


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# A) Create Grid Layer for Clustering ----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# The points are the same across all attribute datasets
# Hence, we could take any of the attribute datasets to establish the grid cell network 
# We take soil 


### 1) Lay a grid network over ROI ----
# Before that: transform the projection 
points_proj <- st_transform(soil_2, crs = 32637)  # UTM Zone 37N (Kenya)
# points_proj <- st_transform(soil_oaf, crs = 32637)  # UTM Zone 37N (Kenya)

# Lay the grids
grid <- st_make_grid(points_proj, cellsize = 2000, square = TRUE)
grid_sf <- st_sf(geometry = grid)
# Transform back to previous projection
grid_sf <- st_transform(grid_sf, crs = st_crs(soil_2))  # Back to original CRS
# grid_sf <- st_transform(grid_sf, crs = st_crs(soil_oaf))  # Back to original CRS


### 2) Filter to ROI grid cells ----

# Add a grid cell ID 
grid_sf$grid_id <- 1:nrow(grid_sf)  # Create a unique ID for each grid cell

# Join points to grid cells (this gives each point a grid cell)
points_joined <- st_join(soil_2, grid_sf, join = st_within)
# points_joined <- st_join(soil_oaf, grid_sf, join = st_within)

# Count the points within each grid cell
grid_with_points <- points_joined %>%
  group_by(grid_id) %>%  # Group by grid cell ID
  summarise(point_count = n()) %>%  # Count the number of points in each grid cell
  filter(point_count > 0)  # Keep only grid cells with at least one point

# Join the grid cells with points count back to the original grid_sf and filter for at least 1 point per grid cell
grid_with_points_geom <- grid_sf %>%
  filter(grid_id %in% grid_with_points$grid_id) %>% 
  left_join(st_drop_geometry(grid_with_points), by = "grid_id") %>%  # Use left_join instead of st_join
  mutate(point_bin = cut(point_count,
                         breaks = c(1, 5, 10, 15, 20, 25, Inf),  # Changed 1 to 0
                         labels = c("1–5", "5–10", "10–15", "15–20", "20–25", ">25"),
                         right = FALSE))  # Changed to right = FALSE

# Check for NA values
sum(is.na(grid_with_points_geom$point_bin))

# Sum stats on points per grid cell: 
summary(grid_with_points$point_count)

### 3) Visual check ----

custom_colors <- c("yellow", "orange", "red", "purple", "blue", "black")

# Plot by count category
ggplot() +
  geom_sf(data = grid_with_points_geom, aes(fill = point_bin), color = "NA", size = 0.2) +
  scale_fill_manual(values = custom_colors, name = "Point Count") +
  geom_sf(data = admin3_clipped, fill = NA, color = "lightgrey", linewidth = 0.6) +
  geom_sf(data = admin1_clipped, fill = NA, color = "black", linewidth = 0.9) +
  theme_minimal() +
  labs(title = "Grid Cells by Point Count Category")

print(length(grid_with_points_geom$grid_id))

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# B) Clustering ----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Define functions for clustering (four separate functions):

### 1.1) Preprocessing function ----
# calculates distance matrices that are needed for the clustering
# (Done in a separate function to prevent recalculating distance matrices for every k)

prepare_clustering_data <- function(keyword, df, func=NULL) {
  # 1) Join spatial information to grids 
  # Ensure proper identification of variables: 
  colnames(df)[colnames(df) != "geometry"] <- paste0(keyword, "_", colnames(df)[colnames(df) != "geometry"])
  
  # Ensure that attribute dataset is a sf object: 
  if (!inherits(df, "sf")) {
    stop("Input data is not a spatial (sf) object.")
  } else {
    message("Input data is a spatial object, continue processing.")
  }
  
  # Spatial join of attribute (observations) with grid cells
  attr_idgrid <- df %>% 
    st_join(grid_with_points_geom, join = st_within) %>%
    select(-point_count, -point_bin)
  
  # Group by grid cell and calculate the mean for each variable
  attr_avggrid <- attr_idgrid %>%
    group_by(grid_id) %>%  # Group by grid cell ID
    summarise(across(starts_with(keyword), \(x) mean(x, na.rm = TRUE)))  
  
  # Join back the aggregated data to the grid geometry
  attr_avggrid_df <- st_drop_geometry(attr_avggrid) # drop geometry structure first to be able to join on grid id
  grid_avg_attr <- grid_with_points_geom %>%
    left_join(attr_avggrid_df, by = "grid_id")
  
  # Rename dataset for easier handling: 
  geocluster_attr <- grid_avg_attr
  
  # Extract attributes first
  attr_raw <- geocluster_attr %>%
    st_drop_geometry() %>%
    select(starts_with(keyword))
  
  # Remove grid cells with all missing values
  all_na_rows <- apply(attr_raw, 1, function(x) all(is.na(x)))
  if (sum(all_na_rows) > 0) {
    message(paste("Found", sum(all_na_rows), "grid cells with all missing values. Removing them."))
    geocluster_attr <- geocluster_attr[!all_na_rows, ]
    attr_raw <- attr_raw[!all_na_rows, ]
  }
  
  ### 2) Define dissimilarity matrices (EXPENSIVE OPERATIONS DONE ONCE)
  
  # Dissimilarity of attributes: 
  # Scale attributes
  attr <- scale(attr_raw)

  # lowertail transformation (if func is not NULL)
  if (!is.null(func)) {
    attr <- matrix(func(as.vector(attr)), nrow = nrow(attr), ncol = ncol(attr))
  }
  
  D0 <- dist(attr) # Feature matrix

  # Spatial dissimilarity
  coords <- st_coordinates(st_centroid(geocluster_attr))
  nb <- knn2nb(knearneigh(coords, k = 10))
  lw <- nb2listw(nb, style = "B", zero.policy = TRUE)
  mat <- listw2mat(lw)
  D1 <- as.dist(1-mat) # Geographic distances between the grid centroids

  # Return preprocessing results
  return(list(
    geocluster_attr = geocluster_attr,
    D0 = D0,
    D1 = D1
  ))
  message("Matrices computed.")
}

### 1.2) Clustering function ----
# (uses precomputed distance matrices)
perform_clustering <- function(prep_data, k, alpha) {
  # 3) Clustering grid cells using precomputed distance matrices
  tree <- hclustgeo(prep_data$D0, prep_data$D1, alpha)
  clusters <- cutree(tree, k = k)
  
  # Add to sf object
  geocluster_attr <- prep_data$geocluster_attr
  geocluster_attr$cluster <- as.factor(clusters)
  
  # Assign cluster value to ward and extract additional metrics
  # Ensure both datasets have the same CRS
  if (st_crs(admin3_clipped) != st_crs(geocluster_attr)) {
    geocluster_attr <- st_transform(geocluster_attr, st_crs(admin3_clipped))
  }
  
  # Spatial join to find which grid cells are within which admin3 units
  intersections <- st_intersects(admin3_clipped, geocluster_attr)
  
  # Create a new dataframe with the results
  result <- admin3_clipped %>%
    mutate(
      geocluster_attr_count = sapply(intersections, length),
      dominant_cluster = sapply(intersections, function(indices) {
        if (length(indices) == 0) return(NA)
        clusters <- geocluster_attr$cluster[indices]
        tab <- table(clusters)
        names(tab)[which.max(tab)]
      }),
      dominant_cluster = as.numeric(dominant_cluster),
      dominant_cluster_count = sapply(seq_along(intersections), function(i) {
        indices <- intersections[[i]]
        if (length(indices) == 0) return(0)
        sum(geocluster_attr$cluster[indices] == dominant_cluster[i])
      })
    ) %>%
    mutate(
      dominant_cluster_pct = ifelse(geocluster_attr_count > 0, 
                                    dominant_cluster_count / geocluster_attr_count * 100, 
                                    0),
      dominant_cluster_area = ifelse(geocluster_attr_count > 0, 
                                     dominant_cluster_count * 4000000, 
                                     0),
      dominant_cluster_area_pct = ifelse(geocluster_attr_count > 0,
                                         as.numeric(dominant_cluster_area/admin3_area),
                                         0)
    )
  
  # Return both dataframes as a named list
  return(list(
    ward_clusters = result,
    grid_clusters = geocluster_attr
  ))
}

### 1.3) fct for single k value ----
generate_single_cluster <- function(keyword, df, k = 14, alpha, func=NULL) {
  prep_data <- prepare_clustering_data(keyword, df, func=func)
  return(perform_clustering(prep_data, k, alpha))
}

### 1.4) fct for multiple k values ----
generate_multi_clusters <- function(keyword, df, k_values, alpha, func=NULL) {
  cat("Preparing data (calculating distance matrices)...\n")
  prep_data <- prepare_clustering_data(keyword, df, func=func)
  
  cat("Running clustering for k values:", paste(k_values, collapse = ", "), "\n")
  
  results_list <- lapply(k_values, function(k) {
    cat("Processing k =", k, "\n")  
    perform_clustering(prep_data, k, alpha)
  })
  
  names(results_list) <- paste0("k_", k_values)
  return(results_list)
}

### 1.5) fct to extract cluster list (ward-level) ----
# Based on generate_multi_clusters output
extract_cluster_results <- function(multi_cluster_results) {
  
  # Get the first result to extract name_concat
  first_result <- multi_cluster_results[[1]]$ward_clusters
  
  # Start with name_concat column
  out <- first_result %>%
    st_drop_geometry() %>%
    select(name_concat)
  
  # Add cluster columns for each k value
  for (k_name in names(multi_cluster_results)) {
    cluster_col <- multi_cluster_results[[k_name]]$ward_clusters$dominant_cluster
    out <- out %>%
      mutate(!!paste0("cluster_", gsub("k_", "", k_name)) := cluster_col+4) # we add 4 here as k only starts at k=5 and not at k=1
  }
  
  return(out)
}

### 1.6) Functions for lowertail clustering ---- 
lowertail <- function(seq, threshold=0){
  return(as.numeric(seq<=threshold))
}

modlog <- function(x){
  return(log(x+8))
}

identity <- function(x){
  return(x)
}

lowerlowertail <- function(x) lowertail(x, threshold = -1.036)
# the number above is the z-score associated with the 15th percentile, 
# roughly where the Estefania-Salazar paper got the best results


### 2) Define parameters ---- 
alpha <- 0.1
k <- 14 # for static analysis
k_values <- 5:100 # for dynamic analysis

admin3_clipped <- 
  admin3_clipped %>% 
  select(name_concat, name_1, geometry)

admin3_clipped$admin3_area <- st_area(admin3_clipped)


### 3) Filter input datasets to years of interest ----
chirp_2 <- chirp_2 %>% select(-contains(c("2016", "2017", "2018", "2019", "2020")))
# 129 var (minus geometry)
evi_2 <- evi_2 %>% select(-contains(c("2016", "2017", "2018", "2019", "2020")))
# 106 var (minus geometry, not summing up to 112 due to early cut off, hence some years only have one var in May instead of two)
ndvi_2 <- ndvi_2 %>% select(-contains(c("2016", "2017", "2018", "2019", "2020")))
# 106 var (minus geometry, not summing up to 112 due to early cut off, hence some years only have one var in May instead of two)
temp_2 <- temp_2 %>% select(-contains(c("2016", "2017", "2018", "2019", "2020")))
# 257 var (minus geometry)
all_2 <- all_2 %>% select(-contains(c("2016", "2017", "2018", "2019", "2020")))


### 4) Generate dynamic clusters (k=5,...100) ---- 

### 4.1) apply fct to gen clusters ----
multi_chirps <- generate_multi_clusters("chirps", chirp_2, k_values, alpha)
multi_temp <- generate_multi_clusters("temp", temp_2, k_values, alpha)
multi_ndvi <- generate_multi_clusters("ndvi", ndvi_2, k_values, alpha)
multi_evi <- generate_multi_clusters("evi", evi_2, k_values, alpha)
multi_soil <- generate_multi_clusters("soil", soil_2, k_values, alpha)
multi_all <- generate_multi_clusters("all", all_2, k_values, alpha)

multi_chirps_log <- generate_multi_clusters("chirps", chirp_2, k_values, alpha, func=modlog)
multi_temp_log <- generate_multi_clusters("temp", temp_2, k_values, alpha, func=modlog)
multi_ndvi_log  <- generate_multi_clusters("ndvi", ndvi_2, k_values, alpha, func=modlog)
multi_evi_log  <- generate_multi_clusters("evi", evi_2, k_values, alpha, func=modlog)
multi_soil_log  <- generate_multi_clusters("soil", soil_2, k_values, alpha, func=modlog)
multi_all_log  <- generate_multi_clusters("all", all_2, k_values, alpha, func=modlog)

multi_chirps_lowertail <- generate_multi_clusters("chirps", chirp_2, k_values, alpha, func=lowertail)
multi_temp_lowertail <- generate_multi_clusters("temp", temp_2, k_values, alpha, func=lowertail)
multi_ndvi_lowertail  <- generate_multi_clusters("ndvi", ndvi_2, k_values, alpha, func=lowertail)
multi_evi_lowertail  <- generate_multi_clusters("evi", evi_2, k_values, alpha, func=lowertail)
multi_soil_lowertail  <- generate_multi_clusters("soil", soil_2, k_values, alpha, func=lowertail)
multi_all_lowertail  <- generate_multi_clusters("all", all_2, k_values, alpha, func=lowertail)

multi_chirps_lowerlowertail <- generate_multi_clusters("chirps", chirp_2, k_values, alpha, func=lowerlowertail)
multi_temp_lowerlowertail <- generate_multi_clusters("temp", temp_2, k_values, alpha, func=lowerlowertail)
multi_ndvi_lowerlowertail  <- generate_multi_clusters("ndvi", ndvi_2, k_values, alpha, func=lowerlowertail)
multi_evi_lowerlowertail  <- generate_multi_clusters("evi", evi_2, k_values, alpha, func=lowerlowertail)
multi_soil_lowerlowertail  <- generate_multi_clusters("soil", soil_2, k_values, alpha, func=lowerlowertail)
multi_all_lowerlowertail  <- generate_multi_clusters("all", all_2, k_values, alpha, func=lowerlowertail)


### 4.2) apply fct to gen clean df for dynamicCE analysis ----
# For the dynamic zone count, we need a df that has the wards as observations 
# and the cluster assignments depending on the k value as columns 
# We write a function for that and apply that to the previously generated multi_attr datasets

dynamic_chirps <- extract_cluster_results(multi_chirps)
dynamic_temp <- extract_cluster_results(multi_temp)
dynamic_ndvi <- extract_cluster_results(multi_ndvi)
dynamic_evi <- extract_cluster_results(multi_evi)
dynamic_soil <- extract_cluster_results(multi_soil)
dynamic_all <- extract_cluster_results(multi_all)

#For Lowertail clusters: 
dynamic_chirps_log <- extract_cluster_results(multi_chirps_log)
dynamic_temp_log <- extract_cluster_results(multi_temp_log)
dynamic_ndvi_log <- extract_cluster_results(multi_ndvi_log)
dynamic_evi_log <- extract_cluster_results(multi_evi_log)
dynamic_soil_log <- extract_cluster_results(multi_soil_log)
dynamic_all_log <- extract_cluster_results(multi_all_log)

dynamic_chirps_lowertail <- extract_cluster_results(multi_chirps_lowertail)
dynamic_temp_lowertail <- extract_cluster_results(multi_temp_lowertail)
dynamic_ndvi_lowertail <- extract_cluster_results(multi_ndvi_lowertail)
dynamic_evi_lowertail <- extract_cluster_results(multi_evi_lowertail)
dynamic_soil_lowertail <- extract_cluster_results(multi_soil_lowertail)
dynamic_all_lowertail <- extract_cluster_results(multi_all_lowertail)

dynamic_chirps_lowerlowertail <- extract_cluster_results(multi_chirps_lowerlowertail)
dynamic_temp_lowerlowertail <- extract_cluster_results(multi_temp_lowerlowertail)
dynamic_ndvi_lowerlowertail <- extract_cluster_results(multi_ndvi_lowerlowertail)
dynamic_evi_lowerlowertail <- extract_cluster_results(multi_evi_lowerlowertail)
dynamic_soil_lowerlowertail <- extract_cluster_results(multi_soil_lowerlowertail)
dynamic_all_lowerlowertail <- extract_cluster_results(multi_all_lowerlowertail)

### 5) Extract k=14 and merge to oaf ----
# As there are 14 admin1 zones in the study region, we use the insurance zone design 
# with 14 zones (k=14) as a common visualization and assessment unit. Hence, 
# there is value in having these insurance zones included in the base dataset. 

chirps_tomerge <- dynamic_chirps %>% mutate(c.chirp = cluster_14) %>% select(name_concat, c.chirp)
ndvi_tomerge <- dynamic_ndvi %>%mutate(c.ndvi = cluster_14) %>% select(name_concat, c.ndvi)
evi_tomerge <- dynamic_evi %>% mutate(c.evi = cluster_14) %>% select(name_concat, c.evi)
temp_tomerge <- dynamic_temp %>% mutate(c.temp = cluster_14) %>% select(name_concat, c.temp)
soil_tomerge <- dynamic_soil %>% mutate(c.soil = cluster_14) %>% select(name_concat, c.soil)
all_tomerge <- dynamic_all %>% mutate(c.all = cluster_14) %>% select(name_concat, c.all)

oaf <- oaf %>% 
    select(-starts_with("c.")) %>% # delete variables first in case the dataset was already saved with cluster variables before
    left_join(chirps_tomerge, by = "name_concat") %>%
    left_join(ndvi_tomerge, by = "name_concat") %>%
    left_join(evi_tomerge, by = "name_concat") %>%
    left_join(temp_tomerge, by = "name_concat") %>%
    left_join(soil_tomerge, by = "name_concat") %>%
    left_join(all_tomerge, by = "name_concat")

# In addition, we extract them separately to have the grid cell assignment and the 
# ward assignment in one R object to be able to visualize these in section E of 
# this script. This separate function ensure to keep both, grid cell and wards assignments. 
cluster_k14_chirps <- generate_single_cluster("chirp", chirp_2, k, alpha)
cluster_k14_temp <- generate_single_cluster("temp",temp_2, k , alpha)
cluster_k14_ndvi <- generate_single_cluster("ndvi", ndvi_2, k , alpha)
cluster_k14_evi <- generate_single_cluster("evi", evi_2, k, alpha)
cluster_k14_soil <- generate_single_cluster("soil", soil_2, k, alpha)
cluster_k14_all <- generate_single_cluster("all", all_2, k, alpha)

# Assess share of wards that has only one constituent grid cell category (footnote 9 in Section 3.1.2)
wards <- cluster_k14_chirps$ward_clusters
grids <- cluster_k14_chirps$grid_clusters

wards <- wards %>%
  mutate(binary = 0)%>%
  mutate(binary = ifelse(geocluster_attr_count == dominant_cluster_count, 1, binary)) %>%
  mutate(share = (dominant_cluster_count/geocluster_attr_count))
table(wards$binary)
summary(wards$share)


### 6) Generate yield clusters ----
#Create name_concat for admin2: 
oaf$name_concat2 <- paste0(oaf$name_1, "_", oaf$name_2)
length(unique(oaf$name_concat2)) 
admin3_clipped$name_concat2 <- paste0(admin3_clipped$name_1, "_", admin3_clipped$name_2)
length(unique(admin3_clipped$name_concat2)) # 101 distinct subcounties

# Ensure data is prepared as in welfare analysis 
# Winsorize maximum yield (roughly 8t is the maximum realistic yield):
quantiles <- quantile(oaf$yield, c(0.01, 0.05, 0.1, 0.25, 0.5,0.75,0.9,0.95, 0.99))
quantiles # 99% equals roughly 8t, hence, we use 99% 
oaf_sorted <- arrange(oaf, yield)
summary(oaf$yield)

oaf$yield <- pmin(oaf$yield, quantile(oaf$yield, 0.99))
# Winsorize minimum yield to same percentile (1%):
oaf$yield <- pmax(oaf$yield, quantile(oaf$yield, 0.01))

summary(oaf$yield)

# Base yield-oracle clusters on wards
# (Instead of grid cells we have to use wards for aggregation and clustering as 
# grid cells do not have missing data in many years. Only 46 grid cells with 
# yield records in all years, hence, rather use 184 wards as used for welfare analysis)

# Multiple actions:  
# - Keep only wards with yield data in every year
# - Select relevant variables: ward, year, yield 
# - Take average by ward and year
# - Transform into wide format 
# - Merge admin data back to it 
all_years <- unique(oaf$year)

oaf_wide <- 
  oaf %>%
  st_drop_geometry()%>% #required, otherwise the pivot wide transformation doesn't work
  group_by(name_concat) %>%
  filter(length(unique(year)) == length(all_years)) %>% # drops 7328k observations 
  ungroup()

#Create name_concat for admin2: 
oaf_wide$name_concat2 <- paste0(oaf_wide$name_1, "_", oaf_wide$name_2)
length(unique(oaf_wide$name_concat2)) # 67 after filtering

oaf_wide <- oaf_wide %>%
  select(name_concat, year, yield)%>%
  group_by(name_concat, year) %>%
  summarize(avg_yield = mean(yield), .groups = 'drop')%>%
  pivot_wider(
    id_cols = name_concat,  
    names_from = year,
    values_from = avg_yield,
    names_prefix = "yield_"
  )%>%
  left_join(admin3_clipped, by="name_concat")%>%
  st_as_sf() %>%
  select(-name_1)#, -admin3_area)

# Prep data for clustering 
# - drop variables not needed for clustering
# - scale variables 
yield_clustering <- oaf_wide %>%
  st_drop_geometry()%>%
  select(name_concat, starts_with("yield"))%>%
  column_to_rownames("name_concat") %>%  # Move name_concat to rownames
  scale() 

# Feature matrix
D0 <- dist(yield_clustering) 

# Calc spatial dissimilarity matrix
coords <- st_coordinates(st_centroid(oaf_wide))
nb <- knn2nb(knearneigh(coords, k = 10))
lw <- nb2listw(nb, style = "B", zero.policy = TRUE)
mat <- listw2mat(lw) 
D1 <- as.dist(1-mat) # Geographic distances between the grid centroids

# Perform clustering hard coded: 
tree <- hclustgeo(D0, D1, alpha)
k = 14
clusters_14 <- cutree(tree, k = k)
table(clusters_14) # check to see how many observations, Should be length k with all > 0

# We choose to design 67 clusters and not 100 as these insurance zones based on 
# the observed yields are not spanning the entire study region as the other data-driven
# insurance zones would. The entire study region would have 100 admin2 zones, but 
# the area for which the yield data is valid, only spans 67 admin2 zones.
k = 67
clusters_67 <- cutree(tree, k = k)
table(clusters_67)

oaf_wide$c.yield <- as.numeric(clusters_14)
oaf_wide$c67.yield <- as.numeric(clusters_67)

# Filter datasets for ward and cluster assignment only (to be used for merging to final dataset)
oaf_wide_wcluster <- 
  oaf_wide %>% 
  st_drop_geometry()%>%
  select(name_concat,c.yield, c67.yield)

oaf <- oaf %>% 
  select(-starts_with("c.yield"), -starts_with("c67.yield"))%>% # Clean in case some other variables had been merged before
  left_join(oaf_wide_wcluster, by = "name_concat")


### 7) Create clusters required for figure A.5 in appendix ----
# Resulting dataset is needed in script 04_calcDescriptives where Fig. A.5 is created

### 7.1) Data prep for clustering based on filtered oaf data (13k) ----

# Filter oaf data to wards we want to keep:
# 1 Prep oaf data to filter by unique id:
# Check that yield data are filtered (exclude yield=0) and winsorized
sum(oaf$yield ==0)
summary(oaf$yield)

# Filter to keep only name_concat values that appear in all years
all_years <- unique(oaf$year)
oaf_filtered <- oaf %>%
  group_by(name_concat) %>%
  filter(length(unique(year)) == length(all_years)) %>% # drops 7k observations
  ungroup()
obstokeep <- oaf_filtered%>% select(unique_id) %>% st_drop_geometry()

# Left join based on unique_id:
chirp_oaf <- obstokeep %>% left_join(chirp_oaf, by="unique_id")%>% select(-unique_id)%>% st_as_sf()

# Filter input datasets to years of interest:
chirp_oaf <- chirp_oaf %>% select(-contains(c("2016", "2017", "2018", "2019", "2020")))

# Check for missing values 
anyNA(chirp_oaf)
sum(is.na(chirp_oaf))
colSums(is.na(chirp_oaf))
# Clean out empty rows: 
chirp_oaf_clean <- chirp_oaf %>%
  filter(!if_all(-geometry, is.na))

# Check again: 
anyNA(chirp_oaf_clean)

### 7.2) Apply fnct ----
k = 14
cluster_k14_oaf_chirps <- generate_single_cluster("chirp", chirp_oaf_clean, k, alpha)

### 7.3) Merge for separate analysis based only on oaf data----
oaf_woafdataclusters <- oaf_filtered %>%
  select(-starts_with("c."), -starts_with("c67")) %>% # delete variables first in case the dataset was already saved with cluster variables before
  left_join(cluster_k14_oaf_chirps$ward_clusters %>% st_drop_geometry() %>%
              select(name_concat, c.chirp = dominant_cluster),
            by = "name_concat") 

chirps <- cluster_k14_oaf_chirps$ward_clusters
sum(is.na(chirps$dominant_cluster))
491-157 # number of wards assigned to an insurance zone

### 8) Anonymize OAF data ----
oaf_anonym <- oaf %>% 
  select(-id_latlonyear) %>% 
  st_drop_geometry()

oaf_woafdataclusters_anonym <- oaf_woafdataclusters %>% 
  select(-id_latlonyear) %>% 
  st_drop_geometry()
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# C) Save datasets ----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# OAF data: 
save(oaf, file = paste0(file_path,"/datasets/oaf.RData"))
save(oaf_anonym, file = paste0(file_path, "/datasets/oaf_anonym.RData"))
save(oaf_woafdataclusters, file = paste0(file_path, "/datasets/oaf_woafdataclusters.RData"))
save(oaf_woafdataclusters_anonym, file = paste0(file_path, "/datasets/oaf_woafdataclusters_anonym.RData"))

# Save dataframes that contain clusters for dynamic CE analysis:
save(dynamic_chirps, file = paste0(file_path,"/datasets/dynamic_chirps.RData"))
save(dynamic_temp, file = paste0(file_path,"/datasets/dynamic_temp.RData"))
save(dynamic_ndvi, file = paste0(file_path,"/datasets/dynamic_ndvi.RData"))
save(dynamic_evi, file = paste0(file_path,"/datasets/dynamic_evi.RData"))
save(dynamic_soil, file = paste0(file_path,"/datasets/dynamic_soil.RData"))
save(dynamic_all, file = paste0(file_path,"/datasets/dynamic_all.RData"))

save(dynamic_chirps_log, file = paste0(file_path,"/datasets/dynamic_chirps_log.RData"))
save(dynamic_temp_log, file = paste0(file_path,"/datasets/dynamic_temp_log.RData"))
save(dynamic_ndvi_log, file = paste0(file_path,"/datasets/dynamic_ndvi_log.RData"))
save(dynamic_evi_log, file = paste0(file_path,"/datasets/dynamic_evi_log.RData"))
save(dynamic_soil_log, file = paste0(file_path,"/datasets/dynamic_soil_log.RData"))
save(dynamic_all_log, file = paste0(file_path,"/datasets/dynamic_all_log.RData"))

save(dynamic_chirps_lowertail, file = paste0(file_path,"/datasets/dynamic_chirps_lowertail.RData"))
save(dynamic_temp_lowertail, file = paste0(file_path,"/datasets/dynamic_temp_lowertail.RData"))
save(dynamic_ndvi_lowertail, file = paste0(file_path,"/datasets/dynamic_ndvi_lowertail.RData"))
save(dynamic_evi_lowertail, file = paste0(file_path,"/datasets/dynamic_evi_lowertail.RData"))
save(dynamic_soil_lowertail, file = paste0(file_path,"/datasets/dynamic_soil_lowertail.RData"))
save(dynamic_all_lowertail, file = paste0(file_path,"/datasets/dynamic_all_lowertail.RData"))

save(dynamic_chirps_lowerlowertail, file = paste0(file_path,"/datasets/dynamic_chirps_lowerlowertail.RData"))
save(dynamic_temp_lowerlowertail, file = paste0(file_path,"/datasets/dynamic_temp_lowerlowertail.RData"))
save(dynamic_ndvi_lowerlowertail, file = paste0(file_path,"/datasets/dynamic_ndvi_lowerlowertail.RData"))
save(dynamic_evi_lowerlowertail, file = paste0(file_path,"/datasets/dynamic_evi_lowerlowertail.RData"))
save(dynamic_soil_lowerlowertail, file = paste0(file_path,"/datasets/dynamic_soil_lowerlowertail.RData"))
save(dynamic_all_lowerlowertail, file = paste0(file_path,"/datasets/dynamic_all_lowerlowertail.RData"))


# Save dataframes that contain dataframe lists to assess cluster stability:
save(cluster_k14_chirps, file = paste0(file_path,"/datasets/cluster_k14_chirps.RData"))
save(cluster_k14_temp, file = paste0(file_path,"/datasets/cluster_k14_temp.RData"))
save(cluster_k14_ndvi, file = paste0(file_path,"/datasets/cluster_k14_ndvi.RData"))
save(cluster_k14_evi, file = paste0(file_path,"/datasets/cluster_k14_evi.RData"))
save(cluster_k14_soil, file = paste0(file_path,"/datasets/cluster_k14_soil.RData"))
save(cluster_k14_all, file = paste0(file_path,"/datasets/cluster_k14_all.RData"))


##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# D) Cluster stability assessment ----
##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
### 1) Pairwise comparison ----
# Generating heatmap presented in appendix A.3.2 

### 1.1) Create matrix of pairwise NMI values ----
# Define the variables to compare
cluster_vars <- c("name_1","c.chirp", "c.evi", "c.ndvi","c.soil", "c.temp", "c.all")

admin3_clusters <- admin3_clipped %>% st_drop_geometry() %>% select(name_concat, name_1) %>%
  left_join(chirps_tomerge, by="name_concat") %>% 
  left_join(ndvi_tomerge, by="name_concat") %>% 
  left_join(evi_tomerge, by="name_concat") %>% 
  left_join(temp_tomerge, by="name_concat") %>% 
  left_join(soil_tomerge, by="name_concat") %>% 
  left_join(all_tomerge, by="name_concat") %>%
  filter(!if_all(starts_with("c."), is.na))
  

# Create a matrix to store the NMI values
nmi_matrix <- matrix(NA, nrow = length(cluster_vars), ncol = length(cluster_vars),
                     dimnames = list(cluster_vars, cluster_vars))

# Loop over all pairs 
for (i in seq_along(cluster_vars)) {
  for (j in seq_along(cluster_vars)) {
    nmi_matrix[i, j] <- NMI(admin3_clusters[[cluster_vars[i]]],
                            admin3_clusters[[cluster_vars[j]]])
  }
}

# Convert to a data frame for easier viewing 
nmi_df <- as.data.frame(nmi_matrix)
# Add rownames as a column
nmi_df <- cbind(Cluster1 = rownames(nmi_df), nmi_df)
print(nmi_df)

### 1.2) Visualize matrix of NMI values ----
nmi_melted <- nmi_df %>% pivot_longer(cols = 2:8,names_to="cluster_2", values_to="nmi") %>% rename(cluster_1 = Cluster1)
# Rename cluster variables for visualization:
nmi_melted$cluster_1 <- factor(nmi_melted$cluster_1,
                               levels = c("c.all", "c.chirp", "c.evi", "c.ndvi", "c.temp", "c.soil","name_1"),
                               labels = c("All", "Precipitation", "EVI", "NDVI", "Temperature", "Soil", "Admin 1"))

nmi_melted$cluster_2 <- factor(nmi_melted$cluster_2,
                               levels = c("c.all", "c.chirp", "c.evi", "c.ndvi", "c.temp", "c.soil", "name_1"),
                               labels = c("All", "Precipitation", "EVI", "NDVI", "Temperature", "Soil", "Admin 1"))

nmi_upper <- nmi_melted %>%
  filter(as.numeric(cluster_1) <= as.numeric(cluster_2))
nmi_upper$cluster_1 <- factor(nmi_upper$cluster_1, levels = rev(levels(nmi_melted$cluster_2)))

nmi_plot <- ggplot(nmi_upper, aes(x = cluster_1, y = cluster_2, fill = nmi)) +
  geom_tile(color = "white") +
  geom_text(aes(label = round(nmi, 2)), size = 3) +
  scale_fill_gradient2(low = "#0072B2", mid= "white", high = "red", midpoint = 0.5, limits = c(0, 1)) +
  scale_x_discrete(position = "top") +
  theme_minimal(base_size = 12) +
  labs(title = "Normalized Mutual Information (NMI) between Cluster Assignments",
       x = "", y = "", fill = "NMI") +
  # theme(axis.text.x = element_text(angle = 45, hjust = 1))
  theme(
    title = element_blank(),
    # plot.title = element_text(hjust = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 0, vjust = 0.5),  # vjust controls vertical placement
    plot.margin = margin(t = 20, r = 10, b = 10, l = 10),              # adds space above plot
   )
nmi_plot

# 1.3) Save Fig.A.3 ----
ggsave(paste0(saving, "/nmi_plot.png"), plot= nmi_plot, width=5, height= 2.5)


### 2) Leave-one-year-out approach ----

### 2.1) Generate reduced datasets and clusters ----
k = 14

# Define years to be left out sequentially:
years <- 2002:2015

# Function to generate dataset with one year left out and clusters based on these reduced datasets:
minus_one_attr <- function(years, attr_data, keyword, k) {
  cluster_list <- lapply(years, function(year_out) {
    # Drop years after crop cuts and the year to be left out
    years_to_exclude <- c("2016", "2017", "2018", "2019", "2020", year_out)
    
    # Subset the data based on the years that we want to leave out
    attr_subset <- attr_data %>%
       select(
         -matches(paste0("(", paste(years_to_exclude, collapse = "|"), ")"))
       )     
    
    # Run the clustering as defined before (in the func under 7.2)
    cluster <- generate_single_cluster(keyword, attr_subset, k, alpha=alpha)
    
    # Extract only name_concat and dominant_cluster
    clust_result <- cluster$ward_clusters %>%
      st_drop_geometry() %>%
      select(name_concat, dominant_cluster)
    
    # Rename the dominant_cluster column
    clust_result <- clust_result %>%
      rename(!!paste0("Cluster_", year_out) := dominant_cluster) 
    
    return(clust_result)
  })
  
  # Reduce to a single wide dataframe by joining on name_concat
  cluster_wide <- reduce(cluster_list, left_join, by = "name_concat")
  
  return(cluster_wide)
}
clusters_chirpsall <- minus_one_attr(years, chirp_2, "chirp", k)
clusters_ndviall <- minus_one_attr(years, ndvi_2, "ndvi", k)
clusters_eviall <- minus_one_attr(years, evi_2, "evi", k)
clusters_tempall <- minus_one_attr(years, temp_2, "temp", k)
clusters_allall <- minus_one_attr(years, all_2, "all", k)

# Join with clusters based on full dataset: 
oaf_ward_chirpclusters <- oaf_ward%>% 
  select(name_concat, cluster = c.chirp) %>% 
  left_join(clusters_chirpsall, by="name_concat")

oaf_ward_tempclusters <- oaf_ward%>% 
  select(name_concat, cluster = c.temp) %>% 
  left_join(clusters_tempall, by="name_concat")

oaf_ward_ndviclusters <- oaf_ward%>% 
  select(name_concat, cluster = c.ndvi) %>% 
  left_join(clusters_ndviall, by="name_concat")

oaf_ward_eviclusters <- oaf_ward%>% 
  select(name_concat, cluster = c.evi) %>% 
  left_join(clusters_eviall, by="name_concat")

oaf_ward_allclusters <- oaf_ward%>% 
  select(name_concat, cluster = c.all) %>% 
  left_join(clusters_allall, by="name_concat")

### 2.2 Assess NMIs ----
# We compare the clusters based on the full dataset to those clusters based on the data for which we left out one year

# Get all relevant year columns
cluster_cols <- paste0("Cluster_", years)

# Define a helper function to compute NMI for a given dataset
compute_nmi <- function(df, prefix) {
  nmis <- sapply(cluster_cols, function(col) {
    NMI(df[[col]], df$cluster)
  })
  names(nmis) <- paste0(prefix, "_", cluster_cols)
  return(nmis)
}

# Compute NMIs for all datasets with unique prefixes
nmi_temp <- compute_nmi(oaf_ward_tempclusters, "temp")
nmi_ndvi <- compute_nmi(oaf_ward_ndviclusters, "ndvi")
nmi_chirps <- compute_nmi(oaf_ward_chirpclusters, "chirps")
nmi_evi <- compute_nmi(oaf_ward_eviclusters, "evi")
nmi_all <- compute_nmi(oaf_ward_allclusters, "all")

### 2.3) Combine in table and export ----
# Combine all into a single named vector
nmi_bind <- c(nmi_temp, nmi_ndvi, nmi_chirps, nmi_evi, nmi_all)

#Export to Excel
nmi_df <- data.frame(
  Comparison = names(nmi_bind),
  NMI = as.numeric(nmi_bind),
  row.names = NULL
)

# Write to Excel
write_xlsx(nmi_df, paste0(saving, "/nmi_results_leaveoneout.xlsx"))

##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# E) Visualization of clusters ----
##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# This section creates Figure 3 in the manuscript
# and Fig. A.4 (Overview of clusters) in the appendix. 

### 1) Preparations ----
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
  "#AEB6BF"   # Taupe
)

# Function to create combined cluster visualization plots
create_cluster_plots <- function(cluster_data, keyword, is_chirps=FALSE) {
  
  # Extract grid and ward clusters from the cluster data
  grids <- cluster_data$grid_clusters
  ward <- cluster_data$ward_clusters
  n <- nrow(grids)
  
  # Create grid-level plot (no legend)
  vis_grid <- ggplot() +
    geom_sf(data = grids, aes(fill = as.factor(cluster)), color = NA, size = 0) +
    geom_sf(data = admin3_clipped, fill = NA, color = "grey35", linewidth = 0.6) +
    geom_sf(data = admin1_clipped, fill = NA, color = "black", linewidth = 0.9) +
    scale_fill_manual(values = cluster_colors, guide = "none") +
    theme_minimal()+ 
    theme(plot.margin = margin(t = -15, r = 5, b = 5, l = 5))
  
  
  if(is_chirps){
  # Only for CHIRPS (as that is the plot that we are displaying in the manuscript and not in the overview in the appendix)
  # Create ward-level plot (with legend for extraction of CHIRPS graph which is separate from overview and hence requires its legend)
  vis_ward_withlegend <- ggplot() +
    geom_sf(data = ward, aes(fill = as.factor(dominant_cluster)), color = NA, size = 0) +
    geom_sf(data = admin3_clipped, fill = NA, color = "grey35", linewidth = 0.6) +
    geom_sf(data = admin1_clipped, fill = NA, color = "black", linewidth = 0.9) +
    scale_fill_manual(values = cluster_colors, name = "Cluster/\nZone") +
    theme_minimal()+
    theme(plot.margin = margin(t = -15, r = 5, b = 5, l = 5))
  
  # Extract legend (to separate the legend in the CHIRPS graph)
  legend <- gtable_filter(ggplot_gtable(ggplot_build(vis_ward_withlegend)), "guide-box")
  } 
  
  # Create ward-level plot without legend (for all other plots as they share the same legend as the yield cluster graph)
  vis_ward <- ggplot() +
    geom_sf(data = ward, aes(fill = as.factor(dominant_cluster)), color = NA, size = 0) +
    geom_sf(data = admin3_clipped, fill = NA, color = "grey35", linewidth = 0.6) +
    geom_sf(data = admin1_clipped, fill = NA, color = "black", linewidth = 0.9) +
    scale_fill_manual(values = cluster_colors, guide = "none") +
    theme_minimal() +
    theme(plot.margin = margin(t = -15, r = 5, b = 5, l = 5))
  
    # Create subtitle grobs
  subtitle_a <- textGrob("a) Clustered grid cells\n", 
                         gp = gpar(fontsize = 12), 
                         hjust = 0, x = 0.03)
  subtitle_b <- textGrob("b) Insurance zones\n(= clustered grid cells aggregated at ward level)", 
                         gp = gpar(fontsize = 12), 
                         hjust = 0, x = 0.03)
  
  # Define main heading
  heading <- paste0("Based on ", keyword)
  
  # Create layout: heading, subtitles, maps 
  if(is_chirps){
  plot <- grid.arrange(
    subtitle_a, subtitle_b, nullGrob(),  # Third column for legend space
    vis_grid, vis_ward, legend,
    ncol = 3,
    nrow = 2,
    widths = c(1, 1, 0.3),  # Equal map widths, small space for legend
    heights = c(1, 10))
    # top = textGrob(heading, gp = gpar(fontsize = 14, fontface = "bold"))
  }
  else
  {
  plot <- grid.arrange(
    subtitle_a, subtitle_b,
    vis_grid, vis_ward,
    ncol = 2,
    nrow = 2,
    widths = c(1, 1),
    heights = c(1, 10),
    top = textGrob(heading, gp = gpar(fontsize = 14, fontface = "bold"))
  )
  }
}

### 2) Generate and save Fig3 ----
plot_combined_chirps <- create_cluster_plots(cluster_k14_chirps, "precipitation", is_chirps = TRUE) #separate code!!

ggsave(paste0(saving,"/chirps_clusters.png"), plot= plot_combined_chirps, width=7.5, height= 6)

### 3) Generate separate plots required for Fig.A.4 ----
plot_combined_evi <- create_cluster_plots(cluster_k14_evi, "EVI", is_chirps=FALSE)
plot_combined_ndvi <- create_cluster_plots(cluster_k14_ndvi, "NDVI", is_chirps=FALSE)
plot_combined_temp <- create_cluster_plots(cluster_k14_temp, "temperature", is_chirps=FALSE)
plot_combined_soil <- create_cluster_plots(cluster_k14_soil, "soil", is_chirps=FALSE)
plot_combined_all <- create_cluster_plots(cluster_k14_all, "all indicators combined", is_chirps=FALSE)

# Plot of yield insurance zones: 
# Create ward-level plot (with legend for extraction)
# First, add sf features back to yield cluster dataset: 
oaf_wide <- 
  oaf_wide %>%
  left_join(admin3_clipped%>%st_drop_geometry(), by="name_concat")%>%
  st_sf()

# # Ward level yield-based zones without legend: 
vis_ward_yield_nolegend <- ggplot() +
  geom_sf(data = oaf_wide, aes(fill = as.factor(c.yield)), color = NA, size = 0) +
  geom_sf(data = admin3_clipped, fill = NA, color = "grey35", linewidth = 0.6) +
  geom_sf(data = admin1_clipped, fill = NA, color = "black", linewidth = 0.9) +
  scale_fill_manual(values = cluster_colors, guide = "none") +
  theme_minimal()+
  labs(title="Insurance zones based on yields",
       subtitle = "(Clustering done at ward level)\n"
  ) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 12)  # Optional: set subtitle size
  )

# Create ward-level plot (with legend for extraction of NDVI graph to ensure that the legend shows NA as well)
wards <- cluster_k14_ndvi$ward_clusters
legend_aux <- ggplot() +
  geom_sf(data = wards, aes(fill = as.factor(dominant_cluster)), color = NA, size = 0) +
  geom_sf(data = admin3_clipped, fill = NA, color = "grey35", linewidth = 0.6) +
  geom_sf(data = admin1_clipped, fill = NA, color = "black", linewidth = 0.9) +
  scale_fill_manual(values = cluster_colors, name = "Cluster/\nZone") +
  theme_minimal()

# Extract legend from NDVI graph
legend <- gtable_filter(ggplot_gtable(ggplot_build(legend_aux)), "guide-box")


### 4) Combine separate graphs into overview ----
# Create empty plots as spacer
spacer <- rectGrob(gp = gpar(col = "white", fill = "white"))

# Combine the yield-based clusters and the legend in one feature that can then be included in the lower right: 
lower_right_box <- arrangeGrob(vis_ward_yield_nolegend, spacer, legend, 
                               nrow = 1, ncol = 4, 
                               widths = c(0.7, 0.03, 0.2, 0.07))

# Overview of all (except for CHIRPS) with spacers: 
overview <- grid.arrange(
  plot_combined_evi, spacer, plot_combined_ndvi, 
  spacer, spacer, spacer,
  plot_combined_temp, spacer, plot_combined_soil,
  spacer, spacer, spacer,
  # plot_combined_all, spacer, vis_ward_yield,
  plot_combined_all, spacer, lower_right_box,
  ncol = 3,  # Now 3 columns to accommodate spacers
  nrow = 5,  # Now 5 rows to accommodate spacers
  widths = c(1, 0.1, 1),    # Middle column is narrow spacer
  heights = c(1, 0.1, 1, 0.1, 1)  # Spacer rows
)
overview 

### 5) Save ----
ggsave(paste0(saving,"/overview_clusters.png"), plot= overview, width=15, height= 20)

##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# F) Visualization of grid cells ---- 
##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# This section generates Fig.A.2 in the appendix of the manuscript. 
# It uses the soil_oaf and the soil_50k dataset which have not been used in this script before. 

### 1) Prep data ----
# Drop duplicates in terms of location: 
soil_oaf <- soil_oaf %>% distinct(geometry, .keep_all = TRUE)
soil_50k <- soil_50k %>% distinct(geometry, keep_all = TRUE)

### 2) Lay grid network over the points----
# (Same procedure as for soil_2 in the beginning of this script)

# Before that: transform the projection 
points_proj_oaf <- st_transform(soil_oaf, crs = 32637)  # UTM Zone 37N (Kenya)
points_proj_50k <- st_transform(soil_50k, crs = 32637)  # UTM Zone 37N (Kenya)

# Lay the grids
grid_oaf <- st_make_grid(points_proj, cellsize = 2000, square = TRUE)
grid_50k <- st_make_grid(points_proj, cellsize = 2000, square = TRUE)

# Transform into sf object 
grid_sf_oaf <- st_sf(geometry = grid_oaf)
grid_sf_50k <- st_sf(geometry = grid_50k)

# Transform back to previous projection
grid_sf_oaf <- st_transform(grid_sf_oaf, crs = st_crs(soil_oaf))  # Back to original CRS
grid_sf_50k <- st_transform(grid_sf_50k, crs = st_crs(soil_50k))  # Back to original CRS

# Add a grid cell ID
grid_sf_oaf$grid_id <- 1:nrow(grid_sf_oaf)  # Create a unique ID for each grid cell
grid_sf_50k$grid_id <- 1:nrow(grid_sf_50k)  # Create a unique ID for each grid cell

# Join points to grid cells (this gives each point a grid cell)
points_joined_oaf <- st_join(soil_oaf, grid_sf_oaf, join = st_within)
points_joined_50k <- st_join(soil_50k, grid_sf_50k, join = st_within)

# Count the points within each grid cell
grid_with_points_oaf <- points_joined_oaf %>%
  group_by(grid_id) %>%  # Group by grid cell ID
  summarise(point_count = n()) %>%  # Count the number of points in each grid cell
  filter(point_count > 0)  # Keep only grid cells with at least one point
grid_with_points_50k <- points_joined_50k %>%
  group_by(grid_id) %>%  # Group by grid cell ID
  summarise(point_count = n()) %>%  # Count the number of points in each grid cell
  filter(point_count > 0)  # Keep only grid cells with at least one point

# Join the grid cells with points count back to the original grid_sf and filter for at least 1 point per grid cell
grid_with_points_geom_oaf <- grid_sf_oaf %>%
  filter(grid_id %in% grid_with_points_oaf$grid_id) %>% 
  left_join(st_drop_geometry(grid_with_points_oaf), by = "grid_id") %>%  # Use left_join instead of st_join
  mutate(point_bin = cut(point_count,
                         breaks = c(1, 5, 10, 15, 20, 25, Inf),  # Changed 1 to 0
                         labels = c("1–5", "5–10", "10–15", "15–20", "20–25", ">25"),
                         right = FALSE))  # Changed to right = FALSE
grid_with_points_geom_50k <- grid_sf_50k %>%
  filter(grid_id %in% grid_with_points_50k$grid_id) %>% 
  left_join(st_drop_geometry(grid_with_points_50k), by = "grid_id") %>%  # Use left_join instead of st_join
  mutate(point_bin = cut(point_count,
                         breaks = c(1, 5, 10, 15, 20, 25, Inf),  # Changed 1 to 0
                         labels = c("1–5", "5–10", "10–15", "15–20", "20–25", ">25"),
                         right = FALSE))  # Changed to right = FALSE

# Check for NA values
sum(is.na(grid_with_points_geom_oaf$point_bin))
sum(is.na(grid_with_points_geom_50k$point_bin))

# Sum stats on points per grid cell: 
summary(grid_with_points$point_count) #based on both (oaf and 50k)
summary(grid_with_points_oaf$point_count) #based on oaf
summary(grid_with_points_50k$point_count) #based on 50k


### 3) Build separate figures ----
# custom_colors <- c("yellow", "orange", "red", "purple", "blue", "black")
custom_colors <- c("#FFFF00", "#FFCC00", "#FF9900", "#FF6600", "#CC3300", "#800000")

# Plot by count category
# Both (oaf and 50k) together first: 
fig_both_withlegend <- ggplot() +
  geom_sf(data = grid_with_points_geom, aes(fill = point_bin), color = "NA", size = 0.2) +
  scale_fill_manual(values = custom_colors, name = "Point Count") +
  geom_sf(data = admin1_clipped, fill = NA, color = "black") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0, margin = margin(0, 0, 10, 0)))+
  # theme(legend.position = "none")+
  labs(title = "c) Coloured by total number of data extraction\nlocations per grid cell")
fig_both <- fig_both_withlegend + theme(legend.position = "none")
# Then oaf: 
fig_oaf <- ggplot() +
  geom_sf(data = grid_with_points_geom_oaf, aes(fill = point_bin), color = "NA", size = 0.2) +
  scale_fill_manual(values = custom_colors, name = "Point Count") +
  geom_sf(data = admin1_clipped, fill = NA, color = "black") +
  theme_minimal() +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0, margin = margin(0, 0, 10, 0)))+
  labs(title = "a) Coloured by number of crop cut locations\nper grid cell")
# Then 50k:
fig_50k <- ggplot() +
  geom_sf(data = grid_with_points_geom_50k, aes(fill = point_bin), color = "NA", size = 0.2) +
  scale_fill_manual(values = custom_colors, name = "Point Count") +
  geom_sf(data = admin1_clipped, fill = NA, color = "black") +
  theme_minimal() +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0, margin = margin(0, 0, 10, 0)))+
  labs(title = "b) Coloured by number of randomly placed locations\non maize cropland per grid cell")

### 4) Layout and combination of graphs ----
# Extract legend from one graph: 
legend <- cowplot::get_legend(fig_both_withlegend)

# Define heading 
heading <- "Distribution of data extraction locations across study region"
# Point of graphic is to show that a single data source would have not allowed full 
# spatial coverage but jointly (oaf and 50k points) the ROI is mostly covered

# Create title grob
title_grob <- textGrob(heading, gp = gpar(fontsize = 14, fontface = "bold"))

# Arrange everything
gridcells <- grid.arrange(
  title_grob,                    # Title on top
  fig_oaf, fig_50k, fig_both, legend,  # Plots and legend in one row
  ncol = 4, nrow = 2,           # 4 columns, 2 rows
  widths = c(1, 1, 1, 0.3),     # Relative widths for the plots row
  heights = c(0.1, 1),          # Small height for title, full for plots
  layout_matrix = rbind(c(1, 1, 1, 1),    # Title spans all columns
                        c(2, 3, 4, 5))     # Plots in second row
)

### 5) Save graph ----
ggsave(paste0(saving,"/gridcells.png"), plot= gridcells, width=16, height= 7)

##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# G) Descriptives on obs per cluster ----
##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
oaf %>%
  st_drop_geometry() %>% 
  summarise(across(starts_with("c."), list(
    n_non_na = ~ sum(!is.na(.)),
    n_unique = ~ n_distinct(., na.rm = TRUE)
  ))) %>%
  pivot_longer(everything(), names_to = "name", values_to = "value") %>%
  separate(name, into = c("cluster_variable", "stat"), sep = "_(?=n_)") %>%
  pivot_wider(names_from = stat, values_from = value)

oaf %>%
  st_drop_geometry() %>% 
  select(starts_with("c."), yield) %>%
  pivot_longer(cols = starts_with("c."), names_to = "index", values_to = "cluster") %>%
  mutate(cluster = as.integer(cluster)) %>%
  count(index, cluster) %>%
  complete(index, cluster = 1:14, fill = list(n = 0)) %>%
  pivot_wider(names_from = cluster, values_from = n, names_sort = TRUE) %>%
  arrange(index)

