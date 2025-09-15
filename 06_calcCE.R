#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Author: Ella Kirchner, Elinor Benami, Andrew Hobbs
# Project name: RS Get in the Zone
# Date Last Updated: # Fri Sep 12 10:01:49 2025
# Purpose: Simulate various insurance schemes and assess and visualize results
# Input:  all to be loaded from the RProj
          # - OAF data 
          # - Admin data on ward level (clipped) 
          # - Insurance zone assignments on ward level (one .RData for every agri-environmental indicator)
          # - Insurance zone assignments using lowertail clustering 
# Output Files: 
          # - Fig. 4 (Static comparison)
          # - Fig. 5 (Heatmaps)
          # - Fig. 6 (Slices of the heatmap) 
          # - Fig. 7 (Variance)
          # - Fig. A.7 in appendix (Dynamic zone count results, ignoring design risk)
          # - Fig. A.8 in appendix (Dynamic zone count results based on lowertail clustering)
          # - Fig. A.9 in appendix (Variance development across insurance zones of all agri-environmental inputs)
          # - Fig. A.10 in appendix (Yield boxplots) 
          # - Fig. A.11 in appendix (Dot-Whiskerplots)
# ReadMe: 1) Adapt file path to work in R-project, then loading of datasets should work 
#         2) Adapt file path for saving (line 38)
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# %%%%%%%%%%%%%%%%%%%%%%%%%%%
# Set Up 
# %%%%%%%%%%%%%%%%%%%%%%%%%%%

# 1. Clean out workspace
rm(list = ls(all = TRUE))

# 2. Set Up File path & directory: 
file_path<- "put_your_file_path_here" # File path to project folder

#setwd(file_path)
getwd()

overleaf_path <- "put_your_overleaf_path_or_other_saving_location_here"

# 3. Set Seed for replication 
set.seed(123456789)

# 4. Set date 
date_today <- as.character(Sys.Date())


# %%%%%%%%%%%%%%%%%%%%%%%%%%%
# Load libraries 
# %%%%%%%%%%%%%%%%%%%%%%%%%%%
library(dplyr)
library(sp)
library(sf) 
library(terra)
library(tidyverse)
library(cowplot)
library(dtplyr)
library(janitor)
library(lubridate)
library(tidyr)
library(tidylog)
library(summarytools)
library(PerformanceAnalytics)
library(RColorBrewer)
library(palmerpenguins)
library(gridExtra)
library(grid)
library(ggpubr)
library(ncdf4)
library(chron)
library(scales)
library(patchwork)
library(furrr)

source("ggsave_latex.R")

# %%%%%%%%%%%%%%%%%%%%%%%%%%%
# General visualization settings 
# %%%%%%%%%%%%%%%%%%%%%%%%%%%
clean_chart_clutter_explore <- 
  theme(    
    panel.grid.major = element_blank(),      # Remove panel grid lines
    panel.grid.minor = element_blank(),      # Remove panel grid lines
    panel.background = element_blank(),      # Remove panel background
    axis.line = element_line(colour = "grey"),       # Add axis line
    axis.title.y = element_text(size = 14, face = "bold", angle = 0, vjust = 0.5), # Rotate y axis so don't have to crank head
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.text.x = element_text(size = 14, face = "bold"),
    axis.text.y = element_text(size = 14, face = "bold"),
    plot.title = element_text(size = 15, face = "bold"),
    plot.subtitle = element_text(size = 14),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14, face = "bold"),
    legend.position = "bottom"
  )

# %%%%%%%%%%%%%%%%%%%%%%%%%%%
# 1 Load Data ---------
# %%%%%%%%%%%%%%%%%%%%%%%%%%%

# Admin data: 
load("datasets/admin3_clipped.RData")

# OAF data anonymized (=without geolocation): 
load("datasets/oaf_anonym.Rdata")
oaf <- oaf_anonym

# Dynamic clusters: 
load("datasets/dynamic_chirps.RData")
load("datasets/dynamic_temp.RData")
load("datasets/dynamic_ndvi.RData")
load("datasets/dynamic_evi.RData")
load("datasets/dynamic_soil.RData")
load("datasets/dynamic_all.RData")


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 2 Yield Data Prep ----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Winsorize maximum yield (roughly 8t is the maximum realistic yield):
quantiles <- quantile(oaf$yield, c(0.01, 0.05, 0.1, 0.25, 0.5,0.75,0.9,0.95, 0.99))
quantiles # 99% equals roughly 8t, hence, we use 99% 
oaf_sorted <- arrange(oaf, yield)
summary(oaf$yield)
oaf$yield <- pmin(oaf$yield, quantile(oaf$yield, 0.99))
# Winsorize minimum yield to same percentile (1%):
oaf$yield <- pmax(oaf$yield, quantile(oaf$yield, 0.01))
summary(oaf$yield)


# Filter to keep only name_concat values that appear in all years
all_years <- unique(oaf$year)

## This operation filters out observations, and in the process, some zones. 
oaf <- oaf %>%
  group_by(name_concat) %>%
  filter(length(unique(year)) == length(all_years)) %>% # drops 7k observations 
  ungroup()

oaf$name_concat2 <- paste0(oaf$name_1, "_", oaf$name_2)

length(unique(oaf$name_concat2)) # 67 subcounties after filtering 
length(unique(oaf$name_concat))  # 184 wards after filtering

# %%%%%%%%%%%%%%%%%%%%%%%%%%%
# 3 Set Params & fcts ----- 
# %%%%%%%%%%%%%%%%%%%%%%%%%%%

# 3.1 Parameters ----
price_per_kg = 40
kes_per_usd =  105
price_kg_usd = price_per_kg/kes_per_usd
price_t_usd=price_kg_usd*1000

# Depending on this value, the trigger changes (always in relation to zone avg yield: mean * trigger_val)
trigger_val = 1   
# Markup assumed to be zero for the analysis:
markup = 0.0
off_farm = (3000/kes_per_usd) # one median monthly income as off-farm income

n_wards <- length(unique(oaf$name_concat))
n_wards

total_people <- n_wards * 14 #gets to approximately 2.5k policyholders per year 
num_cc <- 32
minimum_iwi = 1 # used in case the income with insurance turns negative for any of the used indicators or scenarios

n_groups <- 14 # for current zoning
people_in_group <- total_people/n_groups 

rho = 1.5
maxyield = 8

# 3.2 Basic functions ----
# Utility function
utility <- function(c, rho) if(rho == 1) log(c) else (c ^ (1 - rho)) / (1 - rho)

# Certainty_equivalent
cert_equiv = function(expected_utility, rho) {
  if(rho == 1){
    exp(expected_utility)} else{
      ((1 - rho) * expected_utility) ^ (1 / (1 - rho))
    }
}

# 3.3 Functions to determine ins value ----

# Function on how to determine zone averages 
sample_mean <- function(vec, n){
  if(n<length(vec)){
    samp <- sample(vec, size=n)
  } else {
    samp <- vec
  }
  return(mean(samp))
}

# Function to generate income
gendf_income_wwo_ins <- function(data, grouping, num_cc=NULL, total_cuts=NULL, design_risk=TRUE){
  # zone
  #print(grouping)# to check if the grouping works: it works!
  n_groups <- data %>% pull(!!grouping) %>% n_distinct()
  people_in_group <- total_people/n_groups
  
  if(!is.null(total_cuts)){
    # calculate the number of crop cuts from the number of zones, and override num_cc with that
    num_cc <- floor(total_cuts/n_groups)
  } 
  
  if(!is.null(num_cc) && design_risk){ # specific costs, design risk
    meanfunc <- function(x) sample_mean(x, n=num_cc)
    zone_cost <- num_cc * cost_cc
  } else if (!is.null(num_cc) && !design_risk){ # specific costs, no design risk
    meanfunc <- function(x) mean(x)
    zone_cost <- num_cc * cost_cc
  }  else {
    meanfunc <- function(x) mean(x)
  }
  
  print(paste('People', people_in_group))
  print(paste('Zone cost', zone_cost))
  
  data  %>%
    mutate(y_it = yield, 
           c_it = off_farm + (y_it * price_t_usd)) %>%  
    group_by(!!grouping) %>% # items that vary only across all years
    mutate(
      mean_yield = meanfunc(y_it),
      trigger = mean_yield * trigger_val,
      group_name = !!grouping) %>%
    group_by(!!grouping, year) %>%
    mutate(y_gt = mean(y_it),
           payout = ifelse(y_gt < trigger,  # no grouping by year needed here 
                           (trigger - y_gt)*price_t_usd,0)) %>%     # avg yield per year
    ungroup() %>%
    mutate(
      premium = mean(payout) + mean(payout)*markup,
      c_it_ins = c_it + payout - (premium + (zone_cost/people_in_group)),
      
      #Translate calculated incomes into individual utility values
      util = utility(c_it, rho),
      util_ins = utility(c_it_ins, rho),
            # Transform the individual utility values into expected utility and use that expected utility to calc CE
      ce = cert_equiv(mean(util, na.rm = FALSE), rho),
      ce_ins = cert_equiv(mean(util_ins, na.rm = FALSE), rho),
            # Take the difference of the CE values 
      diff_ce= ce_ins - ce)
}

# Supporting functions to calc income
# to ensure quasi-quotation works appropriately
income_clust <- function(data, cluster_type, num_cc=NULL, total_cuts=NULL, design_risk=design_risk){
    cluster_type <- enquo(cluster_type)
    
    out <-
      gendf_income_wwo_ins(data, cluster_type, num_cc=num_cc, total_cuts=total_cuts, design_risk=design_risk) %>%
      dplyr::select(name_concat, unique_id, year, y_it, c_it, mean_yield, trigger, payout, premium, c_it_ins, group_name, util, util_ins, ce, ce_ins, diff_ce)%>%
      st_drop_geometry()
  }

# Modified merge function to join individual data with ward-level clusters
merge_yields_convert_income <- function(clustering_data, cluster_prefix = "cluster", num_cc=NULL, total_cuts=NULL, design_risk=TRUE){
  merged <-
    oaf %>%
    dplyr::select(-contains("c.")) %>% #remove previous clusters from the data, start fresh! 
    dplyr::left_join(clustering_data, by = "name_concat") %>%  # Join individual data with ward-level clusters
    dplyr::select(name_concat, unique_id, year, yield, contains(cluster_prefix))%>%
    st_drop_geometry()
  
  out <-
    merged %>%
    dplyr::mutate(across(starts_with(cluster_prefix), ~income_clust(merged, .x, num_cc=num_cc, total_cuts=total_cuts, design_risk=design_risk)))
}

# # 3.4 Func to generate plot ----
# # Function generate a plot for one specific scenario of assumed cost and # of crop cuts (but variable number of insurance zones) which can later be combined as required (e.g. for overviews in the appendix etc.)
get_cluster_diff_ce <- function(data_input) {
  do.call(rbind, lapply(5:100, function(i) {
    cluster_name <- paste0("cluster_", i)
    if (cluster_name %in% names(data_input)) {
      data.frame(cluster = i, diff_ce = data_input[[cluster_name]]$diff_ce[1])
    }
  }))
}

chart_ce_cluster <- function(cluster_data, indicator, lower, upper) {
  ggplot(cluster_data, aes(x = cluster, y = diff_ce)) +
    # Invisible point to anchor 0,0 without affecting lines/smoother
    geom_point(data = tibble(cluster = 0, diff_ce = 0), aes(x = cluster, y = diff_ce), alpha = 0) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray30", linewidth = 0.5) +  
    # geom_smooth() +
    geom_line() +
    scale_y_continuous(limits = c(lower, upper), expand = c(0, 0), 
                       breaks = seq(lower, upper, by = 2),
                       minor_breaks = seq(lower, upper, by = 1)
    ) +
    scale_x_continuous(limits = c(0, 105),  # <- force extension beyond last data point
                       breaks = seq(0, 105, by = 10),        
                       minor_breaks = seq(0, 105, by = 5),  
                       expand = c(0, 0)                       ) +
    labs(title = paste0(indicator), #, "-based Clusters
         x = "\nNumber of Insurance Zones",
         y = "USD") +
    clean_chart_clutter_explore +
    guides(
      x = guide_axis(minor.ticks = TRUE),
      y = guide_axis(minor.ticks = TRUE)
    )
}

# %%%%%%%%%%%%%%%%%%%%%%%%%%%
# 4 Application of Functions ----
# %%%%%%%%%%%%%%%%%%%%%%%%%%%
# Per scenario (combi of cost and # of crop cuts), we calculate the insurance value 
# and the corresponding graph as each scenario requires to redefine the parameters

# 4.1 Zero cost scenario ----

# Set parameters 
cost_cc = 0
num_cc = 32 
zone_cost = cost_cc*num_cc

# Calculate incomes
ndvi_income_varclust_0 <- merge_yields_convert_income(dynamic_ndvi, "cluster")
evi_income_varclust_0 <- merge_yields_convert_income(dynamic_evi, "cluster")
temp_income_varclust_0 <- merge_yields_convert_income(dynamic_temp, "cluster")
soil_income_varclust_0 <- merge_yields_convert_income(dynamic_soil, "cluster")
chirp_income_varclust_0 <- merge_yields_convert_income(dynamic_chirps, "cluster")
combined_income_varclust_0 <- merge_yields_convert_income(dynamic_all, "cluster")

# Create graph (lower, upper refer to y-axis limits)
lower = -6
upper = 12
ndvigraph_0 <- chart_ce_cluster(get_cluster_diff_ce(ndvi_income_varclust_0), "NDVI", lower, upper)
chirpgraph_0 <- chart_ce_cluster(get_cluster_diff_ce(chirp_income_varclust_0), "Seasonal Precipitation", lower, upper)
evigraph_0 <- chart_ce_cluster(get_cluster_diff_ce(evi_income_varclust_0), "EVI", lower, upper)
soilgraph_0 <- chart_ce_cluster(get_cluster_diff_ce(soil_income_varclust_0), "Soil", lower, upper)
tempgraph_0 <- chart_ce_cluster(get_cluster_diff_ce(temp_income_varclust_0), "Temperature", lower, upper)
combinedgraph_0 <- chart_ce_cluster(get_cluster_diff_ce(combined_income_varclust_0), "Combined indicators", lower, upper)


# 4.2 Medium cost scenario ----
# Set parameters 
cost_cc= 10
num_cc= 32 
zone_cost = cost_cc*num_cc

# Calculate incomes 
ndvi_income_varclust_10 <- merge_yields_convert_income(dynamic_ndvi)
chirp_income_varclust_10 <- merge_yields_convert_income(dynamic_chirps)
evi_income_varclust_10 <- merge_yields_convert_income(dynamic_evi)
soil_income_varclust_10 <- merge_yields_convert_income(dynamic_soil)
temp_income_varclust_10 <- merge_yields_convert_income(dynamic_temp)
combined_income_varclust_10 <- merge_yields_convert_income(dynamic_all)

# Create graph 
# lower= -10 
# upper = 10
ndvigraph_10 <- chart_ce_cluster(get_cluster_diff_ce(ndvi_income_varclust_10), "NDVI", lower, upper)
chirpgraph_10 <- chart_ce_cluster(get_cluster_diff_ce(chirp_income_varclust_10), "Seasonal Precipitation", lower, upper)
evigraph_10 <- chart_ce_cluster(get_cluster_diff_ce(evi_income_varclust_10), "EVI", lower, upper)
soilgraph_10 <- chart_ce_cluster(get_cluster_diff_ce(soil_income_varclust_10), "Soil", lower, upper)
tempgraph_10 <- chart_ce_cluster(get_cluster_diff_ce(temp_income_varclust_10), "Temperature", lower, upper)
combinedgraph_10 <- chart_ce_cluster(get_cluster_diff_ce(combined_income_varclust_10), "Combined indicators", lower, upper)


# 4.3 High cost scenario ----

# Set parameters 
cost_cc= 36
num_cc= 32 
zone_cost = cost_cc*num_cc

# Calculate incomes 
ndvi_income_varclust_36 <- merge_yields_convert_income(dynamic_ndvi)
chirp_income_varclust_36 <- merge_yields_convert_income(dynamic_chirps)
evi_income_varclust_36 <- merge_yields_convert_income(dynamic_evi)
soil_income_varclust_36 <- merge_yields_convert_income(dynamic_soil)
temp_income_varclust_36 <- merge_yields_convert_income(dynamic_temp)
combined_income_varclust_36 <- merge_yields_convert_income(dynamic_all)

# Create graph 
lower= -10
# upper = 20
ndvigraph_36 <- chart_ce_cluster(get_cluster_diff_ce(ndvi_income_varclust_36), "NDVI", lower, upper)
chirpgraph_36 <- chart_ce_cluster(get_cluster_diff_ce(chirp_income_varclust_36), "Seasonal Precipitation", lower, upper)
evigraph_36 <- chart_ce_cluster(get_cluster_diff_ce(evi_income_varclust_36), "EVI", lower, upper)
soilgraph_36 <- chart_ce_cluster(get_cluster_diff_ce(soil_income_varclust_36), "Soil", lower, upper)
tempgraph_36 <- chart_ce_cluster(get_cluster_diff_ce(temp_income_varclust_36), "Temperature", lower, upper)
combinedgraph_36 <- chart_ce_cluster(get_cluster_diff_ce(combined_income_varclust_36), "Combined indicators", lower, upper)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 5  Computing static results (Fig. 4) ----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# This section generates Figure 4 in the manuscript.

# 5.1 Extract values for static comparison ----
# We extract the values from the previous "dynamic" calculation. 

# 5.1.1 For agri-environmental inputs:
# Get the actual dataset names:
datasets <- list(
  ndvi_0 = ndvi_income_varclust_0,
  evi_0 = evi_income_varclust_0,
  temp_0 = temp_income_varclust_0,
  soil_0 = soil_income_varclust_0,
  all_0 = combined_income_varclust_0,
  chirps_0 = chirp_income_varclust_0,
  ndvi_10 = ndvi_income_varclust_10,
  evi_10 = evi_income_varclust_10,
  temp_10 = temp_income_varclust_10,
  soil_10 = soil_income_varclust_10,
  all_10 = combined_income_varclust_10,
  chirps_10 = chirp_income_varclust_10,
  ndvi_36 = ndvi_income_varclust_36,
  evi_36 = evi_income_varclust_36,
  temp_36 = temp_income_varclust_36,
  soil_36 = soil_income_varclust_36,
  all_36 = combined_income_varclust_36,
  chirps_36 = chirp_income_varclust_36
)

# Function to extract dataframe to build static graph: 
extract_specific_diff_ce <- function(datasets, cluster_number) {
  # Initialize vectors to store results
  dataset_names <- character()
  suffixes <- character()
  diff_ce_values <- numeric()
  
  # Create cluster name
  cluster_name <- paste0("cluster_", cluster_number)
  
  # Loop through each dataset
  for (dataset_key in names(datasets)) {
    # Split dataset name to get base name and suffix
    parts <- strsplit(dataset_key, "_")[[1]]
    base_name <- parts[1]
    suffix <- parts[2]
    
    # Get cluster data
    cluster_data <- datasets[[dataset_key]][[cluster_name]]
    diff_ce_value <- cluster_data$diff_ce[1]
    
    # Store the values
    dataset_names <- c(dataset_names, base_name)
    suffixes <- c(suffixes, suffix)
    diff_ce_values <- c(diff_ce_values, diff_ce_value)
  }
  
  # Create the result dataframe in long format
  result <- data.frame(
    dataset = dataset_names,
    cost = suffixes,
    diff_ce = diff_ce_values,
    stringsAsFactors = FALSE
  )
  
  return(result)
}


# Extract data for different cluster numbers (admin1 = 14, admin2 = 100)
diff_ce_14 <- extract_specific_diff_ce(datasets, 14)
diff_ce_100<- extract_specific_diff_ce(datasets, 100)

# 5.1.2 For admin and yield insurance zones 
# Calculate values for admin1, admin2, and yields separately as they are not included in the income_varcluster dataframes

# Function for that (based on function from dynamic analysis)
calculate_oaf_diff_ce <- function(oaf, grouping_variable, cost_scenarios) {
  # Initialize results dataframe
  results <- data.frame(
    dataset = character(),
    cost = character(),
    diff_ce = numeric(),
    stringsAsFactors = FALSE
  )
  
  # Calculate basic parameters
  n_groups <- oaf[[grouping_variable]] %>% n_distinct()
  people_in_group <- total_people / n_groups
  num_cc <- 32
  
  # Loop through each cost scenario
  for (cost_cc in cost_scenarios) {
    zone_cost <- num_cc * cost_cc
    
    # Calculate diff_ce for this cost scenario
    diff_ce_result <- oaf %>%
      mutate(y_it = yield, 
             c_it = off_farm + (y_it * price_t_usd)) %>%  
      group_by(!!sym(grouping_variable), year) %>%
      mutate(y_gt = mean(y_it)) %>%     # avg yield per year
      group_by(!!sym(grouping_variable)) %>% # items that vary only across all years
      mutate(
        mean_yield = mean(y_gt),
        trigger = mean_yield * trigger_val,
        payout = ifelse(y_gt < trigger,  # no grouping by year needed here 
                        (trigger - y_gt)*price_t_usd, 0),
        
        group_name = !!sym(grouping_variable)) %>%
      ungroup() %>%
      mutate(
        premium = mean(payout) + mean(payout)*markup,
        c_it_ins = c_it + payout - (premium + (zone_cost/people_in_group)),
        
        #Translate calculated incomes into individual utility values
        util = utility(c_it, rho),
        util_ins = utility(c_it_ins, rho),
        
        # Transform the individual utility values into expected utility and use that expected utility to calc CE
        ce = cert_equiv(mean(util, na.rm = FALSE), rho),
        ce_ins = cert_equiv(mean(util_ins, na.rm = FALSE), rho),
        
        # Take the difference of the CE values 
        diff_ce = ce_ins - ce) %>%
      
      # Get the unique diff_ce value for this scenario
      pull(diff_ce) %>%
      first()
   
    # Add to results
    results <- rbind(results, data.frame(
      dataset = grouping_variable,
      cost = as.character(cost_cc),
      diff_ce = diff_ce_result,
      stringsAsFactors = FALSE
    ))
  }
  
  return(results)
}

# Calculate diff_ce for all three variables and cost scenarios
cost_scenarios <- c(0, 10, 36)
variables <- c("name_1", "name_2", "c.yield", "c67.yield")

# Initialize combined results
oaf_results <- data.frame(
  dataset = character(),
  suffix = character(),
  diff_ce = numeric(),
  stringsAsFactors = FALSE
)

# Calculate for each variable
for (var in variables) {
  var_results <- calculate_oaf_diff_ce(oaf, var, cost_scenarios)
  oaf_results <- rbind(oaf_results, var_results)
}

# 5.1.3 Combine for static comparison 
# Combine with the previously generated dataset on other indicators
diff_ce_14 <- rbind(diff_ce_14, oaf_results)
diff_ce_100 <- rbind(diff_ce_100, oaf_results)

# 5.2 Design static barplot (Fig 4) ----

# 5.2.1 Rename & Reorder variables 

# Rename the values in the datasets as we want to see them in the graph: 
var_labels <- c(
  name_1 = "Admin 1",
  name_2 = "Admin 2",
  chirps  = "Seasonal\nPrecip",
  evi    = "EVI",
  ndvi   = "NDVI",
  soil   = "Soil",
  temp   = "Temp",
  all    = "All\nindicators"
)

# Rename datasets (and filter to keep only the assessment for the relevant number of insurance zones) 
diff_ce_admin1 <- diff_ce_14 %>%
  mutate(dataset = recode(dataset, !!!var_labels)) %>%
  mutate(dataset = recode(dataset, "c.yield" = "Yield\n(Oracle)"))%>%
  filter(dataset != "c67.yield")

diff_ce_admin2 <- diff_ce_100 %>%
  mutate(dataset = recode(dataset, !!!var_labels)) %>%
  mutate(dataset = recode(dataset, "c67.yield" = "Yield\n(Oracle)")) %>%
  filter(dataset != "c.yield")

# Fixed inputs
predefined_inputs_admin1 <- c("Yield\n(Oracle)", "Admin 1")
predefined_inputs_admin2 <- c("Yield\n(Oracle)", "Admin 2")

# Create factor order
custom_order_admin1 <- c(
  predefined_inputs_admin1,
  diff_ce_admin1 %>%
    group_by(dataset) %>%
    summarize(val = mean(diff_ce, na.rm = FALSE), .groups = "drop") %>%
    filter(!dataset %in% predefined_inputs_admin1) %>%
    arrange(desc(val)) %>%
    pull(dataset)
)
custom_order_admin2 <- c(
  predefined_inputs_admin2,
  diff_ce_admin2 %>%
    group_by(dataset) %>%
    summarize(val = mean(diff_ce, na.rm = FALSE), .groups = "drop") %>%
    filter(!dataset %in% predefined_inputs_admin2) %>%
    arrange(desc(val)) %>%
    pull(dataset)
)

# Apply custom orders: 
diff_ce_admin1 <- diff_ce_admin1 %>%
  mutate(dataset = factor(dataset, levels = custom_order_admin1))
diff_ce_admin2 <- diff_ce_admin2 %>%
  mutate(dataset = factor(dataset, levels = custom_order_admin2))

# 5.2.2 Design graph 

# Determine the position index of "Yield\n(Oracle)" in custom_order
yield_position <- which(custom_order_admin1 == "Yield\n(Oracle)")

# Add 0.5 to position the line between yield and the next variable
vline_pos <- yield_position + 1.5

# Visualization for paper
diff_ce_barplot_paper_admin1 <-
  diff_ce_admin1 %>%
  filter(dataset != "Admin 2") %>%
  ggplot(aes(fill = factor(cost), y = diff_ce, x = dataset)) +
  geom_bar(position = "dodge", stat = "identity") +
  geom_vline(xintercept = vline_pos, linetype = "dashed", color = "black", linewidth = 0.6) +
  annotate("text", x = yield_position + 0.5, y = 11.9, label = "Benchmark", size = 4.5, fontface = "italic") +
  annotate("text", x = yield_position + 4.5, y = 11.9, label = "Agri-environmental inputs",
           size = 4.5, fontface = "italic") +
  labs(y = "USD ",
       x = "\nInput for Defining Zones",
       fill = "Costs per cropcut (USD)",
       title = "Risk Reduction Value of Insurance Relative to No Insurance",
       subtitle ="Assuming 14 insurance zones with 32 crop cuts per zone and 3 cost scenarios") +
  scale_y_continuous(
    breaks = seq(-4, 12, by = 2),
    minor_breaks = seq(-4, 12, by = 1),
    limits = c(-4, 12.5),
  ) +
  scale_fill_brewer(palette= "PuBu", direction = -1) +
  clean_chart_clutter_explore +
  guides(
    y = guide_axis(minor.ticks = TRUE)
  )

diff_ce_barplot_paper_admin2 <-
  diff_ce_admin2 %>%
  filter(dataset != "Admin 1") %>%
  ggplot(aes(fill = factor(cost), y = diff_ce, x = dataset)) +
  geom_bar(position = "dodge", stat = "identity") +
  geom_vline(xintercept = vline_pos, linetype = "dashed", color = "black", linewidth = 0.6) +
  annotate("text", x = yield_position + 0.5, y = 18, label = "Benchmark", size = 4.5, fontface = "italic") +
  annotate("text", x = yield_position + 4.5, y = 18, label = "Agri-environmental inputs",
           size = 4.5, fontface = "italic") +
  labs(y = "USD ",
       x = "\nInput for Defining Zones",
       fill = "Costs per cropcut (USD)",
       title = "Risk Reduction Value of Insurance Relative to No Insurance",
       subtitle ="Assuming 100 insurance zones with 32 crop cuts per zone and 3 cost scenarios") +
  scale_y_continuous(
    breaks = seq(-30, 20, by = 10),
    minor_breaks = seq(-30, 20, by = 5),
    limits = c(-33, 20),
  ) +
  scale_fill_brewer(palette= "PuBu", direction = -1) +
  clean_chart_clutter_explore +
  guides(
    y = guide_axis(minor.ticks = TRUE)
  )

## Combine the two static charts together with a shared legend
# Suppress individual legends
p1 <- diff_ce_barplot_paper_admin1 + theme(legend.position = "none", axis.title.x = element_blank())
p2 <- diff_ce_barplot_paper_admin2 + theme(legend.position = "none",
                                           plot.title = element_blank())

# Create legend using one of the plots
legend <- get_legend(diff_ce_barplot_paper_admin1) 
final_plot <- (p1 / plot_spacer() / p2) / legend +
  plot_layout(heights = c(1.1, 0.1, 1.2, 0.1))  # Adjust spacer height as needed
final_plot

# 5.3 Save Fig. 4 ----
ggsave.latex(final_plot,
             filename = overleaf_path(paste0("/figs/diff_ce_barplot_admin1_admin2_comparison_", date_today, ".png")),
             label = "fig:diff_ce_barplot_admin1_admin2_comparison", 
             width = 13, height = 9, unit = "in"
)

# %%%%%%%%%%%%%%%%%%%%%%%%%%%
# 6 Dynamic results plots (Fig. A.7) ----
# %%%%%%%%%%%%%%%%%%%%%%%%%%%
# This section generates the graphic presented in Fig.A.7 in the appendix to the manuscript. 

# 6.1 Data prep ----
# Extract all the highlight points first  
n_admin1 = n_distinct(admin3_clipped$name_1) # 14
n_admin2 = n_distinct(admin3_clipped$name_2) #100
a1_zero <- diff_ce_admin1$diff_ce[diff_ce_admin1$dataset == "Admin 1" & diff_ce_admin1$cost == 0][[1]]
a1_ten <- diff_ce_admin1$diff_ce[diff_ce_admin1$dataset == "Admin 1" & diff_ce_admin1$cost == 10][[1]]
a1_36 <- diff_ce_admin1$diff_ce[diff_ce_admin1$dataset == "Admin 1" & diff_ce_admin1$cost == 36][[1]]
a2_zero <- diff_ce_admin2$diff_ce[diff_ce_admin2$dataset == "Admin 2" & diff_ce_admin2$cost == 0][[1]]
a2_ten <- diff_ce_admin2$diff_ce[diff_ce_admin2$dataset == "Admin 2" & diff_ce_admin2$cost == 10][[1]]
a2_36 <- diff_ce_admin2$diff_ce[diff_ce_admin2$dataset == "Admin 2" & diff_ce_admin2$cost == 36][[1]]

# 6.2 Create cluster plot sets for varying inputs ----
make_cluster_plot_set <- function(suffix,
                                  cluster_label = suffix,
                                  cropcut_cost = 0,
                                  highlight_coords = list(x = c(n_admin1, n_admin2), y = c(a1_zero, a2_zero)),
                                  y_limits = c(-20, 20),
                                  return_graphs_only = FALSE) {
  
  # 1. Generate individual plots
  datasets <- list(
    NDVI = paste0("ndvi_income_varclust_", suffix),
    Precipitation = paste0("chirp_income_varclust_", suffix),
    EVI = paste0("evi_income_varclust_", suffix),
    Soil = paste0("soil_income_varclust_", suffix),
    Temperature = paste0("temp_income_varclust_", suffix),
    `Combined indicators` = paste0("combined_income_varclust_", suffix))
  
  graph_list <- lapply(names(datasets), function(name) {
    varname <- datasets[[name]]
    chart_ce_cluster(get_cluster_diff_ce(get(varname)), name, lower, upper)
  })
  names(graph_list) <- names(datasets)
  
  # 2. Highlight points and themes
  highlight_df <- tibble(x = highlight_coords$x, y = highlight_coords$y, point_type = c("Admin 1", "Admin 2"))
  highlight_points <- geom_point(data = highlight_df, aes(x, y, color = point_type), size = 3)
  highlight_colors <- scale_color_manual("Predefined Boundaries", values = c("Admin 1" = "purple", "Admin 2" = "orange"))
  common_scale <- list(
    scale_y_continuous(limits = y_limits, oob = scales::oob_squish),
    scale_x_continuous(breaks = seq(0, 100, by = 20)))
  remove_x_axis <- theme(axis.ticks.x = element_blank(), axis.title.x = element_blank())
  
  common_theme <- theme_minimal() + rremove("ylab") +
    theme(
      axis.text = element_text(size = 11),
      axis.title = element_text(size = 11),
      legend.text = element_text(size = 11),
      legend.title = element_text(size = 11, face = "bold"))
  
  # 3. Extract legend
  legend_plot <- graph_list[["NDVI"]] + common_theme + common_scale + highlight_points + highlight_colors + theme(legend.position = "bottom")
  legend_obj <- get_legend(legend_plot)
  
  # 4. Style individual plots
  graph_list_nolegend <- lapply(graph_list, function(p) {
    p +
      common_theme +
      common_scale +
      highlight_points +
      highlight_colors +
      geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
      remove_x_axis +
      theme(legend.position = "none")
  })
  
  if (return_graphs_only) {return(graph_list_nolegend)}
  
  plot_grid <- do.call(ggarrange, c(graph_list_nolegend, list(ncol = 2, nrow = 3)))
  
  # 5. Title and subtitle
  main_title <- text_grob("Risk Reduction Value by Number of Insurance Zones and Input Data", size = 13, face = "bold")
  subtitle <- text_grob(paste0("(USD ", cropcut_cost, "/cropcut, ", cluster_label, " clusters)"), size = 12, face = "plain")
  title_grob <- arrangeGrob(
    grobs = list(main_title, subtitle),
    ncol = 1,
    heights = unit(c(1.3, 1), "lines")
  )
  
  # 6. Final figure
  fig_export <- grid.arrange(
    title_grob,
    arrangeGrob(
      plot_grid,
      left = text_grob("USD", rot = 90, vjust = 1, size = 13, face = "bold"),
      bottom = text_grob("Number of Zones", size = 13, face = "bold")
    ),
    legend_obj,
    ncol = 1,
    heights = unit(c(1.8, 10, 1.2), "cm")
  )
  
  # 7. Save with ggsave.latex
  ggsave.latex(fig_export,
               filename = overleaf_path(paste0("/figs/variable_clusters_", suffix, "_", date_today, ".png")),
               label = paste0("fig:variable_clusters_paper_", suffix)
               )
  print(fig_export)  # <--- display the plot when the function runs
  invisible(fig_export)
}

# 6.3 Combine them together -----
make_combined_cluster_panel_from_specs <- function(plot_specs, suffix_labels, legend_source = "ndvigraph_0") {
  
  # 1. Build list of plots per spec
  all_plots <- lapply(plot_specs, function(spec) {
    make_cluster_plot_set(
      suffix = spec$suffix,
      cropcut_cost = spec$cropcut_cost,
      highlight_coords = list(x = c(14, n_admin2), y = spec$highlight_y),
      y_limits = spec$ylims,
      return_graphs_only = TRUE
    )
  })
  
  # 2. Add row group labels using text grobs
  row_titles <- lapply(suffix_labels, function(lbl) text_grob(lbl, rot = 90, face = "bold", size = 12))
  
  # 3. Combine each row label with its 6 plots
  combined_rows <- Map(function(label, plots) {
    c(list(label), plots)
  }, row_titles, all_plots)
  
  # 4. Flatten into full list for plotting
  plot_matrix <- do.call(c, combined_rows)
  
  # 5. Legend (using any NDVI plot with a highlight)
  first_spec <- plot_specs[[1]]
  legend_plot <- get(legend_source) +
    theme_minimal() +
    scale_y_continuous(limits = first_spec$ylims, oob = scales::oob_squish) +
    scale_x_continuous(breaks = seq(0, 100, by = 20)) +
    scale_color_manual("Predefined Boundaries", values = c("Admin 1" = "purple", "Admin 2" = "orange")) +
    geom_point(
      data = tibble(x = c(14, n_admin2), y = first_spec$highlight_y, point_type = c("Admin 1", "Admin 2")),
      aes(x, y, color = point_type), size = 3
    ) +
    theme(legend.position = "bottom")
  
  legend_only <- get_legend(legend_plot)
  
  # 6. Compose panel
  grid_plot <- wrap_plots(plot_matrix, ncol = 7, nrow = length(plot_specs),
                          widths = c(0.6, rep(1, 6)))
  
  combined_output <- (
    grid_plot / legend_only +
      plot_layout(heights = c(13, 1.5))
  ) +
    plot_annotation(
      title = "Risk reduction value by crop cut cost, clustering input, and number of zones",
      theme = theme(
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
        plot.margin = ggplot2::margin(2, 2, 2, 2, "mm")
      )
    )
  
  return(combined_output)
}

plot_specs <- list(
  list(suffix = "0",  cropcut_cost = 0,  ylims = c(-5, 15),  highlight_y = c(a1_zero, a2_zero)),
  list(suffix = "10", cropcut_cost = 10, ylims = c(-5, 15),  highlight_y = c(a1_ten, a2_ten)),
  list(suffix = "36", cropcut_cost = 36, ylims = c(-60, 15), highlight_y = c(a1_36, a2_36))
)

suffix_labels <- c("Low Cost\n(0USD/crop cut)", "Med Cost\n(10USD/crop cut)", "High Cost\n(36USD/crop cut)")

combined_dynamic_plot <- make_combined_cluster_panel_from_specs(plot_specs, suffix_labels)
combined_dynamic_plot

# 6.4 Save plot ----
ggsave.latex(combined_dynamic_plot,
             filename = overleaf_path(paste0("/figs/combined_cluster_panel_", date_today, ".png")),
             label = "fig:combined_cluster_panel",
             width = 13, height = 7, units = "in"
)

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 7 Dynamic results with lowertail clustering (Appendix A.6)----
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# This section generates the results presented and discussed in Appendix A.6. 
# Specifically it generates Fig.A.8. 

# 7.1 Functionalize data extraction for lowertail ---- 
# Helper function to load required datasets
load_dynamic_rdata_files <- function(suffix, env = .GlobalEnv) {
  file_stubs <- c("dynamic_chirps", "dynamic_temp", "dynamic_ndvi", "dynamic_evi", "dynamic_soil", "dynamic_all")
  paths <- paste0("datasets/", file_stubs, "_", suffix, ".RData")
  purrr::walk(paths, load, envir = env)
  loaded_objects <- paste0(file_stubs, "_", suffix)
  print(loaded_objects)
}

# 7.2 Lowertail clusters  ----

# Set parameters 
cost_cc = 0
lower = -5 
upper = 15
num_cc = 32
zone_cost = cost_cc*num_cc

# Load dataset with lowertail clusters 
load_dynamic_rdata_files("lowertail")

# Calculate incomes for all clustering variables
ndvi_income_varclust_lowertail <- merge_yields_convert_income(dynamic_ndvi_lowertail)
evi_income_varclust_lowertail <- merge_yields_convert_income(dynamic_evi_lowertail)
temp_income_varclust_lowertail <- merge_yields_convert_income(dynamic_temp_lowertail)
soil_income_varclust_lowertail <- merge_yields_convert_income(dynamic_soil_lowertail)
chirp_income_varclust_lowertail <- merge_yields_convert_income(dynamic_chirps_lowertail)
combined_income_varclust_lowertail <- merge_yields_convert_income(dynamic_all_lowertail)

# Generate sample graphics
ndvigraph_lowertail <- chart_ce_cluster(get_cluster_diff_ce(ndvi_income_varclust_lowertail), "NDVI", lower, upper)
chirpgraph_lowertail <- chart_ce_cluster(get_cluster_diff_ce(chirp_income_varclust_lowertail), "Precipitation", lower, upper)
evigraph_lowertail <- chart_ce_cluster(get_cluster_diff_ce(evi_income_varclust_lowertail), "EVI", lower, upper)
soilgraph_lowertail <- chart_ce_cluster(get_cluster_diff_ce(soil_income_varclust_lowertail), "Soil", lower, upper)
tempgraph_lowertail <- chart_ce_cluster(get_cluster_diff_ce(temp_income_varclust_lowertail), "Temperature", lower, upper)
combinedgraph_lowertail <- chart_ce_cluster(get_cluster_diff_ce(combined_income_varclust_lowertail), "Combined indicators", lower, upper)

make_cluster_plot_set(suffix = "lowertail")

# 7.3 Lowerlowertail clusters ----

# No change in parameters

# Load dataset with lowerlowertail clusters 
load_dynamic_rdata_files("lowerlowertail")

# Calculate incomes for all clustering variables
ndvi_income_varclust_lowerlowertail <- merge_yields_convert_income(dynamic_ndvi_lowerlowertail)
evi_income_varclust_lowerlowertail <- merge_yields_convert_income(dynamic_evi_lowerlowertail)
temp_income_varclust_lowerlowertail <- merge_yields_convert_income(dynamic_temp_lowerlowertail)
soil_income_varclust_lowerlowertail <- merge_yields_convert_income(dynamic_soil_lowerlowertail)
chirp_income_varclust_lowerlowertail <- merge_yields_convert_income(dynamic_chirps_lowerlowertail)
combined_income_varclust_lowerlowertail <- merge_yields_convert_income(dynamic_all_lowerlowertail)

# Generate sample graphics
ndvigraph_lowerlowertail <- chart_ce_cluster(get_cluster_diff_ce(ndvi_income_varclust_lowerlowertail), "NDVI", lower, upper)
chirpgraph_lowerlowertail <- chart_ce_cluster(get_cluster_diff_ce(chirp_income_varclust_lowerlowertail), "Precipitation", lower, upper)
evigraph_lowerlowertail <- chart_ce_cluster(get_cluster_diff_ce(evi_income_varclust_lowerlowertail), "EVI", lower, upper)
soilgraph_lowerlowertail <- chart_ce_cluster(get_cluster_diff_ce(soil_income_varclust_lowerlowertail), "Soil", lower, upper)
tempgraph_lowerlowertail <- chart_ce_cluster(get_cluster_diff_ce(temp_income_varclust_lowerlowertail), "Temperature", lower, upper)
combinedgraph_lowerlowertail <- chart_ce_cluster(get_cluster_diff_ce(combined_income_varclust_lowerlowertail), "Combined indicators", lower, upper)

make_cluster_plot_set(suffix = "lowerlowertail")

# 7.4 log clusters  -----

# No change in parameters

# Load dataset with log clusters 
load_dynamic_rdata_files("log")

# Calculate incomes for all clustering variables
ndvi_income_varclust_log <- merge_yields_convert_income(dynamic_ndvi_log)
evi_income_varclust_log <- merge_yields_convert_income(dynamic_evi_log)
temp_income_varclust_log <- merge_yields_convert_income(dynamic_temp_log)
soil_income_varclust_log <- merge_yields_convert_income(dynamic_soil_log)
chirp_income_varclust_log <- merge_yields_convert_income(dynamic_chirps_log)
combined_income_varclust_log <- merge_yields_convert_income(dynamic_all_log)

# Generate sample graphics
ndvigraph_log <- chart_ce_cluster(get_cluster_diff_ce(ndvi_income_varclust_log), "NDVI", lower, upper)
chirpgraph_log <- chart_ce_cluster(get_cluster_diff_ce(chirp_income_varclust_log), "Precipitation", lower, upper)
evigraph_log <- chart_ce_cluster(get_cluster_diff_ce(evi_income_varclust_log), "EVI", lower, upper)
soilgraph_log <- chart_ce_cluster(get_cluster_diff_ce(soil_income_varclust_log), "Soil", lower, upper)
tempgraph_log <- chart_ce_cluster(get_cluster_diff_ce(temp_income_varclust_log), "Temperature", lower, upper)
combinedgraph_log <- chart_ce_cluster(get_cluster_diff_ce(combined_income_varclust_log), "Combined indicators", lower, upper)

make_cluster_plot_set(suffix = "log")

# 7.5 Join them all together ----

# Generate legend-free plots for each suffix
suffixes <- c("lowertail", "lowerlowertail", "log")
suffix_labels <- c("Lowertail", "Lower-lowertail", "Log")

all_plots <- lapply(suffixes, function(s) make_cluster_plot_set(s, 
                                                                highlight_coords = list(x = c(n_admin1, n_admin2), 
                                                                                        y = c(a1_zero, a2_zero)),
                                                                y_limits = c(-5, 15),
                                                                return_graphs_only = TRUE))

# Add row group labels using text grobs
row_titles <- lapply(suffix_labels, function(lbl) text_grob(lbl, rot = 90, face = "bold", size = 12))

# Add spacing row grobs to match columns
blank_col <- replicate(6, NULL, simplify = FALSE)  # one row of 6 blank cells

# Bind each row label and its 6 plots
combined_rows <- Map(function(label, plots) {
  c(list(label), plots)
}, row_titles, all_plots)

# Flatten into full list for ggarrange
plot_matrix <- do.call(c, combined_rows)

# Legend
legend_plot <- ggplot() +
  geom_point(aes(x = 1, y = 1, color = "Admin 1"), size = 3) +
  geom_point(aes(x = 2, y = 2, color = "Admin 2"), size = 3) +
  scale_color_manual("Predefined Boundaries", 
                     values = c("Admin 1" = "purple", "Admin 2" = "orange")) +
  theme_void()+
  theme(
    legend.direction = "horizontal",
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 12)
  )
legend_only <- cowplot::get_legend(legend_plot)

# Combine them all
grid_plot <- wrap_plots(plot_matrix, 
                        ncol = 7, nrow = 3,
                        widths = c(0.1, rep(1, 6)))

combined_output_dynamic_plot <- (
  grid_plot / legend_only + 
    plot_layout(heights = c(13, 1.5))
) + 
  plot_annotation(
    title = "Risk Reduction by Cluster Strategy and Input",
    theme = theme(
      plot.title = element_text(
        size = 16,
        face = "bold",
        hjust = 0.5
      ),
      plot.margin = ggplot2::margin(2,2,2,2, "mm")
    )
  )

print(combined_output_dynamic_plot)

# 7.6 Save ----
ggsave.latex(combined_output_dynamic_plot,
             filename = overleaf_path(paste0("/figs/combined_dynamic_lowertaillog_", date_today, ".png")),
             label = "fig:combined_dynamic_lowertaillog", 
             width = 13, 
             height = 7, 
             unit = "in")

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 8 Variance analysis ----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# This section generates Fig. 7 of the manuscript and Fig.A.9 of the appendix to the manuscript.

# 8.1 Define the datasets ----
datasets <- list(
  "evi_clusters" = evi_income_varclust_0,
  "ndvi_clusters" = ndvi_income_varclust_0,
  "temp_clusters" = temp_income_varclust_0,  # Replace with your actual dataset names
  "chirp_clusters" = chirp_income_varclust_0,
  "soil_clusters" = soil_income_varclust_0,
  "all_clusters" = combined_income_varclust_0
)

# 8.2 Func to process a single dataset ----
process_dataset <- function(var_data, dataset_name) {
  # Define existing clusters
  cluster_names <- paste0("cluster_", 5:100)
  existing_clusters <- cluster_names[cluster_names %in% names(var_data)]
  
  # Calculate group-level statistics for each cluster
  group_stats_by_cluster <- lapply(existing_clusters, function(cluster_name) {
    cluster_data <- var_data[[cluster_name]]
    cluster_num <- as.numeric(gsub("cluster_", "", cluster_name))
    
    group_stats <- cluster_data %>%
      group_by(group_name) %>%
      summarise(
        n_obs = n(),
        min_yield = min(y_it, na.rm = TRUE),
        max_yield = max(y_it, na.rm = TRUE),
        variance_yield = var(y_it, na.rm = TRUE) * ((n() - 1) / n()),
        mean_yield = mean(y_it, na.rm = TRUE),
        sd_yield = sd(y_it, na.rm=TRUE),
        cc_group =  (1.96 * sd_yield / (0.1*mean_yield))^2,
        .groups = 'drop'
      ) %>%
      mutate(cluster = cluster_name, cluster_num = cluster_num)
    
    return(group_stats)
  }) %>%
    bind_rows()
  
  # Calculate cluster-level summary statistics
  cluster_summary <- group_stats_by_cluster %>%
    group_by(cluster, cluster_num) %>%
    summarise(
      mean_sd = mean(sd_yield, na.rm=TRUE),
      cc_sum = sum(cc_group,na.rm=TRUE),
      mean_variance = mean(variance_yield, na.rm = TRUE),
      min_variance = min(variance_yield, na.rm = TRUE),
      max_variance = max(variance_yield, na.rm = TRUE),
      mean_obs_per_group = mean(n_obs, na.rm = TRUE),
      min_obs_per_group = min(n_obs, na.rm = TRUE),
      max_obs_per_group = max(n_obs, na.rm = TRUE),
      n_groups = n(),
      .groups = 'drop'
    ) %>%
    arrange(cluster_num) %>%
    mutate(dataset = dataset_name)
  
  return(list(
    group_stats = group_stats_by_cluster,
    cluster_summary = cluster_summary
  ))
}

# 8.3 Map func across all datasets ----
all_results <- map_dfr(names(datasets), function(name) {
  process_dataset(datasets[[name]], name)$cluster_summary
})

# 8.4 Format resulting dataframe ----
# Rename the clustering variables for higher clarity in final graph: 
all_results <- all_results %>% 
  mutate(dataset = recode(dataset, 
                          evi_clusters = "EVI", 
                          ndvi_clusters = "NDVI", 
                          chirp_clusters = "Precipitation", 
                          soil_clusters = "Soil", 
                          temp_clusters = "Temperature", 
                          all_clusters = "All indicators combined")) 
 
insurancezones <- list("EVI", "NDVI", "Precipitation", "Soil", "Temperature", "All indicators combined")

# 8.5 Prepare plotting ----

# Define colors for different datasets 
colors <- c(
  "EVI" = "darkgreen",               # Green for vegetation (EVI)
  "NDVI" = "#66bd63",             # Lighter green for NDVI
  "Precipitation" = "#4575b4",    # Blue for rain
  "Soil" = "#a6611a",
  "Temperature" = "#d73027",
  "All indicators combined" = "#54278f"  # Purple as a composite/neutral
)
names(colors) <- insurancezones

# Define variables for plotting
all_results$dataset <- factor(all_results$dataset,
                              levels = c("EVI", "NDVI", "Precipitation", "Temperature", "Soil", "All indicators combined"))

# Common linetype and theme
linetypes <- c("Minimum" = "dotted", "Maximum" = "dashed", "Mean" = "solid")

# Define formatting guide for plot
clean_theme <- theme_minimal(base_size = 12) +
  theme(
    panel.grid.major = element_line(color = "grey85"),
    panel.grid.minor = element_line(color = "grey93"),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 11.5),
    axis.title.x = element_text(size = 11),
    axis.title.y = element_text(size = 11)
  )

# Define cpmmon axis and legend
scale_guides_common <- list(
  scale_color_manual( values = colors, name = "Insurance Zones"),  # Groups: EVI, NDVI, etc.
  scale_linetype_manual(values = linetypes,name = "Statistic") # Groups: Mean, Max, Min
)

# Design plot
make_panel_plot <- function(data, min_var, mean_var, max_var, y_label, title_text, y_max = NULL) {
  ymax_val <- if (is.null(y_max)) {
    # default: compute max across all lines
    max(dplyr::pull(data, {{ min_var }}),
        dplyr::pull(data, {{ mean_var }}),
        dplyr::pull(data, {{ max_var }}),
        na.rm = TRUE) * 1.05  # add 5% headroom
  } else {
    y_max
  }
  
  ggplot(data, aes(x = cluster_num)) +
    geom_line(aes(y = {{ min_var }}, color = dataset, linetype = "Minimum"),
              linewidth = 0.75, alpha = 0.9) +
    geom_line(aes(y = {{ max_var }}, color = dataset, linetype = "Maximum"),
              linewidth = 0.8, alpha = 0.9) +
    geom_line(aes(y = {{ mean_var }}, color = dataset, linetype = "Mean"),
              linewidth = 1, alpha = 0.9) +
    labs(title = title_text, y = y_label, x = NULL) +
    scale_color_manual(values = colors, name = "Insurance Zones") +
    scale_linetype_manual(values = linetypes, name = "Statistic") +
    clean_theme +
    scale_y_continuous(
      expand = c(0, 0),
      limits = c(0, ymax_val)
    ) +
    scale_x_continuous(
      expand = c(0, 0),
      limits = c(0, 100),
      breaks = seq(0, 100, by = 10),
      minor_breaks = seq(0, 100, by = 5)
    ) +
    theme(legend.position = "none")
}

# 8.6 Generate Fig.7 ----
# Filter to just the top performers (soil and all)
temp_all <- all_results %>% filter(
  dataset %in% c("All indicators combined", "Temperature"))

# First two use the helper
variance_plot <- make_panel_plot(temp_all, 
                                 min_var = min_variance,
                                 mean_var = mean_variance,
                                 max_var = max_variance,
                                 y_label = "",
                                 title_text = "Variance by Number of Zones")

observations_plot <- make_panel_plot(temp_all, #all_results,
                                     min_var = min_obs_per_group,
                                     mean_var = mean_obs_per_group,
                                     max_var = max_obs_per_group,
                                     y_label = "",
                                     title_text = "Observations per Zone")

crop_cuts_plot <- ggplot(temp_all, aes(x = cluster_num, y = cc_sum, color = dataset)) +
  geom_line(linewidth = 1) +
  labs(title = "Total Crop Cuts Needed to Estimate Mean Yields with a (+/-) 10% Margin of Error",
       y = "",
       x = "Number of Insurance Zones") +
  clean_theme +
  theme(legend.position = "none") +
  scale_y_continuous(expand = c(0,0), 
                     limits = c(0,6000)) +
  scale_x_continuous(
    expand = c(0, 0),
    limits = c(0, 100),
    breaks = seq(0, 100, by = 10),
    minor_breaks = seq(0, 100, by = 5)) +
  scale_color_manual(values = colors, name = "Insurance Zones")

# Create combined plot with shared legend 
legend_plot_var <- 
  observations_plot +
  guides(
    color = guide_legend(title = "Zone\nInput", direction = "horizontal", nrow = 3, order = 1),
    linetype = guide_legend(title = "Statistic", direction = "horizontal", nrow = 3, order = 2)
  ) +
  theme_void() +
  theme(
    legend.position = "right",  # TEMPORARY so `get_legend()` captures it
    legend.box = "horizontal",
    legend.box.just = "center",
    legend.title = element_text(size = 11, face = "bold"),
    legend.text = element_text(size = 10),
    legend.spacing.x = unit(1, "cm"), # This line controls horizontal spacing
    legend.spacing.y = unit(0.1, "cm"),
    legend.key.height = unit(0.6, "lines") # tighter rows
    )

legend_var <- cowplot::get_legend(legend_plot_var)

# Combine plots using patchwork + cowplot for legend
final_plot <- patchwork::wrap_plots(variance_plot, 
                                    #observations_plot,
                                    crop_cuts_plot, ncol = 1) +
  plot_layout(guides = "collect") +
  plot_annotation(theme = theme(legend.position = "none"))

# Display full graphic with legend
display_final <- cowplot::plot_grid(final_plot, legend_var, ncol = 1, rel_heights = c(1, 0.2))
display_final

# Save/export the final plot
ggsave.latex(display_final,
             filename = overleaf_path(paste0("/figs/variance_analysis_", date_today, ".png")),
             label = "fig:variance_analysis", 
             width = 7, height=5.5, unit = "in"
)

# 8.7 Generate Fig.A.9 ----
# Same approach again but with a different dataset (all_results instead of temp_all)
variance_plot <- make_panel_plot(all_results, 
                                 min_var = min_variance,
                                 mean_var = mean_variance,
                                 max_var = max_variance,
                                 y_label = "",
                                 title_text = "Variance by Number of Zones")

observations_plot <- make_panel_plot(all_results,
                                     min_var = min_obs_per_group,
                                     mean_var = mean_obs_per_group,
                                     max_var = max_obs_per_group,
                                     y_label = "",
                                     title_text = "Observations per Zone")


crop_cuts_plot <- ggplot(all_results, aes(x = cluster_num, y = cc_sum, color = dataset)) +
  geom_line(linewidth = 1) +
  labs(title = "Total Crop Cuts Needed to Estimate Mean Yields with a (+/-) 10% Margin of Error",
       y = "",
       x = "Number of Insurance Zones") +
  clean_theme +
  theme(legend.position = "none") +
  scale_y_continuous(expand = c(0,0), 
                     limits = c(0,6000)) +
  scale_x_continuous(
    expand = c(0, 0),
    limits = c(0, 100),
    breaks = seq(0, 100, by = 10),
    minor_breaks = seq(0, 100, by = 5)) +
  scale_color_manual(values = colors, name = "Insurance Zones")

# Create combined plot with Shared Legend 
legend_plot_var <- 
  observations_plot +
  guides(
    color = guide_legend(title = "Zone\nInput", direction = "horizontal", nrow = 3, order = 1),
    linetype = guide_legend(title = "Statistic", direction = "horizontal", nrow = 3, order = 2)
  ) +
  theme_void() +
  theme(
    legend.position = "right",  # TEMPORARY so `get_legend()` captures it
    legend.box = "horizontal",
    legend.box.just = "center",
    legend.title = element_text(size = 11, face = "bold"),
    legend.text = element_text(size = 10),
    legend.spacing.x = unit(1, "cm"), # This line controls horizontal spacing
    legend.spacing.y = unit(0.1, "cm"),
    legend.key.height = unit(0.6, "lines") # tighter rows
  )

legend_var <- cowplot::get_legend(legend_plot_var)

# Combine plots using patchwork + cowplot for legend
final_plot <- patchwork::wrap_plots(variance_plot, 
                                    crop_cuts_plot, ncol = 1) +
  plot_layout(guides = "collect") +
  plot_annotation(theme = theme(legend.position = "none"))

# Display full graphic with legend
display_final <- cowplot::plot_grid(final_plot, legend_var, ncol = 1, rel_heights = c(1, 0.2))
display_final

# Save/export the final plot
ggsave.latex(display_final,
             filename = overleaf_path(paste0("/figs/app_variance_analysis_", date_today, ".png")),
             label = "fig:variance_analysis",
             width = 7, height=5.5, unit = "in"
)




#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 9 Heatmaps on insurance value (Fig.5)  ----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#This section generates the heatmaps presented in Fig. 5 of the manuscript. 

# # 9.0 Define funcCalculate data for heatmaps ----
# # This step is computationally very expensive. 
# # Hence, we saved the results. It is therefore commented out at this point. 
# 
# future::plan(multisession, workers = parallel::detectCores() - 2)  # leave 2 cores free
# 
# gc() # garbage collection to free up memory
# 
# cropify <- function(num_cc, design_risk=TRUE){
#   clusters <- merge_yields_convert_income(dynamic_all, num_cc=num_cc, design_risk=design_risk)
#   results <- get_cluster_diff_ce(clusters)
#   return(results)
# }
# 
# multicropify <- function(num_cc, reps=24, design_risk=TRUE){
#   results <- future_map_dfr(1:reps, ~ cropify(num_cc=num_cc, design_risk=design_risk))
#   
#   # Now average diff_ce by cluster
#   average_results <- results %>%
#     group_by(cluster) %>%
#     summarise(diff_ce = mean(diff_ce), .groups = "drop")
#   
#   # View the result
#   return(average_results)
# }
# 
# # Set general params:
# cut_min = 1
# cut_delta = 1
# 
# # a) Zero cost scenario
# zone_cost = 0 #TODO: does this need to be changed?
# cost_cc= 0 # assuming USD 0/crop cut for this heatmap
# cuts_per_zone <- seq(cut_min, 100, by=cut_delta)
# 
# results <- lapply(cuts_per_zone, multicropify)
# 
# combined <- bind_rows(results, .id = "cuts") %>% mutate(diff_ce_capped = ifelse(diff_ce < -5, -5, diff_ce)) # cap diff_ce at -5 for better visualization
# combined$cuts <- as.numeric(combined$cuts) * cut_delta + cut_min
# 
# # b) Medium cost scenario
# cost_cc= 10 # assuming USD 10/crop cut for this heatmap
# results10 <- lapply(cuts_per_zone, multicropify)
# 
# combined10 <- bind_rows(results10, .id = "cuts") %>% mutate(diff_ce_capped = ifelse(diff_ce < -5, -5, diff_ce)) # cap diff_ce at -5 for better visualization
# combined10$cuts <- as.numeric(combined10$cuts) * cut_delta + cut_min
# 
# # c) Medium cost scenario IGNORING design risk
# cost_cc= 10
# results10nd <- lapply(cuts_per_zone, multicropify, design_risk=FALSE)
# 
# combined10nd <- bind_rows(results10nd, .id = "cuts") %>% mutate(diff_ce_capped = ifelse(diff_ce < -5, -5, diff_ce)) # cap diff_ce at -5 for better visualization
# combined10nd$cuts <- as.numeric(combined10nd$cuts) * cut_delta + cut_min
# 
# # d) High cost scenario (not presented in manuscript)
# cost_cc= 36 # assuming USD 36/crop cut for this heatmap
# results36 <- lapply(cuts_per_zone, multicropify)
# 
# combined36 <- bind_rows(results36, .id = "cuts") %>% mutate(diff_ce_capped = ifelse(diff_ce < -5, -5, diff_ce)) # cap diff_ce at -5 for better visualization
# combined36$cuts <- as.numeric(combined36$cuts) * cut_delta + cut_min
# 
# # Save data
# write.csv(combined, file = "datasets/heatmap_cost0.csv")
# write.csv(combined10, file = "datasets/heatmap_cost10.csv")
# write.csv(combined10nd, file = "datasets/heatmap_cost10_nodesignrisk.csv")

# 9.1 Load the data for the heatmaps ----
combined <- read.csv("datasets/heatmap_cost0.csv")
combined10 <- read.csv("datasets/heatmap_cost10.csv")
combined10nd <- read.csv("datasets/heatmap_cost10_nodesignrisk.csv")

# 9.2 Generate separate heatmaps ----

# Define function
make_heatmap_plot <- function(data, subtitle = ""){
  plot <- data %>%
    ggplot(aes(y = cuts, x = cluster, fill = diff_ce)) +
    geom_tile(color = "white") +
    scale_fill_gradient2(
      low = "#ef8a62", mid = "#f7f7f7", high = "dark blue",
      midpoint = 0,
      limits = c(-10, 10),   # squished range
      oob = squish,                    # keeps out-of-bound values within color scale
      name = "Diff.\nin CE",
      labels = label_number()
    ) +
    scale_x_continuous(
      breaks = scales::pretty_breaks(),
      limits = c(0, 101),
      minor_breaks = seq(0, 101, by = 10),
      expand = c(0, 0)
      ) +
    scale_y_continuous(
      breaks = scales::pretty_breaks(),
      limits = c(0, 101),
      minor_breaks = seq(0, 101, by = 10),
      expand = c(0, 0)
    ) +
    coord_fixed(xlim = c(0, 105), ylim = c(0, 101)) +
    labs(y = "Crop Cuts\nper Zone",
         x = "Number of Zones",
         #title = "Net value by number of crop cuts and zones",
         subtitle = subtitle) +
    guides(
      x = guide_axis(minor.ticks = TRUE),
      y = guide_axis(minor.ticks = TRUE)
    )  +
    clean_chart_clutter_explore

  return(plot)
}

# 9.2.1 Heatmap for zero cost scenario (panel a)
heatmap_plot <- 
  make_heatmap_plot(combined,
                    subtitle="a) USD 0/crop cut") +
                    theme(legend.position = "none", 
        axis.title.y = element_text(angle = 0, vjust = 0.5))
heatmap_plot

# 9.2.2 Heatmap for medium cost scenario (panel b)
heatmap_plot10 <- 
  make_heatmap_plot(combined10, 
                    subtitle= "b) USD 10/crop cut") +
  theme(legend.position = "none", axis.title.y = element_blank())
heatmap_plot10 


# 9.2.3 Heatmap for medium cost scenario IGNORING design risk (panel c)
heatmap_plot10nd_wlegend <- 
  make_heatmap_plot(combined10nd, 
                    # subtitle = "" ) 
                    subtitle= "c) USD 10/crop cut, no design risk") +
  theme(
    legend.position = "bottom",
    legend.justification = "center",         # Align whole legend box to the left
    legend.title.align = 0,                # Left-align the legend title text
    legend.box.just = "left",              # Align colorbar + text box to the left
    legend.margin = margin(t = 0, b = 0),  # Tighter vertical spacing if needed
    legend.text = element_text(margin = margin(l = 10, r = 4)),  # Nudges text closer to scale bar
    legend.title = element_text(margin = margin(r = 15))         # Nudges title leftward
  )
  
legend <- get_legend(heatmap_plot10nd_wlegend)
legend 

heatmap_plot10nd <- 
  heatmap_plot10nd_wlegend + 
  theme(legend.position = "none", axis.title.y = element_blank())

heatmap_plot10nd

# 9.3 Combine the three heatmaps ----
main_row <- heatmap_plot + heatmap_plot10 + heatmap_plot10nd +
  plot_layout(
    ncol = 3,
    widths = c(1, 1, 1),
    guides = "collect"
  ) &
  theme(
    panel.spacing = unit(0.25, "lines"),   # try smaller values like 0.25
    plot.margin = margin(2, 2, 2, 2)       # reduce outer whitespace
  )

# Combine with legend
combined_plot <- main_row / legend +
  plot_layout(heights = c(1, 0.12))  # adjust legend row height

combined_plot

# 9.4 Save Fig.5 ----
ggsave.latex(combined_plot,
             filename = overleaf_path(paste0("/figs/combined_heatmap", date_today,".png")),
             label = "fig:combined_heatmap" #,
             # width = 11,     
             # height = 3.5 
)

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 10 Heatmap slices (Fig.6)----
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# This section extracts and visualizes slices of the heatmap 
# from panel b in Fig 5. These slices are presented in Fig.6 of the manuscript.

# Merge datasets to be able to depict values with and without design risk:
combined10nd$designrisk <- "No"
combined10$designrisk <- "Yes"
combined10all <- rbind(combined10, combined10nd)

# Function to plot the values
plot_diff_ce <- function(data, x_var, x_label, plot_title=NULL) {
  data %>% 
    ggplot(aes(x = {{ x_var }},
               y = diff_ce,
               color = factor(designrisk),
               group = designrisk)) +      # ensures one line per level
    geom_line() +
    scale_color_manual(values = c("#004488", "#DDAA33")
    ) +
    labs(
      title = plot_title,
      x = x_label,                         # x‑axis label
      y = "Diff. in CE\nin USD",                           # y‑axis label
      color = "Design Risk"                # legend title
    ) +
    
    coord_cartesian(ylim = c(-10, 10)) +      # show only 0–10 on y‑axis
    theme_minimal()
}

# Generating slices for best performing cases: 80 crop cuts and 6 zones (according to panel b in Fig 5)
# Generate slice graphic for 80 cropcuts 
cuts80_plot <- combined10all %>%
  filter(cuts==80) %>%
  plot_diff_ce(cluster, x_label="Zones", plot_title="b) across varying number of zones\nat 80 crop cuts per zone")
cuts80_plot <- cuts80_plot + theme(axis.title.y = element_blank())
cuts80_plot

# Generate slice graphic for 6 zones
zones6_plot <- combined10all %>%
  filter(cluster==6) %>%
  plot_diff_ce(cuts, x_label="Cuts per Zone", plot_title="a) across varying number of crop cuts\nat 6 insurance zones ")
zones6_plot <- zones6_plot + theme(legend.position = "none", axis.title.y = element_text(angle = 0, vjust = 0.5))
zones6_plot

# Combine them
combo_plot <- zones6_plot + cuts80_plot +plot_layout(ncol = 2) + plot_annotation(title = "Development of Net Risk Reduction Value ", theme = theme(plot.title = element_text(hjust = 0.5))) 
combo_plot

# Save 
ggsave.latex(combo_plot,
             filename = overleaf_path(paste0("/figs/heatmap_slices_", date_today,".png")),
             label = "fig:heatmap_slices",
             height = 4, 
             width = 10)

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 11  Examine relationship between clusters and yields (appendix A.9) ----
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# 11.1 Yield boxplots by zone (Fig.A.10) ----
# (k=14, separately for each agri-environmental input)

# extract a dataset from oaf that has yield, name_concat2, then all the `c.` columns
oaf_slim <- 
  oaf %>%
  select(yield, name_1, starts_with("c.")) %>%
  sf::st_drop_geometry()

# Reshape all `c.` columns into long format for plotting
oaf_slim_long <- 
  oaf_slim %>% 
  mutate(c.admin1 = as.numeric(factor(name_1))) %>% 
  sf::st_drop_geometry() %>% 
  pivot_longer(cols = starts_with("c."), names_to = "index", values_to = "cluster") %>%
  mutate(
    cluster = as.integer(cluster),
    index = as.factor(index)
  )

oaf_slim_wide <- 
  oaf_slim_long %>%
  count(index, cluster) %>%
  complete(index, cluster = 1:14, fill = list(n = 0)) %>%
  pivot_wider(names_from = cluster, values_from = n, names_sort = TRUE) %>%
  arrange(index)

# Function to make the plot
make_cluster_plot <- function(index_name, df, show_x_label = TRUE) {
  df_filtered <- df %>% filter(index == index_name)
  
  # Get display name, fallback to original if not found
  # Use case_when for display name
  display_name <- case_when(
    index_name == "c.chirp" ~ "Precipitation",
    index_name == "c.admin1" ~ "Admin 1", 
    index_name == "c.soil" ~ "Soil",
    index_name == "c.evi" ~ "EVI",
    index_name == "c.ndvi" ~ "NDVI",
    index_name == "c.temp" ~ "Temperature",
    index_name == "c.yield" ~ "Yields",
    index_name == "c.all" ~ "Combined indicators",
    TRUE ~ index_name  # fallback to original name
  )
  
  # Re-order clusters by mean yield and assign new labels
  cluster_summary <- df_filtered %>%
    group_by(cluster) %>%
    summarise(
      mean_yield = mean(yield, na.rm = TRUE),
      sd_yield = sd(yield, na.rm = TRUE),
      count = sum(!is.na(yield)),
      .groups = "drop"
    ) %>%
    arrange(mean_yield) %>%
    mutate(
      new_cluster = row_number(),
      label = paste0(new_cluster, "\n(", count, ")\n[", round(sd_yield, 1), "]")
    )
  
  # Join back to full data
  df_plot <- df_filtered %>%
    left_join(cluster_summary, by = "cluster") %>%
    mutate(label = factor(label, levels = cluster_summary$label))
  
  # Compute cluster means for overlay
  cluster_means <- df_plot %>%
    group_by(label) %>%
    summarise(mean_yield = mean(yield, na.rm = TRUE), .groups = "drop")
  
  # Plot
  ggplot(df_plot, aes(x = label, y = yield)) +
    geom_boxplot(outlier.shape = NA, color = "black", fill = NA) +
    geom_line(data = cluster_means, aes(x = label, y = mean_yield, group = 1),
              inherit.aes = FALSE, color = "blue", linewidth = 1) +
    labs(
      title = display_name,
      x = if (show_x_label) "\n Reindexed Zone, ordered by zonal mean yield \n(Obs in Zone) \n[SD]" else NULL,
      y = "Yield (t/ha)"
    ) +
    theme_minimal(base_size = 13)
  }

index_list <- unique(oaf_slim_long$index); index_list

plots <- map2(index_list, seq_along(index_list), ~ 
                make_cluster_plot(.x, oaf_slim_long, show_x_label = .y > 6)
)

nrows <- ceiling(length(index_list) / 2)
rel_heights <- c(rep(1, nrows - 1), 1.4)

output_boxplots <- plot_grid(plotlist = plots, ncol = 2, #align = "hv",   
          axis = "tblr", 
          greedy = TRUE,
          rel_heights = rel_heights)

output_boxplots

# Save Fig.A.10
ggsave.latex(output_boxplots,
             filename = overleaf_path(paste0("/figs/cluster_boxplots_", date_today, ".png")),
             width = 13, 
             height = 10, 
             unit = "in")

# 11.2 Dot-whisker diagnostic (Fig.A.11) ---- 
colors <- c(
  "EVI" = "darkgreen",
  "NDVI" = "#66bd63",
  "Precipitation" = "#4575b4",
  "Soil" = "#8c510a",
  "Temperature" = "#d73027",
  "All indicators combined" = "#756bb1",
  "Admin 1 zones" = "#ff7f00",     # Orange
  "Yield-based zones" = "#f781bf"         # Pink/magenta
)

label_map <- c(
  "c.evi" = "EVI",
  "c.ndvi" = "NDVI",
  "c.chirp" = "Precipitation",
  "c.soil" = "Soil",
  "c.temp" = "Temperature",
  "c.all" = "All indicators combined",
  "c.admin1" = "Admin 1 zones",
  "c.yield" = "Yield-based zones"
)

cluster_stats <- 
  oaf_slim_long %>%
  group_by(index, cluster) %>%
  summarise(
    mean_yield = mean(yield, na.rm = TRUE),
    sd_yield = sd(yield, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  group_by(index) %>%
  arrange(mean_yield) %>%
  mutate(reindexed_cluster = row_number()) %>%
  ungroup()

cluster_stats <- cluster_stats %>%
  group_by(reindexed_cluster) %>%
  arrange(sd_yield, .by_group = TRUE) %>%
  mutate(index_order = row_number()) %>%
  ungroup() %>%
  mutate(index_label = recode(index, !!!label_map))

dot_whisker <- 
  cluster_stats %>% 
  ggplot(aes(x = reindexed_cluster, y = mean_yield, color = index_label)) +
  geom_point(aes(group = index_order), position = position_dodge(width = 0.6), size = 2) +
  geom_errorbar(aes(ymin = mean_yield - sd_yield,
                    ymax = mean_yield + sd_yield,
                    group = index_order),
                position = position_dodge(width = 0.6), width = 0.3) +
  scale_color_manual(values = colors) +
  scale_x_continuous(breaks = seq(min(cluster_stats$reindexed_cluster),
                                  max(cluster_stats$reindexed_cluster), by = 1)) +
  scale_y_continuous(breaks = seq(floor(min(cluster_stats$mean_yield)),
                                  ceiling(max(cluster_stats$mean_yield) + 3), by = 1)) +
  labs(
    x = "\n Zones re-indexed by mean yield, inputs sorted by SD",
    y = "Yield",
    title = "Zonal Mean Yields ± SD by Input Data Source",
    color = "Zone Definition Method"
  ) +
  # theme_minimal(base_size = 13) +
  clean_chart_clutter_explore + 
  theme(legend.position = "bottom", 
        panel.grid.major = element_line(color = "grey80", linewidth = 0.3),
       panel.grid.minor = element_line(color = "grey90", linewidth = 0.2))

dot_whisker 

# Save Fig.A.11
ggsave.latex(dot_whisker,
             filename = overleaf_path(paste0("/figs/cluster_dot_whisker_", date_today, ".png")),
             label = "fig:cluster_dot_whisker")