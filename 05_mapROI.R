#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Author: Ella Kirchner 
# Project name: RS Get in the Zone
# Date Last Updated: # Fri Sep 12 11:34:26 2025 
# Purpose: Create map of study region
# Input files: 
      # - GADM data for Kenya and surrounding countries (loaded via the script)
      # - Shapefile of global oceans 
      # - OAF data 
# Output files: 
      # - Fig. 1: Map of study region
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


# %%%%%%%%%%%%%%%%%%%%%%%%%%%
# Set Up 
# %%%%%%%%%%%%%%%%%%%%%%%%%%%

# 1. Clean out workspace
rm(list = ls(all = TRUE))

# 2. Set Up File path & directory: 
file_path <-  "put_your_file_path_here" # File path to project folder
setwd(file_path)

saving <- paste0(file_path,"/figs")

# 3. Set Seed for replication 
set.seed(123456789)

# %%%%%%%%%%%%%%%%%%%%%%%%%%%
# Load libraries 
# %%%%%%%%%%%%%%%%%%%%%%%%%%%
library(dplyr)
library(tidyverse)
library(ggplot2)
library(sp)
library(sf) 
library(terra)
library(gridExtra)
library(ggpubr)
library(ggspatial)
library(ggpattern)
library(geodata)
library(grid)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 0. Preparations: ----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

remove_chart_clutter <- theme(
  panel.background = element_rect(fill="white", color="black"),
  panel.grid= element_blank(),
  panel.grid.major= element_blank(),
  panel.grid.minor= element_blank(),
  legend.title = element_text(size=16),
  legend.text = element_text(size=14),
  legend.position = "right", 
  axis.title.x = element_blank(),
  axis.title.y = element_blank(),
  axis.text = element_text(size=14),
  legend.direction= "vertical", 
  legend.box = "vertical"
)

# Create list of 14 admin1 districts (counties) included in analysis:
districts_to_keep <- c("Bungoma", "Busia", "Homa Bay", "Kakamega", 
                       "Kericho", "Kisumu", "Kisii", "Migori", "Nandi", "Nyamira", 
                       "Siaya", "Trans Nzoia", "Uasin Gishu", "Vihiga")

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#1. Load data: ----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
kenya_0 <- st_as_sf(gadm("GADM", level=0, country="KEN"))
kenya_1 <- st_as_sf(gadm("GADM", level=1, country="KEN"))
kenya_2 <- st_as_sf(gadm("GADM", level=2, country="KEN"))
kenya_3 <- st_as_sf(gadm("GADM", level=3, country="KEN"))

# Load data for surrounding countries and ocean: 
uganda_0 <- st_as_sf(gadm("GADM", level=0, country="UGA"))
tanzania_0 <- st_as_sf(gadm("GADM", level=0, country="TZA"))
ethiopia_0 <- st_as_sf(gadm("GADM", level=0, country="ETH"))
somalia_0 <- st_as_sf(gadm("GADM", level=0, country="SOM"))
sudan_0 <- st_as_sf(gadm("GADM", level=0, country="SSD"))
ocean <- read_sf(paste0(file_path, "/datasets/water-polygons-split-4326/water_polygons.shp"))

# # Oaf Data:
# load("datasets/oaf.Rdata")

# Anonymized oaf data instead (without exact geolocations):
load("datasets/oaf_anonym.Rdata")
oaf <- oaf_anonym

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 2. Country map: ----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Prepare data:
kenya_1 <- kenya_1 %>%
  mutate(
    included = ifelse(NAME_1 %in% districts_to_keep, 1,0)
  )
roi <- st_as_sf(subset(kenya_1, NAME_1 %in% districts_to_keep))

# Define extent of country map: 
xmin=32
xmax=43.5
ymin=-5.05
ymax=6

# Crop the surrounding countries and the ocean to the area that we want to map: 
uganda_cropped <- st_crop(uganda_0, xmin = xmin, xmax = xmax, ymin = ymin, ymax =ymax)
tanzania_cropped <- st_crop(tanzania_0, xmin = xmin, xmax = xmax, ymin = ymin, ymax =ymax)
ethiopia_cropped <- st_crop(ethiopia_0, xmin = xmin, xmax = xmax, ymin = ymin, ymax =ymax)
somalia_cropped <- st_crop(somalia_0, xmin = xmin, xmax = xmax, ymin = ymin, ymax =ymax)
sudan_cropped <- st_crop(sudan_0, xmin = xmin, xmax = xmax, ymin = ymin, ymax =ymax)
ocean_cropped <- ocean %>%
  filter(x>32) %>%
  filter(x<43.5) %>%
  filter(y<6) %>%
  filter(y>-7)
ocean_cropped2 <- st_union(ocean_cropped)
ocean_cropped3 <- st_crop(ocean_cropped2, xmin = xmin, xmax = xmax, ymin = ymin, ymax =ymax)

# Location of ROI: 
location_roi <- ggplot(kenya_1) + 
  geom_sf(fill = "white")+
  theme(plot.background = element_blank())+
  geom_sf(data=kenya_1, aes(fill=ifelse(NAME_1 %in% districts_to_keep, "Study area", "Not included in study")), alpha=0.9)+
  geom_sf(data=ocean_cropped3, fill="lightgrey", color="grey")+
  geom_sf(data=uganda_cropped, fill="white", color="grey")+
  geom_sf(data=tanzania_cropped, fill="white", color="grey")+
  geom_sf(data=ethiopia_cropped, fill="white", color="grey")+
  geom_sf(data=somalia_cropped, fill="white", color="grey")+
  geom_sf(data=sudan_cropped, fill="white", color="grey") +
  geom_sf(data=kenya_0, fill="NA", color="black", size=2, alpha=.3)+
  geom_text(data = data.frame(lon = c(33.5, 33, 38, 42.2, 35), lat = c(2, 5, 5.5, 2, -4), label = c("Uganda", "South\nSudan", "Ethiopia", "Somalia", "Tanzania")), aes(x = lon, y = lat, label = label), size = 5, color = "grey") +
  scale_fill_manual(values = c("Study area" = "grey", "Not included in study" = "white")) +
  annotate("rect", xmin=33.8, xmax=36, ymin=-1.6, ymax=1.6, color="red", fill=NA)+
  # geom_rect(aes(xmin=33.8, xmax=36, ymin=-1.6, ymax=1.6), color="red", fill="NA")+
  remove_chart_clutter +
  theme(
    legend.position="none", 
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks = element_blank(),
    panel.background = element_rect(fill="white", color="NA")  
  )+  
  annotate("label", x=38, y= 0.5, label="Kenya", size=14/.pt, color="black")

location_roi

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 3. Map of the ROI (region of interest): ----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Filter oaf data for wards of interest:
all_years <- unique(oaf$year)

oaf_wallyears <- oaf %>%
  group_by(name_concat) %>%
  filter(length(unique(year)) == length(all_years)) %>% # drops 7k observations 
  ungroup()%>%
  select(name_concat, name_1, name_3, starts_with("c")) %>%
  # mutate(geometry = NULL ) %>%
  as_tibble() %>%
  distinct()

# Filter Kenya data to ROI:
kenya_1f <- filter(kenya_1, NAME_1 %in% districts_to_keep)
kenya_2f <- filter(kenya_2, NAME_1 %in% districts_to_keep)
kenya_3f <- filter(kenya_3, NAME_1 %in% districts_to_keep)
 
#Add dummies on whether name_concat is in study or not:
wards_available_allyears <- as.list(oaf_wallyears$name_concat)
kenya_3f <- kenya_3f %>%
  mutate(
    name_concat= paste(NAME_1, NAME_2, NAME_3, sep="_"),
    ward_available_allyears =ifelse(name_concat %in% wards_available_allyears, 1,0),
    ward_available_allyears_letter=ifelse(name_concat %in% wards_available_allyears, "A", "B"),
  )

# Map of study region annotated with real "labels": 
studyregion <- ggplot()+ 
  geom_sf(data=kenya_3f, aes(fill=ward_available_allyears_letter), linewidth=0.5)+
  scale_fill_manual(
    values=c(A='lightgrey',B='white'),
    labels = c("Ward with at least one crop cut per year (2016-20)", "Wards with gaps in crop cut coverage across years"),
  )+
  geom_sf(data=kenya_3f, aes(color="Ward (administrative level 3)"), fill=NA, linewidth=0.6)+
  geom_sf(data=kenya_2f, aes(color="Subcounty (administrative level 2)"), fill=NA, linewidth=0.7)+
  geom_sf(data=kenya_1f, aes(color="County (administrative level 1)"), fill=NA, linewidth=1.4) +
  scale_color_manual(
    values = c("Ward (administrative level 3)" = "darkgrey", "Subcounty (administrative level 2)" = "black", "County (administrative level 1)" = "black"),
    name = "Administrative Boundaries"
  ) +
  guides(
    fill = guide_legend(override.aes = list(color = NA)),
    color = guide_legend(
      override.aes = list(
        fill = "white", 
        linewidth = c("Admin1" = 1.4, "Admin2" = 0.7, "Admin3" = 0.6)
      )))+
  ggspatial::annotation_north_arrow(location='br', height=unit(1,"cm"), width = unit(1,"cm"))+
  ggspatial::annotation_scale(style='ticks', unit_category='metric', location='bl', text_cex=1.1)+
  labs(
    fill="Availability of crop cut data in wards"
  )+
  theme(
    text=element_text(size=12)
  )+
  remove_chart_clutter+
  annotate("label", x=34.9, y=1.1, label= "Trans Nzoia", size=14/.pt)+
  annotate("label", x=35.3, y=0.5, label= "Uasin\nGishu", size=14/.pt)+
  annotate("label", x=34.8, y=0.7, label="Bungoma", size=14/.pt)+
  annotate("label", x=34.24, y=0.38, label="Busia", size=14/.pt)+
  annotate("label", x=34.67, y=0.33, label="Kakagema", size=14/.pt)+
  annotate("label", x=35.1, y=0.29, label="Nandi", size=14/.pt)+
  annotate("label", x=34.2, y=-0.1, label="Siaya", size=14/.pt)+
  annotate("label", x=34.3, y=-0.6, label="Homa Bay", size=14/.pt)+
  annotate("label", x=34.22, y=-0.96, label="Migori", size=14/.pt)+
  annotate("label", x=34.75, y=-0.80, label="Kisii", size=14/.pt)+
  annotate("label", x=34.97, y=-0.65, label="Nya-\nmira", size=14/.pt)+
  annotate("label", x=35.49, y=-0.18, label="Kericho", size=14/.pt)+
  annotate("label", x=34.85, y=-0.18, label="Kisumu", size=14/.pt)+
  annotate("label", x=34.71, y= 0.07, label="Vihiga", size=14/.pt)+
  coord_sf(ylim=c(-1.5,1.3),xlim=c(33.9,35.7), clip="off")
studyregion

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 4. Combine both maps: ----

small <- ggplotGrob(location_roi)

map <- studyregion + 
  annotation_custom(
    grob = small, 
    xmin=35.65, 
    xmax=36.95, 
    ymin=-1.4, 
    ymax=-0.4)+
  annotation_custom(
    grob = text_grob("Data sources: GADM 4.1\nProduced by Ella Kirchner\nLast updated: June 2025\n ", hjust=0, just="left", size=14),
    xmin=35.62, 
    xmax=36.1, 
    ymin=-1.6, 
    ymax=-1.5)
map 

# Save Fig.1
ggsave(paste0(saving,"/studyregiontest.png"), plot= map, width=14, height= 12)


