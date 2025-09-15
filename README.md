# GetInTheZone_JDE
Replication package for the scientific journal article: *Kirchner et al. (2025). Get in the Zone: the Risk-Adjusted Welfare Effects of Data-Driven vs. Administrative Borders for Index Insurance Zones. Forthcoming in Journal of Development Economics.*

*Authors/Team:* Ella Kirchner, Elinor Benami, Andrew Hobbs, Michael Carter, Zhenong Jin 

*Abstract/Project Description:* Agricultural index insurance aims to protect producers from widespread shocks within defined areas, often based on administrative boundaries. However, these boundaries may poorly reflect yield variation, imposing costs on both policyholders and the public. Advances in geospatial data accessibility and clustering techniques offer potential to define more homogeneous insurance zones. This study evaluates the impact of redrawing zone boundaries based on observed agri-environmental conditions in western Kenya. Using over 13,000 crop cut observations with satellite-based estimates on growing conditions, we assess how data-driven zones affect producer welfare in a simulated area-yield index insurance program. While some data-driven zones modestly improve risk reduction with a fixed number of zones, greater flexibility in zone number and field sampling intensity offers more potential to balance performance and cost. We present a conceptual model and empirical simulations to characterize how zone design and data collection intensity jointly shape insurance outcomes in resource-constrained settings.


**Code and Data Repository**

This repository contains the instructions and scripts used for analysis. The table below summarizes each script and its purpose.

| Script          | Purpose                                                                                                                                         |
|-----------------|-------------------------------------------------------------------------------------------------------------------------------------------------|
| **00_extract_data_gee**      | Extract agri-environmental input data using Google Earth Engine                                                                                 |
| **01_prepData** | Generate all input datasets required for analysis                                                                                               |
| **02_predictYields** | Assess relationship between yields and agri-environmental spatial variables (compute ML models)                                           |
| **03_mergeCluster**  | Perform clustering, transform clustering into insurance zones, visualize process, assess cluster stability, and add/save insurance zone assignment |
| **04_calcDescriptives** | Generate descriptive statistics on input data (crop cut and spatial data), compute specific numbers mentioned in the manuscript, and figures for the appendix |
| **05_mapROI**   | Create map of the study region                                                                                                                  |
| **06_calcCE**   | Simulate various insurance scenarios and assess and visualize results                                                                           |


**Data Privacy Restrictions**

The analysis draws on data provided by an external partner, as detailed in the article. While the data are georeferenced, the exact locations cannot be shared for privacy reasons. As a result, steps that depend on those coordinates—such as spatial data extraction, clustering, and some visualizations—cannot be fully reproduced from these files. However, for transparency, we provide the full analysis code, and anonymized datasets that allow the replication of the insurance simulations and assessments can be provided upon reasonable request. The original dataset may be available upon agreement with the data provider. 


**Instructions**

Implementing this replication package requires access to Google Earth Engine and the statistical software R.  
1. Script 00_gee: This script is a collection of javascript code that extracts the agri-environmental input variables of interest. The code needs to be pasted into the Google Earth Engine Code Editor and run section by section as detailed in the script. The code can be partially replicated without the original, georeferenced dataset. The script has three major objectives: 1) to generate 50k random locations on maize cropland (replicable), 2) to extract agri-environmental data from the 50k locations (also replicable, except for the soil data), and 3) to extract the agri-environmental data associated with the locations of the input dataset (not replicable without the original dataset). Prerequisite to running step 1) is the generation of a Google Earth Engine asset of the study region for which the shapefile is provided in the data folder corresponding to the name of the script. 
2. Scripts 01 - 06: These  R scripts transform and clean the data, simulate and evaluate each insurance scenario  insurance scheme and generate visualizations. 
    - Script 01, 02, and 03 are provided as reference but are not possible to replicate without the original dataset including its locations. 
    - In script 03, we save the data set without the locations so that it can be used in later scripts.
    - Script 04, 05 and 06 can be replicated using the anonymized data (i.e., original data without locations). Minor parts rely on the dataset with the location, which are therefore commented out but still provided for full transparency. 


The R scripts in this repository are designed to be run sequentially, not as standalone files. Script 01 generates datasets used by scripts 02 and 03. Script 03 then produces datasets required for scripts 04, 05, and 06. Upon request, the author team can provide anonymized data files needed to run scripts 04–06.


The table below shows which manuscript and appendix figures/tables can be replicated with the anonymized dataset, and where to find the code to generate them.


| Fig./Tab. | Short description of Fig./Tab. | R Script | Replicable based on anonymized data? |
|-----------|--------------------------------|----------|-------------------------------------|
|*Manuscript*| | | |
| Fig 1     | Map of study region            | 05_mapROI | Yes                                 | 
| Fig 2     | Flowchart outlining the data extraction and clustering (zone-development) protocols | No code, authors own visualization | NA |
| Fig 3     | Example for zone development based on precipitation clusters | 03_mergeCluster | No |
| Fig 4     | Net risk reduction values of insurance relative to no insurance under different (static) zone definitions | 06_calcCE | Yes |
| Fig 5     | (Heatmaps) Net Value by Number of Crop Cuts per Zone and Number of Zones | 06_calcCE | Yes |
| Fig 6     | Slices of heatmaps             | 06_calcCE | Yes |
| Fig 7     | Range of insurance zone variance for two key indicators | 06_calcCE | Yes |
| Table 1   | Descriptive statistics on retained crop cuts | 04_calcDescriptives | Yes | 
| Table 2   | Descriptive statistics on agri-environmental conditions extracted at crop cut locations | 04_calcDescriptives | No |
|*Appendix*| | | |
| Fig A.1   | Distribution of crop cut yields | 04_calcDescriptives | Yes |
| Fig A.2   | Location of data extraction points grouped by grid cell | 03_mergeCluster | No |
| Fig A.3   | Pairwise NMI across insurance zones based on agri-environmental clustering inputs | 03_mergeCluster | No |
| Fig A.4   | Illustration of grid cell clusters and insurance zones by input | 03_mergeCluster | No |
| Fig A.5   | Zoning approaches and alternatives | 04_calcDescriptives | Partially: 1st map not, but 2nd and 3rd map of figure are replicable |
| Fig A.6   | Number of zones considered in analysis vs. number of zones designed for entire study region | 04_calcDescriptives | Yes |
| Fig A.7   | Risk reduction value for 5-100 zones | 06_calcCE | Yes |
| Fig A.8   | Risk reduction value for 5-100 zones assigning higher weight to lower tail of the yield distribution | 06_calcCE | Yes |
| Fig A.9   | Range of insurance zone variance for all indicators | 06_calcCE | Yes |
| Fig A.10  | Boxplots of yields by zone and input | 06_calcCE | Yes |
| Fig A.11  | Dot-Whisker plots on yields by zone and input | 06_calcCE | Yes |
| Tab A.1   | Insurance zone stability by year left out | 03_mergeCluster | No |
| Tab A.2   | Distribution of crop cut locations across data-driven and administrative insurance zones assuming 14 insurance zones spanning the entire study region | 04_calcDescriptives | Yes |
| Tab A.3   | Model performance results | 02_predictYields | No |
