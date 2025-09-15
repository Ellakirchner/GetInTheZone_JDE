#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Author: Andrew Hobbs 
# Project name: RS Get in the Zone
# Date Last Updated: # Jul 3 2025 (by Andrew)
# Purpose: Assess relation between yields and agri-environmental spatial variables
# Input: 
        # - raw oaf data 
        # - agri-environmental input data for all indicators (later filtered to relevant years)
# Output Files: Table on measures of fit of four models
# ReadMe: 1) Adapt filepath in line 21
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# %%%%%%%%%%%%%%%%%%%%%%%%%%%
# Set Up ----
# %%%%%%%%%%%%%%%%%%%%%%%%%%%

# Clean out workspace
rm(list = ls(all = TRUE))

# Set Up File path & directory: 
file_path <- "put_your_file_path_here" #File path to project folder

#setwd(file_path)
getwd()

setwd(file_path)

# Set Seed for replication 
set.seed(123456789)

# %%%%%%%%%%%%%%%%%%%%%%%%%%%
# Load libraries ----
# %%%%%%%%%%%%%%%%%%%%%%%%%%%
library(dplyr)
library(tidyr)
library(stringr)
library(caret)
library(dplyr)
library(stringr)
library(randomForest)
library(xgboost)
library(zoo)
library(doParallel)
library(sf)  
library(purrr)
library(xtable)


# %%%%%%%%%%%%%%%%%%%%%%%%%%%
# Load data ----
# %%%%%%%%%%%%%%%%%%%%%%%%%%%
load("datasets/oaf_raw.RData")
yields <- oaf_raw
load('datasets/all_oaf.Rdata') 

#%%%%%%%%%%%%%%%%%%%%%%%%%
# 1 Prepare data ----
#%%%%%%%%%%%%%%%%%%%%%%%%%

# 1.1 - Prep agri-envir. var ----
# Define relevant years first 
# (We are interested in keeping those years for which we have yield data 
# to compare the association between the agri-environmental var and the 
# observed yields)
years <- c("2016", "2017", "2018", "2019", "2020")

# Long-format concatenation
combined_data <- map_dfr(years, function(kw) {
  all_oaf %>%
    select(unique_id, contains(kw) | contains("soil")) %>%
    rename_with(~ sub(kw, "", .), contains(kw)) %>%
    mutate(year = kw)
}) %>% 
# drop geometry  
st_drop_geometry()

combined_data <- combined_data %>% merge(yields, on=c("unique_id", "year"))

# combine similar columns and drop other columns that we don't have data every year for
combined_data <- combined_data %>%
  mutate(`ndvi_-05-17` = coalesce(`ndvi_-05-17`, `ndvi_-05-16`),
         `ndvi_-06-02` = coalesce(`ndvi_-06-02`, `ndvi_-06-01`),
         `ndvi_-06-18` = coalesce(`ndvi_-06-18`, `ndvi_-06-17`),
         `ndvi_-07-04` = coalesce(`ndvi_-07-04`, `ndvi_-07-03`),
         `ndvi_-07-20` = coalesce(`ndvi_-07-20`, `ndvi_-07-19`),
         `ndvi_-08-05` = coalesce(`ndvi_-08-05`, `ndvi_-08-04`),
         `ndvi_-08-21` = coalesce(`ndvi_-08-21`, `ndvi_-08-20`)) %>%
  mutate(`evi_-05-17` = coalesce(`evi_-05-17`, `evi_-05-16`),
         `evi_-06-02` = coalesce(`evi_-06-02`, `evi_-06-01`),
         `evi_-06-18` = coalesce(`evi_-06-18`, `evi_-06-17`),
         `evi_-07-04` = coalesce(`evi_-07-04`, `evi_-07-03`),
         `evi_-07-20` = coalesce(`evi_-07-20`, `evi_-07-19`),
         `evi_-08-05` = coalesce(`evi_-08-05`, `evi_-08-04`),
         `evi_-08-21` = coalesce(`evi_-08-21`, `evi_-08-20`)) %>% 
  # drop remaining columns with NA values
  select(where(~ all(!is.na(.))))


# 2 Set up parallel backend ----
# Check number of cores on your machine and set it to two less than your available 
num_cores <- detectCores()
cl <- makeCluster(num_cores)
registerDoParallel(cl)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 3 Prepare models ----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# 3.1 - Stratified train/test split ----
train_indices <- combined_data %>%
  group_by(year) %>%
  group_split() %>%
  lapply(function(df) sample(nrow(df), size = floor(0.7 * nrow(df)))) %>%
  unlist() %>%
  as.integer()

train_data <- combined_data[train_indices, ]
test_data  <- combined_data[-train_indices, ]

# 3.2 - Feature selection ----
feature_cols <- names(combined_data)[str_detect(names(combined_data), "temp|chirp|ndvi|evi|soil")]
target_col <- "yield_kg_ph"

train_data <- train_data %>% select(all_of(c(feature_cols, target_col))) %>% na.omit()
test_data  <- test_data %>% select(all_of(c(feature_cols, target_col))) %>% na.omit()

x_train <- train_data[, feature_cols]
y_train <- train_data[[target_col]]
x_test  <- test_data[, feature_cols]
y_test  <- test_data[[target_col]]

# 3.3 - Training control ----
ctrl <- trainControl(method = "cv", number = 5, allowParallel = TRUE, verboseIter = TRUE)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 4 Train models ----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 4.1 - Linear Regression  ----
lm_model <- train(x = x_train, y = y_train, method = "lm", trControl = ctrl)

# 4.2 - Elastic Net ----
glmnet_grid <- expand.grid(
  alpha = seq(0, 1, by = 0.05),  
  lambda = 10^seq(-5, 0, length.out = 50)  
)

glmnet_model <- train(
  x = x_train,
  y = y_train,
  method = "glmnet",
  tuneGrid = glmnet_grid,
  trControl = ctrl,
  preProcess = c("center", "scale")  
)

# 4.3 - Random Forest ----
rf_grid <- expand.grid(
  mtry = unique(floor(seq(1, length(feature_cols), length.out = 10)))
)

rf_model <- train(
  x = x_train,
  y = y_train,
  method = "rf",
  tuneGrid = rf_grid,
  trControl = ctrl,
  ntree = 750,  
  importance = TRUE
)

# 4.4 - XGBoost  ----
xgb_grid <- expand.grid(
  nrounds = c(100, 200, 300),
  max_depth = c(3, 6, 9),
  eta = c(0.01, 0.1, 0.3),
  gamma = c(0, 1, 5),
  colsample_bytree = c(0.6, 0.8, 1),
  min_child_weight = c(1, 5),
  subsample = c(0.6, 0.8, 1)
)

xgb_model <- train(
  x = as.matrix(x_train),
  y = y_train,
  method = "xgbTree",
  trControl = ctrl,
  tuneGrid = xgb_grid,
  verbosity = 0
)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 5 Evaluate models ---- 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%

lm_preds  <- predict(lm_model,  newdata = x_test)
glmnet_preds <- predict(glmnet_model, newdata = x_test)
rf_preds  <- predict(rf_model,  newdata = x_test)
xgb_preds <- predict(xgb_model, newdata = as.matrix(x_test))

results <- data.frame(
  Model = c("Linear Regression", "Elastic Net", "Random Forest", "XGBoost"),
  R2 = c(
    R2(lm_preds, y_test),
    R2(glmnet_preds, y_test),
    R2(rf_preds, y_test),
    R2(xgb_preds, y_test)
  )
)

print(results)

# Show best parameters
print(rf_model$bestTune)
print(xgb_model$bestTune)

# Clean up parallel backend
stopCluster(cl)
registerDoSEQ()

xtable(results, digits=4, 
       caption = "Model Performance Results",
       label = "tab:model_performance_results",
       align = "lcc") %>%
  print(include.rownames = FALSE, 
        booktabs = TRUE, 
        floating = TRUE, 
        size = "small", 
        caption.placement = "top") 

#%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 6 Save models ----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Create a folder if it doesn't exist
dir.create("models", showWarnings = FALSE, recursive = TRUE)

# Save each model
saveRDS(lm_model,    file = "models/lm_model.rds")
saveRDS(glmnet_model,file = "models/glmnet_model.rds")
saveRDS(rf_model,    file = "models/rf_model.rds")
saveRDS(xgb_model,   file = "models/xgb_model.rds")