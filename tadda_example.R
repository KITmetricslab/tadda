### Code for empirical example illustrating properties of forecasts for log changes of armed conflict fatalities using TADDA variants, MAE or MSE
# Lotta Rüter
# lotta.rueter@kit.edu


# load packages
# library(dplyr) # included in tidyverse
library(tidyverse) # required for "reduce" function
library(rlist)

# set working directory
current_path = rstudioapi::getActiveDocumentContext()$path # get path of this file
setwd(dirname(current_path))

# get scoring functions
source("Illustration/functions.R")
source("bayes_acts_functions.R")

# read in data
data_fatalities <- read.csv(paste("Data/data_fatalities.csv"))[,-1]

# set parameters & choose task windows
epsilon <- 0.048
window_length <- 6 # use data of the past 6 months for predictions
month_id_task1 <- 490:495 # true future of challenge (10/20-03/21)
month_id_task2 <- 445:480 # test set (01/17-12/19)
window_months <- 408:495 # data contains obs. for month_id 400 onwards, 7-step ahead forecast exists for 408 onwards

# TADDA-scores that were used in paper: TADDA1 = TADDA1_L1, TADDA2 = TADDA2_L1, each with epsilon = 0.048
BA_names <- c("BA_SE", "BA_TADDA1_L1", "BA_TADDA2_L1")

# initialise objects
predictions_uninformed <- predictions_informed <- predictions_combined <- list()

# compute predictions for each country, score and prediction horizon
for (country in 1:length(unique(data_fatalities$country_name))) {
  current_country_name <- unique(data_fatalities$country_name)[country]; print(current_country_name)
  data_country <- data_fatalities %>% filter(country_name == current_country_name) %>% filter(month_id %in% window_months)
  
  # case 1 (uninformed): use log change distribution directly
  list_uninformed <- lapply(1:7, function(step_ahead) bayes_acts_predictions(variant = "uninformed", s = step_ahead))
  predictions_uninformed[[country]] <- cbind("country_name" = current_country_name, reduce(list_uninformed, full_join, by = "month_id"))
  
  # case 2 (informed): construct the log change distribution based on current observation and past fatalities
  list_informed <- lapply(1:7, function(step_ahead) bayes_acts_predictions(variant = "informed", s = step_ahead))
  predictions_informed[[country]] <- cbind("country_name" = current_country_name, reduce(list_informed, full_join, by = "month_id"))
}

# combine lists to data frame and select months that are relevant for evaluation
data_predictions_uninformed <- merge(data_fatalities, cbind(bind_rows(predictions_uninformed), "No_Change" = 0),
                                     by = c("country_name", "month_id")) %>% filter(month_id %in% 445:495) # 445:495 comprises both evaluation sets
data_predictions_informed <- merge(data_fatalities, cbind(bind_rows(predictions_informed), "No_Change" = 0),
                                     by = c("country_name", "month_id")) %>% filter(month_id %in% 445:495) # 445:495 comprises both evaluation sets

# compute losses for each prediction horizon and loss function
SE_loss_uninformed_df <- loss_df(data_predictions_uninformed, "SE")
TADDA1_L1_loss_uninformed_df <- loss_df(data_predictions_uninformed, "TADDA1_L1")
TADDA2_L1_loss_uninformed_df <- loss_df(data_predictions_uninformed, "TADDA2_L1")

SE_loss_informed_df <- loss_df(data_predictions_informed, "SE")
TADDA1_L1_loss_informed_df <- loss_df(data_predictions_informed, "TADDA1_L1")
TADDA2_L1_loss_informed_df <- loss_df(data_predictions_informed, "TADDA2_L1")

# summarise losses
mean_loss_uninformed_task2 <- mean_loss(data_predictions_uninformed, month_id_task2)
mean_loss_informed_task2 <- mean_loss(data_predictions_informed, month_id_task2)
write.csv(round(mean_loss_uninformed_task2, 3), "mean_loss_uninformed_task2.csv")
write.csv(round(mean_loss_informed_task2, 3), "mean_loss_informed_task2.csv")
