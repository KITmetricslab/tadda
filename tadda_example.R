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
data_fatalities <- read.csv(paste("Data/fatalities.csv"))[,-1]
country_names <- unique(data_fatalities$country_name)

# set parameters & choose task windows
epsilon <- 0.048
window_length <- 10 # use data of the past 9 months for predictions

# task 1
pred_month_id_task1 <- 490:495 # true future of challenge (10/20-03/21)
train_month_id_task1 <- 445:488

# task 2
pred_month_id_task2 <- 445:480 # test set (01/17-12/19)
train_month_id_task2 <- 409:444 

# overall data window
# data generally contains obs. for month_id 121 onwards, 7-step ahead forecast exists for 128 onwards; 
# however, data for South Sudan (incl. 7 step-ahead forecast) is only available from 386 onwards so we choose 386 as the starting month
window_months <- 386:495 

# TADDA-score that was used in paper: TADDA1 = TADDA1_L1 with epsilon = 0.048
BA_names <- c("BA_SE", "BA_TADDA1", "BA_TADDA2")

# initialise objects
predictions <- list()

# compute predictions for each country, score and prediction horizon
for (country in 1:length(country_names)) {
  current_country_name <- country_names[country]; print(current_country_name)
  data_country <- data_fatalities %>% filter(country_name == current_country_name) %>% filter(month_id %in% window_months)
  
  # construct the log change distribution based on current observation and past fatalities
  predictions_ls <- lapply(1:7, function(step_ahead) bayes_acts_predictions(s = step_ahead))
  predictions[[country]] <- cbind("country_name" = current_country_name, reduce(predictions_ls, full_join, by = "month_id"))
}

# combine list to data frame and select months that are relevant for evaluation
# predictions are available for (window_months[1] + window_length + s) : window
data_predictions <- merge(data_fatalities, cbind(bind_rows(predictions), "No_Change" = 0), by = c("country_name", "month_id")) %>% 
  select(., c("country_name", "month_id", starts_with("log_change"), starts_with("BA_SE"), starts_with("BA_TADDA1"), starts_with("BA_TADDA2"), "No_Change"))
# write.csv(data_predictions, paste("predictions_w", window_length, ".csv", sep = ""))

# compute losses for each prediction horizon and loss function
loss_all <- compute_losses(data_predictions)
# write.csv(loss_all, paste("losses_w", window_length, ".csv", sep = ""))

# summarise losses
mean_loss_task2 <- mean_loss(loss_all, pred_month_id_task2); mean_loss_task2 # task 2
# write.csv(mean_loss_task2, paste("mean_loss_task2_w", window_length, ".csv", sep = ""))


