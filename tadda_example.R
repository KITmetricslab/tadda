### Code for empirical example in Section 5
# Lotta Rüter
# lotta.rueter@kit.edu

### Part 1: Determines optimal window length based on averages losses on training set of task 2
# Produces "results/average_scores_for_different_window_lengths.csv" which shows that w=5 would be optimal for minimizing TADDA1 via TADDA1_OPF

### Part 2: Yields results for chosen window length w = 9 for all of Africa
# Computes the predictions for the log-changes in fatalities and corresponding losses for each African country, month in 395:495, OPF and lead time s=2,...,7, see "results/individual_predictions_w9.csv" and "results/individual_losses_w9.csv"
# Produces the central results presented in Table 2, see "results/average_scores_w9.csv"
# Produces the basis of Table 3, see "results/empirical_quantiles_w9.csv"

# ---------------
### Part 0: Initialization
# ---------------
# load packages
library(tidyverse) # includes dplyr, required for "reduce" function
library(rlist) # functions for manipulation of lists

# set working directory
current_path = rstudioapi::getActiveDocumentContext()$path # get path of this file
setwd(dirname(current_path))

# get scoring functions
source("Theory/functions.R")
source("Theory/bayes_acts_functions.R")

# read in data
data_fatalities <- read.csv(paste("Data/fatalities.csv"))[,-1]
country_names <- unique(data_fatalities$country_name)

# set parameters & choose task windows
epsilon <- 0.048

# task 2
train_month_id_task2 <- 409:444 # training set (01/14-12/16)
pred_month_id_task2 <- 445:480 # test set (01/17-12/19)

# overall data window
# data generally contains obs. for month_id 121 onwards, 7-step ahead forecast exists for 128 onwards; 
# however, data for South Sudan (incl. 7 step-ahead forecast) is only available from 386 onwards so we choose 386 as the starting month
window_months <- 386:495 

# TADDA-score that was used in paper: TADDA1 = TADDA1_L1 with epsilon = 0.048
OPF_names <- c("OPF_SE", "OPF_TADDA1")

# initialise objects
predictions <- list()

# mean scores of VIEWS competition ensemble (from Table 2 in Vesco et al., 2022)
MSE_ensemble_cm <- c(.504, .551, .579, .548, .573, .599, mean(c(.504, .551, .579, .548, .573, .599)))
TADDA_ensemble_cm <- c(.371, .379, .394, .381, .386, .400, mean(c(.371, .379, .394, .381, .386, .400)))
# ---------------


# ---------------
### Part 1: Determine optimal window length based on averages losses on training set of task 2
# ---------------
loss_means_df_train <- loss_means_df_pred <- data.frame(matrix(NA, ncol = 6, nrow = 17))
colnames(loss_means_df_train) <- colnames(loss_means_df_pred) <- c("MSE_OPF_SE", "MSE_OPF_TADDA1", "MSE_No_Change",
                                                                   "MTADDA1_OPF_SE", "MTADDA1_OPF_TADDA1", "MTADDA1_No_Change")

for(window_length in 1:17) { # 17 is highest possible window length, since 409-17-7+1 = 386, the first month for which observations of all countries are available
  print(paste("window_length =", window_length))
  
  # compute predictions for each country, score and prediction horizon
  for (country in 1:length(country_names)) {
    current_country_name <- country_names[country] # print(current_country_name)
    data_country <- data_fatalities %>% filter(country_name == current_country_name) %>% filter(month_id %in% window_months)
    
    # construct the log change distribution based on current observation and past fatalities
    predictions_ls <- lapply(1:7, function(step_ahead) bayes_acts_predictions(s = step_ahead))
    predictions[[country]] <- cbind("country_name" = current_country_name, reduce(predictions_ls, full_join, by = "month_id"))
  }
  
  # combine list to data frame and select months that are relevant for evaluation
  # predictions are available for (window_months[1] + window_length + s) : window
  data_predictions <- merge(data_fatalities, cbind(bind_rows(predictions), "No_Change" = 0), by = c("country_name", "month_id")) %>% 
    select(., c("country_name", "month_id", starts_with("log_change"), starts_with("OPF_SE"), starts_with("OPF_TADDA1"), "No_Change"))
  
  # compute losses for each prediction horizon and loss function
  loss_all <- compute_losses(data_predictions)
  
  # summarize losses
  mean_loss_task2_train <- mean_loss(loss_all, train_month_id_task2) # task 2
  mean_loss_task2_pred <- mean_loss(loss_all, pred_month_id_task2) # task 2
  
  loss_means_df_train[window_length,] <- mean_loss_task2_train[7,] # stores average loss across window
  loss_means_df_pred[window_length,] <- mean_loss_task2_pred[7,] # stores average loss across window
}
round(loss_means_df_pred[1:17,], 3)
write.csv(round(loss_means_df_train[1:17,], 3), "average_scores_for_different_window_lengths.csv")
# ---------------


# ---------------
### Part 2: Results for chosen window w = 9
# ---------------
window_length <- 9 # use data of the past 9 months for predictions

# compute predictions for each country, score and prediction horizon
for (country in 1:length(country_names)) {
  current_country_name <- country_names[country] # print(current_country_name)
  data_country <- data_fatalities %>% filter(country_name == current_country_name) %>% filter(month_id %in% window_months)
  
  # construct the log change distribution based on current observation and past fatalities
  predictions_ls <- lapply(1:7, function(step_ahead) bayes_acts_predictions(s = step_ahead))
  predictions[[country]] <- cbind("country_name" = current_country_name, reduce(predictions_ls, full_join, by = "month_id"))
}

# combine list to data frame and select months that are relevant for evaluation
# predictions are available for (window_months[1] + window_length + s) : window
data_predictions <- merge(data_fatalities, cbind(bind_rows(predictions), "No_Change" = 0), by = c("country_name", "month_id")) %>% 
  select(., c("country_name", "month_id", starts_with("log_change"), starts_with("OPF_SE"), starts_with("OPF_TADDA1"), "No_Change"))
write.csv(data_predictions, paste("individual_predictions_w", window_length, ".csv", sep = ""))

## 2.1 Summary statistics
data_OPF_TADDA <- unlist(data_predictions %>% filter(month_id %in% pred_month_id_task2) %>%
                           select(., c(starts_with("OPF_TADDA1"))))
data_OPF_SE <- unlist(data_predictions %>% filter(month_id %in% pred_month_id_task2) %>%
                        select(., c(starts_with("OPF_SE"))))
data_log_change <- unlist(data_predictions %>% filter(month_id %in% pred_month_id_task2) %>%
                            select(., c(starts_with("log_change"))))

empirical_quantiles <- rbind(round(quantile(data_OPF_SE, (c(1:4, 15:20))/20), 3),
                             round(quantile(data_OPF_TADDA, (c(1:4, 15:20))/20), 3),
                             round(quantile(data_log_change, (c(1:4, 15:20))/20), 3))
rownames(empirical_quantiles) <- c(OPF_names, "True_log_changes")
write.csv2(empirical_quantiles, paste("empirical_quantiles_w", window_length, ".csv", sep = ""))

# the following characteristics are computed across all forecasting horizons s=2,...,7
mean(data_OPF_SE) # average mean forecast
sum(data_OPF_SE==0)/length(data_OPF_SE) # percentage of mean forecasts that are zero
sum(abs(data_OPF_SE)>epsilon)/length(data_OPF_SE) # percentage of mean forecasts outside the epsilon interval

mean(data_OPF_TADDA) # average OPF_TADDA forecast
sum(data_OPF_TADDA==0)/length(data_OPF_TADDA) # percentage of OPF_TADDA forecasts that are zero
sum(abs(data_OPF_TADDA)>epsilon)/length(data_OPF_TADDA) # percentage of OPF_TADDA forecasts outside the epsilon interval

mean(data_log_change) # average true log change
sum(data_log_change==0)/length(data_log_change) # percentage of true log changes that are zero
sum(abs(data_log_change)>epsilon)/length(data_log_change) # percentage of true log_changes outside the epsilon interval

## 2.2
# compute losses for each prediction horizon and loss function
loss_all <- compute_losses(data_predictions)
write.csv(loss_all, paste("individual_losses_w", window_length, ".csv", sep = ""))

# summarise losses
mean_loss_task2_pred <- mean_loss(loss_all, pred_month_id_task2); mean_loss_task2_pred # task 2
mean_loss_task2_pred_incl_ensembles <- data.frame(mean_loss_task2_pred[,1:2], MSE_ensemble_cm, mean_loss_task2_pred[,3:5], TADDA_ensemble_cm, "MTADDA1_No_Change" = mean_loss_task2_pred[,6])
write.csv(round(mean_loss_task2_pred_incl_ensembles, 3), paste("average_scores_w", window_length, ".csv", sep = ""))
