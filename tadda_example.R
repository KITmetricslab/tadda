### Code for empirical example illustrating properties of forecasts for log changes of armed conflict fatalities using TADDA variants, MAE or MSE
# Lotta Rüter
# lotta.rueter@kit.edu

library(dplyr)

# set working directory
current_path = rstudioapi::getActiveDocumentContext()$path # get path of this file
setwd(dirname(current_path))

# get scoring functions
source("Illustration/functions.R")
source("bayes_acts_functions.R")

# read in data
country_name_selected <- c("Mozambique", "Sudan", "Congo, DRC", "Mali", "Nigeria", "Somalia")
i <- 2
data_fatalities <- read.csv2(paste("Data/fatalities_", country_name_selected[i], ".csv", sep = ""))

epsilon <- 0.048
window_length <- 24 # use data of the past 2 years for predictions

# theoretical bayes acts of last 24 observations
# 1 step-ahead forecast
BA_analytical_uninformed <- BA_analytical_informed <- data.frame(matrix(nrow = 0, ncol = 9))
colnames(BA_analytical_uninformed) <- colnames(BA_analytical_informed) <- 
  c("month_id",
    "BA_AE", "BA_TADDA_L1", "BA_TADDA1_L1", "BA_TADDA2_L1",
    "BA_SE", "BA_TADDA_L2", "BA_TADDA1_L2", "BA_TADDA2_L2")

# case 1: using the log change distribution of the window directly (uninformed)
for (i in 1:(length(data_fatalities$fatalities)-window_length-1)) {
  window_begin <- i+1 # first order log-change is only available from t=2 onwards
  window_end <- i+window_length
  log_change_distribution <- data_fatalities$log_change_s1[window_begin:window_end]
  BA_analytical_uninformed[i,] <- c(data_fatalities$month_id[window_end+1], # use bayes acts of past 23 log changes as prediction for observation t+1
                                    compute_bayes_acts(log_change_distribution, epsilon))
}

# case 2: constructing the log change distribution based on current observation and past fatalities (informed)
for (i in 1:(length(data_fatalities$fatalities)-window_length-1)) {
  window_begin <- i
  window_end <- i+window_length
  fatalities_window <- data_fatalities$fatalities[window_begin:(window_end-1)]
  log_change_distribution <- log(fatalities_window+1) - log(data_fatalities$fatalities[window_end]+1) # use last observed value as reference point for log change distribution
  BA_analytical_informed[i,] <- c(data_fatalities$month_id[window_end+1], # prediction for t+1
                                  compute_bayes_acts(log_change_distribution, epsilon))
}

# compute mean losses
y_true <- unlist(data_fatalities %>% filter(month_id %in% BA_analytical_informed$month_id) %>% select(log_change_s1))

loss_tables_uninformed <- loss_tables(BA_analytical_uninformed[,-1], y_true, epsilon)
loss_tables_informed <- loss_tables(BA_analytical_informed[,-1], y_true, epsilon)

mean_loss_table_uninformed <- mean_loss_table(BA_analytical_uninformed[,-1], y_true, epsilon); mean_loss_table_uninformed
mean_loss_table_informed <- mean_loss_table(BA_analytical_informed[,-1], y_true, epsilon); mean_loss_table_informed

## Plots where we can ain't see nothin' yet -> choose different plotting window
## TADDA-scores that were used in paper: TADDA1 = TADDA1_L1, TADDA2 = TADDA2_L1, each with epsilon = 0.048
plot(BA_analytical_informed$month_id, BA_analytical_informed$BA_SE, type = "l", col = "red", xlab = "Month ID", ylab = "log change")
lines(BA_analytical_informed$month_id, BA_analytical_informed$BA_TADDA1_L1, type = "l", col = "blue")
lines(BA_analytical_informed$month_id, BA_analytical_informed$BA_TADDA2_L1, type = "l", col = "green")
lines(BA_analytical_informed$month_id, y_true, type = "l", col = "black")

# plot for SE loss given different predictions optimised for SE, TADDA1_L1 and TADDA2_L1
plot(BA_analytical_informed$month_id, loss_tables_uninformed$SE$BA_SE, type = "l", col = "red", xlab = "Month ID", ylab = "loss")
lines(BA_analytical_informed$month_id, loss_tables_uninformed$SE$BA_TADDA1_L1, type = "l", col = "blue", xlab = "Month ID", ylab = "loss")
lines(BA_analytical_informed$month_id, loss_tables_uninformed$SE$BA_TADDA2_L1, type = "l", col = "green", xlab = "Month ID", ylab = "loss")

# plot for TADDA1_L1 loss given different predictions optimised for SE, TADDA1_L1 and TADDA2_L1
plot(BA_analytical_informed$month_id, loss_tables_uninformed$TADDA1_L1$BA_SE, type = "l", col = "red", xlab = "Month ID", ylab = "loss")
lines(BA_analytical_informed$month_id, loss_tables_uninformed$TADDA1_L1$BA_TADDA1_L1, type = "l", col = "blue", xlab = "Month ID", ylab = "loss")
lines(BA_analytical_informed$month_id, loss_tables_uninformed$TADDA1_L1$BA_TADDA2_L1, type = "l", col = "green", xlab = "Month ID", ylab = "loss")

# plot for TADDA2_L1 loss given different predictions optimised for SE, TADDA1_L1 and TADDA2_L1
plot(BA_analytical_informed$month_id, loss_tables_uninformed$TADDA2_L1$BA_SE, type = "l", col = "red", xlab = "Month ID", ylab = "loss")
lines(BA_analytical_informed$month_id, loss_tables_uninformed$TADDA2_L1$BA_TADDA1_L1, type = "l", col = "blue", xlab = "Month ID", ylab = "loss")
lines(BA_analytical_informed$month_id, loss_tables_uninformed$TADDA2_L1$BA_TADDA2_L1, type = "l", col = "green", xlab = "Month ID", ylab = "loss")

## To do:
# create nice plots
# generalize for multi-step ahead forecast
# check true future






