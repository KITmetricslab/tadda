
scores <- c(AE, SE, TADDA_L1, TADDA_L2, TADDA_L1_v1, TADDA_L1_v2, TADDA_L2_v1, TADDA_L2_v2)

BA <- data.frame()

BA <- c()
# use rolling windows over last 2 years, i.e. 24 months
for (i in 1:(length(data_fatalities$fatalities)-24)) {
  fatalities_window <- data_fatalities$fatalities[i:(i+23)]
  # for (s in 1:length(scores)){
  log_change_distribution <- log(data_fatalities$fatalities[i+23]+1)-log(fatalities_window+1)
  BA[i] <- get_bayes_acts(log_change_distribution, epsilon, grid_y, TADDA_L2_v2)
  # }
}

# use rolling windows over last 3 years, i.e. 36 months
for (i in 1:(length(data_fatalities$fatalities)-37)) {
  fatalities_window <- data_fatalities$fatalities[i:(i+35)]
  log_change_distribution <- log(data_fatalities$fatalities[(i+36)]+1)-log(fatalities_window+1)
  # BA[i,4] <-BA[i,3] <- BA[i,2] <- BA[i,1]<- NA
  BA[i,1] <- get_bayes_acts_wo_epsilon(log_change_distribution, grid_y, AE)
  BA[i,2] <- get_bayes_acts_wo_epsilon(log_change_distribution, grid_y, SE)
  BA[i,3] <- get_bayes_acts_wo_epsilon(log_change_distribution, grid_y, TADDA_L1)
  BA[i,4] <- get_bayes_acts_wo_epsilon(log_change_distribution, grid_y, TADDA_L2)
  BA[i,5] <- get_bayes_acts(log_change_distribution, epsilon, grid_y, TADDA_L1_v1)
  BA[i,6] <- get_bayes_acts(log_change_distribution, epsilon, grid_y, TADDA_L1_v2)
  BA[i,7] <- get_bayes_acts(log_change_distribution, epsilon, grid_y, TADDA_L2_v1)
  BA[i,8] <- get_bayes_acts(log_change_distribution, epsilon, grid_y, TADDA_L2_v2)
}

# use empirical cdfs of fatalities to compute distribution and hence Bayes acts for TADDA & MSE for each step ahead at each time t
# theoretical points (derived from Bayes acts)
# numerical points (derived via score optimisation)

# mean for each score and forecasting horizon on training set

# score value on test set per score, forecast value per score, and horizon
