### Code for empirical example illustrating properties of forecasts for log changes of armed conflict fatalities using TADDA variants, MAE or MSE
# Lotta Rüter
# lotta.rueter@kit.edu

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

# theoretical bayes acts of last 24 observations
# 1 step-ahead forecast

# case 1: using the log change distribution of the window directly 
BA_analytical_uninformed <- BA_analytical_informed <- data.frame(matrix(nrow = 0, ncol = 9))
colnames(BA_analytical_uninformed) <- colnames(BA_analytical_informed) <- 
  c("month_id",
    "BA_AE", "BA_TADDA_L1", "BA_TADDA1_L1", "BA_TADDA2_L1",
    "BA_SE", "BA_TADDA_L2", "BA_TADDA1_L2", "BA_TADDA2_L2")

for (i in 1:(length(data_fatalities$fatalities)-24)) {
  log_change_distribution <- data_fatalities$log_change_s1[(i+1):(i+24)]
  BA_analytical_uninformed[i,] <- c(data_fatalities$month_id[i+24], compute_bayes_acts(log_change_distribution, epsilon))
}


# for (i in 1:(length(data_fatalities$fatalities)-25)) {
#   log_change_distribution <- data_fatalities$log_change_s1[(i+1):(i+24)]
#   BA_analytical_uninformed[i,1] <- BA_AE(log_change_distribution)
#   BA_analytical_uninformed[i,2] <- BA_TADDA_L1(log_change_distribution)
#   BA_analytical_uninformed[i,3] <- BA_TADDA1_L1(log_change_distribution, epsilon)
#   BA_analytical_uninformed[i,4] <- BA_TADDA2_L1(log_change_distribution, epsilon)
#   BA_analytical_uninformed[i,5] <- BA_SE(log_change_distribution)
#   BA_analytical_uninformed[i,6] <- BA_TADDA_L2(log_change_distribution)
#   BA_analytical_uninformed[i,7] <- BA_TADDA1_L2(log_change_distribution, epsilon)
#   BA_analytical_uninformed[i,8] <- BA_TADDA2_L2(log_change_distribution, epsilon)
# }

# case 2: constructing the log change distribution based on current observation and past fatalities
for (i in 1:(length(data_fatalities$fatalities)-24)) {
  fatalities_window <- data_fatalities$fatalities[i:(i+23)]
  log_change_distribution <- log(data_fatalities$fatalities[(i+24)]+1)-log(fatalities_window+1)
  BA_analytical_informed[i,] <- c(data_fatalities$month_id[i+24], compute_bayes_acts(log_change_distribution, epsilon))
}


for (i in 1:(length(data_fatalities$fatalities)-24)) {
  fatalities_window <- data_fatalities$fatalities[i:(i+23)]
  log_change_distribution <- log(data_fatalities$fatalities[(i+24)]+1)-log(fatalities_window+1)
  BA_analytical_informed[i,1] <- BA_AE(log_change_distribution)
  BA_analytical_informed[i,2] <- BA_TADDA_L1(log_change_distribution)
  BA_analytical_informed[i,3] <- BA_TADDA1_L1(log_change_distribution, epsilon)
  BA_analytical_informed[i,4] <- BA_TADDA2_L1(log_change_distribution, epsilon)
  BA_analytical_informed[i,5] <- BA_SE(log_change_distribution)
  BA_analytical_informed[i,6] <- BA_TADDA_L2(log_change_distribution)
  BA_analytical_informed[i,7] <- BA_TADDA1_L2(log_change_distribution, epsilon)
  BA_analytical_informed[i,8] <- BA_TADDA2_L2(log_change_distribution, epsilon)
}




# CHECK AGAIN OR IMPROVE / CORRECT / WHATEVER! -> numerical version

grid_y <- seq(-10, 10, 0.001) # HOW TO CHOOSE GRID OVER Y?
BA <- data.frame()

BA <- c()
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
