### Code for data preparation of empirical example

# Extracts time series of country month fatalities due to state based conflict "ged_best_sb" for each country in Africa
# Computes true s=1,...,7 step ahead log-changes for each country and month
# Saves the results in "fatalities.csv"

# Lotta Rüter
# lotta.rueter@kit.edu

# load packages
library(arrow) # parquet data
library(dplyr) # data manipulation
library(rlist) # functions for manipulation of lists

# define function for computing the log change
compute_log_change <- function(vector, step_ahead) {
  n <- length(vector)
  log_change <- c(rep(NA, step_ahead), log(vector[(step_ahead+1):n] + 1) - log(vector[1:(n-step_ahead)] + 1))
  log_change
}

# set working directory
current_path = rstudioapi::getActiveDocumentContext()$path # get the path of your current open file
setwd(dirname(current_path))

# read in & extract relevant data of all African countries
skeleton_cm_africa <- read_parquet('skeleton_cm_africa.parquet') # contains data on only African countries
ged_cm_postpatch <- read_parquet('ged_cm_postpatch.parquet') %>% 
  filter(country_id %in% skeleton_cm_africa$country_id) %>%
  select(country_id, month_id, ged_best_sb) %>%
  rename("fatalities" = "ged_best_sb") %>%
  arrange(month_id) %>%
  arrange(country_id)

months_used <- 121:495

month_country_information <-
  skeleton_cm_africa %>% filter(month_id %in% months_used) %>% distinct %>% select(-in_africa)

ged_cm_postpatch <- ged_cm_postpatch %>%
  filter(month_id %in% months_used)

data_merged <- merge(month_country_information, ged_cm_postpatch, by = c("country_id", "month_id"))

n <- nrow(data_merged)
n_month <- length(months_used)

data_fatalities <- data.frame(matrix(ncol = 13, nrow = n))
colnames(data_fatalities) <- c("country_name", "country_id", "month_id", "month", "year", "fatalities", "log_change_s1", "log_change_s2", "log_change_s3", "log_change_s4", "log_change_s5", "log_change_s6", "log_change_s7")
data_fatalities[,1:6] <- data_merged[,c("country_name", "country_id", "month_id", "month", "year", "fatalities")]

# compute log-changes
for (i in 1:length(unique(data_merged$country_name))) {
  country <- unique(data_merged$country_name)[i]
  data_country <- data_merged %>% filter(country_name == country)
  log_changes <- lapply(1:7, function(s) {compute_log_change(data_country$fatalities, s)})
  data_fatalities[which(data_fatalities$country_name == country), 7:13] <- list.cbind(log_changes)
}

data_fatalities <- data_fatalities # %>% filter(month_id >= 426) # observation 445-12-7 required for log-change distributions of s7-step ahead forecast of uninformed variant
write.csv(data_fatalities, "fatalities.csv")
