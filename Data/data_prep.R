# load packages
library(arrow) # parquet data
library(dplyr) # data manipulation

# define function for computing the log change
compute_log_change <- function(vector, step_ahead) {
  n <- length(vector)
  log_change <- log(vector[(step_ahead+1):n] + 1) - log(vector[1:(n-step_ahead)] + 1)
  log_change <- c(rep(NA, step_ahead), log_change)
  return(log_change)
}

# set working directory
current_path = rstudioapi::getActiveDocumentContext()$path # get the path of your current open file
setwd(dirname(current_path))

# read in & extract relevant data
country_name_selected <- c("Mozambique", "Sudan", "Congo, DRC", "Mali", "Nigeria", "Somalia")
skeleton_cm_africa <- read_parquet('skeleton_cm_africa.parquet') %>%
  filter(country_name %in% country_name_selected)

country_id_selected <- unique((skeleton_cm_africa %>%
                                 filter(country_name %in% country_name_selected))$country_id)

ged_cm_postpatch <- read_parquet('ged_cm_postpatch.parquet') %>% 
  mutate("fatalities" = rowSums(select(., starts_with("ged_best")))) %>%
  filter(country_id %in% country_id_selected) %>%
  select(country_id, month_id, fatalities)

# only include observations with fatalities != NA
# discard last 3 months since they weren't included in the requested true future horizon of the competition
months_available <- head(unique(ged_cm_postpatch$month_id[which(!is.na(ged_cm_postpatch$fatalities))]),-3)

month_country_information <-
  skeleton_cm_africa %>% filter(month_id %in% months_available) %>% distinct %>% select(-in_africa)

ged_cm_postpatch <- ged_cm_postpatch %>%
  filter(month_id %in% months_available)

data_fatalities <- merge(month_country_information, ged_cm_postpatch, by = c("country_id", "month_id"))

fatalities <- list()

for (i in 1:length(country_name_selected)) {
  fatalities[[i]] <-
    data_fatalities[which(data_fatalities$country_name == country_name_selected[i]),
                    c("month_id", "month", "year", "fatalities")]
  
  fatalities_log_changes <- lapply(1:7, function(s_ahead) {compute_log_change(fatalities[[i]]$fatalities, s_ahead)})
  names(fatalities_log_changes) <- c("log_change_s1", "log_change_s2", "log_change_s3", "log_change_s4", "log_change_s5", "log_change_s6", "log_change_s7")
  
  fatalities[[i]] <- cbind(fatalities[[i]], as.data.frame(do.call(cbind, fatalities_log_changes)))
  write.csv2(fatalities[[i]], paste("fatalities_", country_name_selected[i], ".csv", sep = "")) # sep = ";"
}

