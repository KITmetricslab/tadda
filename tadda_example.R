# install and load packages
install.packages("arrow")
library(arrow) # parquet data
library(dplyr) # data manipulation

# set working directory
current_path = rstudioapi::getActiveDocumentContext()$path # get the path of your current open file
setwd(dirname(current_path))

# read in data
ged_cm_prepatch <- read_parquet('ged_cm_prepatch.parquet') %>% # ? noch nicht aufbereitet
  mutate("fatalities" = rowSums(select(., starts_with("ged_best"))))
ged_cm_postpatch <- read_parquet('ged_cm_postpatch.parquet') %>% # ? aufbereitete Daten?!
  mutate("fatalities" = rowSums(select(., starts_with("ged_best"))))
skeleton_cm_africa <- read_parquet('skeleton_cm_africa.parquet')

log_change_fatalities_postpatch_selected <- list()

par(mfrow=c(6,2)) 

for (i in 1:6){
  # selected countries
  country_name_selected <- c("Mozambique", "Sudan", "Congo, DRC", "Mali", "Nigeria", "Somalia") [i]
  
  # extract data of selected country
  country_id_selected <- unique((skeleton_cm_africa %>%
                                   filter(country_name %in% country_name_selected))$country_id)
  
  ged_cm_prepatch_selected <- ged_cm_prepatch %>%
    filter(country_id %in% country_id_selected)
  ged_cm_postpatch_selected <- ged_cm_postpatch %>%
    filter(country_id %in% country_id_selected)
  n <- length(ged_cm_postpatch_selected$fatalities)
  
  # beginnt erst bei 108
  # endet bei 498
  
  log_change_fatalities_postpatch_selected[[i]] <- log( (ged_cm_postpatch_selected$fatalities[2:n] + 1) / (ged_cm_postpatch_selected$fatalities[1:n-1] + 1 ))
  plot(ged_cm_postpatch_selected$month_id, ged_cm_postpatch_selected$fatalities, xlab = "month ID", ylab = "fatalities per month", xlim = c(100,500), main = country_name_selected)
  plot(ged_cm_postpatch_selected$month_id[2:n], log_change_fatalities_postpatch_selected[[i]], xlab = "month ID", ylab = "log change", xlim = c(100,500), main = country_name_selected)
}
