### Functions for computing the optimal point forecasts (OPFs)
# of different scoring functions (AE, SE and variants of TADDA)
# as well as a summary function for a nicer represenation of the results

# Lotta Rüter
# lotta.rueter@kit.edu

# ----
### Optimal point forecasts (OPFs)
# ----
type_no <- 1 # calculation type for quantile estimation

## L1 (Absolute Distance) Optimal Point Forecasts
OPF_AE <- function(distribution_Y) {
  median(distribution_Y, type = type_no)
}

OPF_TADDA1_L1 <- function(distribution_Y, epsilon){
  pi_minus <- mean(distribution_Y < -epsilon)
  pi_plus <- mean(distribution_Y > epsilon)
  m <- median(distribution_Y, type = type_no)
  
  if(pi_minus >= 0.5*(1 + pi_plus)){
    return(quantile(distribution_Y, 0.5*(1 + pi_plus), type = type_no))
  }
  
  if(pi_minus > 0.5 & pi_minus < 0.5*(1 + pi_plus)){
    return(-epsilon)
  }
  
  if(pi_minus <= 0.5 & pi_plus <= 0.5){
    return(m)
  }
  
  if(pi_plus > 0.5 & pi_plus <= 0.5*(1 + pi_minus)){
    return(epsilon)
  }
  
  if(pi_plus > 0.5*(1 + pi_minus)){
    return(quantile(distribution_Y, 0.5*(1 - pi_minus), type = type_no))
  }
}

OPF_TADDA2_L1 <- function(distribution_Y, epsilon){
  pi_minus <- mean(distribution_Y < -epsilon)
  pi_plus <- mean(distribution_Y > epsilon)
  Pr_minus_epsilon <- mean(distribution_Y == epsilon)
  m <- median(distribution_Y, type = type_no)
  
  if(pi_minus >= 2/3){
    return(quantile(distribution_Y, 0.5*(2 - pi_minus), type = type_no))
  }
  
  if(m < -epsilon & pi_minus < 2/3){
    return(-epsilon)
  }
  
  if(-epsilon <= m & m <= epsilon & pi_minus >= (1 + pi_plus - 2*Pr_minus_epsilon)/3){
    return(-epsilon)
  }
  
  if(-epsilon <= m & m <= epsilon & pi_plus < (1 + pi_minus)/3 & pi_minus < (1 + pi_plus - 2*Pr_minus_epsilon)/3){
    return(quantile(distribution_Y, 0.5*(1 - pi_minus + pi_plus), type = type_no))
  }
  
  if(-epsilon <= m & m <= epsilon & -epsilon < m & m < epsilon & pi_plus >= (1 + pi_minus)/3){
    return(epsilon)
  }
  
  if(m > epsilon & pi_plus <= 2/3){
    return(epsilon)
  }
  
  if(pi_plus > 2/3){
    return(quantile(distribution_Y, 0.5*pi_plus), type = type_no)
  }
}

## L2 (Squared Distance) Optimal Point Forecasts
OPF_SE <- function(distribution_Y) {
  mean(distribution_Y)
}

OPF_TADDA1_L2 <- function(distribution_Y, epsilon) {
  pi_minus <- mean(distribution_Y < -epsilon)
  pi_plus <- mean(distribution_Y > epsilon)
  mu <- mean(distribution_Y)
  
  if(mu < -epsilon){
    return(mu*(1/(1 + pi_plus)) - epsilon*(pi_plus/(1 + pi_plus)))
  }
  
  if(mu >= -epsilon & mu <= epsilon){
    return(mu)
  }
  
  if(mu > epsilon){
    return(mu*(1/(1 + pi_minus)) + epsilon*(pi_minus/(1 + pi_minus)))
  }
}

# ----
### Summary Functions
# ----

# bayes acts of SE, TADDA1_L1, TADDA2_L1
bayes_acts <- function(distribution_Y, epsilon) {
  data.frame(
    "OPF_SE" = OPF_SE(distribution_Y),
    "OPF_TADDA1" = OPF_TADDA1_L1(distribution_Y, epsilon)
  )
}

# Predictions for rolling windows
bayes_acts_predictions <- function(data_c = data_country, s = step_ahead, w_length = window_length, eps = epsilon) {
  predictions <- lapply(1:(nrow(data_c)-w_length-s+1), function(w_begin) {
    w_end <- w_begin + w_length - 1
    fatalities_w <- data_c$fatalities[w_begin:w_end]
    log_change_distribution <- log(fatalities_w + 1) - log(tail(fatalities_w, 1) + 1) # use last observed value as reference point for log change distribution
    preds <- data.frame("month_id" = data_c$month_id[w_end + s], bayes_acts(log_change_distribution, eps))
    colnames(preds) <- c("month_id", paste(colnames(preds)[2:ncol(preds)], "_s", s, sep = ""))
    rownames(preds) <- NULL
    preds
  })
  bind_rows(predictions)
}

# Computes losses with loss function "loss" for all forecasts of forecasting horizon s
compute_loss_per_s_and_loss_fun <- function(s, loss, data_all, eps = epsilon, bayes_act_names = OPF_names) {
  true_data_colname <- paste("log_change_s", s, sep = "")
  predictions_colnames <- c(paste(bayes_act_names, s, sep = "_s"), "No_Change")
  if(loss == "SE") {
    df_losses <- SE(data_all[, predictions_colnames], data_all[, true_data_colname])
  } else if(loss == "TADDA1") {
    df_losses <- list.cbind(lapply(predictions_colnames, function(colname) TADDA_L1_v1(data_all[, colname], data_all[, true_data_colname], eps)))
  } else { stop("'loss' must bei either 'SE', 'TADDA1'")  }
  colnames(df_losses) <- predictions_colnames
  colnames(df_losses)[ncol(df_losses)] <- paste(colnames(df_losses)[ncol(df_losses)], "_s", s, sep = "")
  df_losses <- data.frame("country_name" = data_all$country_name, "month_id" = data_all$month_id, df_losses)
  df_losses
}

compute_losses <- function(data_predictions) {
  loss_SE_ls <- lapply(1:7, function(s) compute_loss_per_s_and_loss_fun(s, loss = "SE", data_all = data_predictions))
  loss_TADDA1_ls <- lapply(1:7, function(s) compute_loss_per_s_and_loss_fun(s, loss = "TADDA1", data_all = data_predictions))
  
  loss_SE <- list.cbind(loss_SE_ls)[, !duplicated(colnames(list.cbind(loss_SE_ls)))]
  loss_TADDA1 <- list.cbind(loss_TADDA1_ls)[, !duplicated(colnames(list.cbind(loss_TADDA1_ls)))]
  
  colnames(loss_SE)[3:length(colnames(loss_SE))] <- paste("loss_SE_", colnames(loss_SE)[3:length(colnames(loss_SE))], sep = "")
  colnames(loss_TADDA1)[3:length(colnames(loss_TADDA1))] <- paste("loss_TADDA1_", colnames(loss_TADDA1)[3:length(colnames(loss_TADDA1))], sep = "")
  
  loss_SE_ordered <- loss_SE %>% select(., "country_name", "month_id", starts_with("loss_SE_OPF_SE"), starts_with("loss_SE_OPF_TADDA1"), starts_with("loss_SE_No_Change"))
  loss_TADDA1_ordered <- loss_TADDA1 %>% select(., "country_name", "month_id", starts_with("loss_TADDA1_OPF_SE"), starts_with("loss_TADDA1_OPF_TADDA1"), starts_with("loss_TADDA1_No_Change"))
  
  merge(loss_SE_ordered, loss_TADDA1_ordered, by = c("country_name", "month_id"))
}

mean_loss <- function(loss, month_id_task) {
  mloss <- colMeans(loss[which(loss$month_id %in% month_id_task), ][-c(1,2)])
  mloss_ordered <- matrix(mloss, nrow = 7, byrow = FALSE)
  mloss_ordered <- mloss_ordered[-1,]
  mloss_ordered <- rbind(mloss_ordered, colMeans(mloss_ordered))
  
  colnames(mloss_ordered) <- c("MSE_OPF_SE", "MSE_OPF_TADDA1", "MSE_No_Change",
                               "MTADDA1_OPF_SE", "MTADDA1_OPF_TADDA1", "MTADDA1_No_Change")
  
  rownames(mloss_ordered) <- c(paste("s", 2:7, sep = ""), "colmean")
  
  # round(mloss_ordered, 4)
  mloss_ordered
}
