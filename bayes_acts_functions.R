### Functions for Bayes Acts of different scoring functions of AE, SE and TADDA
# Lotta Rüter
# lotta.rueter@kit.edu

# ----
### BAYES ACTS
# ----
## L1 Bayes Acts
BA_AE <- function(distribution_Y) {
  median(distribution_Y, type = 1)
}

BA_TADDA_L1 <- function(distribution_Y) {
  m <- median(distribution_Y, type = 1)
  pi <- sum(sign(m) != sign(distribution_Y))/length(distribution_Y)
  ifelse(pi<1/3, quantile(distribution_Y, 0.5 * (1 - sign(m)*pi)), 0)
}

BA_TADDA1_L1 <- function(distribution_Y, epsilon) {
  m <- median(distribution_Y, type = 1)
  
  if(m < (-epsilon)) {
    F_minus_epsilon <- ecdf(distribution_Y)(-epsilon)
    F_epsilon <- ecdf(distribution_Y)(epsilon)
    critical_value <- F_minus_epsilon / (2 - F_epsilon)
    if (critical_value > 0.5) {
      return(quantile(distribution_Y, 0.5*(2-F_epsilon), type = 1))
    } else {
      return(-epsilon)
    }
  }
  
  else if(m <= epsilon) {
    return(m)
  }
  
  else {
    F_minus_epsilon <- ecdf(distribution_Y)(-epsilon)
    F_epsilon <- ecdf(distribution_Y)(epsilon)
    critical_value <- (F_epsilon - F_minus_epsilon) / (1 + F_minus_epsilon)
    if (critical_value >= 0.5) {
      return(epsilon)
    } else {
      return(quantile(distribution_Y, 0.5*(1-F_minus_epsilon), type = 1))
    }
  }
}

BA_TADDA2_L1 <- function(distribution_Y, epsilon) {
  m <- median(distribution_Y, type = 1)
  
  # case m < -epsilon
  if(m < (-epsilon)) {
    F_minus_epsilon <- ecdf(distribution_Y)(-epsilon)
    critical_value <- F_minus_epsilon / (2 - F_minus_epsilon)
    if (critical_value > 0.5) {
      return(quantile(distribution_Y, 0.5*(2-F_minus_epsilon), type = 1))
    } else {
      return(-epsilon)
    }
  }
  
  else if(m <= epsilon) {
    return(m)
  }
  
  else {
    F_epsilon <- ecdf(distribution_Y)(epsilon)
    critical_value <- 2 * F_epsilon / (1 + F_epsilon)
    if (critical_value >= 0.5) {
      return(epsilon)
    } else {
      return(quantile(distribution_Y, 0.5*(1 - F_epsilon), type = 1))
    }
  }
}

## L2 Bayes Acts
BA_SE <- function(distribution_Y) {
  mean(distribution_Y)
}

BA_TADDA_L2 <- function(distribution_Y) {
  mu <- mean(distribution_Y)
  pi <- sum(sign(mu) != sign(distribution_Y))/length(distribution_Y)
  mu / (1+pi)
}

BA_TADDA1_L2 <- function(distribution_Y, epsilon) {
  mu <- mean(distribution_Y)
  
  if(mu < (-epsilon)) {
    F_epsilon <- ecdf(distribution_Y)(epsilon)
    denominator <- 2 - F_epsilon
    return(mu/denominator
           - epsilon * (1-F_epsilon)/denominator)
  }
  
  else if(mu <= epsilon){
    return(mu)
  }
  
  else {
    F_minus_epsilon <- ecdf(distribution_Y)(-epsilon)
    denominator <- 1 + F_minus_epsilon
    return(mu/denominator
           + epsilon * F_minus_epsilon/denominator)
  }
}

BA_TADDA2_L2 <- function(distribution_Y, epsilon) {
  mu <- mean(distribution_Y)
  
  if(mu < (-epsilon)) {
    F_minus_epsilon <- ecdf(distribution_Y)(-epsilon)
    F_epsilon <- ecdf(distribution_Y)(epsilon)
    denominator <- 2 - F_minus_epsilon
    return(mu/denominator
           - epsilon * (F_epsilon-F_minus_epsilon)/denominator
           + epsilon * (1-F_epsilon)/denominator)
  }
  
  else if(mu <= epsilon){
    return(mu)
  }
  
  else {
    F_minus_epsilon <- ecdf(distribution_Y)(-epsilon)
    F_epsilon <- ecdf(distribution_Y)(epsilon)
    denominator <- 1 + F_epsilon
    
    return(mu/denominator
           - epsilon * F_minus_epsilon/denominator
           + epsilon * (F_epsilon-F_minus_epsilon)/denominator)
  }
}

# ----
### SUMMARY FUNCTIONS
# ----

# bayes acts of SE, TADDA1_L1, TADDA2_L1
bayes_acts <- function(distribution_Y, epsilon) {
  data.frame(
    "BA_SE" = BA_SE(distribution_Y),
    "BA_TADDA1" = BA_TADDA1_L1(distribution_Y, epsilon),
    "BA_TADDA2" = BA_TADDA2_L1(distribution_Y, epsilon)
  )
}

# predictions for rolling windows
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

# computes losses with loss function "loss" for all forecasts of forecasting horizon s
compute_loss_per_s_and_loss_fun <- function(s, loss, data_all, eps = epsilon, bayes_acts = BA_names) {
  true_data_colname <- paste("log_change_s", s, sep = "")
  predictions_colnames <- c(paste(bayes_acts, s, sep = "_s"), "No_Change")
  if(loss == "SE") {
    df_losses <- SE(data_all[, predictions_colnames], data_all[, true_data_colname])
  } else if(loss == "TADDA1") {
    df_losses <- list.cbind(lapply(predictions_colnames, function(colname) TADDA_L1_v1(data_all[, colname], data_all[, true_data_colname], eps)))
  } else if(loss == "TADDA2") {
    df_losses <- list.cbind(lapply(predictions_colnames, function(colname) TADDA_L1_v2(data_all[, colname], data_all[, true_data_colname], eps)))
  } else { stop("'loss' must bei either 'SE', 'TADDA1' or 'TADDA2'")  }
  colnames(df_losses) <- predictions_colnames
  colnames(df_losses)[ncol(df_losses)] <- paste(colnames(df_losses)[ncol(df_losses)], "_s", s, sep = "")
  df_losses <- data.frame("country_name" = data_all$country_name, "month_id" = data_all$month_id, df_losses)
  df_losses
}

compute_losses <- function(data_predictions) {
  loss_SE_ls <- lapply(1:7, function(s) compute_loss_per_s_and_loss_fun(s, loss = "SE", data_all = data_predictions))
  loss_TADDA1_ls <- lapply(1:7, function(s) compute_loss_per_s_and_loss_fun(s, loss = "TADDA1", data_all = data_predictions))
  loss_TADDA2_ls <- lapply(1:7, function(s) compute_loss_per_s_and_loss_fun(s, loss = "TADDA2", data_all = data_predictions))
  
  loss_SE <- list.cbind(loss_SE_ls)[, !duplicated(colnames(list.cbind(loss_SE_ls)))]
  loss_TADDA1 <- list.cbind(loss_TADDA1_ls)[, !duplicated(colnames(list.cbind(loss_TADDA1_ls)))]
  loss_TADDA2 <- list.cbind(loss_TADDA2_ls)[, !duplicated(colnames(list.cbind(loss_TADDA2_ls)))]
  
  colnames(loss_SE)[3:length(colnames(loss_SE))] <- paste("loss_SE_", colnames(loss_SE)[3:length(colnames(loss_SE))], sep = "")
  colnames(loss_TADDA1)[3:length(colnames(loss_TADDA1))] <- paste("loss_TADDA1_", colnames(loss_TADDA1)[3:length(colnames(loss_TADDA1))], sep = "")
  colnames(loss_TADDA2)[3:length(colnames(loss_TADDA2))] <- paste("loss_TADDA2_", colnames(loss_TADDA2)[3:length(colnames(loss_TADDA2))], sep = "")
  
  loss_SE_ordered <- loss_SE %>% select(., "country_name", "month_id", starts_with("loss_SE_BA_SE"), starts_with("loss_SE_BA_TADDA1"), starts_with("loss_SE_BA_TADDA2"), starts_with("loss_SE_No_Change"))
  loss_TADDA1_ordered <- loss_TADDA1 %>% select(., "country_name", "month_id", starts_with("loss_TADDA1_BA_SE"), starts_with("loss_TADDA1_BA_TADDA1"), starts_with("loss_TADDA1_BA_TADDA2"), starts_with("loss_TADDA1_No_Change"))
  loss_TADDA2_ordered <- loss_TADDA2 %>% select(., "country_name", "month_id", starts_with("loss_TADDA2_BA_SE"), starts_with("loss_TADDA2_BA_TADDA1"), starts_with("loss_TADDA2_BA_TADDA2"), starts_with("loss_TADDA2_No_Change"))
  
  merge(loss_SE_ordered, merge(loss_TADDA1_ordered, loss_TADDA2_ordered, by = c("country_name", "month_id")), by = c("country_name", "month_id"))
}

loss <- loss_all
month_id_task <- pred_month_id_task2
length(mloss)

mean_loss <- function(loss, month_id_task) {
  mloss <- colMeans(loss[which(loss$month_id %in% month_id_task), ][-c(1,2)])
  mloss_ordered <- matrix(mloss, nrow = 7, byrow = FALSE)
  mloss_ordered <- mloss_ordered[-1,]
  mloss_ordered <- rbind(mloss_ordered, colMeans(mloss_ordered))
  
  colnames(mloss_ordered) <- c("MSE_BA_SE", "MSE_BA_TADDA1", "MSE_BA_TADDA2", "MSE_No_Change",
                               "MTADDA1_BA_SE", "MTADDA1_BA_TADDA1", "MTADDA1_BA_TADDA2", "MTADDA1_No_Change",
                               "MTADDA2_BA_SE", "MTADDA2_BA_TADDA1", "MTADDA2_BA_TADDA2", "MTADDA2_No_Change")
  rownames(mloss_ordered) <- c(paste("s", 2:7, sep = ""), "colmean")
  
  round(mloss_ordered, 3)
}

