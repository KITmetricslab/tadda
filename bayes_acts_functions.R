### Functions for Bayes Acts of different scoring functions of AE, SE and TADDA
# Lotta Rüter
# lotta.rueter@kit.edu

# ----
### BAYES ACTS
# ----
## L1 Bayes Acts
BA_AE <- function(distribution_Y) {
  median(distribution_Y, type = 4)
}

BA_TADDA_L1 <- function(distribution_Y) {
  m <- median(distribution_Y, type = 4)
  pi <- sum(sign(m) != sign(distribution_Y))/length(distribution_Y)
  ifelse(pi<1/3, quantile(distribution_Y, 0.5 * (1 - sign(m)*pi), type = 4), 0)
}

BA_TADDA1_L1 <- function(distribution_Y, epsilon) {
  m <- median(distribution_Y, type = 4)
  
  if(m < (-epsilon)) {
    F_minus_epsilon <- ecdf(distribution_Y)(-epsilon)
    F_epsilon <- ecdf(distribution_Y)(epsilon)
    critical_value <- F_minus_epsilon / (2 - F_epsilon)
    if (critical_value > 0.5) {
      return(quantile(distribution_Y, 0.5*(2-F_epsilon), type = 4))
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
      return(quantile(distribution_Y, 0.5*(1-F_minus_epsilon), type = 4))
    }
  }
}

BA_TADDA2_L1 <- function(distribution_Y, epsilon) {
  m <- median(distribution_Y, type = 4)
  
  # case m < -epsilon
  if(m < (-epsilon)) {
    F_minus_epsilon <- ecdf(distribution_Y)(-epsilon)
    critical_value <- F_minus_epsilon / (2 - F_minus_epsilon)
    if (critical_value > 0.5) {
      return(quantile(distribution_Y, 0.5*(2-F_minus_epsilon), type = 4))
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
      return(quantile(distribution_Y, 0.5*(1 - F_epsilon), type = 4))
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
    "BA_TADDA1_L1" = BA_TADDA1_L1(distribution_Y, epsilon),
    "BA_TADDA2_L1" = BA_TADDA2_L1(distribution_Y, epsilon)
  )
}

# predictions for rolling windows
bayes_acts_predictions <- function(variant = "informed", data_c = data_country, s = step_ahead, w_length = window_length, eps = epsilon) {
  if(variant == "uninformed") {
    col_name <- paste("log_change_s", s, sep = "")
    predictions <- lapply(1:(nrow(data_c)-w_length), function(w_begin) {
      log_change_distribution <- data_c[w_begin:(w_begin + w_length - 1), col_name]
      preds <- data.frame("month_id" = data_c$month_id[w_begin + w_length + s - 1], bayes_acts(log_change_distribution, eps))
      colnames(preds) <- c("month_id", paste(colnames(preds)[2:ncol(preds)], "_s", s, sep = ""))
      rownames(preds) <- NULL
      preds
      })
  } else if(variant == "informed") {
    predictions <- lapply(1:(nrow(data_c)-w_length-s), function(w_begin) {
      w_end <- w_begin + w_length - 1
      fatalities_w <- data_c$fatalities[w_begin:w_end]
      log_change_distribution <- log(fatalities_w + 1) - log(tail(fatalities_w, 1) + 1) # use last observed value as reference point for log change distribution
      preds <- data.frame("month_id" = data_c$month_id[w_end + s], bayes_acts(log_change_distribution, eps))
      colnames(preds) <- c("month_id", paste(colnames(preds)[2:ncol(preds)], "_s", s, sep = ""))
      rownames(preds) <- NULL
      preds
    }) 
  } else if(variant == "combined") { # does not work very well -> discard it
      predictions <- lapply(1:(nrow(data_c)-w_length-s), function(w_begin) {
        combinations <- combn(data_c$fatalities[w_begin:(w_begin + w_length)], 2)# , data_c$fatalities[w_begin:(w_begin + w_length)])
        log_change_distribution <- log(combinations[1,] + 1) - log(combinations[2,] + 1)
        preds <- data.frame("month_id" = data_c$month_id[w_begin + w_length + s], bayes_acts(log_change_distribution, eps))
        colnames(preds) <- c("month_id", paste(colnames(preds)[2:ncol(preds)], "_s", s, sep = ""))
        rownames(preds) <- NULL
        preds
      })
  } else { stop("variant must bei either 'informed', 'uninformed' or 'combined'") }
  
  bind_rows(predictions)
}

compute_losses <- function(s, loss, data_all, eps = epsilon) {
  true_data_colname <- paste("log_change_s", s, sep = "")
  predictions_colnames <- c(paste(BA_names, s, sep = "_s"), "No_Change")
  if(loss == "SE") {
    df_losses <- SE(data_all[, predictions_colnames], data_all[, true_data_colname])
  } else if(loss == "TADDA1_L1") {
    df_losses <- list.cbind(lapply(predictions_colnames, function(colname) TADDA_L1_v1(data_all[, colname], data_all[, true_data_colname], eps)))
  } else if(loss == "TADDA2_L1") {
    df_losses <- list.cbind(lapply(predictions_colnames, function(colname) TADDA_L1_v2(data_all[, colname], data_all[, true_data_colname], eps)))
  } else { stop("'loss' must bei either 'SE', 'TADDA1_L1' or 'TADDA2_L1'")  }
  colnames(df_losses) <- predictions_colnames
  colnames(df_losses)[ncol(df_losses)] <- paste(colnames(df_losses)[ncol(df_losses)], "_s", s, sep = "")
  df_losses <- cbind("month_id" = data_all$month_id, df_losses)
  df_losses
}

loss_df <- function(data_predictions, loss_fun) {
  loss <- lapply(1:7, function(s) compute_losses(s, loss = loss_fun, data_all = data_predictions))
  list.cbind(loss)[, !duplicated(colnames(list.cbind(loss)))]
}

mean_loss <- function(data_predictions, month_id_task) {
  SE_loss_df <- loss_df(data_predictions, "SE")
  TADDA1_L1_loss_df <- loss_df(data_predictions, "TADDA1_L1")
  TADDA2_L1_loss_df <- loss_df(data_predictions, "TADDA2_L1")
  
  mean_loss_SE <- colMeans(SE_loss_df[which(SE_loss_df$month_id %in% month_id_task), ])
  mean_loss_TADDA1_L1 <- colMeans(TADDA1_L1_loss_df[which(TADDA1_L1_loss_df[,1] %in% month_id_task), ])
  mean_loss_TADDA2_L1 <- colMeans(TADDA2_L1_loss_df[which(TADDA2_L1_loss_df[,1] %in% month_id_task), ])
  
  rbind(mean_loss_SE, mean_loss_TADDA1_L1, mean_loss_TADDA2_L1)
}
