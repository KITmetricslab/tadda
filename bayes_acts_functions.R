### Functions for Bayes Acts of different scoring functions of AE, SE and TADDA
# Lotta Rüter
# lotta.rueter@kit.edu


## L1 Bayes Acts
BA_AE <- function(distribution_Y) {
  median(distribution_Y)
}

BA_TADDA_L1 <- function(distribution_Y) {
  m <- median(distribution_Y)
  pi <- sum(sign(m) != sign(distribution_Y))/length(distribution_Y)
  ifelse(pi<1/3, quantile(distribution_Y, 0.5 * (1 - sign(m)*pi)), 0)
}

BA_TADDA1_L1 <- function(distribution_Y, epsilon) {
  m <- median(distribution_Y)
  
  if(m < (-epsilon)) {
    F_minus_epsilon <- ecdf(distribution_Y)(-epsilon)
    F_epsilon <- ecdf(distribution_Y)(epsilon)
    critical_value <- F_minus_epsilon / (2 - F_epsilon)
    if (critical_value > 0.5) {
      return(quantile(distribution_Y, 0.5*(2-F_epsilon)))
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
      return(quantile(distribution_Y, 0.5*(1-F_minus_epsilon)))
    }
  }
}

BA_TADDA2_L1 <- function(distribution_Y, epsilon) {
  m <- median(distribution_Y)
  
  # case m < -epsilon
  if(m < (-epsilon)) {
    F_minus_epsilon <- ecdf(distribution_Y)(-epsilon)
    critical_value <- F_minus_epsilon / (2 - F_minus_epsilon)
    if (critical_value > 0.5) {
      return(quantile(distribution_Y, 0.5*(2-F_minus_epsilon)))
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
      return(quantile(distribution_Y, 0.5*(1 - F_epsilon)))
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

# table of all bayes acts
compute_bayes_acts <- function(distribution_Y, epsilon) {
  c(
    BA_AE(distribution_Y),
    BA_TADDA_L1(distribution_Y),
    BA_TADDA1_L1(distribution_Y, epsilon),
    BA_TADDA2_L1(distribution_Y, epsilon),
    BA_SE(distribution_Y),
    BA_TADDA_L2(distribution_Y),
    BA_TADDA1_L2(distribution_Y, epsilon),
    BA_TADDA2_L2(distribution_Y, epsilon)
  )
}

# # scores as loss functions
# loss <- function(y_hat, y_true, score) {
#   loss <- data.frame(matrix(nrow = 0, ncol = ncol(y_hat)))
#   colnames(loss) <- colnames(y_hat)
#   for (i in 1:nrow(y_hat)) {
#     for (j in 1:ncol(y_hat)) {
#       loss[i,j] <- score(y_hat[i,j], y_true[i])
#     }
#   }
#   loss
# }

# scores as loss functions
loss <- function(y_hat, y_true, score) {
  loss <- data.frame(matrix(nrow = 0, ncol = ncol(y_hat)))
  colnames(loss) <- colnames(y_hat)
  for (i in 1:nrow(y_hat)) {
    for (j in 1:ncol(y_hat)) {
      loss[i,j] <- score(y_hat[i,j], y_true[i])
    }
  }
  loss
}

# scores as loss functions with epsilon
loss_epsilon <- function(y_hat, y_true, epsilon, score) {
  loss <- data.frame(matrix(nrow = 0, ncol = ncol(y_hat)))
  colnames(loss) <- colnames(y_hat)
  for (i in 1:nrow(y_hat)) {
    for (j in 1:ncol(y_hat)) {
      loss[i,j] <- score(y_hat[i,j], y_true[i], epsilon)
    }
  }
  loss
}

# loss tables, detailed
loss_tables <- function(y_hat, y_true, epsilon) {
  y_hat <- data.frame(y_hat, "No_Change" = 0) # add no-change baseline
  loss_AE <- loss(y_hat, y_true, AE)
  loss_TADDA_L1 <- loss(y_hat, y_true, TADDA_L1)
  loss_TADDA1_L1 <- loss_epsilon(y_hat, y_true, epsilon, TADDA_L1_v1)
  loss_TADDA2_L1 <- loss_epsilon(y_hat, y_true, epsilon, TADDA_L1_v2)
  
  loss_SE <- loss(y_hat, y_true, SE)
  loss_TADDA_L2 <- loss(y_hat, y_true, TADDA_L2)
  loss_TADDA1_L2 <- loss_epsilon(y_hat, y_true, epsilon, TADDA_L2_v1)
  loss_TADDA2_L2 <- loss_epsilon(y_hat, y_true, epsilon, TADDA_L2_v2)
  
  loss_tables_list <- list(loss_AE, loss_TADDA_L1, loss_TADDA1_L1, loss_TADDA2_L1,
                           loss_SE, loss_TADDA_L2, loss_TADDA1_L2, loss_TADDA2_L2) 
  names(loss_tables_list) <- c("AE", "TADDA_L1", "TADDA1_L1", "TADDA2_L1", "SE", "TADDA_L2", "TADDA1_L2", "TADDA2_L2")
  loss_tables_list
}

# overview loss table
mean_loss_table <- function(y_hat, y_true, epsilon) {
  y_hat <- data.frame(y_hat, "No_Change" = 0) # add no-change baseline
  loss_AE <- loss(y_hat, y_true, AE)
  loss_TADDA_L1 <- loss(y_hat, y_true, TADDA_L1)
  loss_TADDA1_L1 <- loss_epsilon(y_hat, y_true, epsilon, TADDA_L1_v1)
  loss_TADDA2_L1 <- loss_epsilon(y_hat, y_true, epsilon, TADDA_L1_v2)
  
  loss_SE <- loss(y_hat, y_true, SE)
  loss_TADDA_L2 <- loss(y_hat, y_true, TADDA_L2)
  loss_TADDA1_L2 <- loss_epsilon(y_hat, y_true, epsilon, TADDA_L2_v1)
  loss_TADDA2_L2 <- loss_epsilon(y_hat, y_true, epsilon, TADDA_L2_v2)
  
  loss_table <- rbind(colMeans(loss_AE), colMeans(loss_TADDA_L1), colMeans(loss_TADDA1_L1), colMeans(loss_TADDA2_L1),
                      colMeans(loss_SE), colMeans(loss_TADDA_L2), colMeans(loss_TADDA1_L2), colMeans(loss_TADDA2_L2))
  
  row_minima <- apply(loss_table, 1, min)
  minimizer <- c()
  
  for (r in 1:nrow(loss_table)) {
    minimizer[r] <- paste(colnames(y_hat)[which(loss_table[r,]%in%row_minima[r])], collapse=', ')
  }
  
  loss_table <- data.frame(loss_table, "Minimizer" = minimizer)
  rownames(loss_table) <- c("MAE", "MTADDA_L1", "MTADDA1_L1", "MTADDA2_L1", "MSE", "MTADDA_L2", "MTADDA1_L2", "MTADDA2_L2")
  
  loss_table
}
