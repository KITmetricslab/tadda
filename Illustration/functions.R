### Some functions used to generate illustrations and simulation examples.
# Johannes Bracher
# johannes.bracher@kit.edu

# absolute error
AE <- function(y_hat, y){
  abs(y_hat - y)
}

# TADDA score with L1 loss, epsilon = 0
TADDA_L1 <- function(y_hat, y){
  abs(y_hat - y) + (sign(y_hat) != sign(y))*abs(y_hat)
}

# TADDA score with L1 loss, epsilon > 0, version 1
TADDA_L1_v1 <- function(y_hat, y, epsilon){
  ae <- abs(y_hat - y)
  penalty1 <- (y_hat > epsilon & y < -epsilon) * abs(y_hat - epsilon)
  penalty2 <- (y_hat < -epsilon & y > epsilon) * abs(y_hat + epsilon)
  ae + penalty1 + penalty2
}

# TADDA score with L1 loss, epsilon > 0, version 2
TADDA_L1_v2 <- function(y_hat, y, epsilon){
  ae <- abs(y_hat - y)
  penalty1 <- ((y_hat < epsilon & y > epsilon) | 
                 (y_hat > epsilon & y >= -epsilon & y <= epsilon)) * abs(y_hat - epsilon)
  penalty2 <- ((y_hat > -epsilon & y < -epsilon) | 
                 (y_hat < -epsilon & y >= -epsilon & y <= epsilon)) * abs(y_hat + epsilon)
  ae + penalty1 + penalty2
}

# squared error
SE <- function(y_hat, y){
  (y_hat - y)^2
}

# TADDA score with L2 loss, epsilon = 0
TADDA_L2 <- function(y_hat, y){
  (y_hat - y)^2 + (sign(y_hat) != sign(y))*(y_hat)^2
}

# TADDA score with L2 loss, epsilon > 0, version 1
TADDA_L2_v1 <- function(y_hat, y, epsilon){
  se <- (y_hat - y)^2
  penalty1 <- (y_hat > epsilon & y < -epsilon) * (y_hat - epsilon)^2
  penalty2 <- (y_hat < -epsilon & y > epsilon) * (y_hat + epsilon)^2
  se + penalty1 + penalty2
}

# TADDA score with L2 loss, epsilon > 0, version 2
TADDA_L2_v2 <- function(y_hat, y, epsilon){
  se <- (y_hat - y)^2
  penalty1 <- ((y_hat < epsilon & y > epsilon) | 
                 (y_hat > epsilon & y >= -epsilon & y <= epsilon)) * (y_hat - epsilon)^2
  penalty2 <- ((y_hat > -epsilon & y < -epsilon) | 
                 (y_hat < -epsilon & y >= -epsilon & y <= epsilon)) * (y_hat + epsilon)^2
  se + penalty1 + penalty2
}

# compute the Bayes act from grid_y for a given score under a given forecast distribution (samples_y);
# repeated for a set of values epsilon
get_bayes_acts <- function(grid_y, grid_epsilon, samples_y, score){
  bayes_acts <- numeric(length(grid_epsilon))
  for(i in seq_along(grid_epsilon)){
    scores_temp <- sapply(grid_y, function(y_hat) mean(score(y_hat, samples_y, grid_epsilon[i])))
    bayes_acts[i] <- grid_y[which.min(scores_temp)]
  }
  bayes_acts
}

# compute the Bayes act from grid_y for a given score under a given forecast distribution (samples_y);
# without epsilon for AE, SE, TADDA_L1 and TADDA_L2
get_bayes_acts_wo_epsilon <- function(grid_y, samples_y, score){
  scores_temp <- sapply(grid_y, function(y_hat) mean(score(y_hat, samples_y)))
  bayes_act <- grid_y[which.min(scores_temp)]
  bayes_act
}

# helper function to add text in a box to a plot:
text_in_box <- function(x, y, txt, col, cex = 1){
  legend(x, y, txt,
         xjust = 0.5,      # 0.5 means center adjusted
         yjust = 0.5,      # 0.5 means center adjusted
         x.intersp = -0.5, # adjust character interspacing as you like to effect box width
         y.intersp = 0.1,  # adjust character interspacing to effect box height
         adj = c(0, 0.5),
         box.col = col,
         text.col = col,
         bg = "white",
         cex = cex)
}
