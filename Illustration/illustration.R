### Some illustrations and simulation examples on the TADDA scores
# Johannes Bracher
# johannes.bracher@kit.edu

# setwd("/home/johannes/Documents/Ideas/armed_conflicts/tadda/Illustration")

# library for skew normal
library(sn)
library(plotrix)

current_path <- rstudioapi::getActiveDocumentContext()$path # get path of this file
setwd(dirname(current_path))

source("functions.R")

# parameters of skew normal:
xi <- -0.15
omega <- 0.4
alpha <- 8

# tolerance value:
epsilon <- 0.048

# generate samples:
set.seed(123)
samples_y <- rsn(100000, xi = xi, omega = omega, alpha = alpha)


# compute mean and median:
(mu <- mean(samples_y))
(med <- median(samples_y))
grid_y <- seq(from = 0, to = 0.2, by = 0.001)
(ba_tadda1 <- get_bayes_acts(grid_y = grid_y, grid_epsilon = epsilon, samples_y = samples_y, score = TADDA_L1_v1))

#####################
### plot density function:

# grid:
grid_y_density <- seq(from = -0.4, to = 1.2, by = 0.001)

# evaluate density:
f <- dsn(x = grid_y_density, xi = xi, omega = omega, alpha = alpha)


pdf("figures/illustration.pdf", width = 6, height = 3)
par(mar = c(4.2, 4.2, 0.5, 0.5), las = 1)
plot(grid_y_density, f, type = "l", ylab = "f(y)", xlab = "y", ylim = c(0, 2))
# rect(-epsilon, 0, epsilon, c(-1, 2.5), col = "grey97", border = NA)
# abline(v = epsilon, lty = 3, col = "grey")
# abline(v = -epsilon, lty = 3, col = "grey")

box()
lines(grid_y_density, f)

abline(v = mu, lty = 2, col = "black")
abline(v = med, lty = 3, col = "black")
abline(v = ba_tadda1, lty = 4, col = "red")

legend("topright", legend = c(paste0("Mean: ", round(mu, 3)), # hard-coded, modify as needed
                              paste0("Median: ", round(med, 3)),
                              expression(BA~under~TADDA[0.048]: 0.06)),
                              # expression(BA~under~TADDA[0]: 0.019)),
       lty = c(2, 3, 4, 5), bty = "n", col = c("black", "black", "red"), cex = 0.9)
dev.off()



#####################
### plot expected scores as function of y_hat:

epsilon <- 0.048

grid_y_hat <- seq(from = -0.1, to = 0.2, by = 0.001)
average_scores_ae <-
  average_scores_se <-
  average_scores_tadda_l1_0 <-
  average_scores_tadda_l1_v1_epsilon <-
  average_scores_tadda_l1_v2_epsilon <-
  average_scores_tadda_l2_0 <-
  average_scores_tadda_l2_v1_epsilon <-
  average_scores_tadda_l2_v2_epsilon <-
  numeric(length(grid_y_hat))

for(i in seq_along(grid_y_hat)){
  y_hat_temp <- grid_y_hat[i]
  average_scores_ae[i] <- mean(abs(y_hat_temp - samples_y))
  average_scores_se[i] <- mean((y_hat_temp - samples_y)^2)
  
  average_scores_tadda_l1_0[i] <- mean(TADDA_L1(y_hat_temp, samples_y))
  average_scores_tadda_l1_v1_epsilon[i] <- mean(TADDA_L1_v1(y_hat_temp, samples_y, epsilon = epsilon))
  average_scores_tadda_l1_v2_epsilon[i] <- mean(TADDA_L1_v2(y_hat_temp, samples_y, epsilon = epsilon))
  
  average_scores_tadda_l2_0[i] <- mean(TADDA_L2(y_hat_temp, samples_y))
  average_scores_tadda_l2_v1_epsilon[i] <- mean(TADDA_L2_v1(y_hat_temp, samples_y, epsilon = epsilon))
  average_scores_tadda_l2_v2_epsilon[i] <- mean(TADDA_L2_v2(y_hat_temp, samples_y, epsilon = epsilon))
}

pdf("figures/expected_scores.pdf", width = 7, height = 2.6)
layout(matrix(1:3, nrow = 1), widths = c(2, 2, 1))
par(las = 1, mar = c(4.2, 4.2, 2.5, 2))
yl <- c(0.1, 0.4)
cex.pt <- 2

# L1
plot(grid_y_hat, average_scores_tadda_l1_0, type = "l", col = "black", xlab = expression(hat(y)),
     ylab = "expected score", xlim = c(-0.1, 0.2), ylim = yl)
mtext("(a) L1", side = 3, line = 0.5)


rect(-epsilon, 0, epsilon, 1.2*yl[2], col = "grey97", border = NA)

abline(v = med, col = "darkgrey")
text_in_box(med, 0.36, "m", col = "darkgrey")

# abline(v = mu,  col = "darkgrey")
# text_in_box(mu, 0.36, expression(mu), col = "darkgrey")

abline(v = 0, col = "darkgrey", lty = "dotted")
text_in_box(0, 0.9*yl[2], 0, col = "darkgrey")

abline(v = epsilon,  col = "darkgrey", lty = "dotted")
text_in_box(epsilon, 0.9*yl[2], expression(epsilon), col = "darkgrey")

abline(v = -epsilon,  col = "darkgrey", lty = "dotted")
text_in_box(-epsilon, 0.9*yl[2], expression(-epsilon), col = "darkgrey")

lines(grid_y_hat, average_scores_tadda_l1_0, type = "l", col = "black", lty = "solid")
lines(grid_y_hat, average_scores_tadda_l1_v1_epsilon, type = "l", col = "red", lty = "dashed")
lines(grid_y_hat, average_scores_tadda_l1_v2_epsilon, type = "l", col = "blue", lty = "twodash")
lines(grid_y_hat, average_scores_ae, type = "l", col = "darkgreen", lty = "dotted")

points(grid_y_hat[which.min(average_scores_tadda_l1_0)], min(average_scores_tadda_l1_0),
       col = "black", pch = 18, cex = cex.pt)
points(grid_y_hat[which.min(average_scores_tadda_l1_v1_epsilon)], min(average_scores_tadda_l1_v1_epsilon),
       col = "red", pch = 18, cex = cex.pt)
points(grid_y_hat[which.min(average_scores_tadda_l1_v2_epsilon)], min(average_scores_tadda_l1_v2_epsilon),
       col = "blue", pch = 18, cex = cex.pt)
points(grid_y_hat[which.min(average_scores_ae)], min(average_scores_ae),
       col = "darkgreen", pch = 18, cex = cex.pt)

box()

# L2
yl2 <- c(0.04, 0.16)
plot(grid_y_hat, average_scores_tadda_l2_0, type = "l", col = "black", xlab = expression(hat(y)),
     ylab = "expected score", xlim = c(-0.1, 0.2), ylim = yl2)
mtext("(b) L2", side = 3, line = 0.5)

rect(-epsilon, 0, epsilon, 1.2*yl2[2], col = "grey97", border = NA)

# abline(v = med, col = "darkgrey")
# text_in_box(med, 0.8*yl2[2], "m", col = "darkgrey")

abline(v = mu,  col = "darkgrey")
text_in_box(mu, 0.9*yl2[2], expression(mu), col = "darkgrey")

abline(v = 0, col = "darkgrey", lty = "dotted")
text_in_box(0, 0.9*yl2[2], 0, col = "darkgrey")

abline(v = epsilon,  col = "darkgrey", lty = "dotted")
text_in_box(epsilon, 0.9*yl2[2], expression(epsilon), col = "darkgrey")

abline(v = -epsilon,  col = "darkgrey", lty = "dotted")
text_in_box(-epsilon, 0.9*yl2[2], expression(-epsilon), col = "darkgrey")

lines(grid_y_hat, average_scores_tadda_l2_0, type = "l", col = "black", lty = "solid")
lines(grid_y_hat, average_scores_tadda_l2_v1_epsilon, type = "l", col = "red", lty = "dashed")
lines(grid_y_hat, average_scores_tadda_l2_v2_epsilon, type = "l", col = "blue", lty = "twodash")
lines(grid_y_hat, average_scores_se, type = "l", col = "purple", lty = "dotted")

points(grid_y_hat[which.min(average_scores_tadda_l2_0)], min(average_scores_tadda_l2_0),
       col = "black", pch = 18, cex = cex.pt)
points(grid_y_hat[which.min(average_scores_tadda_l2_v1_epsilon)], min(average_scores_tadda_l2_v1_epsilon),
       col = "red", pch = 18, cex = cex.pt)
points(grid_y_hat[which.min(average_scores_tadda_l2_v2_epsilon)], min(average_scores_tadda_l2_v2_epsilon),
       col = "blue", pch = 18, cex = cex.pt)
points(grid_y_hat[which.min(average_scores_se)], min(average_scores_se),
       col = "purple", pch = 18, cex = cex.pt)

box()

# Legend
par(mar = c(0, 0, 0, 0))
plot(NULL, xlim = 0:1, ylim = 0:1, xlab = "", ylab = "", axes = FALSE)
legend("center", legend = c(expression(AE(hat(y), y)),
                            expression(SE(hat(y), y)),
                            expression({TADDA[0]}(hat(y), y)),
                            expression({TADDA1[epsilon]}(hat(y), y)),
                            expression({TADDA2[epsilon]}(hat(y), y))),
       bty = "n", cex = 1, lty = c("dotted", "dotdash", "solid", "dashed", "twodash"), 
       col = c("darkgreen", "purple", "black", "red", "blue"),
       pch = 18, pt.cex = 2,
       y.intersp = 1.5)
dev.off()


#####################
### Table with expected scores of Bayes acts:
expected_scores <- matrix(ncol = 9, nrow = 9)
colnames(expected_scores) <- c("value", "AE", "SE", 
                               "TADDA_L1_0", "TADDA_L1_v1", "TADDA_L1_v2",
                               "TADDA_L2_0", "TADDA_L2_v1", "TADDA_L2_v2")
rownames(expected_scores) <- c("median", "mean",
                               "BA TADDA_L1_0", "BA TADDA_L1_v1", "BA TADDA_L1_v2",
                               "BA TADDA_L2_0", "BA TADDA_L2_v1", "BA TADDA_L2_v2",
                               "zero")

ba_tadda_l1_v1_epsilon <- grid_y_hat[which.min(average_scores_tadda_l1_v1_epsilon)]
ba_tadda_l1_v2_epsilon <- grid_y_hat[which.min(average_scores_tadda_l1_v2_epsilon)]
ba_tadda_l2_v1_epsilon <- grid_y_hat[which.min(average_scores_tadda_l2_v1_epsilon)]
ba_tadda_l2_v2_epsilon <- grid_y_hat[which.min(average_scores_tadda_l2_v2_epsilon)]


expected_scores[, "value"] <- c(med, mu, 
                                med_modified, ba_tadda_l1_v1_epsilon, ba_tadda_l1_v2_epsilon,
                                mu_modified, ba_tadda_l2_v1_epsilon, ba_tadda_l2_v2_epsilon,
                                0)

expected_scores[, "AE"] <- sapply(expected_scores[, "value"], 
                                  function(y_hat) mean(abs(y_hat - samples_y)))
expected_scores[, "SE"] <- sapply(expected_scores[, "value"], 
                                  function(y_hat) mean((y_hat - samples_y)^2))

expected_scores[, "TADDA_L1_0"] <- sapply(expected_scores[, "value"], 
                                          function(x) mean(TADDA_L1(x, y = samples_y)))
expected_scores[, "TADDA_L1_v1"] <- sapply(expected_scores[, "value"],
                                           function(x) mean(TADDA_L1_v1(x, y = samples_y, epsilon = epsilon)))
expected_scores[, "TADDA_L1_v2"] <- sapply(expected_scores[, "value"],
                                           function(x) mean(TADDA_L1_v2(x, y = samples_y, epsilon = epsilon)))

expected_scores[, "TADDA_L2_0"] <- sapply(expected_scores[, "value"], 
                                          function(x) mean(TADDA_L2(x, y = samples_y)))
expected_scores[, "TADDA_L2_v1"] <- sapply(expected_scores[, "value"],
                                           function(x) mean(TADDA_L2_v1(x, y = samples_y, epsilon = epsilon)))
expected_scores[, "TADDA_L2_v2"] <- sapply(expected_scores[, "value"],
                                           function(x) mean(TADDA_L2_v2(x, y = samples_y, epsilon = epsilon)))

expected_scores_small <- rbind(expected_scores[c("median", "mean"), ],
                               BA = NA,
                               zero = expected_scores["zero", ])
expected_scores_small <- expected_scores_small[, c("AE", "SE", 
                                                   "TADDA_L1_v1")]
expected_scores_small["BA", ] <- diag(expected_scores[-nrow(expected_scores), -1])

library(xtable)
# scores:
xtable(expected_scores_small, 3, digits = 3)

# Bayes acts:
xtable(matrix(expected_scores[1:8, "value"], nrow = 1), digits = 3)

#########################
### plotting Bayes act as a function of epsilon:

# grif of y-values to check for Bayes act:
grid_y <- seq(from = 0.01, to = 0.3, by = 0.002)

# generate a shorter set of samples to speed up computation:
samples_y_short <- rsn(30000, xi = xi, omega = omega, alpha = alpha)

# values of epsilon for which to compute Bayes acts:
grid_epsilon <- (seq(from = sqrt(0.01), to = sqrt(0.5), length.out = 100))^2

# compute Bayes acts under different scores:
bayes_acts_L1_v1 <- get_bayes_acts(grid_y = grid_y, grid_epsilon = grid_epsilon,
                                   samples_y = samples_y_short, score = TADDA_L1_v1)
bayes_acts_L1_v2 <- get_bayes_acts(grid_y = grid_y, grid_epsilon = grid_epsilon,
                                   samples_y = samples_y_short, score = TADDA_L1_v2)

bayes_acts_L2_v1 <- get_bayes_acts(grid_y = grid_y, grid_epsilon = grid_epsilon,
                                   samples_y = samples_y_short, score = TADDA_L2_v1)
bayes_acts_L2_v2 <- get_bayes_acts(grid_y = grid_y, grid_epsilon = grid_epsilon,
                                   samples_y = samples_y_short, score = TADDA_L2_v2)


# Plot:
pdf("figures/ba_epsilon.pdf", width = 7, height = 2.6)

# tructure plot region:
layout(matrix(1:3, nrow = 1), widths = c(2, 2, 1))
par(las = 1, mar = c(4.2, 4.2, 2.5, 2))
# manually preparing square-root axes...
x_labels <- c(0.01, 0.1, 0.25, 0.5, 1, 2, 4)
y_labels <- c(0.01, 0.1, 0.25, 0.5, 1)

# L1:
plot(sqrt(grid_epsilon), sqrt(bayes_acts_L1_v1), type = "l", col = "red", lty = "dashed",
     xlab = expression(epsilon), ylab = "Bayes act", ylim = sqrt(c(0.01, 0.3)), axes = FALSE)
abline(0:1, col = "darkgrey")
mtext("(a) L1", side = 3, line = 0.5)
axis(1, at = sqrt(x_labels), labels = x_labels)
axis(2, at = sqrt(x_labels), labels = x_labels)

lines(sqrt(grid_epsilon), sqrt(bayes_acts_L1_v2), col = "blue", lty = "twodash")
lines(sqrt(grid_epsilon), sqrt(bayes_acts_L1_v1), type = "l", col = "red", lty = "dashed")
abline(h = sqrt(med), col = "darkgrey")
abline(v = sqrt(med), col = "darkgrey")

text_in_box(0.25, sqrt(med), "m", col = "darkgrey")
text_in_box(sqrt(med), 0.25, "m", col = "darkgrey")

text_in_box(sqrt(0.28), sqrt(0.28), expression(epsilon), col = "darkgrey")

box()

# L2:
plot(sqrt(grid_epsilon), sqrt(bayes_acts_L2_v1), type = "l", col = "red", lty = "dashed",
     xlab = expression(epsilon), ylab = "Bayes act", ylim = sqrt(c(0.01, 0.3)), axes = FALSE)
abline(0:1, col = "darkgrey")
mtext("(b) L2", side = 3, line = 0.5)
axis(1, at = sqrt(x_labels), labels = x_labels)
axis(2, at = sqrt(x_labels), labels = x_labels)

lines(sqrt(grid_epsilon), sqrt(bayes_acts_L2_v2), col = "blue", lty = "twodash")
lines(sqrt(grid_epsilon), sqrt(bayes_acts_L2_v1), type = "l", col = "red", lty = "dashed")
abline(h = sqrt(mu), col = "darkgrey")
abline(v = sqrt(mu), col = "darkgrey")
text_in_box(0.25, sqrt(mu), expression(mu), col = "darkgrey")
text_in_box(sqrt(mu), 0.25, expression(mu), col = "darkgrey")
text_in_box(sqrt(0.28), sqrt(0.28), expression(epsilon), col = "darkgrey")
box()

# Legend:
par(mar = c(0, 0, 0, 0))
plot(NULL, xlim = 0:1, ylim = 0:1, xlab = "", ylab = "", axes = FALSE)
legend("center", legend = c(expression({TADDA1[epsilon]}(hat(y), y)),
                            expression({TADDA2[epsilon]}(hat(y), y))),
       bty = "n", cex = 1, lty = c("dashed", "twodash"), 
       col = c("red", "blue"),
       y.intersp = 1.5)
dev.off()



#################
# Illustration of TADDA scores:

epsilon <- 0.048

### L1
pdf("figures/curves_scores_L1.pdf", width = 8, height = 4)
layout(matrix(c(1, 2, 5,
                3, 4, 5), byrow = TRUE, ncol = 3), widths = c(2, 2, 1))
par(mar = c(4.2, 4, 3, 0.5), las = 1)
xlim <- c(-0.15, 0.15)
ylim <- c(0, 0.4)
ylab <- "score" # expression({TADDA^(1)}(hat(y), y))
shift <- 0.01

# with fixed y, as a function of y_hat

# y in [-epsilon, epsilon]
grid_y <- seq(from = -0.15, to = 0.25, by = 0.005)
y_true <- 0.03
grid_TADDA_a <- TADDA_L1(grid_y, y_true)
grid_TADDA_v1_a <- TADDA_L1_v1(grid_y, y_true, epsilon = epsilon)
grid_TADDA_v2_a <- TADDA_L1_v2(grid_y, y_true, epsilon = epsilon)
grid_ae_a <- abs(grid_y - y_true)

plot(grid_y, grid_TADDA_v1_a, xlab = expression(hat(y)),
     ylab = ylab, ylim = ylim, xlim = xlim, col = "white")
mtext(expression(
  TADDA~with~y%in%paste("(", -epsilon, ",", epsilon, ")"~", varying"~hat(y))
), 3, cex = 0.9, line = 0.3)

rect(-epsilon, -1, epsilon, 11, col = "grey97", border = NA)
abline(h = 0:10, col = "grey95")
abline(v = -6:8, col = "grey95")
box()
# lines(grid_y, grid_TADDA_a + shift, col = "red", lty = 5)
lines(grid_y, grid_TADDA_v1_a, col = "red", lty = 4)
lines(grid_y, grid_ae_a - shift, col = "black", lty = 3)
# lines(grid_y, grid_TADDA_v2_a + shift, col = "blue", lty = "twodash")

abline(v = 0, col = "darkgrey", lty = "dotted")
abline(v = y_true, lty = "solid", col = "darkgrey")
abline(v = c(epsilon, -epsilon), lty = "dotted", col = "darkgrey")

text_in_box(0, 0.9*ylim[2], expression(0), col = "darkgrey")
text_in_box(epsilon, 0.9*ylim[2], expression(epsilon), col = "darkgrey")
text_in_box(-epsilon,  0.9*ylim[2], expression(-epsilon), col = "darkgrey")
text_in_box(y_true,  0.5*ylim[2], expression(y), col = "darkgrey")


# y > epsilon
y_true <- 0.1
grid_TADDA_b <- TADDA_L1(grid_y, y_true)
grid_TADDA_v1_b <- TADDA_L1_v1(grid_y, y_true, epsilon = epsilon)
grid_TADDA_v2_b <- TADDA_L1_v2(grid_y, y_true, epsilon = epsilon)
grid_ae_b <- abs(grid_y - y_true)


plot(grid_y, grid_TADDA_v1_b + shift, type = "l", xlab = expression(hat(y)),
     ylab = "", ylim = ylim, xlim = xlim, col = "white")
mtext(expression(TADDA~with~y > epsilon~", varying"~hat(y)), 3, cex = 0.9, line = 0.3)
rect(-epsilon, -1, epsilon, 11, col = "grey97", border = NA)
abline(h = 0:10, col = "grey95")
abline(v = -6:8, col = "grey95")
box()
# lines(grid_y, grid_TADDA_b + shift, col = "red", lty = 5)
lines(grid_y, grid_TADDA_v1_b, col = "red", lty = 4)
lines(grid_y, grid_ae_b - shift, col = "black", lty = 3)

# lines(grid_y, grid_TADDA_v2_b + shift, col = "blue", lty = "twodash")

abline(v = 0, col = "darkgrey", lty = "dotted")
abline(v = y_true, lty = "solid", col = "darkgrey")
abline(v = c(epsilon, -epsilon), lty = "dotted", col = "darkgrey")

text_in_box(0, 0.9*ylim[2], expression(0), col = "darkgrey")
text_in_box(epsilon, 0.9*ylim[2], expression(epsilon), col = "darkgrey")
text_in_box(-epsilon,  0.9*ylim[2], expression(-epsilon), col = "darkgrey")
text_in_box(y_true,  0.5*ylim[2], expression(y), col = "darkgrey")


# with fixed y_hat, as a function of y_hat

# y_hat in [-epsilon, epsilon]
y_hat <- 0.03
grid_TADDA_c <- TADDA_L1(y_hat, grid_y)
grid_TADDA_v1_c <- TADDA_L1_v1(y_hat, grid_y, epsilon = epsilon)
# grid_TADDA_v2_c <- TADDA_L1_v2(y_hat, grid_y, epsilon = epsilon)
grid_ae_c <- abs(y_hat - grid_y)

plot(grid_y, grid_TADDA_v1_c + shift, type = "l", xlab = expression(y),
     ylab = ylab, ylim = ylim, xlim = xlim, col = "white")
mtext(
  expression(TADDA~with~hat(y)%in%paste("(", -epsilon, ",", epsilon, ")"~", varying"~y)),
  3, cex = 0.9, line = 0.3)

rect(-epsilon, -1, epsilon, 11, col = "grey97", border = NA)
abline(h = 0:10, col = "grey95")
abline(v = -6:8, col = "grey95")
box()
# lines(grid_y, grid_TADDA_c + shift, col = "red", lty = 5)
lines(grid_y, grid_TADDA_v1_c, col = "red", lty = 4)
lines(grid_y, grid_ae_c - shift, col = "black", lty = 3)

# lines(grid_y, grid_TADDA_v2_c - shift, col = "blue", lty = "twodash")

abline(v = 0, col = "darkgrey", lty = "dotted")
abline(v = y_hat, lty = "solid", col = "darkgrey")
abline(v = c(epsilon, -epsilon), lty = "dotted", col = "darkgrey")

text_in_box(0, 0.9*ylim[2], expression(0), col = "darkgrey")
text_in_box(epsilon, 0.9*ylim[2], expression(epsilon), col = "darkgrey")
text_in_box(-epsilon,  0.9*ylim[2], expression(-epsilon), col = "darkgrey")
text_in_box(y_hat,  0.5*ylim[2], expression(hat(y)), col = "darkgrey")


# y > epsilon
y_hat <- 0.1
grid_TADDA_d <- TADDA_L1(y_hat, grid_y)
grid_TADDA_v1_d <- TADDA_L1_v1(y_hat, grid_y, epsilon = epsilon)
grid_TADDA_v2_d <- TADDA_L1_v2(y_hat, grid_y, epsilon = epsilon)
grid_ae_d <- abs(y_hat - grid_y)


plot(grid_y, grid_TADDA_v1_d + shift, type = "l", xlab = expression(y),
     ylab = "", ylim = ylim, xlim = xlim, col = "white")
mtext(expression(TADDA~with~hat(y) > epsilon~", varying"~y), 3, cex = 0.9, line = 0.3)
rect(-epsilon, -1, epsilon, 11, col = "grey97", border = NA)
abline(h = 0:10, col = "grey95")
abline(v = -6:8, col = "grey95")
box()
# lines(grid_y, grid_TADDA_d + shift, col = "red", lty = 5)
lines(grid_y, grid_TADDA_v1_d, col = "red", lty = 4)
lines(grid_y, grid_ae_d - shift, col = "black", lty = 3)

abline(v = 0, col = "darkgrey", lty = "dotted")
abline(v = y_true, lty = "solid", col = "darkgrey")
abline(v = c(epsilon, -epsilon), lty = "dotted", col = "darkgrey")

text_in_box(0, 0.9*ylim[2], expression(0), col = "darkgrey")
text_in_box(epsilon, 0.9*ylim[2], expression(epsilon), col = "darkgrey")
text_in_box(-epsilon,  0.9*ylim[2], expression(-epsilon), col = "darkgrey")
text_in_box(y_true,  0.5*ylim[2], expression(hat(y)), col = "darkgrey")

par(mar = c(0, 0, 0, 0))
plot(NULL, xlim = 0:1, ylim = 0:1, xlab = "", ylab = "", axes = FALSE)
legend("center", legend = c(# expression({TADDA[0]}(hat(y), y)),
                            expression({TADDA[epsilon]}(hat(y), y)),
                            expression(AE(hat(y), y))),
       bty = "n", cex = 1.2, lty = 5:3, col = c("red", "black"),
       y.intersp = 1.5)
dev.off()



### L2

pdf("figures/curves_scores_L2.pdf", width = 8, height = 4)
layout(matrix(c(1, 2, 5,
                3, 4, 5), byrow = TRUE, ncol = 3), widths = c(2, 2, 1))
par(mar = c(4.2, 4, 3, 0.5), las = 1)
xlim <- c(-0.1, 0.2)
ylim <- c(0, 0.1)
ylab <- "score" # expression({TADDA^(2)}(hat(y), y))
shift <- 0

# with fixed y, as a function of y_hat

# y in [-epsilon, epsilon]
y_true <- 0.03
grid_TADDA_L2_a <- TADDA_L2(grid_y, y_true)
grid_TADDA_L2_v1_a <- TADDA_L2_v1(grid_y, y_true, epsilon = epsilon)
grid_TADDA_L2_v2_a <- TADDA_L2_v2(grid_y, y_true, epsilon = epsilon)

plot(grid_y, grid_TADDA_L2_a, type = "l", xlab = expression(hat(y)),
     ylab = ylab, ylim = ylim, xlim = xlim)
mtext(expression(
  TADDA^(2)~with~y%in%paste("(", -epsilon, ",", epsilon, ")"~", varying"~hat(y))
), 3, cex = 0.9, line = 0.3)

rect(-epsilon, -1, epsilon, 11, col = "grey97", border = NA)
abline(h = 0:10, col = "grey95")
abline(v = -6:8, col = "grey95")
box()
lines(grid_y, grid_TADDA_L2_a)
lines(grid_y, grid_TADDA_L2_v1_a - shift, col = "red", lty = "dashed")
lines(grid_y, grid_TADDA_L2_v2_a + shift, col = "blue", lty = "twodash")

abline(v = 0, col = "darkgrey", lty = "dotted")
abline(v = y_true, lty = "solid", col = "darkgrey")
abline(v = c(epsilon, -epsilon), lty = "dotted", col = "darkgrey")

text_in_box(0, 0.9*ylim[2], expression(0), col = "darkgrey")
text_in_box(epsilon, 0.9*ylim[2], expression(epsilon), col = "darkgrey")
text_in_box(-epsilon,  0.9*ylim[2], expression(-epsilon), col = "darkgrey")
text_in_box(y_true,  0.5*ylim[2], expression(y), col = "darkgrey")


# y > epsilon
y_true <- 0.1
grid_TADDA_L2_b <- TADDA_L2(grid_y, y_true)
grid_TADDA_L2_v1_b <- TADDA_L2_v1(grid_y, y_true, epsilon = epsilon)
grid_TADDA_L2_v2_b <- TADDA_L2_v2(grid_y, y_true, epsilon = epsilon)

plot(grid_y, grid_TADDA_L2_b, type = "l", xlab = expression(hat(y)),
     ylab = "", ylim = ylim, xlim = xlim)
mtext(expression(TADDA^(2)~with~y > epsilon~", varying"~hat(y)), 3, cex = 0.9, line = 0.3)
rect(-epsilon, -1, epsilon, 11, col = "grey97", border = NA)
abline(h = 0:10, col = "grey95")
abline(v = -6:8, col = "grey95")
box()
lines(grid_y, grid_TADDA_L2_b)
lines(grid_y, grid_TADDA_L2_v1_b - shift, col = "red", lty = "dashed")
lines(grid_y, grid_TADDA_L2_v2_b + shift, col = "blue", lty = "twodash")

abline(v = 0, col = "darkgrey", lty = "dotted")
abline(v = y_true, lty = "solid", col = "darkgrey")
abline(v = c(epsilon, -epsilon), lty = "dotted", col = "darkgrey")

text_in_box(0, 0.9*ylim[2], expression(0), col = "darkgrey")
text_in_box(epsilon, 0.9*ylim[2], expression(epsilon), col = "darkgrey")
text_in_box(-epsilon,  0.9*ylim[2], expression(-epsilon), col = "darkgrey")
text_in_box(y_true,  0.5*ylim[2], expression(y), col = "darkgrey")



# with fixed y_hat, as a function of y_hat

# y_hat in [-epsilon, epsilon]
y_hat <- 0.03
grid_TADDA_L2_c <- TADDA_L2(y_hat, grid_y)
grid_TADDA_L2_v1_c <- TADDA_L2_v1(y_hat, grid_y, epsilon = epsilon)
grid_TADDA_L2_v2_c <- TADDA_L2_v2(y_hat, grid_y, epsilon = epsilon)

plot(grid_y, grid_TADDA_L2_c, type = "l", xlab = expression(y),
     ylab = ylab, ylim = ylim, xlim = xlim)
mtext(
  expression(TADDA^(2)~with~hat(y)%in%paste("(", -epsilon, ",", epsilon, ")"~", varying"~y)),
  3, cex = 0.9, line = 0.3)

rect(-epsilon, -1, epsilon, 11, col = "grey97", border = NA)
abline(h = 0:10, col = "grey95")
abline(v = -6:8, col = "grey95")
box()
lines(grid_y, grid_TADDA_L2_c)
lines(grid_y, grid_TADDA_L2_v1_c - shift, col = "red", lty = "dashed")
lines(grid_y, grid_TADDA_L2_v2_c + shift, col = "blue", lty = "twodash")

abline(v = 0, col = "darkgrey", lty = "dotted")
abline(v = y_hat, lty = "solid", col = "darkgrey")
abline(v = c(epsilon, -epsilon), lty = "dotted", col = "darkgrey")

text_in_box(0, 0.9*ylim[2], expression(0), col = "darkgrey")
text_in_box(epsilon, 0.9*ylim[2], expression(epsilon), col = "darkgrey")
text_in_box(-epsilon,  0.9*ylim[2], expression(-epsilon), col = "darkgrey")
text_in_box(y_hat,  0.5*ylim[2], expression(hat(y)), col = "darkgrey")


# y > epsilon
y_hat <- 0.1
grid_TADDA_L2_d <- TADDA_L2(y_hat, grid_y)
grid_TADDA_L2_v1_d <- TADDA_L2_v1(y_hat, grid_y, epsilon = epsilon)
grid_TADDA_L2_v2_d <- TADDA_L2_v2(y_hat, grid_y, epsilon = epsilon)

plot(grid_y, grid_TADDA_L2_d, type = "l", xlab = expression(y),
     ylab = "", ylim = ylim, xlim = xlim)
mtext(expression(TADDA^(2)~with~hat(y) > epsilon~", varying"~y), 3, cex = 0.9, line = 0.3)
rect(-epsilon, -1, epsilon, 11, col = "grey97", border = NA)
abline(h = 0:10, col = "grey95")
abline(v = -6:8, col = "grey95")
box()
lines(grid_y, grid_TADDA_L2_d)
lines(grid_y, grid_TADDA_L2_v1_d - shift, col = "red", lty = "dashed")
lines(grid_y, grid_TADDA_L2_v2_d + shift, col = "blue", lty = "twodash")

abline(v = 0, col = "darkgrey", lty = "dotted")
abline(v = y_true, lty = "solid", col = "darkgrey")
abline(v = c(epsilon, -epsilon), lty = "dotted", col = "darkgrey")

text_in_box(0, 0.9*ylim[2], expression(0), col = "darkgrey")
text_in_box(epsilon, 0.9*ylim[2], expression(epsilon), col = "darkgrey")
text_in_box(-epsilon,  0.9*ylim[2], expression(-epsilon), col = "darkgrey")
text_in_box(y_true,  0.5*ylim[2], expression(hat(y)), col = "darkgrey")

par(mar = c(0, 0, 0, 0))
plot(NULL, xlim = 0:1, ylim = 0:1, xlab = "", ylab = "", axes = FALSE)
legend("center", legend = c(expression({TADDA[0]^(2)}(hat(y), y)),
                            expression({TADDA1[epsilon]^(2)}(hat(y), y)),
                            expression({TADDA2[epsilon]^"(2)"}(hat(y), y))),
       bty = "n", cex = 1.2, lty = c("solid", "dashed", "twodash"), col = c("black", "red", "blue"),
       y.intersp = 1.5)
dev.off()


#####################
# Illustration of Fabian's score


x <- seq(from = -1, to = 1, by = 0.001)
d <- 0.05
yl <- c(0, 4)
pi <- 0.2
epsilon <- 0.048
shift <- 0.05

new_score <- function(x, y, d, pi){
  (x - y)^2 + pi*(log((1 + exp(y/d))/(1 + exp(x/d))) - (y - x)/(1 + exp(-x/d))/d)
}

pdf(paste0("figures/new_score_", d, ".pdf"), width = 7, height = 2.7)
layout(matrix(c(1, 2, 3, 4, 4, 4), nrow = 2, byrow = TRUE), heights = c(0.85, 0.15))
par(mar = c(4, 4.4, 1, 1))
y <- 0
sc <- new_score(x, y, d, pi)
se <- (x - y)^2
TA2 <- TADDA_L2_v1(x, y, epsilon = epsilon)

plot(x, sc + shift, type = "l", xlab = expression(hat(y)), ylab = expression(s(hat(y), y)), ylim = yl,
     col = "orange", lty = 6)
lines(x, se - shift, col = "black", lty = 3)
lines(x, TA2, col = "blue", lty = 4)
abline(v = y, col = "darkgrey")
text_in_box(y, 3, "y", col = "darkgrey")
text_in_box(0, 2, "0", col = "darkgrey")


y <- 0.2
sc <- new_score(x, y, d, pi)
se <- (x - y)^2
TA2 <- TADDA_L2_v1(x, y, epsilon = epsilon)

plot(x, sc + shift, type = "l", xlab = expression(hat(y)), ylab = expression(s(hat(y), y)), ylim = yl,
     col = "orange", lty = 6)
lines(x, se - shift, col = "black", lty = 3)
lines(x, TA2, col = "blue", lty = 4)
abline(v = y, col = "darkgrey")
text_in_box(y, 3, "y", col = "darkgrey")

abline(v = 0, col = "darkgrey", lty = 3)
text_in_box(0, 2, "0", col = "darkgrey")


y <- 0.8
sc <- new_score(x, y, d, pi)
se <- (x - y)^2
TA2 <- TADDA_L2_v1(x, y, epsilon = epsilon)

plot(x, sc + shift, type = "l", xlab = expression(hat(y)), ylab = expression(s(hat(y), y)), ylim = yl,
     col = "orange", lty = 6)
lines(x, se - shift, col = "black", lty = 3)
lines(x, TA2, col = "blue", lty = 4)
abline(v = y, col = "darkgrey")
text_in_box(y, 3, "y", col = "darkgrey")

abline(v = 0, col = "darkgrey", lty = 3)
text_in_box(0, 2, "0", col = "darkgrey")

par(mar = c(0, 0, 0, 0))
plot(NULL, xlab = "", ylab = "", axes = FALSE, xlim = 0:1, ylim = 0:1)
legend("center", legend = c("suggested score", expression(TADDA[0.048]~with~L2~distance), "SE"),
       col = c("orange", "blue", "black"), lty = c(6, 4, 3), ncol = 3, bty = "n")

dev.off()

# test new score:
samples_y <- rsn(100000, xi = xi, omega = omega, alpha = alpha)
mu <- mean(samples_y)

new_sc <- numeric(length(x))
for(i in seq_along(x)){
  new_sc[i] <- mean(new_score(x = x[i], y = samples_y, d = 0.1, pi = 1))
}
plot(x, new_sc, type = "l")
abline(v = mu)
x[which.min(new_sc)]
