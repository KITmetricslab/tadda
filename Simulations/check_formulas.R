### Checking the formulas for Bayes acts against simulations

# Here we obtain OPFs / Bayes acts numerically and using the analytical formulas
# The agreement across different values of epsilon confirms the correctness
# of the formulas.

# note on nomenclature: this script uses the abbreviation "ba" ("Bayes act") 
# rather than "opf" ("optimal point forecast")

# Johannes Bracher
# johannes.bracher@kit.edu

library(sn) # library for skew normal

# set path to current directory (can also be done manually)
current_path <- rstudioapi::getActiveDocumentContext()$path # get path of this file
setwd(dirname(current_path))

# get helper functions
source("functions.R")
source("../Empirical Example/bayes_acts_functions.R")

# parameters of skew normal:
xi <- -0.15
omega <- 0.4
alpha <- 8

# generate samples from the specified skew normal (shorter than usually to speed up computations):
set.seed(123)
samples_y_short <- rsn(30000, xi = xi, omega = omega, alpha = alpha)
med <- median(samples_y_short)
mu <- mean(samples_y_short)

# grid of y-values to check for Bayes act numerically (grid search):
grid_y <- seq(from = 0.01, to = 0.3, by = 0.002)

# values of epsilon for which to compute Bayes acts:
grid_epsilon <- (seq(from = sqrt(0.01), to = sqrt(0.5), length.out = 50))^2


# compute Bayes acts *numerically* under different scores:
ba_num_tadda1_l1 <- get_bayes_acts(grid_y = grid_y, grid_epsilon = grid_epsilon,
                                   samples_y = samples_y_short, score = TADDA_L1_v1)
ba_num_tadda2_l1 <- get_bayes_acts(grid_y = grid_y, grid_epsilon = grid_epsilon,
                                  samples_y = samples_y_short, score = TADDA_L1_v2)
ba_num_tadda1_l2 <- get_bayes_acts(grid_y = grid_y, grid_epsilon = grid_epsilon,
                                  samples_y = samples_y_short, score = TADDA_L2_v1)


# compute Bayes acts *analytically* under different scores:
# TADDA1, L1, version used in manuscript
ba_ana_tadda1_l1 <- sapply(grid_epsilon, OPF_TADDA1_L1, distribution_Y = samples_y_short)
ba_ana_tadda2_l1 <- sapply(grid_epsilon, OPF_TADDA2_L1, distribution_Y = samples_y_short)
ba_ana_tadda1_l2 <- sapply(grid_epsilon, OPF_TADDA1_L2, distribution_Y = samples_y_short)


# Plot to compare curves:
pdf("Figures/ba_numerical_vs_analytical.pdf", width = 7, height = 2.6)

# structure plot region:
layout(matrix(1:4, nrow = 1), widths = c(2, 2, 2, 1))
par(las = 1, mar = c(4.2, 4.2, 2.5, 2))
# manually preparing square-root axes...
x_labels <- c(0.01, 0.1, 0.25, 0.5, 1, 2, 4)
y_labels <- c(0.01, 0.1, 0.25, 0.5, 1)

# TADDA1, L1:
plot(sqrt(grid_epsilon), sqrt(ba_num_tadda1_l1), type = "l", col = "red", lty = "dashed",
     xlab = expression(epsilon), ylab = "Bayes act", ylim = sqrt(c(0.01, 0.3)), axes = FALSE)
abline(0:1, col = "darkgrey")
mtext("(a) TADDA1, L1", side = 3, line = 0.5)
axis(1, at = sqrt(x_labels), labels = x_labels)
axis(2, at = sqrt(x_labels), labels = x_labels)


abline(h = sqrt(med), col = "darkgrey")
abline(v = sqrt(med), col = "darkgrey")

text_in_box(0.25, sqrt(med), "m", col = "darkgrey")
text_in_box(sqrt(med), 0.25, "m", col = "darkgrey")

text_in_box(sqrt(0.28), sqrt(0.28), expression(epsilon), col = "darkgrey")

lines(sqrt(grid_epsilon), sqrt(ba_ana_tadda1_l1), col = "blue", lty = "solid")
lines(sqrt(grid_epsilon), sqrt(ba_num_tadda1_l1), col = "red", lty = "dashed")


box()


# TADDA2, L1:
plot(sqrt(grid_epsilon), sqrt(ba_num_tadda2_l1), type = "l", col = "red", lty = "dashed",
     xlab = expression(epsilon), ylab = "Bayes act", ylim = sqrt(c(0.01, 0.3)), axes = FALSE)
abline(0:1, col = "darkgrey")
mtext("(b) TADDA2, L1", side = 3, line = 0.5)
axis(1, at = sqrt(x_labels), labels = x_labels)
axis(2, at = sqrt(x_labels), labels = x_labels)


abline(h = sqrt(med), col = "darkgrey")
abline(v = sqrt(med), col = "darkgrey")

text_in_box(0.25, sqrt(med), "m", col = "darkgrey")
text_in_box(sqrt(med), 0.25, "m", col = "darkgrey")

text_in_box(sqrt(0.28), sqrt(0.28), expression(epsilon), col = "darkgrey")

lines(sqrt(grid_epsilon), sqrt(ba_ana_tadda2_l1), col = "blue", lty = "solid")
lines(sqrt(grid_epsilon), sqrt(ba_num_tadda2_l1), col = "red", lty = "dashed")


box()

# TADDA1, L2:
plot(sqrt(grid_epsilon), sqrt(ba_num_tadda1_l2), type = "l", col = "red", lty = "dashed",
     xlab = expression(epsilon), ylab = "Bayes act", ylim = sqrt(c(0.01, 0.3)), axes = FALSE)
abline(0:1, col = "darkgrey")
mtext("(c) TADDA1, L2", side = 3, line = 0.5)
axis(1, at = sqrt(x_labels), labels = x_labels)
axis(2, at = sqrt(x_labels), labels = x_labels)

abline(h = sqrt(mu), col = "darkgrey")
abline(v = sqrt(mu), col = "darkgrey")
text_in_box(0.25, sqrt(mu), expression(mu), col = "darkgrey")
text_in_box(sqrt(mu), 0.25, expression(mu), col = "darkgrey")
text_in_box(sqrt(0.28), sqrt(0.28), expression(epsilon), col = "darkgrey")
box()

lines(sqrt(grid_epsilon), sqrt(ba_ana_tadda1_l2), col = "blue", lty = "solid")
lines(sqrt(grid_epsilon), sqrt(ba_num_tadda1_l2), type = "l", col = "red", lty = "dashed")

# Legend:
par(mar = c(0, 0, 0, 0))
plot(NULL, xlim = 0:1, ylim = 0:1, xlab = "", ylab = "", axes = FALSE)
legend("center", legend = c("analytical", "numerical"),
       bty = "n", cex = 1, lty = c("solid", "dashed"), 
       col = c("red", "blue"),
       y.intersp = 1.5)

dev.off()