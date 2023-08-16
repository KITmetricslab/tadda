### Some illustrations and simulation examples on the TADDA scores

# This script generates FIGURE 1, FIGURE 2, TABLE 1, SUPPLEMENTARY FIGURE S6

# note on nomenclature: this script uses the abbreviation "ba" ("Bayes act") 
# rather than "opf" ("optimal point forecast")

# Johannes Bracher
# johannes.bracher@kit.edu

library(sn) # library for skew normal
library(plotrix) # library needed for some plotting functionality
library(xtable) # required for tex tables


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

# tolerance value:
epsilon <- 0.048

# generate samples from the specified skew normal:
set.seed(123)
samples_y <- rsn(100000, xi = xi, omega = omega, alpha = alpha)

# compute mean and median of the skew normal, i.e., OPFs under squared and absolute error:
(mu <- mean(samples_y))
(med <- median(samples_y))

# compute TADDA Bayes acts:
grid_y <- seq(from = 0, to = 0.2, by = 0.001) # requires a grid of values to search
# compute numerically:
(ba_tadda1 <- get_bayes_acts(grid_y = grid_y, grid_epsilon = epsilon, samples_y = samples_y, score = TADDA_L1_v1))
# for comparison: analytical solution:
OPF_TADDA1_L1(samples_y, epsilon)
# compute numerically:
(ba_tadda2 <- get_bayes_acts(grid_y = grid_y, grid_epsilon = epsilon, samples_y = samples_y, score = TADDA_L1_v2))
# for comparison: analytical solution:
OPF_TADDA2_L1(samples_y, epsilon)
# analytical and numerical solutions agree

#####################
### FIGURE 2 from the manuscript:

# grid of values for y (shown on the x-axis of Fig 2):
grid_y_density <- c(seq(from = -0.5, to = 1.1, by = 0.001))

# evaluate density:
f <- dsn(x = grid_y_density, xi = xi, omega = omega, alpha = alpha)


# compute expected scores as a function of y_hat:
# grid: more dense where minima occur, less dense for rest
grid_y_hat <- c(seq(from = -0.5, to = 0, by = 0.01),
                seq(from = 0.001, to = 0.3, by = 0.001),
                seq(from = 0.31, to = 1.1, by = 0.01))

#empty vectors to store expected scores:
average_scores_ae <-
  average_scores_se <-
  average_scores_tadda_l1_0 <-
  average_scores_tadda_l1_v1_epsilon <-
  average_scores_tadda_l1_v2_epsilon <-
  average_scores_tadda_l2_0 <-
  average_scores_tadda_l2_v1_epsilon <-
  average_scores_tadda_l2_v2_epsilon <-
  numeric(length(grid_y_hat))

# fill these vectors:
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

# plot:
pdf("Figures/illustration.pdf", width = 8.5, height = 3)
# structure plot area:
layout(matrix(1:2, ncol = 2), widths = c(0.68, 0.32))
par(mar = c(4.2, 4.2, 0.5, 4.5), las = 1)

# plot density:
plot(grid_y_density, f, type = "l", ylab = "f(y)", xlab = "", ylim = c(0, 2))

# custom x-axis labelling:
mtext(expression(y), side = 1, line = 2.5, at = 0.2)
mtext(expression(hat(y)), side = 1, line = 2.5, at = 0.3, col = "darkgrey")

# add second axis (right):
axis(4, at = 0:4/2, labels = paste0(0:4, ".0"), col.axis = "darkgrey", col.ticks = "darkgrey")
par(las = 0)
mtext(side = 4, at = 1, line = 3, "expected score", col = "darkgrey")
par(las = 1)
# # show tolerance region (commented out as plot gets busy)
# rect(-epsilon, 0, epsilon, c(-1, 2.5), col = "grey97", border = NA)
# abline(v = epsilon, lty = 3, col = "grey")
# abline(v = -epsilon, lty = 3, col = "grey")

box()

# light lines with expected scores:
lines(grid_y_hat, 2*average_scores_tadda_l1_v1_epsilon, type = "l", col = rgb(1, 0, 0, 0.25), lty = 4)
lines(grid_y_hat, 2*average_scores_ae, type = "l", col = "darkgrey", lty = "dotted")
lines(grid_y_hat, 2*average_scores_se, type = "l", col = "darkgrey", lty = "dashed")

# add density on top:
lines(grid_y_density, f)

# vertical lines indicating OPFs:
abline(v = mu, lty = 2, col = "black")
abline(v = med, lty = 3, col = "black")
abline(v = ba_tadda1, lty = 4, col = "red")

# extra panel for legend:
par(mar = c(0, 0, 0, 0), las = 1)
plot(NULL, xlim = 0:1, ylim = 0:1, axes = FALSE, xlab = "", ylab = "")
legend("center", legend = c("Functionals of F:",
                            paste0("Mean: ", round(mu, 3)), # hard-coded, modify as needed
                              paste0("Median: ", round(med, 3)),
                              expression(OPF~TADDA[0.048]: 0.06),
                              "",
                              expression(Expected~score~ of~hat(y)~under~"F:"),
                              "SE",
                              "AE",
                              expression(TADDA[0.048])),
                              # expression(BA~under~TADDA[0]: 0.019)),
       lty = c(NA, 2, 3, 4, NA, NA, 2, 3, 4), # bty = "n",
       col = c(NA, "black", "black", "red", NA,
               NA, "darkgrey", "darkgrey", rgb(1, 0, 0, 0.25)), cex = 0.9, bty = "n")
dev.off()



#####################
### additional figure ultimately not used in manuscript: plot expected scores as function of y_hat:

pdf("Figures/expected_scores.pdf", width = 7, height = 2.6)
# structure plot area:
layout(matrix(1:3, nrow = 1), widths = c(2, 2, 1))
par(las = 1, mar = c(4.2, 4.2, 2.5, 2))
yl <- c(0.1, 0.4)
cex.pt <- 2

# L1
plot(grid_y_hat, average_scores_tadda_l1_0, type = "l", col = "black", xlab = expression(hat(y)),
     ylab = "expected score", xlim = c(-0.1, 0.2), ylim = yl)
mtext("(a) L1", side = 3, line = 0.5)

# highlight tolerance area:
rect(-epsilon, 0, epsilon, 1.2*yl[2], col = "grey97", border = NA)

abline(v = 0, col = "darkgrey", lty = "dotted")
text_in_box(0, 0.9*yl[2], 0, col = "darkgrey")

abline(v = epsilon,  col = "darkgrey", lty = "dotted")
text_in_box(epsilon, 0.9*yl[2], expression(epsilon), col = "darkgrey")

abline(v = -epsilon,  col = "darkgrey", lty = "dotted")
text_in_box(-epsilon, 0.9*yl[2], expression(-epsilon), col = "darkgrey")

# highlight median:
abline(v = med, col = "darkgrey")
text_in_box(med, 0.36, "m", col = "darkgrey")

# add curves:
lines(grid_y_hat, average_scores_tadda_l1_0, type = "l", col = "black", lty = "solid")
lines(grid_y_hat, average_scores_tadda_l1_v1_epsilon, type = "l", col = "red", lty = "dashed")
lines(grid_y_hat, average_scores_tadda_l1_v2_epsilon, type = "l", col = "blue", lty = "twodash")
lines(grid_y_hat, average_scores_ae, type = "l", col = "darkgreen", lty = "dotted")

# add points for minima:
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

# highlight tolerance area:
rect(-epsilon, 0, epsilon, 1.2*yl2[2], col = "grey97", border = NA)

abline(v = 0, col = "darkgrey", lty = "dotted")
text_in_box(0, 0.9*yl2[2], 0, col = "darkgrey")

abline(v = epsilon,  col = "darkgrey", lty = "dotted")
text_in_box(epsilon, 0.9*yl2[2], expression(epsilon), col = "darkgrey")

abline(v = -epsilon,  col = "darkgrey", lty = "dotted")
text_in_box(-epsilon, 0.9*yl2[2], expression(-epsilon), col = "darkgrey")

# highligh mu:
abline(v = mu,  col = "darkgrey")
text_in_box(mu, 0.9*yl2[2], expression(mu), col = "darkgrey")

# add curves:
lines(grid_y_hat, average_scores_tadda_l2_0, type = "l", col = "black", lty = "solid")
lines(grid_y_hat, average_scores_tadda_l2_v1_epsilon, type = "l", col = "red", lty = "dashed")
lines(grid_y_hat, average_scores_tadda_l2_v2_epsilon, type = "l", col = "blue", lty = "twodash")
lines(grid_y_hat, average_scores_se, type = "l", col = "purple", lty = "dotted")

# add points for minima:
points(grid_y_hat[which.min(average_scores_tadda_l2_0)], min(average_scores_tadda_l2_0),
       col = "black", pch = 18, cex = cex.pt)
points(grid_y_hat[which.min(average_scores_tadda_l2_v1_epsilon)], min(average_scores_tadda_l2_v1_epsilon),
       col = "red", pch = 18, cex = cex.pt)
points(grid_y_hat[which.min(average_scores_tadda_l2_v2_epsilon)], min(average_scores_tadda_l2_v2_epsilon),
       col = "blue", pch = 18, cex = cex.pt)
points(grid_y_hat[which.min(average_scores_se)], min(average_scores_se),
       col = "purple", pch = 18, cex = cex.pt)

box()

# Legend in separate panel:
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
### TABLE 1: Table with expected scores of Bayes acts:

# matrix tro store expected scores:
expected_scores <- matrix(ncol = 9, nrow = 9)
colnames(expected_scores) <- c("value", "AE", "SE", 
                               "TADDA_L1_0", "TADDA_L1_v1", "TADDA_L1_v2",
                               "TADDA_L2_0", "TADDA_L2_v1", "TADDA_L2_v2")
rownames(expected_scores) <- c("median", "mean",
                               "BA TADDA_L1_0", "BA TADDA_L1_v1", "BA TADDA_L1_v2",
                               "BA TADDA_L2_0", "BA TADDA_L2_v1", "BA TADDA_L2_v2",
                               "zero")

# extract OPFs from previously cimputed vectors:
ba_tadda_l1_v1_epsilon <- grid_y_hat[which.min(average_scores_tadda_l1_v1_epsilon)]
ba_tadda_l1_v2_epsilon <- grid_y_hat[which.min(average_scores_tadda_l1_v2_epsilon)]
ba_tadda_l2_v1_epsilon <- grid_y_hat[which.min(average_scores_tadda_l2_v1_epsilon)]
ba_tadda_l2_v2_epsilon <- grid_y_hat[which.min(average_scores_tadda_l2_v2_epsilon)]

# add OPFs to matrix:
expected_scores[, "value"] <- c(med, mu, 
                                ba_tadda1, ba_tadda_l1_v1_epsilon, ba_tadda_l1_v2_epsilon,
                                ba_tadda2, ba_tadda_l2_v1_epsilon, ba_tadda_l2_v2_epsilon,
                                0)

expected_scores[, "AE"] <- sapply(expected_scores[, "value"], 
                                  function(y_hat) mean(abs(y_hat - samples_y)))
expected_scores[, "SE"] <- sapply(expected_scores[, "value"], 
                                  function(y_hat) mean((y_hat - samples_y)^2))

# add expected scores to matrix:
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

# use only a subset of the generated matrix for the manuscript:
expected_scores_small <- rbind(expected_scores[c("median", "mean", "BA TADDA_L1_v1", "BA TADDA_L1_v2"), ],
                               zero = expected_scores["zero", ])
expected_scores_small <- expected_scores_small[, c("AE", "SE", "TADDA_L1_v1", "TADDA_L1_v2")]

# generate latex code: expected scores
xtable(expected_scores_small, 3, digits = 3)

# OPFs:
xtable(matrix(expected_scores[1:8, "value"], nrow = 1), digits = 3)



#########################
### additional figure ultimately not used in manuscript: plotting Bayes act as a function of epsilon:

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
pdf("Figures/ba_epsilon.pdf", width = 7, height = 2.6)

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
# FIGURE 1: Illustration of TADDA1 score:

# set epsilon
epsilon <- 0.048

### L1, version TADDA2
pdf("Figures/curves_scores_L1.pdf", width = 8, height = 4)
# structure plot area:
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
  TADDA~with~y%in%paste("[", -epsilon, ",", epsilon, "]"~", varying"~hat(y))
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
grid_TADDA_v2_c <- TADDA_L1_v2(y_hat, grid_y, epsilon = epsilon)
grid_ae_c <- abs(y_hat - grid_y)

plot(grid_y, grid_TADDA_v1_c + shift, type = "l", xlab = expression(y),
     ylab = ylab, ylim = ylim, xlim = xlim, col = "white")
mtext(
  expression(TADDA~with~hat(y)%in%paste("[", -epsilon, ",", epsilon, "]"~", varying"~y)),
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

# legend in separate panel:
par(mar = c(0, 0, 0, 0))
plot(NULL, xlim = 0:1, ylim = 0:1, xlab = "", ylab = "", axes = FALSE)
legend("center", legend = c(# expression({TADDA[0]}(hat(y), y)),
  expression({TADDA2[epsilon]}(hat(y), y)),
  expression(AE(hat(y), y))),
  bty = "n", cex = 1.2, lty = 5:3, col = c("red", "black"),
  y.intersp = 1.5)
dev.off()



#################
# SUPPLEMENTARY FIGURE S6: Illustration of TADDA2 score:

# set epsilon
epsilon <- 0.048

### L1, version TADDA2
pdf("Figures/illustration_TADDA2.pdf", width = 8, height = 4)
# structure plot area:
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

plot(grid_y, grid_TADDA_v2_a, xlab = expression(hat(y)),
     ylab = ylab, ylim = ylim, xlim = xlim, col = "white")
mtext(expression(
  TADDA~with~y%in%paste("[", -epsilon, ",", epsilon, "]"~", varying"~hat(y))
), 3, cex = 0.9, line = 0.3)

rect(-epsilon, -1, epsilon, 11, col = "grey97", border = NA)
abline(h = 0:10, col = "grey95")
abline(v = -6:8, col = "grey95")
box()
# lines(grid_y, grid_TADDA_a + shift, col = "red", lty = 5)
lines(grid_y, grid_TADDA_v2_a, col = "red", lty = 4)
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


plot(grid_y, grid_TADDA_v2_b + shift, type = "l", xlab = expression(hat(y)),
     ylab = "", ylim = ylim, xlim = xlim, col = "white")
mtext(expression(TADDA~with~y > epsilon~", varying"~hat(y)), 3, cex = 0.9, line = 0.3)
rect(-epsilon, -1, epsilon, 11, col = "grey97", border = NA)
abline(h = 0:10, col = "grey95")
abline(v = -6:8, col = "grey95")
box()
# lines(grid_y, grid_TADDA_b + shift, col = "red", lty = 5)
lines(grid_y, grid_TADDA_v2_b, col = "red", lty = 4)
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
grid_TADDA_v2_c <- TADDA_L1_v2(y_hat, grid_y, epsilon = epsilon)
grid_ae_c <- abs(y_hat - grid_y)

plot(grid_y, grid_TADDA_v2_c + shift, type = "l", xlab = expression(y),
     ylab = ylab, ylim = ylim, xlim = xlim, col = "white")
mtext(
  expression(TADDA~with~hat(y)%in%paste("[", -epsilon, ",", epsilon, "]"~", varying"~y)),
  3, cex = 0.9, line = 0.3)

rect(-epsilon, -1, epsilon, 11, col = "grey97", border = NA)
abline(h = 0:10, col = "grey95")
abline(v = -6:8, col = "grey95")
box()
# lines(grid_y, grid_TADDA_c + shift, col = "red", lty = 5)
lines(grid_y, grid_TADDA_v2_c, col = "red", lty = 4)
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


plot(grid_y, grid_TADDA_v2_d + shift, type = "l", xlab = expression(y),
     ylab = "", ylim = ylim, xlim = xlim, col = "white")
mtext(expression(TADDA~with~hat(y) > epsilon~", varying"~y), 3, cex = 0.9, line = 0.3)
rect(-epsilon, -1, epsilon, 11, col = "grey97", border = NA)
abline(h = 0:10, col = "grey95")
abline(v = -6:8, col = "grey95")
box()
# lines(grid_y, grid_TADDA_d + shift, col = "red", lty = 5)
lines(grid_y, grid_TADDA_v2_d, col = "red", lty = 4)
lines(grid_y, grid_ae_d - shift, col = "black", lty = 3)

abline(v = 0, col = "darkgrey", lty = "dotted")
abline(v = y_true, lty = "solid", col = "darkgrey")
abline(v = c(epsilon, -epsilon), lty = "dotted", col = "darkgrey")

text_in_box(0, 0.9*ylim[2], expression(0), col = "darkgrey")
text_in_box(epsilon, 0.9*ylim[2], expression(epsilon), col = "darkgrey")
text_in_box(-epsilon,  0.9*ylim[2], expression(-epsilon), col = "darkgrey")
text_in_box(y_true,  0.5*ylim[2], expression(hat(y)), col = "darkgrey")

# legend in separate panel:
par(mar = c(0, 0, 0, 0))
plot(NULL, xlim = 0:1, ylim = 0:1, xlab = "", ylab = "", axes = FALSE)
legend("center", legend = c(# expression({TADDA[0]}(hat(y), y)),
                            expression({TADDA2[epsilon]}(hat(y), y)),
                            expression(AE(hat(y), y))),
       bty = "n", cex = 1.2, lty = 5:3, col = c("red", "black"),
       y.intersp = 1.5)
dev.off()



### additional figure ultimately not used: illustration of L2 scores

pdf("Figures/curves_scores_L2.pdf", width = 8, height = 4)
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
  TADDA^(2)~with~y%in%paste("[", -epsilon, ",", epsilon, "]"~", varying"~hat(y))
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
  expression(TADDA^(2)~with~hat(y)%in%paste("[", -epsilon, ",", epsilon, "]"~", varying"~y)),
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