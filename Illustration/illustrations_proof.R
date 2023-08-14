### Illustrations for proofs of OPF formulations

# This script generates SUPPLEMENTARY FIGURES S4 and S5

# Johannes Bracher
# johannes.bracher@kit.edu

library(sn)

current_path <- rstudioapi::getActiveDocumentContext()$path # get path of this file
setwd(dirname(current_path))

source("functions.R")


#################################################
# SUPPLEMENTARY FIGURE S4: TADDA_epsilon, m > epsilon

# epsilon:
epsilon <- 0.5 # using larger epsilon for illustration

# parameters of skew normal:
xi <- -1.7
omega <- 5
alpha <- 8

# generate samples:
set.seed(123)
samples_y <- rsn(100000, xi = xi, omega = omega, alpha = alpha)

# prepare plotting density:
grid_y <- seq(from = -5, to = 10, by = 0.01)
F <- psn(x = grid_y, xi = xi, omega = omega, alpha = alpha)
# compute mean and median:
(mu <- mean(samples_y))
(med <- median(samples_y))

# evaluate density at epsilon and minus epsilon
F_minus_epsilon <-  mean(samples_y < -epsilon)
F_epsilon <-  mean(samples_y < epsilon)

# inflate with epsilons:
n_epsilon <- sum(samples_y < -epsilon)
samples_y_epsilon <- c(samples_y, rep(epsilon, n_epsilon))

# prepare plotting density of epsilon-inflated distribution:
grid_p <- 1:100/100
quantiles_y_epsilon <- quantile(samples_y_epsilon, p = grid_p)

# obtain median and mean of inflated distribution:
(med_modified <- max((quantile(samples_y, 0.5*(1 - sign(med)*F_minus_epsilon))), 0))
(mu_modified <- mean(samples_y)/(1 + pi))


# Plot:
pdf("figures/F_vs_G_epsilon.pdf", width = 6, height = 3.5)
par(las = 1, mar = c(4.2, 4.2, 0.5, 0.5))
plot(grid_y, F, type = "l", col = "black", xlab = "y", ylab = "cumulative density", ylim = 0:1)

# highlight tolerance region:
rect(-epsilon, 0, epsilon, 1.2, col = "grey97", border = NA)

# adding vertical / horizontal lines and labelling:
abline(v = -epsilon, col = "lightgrey")
abline(v = epsilon, col = "lightgrey")
abline(h = 0.5, col = "lightgrey", lty  ="dashed")
abline(h = F_epsilon/(1 + F_minus_epsilon), col = "lightgrey")
abline(h = (F_epsilon + F_minus_epsilon)/(1 + F_minus_epsilon), col = "lightgrey")
abline(v = med, col = "black", lty = "dashed")
abline(v = med_modified, col = "red", lty  ="dashed")

# adding CDFs:
lines(grid_y, F)
lines(quantiles_y_epsilon, grid_p, col = "red")

# some more labelling:
text_in_box(-3, 0.5, "0.5", col = "lightgrey")
text_in_box(epsilon, 0.7, expression(epsilon), col = "lightgrey")
text_in_box(-epsilon, 0.7, expression(-epsilon), col = "lightgrey")
text_in_box(med_modified, 0.95, "OPF: median of Z", col = "red")
text_in_box(med + 1.25, 0.05, "median of Y", col = "black")
text_in_box(6, (F_epsilon + F_minus_epsilon)/(1 + F_minus_epsilon), 
            expression((F(epsilon) + pi["-"])/(1 + pi["-"])), col = "lightgrey", cex = 0.8)
text_in_box(6, F_epsilon/(1 + F_minus_epsilon), 
            expression(F(epsilon)/(1 + pi["-"])), col = "lightgrey", cex = 0.8)

legend("topleft", legend = c("CDF of Y", "CDF of Z"), col = c("black", "red"), lty = 1, bty = "n")

box()

dev.off()



##########################################
# SUPPLEMENTARY FIGURE S5: TADDA_epsilon, m > epsilon

# epsilon:
epsilon <- 0.5 # larger epsilon for illustration:

# parameters of skew normal:
xi <- 1.7
omega <- 5
alpha <- -8

# generate samples:
set.seed(123)
samples_y <- rsn(100000, xi = xi, omega = omega, alpha = alpha)

# prepare plotting density:
grid_y <- seq(from = -10, to = 5, by = 0.01)
F <- psn(x = grid_y, xi = xi, omega = omega, alpha = alpha)
# compute mean and median:
(mu <- mean(samples_y))
(med <- median(samples_y))

# evaluate density at epsilon and minus epsilon
F_minus_epsilon <-  mean(samples_y < -epsilon)
F_epsilon <-  mean(samples_y < epsilon)

# inflate with -epsilons:
n_minus_epsilon <- sum(samples_y > epsilon)
samples_y_minus_epsilon <- c(samples_y, rep(-epsilon, n_minus_epsilon))

# prepare plotting density of -epsilon-inflated distribution:
grid_p <- 1:100/100
quantiles_y_minus_epsilon <- quantile(samples_y_minus_epsilon, p = grid_p)

# obtain median and mean of inflated distribution:
(med_modified <- quantile(samples_y, 0.5*(2 - F_epsilon)))
(mu_modified <- mean(samples_y)/(1 + pi))

# Plot:

pdf("figures/F_vs_G_minus_epsilon.pdf", width = 6, height = 3.5)
par(las = 1, mar = c(4.2, 4.2, 0.5, 0.5))
plot(grid_y, F, type = "l", col = "black", xlab = "y", ylab = "cumulative density", ylim = 0:1)

# highlight tolerance region:
rect(-epsilon, 0, epsilon, 1.2, col = "grey97", border = NA)

# various vertical and horizontal lines with labels:
abline(v = -epsilon, col = "lightgrey")
abline(v = epsilon, col = "lightgrey")
abline(h = 0.5, col = "lightgrey", lty = "dashed")
abline(h = F_minus_epsilon/(2 - F_epsilon), col = "lightgrey")
abline(h = (1 - F_epsilon + F_minus_epsilon)/(2 - F_epsilon), col = "lightgrey")
abline(v = med, col = "black", lty = "dashed")
abline(v = med_modified, col = "red", lty  ="dashed")

# add CDFs:
lines(quantiles_y_minus_epsilon, grid_p, col = "red")
lines(grid_y, F)

# some more labelling:
text_in_box(3, 0.5, "0.5", col = "lightgrey")
text_in_box(epsilon, 0.3, expression(epsilon), col = "lightgrey")
text_in_box(-epsilon, 0.3, expression(-epsilon), col = "lightgrey")
text_in_box(-6, F_minus_epsilon/(2 - F_epsilon), 
            expression(F(-epsilon)/(1 + pi["+"])), col = "lightgrey", cex = 0.8)
text_in_box(-5.5, (1 - F_epsilon + F_minus_epsilon)/(2 - F_epsilon),
            expression((F(-epsilon) + pi["+"])/(1 + pi["+"])), col = "lightgrey", cex = 0.8)
text_in_box(med_modified + 1.6, 0.15, "BA: median of Z", col = "red")

text_in_box(med, 0.95, "median of Y", col = "black")


# legend:
legend("topleft", legend = c("CDF of Y", "CDF of Z"), col = c("black", "red"), lty = 1, bty = "n")

box()

dev.off()
