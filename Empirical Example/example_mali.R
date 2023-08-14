### Code for empirical example illustrating properties of forecasts for log changes of armed conflict fatalities using TADDA variants, MAE or MSE

# This generates FIGURE 3

# Johannes Bracher
# johannes.bracher@kit.edu


# load packages
# library(dplyr) # included in tidyverse
library(tidyverse) # required for "reduce" function
library(rlist) # for data manipulation
library(purrr) # for data manipulation

# set working directory (can also be done manually)
current_path = rstudioapi::getActiveDocumentContext()$path # get path of this file
setwd(dirname(current_path))

# get scoring functions
source("../bayes_acts_functions.R")

# read in data
data_fatalities <- read.csv(paste("../Data/fatalities.csv"))[,-1]
country_names <- unique(data_fatalities$country_name)

# define week when forecast is issued:
end <- 470
# data subset to plot:
mali <- subset(data_fatalities, country_name == "Mali" & month_id %in% c(445:end))
# time on calendar scale:
mali$time <- mali$year + mali$month/12
# extract last value (needed in computations:)
last_value <- tail(mali$fatalities, 1)
yl <- c(5, 115)


# Plot:
pdf("figures/example_mali.pdf", width = 9, height = 3.5)

# structure plot area:
layout(matrix(1:2, ncol = 2), widths = c(0.6, 0.4))
par(las = 1, mar = c(2.5, 4, 2.5, 5))


# time series plot:
plot(mali$time, mali$fatalities, type = "h", xlim = c(2017, 2019.5), axes = FALSE,
     xlab = "", ylab = "fatalities", ylim = yl,
     main = "", log = "y", col = "white")
mtext("Time series of fatalities", line = 1)
# add axes manually:
# left:
axis(1, at = c(2017 + c(1, 4, 7, 10, 13, 16, 19, 22, 25, 28, 31)/12),
     labels = c("Jan 17", "", "Jul 17", "", "Jan 18", "", "Jul 18", "", "Jan 19", "", "Jul 19"))
axis(2, at = c(5, 10, 20, 50, 100))
# right:
ylabs_right <- -4:2
at_right <- last_value*exp(ylabs_right)
axis(4, at = at_right, labels = ylabs_right)
mtext("log-difference to observation Feb 2019", side = 4, las = 0, line = 2.5)
box()

# highlightinh and labelling:
now <- tail(mali$time, 1)
abline(v = now + 0.2/12, lty = 3)
rect(now - 8.2/12, 0.01, now + 0.2/12, 150, col = rgb(0.2, 0.2, 0.2, 0.2), border = NA)
text(2018.8, 90, "last 9\n observations", cex = 0.75)
text(now + 1.3/12, 25, "last available obsrevation: Feb 2019", srt = 90, cex = 0.75)
abline(h = last_value, col = "darkgrey", lty = "solid")

# add points for last obsrevations:
points(tail(mali$time, 9), tail(mali$fatalities, 9), pch = 16, cex = 0.75)
points(tail(mali$time, 1), tail(mali$fatalities, 1), pch = 4, cex = 1.5)
lines(mali$time, mali$fatalities, type = "h")

# Boxplot:
par(mar = c(2.5, 4, 2.5, 1))
reference_obs <- tail(mali$fatalities, 9)
last_obs <- tail(mali$fatalities, 1)
reference_log_changes <- log(reference_obs + 1) - log(last_obs + 1)
boxplot(reference_log_changes, col = "white", border = "grey", ylim = log(yl/last_value),
        ylab = "log change", main = "", axes = FALSE)
axis(2, at = c(-2, -1, 0, 1))
box()
mtext(side = 3, "Forecast distribution\n for log-change")

# compute OPFs:
mu <- mean(reference_log_changes)
TADDA_BA <- OPF_TADDA1_L1(reference_log_changes, epsilon = 0.048)
# add to plot:
abline(h = TADDA_BA, col = "red", lty = 4)
abline(h = mu, col = "black", lty = 2)
abline(h = 0, col = "darkgrey", lty = "solid")

# boxplot on top:
boxplot(reference_log_changes, col = "white", border = "grey", add = TRUE, axes = FALSE)
# add points for individual values:
points(rep(1, length(reference_log_changes)), reference_log_changes, pch = 16, cex = 0.75)
points(1, 0, pch = 4, cex = 1.2)

# legend:
legend("bottom", legend = c(paste("Mean:", round(mu, 2)),
                            expression(OPF~TADDA[0.048]: -0.048)),
       cex = 0.75, bty = "n", lty = c(2, 4), col = c("black", "red"))
dev.off()
