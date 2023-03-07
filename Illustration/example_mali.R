### Code for empirical example illustrating properties of forecasts for log changes of armed conflict fatalities using TADDA variants, MAE or MSE
# Johannes Bracher
# johannes.bracher@kit.edu


# load packages
# library(dplyr) # included in tidyverse
library(tidyverse) # required for "reduce" function
library(rlist)
library(purrr)

# set working directory
current_path = rstudioapi::getActiveDocumentContext()$path # get path of this file
setwd(dirname(current_path))

# get scoring functions
source("functions.R")
source("../bayes_acts_functions.R")

# read in data
data_fatalities <- read.csv(paste("../Data/fatalities.csv"))[,-1]
country_names <- unique(data_fatalities$country_name)

end <- 470
mali <- subset(data_fatalities, country_name == "Mali" & month_id %in% c(445:end))
mali$time <- mali$year + mali$month/12
last_value <- tail(mali$fatalities, 1)
yl <- last_value*exp(c(-2, 1))


# Plot:
pdf("figures/example_mali.pdf", width = 9, height = 3.5)

# structure plot area:
layout(matrix(1:2, ncol = 2), widths = c(0.6, 0.4))
par(las = 1, mar = c(2.5, 4, 2.5, 5))


# time series plot:
plot(mali$time, mali$fatalities, type = "h", xlim = c(2017, 2019.5), axes = FALSE,
     xlab = "", ylab = "fatalities (log-scale)", ylim = yl,
     main = "", log = "y")
mtext("Time series of fatalities", line = 1)
axis(1, at = c(2017 + c(1, 4, 7, 10, 13, 16, 19, 22, 25, 28, 31)/12),
     labels = c("Jan 17", "", "Jul 17", "", "Jan 18", "", "Jul 18", "", "Jan 19", "", "Jul 19"))
axis(2)
ylabs_right <- -4:2
at_right <- last_value*exp(ylabs_right)
axis(4, at = at_right, labels = ylabs_right)
mtext("log-difference to last observation", side = 4, las = 0, line = 2)
box()

now <- tail(mali$time, 1)
abline(v = now + 0.5/12, lty = 3)
rect(now - 8.5/12, 0.01, now + 0.5/12, 150, col = rgb(0.2, 0.2, 0.2, 0.2), border = NA)
text(2018.8, 90, "last 9\n observations", cex = 0.75)
text(now + 1.2/12, 45, "time of forecast", srt = 90, cex = 0.75)
points(tail(mali$time, 9), tail(mali$fatalities, 9), pch = 16, cex = 0.75)


# Boxplot:
par(mar = c(2.5, 4, 2.5, 1))
reference_obs <- tail(mali$fatalities, 9)
last_obs <- tail(mali$fatalities, 1)
reference_log_changes <- log(reference_obs + 1) - log(last_obs + 1)
boxplot(reference_log_changes, col = "white", border = "grey", ylim = log(yl/last_value),
        ylab = "log change", main = "")
mtext(side = 3, "Forecast distribution\n for log-change")

mu <- mean(reference_log_changes)
TADDA_BA <- BA_TADDA1_L1(reference_log_changes, epsilon = 0.048)

abline(h = TADDA_BA, col = "red", lty = 4)
abline(h = mu, col = "black", lty = 2)

boxplot(reference_log_changes, col = "white", border = "grey", add = TRUE)
points(rep(1, length(reference_log_changes)), reference_log_changes, pch = 16, cex = 0.75)

legend("bottom", legend = c(paste("Mean:", round(mu, 2)),
                            expression(BA~under~TADDA[0.048]: 0)),
       cex = 0.75, bty = "n", lty = c(2, 4), col = c("black", "red"))
dev.off()
