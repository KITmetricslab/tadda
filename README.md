# Code for the paper _Direction Augmentation in the Evaluation of Armed Conflict Predictions_

by Johannes Bracher, Lotta Rüter, Fabian Krüger, Sebastian Lerch and Melanie Schienle published in _International Interactions_, preprint available [here](https://arxiv.org/abs/2304.12108).

## Notes

+ The contents of this repository are under a Creative Commons Attribution-NonCommercial 4.0 International Public License.
+ Please note that we use the terms optimal point forecast (OPF) and Bayes act (BA) synonymously in this code.
+ All computations were performed with R version 4.3.0 (2023-04-21)
+ Required packages: `tidyverse`, `rlist`, `purrr`, `sn`, `plotrix`, `xtable`

## Contents

### `Empirical Example` – data and code to reproduce all empirical results
+ `tadda_example.R` reproduces the empirical results from Section 5, generates the files in the `Results` folder, computes the optimal window size *w* and generates Tables 2 and 3.
+ `example_mali.R` generates Figure 3.
+ `bayes_acts_functions.R` contains functions for computing the optimal point forecasts (OPFs) of different scoring functions (AE, SE and variants of TADDA) as well as a summary function for a nicer representation of the results.

#### `/Data`
+ `data_prep.R` (i) extracts time series of country month fatalities due to state based conflict `ged_best_sb` for each country in the data set and (ii) computes true *s = 1, ..., 7* step ahead log-changes for each country and month.
+ `fatalities.csv` contains results of `data_prep.R`.
+ `ged_cm_postpatch.parquet` and `skeleton_cm_africa.parquet` contain the data in their original format. These were retrieved from https://github.com/UppsalaConflictDataProgram/views_competition/tree/main/data (published under Creative Commons Attribution-NonCommercial 4.0 International Public License).

#### `/Results`
+ `average_scores_for_different_window_lengths.csv` contains average scores by forecast horizon for different values of *w*; shows that *w = 5* would be optimal for minimizing TADDA1 via `TADDA1_OPF`.
+ `individual_predictions_w9.csv` and `individual_losses_w9.csv` contain the predictions using *w = 9* for the log-changes in fatalities and corresponding losses for each African country, month in 395:495, lead time *s = 2, ..., 7* and scoring function / OPF.
+ `average_scores_w9.csv` contains the central results presented in Table 2 (i.e., with *w = 9*).
+ `empirical_quantiles_w9.csv` is the basis of Table 3.

#### `/Figures`
+ Figure 3: `example_mali.pdf`

### `Simulations` – illustrative figures and small simulation examples
+ `check_formulas.R` compares OPFs / Bayes acts computed using numerical optimization to the analytical results obtained using the formulas from the manuscript. The agreement between the two indicates that our derivations are correct.
+ `functions.R` contains functions (i) to compute different scoring functions (AE, SE and variants of TADDA), (ii) to numerically determine the optimal point forecast (Bayes act) given a distribution and scoring rule as well as (iii) to provide text annotations in a plot.
+ `illustration.R` generates Figure 2, Table 1, and Supplementary Figure S6.
+ `illustrations_proof.R` generates Supplementary Figures S4 and S5.

#### `/Figures`
+ Figure 1: `curves_scores_L1.pdf`
+ Figure 2: `illustration.pdf`
+ Figure S4: `F_vs_G_epsilon.pdf`
+ Figure S5: `F_vs_G_minus_epsilon.pdf`
+ Figure S6: `illustration_TADDA2.pdf`
+ `ba_numerical_vs_analytical.pdf` (plausibility check of analytical results), `curves_scores_L2.pdf` (illustration of the L2 version of TADDA1) and `expected_scores.pdf` (more detailed version of light grey / red curves from Figure 2) further illustrate TADDA1 and TADDA2, not included in the paper.

