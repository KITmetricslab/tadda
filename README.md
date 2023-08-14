# Code for the paper _Direction Augmentation in the Evaluation of Armed Conflict Predictions_

by Johannes Bracher, Lotta Rüter, Fabian Krüger, Sebastian Lerch and Melanie Schienle published in _International Interactions_ (forthcoming)

## Notes

+ The contents of this repository are under a Creative Commons Attribution-NonCommercial 4.0 International Public License.
+ Please note that we use the terms optimal point forecast (OPF) and Bayes act (BA) synonymously in this code.

## Contents

### **_Empirical Example_** – data and code to reproduce all empirical results
+ _tadda_example.R_ contains the code for empirical example in Section 5: computes the optimal window size and generates Table 2 and Table 3
+ _example_mali.R_ generates Figure 3
+ _bayes_acts_functions.R_ contains functions for computing the optimal point forecasts (OPFs) of different scoring functions (AE, SE and variants of TADDA) as well as a summary function for a nicer represenation of the results

**_Data_** – data for the empirical analysis
+ _data_prep.R_
  * extracts time series of country month fatalities due to state based conflict "ged_best_sb" for each country in Africa
  * computes true s = 1, ..., 7 step ahead log-changes for each country and month
+ _fatalities.csv_ contains results of _data_prep.R_
+ _ged_cm_postpatch.parquet_ and _skeleton_cm_africa.parquet_ were retrieved from https://github.com/UppsalaConflictDataProgram/views_competition/tree/main/data (published under Creative Commons Attribution-NonCommercial 4.0 International Public License).

**_Results_** – optimal window size, Table 2 and Table 3
+ _average_scores_for_different_window_lengths.csv_ shows that w=5 would be optimal for minimizing TADDA1 via TADDA1_OPF
+ _individual_predictions_w9.csv_ and _individual_losses_w9.csv_ contain the predictions for the log-changes in fatalities and corresponding losses for each African country, month in 395:495, lead time s = 2, ..., 7 and scoring function / OPF
+ _average_scores_w9.csv_ contains the central results presented in Table 2
+ _empirical_quantiles_w9.csv_ is the basis of Table 3

### **_Simulation_** – illustrative figures and small simulation examples
+ _functions.R_ contains functions to compute different scoring functions (AE, SE and variants of TADDA), to numerically determine the optimal point forecast (Bayes act) given a distribution and to provide text annotation in a plot
+ _illustration.R_ generates Figure 2, Table 1, and Supplementary Figure S6
+ _illustrations_proof.R_ generates Supplementary Figures S4 and S5
