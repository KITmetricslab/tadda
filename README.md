# Code for the paper _Direction Augmentation in the Evaluation of Armed Conflict Predictions_

by Johannes Bracher, Lotta Rüter, Fabian Krüger, Sebastian Lerch and Melanie Schienle published in _International Interactions_ (forthcoming)

## Contents
Please note that we use the terms optimal point forecast (OPF) and bayes act (BA) synonymously in this code.

### **Empirical Example**
contains data and code to reproduce all empirical results
+ The file _bayes_acts_functions.R_ contains functions for computing the optimal point forecasts (OPFs) of different scoring functions (AE, SE and variants of TADDA) as well as a summary function for a nicer represenation of the results
+ The file _example_mali.R_ generates Figure 3
+ The file _tadda_example.R_ contains the code for empirical example in Section 5: it computes the optimal window size and generates Table 2 and Table 3

**Data**
+ The file _data_prep.R_
  * extracts time series of country month fatalities due to state based conflict "ged_best_sb" for each country in Africa
  * computes true s=1,...,7 step ahead log-changes for each country and month
  * saves the results in "fatalities.csv"
+ The files _ged_cm_postpatch.parquet_ and _skeleton_cm_africa.parquet_ were retrieved from https://github.com/UppsalaConflictDataProgram/views_competition/tree/main/data (published under Creative Commons Attribution-NonCommercial 4.0 International Public License).

**Results**


### **Simulation**
contains code to reproduce illustrative figures and small simulation examples


The contents of this repository are likewise under a Creative Commons Attribution-NonCommercial 4.0 International Public License.
