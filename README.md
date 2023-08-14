# Code for the paper _Direction Augmentation in the Evaluation of Armed Conflict Predictions_

by Johannes Bracher, Lotta Rüter, Fabian Krüger, Sebastian Lerch and Melanie Schienle published in _International Interactions_ (forthcoming).

## Contents
### **Empirical Example**
contains data and code to reproduce all empirical results
**Data**
+ The file _data_prep.R_
  * extracts time series of country month fatalities due to state based conflict "ged_best_sb" for each country in Africa
  * Computes true s=1,...,7 step ahead log-changes for each country and month
  * Saves the results in "fatalities.csv"
+ The files _ged_cm_postpatch.parquet_ and _skeleton_cm_africa.parquet_ were retrieved from https://github.com/UppsalaConflictDataProgram/views_competition/tree/main/data (published under Creative Commons Attribution-NonCommercial 4.0 International Public License).


### **Simulation**
contains code to reproduce illustrative figures and small simulation examples.


The contents of this repository are likewise under a Creative Commons Attribution-NonCommercial 4.0 International Public License.
