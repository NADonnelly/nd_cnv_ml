# ND-CNV Machine Learning

This project contains the code used in the ND-CNV machine learning paper (link here)

The data will be stored elsewhere (or here? or on OSF?)

The Main Scripts folder contains the R scripts that do the analysis:

* 1-Data-Preparation.R contains code that takes the raw data from a master spreadsheet and applies a process of data cleaning: only numeric data are selected, missing data codes are harmonised, variables and participants with > 25% missing data are removed (in that order) and the cleaned raw data are saved.

* 2-Descriptives.R contains code that makes the demographic details table

* 3-PLS.R performs principal components analysis and (sparse) Partial Least Squares Discriminant Analysis

* 4-ML.R fits machine learning models to the full set of variables using nested cross validation, determines variable importance with permutation testing, selects the most important variables from the models and re-fits ML models with the reduce variable sets, then fits models with the final sets of variables and final model hyperparameters to the held-out test data

* 5-Variable-Dimensions.R uses exploratory graph analysis to investigate the underlying dimensional structure of the variables selected to be most important to ND-CNV classification.

* 6-ML-Final.R takes the minimal set of variables from the dimensions identified by the EGA analysis and fits the ML classification models using only these variables, and measures performance