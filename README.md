# ND-GC Machine Learning Project

This project contains the code used in our ND-GC variable selection machine learning paper: [preprint](https://www.medrxiv.org/content/10.1101/2022.12.16.22283581v1)

The data used in the study is available via the [IMAGINE ID study](https://imagine-id.org/healthcare-professionals/)

The Main Scripts folder contains the R scripts that do the analysis:

* 1-Data-Preparation.R contains code that takes the raw data from a master spreadsheet and applies a process of data cleaning: only numeric data are selected, some variables are recoded, missing data codes are harmonised, variables with > 90% of responses the same are removed, variables and participants with > 25% missing data are removed (in that order) and highly correlated (>0.8) variables are removed. The cleaned raw data are saved.

* 2-Descriptives.R contains code that makes the demographic details table

* 3-Data-Split.R contains code that performs the initial split of the dataset into training (80%) and test (20%) data, stratified by group, gender and age

* 4-PLS.R performs principal components analysis and (sparse) Partial Least Squares Discriminant Analysis

* 5-ML-All-Variables.R fits machine learning models to the full set of variables using nested cross validation

* 6-Variable-Importance.R determines variable importance with permutation testing and selects the most important variables from the models 

* 7-ML-Selected-Variables.R re-fits ML models with the reduced variable sets

* 8-Model-Evaluation.R fits the best performing models with the final sets of variables and final model hyperparameters to the held-out test data

* 9-Variable-Dimensions.R uses exploratory graph analysis to investigate the underlying dimensional structure of the variables selected to be most important to ND-GC classification.

* 10-ML-Final.R takes the minimal set of variables from the dimensions identified by the EGA analysis and fits the ML classification models using only these variables, and measures performance

* 11-Model-Deployment.R contains code for making models for the accompanying [shiny app](https://nadonnelly.shinyapps.io/cnv_ml_app/)

The 
