
#Introduction ========

#Here we are going to try and apply machine learning approaches to our dataset
#
#Our goal will be twofold:
#
# - To find a model that has the best possible performance
#
# - To compare models fit on different sets of variables for their predictive accuracy
#
# We are going to do a split on our dataset into test and training sets
#
# Then we will build some models using cross validation
#
# Apparently nested CV is the best approach?
#
# It is also necessary to do imputation in each CV fold to get best performance I believe?
#
# Which ML models shall we use?
#
# - Plain old logistic regression
# - Regularised regression with elastic net (glmnet)
# - Linear SVM
# - Something boosty e.g. XGBoost
# - ANN?
# - Ensembles/stacking models?
#
# Target metric for optimisation? I think we will use the AUC as the target


# Load packages ------

library(tidyverse)
library(tidymodels)
library(workflowsets)
library(tidyposterior)
library(finetune)
library(probably)
library(stacks)
library(doFuture)
library(baguette)

tidymodels_prefer()

# Load data ======

DF = read_csv("FilteredData2.csv")
# DF = read_csv("ImputedData2.csv")

DF0 = 
  DF |>
  mutate(group = map_chr(GenotypeCode,~ifelse(.x == 16,"Control","ND-CNV"))) |>
  mutate(group = factor(group, levels = c("Control","ND-CNV"))) |>
  select(-c(IDs,GenotypeCode)) |>
  relocate(group)

#Have as quick look at all variables
# skimr::skim(DF0)


# Prepare variable sets ------- 

vars_glm       = read_csv("selected_variables_glm.csv")
vals_LR_auc    = read_csv("selected_variables_LR_auc.csv")
vars_lasso     = read_csv("selected_variables_lasso.csv")
vars_low_corr  = read_csv("selected_variables_low_corr.csv")
vars_horseshoe = read_csv("selected_variables_horseshoe.csv")
vars_pp_11     = read_csv("selected_variables_projpred_11.csv")
vars_pp_25     = read_csv("selected_variables_projpred_25.csv")


# Compare our variables?
var_set = 
  bind_rows(vars_glm       |> select(variable) |> mutate(type = "glm"),
            vals_LR_auc    |> select(variable,type),
            vars_lasso     |> select(variable = Variable,type),
            vars_low_corr  |> select(variable,type),
            vars_horseshoe |> select(variable,type),
            vars_pp_11     |> select(variable,type),
            vars_pp_25     |> select(variable,type)) |>
  count(variable)


#Make a plot

var_set   |>
  ggplot(aes(x = fct_reorder(variable,n,max), y = n)) +
  geom_col() +
  theme_bw() +
  theme(axis.text.x.bottom = element_text(size = 4, angle = -90))


#Take only variables included in more than 1 method for extraction
vars_more_than_one = 
  var_set |>
  filter(n > 1)


#take the vars with the highest number of inclusions
vars_max = 
  var_set |>
  slice_max(order_by = n)

#Now we have our sets of variables, lets get some data ready...


# Splits and folds =====

# Lets force our ordinal variables to be integers (necessary for imputation)

## work out if a variable is all integers
is_int_var = 
  DF0 |> 
  pivot_longer(-c(group)) |> 
  group_by(name) |> 
  nest() |>
  mutate(all_int = map_lgl(data, ~ .x |>
                             drop_na() |> 
                             mutate(frac = map_lgl(value,~.x%%1 == 0)) |> 
                             summarise(is_int = all(frac)) |>
                             pull())) |>
  select(name,all_int) 

DF1 = 
  DF0 |>
  mutate(across(is_int_var |> filter(all_int == T) |> pull(name), as.integer)) 


set.seed(26052022)

# Prepare our intial data splits
d_split = initial_split(DF1, strata = group,prop = 3/4)

#Get test and training data
d_train = training(d_split)
d_test  = testing(d_split)

#Prepare for 10 fold cross validation
d_folds = vfold_cv(d_train, strata = group, v = 10, repeats = 1)

# Fit same model to the different sets of variables =======

#Lets fit the same model to our difference variable sets

## Define our pre-processing recipe =====

#This should include our imputation
#a reasonable place to start https://recipes.tidymodels.org/reference/step_impute_bag.html 

#Make a list of formulas

form_list = 
  bind_rows(vars_glm            |> select(variable) |> mutate(type = "glm"),
            vals_LR_auc         |> select(variable,type),
            vars_lasso          |> select(variable = Variable,type),
            vars_low_corr       |> select(variable,type),
            vars_horseshoe      |> select(variable,type),
            vars_pp_11          |> select(variable,type),
            vars_pp_25          |> select(variable,type),
            vars_max            |> select(variable) |> mutate(type = "max_vars"),
            vars_more_than_one  |> select(variable) |> mutate(type = "more_than_one")) |>
  group_by(type) |>
  nest() |>
  
  #Make our formulas
  mutate(formulas = map_chr(data,~paste("group ~",.x |> 
                                          pull(variable) |> 
                                          paste(collapse = " + ")))) |>
  
  
  #Make recipes which we will append our formulas into
  mutate(recipes = map(formulas, ~  recipe(formula = as.formula(.x),
                                           data = d_train) %>%
                         # step_normalize(!!!is_int_var |> filter(all_int == F) |> pull(name)) %>%
                         step_impute_bag(all_predictors())))


rm(list = ls(pattern = "vars_"))

## Define a model ----

#Define a model - Elastic Net Regression with glmnet
logistic_reg_spec <- 
  logistic_reg(mode = "classification",
               penalty = tune(), mixture = tune()) %>% 
  set_engine("glmnet")


## Make a workflow ===== 

#Make the workflow set
wf <- 
  workflow_set(
    preproc = form_list$recipes, 
    models = list(logistic_reg = logistic_reg_spec)
  )

wf = 
  wf |>
  mutate(var_type = form_list$type) |>
  mutate(wflow_id = map2_chr(wflow_id,var_type, ~str_replace(.x,"recipe_[0-9]",.y))) |>
  mutate(wflow_id = map_chr(wflow_id,str_remove,"_logistic_reg")) |>
  dplyr::select(-var_type)

## Fit models =====

race_ctrl <-
  control_race(
    save_pred = TRUE,
    parallel_over = "resamples",
    save_workflow = TRUE
  )


library(doParallel)

# Create a cluster object and then register: 
cl <- makePSOCKcluster(14)
registerDoParallel(cl)


#Fit our workflows
grid_results <-
  wf %>%
  workflow_map(
    "tune_race_anova",
    seed = 1503,
    resamples = d_folds,
    control = race_ctrl,
    grid = 25,
    metrics = metric_set(roc_auc)
  )


stopCluster(cl)

#Look at the results
grid_results |>
  rank_results(rank_metric = "roc_auc") %>% 
  filter(.metric == "roc_auc")

#Plot
autoplot(grid_results , metric = "roc_auc")


#We can directly compare with tidyposterior, which takes the best submodel from
#each workflow and then does a stan_glm on the performance metrics
roc_mod <- 
  perf_mod(grid_results, 
           metric = "roc_auc", 
           seed = 1, refresh = 0,
           chains = 4, cores = 8,
           prior = rstanarm::normal(0,1),
           prior_intercept = rstanarm::normal(0,1.5),
           prior_aux = rstanarm::exponential(rate = 1))

autoplot(roc_mod)
autoplot(roc_mod, type = "ROPE", size = 0.025)

#Looks like recipe 5 and recipe 9 are best. These are the horseshoe and the "more than
#one" variables - in general models with more variables do well. The variable sets with 
#fewest included variables do less well; although we are still getting AUCs of > 0.92
#which is very impressive

best_set = 
  roc_mod |>
  tidy() |> 
  group_by(model) |> 
  tidybayes::median_hdi() |> 
  arrange(-posterior)
  slice_max(posterior) |>
  pull(model)


# Fit multiple models =====

#I think we should use the horseshoe variables as these did absolutely best, and then
#also the  max_vars because they seem to have good performance and theres only 7 of them

#Get our recipes that did best
final_recipes = 
    form_list |>
    filter(type %in% c("Horseshoe","max_vars")) |>
    pull(recipes) 
  
#Name the recipes  
names(final_recipes) <-  c("Horseshoe","max_vars") 
  

## Define models we intend to use =======

#Set up a bunch of potential models

#Elastic Net Regression with glmnet
logistic_reg_spec <- 
  logistic_reg(mode = "classification",
               penalty = tune(), mixture = tune()) %>% 
  set_engine("glmnet")

# Radial SVM
svm_r_spec <- 
  svm_rbf(cost = tune(), rbf_sigma = tune()) %>% 
  set_engine("kernlab") %>% 
  set_mode("classification")

#Linear SVM
svm_l_spec = 
  svm_linear(cost = tune()) %>%
  set_engine("kernlab") %>% 
  set_mode("classification")

#Polynomial SVM
svm_p_spec <- 
  svm_poly(cost = tune(), degree = tune()) %>% 
  set_engine("kernlab") %>% 
  set_mode("classification")

#Neural Network with nnet
nnet_spec <- 
  mlp(hidden_units = tune(), penalty = tune(), epochs = tune()) %>% 
  set_engine("nnet", MaxNWts = 2600) %>% 
  set_mode("classification")

#Neural network with keras - this is very slow and needs tensorflow, which really needs a CUDA GPU
#So we will give in a miss
# nnet_ker_spec <-
#   mlp(epochs = tune(), hidden_units = tune(), dropout = tune()) %>%
#   set_mode("classification") %>%
#   # Also set engine-specific `verbose` argument to prevent logging the results:
#   set_engine("keras", verbose = 0)

#MARS
mars_spec <- 
  mars(prod_degree = tune()) %>%  #<- use GCV to choose terms
  set_engine("earth") %>% 
  set_mode("classification")

#KNN
knn_spec <- 
  nearest_neighbor(neighbors = tune(), dist_power = tune(), weight_func = tune()) %>% 
  set_engine("kknn") %>% 
  set_mode("classification")

#CART
cart_spec <- 
  decision_tree(cost_complexity = tune(), min_n = tune()) %>% 
  set_engine("rpart") %>% 
  set_mode("classification")

#Bagging CART? This appears to have problems, possibly with parallel processing. We can try and repeat it elsewhere
bag_cart_spec <-
  bag_tree(mode = "classification",
           min_n = tune()) %>%
  set_engine("C5.0")

#Random Forest
rf_spec <- 
  rand_forest(mtry = tune(), min_n = tune(), trees = 1000) %>% 
  set_engine("ranger") %>% 
  set_mode("classification")

#Boosted Trees
xgb_spec <- 
  parsnip::boost_tree(
    mode = "classification",
    trees = 1000,
    min_n = tune(),
    tree_depth = tune(),
    learn_rate = tune(),
    loss_reduction = tune() ) %>%
  set_engine("xgboost", objective = "binary:logistic")

#lightGBM?


## Make a workflowset ======

wf_par <- 
  workflow_set(
    preproc = final_recipes, 
    models  = list(logistic_reg = logistic_reg_spec,
                   SVM_radial   = svm_r_spec, 
                   SVM_linear   = svm_l_spec,
                   SWM_poly     = svm_p_spec,
                   nnet_n       = nnet_spec,
                   mars         = mars_spec,
                   knn          = knn_spec,
                   cart         = cart_spec,
                   rand_forest  = rf_spec)
  )

#Some models will run in parallel happily, some produce problems

# xgboost seems to be unhappy
# the keras models are unhappy (they really want to run on a GPU)

#These have to be fit not in parallel for some reason
wf_seq <- 
  workflow_set(
    preproc = final_recipes, 
    models  = list(nnet_k       = nnet_ker_spec,     
                   bagged_tree  = bag_cart_spec,
                   xgb_boost    = xgb_spec))
                   

#You can look at the full recipe for any of these e.g.
# wf |> extract_workflow(id = "Horseshoe_nnet_n")


## Do some model fitting
race_ctrl <-
  control_race(
    save_pred = TRUE,
    parallel_over = "everything",
    save_workflow = TRUE
  )


#Set up to run these in parallel
registerDoFuture()
plan(multisession, workers = 16)

#Fit our workflows that work in parallel
wf_par_results <-
  wf_par |>
  workflow_map(
    "tune_race_anova",
    seed = 1503,
    resamples = d_folds,
    control = race_ctrl,
    grid = 25,
    metrics = metric_set(roc_auc)
  )

wf_par_results = 
  wf_par_results |> 
  filter(!str_detect(wflow_id,"xgb"))

#Turn off parallel mode
plan(sequential)


#Fit the workflows that seem to work best sequentially
wf_seq_results <-
  wf_seq |>
  workflow_map(
    "tune_race_anova",
    seed = 1503,
    resamples = d_folds,
    control = race_ctrl,
    grid = 25,
    metrics = metric_set(roc_auc)
  )


#Unite the datasets
wf_results = 
  bind_rows(wf_par_results,wf_seq_results)

#Look at the results
wf_results |>
  rank_results(rank_metric = "roc_auc") %>% 
  filter(.metric == "roc_auc")

#Plot
autoplot(wf_results, metric = "roc_auc")

#Do our own
wf_results |>
  rank_results(rank_metric = "roc_auc") |>
  filter(.metric == "roc_auc") |>
  mutate(var_type = ifelse(str_detect(wflow_id,"Horseshoe"),"Horseshoe","Max Vars")) |>
  ggplot(aes(x = rank,y = mean, ymax = mean+std_err,ymin = mean-std_err,colour = model,shape = var_type)) +
  geom_errorbar() +
  geom_point(size = 2) +
  theme_bw()


#We can directly compare with tidyposterior, which takes the best submodel from
#each workflow and then does a stan_glm on the performance metrics
roc_mod <- 
  perf_mod(wf_results, 
           metric = "roc_auc", 
           seed = 1, refresh = 0,
           chains = 4, cores = 8,
           iter = 4000, warmup = 1000,
           prior = rstanarm::normal(0,1),
           prior_intercept = rstanarm::normal(0,1.5),
           prior_aux = rstanarm::exponential(rate = 1))


#Make our own plot
roc_mod_post = 
  roc_mod |>
  tidy() |>
  separate(model,into = c("vars","mod_type"),sep = "_", extra = "merge") |>
  mutate(mod_type = map_chr(mod_type,str_remove,pattern = "vars_")) |>
  group_by(vars,mod_type) |>
  tidybayes::median_hdi()

roc_mod_post |>
  ggplot(aes(x = forcats::fct_reorder(mod_type,posterior,max),y = posterior,ymin = .lower,ymax = .upper,group = vars)) +
  geom_errorbar(position = position_dodge(width = 0.5)) +
  geom_point(aes(shape = vars),position = position_dodge(width = 0.5)) +
  theme_bw() +
  theme(axis.text.x.bottom = element_text(size = 8, angle = -45)) +
  labs(y = "Model AUC", x = "Model Type",shape = "Variable \nSet")


#Do the ROPE Plot
autoplot(roc_mod, type = "ROPE", size = 0.025)



#they all look kind of similar apart from two which did poorly

#The best model appears to be the linear SVM, so we will go with that
best_results = 
  wf_results |>
  extract_workflow_set_result("Horseshoe_SVM_linear") |>
  select_best(metric = "roc_auc")

#Finalise
best_test_results <- 
  wf_results %>% 
  extract_workflow("Horseshoe_SVM_linear") %>% 
  finalize_workflow(best_results) %>% 
  last_fit(split = d_split)

#Get our final AUC
collect_metrics(best_test_results)

#Get the metrics for the  best version of every model?
best_results_wf = 
  tibble(wflow_id = wf_results$wflow_id) |>
  mutate(res = map(wflow_id,~extract_workflow_set_result(wf_results,id = .x)))  |>
  mutate(best_metric = map(res, select_best,metric = "roc_auc")) |>
  mutate(best_wf  = map(wflow_id, ~extract_workflow(wf_results,id = .x))) |>
  mutate(final_fit = map2(best_wf,best_metric,~finalize_workflow(.x,.y) |> 
                               last_fit(split = d_split)),
         final_metrics = map(final_fit,collect_metrics))


best_results_wf |>
  unnest(final_metrics) |>
  filter(.metric == "roc_auc") |>
  select(wflow_id,.metric,.estimate,.config) |>
  arrange(-.estimate) |>
  knitr::kable()


# Variable importance ------


library(DALEXtra)

set.seed(1803)

#Dalex does not seem to work well with missing data so we will bake our trianing data 
#again using the imputation within each workflow

explainer_svm_poly = 
  explain_tidymodels(best_results_wf$final_fit[[4]] |> 
                     extract_workflow(),
                     data = best_results_wf$final_fit[[4]] |> 
                       extract_recipe() |>
                       bake(d_train) |> 
                       select(-group) |>
                       mutate(across(where(is.integer),as.double)),
                     y = d_train$group |> as.numeric(),
                     label = "svm_poly",
                     verbose = FALSE) 

explainer_svm_l = 
  explain_tidymodels(best_results_wf$final_fit[[3]] |> 
                       extract_workflow(),
                     data = best_results_wf$final_fit[[3]] |> 
                       extract_recipe() |>
                       bake(d_train) |> 
                       select(-group) |>
                       mutate(across(where(is.integer),as.double)),
                     y = d_train$group |> as.numeric(),
                     label = "svm_l",
                     verbose = FALSE) 


explainer_rf = 
  explain_tidymodels(best_results_wf$final_fit[[18]] |> 
                       extract_workflow(),
                     data = best_results_wf$final_fit[[18]] |> 
                       extract_recipe() |>
                       bake(d_train) |> 
                       select(-group) |>
                       mutate(across(where(is.integer),as.double)),
                     y = d_train$group |> as.numeric(),
                     label = "rf",
                     verbose = FALSE) 


#Lets do one with the fewer components



vip_svm_poly <- 
  model_parts(explainer_svm_poly)


vip_svm_l <- 
  model_parts(explainer_svm_l)


vip_rf <- 
  model_parts(explainer_rf)




ggplot_imp <- function(...) {
  
  obj <- list(...)
  
  metric_name <- attr(obj[[1]], "loss_name")
  metric_lab <- paste(metric_name, 
                      "after permutations\n(higher indicates more important)")
  
  full_vip <- bind_rows(obj) %>%
    filter(variable != "_baseline_")
  
  perm_vals <- full_vip %>% 
    filter(variable == "_full_model_") %>% 
    group_by(label) %>% 
    summarise(dropout_loss = mean(dropout_loss))
  
  p <- full_vip %>%
    filter(variable != "_full_model_") %>% 
    mutate(variable = fct_reorder(variable, dropout_loss)) %>%
    ggplot(aes(dropout_loss, variable)) 
  
  if(length(obj) > 1) {
    p <- p + 
      facet_wrap(vars(label),ncol = 1) +
      geom_vline(data = perm_vals, aes(xintercept = dropout_loss, color = label),
                 size = 1.4, lty = 2, alpha = 0.7) +
      geom_boxplot(aes(color = label, fill = label), alpha = 0.2)
  } else {
    p <- p + 
      geom_vline(data = perm_vals, aes(xintercept = dropout_loss),
                 size = 1.4, lty = 2, alpha = 0.7) +
      geom_boxplot(fill = "#91CBD765", alpha = 0.4)
    
  }
  p +
    theme(legend.position = "none") +
    labs(x = metric_lab, 
         y = NULL,  fill = NULL,  color = NULL)
}


#Plot the variable importance thing
ggplot_imp(vip_svm_poly, vip_svm_l,vip_rf)


#They look similar enough.

# Ensembles =======

#Lets do an ensemble for the set of 7 variables because we would most like to increase the performance of models with only those variables-

library(stacks)

d_stack = 
  stacks() |>
  add_candidates(wf_results |> 
                   filter(str_detect(wflow_id,"max")) 
  ) 

set.seed(27052022)

d_stack_blend = 
  blend_predictions(d_stack,
                    metric = metric_set(roc_auc))

d_stack_fit = 
  d_stack_blend |>
  fit_members()


# Predict on new data

member_preds <- 
  d_test %>%
  select(group) %>%
  bind_cols(
    predict(
      d_stack_fit, 
      d_test, 
      members = TRUE,
      type = "prob"
    )
  )

member_metrics = 
  member_preds |> 
  pivot_longer(-group) |>
  filter(!str_detect(name,"CNV")) |>
  group_by(name) |>
  nest() 

member_metrics |>
  mutate(model_auc = map(data,roc_auc, truth = group,value)) |>
  unnest(model_auc)

#So the ensemble prediction is... worse than individual members?
#Suggests the ensemble is overfitting?


# Thresholds =====


#Make a plot of classification performance at difference thresholds


-


best_test_results %>% 
  collect_predictions() |>
  ggplot(aes(`.pred_ND-CNV`,colour = group,fill = group)) +
  geom_density(position = "identity",alpha = 0.5)


  best_test_results %>% 
  collect_predictions() %>% 
  probably::threshold_perf(group,`.pred_ND-CNV`,
                           thresholds = 0.5)


threshold_data = 
  best_test_results %>% 
  collect_predictions() %>% 
  probably::threshold_perf(group,.pred_Control,
                           thresholds = seq(0, 1, by = 0.0025))


threshold_data <- 
  threshold_data %>%
  filter(.metric != "distance") %>%
  mutate(group = case_when(
    .metric == "sens" | .metric == "spec" ~ "1",
    TRUE ~ "2"
  ))

max_j_index_threshold <- 
  threshold_data %>%
  filter(.metric == "j_index") %>%
  filter(.estimate == max(.estimate)) %>%
  pull(.threshold) |>
  median()

ggplot(threshold_data, aes(x = -.threshold, y = .estimate, color = .metric)) +
  geom_line(size = 1.25) +
  theme_bw() +
  # scale_color_viridis_d(end = 0.9) +
  scale_alpha_manual(values = c(.4, 1), guide = "none") +
  geom_vline(xintercept = -max_j_index_threshold, alpha = .6, color = "grey30",
             size = 1.25, lty = 2) +
  labs(
    x = "Threshold",
    y = "Metric Estimate",
    title = "Performance by threshold",
  ) +
  theme(panel.grid = element_blank())




# Prepare data =======

# Fit models to the variables which had the highest coefficients in the stan_glm model
#with the horseshoe/heirarchical shrinkage prior
d_hs = 
  DF0 |>
  select(group, all_of(vars_horseshoe$variable))

#Split the data into training and testing sets
split_hs = initial_split(d_hs,prop = 4/5)
train_hs = training(split_hs)
test_hs  = testing(split_hs)








# Nested resampling =======

#Make data for nested resampling
results <- nested_cv(train_hs, 
                     outside = vfold_cv(v = 5, repeats = 1), 
                     inside = bootstraps(times = 25))

# results$inner_resamples[[2]]


#Begin building up a recipe
library(kernlab)

# `object` will be an `rsplit` object from our `results` tibble
# `cost` is the tuning parameter
svm_auc <- function(object, cost = 1) {

  mod_rec = 
    object |>
    analysis() |>
    recipe(group ~ .) |>
    step_impute_bag(all_predictors()) 
  
  mod_spec <- 
    svm_rbf(mode = "classification", cost = 0.25,
            engine = "kernlab") 
  
  mod_wf = 
    workflow() |>
    add_recipe(mod_rec)
  
  mod_fit = 
    mod_wf %>%
    add_model(mod_spec) %>%
    fit(object |> analysis() )
  
  
  mod_auc <- 
    predict(mod_fit, assessment(object) |> dplyr::select(-group), type = "prob") |> 
    bind_cols(assessment(object) |> dplyr::select(group)) |>
    roc_auc(truth = group, estimate = .pred_Control)
  
  return(mod_auc$.estimate)
  
}

# In some case, we want to parameterize the function over the tuning parameter:
auc_wrapper <- function(cost, object) svm_auc(object, cost)


# `object` will be an `rsplit` object for the bootstrap samples
tune_over_cost <- function(object) {
  tibble(cost = 2 ^ seq(-2, 8, by = 1)) %>% 
    mutate(AUC = map_dbl(cost, auc_wrapper, object = object))
}


# `object` is an `rsplit` object in `results$inner_resamples` 
summarize_tune_results <- function(object) {
  
  # Return row-bound tibble that has the 25 bootstrap results
  map_df(object$splits, tune_over_cost) %>%
    
    # For each value of the tuning parameter, compute the 
    # average AUC which is the inner bootstrap estimate. 
    group_by(cost) %>%
    summarize(mean_AUC = mean(AUC, na.rm = TRUE),
              n = length(AUC),
              .groups = "drop")
}


library(furrr)
plan(multisession)

tuning_results <- 
  future_map(results$inner_resamples, summarize_tune_results) 


#So this has used our bootstrapped samples to tune the cost parameter
#We can look at the results of this


library(scales)

pooled_inner <- 
  tuning_results %>% bind_rows

best_cost <- 
  function(dat) dat[which.max(dat$mean_AUC),]

p <- 
  ggplot(pooled_inner, aes(x = cost, y = mean_AUC)) + 
  scale_x_continuous(trans = 'log2') +
  xlab("SVM Cost") + ylab("Inner AUC")

for (i in 1:length(tuning_results))
  p <- p  +
  geom_line(data = tuning_results[[i]], alpha = .2) +
  geom_point(data = best_cost(tuning_results[[i]]), pch = 16, alpha = 3/4)

p <- p + geom_smooth(data = pooled_inner, se = FALSE)
p



# And then we use this value to fit models to our outer datasets
cost_vals <- 
  tuning_results %>% 
  map_df(best_cost) %>% 
  select(cost)

results <- 
  bind_cols(results, cost_vals) %>% 
  mutate(cost = factor(cost, levels = paste(2 ^ seq(-2, 8, by = 1))))

ggplot(results, aes(x = cost)) + 
  geom_bar() + 
  xlab("SVM Cost") + 
  scale_x_discrete(drop = FALSE)


#Now we use the cost parameter optimised on the inner loops to 
#calculate performance on each of our k folds
results <- 
  results %>% 
  mutate(AUC = map2_dbl(splits, cost, svm_auc))


summary(results$AUC)


not_nested <- 
  map(results$splits, tune_over_cost) %>%
  bind_rows


outer_summary <- not_nested %>% 
  group_by(cost) %>% 
  summarize(outer_AUC = mean(AUC), n = length(AUC))



ggplot(outer_summary, aes(x = cost, y = outer_AUC)) + 
  geom_point() + 
  geom_line() + 
  scale_x_continuous(trans = 'log2') +
  xlab("SVM Cost") + ylab("AUC")




