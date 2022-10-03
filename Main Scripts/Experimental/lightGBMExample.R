
# Penguin Example ======

# Try and fit an example model with lightGBM and tidymodels, including
# variable importance

pacman::p_load(bonsai, tidymodels, modeldata, tidyverse, tidyposterior)

data(penguins)


#Make some data splits
d_split = initial_split(penguins, strata = species,prop = 4/5)

#Get test and training data
d_train = training(d_split)
d_test  = testing(d_split)
d_fold  = rsample::vfold_cv(d_train, v = 5, repeats = 2, strata = species)


p_rec = 
  d_train |>
  recipes::recipe(species ~ .) |>
  step_normalize(all_numeric_predictors()) |>
  step_zv(all_predictors()) 
  
  
bt_mod <- 
  boost_tree(trees = 1000,
             min_n = tune(),
             tree_depth = tune(),
             learn_rate = tune(),
             loss_reduction = tune()) %>%
  set_engine(engine = "lightgbm") %>%
  set_mode(mode = "classification") 

p_params <- 
  bt_mod |>
  extract_parameter_set_dials()

# p_params <- 
#   p_params|>
#   update(`mtry` = finalize(mtry(),d_train))

p_wf <- 
  workflow() |>
  add_recipe(p_rec) |>
  add_model(bt_mod)

p_grid <- 
  dials::grid_max_entropy(
    p_params, 
    size = 30 # set this to a higher number to get better results
    # I don't want to run this all night, so I set it to 30
  )


#Do grid search
lgbm_tuned <- 
  p_wf |>
  tune::tune_grid(
    resamples = d_fold,
    grid      = p_grid,
    metrics = yardstick::metric_set(roc_auc, accuracy, f_meas),
    control = tune::control_grid(verbose = TRUE) 
)

#Get our best set of tuning parameters
lgbm_tuned |>
  show_best(metric = "roc_auc")

#We can make a plot of performance by parameter
lgbm_tuned %>%  
  tune::show_best(metric = "roc_auc",n = 10) %>% 
  tidyr::pivot_longer(min_n:loss_reduction, names_to="variable",values_to="value" ) %>% 
  ggplot(aes(value,mean)) + 
  geom_line(alpha=1/2)+ 
  geom_point()+ 
  facet_wrap(~variable,scales = "free")+
  ggtitle("Best parameters for AUC")


#Get the best parameters
lgbm_best_params <- 
  lgbm_tuned %>%
  tune::select_best("roc_auc")

#Finalize our modelling workflow
p_wf_final <- 
  p_wf %>% 
  finalize_workflow(lgbm_best_params)

#Fit on the split data, testing on the previously unseen test data
p_final_fit =
  p_wf_final |>
  last_fit(d_split)

#Plot our predictions
p_final_fit %>%
  collect_predictions() %>% 
  roc_curve(species, .pred_Adelie:.pred_Gentoo) %>%
  autoplot()

#Tabulate performance
p_final_fit %>%
  collect_predictions() %>%
  roc_auc(species, .pred_Adelie:.pred_Gentoo)

#So the model can get perfect performance...

#And you can get the importance from this model too...
p_final_fit |>
  extract_fit_engine() |>
  lightgbm::lgb.importance()



#But can we do the same thing with our data? =====

# We will use non - nested CV as a simple toy example

## Prepare Data =====

`%nin%` = Negate(`%in%`)

#Load the full dataset
DF = read_csv("FilteredData2.csv")

#Make the class labels
DF0 = 
  DF |>
  mutate(group = map_chr(GenotypeCode,~ifelse(.x == 16,"Control","ND-CNV"))) |>
  mutate(group = factor(group, levels = c("Control","ND-CNV"))) |>
  select(-c(IDs,GenotypeCode)) |>
  relocate(group)

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

#Set a random number seed for reproducibility
set.seed(09062022)


#Set our data splitting ratios and number of folds
split_prop = 3/4

n_outer    = 10
n_inner    = 10
n_fold     = 5


# Prepare our initial data splits
d_split = initial_split(DF1, strata = group,prop = split_prop)

#Get test and training data
d_train = training(d_split)
d_test  = testing(d_split)

rec_impute =
  d_train |>
  recipe(group ~ .) |>
  step_zv(all_predictors()) |>
  step_impute_bag(all_predictors())   

#Prep a model for imputing missing data
rec_impute_prep = 
  rec_impute |>
  prep()

#Get our imputed data back
d_train = 
  rec_impute_prep |>
  juice()

#And make a test set using the same imputation model, but the test data
d_test = 
  rec_impute_prep |>
  bake(new_data = 
         d_test)


d_fold  = rsample::vfold_cv(d_train, v = 5, repeats = 2, strata = group)

#We could maybe rebuild a split object with the imputed data
d_split_impute = 
  make_splits(d_train,d_test)

# Set up model recipe =====
lgbm_rec = 
  d_train |>
  recipes::recipe(group ~ .) |>
  step_zv(all_predictors()) 


lgbm_mod <- 
  boost_tree(trees = 1000,
             min_n = tune(),
             tree_depth = tune(),
             learn_rate = tune(),
             loss_reduction = tune()) %>%
  set_engine(engine = "lightgbm") %>%
  set_mode(mode = "classification") 

lgbm_params <- 
  lgbm_mod |>
  extract_parameter_set_dials()

lgbm_wf <- 
  workflow() |>
  add_recipe(lgbm_rec) |>
  add_model(lgbm_mod)


## Do an initial grid search =====

lgbm_initial_grid <- 
  dials::grid_max_entropy(
    lgbm_params, 
    size = 10 #Low number because we are going to use tune_bayes later
  )


#Do grid search
lgbm_initial_tuned <- 
  lgbm_wf |>
  tune::tune_grid(
    resamples = d_fold,
    grid      = lgbm_initial_grid,
    metrics = yardstick::metric_set(roc_auc, accuracy, f_meas),
    control = tune::control_grid(verbose = TRUE) 
  )

lgbm_initial_tuned |>
  autoplot()


## Bayesian grid search =====

#Now do bayesian grid search to further optimise
set.seed(1403)

lgbm_bo <-
  lgbm_wf %>%
  tune_bayes(
    resamples = d_fold,
    metrics = yardstick::metric_set(roc_auc, accuracy, f_meas),
    initial = lgbm_initial_tuned,
    param_info = lgbm_params,
    iter = 25,
    control = control_bayes(verbose = TRUE)
  )

lgbm_bo |>
  autoplot(type = "performance")


#Get our best set of tuning parameters
lgbm_bo |>
  show_best(metric = "roc_auc")

#Compare with the un-tuned grid search
lgbm_initial_tuned|>
  show_best(metric = "roc_auc")

#We can make a plot of performance by parameter
lgbm_bo %>%  
  tune::show_best(metric = "roc_auc",n = 10) %>% 
  tidyr::pivot_longer(min_n:loss_reduction, names_to="variable",values_to="value" ) %>% 
  ggplot(aes(value,mean)) + 
  geom_line(alpha=1/2)+ 
  geom_point()+ 
  facet_wrap(~variable,scales = "free")+
  ggtitle("Best parameters for AUC")


## Finalise model and predict test data ====

#Get the best parameters
lgbm_best_params <- 
  lgbm_bo %>%
  tune::select_best("roc_auc")

#Finalize our modelling workflow
lgbm_wf_final <- 
  lgbm_wf %>% 
  finalize_workflow(lgbm_best_params)

#Fit on the split data, testing on the previously unseen test data
lgbm_final_fit =
  lgbm_wf_final |>
  last_fit(d_split)

#Plot our predictions
lgbm_final_fit %>%
  collect_predictions() %>% 
  roc_curve(group, .pred_Control) %>%
  autoplot()

#Tabulate performance
lgbm_final_fit %>%
  collect_predictions() %>%
  roc_auc(group, .pred_Control)

lgbm_final_fit %>%
  collect_predictions() %>%
  average_precision(group,.pred_Control )
         
#This model can get extremely good performance


#Make a plot of classification performance at difference thresholds


  
lgbm_final_fit %>% 
  collect_predictions() |>
  ggplot(aes(.pred_Control,colour = group,fill = group)) +
  geom_density(position = "identity",alpha = 0.5)



threshold_data = 
  lgbm_final_fit %>% 
  collect_predictions() %>% 
  probably::threshold_perf(group,.pred_Control,
                           thresholds = seq(0, 1, by = 0.0025)) %>%
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

ggplot(threshold_data, aes(x = .threshold, y = .estimate, color = .metric)) +
  geom_line(size = 1.25) +
  theme_bw() +
  # scale_color_viridis_d(end = 0.9) +
  scale_alpha_manual(values = c(.4, 1), guide = "none") +
  geom_vline(xintercept = max_j_index_threshold, alpha = .6, color = "grey30",
             size = 1.25, lty = 2) +
  labs(
    x = "Threshold",
    y = "Metric Estimate",
    title = "Performance by threshold"  ) +
  theme(panel.grid = element_blank())

## Variable importance ====

#And you can get the importance from this model too...
lgbm_final_fit |>
  extract_fit_engine() |>
  lightgbm::lgb.importance() |>
  as_tibble() |>
  ggplot(aes(y = forcats::fct_reorder(Feature,Gain,max), x = Gain)) +
  geom_col() 


lgbm_vars = 
  lgbm_final_fit |>
  extract_fit_engine() |>
  lightgbm::lgb.importance() |>
  as_tibble() |>
  slice_max(Gain,n = 30) 

lgbm_vars |>
  knitr::kable()



# Compare with another model ======

#Lets find out if we can get similar performance with a different model
#and also use this as a chance to compare variable importance. We are going
#to use a linear SVM here

#Linear SVM
svm_l_mod = 
  svm_linear(cost = tune(),margin = tune()) %>%
  set_engine("kernlab") %>% 
  set_mode("classification")


svm_l_params <- 
  svm_l_mod |>
  extract_parameter_set_dials()

svm_l_wf <- 
  workflow() |>
  add_recipe(lgbm_rec) |>
  add_model(svm_l_mod)


## Do an initial grid search =====

svm_l_initial_grid <- 
  dials::grid_max_entropy(
    svm_l_params, 
    size = 10 #Low number because we are going to use tune_bayes later
  )


#Do grid search
svm_l_initial_tuned <- 
  svm_l_wf |>
  tune::tune_grid(
    resamples = d_fold,
    grid      = svm_l_initial_grid,
    metrics = yardstick::metric_set(roc_auc, accuracy, f_meas),
    control = tune::control_grid(verbose = TRUE) 
  )

svm_l_initial_tuned |>
  autoplot()


## Bayesian grid search =====

#Now do bayesian grid search to further optimise
set.seed(1403)

svm_l_bo <-
  svm_l_wf %>%
  tune_bayes(
    resamples = d_fold,
    metrics = yardstick::metric_set(roc_auc, accuracy, f_meas),
    initial = svm_l_initial_tuned,
    param_info = svm_l_params,
    iter = 25,
    control = control_bayes(verbose = TRUE)
  )

svm_l_bo |>
  autoplot(type = "performance")


#Get our best set of tuning parameters
svm_l_bo |>
  show_best(metric = "roc_auc")

#Compare with the un-tuned grid search
svm_l_initial_tuned|>
  show_best(metric = "roc_auc")

#Pretty similar


#We can make a plot of performance by parameter
svm_l_bo %>%  
  tune::show_best(metric = "roc_auc",n = 10) %>% 
  tidyr::pivot_longer(cost:margin, names_to="variable",values_to="value" ) %>% 
  ggplot(aes(value,mean)) + 
  geom_line(alpha=1/2)+ 
  geom_point()+ 
  facet_wrap(~variable,scales = "free")+
  ggtitle("Best parameters for AUC")


## Finalise model and predict test data ====

#Get the best parameters
svm_l_best_params <- 
  svm_l_bo %>%
  tune::select_best("roc_auc")

#Finalize our modelling workflow
svm_l_wf_final <- 
  svm_l_wf %>% 
  finalize_workflow(svm_l_best_params)

#Fit on the split data, testing on the previously unseen test data
svm_l_final_fit =
  svm_l_wf_final |>
  fit(data = d_train)

#Plot our predictions
d_test |>
  bind_cols(
    predict(svm_l_final_fit, d_test, type = "prob")
  ) %>% 
  roc_curve(group, .pred_Control) %>%
  autoplot()

#Tabulate performance
d_test |>
  bind_cols(
    predict(svm_l_final_fit, d_test, type = "prob")
  ) %>%
  roc_auc(group, .pred_Control)

d_test |>
  bind_cols(
    predict(svm_l_final_fit, d_test, type = "prob")
  )%>%
  average_precision(group,.pred_Control )

#This model can also get extremely good performance - perhaps even better than the lightgbm


#Make a plot of classification performance at difference thresholds



threshold_data = 
  d_test |>
  bind_cols(
    predict(svm_l_final_fit, d_test, type = "prob")
  ) %>% 
  probably::threshold_perf(group,.pred_Control,
                           thresholds = seq(0, 1, by = 0.0025)) %>%
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

ggplot(threshold_data, aes(x = .threshold, y = .estimate, color = .metric)) +
  geom_line(size = 1.25) +
  theme_bw() +
  # scale_color_viridis_d(end = 0.9) +
  scale_alpha_manual(values = c(.4, 1), guide = "none") +
  geom_vline(xintercept = max_j_index_threshold, alpha = .6, color = "grey30",
             size = 1.25, lty = 2) +
  labs(
    x = "Threshold",
    y = "Metric Estimate",
    title = "Performance by threshold"  ) +
  theme(panel.grid = element_blank())



# Lets do another model - lets say nnet =====

nnet_mod <- 
  mlp(hidden_units = tune(), penalty = tune(), epochs = tune()) %>% 
  set_engine("nnet", MaxNWts = 2600) %>% 
  set_mode("classification")

nnet_params <- 
  nnet_mod |>
  extract_parameter_set_dials()

nnet_wf <- 
  workflow() |>
  add_recipe(lgbm_rec) |>
  add_model(nnet_mod)


## Do an initial grid search =====

nnet_initial_grid <- 
  dials::grid_max_entropy(
    nnet_params, 
    size = 10 #Low number because we are going to use tune_bayes later
  )


#Do grid search
nnet_initial_tuned <- 
  nnet_wf |>
  tune::tune_grid(
    resamples = d_fold,
    grid      = nnet_initial_grid,
    metrics = yardstick::metric_set(roc_auc, accuracy, f_meas),
    control = tune::control_grid(verbose = TRUE) 
  )

nnet_initial_tuned |>
  autoplot()


## Bayesian grid search =====

#Now do bayesian grid search to further optimise
set.seed(1403)

nnet_bo <-
  nnet_wf %>%
  tune_bayes(
    resamples = d_fold,
    metrics = yardstick::metric_set(roc_auc, accuracy, f_meas),
    initial = nnet_initial_tuned,
    param_info = nnet_params,
    iter = 25,
    control = control_bayes(verbose = TRUE)
  )

nnet_bo |>
  autoplot(type = "performance")


#Get our best set of tuning parameters
nnet_bo |>
  show_best(metric = "roc_auc")

#Compare with the un-tuned grid search
nnet_initial_tuned|>
  show_best(metric = "roc_auc")

#Pretty similar


#We can make a plot of performance by parameter
nnet_bo %>%  
  tune::show_best(metric = "roc_auc",n = 10) %>% 
  tidyr::pivot_longer(hidden_units:epochs, names_to="variable",values_to="value" ) %>% 
  ggplot(aes(value,mean)) + 
  geom_line(alpha=1/2)+ 
  geom_point()+ 
  facet_wrap(~variable,scales = "free")+
  ggtitle("Best parameters for AUC")


## Finalise model and predict test data ====

#Get the best parameters
nnet_best_params <- 
  nnet_bo %>%
  tune::select_best("roc_auc")

#Finalize our modelling workflow
nnet_wf_final <- 
  nnet_wf %>% 
  finalize_workflow(nnet_best_params)

#Fit on the split data, testing on the previously unseen test data
nnet_final_fit =
  nnet_wf_final |>
  fit(data = d_train)

#Plot our predictions
d_test |>
  bind_cols(
    predict(nnet_final_fit, d_test, type = "prob")
  ) %>% 
  roc_curve(group, .pred_Control) %>%
  autoplot()

#Tabulate performance
d_test |>
  bind_cols(
    predict(nnet_final_fit, d_test, type = "prob")
  ) %>%
  roc_auc(group, .pred_Control)

d_test |>
  bind_cols(
    predict(nnet_final_fit, d_test, type = "prob")
  )%>%
  average_precision(group,.pred_Control )


#This model can also get extremely good performance - perhaps even better than the lightgbm


#Make a plot of classification performance at difference thresholds

threshold_data = 
  d_test |>
  bind_cols(
    predict(nnet_final_fit, d_test, type = "prob")
  ) %>% 
  probably::threshold_perf(group,.pred_Control,
                           thresholds = seq(0, 1, by = 0.0025)) %>%
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

ggplot(threshold_data, aes(x = .threshold, y = .estimate, color = .metric)) +
  geom_line(size = 1.25) +
  theme_bw() +
  # scale_color_viridis_d(end = 0.9) +
  scale_alpha_manual(values = c(.4, 1), guide = "none") +
  geom_vline(xintercept = max_j_index_threshold, alpha = .6, color = "grey30",
             size = 1.25, lty = 2) +
  labs(
    x = "Threshold",
    y = "Metric Estimate",
    title = "Performance by threshold"  ) +
  theme(panel.grid = element_blank())


## Variable importance

#And you can get the importance from this model too...
# nnet_vars =
#   nnet_final_fit |>
#   extract_fit_engine() |>
#   vip::vi_permute(train = d_train,
#                   target = "group",
#                   metric = "auc",
#                   pred_wrapper = predict,
#                   reference_class = "Control")
#   
# 
# nnet_vars |> 
#   arrange(-Importance)|>
#   ggplot(aes(y = forcats::fct_reorder(Variable,Importance,max), x = Importance)) +
#   geom_col() 
# 
# 
# nnet_vars |> 
#   slice_max(Importance,n = 30) |>
#   knitr::kable()


nnet_explainer = 
  nnet_final_fit |>
  extract_fit_parsnip() |>
  explain_tidymodels(data  = d_train |> select(-group),
                     y     = d_train |> 
                       select(group) |> 
                       mutate(group = ifelse(group == "Control",1,0)) |>
                       pull()
                     ) 

nnet_vars =  
  nnet_explainer|> 
  model_parts(type = "variable_importance")

nnet_vars |> 
  plot()



nnet_vars_30 = 
  nnet_vars |>
  as_tibble() |> 
  group_by(variable) |>
  summarise(imp_mu = mean(dropout_loss),
            imp_sd = sd(dropout_loss)) |>
  filter(variable %nin% c("_baseline_","_full_model_")) |>
  slice_max(imp_mu,n = 30)

nnet_vars_30 |>
  knitr::kable()


# Finish with elastic net ======

en_mod <- 
  logistic_reg(mode = "classification",
               penalty = tune(), mixture = tune()) %>% 
  set_engine("glmnet")


en_params <- 
  en_mod |>
  extract_parameter_set_dials()

en_wf <- 
  workflow() |>
  add_recipe(lgbm_rec) |>
  add_model(en_mod)


## Do an initial grid search =====

en_initial_grid <- 
  dials::grid_max_entropy(
    en_params, 
    size = 10 #Low number because we are going to use tune_bayes later
  )


#Do grid search
en_initial_tuned <- 
  en_wf |>
  tune::tune_grid(
    resamples = d_fold,
    grid      = en_initial_grid,
    metrics = yardstick::metric_set(roc_auc, accuracy, f_meas),
    control = tune::control_grid(verbose = TRUE) 
  )

en_initial_tuned |>
  autoplot()


## Bayesian grid search =====

#Now do bayesian grid search to further optimise
set.seed(1403)

en_bo <-
  en_wf %>%
  tune_bayes(
    resamples = d_fold,
    metrics = yardstick::metric_set(roc_auc, accuracy, f_meas),
    initial = en_initial_tuned,
    param_info = en_params,
    iter = 25,
    control = control_bayes(verbose = TRUE)
  )

en_bo |>
  autoplot(type = "performance")


#Get our best set of tuning parameters
en_bo |>
  show_best(metric = "roc_auc")

#Compare with the un-tuned grid search
en_initial_tuned|>
  show_best(metric = "roc_auc")

#Pretty similar


#We can make a plot of performance by parameter
en_bo %>%  
  tune::show_best(metric = "roc_auc",n = 10) %>% 
  tidyr::pivot_longer(penalty:mixture, names_to="variable",values_to="value" ) %>% 
  ggplot(aes(value,mean)) + 
  geom_line(alpha=1/2)+ 
  geom_point()+ 
  facet_wrap(~variable,scales = "free")+
  ggtitle("Best parameters for AUC")


## Finalise model and predict test data ====

#Get the best parameters
en_best_params <- 
  en_bo %>%
  tune::select_best("roc_auc")

#Finalize our modelling workflow
en_wf_final <- 
  en_wf %>% 
  finalize_workflow(en_best_params)

#Fit on the split data, testing on the previously unseen test data
en_final_fit =
  en_wf_final |>
  fit(data = d_train)

#Plot our predictions
d_test |>
  bind_cols(
    predict(en_final_fit, d_test, type = "prob")
  ) %>% 
  roc_curve(group, .pred_Control) %>%
  autoplot()

#Tabulate performance
d_test |>
  bind_cols(
    predict(en_final_fit, d_test, type = "prob")
  ) %>%
  roc_auc(group, .pred_Control)

d_test |>
  bind_cols(
    predict(en_final_fit, d_test, type = "prob")
  )%>%
  average_precision(group,.pred_Control )


#This model can also get extremely good performance - perhaps even better than the lightgbm


#Make a plot of classification performance at difference thresholds

threshold_data = 
  d_test |>
  bind_cols(
    predict(en_final_fit, d_test, type = "prob")
  ) %>% 
  probably::threshold_perf(group,.pred_Control,
                           thresholds = seq(0, 1, by = 0.0025)) %>%
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

ggplot(threshold_data, aes(x = .threshold, y = .estimate, color = .metric)) +
  geom_line(size = 1.25) +
  theme_bw() +
  # scale_color_viridis_d(end = 0.9) +
  scale_alpha_manual(values = c(.4, 1), guide = "none") +
  geom_vline(xintercept = max_j_index_threshold, alpha = .6, color = "grey30",
             size = 1.25, lty = 2) +
  labs(
    x = "Threshold",
    y = "Metric Estimate",
    title = "Performance by threshold"  ) +
  theme(panel.grid = element_blank())


## Variable importance

#And you can get the importance from this model too...
# nnet_vars =
#   nnet_final_fit |>
#   extract_fit_engine() |>
#   vip::vi_permute(train = d_train,
#                   target = "group",
#                   metric = "auc",
#                   pred_wrapper = predict,
#                   reference_class = "Control")
#   
# 
# nnet_vars |> 
#   arrange(-Importance)|>
#   ggplot(aes(y = forcats::fct_reorder(Variable,Importance,max), x = Importance)) +
#   geom_col() 
# 
# 
# nnet_vars |> 
#   slice_max(Importance,n = 30) |>
#   knitr::kable()


en_explainer = 
  en_final_fit |>
  extract_fit_parsnip() |>
  explain_tidymodels(data  = d_train |> select(-group),
                     y     = d_train |> 
                       select(group) |> 
                       mutate(group = ifelse(group == "Control",1,0)) |>
                       pull()
  ) 

en_vars =  
  en_explainer|> 
  model_parts(type = "variable_importance")


en_vars |> 
  plot()


en_vars_30 = 
  en_vars |>
  as_tibble() |> 
  group_by(variable) |>
  summarise(imp_mu = mean(dropout_loss),
            imp_sd = sd(dropout_loss)) |>
  filter(variable %nin% c("_baseline_","_full_model_")) |>
  arrange(-imp_mu) |>
  dplyr::slice(1:30)

en_vars_30 |>
  knitr::kable()



# Compare variable importances ======


# Look at how much the variables from nnet, svm and light gbm overlap

all_vars = 
  bind_rows(
    lgbm_vars|> 
        select(Feature) |>
      rename(variable = Feature) |>
      mutate(model = "light_gbm"),
    svm_l_vars_30 |> 
      select(variable) |>
      mutate(model = "svm"),
    nnet_vars_30 |> 
      select(variable) |>
      mutate(model = "nnet"),
    en_vars_30 |> 
      select(variable) |>
      mutate(model = "en")
    )


all_vars |>
  count(variable) |>
  ggplot(aes(y = forcats::fct_reorder(variable,n,max), x = n)) +
  geom_col() 

#So we can see here that only 8 variables appear in both top 30 lists, suggesting that the two models are
#deriving their good performance from multiple variables - perhaps there is redundancy between variables

#We could try and fit a model using only the variables included in > 1 model
reduced_var_rec = 
  d_train |>
  recipes::recipe(all_vars |>
                    count(variable) |>
                    filter(n > 1) |>
                    pull(variable) |>
                    paste(collapse = " + ") %>%
                    paste("group ~ ", ., sep = "") |>
                    as.formula()
                  ) |>
  step_zv(all_predictors()) 



wf_par <- 
  workflow_set(
    preproc = list(reduced_var  = reduced_var_rec), 
    models  = list(lr.en        = en_mod,
                   SVM.linear   = svm_l_mod,
                   nnet         = nnet_mod,
                   lgbm         = lgbm_mod)
  )


wf_par_results <-
  wf_par |>
  workflow_map(
    "tune_bayes",
    seed = 1503,
    resamples = d_fold,
    control = control_bayes(verbose = TRUE,
                            save_pred = TRUE,
                            parallel_over = "everything",
                            save_workflow = TRUE),
    initial = 20,
    iter = 25,
    metrics = yardstick::metric_set(roc_auc, accuracy, f_meas),
  )



wf_par_results %>% 
  rank_results() %>% 
  filter(.metric == "roc_auc") %>% 
  select(model, .config, roc_auc = mean, rank)

autoplot(
  wf_par_results,
  rank_metric = "roc_auc",  # <- how to order models
  metric = "roc_auc",       # <- which metric to visualize
  select_best = TRUE     # <- one point per workflow
)



roc_mod <- 
  perf_mod(wf_par_results,
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
  ggplot(aes(x = forcats::fct_reorder(mod_type,posterior,max),
             y = posterior,ymin = .lower,ymax = .upper,group = vars,
             colour = mod_type)) +
  geom_errorbar(position = position_dodge(width = 0.5)) +
  geom_point(aes(shape = vars),position = position_dodge(width = 0.5)) +
  theme_bw() +
  theme(axis.text.x.bottom = element_text(size = 8, angle = -45)) +
  labs(y = "Model AUC", x = "Model Type",shape = "Variable \nSet")


#Do the ROPE Plot
autoplot(roc_mod, type = "ROPE", size = 0.025)

# Very similar

best_results <- 
  wf_par_results %>% 
  extract_workflow_set_result("reduced_var_lgbm") %>% 
  select_best(metric = "roc_auc")


boosting_test_results <- 
  wf_par_results %>% 
  extract_workflow("reduced_var_lgbm") %>% 
  finalize_workflow(best_results) %>% 
  last_fit(split = d_split_impute)


threshold_data = 
  boosting_test_results |>
  collect_predictions() |>
  probably::threshold_perf(group,.pred_Control,
                           thresholds = seq(0, 1, by = 0.0025)) %>%
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

ggplot(threshold_data, aes(x = .threshold, y = .estimate, color = .metric)) +
  geom_line(size = 1.25) +
  theme_bw() +
  # scale_color_viridis_d(end = 0.9) +
  scale_alpha_manual(values = c(.4, 1), guide = "none") +
  geom_vline(xintercept = max_j_index_threshold, alpha = .6, color = "grey30",
             size = 1.25, lty = 2) +
  labs(
    x = "Threshold",
    y = "Metric Estimate",
    title = "Performance by threshold"  ) +
  theme(panel.grid = element_blank())


#So, we get pretty much the same performance with our 24 variables, compared to the full dataset. 

#And the lightGBM model is the best, albeit in a very close race


boosting_test_results |>
  extract_fit_engine() |>
  lightgbm::lgb.importance() |>
  as_tibble() |>
  ggplot(aes(y = forcats::fct_reorder(Feature,Gain,max), x = Gain)) +
  geom_col() 

#DALEX doesn't work with light GBM in tidymodels at present so we will use the model-dervied 
#variable importance



# So the next step would be explore these variables interrelationships
