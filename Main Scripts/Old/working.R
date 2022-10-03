
# Model fitting scratchpad =====

plan(sequential)

# XGBoost model specification
xgboost_model <- 
  parsnip::boost_tree(
    mode = "classification",
    trees = 1000,
    min_n = tune(),
    tree_depth = tune(),
    learn_rate = tune(),
    loss_reduction = tune()
  ) %>%
  # set_engine("C5.0")
  set_engine("xgboost", objective = "binary:logistic")


# grid specification
xgboost_params <- 
  dials::parameters(
    min_n(),
    tree_depth(),
    learn_rate(),
    loss_reduction()
  )


xgboost_grid <- 
  dials::grid_max_entropy(
    xgboost_params, 
    size = 60
  )

# knitr::kable(head(xgboost_grid))

xgboost_wf <- 
  workflows::workflow() %>%
  add_model(xgboost_model) %>% 
  add_recipe(final_recipes$max_vars)

# hyperparameter tuning
xgboost_tuned <- 
  tune::tune_grid(
    object = xgboost_wf,
    resamples = d_folds,
    grid = xgboost_grid,
    metrics = yardstick::metric_set(accuracy,roc_auc),
    control = tune::control_grid(verbose = TRUE)
)

xgboost_tuned |> 
  collect_notes()

#So the issue looks like something to do with the parallel processing?

xgboost_tuned %>%
  tune::show_best(metric = "roc_auc") 


xgboost_best_params <- 
  xgboost_tuned %>%
  tune::select_best("roc_auc")

#Can we tune with racing?
xgboost_race <- 
    finetune::tune_race_anova(
      object = xgboost_wf,
      resamples = d_folds,
      grid = xgboost_grid,
      metrics = yardstick::metric_set(accuracy,roc_auc),
      control = finetune::control_race(save_pred = TRUE,
                                       save_workflow = TRUE)
  )

xgboost_race %>%
  tune::show_best(metric = "roc_auc") 

# Lets try the bagged tree thing
bag_cart_spec <-
  bag_tree(mode = "classification",
           min_n = tune()) %>%
  set_engine("C5.0")


# Make a workflow
bag_cart_wf <- 
  workflows::workflow() %>%
  add_model(bag_cart_spec) %>% 
  add_recipe(final_recipes$max_vars)


# grid specification
bag_cart_wf |>
  extract_parameter_set_dials()


bag_cart_params <- 
  dials::parameters(
    min_n())


bag_cart_grid <- 
  dials::grid_max_entropy(
    bag_cart_params, 
    size = 60
  )


#Can we tune with grid search?
bag_cart_tuned <- 
  tune::tune_grid(
    object = bag_cart_wf,
    resamples = d_folds,
    grid = bag_cart_grid,
    metrics = yardstick::metric_set(accuracy,roc_auc),
    control = tune::control_grid(verbose = TRUE)
  )

bag_cart_tuned %>%
  tune::show_best(metric = "roc_auc") 


bag_cart_best_params <- 
  bag_cart_tuned %>%
  tune::select_best("roc_auc")


#Can we tune with racing?
bag_cart_race <- 
  finetune::tune_race_anova(
    object = bag_cart_wf,
    resamples = d_folds,
    grid = 25,
    metrics = yardstick::metric_set(accuracy,roc_auc),
    control = finetune::control_race(save_pred = TRUE,
                                     save_workflow = TRUE)
  )



# Model ensembles ======

library(tidymodels)
library(tidyverse)
library(stacks)

data("tree_frogs")

# subset the data
tree_frogs <- 
  tree_frogs %>%
  select(-c(clutch, latency))



ggplot(tree_frogs) +
  aes(x = treatment, y = age, color = reflex) +
  geom_jitter() +
  labs(y = "Embryo Age (s)", 
       x = "treatment",
       color = "Response")



# some setup: resampling and a basic recipe
set.seed(1)

tree_frogs_split <- initial_split(tree_frogs)
tree_frogs_train <- training(tree_frogs_split)
tree_frogs_test  <- testing(tree_frogs_split)

folds <- rsample::vfold_cv(tree_frogs_train, v = 5)

tree_frogs_rec <- 
  recipe(reflex ~ ., data = tree_frogs_train) %>%
  step_dummy(all_nominal(), -reflex) %>%
  step_zv(all_predictors())

tree_frogs_wflow <- 
  workflow() %>% 
  add_recipe(tree_frogs_rec)


ctrl_grid <- control_stack_grid()



rand_forest_spec <- 
  rand_forest(
    mtry = tune(),
    min_n = tune(),
    trees = 500
  ) %>%
  set_mode("classification") %>%
  set_engine("ranger")

rand_forest_wflow <-
  tree_frogs_wflow %>%
  add_model(rand_forest_spec)

rand_forest_res <- 
  tune_grid(
    object = rand_forest_wflow, 
    resamples = folds, 
    grid = 10,
    control = ctrl_grid
  )


nnet_spec <-
  mlp(hidden_units = tune(), penalty = tune(), epochs = tune()) %>%
  set_mode("classification") %>%
  set_engine("nnet")

nnet_rec <- 
  tree_frogs_rec %>% 
  step_normalize(all_predictors())

nnet_wflow <- 
  tree_frogs_wflow %>%
  add_model(nnet_spec)

nnet_res <-
  tune_grid(
    object = nnet_wflow, 
    resamples = folds, 
    grid = 10,
    control = ctrl_grid
  )



tree_frogs_model_st <- 
  # initialize the stack
  stacks() %>%
  # add candidate members
  add_candidates(rand_forest_res) %>%
  add_candidates(nnet_res) %>%
  # determine how to combine their predictions
  blend_predictions() %>%
  # fit the candidates with nonzero stacking coefficients
  fit_members()

tree_frogs_model_st


autoplot(tree_frogs_model_st) 
autoplot(tree_frogs_model_st, type = "members")
autoplot(tree_frogs_model_st, type = "weights")


collect_parameters(tree_frogs_model_st, "rand_forest_res")


tree_frogs_pred <-
  tree_frogs_test %>%
  bind_cols(predict(tree_frogs_model_st, ., type = "prob"))



yardstick::roc_auc(
  tree_frogs_pred,
  truth = reflex,
  contains(".pred_")
)




tree_frogs_pred <-
  tree_frogs_test %>%
  select(reflex) %>%
  bind_cols(
    predict(
      tree_frogs_model_st,
      tree_frogs_test,
      type = "class",
      members = TRUE
    )
  )


map_dfr(
  setNames(colnames(tree_frogs_pred), colnames(tree_frogs_pred)),
  ~mean(tree_frogs_pred$reflex == pull(tree_frogs_pred, .x))
) %>%
  pivot_longer(c(everything(), -reflex))


# Probability
p_d = 
  crossing(pops  = 10^(seq(2,7,0.1)) |> round(),
           probs = 10^(seq(-5,0,0.5))) |>
  mutate(dist = map2(probs,pops,~rbinom(10000,.y,.x))) |>
  # mutate(dist = map2(dist,pops,~.x/.y)) |>
  mutate(prob_occur = map_int(dist, ~(.x > 0) %>% sum())) |>
  mutate(prob_occur = prob_occur/10000)



p_d |>
  mutate(perc = round(probs * 100,digits = 3)) |>
  ggplot(aes(x = pops,y = prob_occur, group = factor(perc),
             colour = factor(perc))) +
  geom_line() +
  geom_point() +
  labs(x = "Population Size", y = "Probability of > 1 hit",color = "Allele Frequency (%)") +
  scale_x_log10() +
  geom_hline(yintercept = 0.95,lty = 2) +
  geom_vline(xintercept = 10000) +
  theme_minimal()


