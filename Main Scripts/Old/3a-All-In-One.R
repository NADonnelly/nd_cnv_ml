
# Introduction =========

#What we are going to do here is to combine imputation, variable selection
#and model fitting into one cross-validation process (possibly nested if I
#can work that out).

#So - we will start with the data filtered of the non-quantitative variables
# - Then we will do our initial test/training split
# - Then with the training data we will split into nested CV
# - The inner folds will be used for the process of training imputation and
#   variable selection models. I might also pop PLS in.
# - Next we tune machine learning models as models x variable sets
# - Use the outer folds to evaluate model performance

#This does raise problems about comparing models though because if we do
#variable selection within CV loops, how do we compare performance?
#For the actual goals of this study, we are trying to both make a model that 
#makes the best predictions on out-of-sample test data, and select a reduced
#set of variables to make a screening questionnaire which is able to preserve
#much of the accuracy that you get from having access to a full database of
#research items, which are not practical to gather in the clinical setting

#Of these objectives, the second is perhaps most important



# Load packages and data ========

#Lets load up
pacman::p_load(tidyverse,tidymodels,workflowsets,tidyposterior,
               finetune,probably,stacks,doFuture,baguette,furrr,
               rstanarm,projpred, bonsai, learntidymodels, vip,
               ggforce,tabnet,discrim)

tidymodels_prefer()



#Load the full dataset
DF = read_csv("FilteredData2.csv")



# Initial data splits ========

#Make the class labels
DF0 = 
  DF |>
  mutate(group = map_chr(GenotypeCode,~ifelse(.x == 16,"Control","ND-CNV"))) |>
  mutate(group = factor(group, levels = c("Control","ND-CNV"))) |>
  select(-c(IDs,GenotypeCode)) |>
  relocate(group)

#Have as quick look at all variables
# skimr::skim(DF0)


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
split_prop = 4/5

n_outer    = 10
n_inner    = 20
n_fold     = 5


# Prepare our initial data splits
d_split = initial_split(DF1, strata = group,prop = split_prop)

#Get test and training data
d_train = training(d_split)
d_test  = testing(d_split)

#Prepare for (nested) cross validation
d_folds = nested_cv(d_train, 
                    outside = vfold_cv(v = n_fold, repeats = n_outer/n_fold), 
                    inside  = bootstraps(times = n_inner)) 

# So the state of play is that we have:

# - The full dataset with 489 individuals and 233 variables
# - The testing holdout dataset, which we will not be looking at until the very end, with 99 individuals
# - Nested CV datasets; each outer fold contains 390 individuals, split into 351 analysis, 39 assessment
# - Nested CV inner folds contain 351 individuals, 123 assessment, 351 analysis bootstraps

# We might have to start doing this in a messy way and then develop some functions to
# streamline the code. But lets start



# All variables ======

#I would like to begin with a look at the full variable set. We would be doing any analysis on each folder
#but lets look at all the training data together at first



## Fit models =====

#OK so now we are going to do ML to the full variable set. To do nested CV we are going to loop
#I know loops are bad etc but this is a good way to do things quickly. Plus I'm used to MATLAB

#Preallocate an output tibble
d_fold_results = 
  d_folds |>
  select(id,id2) |>
  mutate(.results       = vector(mode = "list",length = nrow(d_folds)),
         .mod_posterior = vector(mode = "list",length = nrow(d_folds)))

#Loop through each outer fold
for(i in 1:nrow(d_folds)){
  
  
  #Get the outer data for this iteration of the loop
  d_outer = d_folds$splits[[i]]
  
  #Get the inner data for this iteration of the loop
  d_inner = d_folds$inner_resamples[[i]]
  
  
  
  #Make pre-processing recipes

  #Set up a recipe without PLS, just all the variables all at once, with imputation
  rec_impute =
    d_outer |>
    analysis() |>
    recipe(group ~ .) |>
    step_zv(all_predictors()) |>
    step_impute_bag(all_predictors())   
  
  
  #Set a recipe for our initial data preparation including imputation and PLS
  rec_impute_pls = 
    d_outer |>
    analysis() |>
    recipe(group ~ .) |>
    step_zv(all_predictors()) |>
    step_impute_bag(all_predictors()) |>
    step_pls(all_predictors(),outcome = vars(group),num_comp = 2)
  
  
  #Set up a recipe without any imputation
  rec_simple = 
    d_outer |>
    analysis() |>
    recipe(group ~ .) |>
    step_zv(all_predictors())
  
  #Set up a recipe without any imputation, with with pls
  rec_simple_pls = 
    d_outer |>
    analysis() |>
    recipe(group ~ .) |>
    step_zv(all_predictors()) |>
    step_pls(all_predictors(),outcome = vars(group),num_comp = 2)
  

  #We are going to set the number of PLS components to 2 because in
  #testing it looks like we always pick 2 when we fit, and having tuning in 
  #a recipe makes things horrid in several places lower down
  
  
  #Now, rather than impute on every single model training step,
  #we are going to do the imputation on this outer fold here
  #I'm not sure if this is imperfect, but it is going to save a 
  #lot of computational time

  #Prep a model for imputing missing data
  rec_impute_prep = 
    rec_impute |>
    prep()
  
  #Get our imputed data back
  d_impute = 
    rec_impute_prep |>
    juice()
  
  
  #We can make bootstraps from the imputed dataset (replaces our existing bootstrapped inner samples)
  d_inner = 
    d_impute |>
    bootstraps(times = n_inner)
  
  #And make a test set for this outer loop using the same imputation model, but the test (assessment) 
  #data
  d_outer_test = 
    rec_impute_prep |>
    bake(new_data = 
           d_outer |> 
           assessment())
  
  #We could maybe rebuild a split object with the imputed data
  d_outer_split = 
    make_splits(d_impute,d_outer_test)
  
  

  
  #Set up a bunch of potential models
  
  
  #Plain old LR - will work horribly
  glm_spec <- 
    logistic_reg(mode = "classification") %>% 
    set_engine("glm")
  
  
  #Bayesian LR with RSTANARM
  glm_bayes_spec <- 
    logistic_reg(mode = "classification") %>% 
    set_engine("stan",  
               seed = 1, refresh = 0,
               chains = 4, cores = 16,
               iter = 4000, warmup = 1000,
               prior = rstanarm::normal(0,1),
               prior_intercept = rstanarm::normal(0,1.5),
               prior_aux = rstanarm::exponential(rate = 1))
  
  
  #Bayesian LR with RSTANARM and the Horseshoe prior - removed for being very slow to fit and not
  #notably different to the Bayesian GLM with and normal(0,1) prior, which is also weakly regularising
 
  # #Make our horseshoe prior
  # n <- nrow(d_impute)
  # D <- ncol(d_impute)-1
  # 
  # #How many variables do you think will be useful predictors? Lets say 30 as a 
  # #starting value
  # p0 <- 30
  # 
  # tau0 <- p0/(D - p0) * 1/sqrt(n)
  # prior_coeff <- hs(global_scale = tau0, slab_scale = 1)
  # 
  # #Right well, lets fit it
  # glm_horseshoe_spec <- 
  #   logistic_reg(mode = "classification") %>% 
  #   set_engine("stan",  
  #              seed = 1, refresh = 0,
  #              chains = 4, cores = 4,
  #              iter = 4000, warmup = 1000,
  #              prior = prior_coeff,
  #              refresh = 1,
  #              adapt_delta = 0.999,
  #              prior_intercept = rstanarm::normal(0,1.5),
  #              prior_aux = rstanarm::exponential(rate = 1))
  
  
  #Elastic Net Regression with glmnet
  logistic_reg_spec <- 
    logistic_reg(mode = "classification",
                 penalty = tune(), 
                 mixture = tune()) %>% 
    set_engine("glmnet")
  
  
  #Linear Discriminant analysis
  lda_spec <- 
    discrim_linear(mode = "classification",
                   engine = "mda",
                   penalty = tune())
  
  
  #Flexible Discriminant Analysis
  fda_spec <- 
    discrim_flexible(prod_degree = tune(),
                     mode = "classification",
                     engine = "earth")
  
  
  #Regularised Discriminant Analysis
  rda_spec <- 
    discrim_regularized(mode = "classification",
                        engine = "klaR",
                        frac_common_cov = tune(),
                        frac_identity   = tune())
  
  
  #Linear SVM
  svm_l_spec = 
    svm_linear(cost   = tune(),
               margin = tune()) %>%
    set_engine("kernlab") %>% 
    set_mode("classification")
  
  
    # Radial SVM
  svm_r_spec <- 
    svm_rbf(cost     = tune(), 
            rbf_sigma = tune()) %>% 
    set_engine("kernlab") %>% 
    set_mode("classification")

  #Polynomial SVM
  svm_p_spec <- 
    svm_poly(cost   = tune(), 
             degree = tune()) %>% 
    set_engine("kernlab") %>% 
    set_mode("classification")
  
  
  #Neural Network with nnet
  nnet_spec <- 
    mlp(hidden_units = tune(), 
        penalty      = tune(), 
        epochs       = tune()) %>% 
    set_engine("nnet", MaxNWts = 2600) %>% 
    set_mode("classification")
  
  
  #Neural Network with tabnet
  tn_spec <- 
    tabnet(epochs = 50, batch_size = 1024) %>%
    set_engine("torch", verbose = TRUE) %>%
    set_mode("classification")
  
  #Of course there are lots of settings to tune in a tabnet model, but I will
  #use the defaults because I am currently too ignorant to do otherwise, and 
  #its a computational beast to fit  
  
  #Tabnet with tuning as suggested by RStudio https://blogs.rstudio.com/ai/posts/2021-02-11-tabnet/
  tn_spec <- 
    tabnet(batch_size = 16384, 
           epochs          = tune(), 
           decision_width  = tune(), 
           attention_width = tune(),
           num_steps       = tune(), 
           penalty         = tune(),
           learn_rate      = tune()) %>%
    set_engine("torch", verbose = TRUE) %>%
    set_mode("classification")
  
  #Neural Network with Brulee/Torch (tabnet also uses Torch as the backend and is
  #specifically meant to be good for tabular data )
  brulee_spec <- 
    mlp(hidden_units = tune(),
        penalty      = tune(),
        epochs       = tune(),
        learn_rate   = tune(),
       # mixture      = tune(),
        activation   = "softmax")   %>% 
    set_engine("brulee") %>% 
    set_mode("classification")%>% 
    translate()
  

  #MARS
  mars_spec <- 
    mars(prod_degree = tune()) %>%  #<- use GCV to choose terms
    set_engine("earth") %>% 
    set_mode("classification")
  
  
  #KNN
  knn_spec <- 
    nearest_neighbor(neighbors   = tune(), 
                     dist_power  = tune(), 
                     weight_func = tune()) %>% 
    set_engine("kknn") %>% 
    set_mode("classification")
  
  
  #Naive Bayes
  bayes_spec <- 
    naive_Bayes(smoothness = tune()) %>%
    set_engine('klaR')
  
  
  #CART
  cart_spec <- 
    decision_tree(cost_complexity = tune(), 
                  min_n = tune()) %>% 
    set_engine("rpart") %>% 
    set_mode("classification")
  
  
  #Bagging CART
  bag_cart_spec <-
    bag_tree(mode = "classification",
             min_n = tune()) %>%
    set_engine("C5.0")
  
  
  #Random Forest
  rf_spec <- 
    rand_forest(trees = 500,min_n = tune()) %>% 
    set_engine("ranger") %>% 
    set_mode("classification")
  
  
  #Conditional random forest
  cond_forest_spec <- 
    rand_forest(trees = 500,min_n = tune()) %>%
    set_engine("partykit") %>%
    set_mode("classification") 
  
  #Note we are not going to tune mtry because for the pls it doesnt make sense
  
  
  #Boosted Trees with XGBoost
  xgb_spec <- 
    parsnip::boost_tree(
      mode = "classification",
      trees = 500,
      min_n = tune(),
      tree_depth = tune(),
      learn_rate = tune(),
      loss_reduction = tune() ) %>%
    set_engine("xgboost", objective = "binary:logistic")  
  
  
  #Boosted Trees with lightgbm
  light_gbm_spec <- 
    parsnip::boost_tree(
      mode = "classification",
      engine = "lightgbm",
      trees = 500,
      min_n = tune(),
      tree_depth = tune(),
      learn_rate = tune(),
      loss_reduction = tune()   )
  
  
  #Make workflow sets
  
  #Some models will run in parallel happily, some produce problems
  
  #These are the ones that fit in parallel
  wf_par <- 
    workflow_set(
      preproc = list(impute.pls   = rec_simple_pls,
                     impute       = rec_simple), 
      models  = list(lr.glm       = glm_spec,
                     lr.en        = logistic_reg_spec,
                     lda          = lda_spec,
                     fda          = fda_spec,
                     rda          = rda_spec,
                     SVM.radial   = svm_r_spec, 
                     SVM.linear   = svm_l_spec,
                     SVM.poly     = svm_p_spec,
                     nnet.n       = nnet_spec,
                     nb           = bayes_spec,
                     mars         = mars_spec,
                     knn          = knn_spec,
                     cart         = cart_spec,
                     rand.forest  = rf_spec,
                     cond.forest  = cond_forest_spec)
    )
  
  
  # #We need to so an extra step to make our forests work with the tuning pls - or we just remove the mtry argument
  # Or we remove the tuning on the pls
  
  #  cond_forest_params = 
  #    wf_par |> 
  #     extract_workflow(id = "impute.pls_cond.forest") |>
  #     extract_parameter_set_dials() |>
  #     update(`mtry` = finalize(mtry(),d_impute))
  #   
  # rand_forest_params = 
  #    wf_par |> 
  #    extract_workflow(id = "impute.pls_rand.forest") |>
  #    extract_parameter_set_dials() |>
  #    update(`mtry` = finalize(mtry(),d_impute))
  # 
  # fda_params = 
  #   wf_par |>
  #   extract_workflow(id = "impute.pls_fda") |>
  #   extract_parameter_set_dials() |>
  #   update(`num_terms ` = finalize(num_terms (),d_impute))
  #    
  # wf_par = 
  #   wf_par |>
  #   option_add(param_info = cond_forest_params, id = "impute.pls_cond.forest") |>
  #   option_add(param_info = rand_forest_params, id = "impute.pls_rand.forest") |>
  #   option_add(param_info = fda_params        , id = "impute.pls_fda")
  
  
  #These have to be fit not in parallel for some reason
  wf_seq <- 
    workflow_set(
      preproc = list(impute.pls   = rec_simple_pls,
                     impute       = rec_simple), 
      models  = list(lr.bayes     = glm_bayes_spec,
                     # lr.horseshoe = glm_horseshoe_spec,
                     bagged.tree  = bag_cart_spec,
                     xgb.boost    = xgb_spec,
                     lgb.boost    = light_gbm_spec,
                     tabnet       = tn_spec,
                     brulee       = brulee_spec))
  
  #You can look at the full recipe for any of these e.g.
  # wf_par |> extract_workflow(id = "impute.pls_lr.en")

  #You can look at the paramters used thus:
     # wf_par |>
     # extract_workflow(id = "impute.pls_rand.forest") |>
     # extract_parameter_set_dials() |>
     #   pull(object)
  
  
  
  ## Do some model fitting

  #Set up controls for the finetune anova race method
  race_ctrl <-
    control_race(
      save_pred = TRUE,
      parallel_over = "everything",
      save_workflow = TRUE
    )
  
  
  #Set up to run these in parallel
  registerDoFuture()
  plan(multisession, workers = 8)
  
  #Fit our workflows that work in parallel
  wf_par_results <-
    wf_par |>
    workflow_map(
      "tune_race_anova",
      seed = 1503,
      resamples = d_inner,
      control = race_ctrl,
      grid = 25,
      metrics = metric_set(roc_auc)
    )
  

  #What happens if we use tunebayes instead? Results look very similar,
  #is the answer to that
  # wf_par_results_bo <-
  #   wf_par |>
  #   workflow_map(
  #     "tune_bayes",
  #     seed = 1503,
  #     resamples = d_inner,
  #     control   = control_bayes(
  #       verbose = TRUE,
  #       save_pred = TRUE,
  #       parallel_over = "everything",
  #       save_workflow = TRUE),
  #     initial = 20,
  #     iter = 20,
  #     metrics = metric_set(roc_auc)
  #   )
  
  
  
  
  
  #Turn off parallel mode
  plan(sequential)

  #Fit the workflows that seem to work best sequentially
  wf_seq_results <-
    wf_seq |>
    workflow_map(
      "tune_race_anova",
      seed = 1503,
      resamples = d_inner,
      control = race_ctrl,
      grid = 25,
      metrics = metric_set(roc_auc),
      verbose = TRUE
    )


  
  #Unite the datasets
  wf_results = 
    bind_rows(wf_par_results,wf_seq_results)
  
  # #Rank models by their resampling results
  # rankings <-
  #   wf_results |>
  #   rank_results(select_best = TRUE)
  # 
  # filter(rankings, rank <= 5) %>% 
  #   dplyr::select(rank, mean, model, wflow_id)  
  
  
  #Plot individual model's fitting characteristics
  # autoplot(wf_results, metric = "roc_auc", id = "impute.pls_SVM.linear")
  
  
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
    arrange(-posterior)

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

  
  #Get the settings for the  best version of every model, fit to
  #the outer fold training data, and then test on the outer folder test
  #data
  
  best_results_wf = 
    tibble(wflow_id = wf_results$wflow_id) |>
    mutate(res         = purrr::map(wflow_id,~workflowsets::extract_workflow_set_result(wf_results,id = .x)))  |>
    mutate(best_params = purrr::map(res, tune::select_best,metric = "roc_auc")) |>
    mutate(best_wf     = purrr::map(wflow_id, ~workflowsets::extract_workflow(wf_results,id = .x))) 
  
  
  best_results_wf = 
    best_results_wf |>

    #We are getting a problem for models were we haven't tuned anything. Lets make
    #a tiny bit of logic about this
    mutate(npar = map_lgl(best_params,~ncol(.x) > 1))   |>
    mutate(best_wf_final  = pmap(list(a = npar, b = best_wf, c = best_params),
                                 function(a,b,c){
                                   if(a){
                                     wf_final = finalize_workflow(b,c)
                                     }else{
                                       wf_final = b
                                       }
                                   return(wf_final)
                                   })) |>
          
    #Fit the best models to the full training data for this outer fold       
    mutate(outer_train_fit    = purrr::map(best_wf_final, last_fit,split = d_outer_split)) 
    
    
    # mutate(outer_train_fit    = purrr::map(best_wf_final,fit,data = d_impute)) |>           
    # 
    # #Apply the best models to the test data for this outer fold (not previously seen by any
    # #part of  the fitting process)
    # mutate(outer_test_metrics = purrr::map(outer_train_fit,~d_outer_test |>
    #                                          bind_cols(
    #                                            predict(.x, d_outer_test, type = "prob")
    #                                            ) %>%
    #                                          roc_auc(group, .pred_Control)))
  
  
  #We would also like to derive variable importance measures in each loop. There are different
  #ways to do this. Its quickest to use model-based variable importance, but we can also use 
  #permutation. Also notably two of the newer models - light gbm and tabnet do not work with
  #the existing vip::vi or DALEXtra::explain frameworks, but do have their own internal variable
  #importance functions which we can use.
  
  
  
  
  #Store the results for this fold
  d_fold_results$.results[[i]] =
    best_results_wf |>
    dplyr::select(-c(res,best_wf)) 
  
  
  #Store the performance model data too
  d_fold_results$.mod_posterior[[i]] = 
    roc_mod_post
  
  #lets do some cleaning up
  rm(list = c("d_impute","d_inner", "d_outer_test", "d_outer",
              "rec_impute","rec_impute_pls","rec_simple","rec_simple_pls",
              "wf_par","wf_seq","wf_results","wf_seq_results","wf_par_results",
              "best_results_wf",
              "roc_mod","roc_mod_post"))
  
  
  # #Tabulate our metrics
  # best_results_wf |>
  #   unnest(outer_test_metrics) |>
  #   filter(.metric == "roc_auc") |>
  #   select(wflow_id,.metric,.estimate) |>
  #   arrange(-.estimate) |>
  #   knitr::kable()
  # 
  # #Or plot I suppose
  # best_results_wf |>
  #   unnest(outer_test_metrics) |>
  #   filter(.metric == "roc_auc") |>
  #   dplyr::select(wflow_id,.metric,.estimate) |>
  #   ggplot(aes(x = forcats::fct_reorder(wflow_id,.estimate,max),
  #              y = .estimate)) +
  #   geom_col() +
  #   theme_bw() +
  #   theme(axis.text.x.bottom = element_text(angle = 45))
  

  
}


#Save this (it takes 12 hours or so to fit!)
write_rds(d_fold_results,"nested_cv_result_all_vars.rds")

d_fold_results = read_rds("nested_cv_result_all_vars.rds")

#Lets plot some results
d_fold_results |> 
  select(-.results) |>
  unnest(.mod_posterior) |>
  ggplot(aes(x = forcats::fct_reorder(mod_type,posterior,max),
             y = posterior,ymin = .lower,ymax = .upper,group = vars,
             colour = mod_type)) +
  geom_errorbar(position = position_dodge(width = 0.5)) +
  geom_point(aes(shape = vars),position = position_dodge(width = 0.5)) +
  theme_bw() +
  theme(axis.text.x.bottom = element_text(size = 8, angle = -45)) +
  labs(y = "Model AUC", x = "Model Type",shape = "Variable \nSet")


#Look at the averaged performance on the outer fold test data
d_fold_results |>
  select(-.mod_posterior) |>
  unnest(.results) |>
  unnest(outer_test_metrics) |>
  select(id,id2,wflow_id,.estimate) |> 
  mutate(id = interaction(id,id2,sep = "_")) |>
  select(-id2) |>
  separate(wflow_id,into = c("vars","mod_type"),sep = "_", extra = "merge") |>
  group_by(vars,mod_type) |>
  summarise(m_p = mean(.estimate)) |>
  arrange(-m_p) |>
  knitr::kable()
  

d_fold_results |>
  select(-.mod_posterior) |>
  unnest(.results) |>
  unnest(outer_test_metrics) |>
  select(id,id2,wflow_id,.estimate) |> 
  mutate(id = interaction(id,id2,sep = "_")) |>
  select(-id2) |>
  separate(wflow_id,into = c("vars","mod_type"),sep = "_", extra = "merge") |>
  group_by(vars,mod_type) |>
  summarise(m_p = mean(.estimate),
            s_p = sd(.estimate)) |>
  filter(!str_detect(vars,"pls")) |>
  ggplot(aes(x = forcats::fct_reorder(mod_type,m_p,max),
             y = m_p,ymin = m_p-s_p,ymax = m_p + s_p,
             group = vars,
             colour = mod_type,
             fill   = mod_type)) +
  geom_point(aes(shape = vars),position = position_dodge(width = 0.5)) +
  geom_linerange(position = position_dodge(width = 0.5)) +
  theme_bw() +
  theme(axis.text.x.bottom = element_text(size = 8, angle = -45)) +
  labs(y = "Model AUC", x = "Model Type",shape = "Variable \nSet")


#Do something like the perf_mod function in tidymodels
cv_results = 
  d_fold_results |>
  select(-.mod_posterior) |>
  unnest(.results) |>
  unnest(outer_test_metrics) |>
  select(id,id2,wflow_id,.estimate) |> 
  mutate(id = interaction(id,id2,sep = "_")) |>
  select(-id2)


#We can directly compare with tidyposterior, which takes the best submodel from
#each workflow and then does a stan_glm on the performance metrics
cv_mod <- 
  stan_glmer(.estimate ~ 0 + wflow_id + (1|id),
              data = cv_results,
              seed = 1, refresh = 0,
              chains = 4, cores = 8,
              iter = 10000, warmup = 2000,
              prior = rstanarm::normal(0,1),
              prior_intercept = rstanarm::normal(0,1.5),
              prior_aux = rstanarm::exponential(rate = 1))

cv_mod |> summary()


#Make our own plot
cv_mod_post = 
  cv_mod |>
  broom.mixed::tidy( conf.int = T, conf.level = 0.95, conf.method = "HPDinterval") |>
  mutate(term = str_remove(term,"wflow_id")) |>
  separate(term,into = c("vars","mod_type"),sep = "_", extra = "merge") 


cv_mod_post |>
  ggplot(aes(x = forcats::fct_reorder(mod_type,estimate,max),
             y = estimate,group = vars,ymin =  conf.low,ymax = conf.high,
             colour = mod_type,
             fill   = mod_type)) +
  geom_errorbar(position = position_dodge(width = 0.5)) +
  geom_point(aes(shape = vars),position = position_dodge(width = 0.5)) +
  theme_bw() +
  theme(axis.text.x.bottom = element_text(size = 8, angle = -45)) +
  labs(y = "Model AUC", x = "Model Type",shape = "Variable \nSet")

#Get rid of PLS, which models do best with all variables (so we can get variable importance)
cv_mod_post |>
  filter(!str_detect(vars,"pls")) |>
  ggplot(aes(x = forcats::fct_reorder(mod_type,estimate,max),
             y = estimate,group = vars,ymin =  conf.low,ymax = conf.high,
             colour = mod_type,
             fill   = mod_type)) +
  geom_errorbar(position = position_dodge(width = 0.5)) +
  geom_point(aes(shape = vars),position = position_dodge(width = 0.5)) +
  theme_bw() +
  theme(axis.text.x.bottom = element_text(size = 8, angle = -45)) +
  labs(y = "Model AUC", x = "Model Type",
       title = "All variable models only")




## Get variable importance ====


#We can use permutation based variable importance. We might be best off using the full 
#training set - but will we need to impute it. Lets do that

rec_impute =
  d_train |>
  recipe(group ~ .) |>
  step_zv(all_predictors()) |>
  step_impute_bag(all_predictors())   

rec_impute_prep = 
  rec_impute |>
  prep()

#Get our imputed data back
d_impute = 
  rec_impute_prep |>
  juice()


#Get variable importance using DALEX
d_fold_vi = 
  d_fold_results |>
  select(-.mod_posterior) |>
  unnest(.results) |>
  mutate(id = interaction(id,id2,sep = "_")) |>
  separate(wflow_id,into = c("vars","mod_type"),sep = "_", extra = "merge") |>
  select(id,vars,mod_type,outer_train_fit) |>
  filter(!str_detect(vars,"pls")) |>
  select(-vars)


d_fold_vi = 
  d_fold_vi|>
  filter(!(str_detect(mod_type,"tab") | str_detect(mod_type,"lgb")))|>
  mutate(v_explainer = map2(outer_train_fit,mod_type,
                            ~extract_fit_parsnip(.x) |> 
                             explain_tidymodels(label = .y,
                                                data  = d_impute |> select(-group),
                                                y     = d_impute |> 
                                                          select(group) |> 
                                                          mutate(group = ifelse(group == "Control",1,0)) |> 
                                                          pull()
                                                ) )) |>
  mutate(v_parts = map(v_explainer,  model_parts, type = "variable_importance"))


# So models won't work like this, so we will need to manually produce variable importances
# - Light GBM
# - Tabnet


#We can do that here






# Then we aggregate everything together



# And save - the permutation test takes a long term
write_rds(d_fold_vi,"nested_cv_variable_importance.rds")



# Then we need to ask ourselves the question - what do we do with this?









explainer_i = 
  explain_tidymodels(
    model = d_fold_vi$outer_train_fit[[i]] |>
      extract_fit_parsnip(),
    data  = d_impute |> select(-group),
    y     = d_impute |> select(group) |>mutate(group = ifelse(group == "Control",1,0)) |> pull(),
    label = d_fold_vi$mod_type[[i]]
  )

set.seed(1803)
vip_i = 
  explainer_i |>
  model_parts(type = "variable_importance")







d_fold_vi = 
  d_fold_results |>
  select(-.mod_posterior) |>
  unnest(.results) |>
  mutate(id = interaction(id,id2,sep = "_")) |>
  separate(wflow_id,into = c("vars","mod_type"),sep = "_", extra = "merge") |>
  select(id,vars,mod_type,outer_train_fit) |>
  filter(!str_detect(vars,"pls")) |>
  select(-vars) |>
  filter(!str_detect(mod_type,"da")) 



d_fold_vi$vi = vector(mode = "list", length = nrow(d_fold_vi))

for(i in 1:nrow(d_fold_vi)){
  
  if(!str_detect(d_fold_vi$mod_type[[i]],"da" )){
    
    d_fold_vi$vi[[i]] = 
      d_fold_vi$outer_train_fit[[i]] |>
      extract_fit_parsnip() |>
      vip::vi_model()
    
  }else{
    
    d_fold_vi$vi[[i]] = 
      d_fold_vi$outer_train_fit[[i]] |>
      extract_fit_parsnip() |>
      vip::vi(method = "permute",
            train = d_impute,
            target = "group",
            metric = "auc",
            pred_wrapper = function(object, newdata) predict(object, newdata),
            reference_class = "Control")    
    
    
    d_fold_vi$outer_train_fit[[i]] |> extract_fit_parsnip() |> caret::varImp()

  }

}

explainer_i = 
  explain_tidymodels(
    model = d_fold_vi$outer_train_fit[[i]] |>
      extract_fit_parsnip(),
    data  = d_impute |> select(-group),
    y     = d_impute |> select(group) |>mutate(group = ifelse(group == "Control",1,0)) |> pull(),
    label = d_fold_vi$mod_type[[i]]
)

set.seed(1803)
vip_i = 
  explainer_i |>
  model_parts(type = "variable_importance")


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
  
  p <- 
    full_vip %>%
    filter(variable != "_full_model_") %>% 
    mutate(variable = fct_reorder(variable, dropout_loss)) %>%
    ggplot(aes(dropout_loss, variable)) 
  
  if(length(obj) > 1) {
    
    p <- p + 
      facet_wrap(vars(label)) +
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


ggplot_imp(vip_i)


full_vip <- 
  bind_rows(vip_i) %>%
  filter(variable != "_baseline_")

perm_vals <- 
  full_vip %>% 
  filter(variable == "_full_model_") %>% 
  group_by(label) %>% 
  summarise(dropout_loss = mean(dropout_loss))



#Get variable importance
e_net_vi = 
  e_net_final_fit %>%
  extract_fit_parsnip() %>%
  vip::vi() %>%
  mutate(
    Importance = abs(Importance),
    Variable = fct_reorder(Variable, Importance)
  ) 

#Plot our important variables by importance
e_net_vi |> 
  filter(Importance > 0) |>
  ggplot(aes(x = Importance, y = Variable, fill = Sign)) +
  geom_col() +
  scale_x_continuous(expand = c(0, 0)) +
  labs(y = NULL)

#Unfortunately as we are doing PLS, we just get the PLS components back... 




# Variable selection =====

# So we can use the inner folds to do variable selection repeatedly as part of our broader ML approach

#Get an example

d_0_1 = d_folds$splits[[1]]
#So which methods for variable selection should we try? These are the methods I have used previously

# - Using a GLM and selecting the variables that are significant (variables fit individually)
# - Using a GLM and using the GLM to predict group with CV, ranking variables by AUC and taking the top 30
# - Using LASSO regression (with CV on the penalty) and selecting all variables with an importance > 0
# - Calculating the polychoric correlation between all variables and selecting only those with an average
#   correlation < 0.5 (to select hopefully the most informative set of variables)
# - Use a penalised Bayesian GLM, with the horseshoe prior set for 30 expected informative variables, and
#   taking the 30 variables with the largest posterior values
# - Extending the Bayesian GLM by using the projective prediction approach and identifying the top 30 variables

# - Other things - Boruta method?

# - Then some derived variable approaches: 
#   - all variables selected by more than one of the above methods
#   - those variables included in the largest number of sets (was 6 the first time)

# We would then do this 25 times to get final sets of variables which can be used for ML fitting?

#I suppose one approach would be to do the bootstrapping thing a very large number of times and select
#only variables with a value significantly different from 0.5. Lets think about that at a later point


## Single Variable GLM Significance =====

#This is straightforward enough for each outer fold. We simply take the training data for this outer folder like this:

nvar = ncol(  d_0_1 |>
                analysis() ) - 1

d_sig = 
  d_0_1 |>
  analysis() |>
  pivot_longer(-c(group),names_to = "variable",values_to = "value") |>
  group_by(variable) |>
  nest()

d_sig = 
  d_sig |>
  mutate(pv = map_dbl(data,~glm(group ~ value,
                                data = .x, 
                                family = binomial(link = "logit")) %>%
                        broom::tidy() |>
                        filter(term == "value") |>
                        pull(p.value))) |>
  mutate(pv_corr = pv*nvar)

vars_sig = 
  d_sig |>
  filter(pv_corr < 0.05) |>
  select(variable,pv_corr) |>
  mutate(type = "single_LR_sig")


## Single Variable GLM AUC =====

# We are going to work out the AUC for each variable independently 
d_auc = 
  d_0_1 |>
  analysis() |>
  pivot_longer(-c(group),names_to = "variable",values_to = "value") |>
  group_by(variable) |>
  nest() |>
  mutate(folds = map(data,bootstraps, strata = group,times = 25))


#Make a tidymodel for a glm
log_model = 
  logistic_reg(mode = "classification",
               engine = "glm") 

log_wflow = 
  workflow() |>
  add_model(log_model) |>
  add_formula(group ~ value)


#Make a control variable for model fitting
use_par <-
  control_resamples(save_pred = TRUE, save_workflow = TRUE, allow_par = TRUE,
                    parallel_over = "resamples")


#Fit bootstraps for classification accuracy to all variables
plan(multisession, workers = 14)

d_auc = 
  d_auc |>
  ungroup() |>
  mutate(var_metrics = furrr::future_map(folds,~log_wflow |> 
                                           fit_resamples(resamples = .x,
                                                         control = use_par) |> 
                                           collect_metrics() |> 
                                           dplyr::select(.metric,mean,std_err) |> 
                                           filter(.metric == "roc_auc")
  ))



#Extract the average AUC for each variable and plot
d_auc |> 
  dplyr::select(variable,var_metrics) |>
  unnest(var_metrics) |>
  arrange(-mean) |>
  ggplot(aes(x = fct_reorder(variable,mean,max), y = mean, ymin = mean-std_err,ymax = mean+std_err)) +
  geom_col() +
  geom_errorbar() +
  coord_cartesian(ylim = c(0.4,1)) +
  theme_bw() +
  geom_hline(yintercept = 0.5) +
  theme(axis.text.x.bottom = element_text(size = 5, angle = -90))


# We could take the top say 30 variables
vars_auc = 
  d_auc |> 
  dplyr::select(variable,var_metrics) |>
  unnest(var_metrics) |>
  ungroup() |>
  slice_max(order_by = mean,n = 30) |>
  mutate(type = "single_LR_auc")



## Bayesian regression =====

#We can do two things - fit a model with rstanarm, and use the projpred thing

#We have to impute the dataset to make this work

d_bayes = 
  d_0_1 |>
  analysis() 

rec_bayes = 
  d_bayes |>
  recipe(group ~ .) %>%
  step_impute_bag(all_predictors()) %>%
  prep(training = d_bayes)

d_bayes = 
  rec_bayes |>
  juice()


set.seed(1234)



#Make our horseshoe prior
n <- nrow(d_bayes)
D <- ncol(d_bayes)-1

#How many variables do you think will be useful predictors? Lets say 30 
p0 <- 30

tau0 <- p0/(D - p0) * 1/sqrt(n)
prior_coeff <- hs(global_scale = tau0, slab_scale = 1)


#Right well, lets fit it
options(mc.cores = 12)

fit_1 = 
  stan_glm(group ~ ., 
           data = d_bayes,
           family = binomial(link = "logit"),
           prior           = prior_coeff, 
           prior_intercept = normal(0, 1.5),
           cores = 12, seed = 12345,
           chains = 4, iter = 4000, warmup = 1000,
           refresh = 1,
           adapt_delta = 0.999)


#Plot our top 30 predictors

left_join(
  fit_1 |>
    coef() |>
    as_tibble(rownames = "variable") |>
    rename(estimate = value),
  fit_1 |>
    posterior_interval(prob = 0.90) |>
    as_tibble(rownames = "variable"),
  by = "variable") |>
  filter(variable != "(Intercept)") |>
  slice_max(abs(estimate), n = 30) |>
  mutate(
    variable = fct_reorder(variable, abs(estimate))
  ) |>
  ggplot(aes(x = estimate, y = variable, xmin = `5%`,xmax = `95%`)) +
  geom_errorbar() +
  geom_point() +
  labs(y = NULL) +
  theme_bw() +
  geom_vline(xintercept = 0,lty = 3)


#We can save the variables we extracted thus
vars_bayes_horseshoe = 
  left_join(
    fit_1 |>
      coef() |>
      as_tibble(rownames = "variable") |>
      rename(estimate = value),
    fit_1 |>
      posterior_interval(prob = 0.90) |>
      as_tibble(rownames = "variable"),
    by = "variable") |>
  filter(variable != "(Intercept)") |>
  slice_max(abs(estimate), n = 30) |>
  mutate(
    variable = fct_reorder(variable, abs(estimate))
  ) |>
  mutate(type = "Horseshoe")



#Do the cross validation variable selection thing - this is very slow but doesn't use
#a lot of CPU power, which is very frustrating
cvvs <- 
  fit_1 |>
  cv_varsel(
    # cv_method = "LOO",
    cv_method = "kfold",
    k = 10,
    validate_search = TRUE,
    nclusters_pred = 30,
    nterms_max = 30,
    seed = 411183
  )



#Plot our results
plot(cvvs, stats = c("auc","acc"), deltas = TRUE)

#How many terms should we take?


#The package itself suggests we take a certain number
modsize_decided <- suggest_size(cvvs)


( soltrms <- solution_terms(cvvs) )
( soltrms_final <- head(soltrms, modsize_decided) )


#Now we can project the reference model onto the final submodel
prj <- 
  project(
    fit_1,
    solution_terms = soltrms_final,
    seed = 15705533
  )


# Get the projected posterior draws
prj_mat <- as.matrix(prj)


#Then process them
prj_drws <- as_draws_matrix(prj_mat)

# In the following call, as.data.frame() is used only because pkgdown
# versions > 1.6.1 don't print the tibble correctly.

prj_drws |>
  posterior::summarize_draws("median", "mad", 
                  function(x) quantile(x, probs = c(0.025, 0.975))) |>
  filter(variable != "(Intercept)") 



#Plot marginal posteriors
mcmc_intervals(prj_mat) +
  ggplot2::coord_cartesian(xlim = c(-1.5, 1.6))

# Do an equivalent of pp_check
prj_predict <- proj_predict(prj, .seed = 762805)


# We can save the variables we identified here

vars_bayes_projpred = 
  prj_drws |>
  posterior::summarize_draws("median", "mad", 
                  function(x) quantile(x, probs = c(0.025, 0.975))) |>
  filter(variable != "(Intercept)") |>
  mutate(type = "projpred")



  ## LASSO =====
  
  
  # Next we tackle the lasso
  #Make a recipe which includes an imputation step (lets not do that for the time being)
  lasso_rec =
    recipe(group ~ .) %>%
    step_impute_bag(all_predictors())
  

  #Specify a basic lasso model
  lasso_spec = 
    logistic_reg(mode = "classification",
                 engine = "glmnet",
                 penalty = tune(), 
                 mixture = 1)
  
  #Make a workflow
  lasso_wf = 
    workflow() %>%
    add_model(lasso_spec) %>%
    add_recipe(lasso_rec)
  

  #Can we speed things up with some finetune?
  ctrl <- control_sim_anneal(verbose = TRUE)
  
  
  #We can resample using our inner sample bootstraps
  lasso_sa = 
    lasso_wf %>% 
    tune_sim_anneal(resamples = d_folds$inner_resamples[[1]], 
                    iter = 20, control = ctrl, 
                    metrics = metric_set(roc_auc))
  
  
  
  #Make a plot of our grid search
  lasso_sa %>%
    collect_metrics() %>%
    ggplot(aes(x = penalty,y = mean,ymin = mean-std_err,ymax = mean+std_err)) +
    geom_point()+
    geom_errorbar()+
    geom_line() +
    facet_wrap(~.metric,ncol = 1) +
    theme_bw() +
    scale_x_log10()
  
  #Get the best penalty term for auc
  highest_auc <- 
    lasso_sa %>%
    select_best("roc_auc")
  
  #Finalise our model
  final_lasso <- 
    finalize_workflow(
      lasso_wf %>% 
        add_model(tune_spec),
      highest_auc)
  
  #Fit the model with the best fit
  lasso_final = 
    final_lasso %>%
    fit(DF1)
  
  ## Get variable importance from our LASSO model ====
  
  #Get variable importance
  lasso_vi = 
    lasso_final %>%
    extract_fit_parsnip() %>%
    vi(lambda = highest_auc$penalty) %>%
    mutate(
      Importance = abs(Importance),
      Variable = fct_reorder(Variable, Importance)
    ) 
  
  #Plot our important variables by importance
  lasso_vi |> 
    filter(Importance > 0) |>
    ggplot(aes(x = Importance, y = Variable, fill = Sign)) +
    geom_col() +
    scale_x_continuous(expand = c(0, 0)) +
    labs(y = NULL)
  
  #This is the auc for the whole set - obvious this will be overfit as we don't have a test and training split
  lasso_final %>%
    augment(rec_2_prep ) |>
    select(group,.pred_class,.pred_Control) |>
    roc_auc(group,.pred_Control)
  

  ## Tabnet ======
  
  
  
  
  # Faff with tabnet
  # 
  # #Now it seems that we can't get tabnet to work with the finalize_workflow function
  # 
  # #Lets try and fit tabnet on its own, without using a workflowset to
  # #we can tune better
  # 
  # 
  # tn_mod <- 
  #   tabnet(batch_size = 1024,
  #          epochs = tune(), 
  #          decision_width = tune(), 
  #          attention_width = tune(),
  #          num_steps = tune(), 
  #          penalty = tune(), 
  #          learn_rate = tune()) %>%
  #   set_engine("torch", verbose = TRUE) %>%
  #   set_mode("classification")
  # 
  # tn_wf <- 
  #   workflow() %>%
  #   add_model(tn_mod) %>%
  #   add_recipe(rec_simple )
  # 
  # #Should we also tune on epochs?
  # 
  # tn_grid <-
  #   tn_wf %>%
  #   extract_parameter_set_dials() %>%
  #   update(
  #     epochs          = epochs(range = c(10,100)),
  #     decision_width  = decision_width(range = c(20, 40)),
  #     attention_width = attention_width(range = c(20, 40)),
  #     num_steps       = num_steps(range = c(4, 6)),
  #     learn_rate      = learn_rate(range = c(-2.5, -1))
  #   ) %>%
  #   grid_max_entropy(size = 20)
  # 
  # 
  # 
  # tn_ctrl <- control_race(verbose_elim = TRUE)
  # set.seed(777)
  # 
  # tn_res <- 
  #   tn_wf %>% 
  #   tune_race_anova(
  #     resamples = d_inner, 
  #     grid = tn_grid,
  #     control = tn_ctrl
  #   )
  # 
  # #Get the best params
  # show_best(tn_res, metric = "roc_auc")
  # 
  # #So this isn't working
  # tn_final = 
  #   tn_wf |>
  #   finalize_workflow(select_best(tn_res, metric = "roc_auc") )
  # 
  # #Nor is this
  # tn_mod_final =
  #   update(tn_mod,
  #          parameters = select_best(tn_res, metric = "roc_auc") )
  # 
  # 
  # #I suppose we can manually implant the parameters
  # tn_best_params = 
  #   select_best(tn_res, metric = "roc_auc") 
  # 
  # tn_mod_final <- 
  #   tabnet(batch_size = 1024, 
  #          decision_width  = !!tn_best_params$decision_width, 
  #          attention_width = !!tn_best_params$attention_width,
  #          num_steps       = !!tn_best_params$num_steps, 
  #          penalty         = !!tn_best_params$penalty, 
  #          learn_rate      = !!tn_best_params$learn_rate,
  #          epochs          = !!tn_best_params$epochs) %>%
  #   set_engine("torch", verbose = TRUE) %>%
  #   set_mode("classification")
  # 
  # tn_wf_final <- 
  #   workflow() %>%
  #   add_model(tn_mod_final) %>%
  #   add_recipe(rec_simple )
  # 
  # tn_fit_final = 
  #   tn_wf_final |>
  #   fit(data  = d_impute)
  # 
  # tn_pred = 
  #   d_outer_test |>
  #   bind_cols(
  #     predict(tn_fit_final, d_outer_test, type = "prob") ) |>
  #   roc_auc(group, .pred_Control)
  
