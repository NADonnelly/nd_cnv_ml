
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
               finetune,probably,doFuture,furrr,
               rstanarm,ggforce,doParallel,patchwork,
               DALEX,DALEXtra,readxl,lubridate, crayon)

tidymodels_prefer()



#Load the full dataset
DF = read_csv("FilteredData2.csv")

`%nin%` = Negate(`%in%`)


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

n_outer    = 20
n_inner    = 20
n_fold     = 5


# Prepare our initial data splits
d_split = initial_split(DF1, strata = group,prop = split_prop)

#Get test and training data
d_train = training(d_split)
d_test  = testing(d_split)

#Get the IDs of the participants in the assessment.test split
d_test_id = d_split |> complement()

#Prepare for (nested) cross validation
d_folds = nested_cv(d_train, 
                    outside = vfold_cv(v = n_fold, repeats = n_outer/n_fold, strata = group), 
                    inside  = bootstraps(times = n_inner, strata = group)) 

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
#I know loops are bad etc but this is a good way to get going quickly. Plus I'm used to MATLAB

#Preallocate an output tibble
d_fold_results = 
  d_folds |>
  select(id,id2) |>
  mutate(.results       = vector(mode = "list",length = nrow(d_folds)),
         .mod_posterior = vector(mode = "list",length = nrow(d_folds)))


fancy  <- combine_styles(make_style("ivory"), 
                         make_style("grey20", bg = TRUE))


#Loop through each outer fold
for(i in 1:nrow(d_folds)){
  
  #Tell the user where we are in our looping
  cat(green(
    'This is loop ' %+%
      blue$underline$bold(i) %+%
      ' of ' %+%
      blue$underline$bold(nrow(d_folds)) %+% '\n'
  ))
  
  #Get the outer data for this iteration of the loop
  d_outer = d_folds$splits[[i]]
  
  #Get the inner data for this iteration of the loop
  d_inner = d_folds$inner_resamples[[i]]
  
  
  
  #Make pre-processing recipes
  cat(fancy("Pre-processing..."), "\n")
  

  #Set up a recipe with all the variables all at once, with imputation
  rec_impute =
    d_outer |>
    analysis() |>
    recipe(group ~ .) |>
    step_zv(all_predictors()) |>
    step_impute_bag(all_predictors())   
  
  #Set up a recipe without any imputation
  rec_simple = 
    d_outer |>
    analysis() |>
    recipe(group ~ .) |>
    step_zv(all_predictors())
  
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
  
  #Set a random number seed for reproducibility
  set.seed(09062022)
  
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
  
  #We now rebuild a split object with the imputed data
  d_outer_split = 
    make_splits(d_impute,d_outer_test)
  
  
  # #Given we have an imbalance between the number of cases and controls, should we be doing something about this?
  # rec_themis = 
  #   d_impute |>
  #   recipe(group ~ .) |>
  #   step_zv(all_predictors()) |>
  #   step_smote(group,over_ratio = 1)
  # 
  
  # rec_themis |>
  #   prep() |>
  #   bake(new_data = NULL) |>
  #   count(group)

  
  #Set up a bunch of potential models. Agreed with Marianne and Valentina:

  # - Penalised Regression
  # - Linear SVM
  # - Random Forests
  # - Neural Network (with nnet for speed)
  
  
  cat(fancy("Fitting Models..."), "\n")
  
  
  #Elastic Net Regression with glmnet
  logistic_reg_spec <- 
    logistic_reg(mode = "classification",
                 penalty = tune(), 
                 mixture = tune()) %>% 
    set_engine("glmnet")
  
  
  #Linear SVM
  svm_l_spec = 
    svm_linear(cost   = tune(),
               margin = tune()) %>%
    set_engine("kernlab") %>% 
    set_mode("classification")
  
  
  #Neural Network with nnet
  nnet_spec <- 
    mlp(hidden_units = tune(), 
        penalty      = tune(), 
        epochs       = tune()) %>% 
    set_engine("nnet", MaxNWts = 2600) %>% 
    set_mode("classification")
  
  
  #Random Forest
  rf_spec <-
    rand_forest(mtry = tune(), min_n = tune(), trees = 1000) %>%
    set_engine("ranger") %>%
    set_mode("classification")

  
  #Make workflow sets
  
  
  #Some models will run in parallel happily
  wf_par <- 
    workflow_set(
      preproc = list(impute       = rec_simple),
      models  = list(lr.en        = logistic_reg_spec,
                     SVM.linear   = svm_l_spec,
                     nnet.n       = nnet_spec,
                     rf           = rf_spec))
  
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
  plan(multisession, workers = 16)
  
  #Fit our workflows that work in parallel
  wf_results <-
    wf_par |>
    workflow_map(
      "tune_race_anova",
      seed = 1503,
      resamples = d_inner,
      control = race_ctrl,
      grid = 30,
      metrics = metric_set(roc_auc,mn_log_loss)
    )
  
  
  #Turn off parallel mode
  plan(sequential)


  # #Rank models by their resampling results
  rankings <-
    wf_results |>
    rank_results(select_best = TRUE)
  
  # autoplot(wf_results)
  # autoplot(wf_results, id ="impute_lr.en")
  
  
  cat(fancy("Calculating Model Performance..."), "\n")


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
  
  
  #Extract posterior credible intervals for the models
  roc_mod_post = 
    roc_mod |>
    tidy() |>
    separate(model,into = c("vars","mod_type"),sep = "_", extra = "merge") |>
    mutate(mod_type = map_chr(mod_type,str_remove,pattern = "vars_")) |>
    group_by(vars,mod_type) |>
    tidybayes::median_hdi() |>
    arrange(-posterior)

  #Make our own plot
  # roc_mod_post |>
  #   ggplot(aes(x = forcats::fct_reorder(mod_type,posterior,max),
  #              y = posterior,ymin = .lower,ymax = .upper,group = vars,
  #              colour = mod_type)) +
  #   geom_errorbar(position = position_dodge(width = 0.5)) +
  #   geom_point(aes(shape = vars),position = position_dodge(width = 0.5)) +
  #   theme_bw() +
  #   theme(axis.text.x.bottom = element_text(size = 8, angle = -45)) +
  #   labs(y = "Model AUC", x = "Model Type",shape = "Variable \nSet")
  # 
  # #Do the ROPE Plot
  # autoplot(roc_mod, type = "ROPE", size = 0.025)

  
  
  #Get the settings for the  best version of every model, fit to
  #the outer fold training data, and then test on the outer folder test
  #data
  best_results_wf = 
    tibble(wflow_id = wf_results$wflow_id) |>
    
    #Here we extract the workflows and the best parameters in each workflow
    mutate(res         = purrr::map(wflow_id,~workflowsets::extract_workflow_set_result(wf_results,id = .x)))  |>
    mutate(best_params = purrr::map(res, tune::select_best,metric = "roc_auc")) |>
    mutate(best_wf     = purrr::map(wflow_id, ~workflowsets::extract_workflow(wf_results,id = .x))) 
  
  
  best_results_wf = 
    best_results_wf |>
    
    #We finalize the workflow by combining the best paramters with the workflows
    mutate(best_wf_final  = map2(best_wf,best_params,~finalize_workflow(.x,.y))) |>

    #Fit the best models to the full training data for this outer fold and apply to the test 
    #data for this outer fold (not previously seen by any part of the fitting process)
    mutate(outer_full_fit = purrr::map(best_wf_final, last_fit,split = d_outer_split)) |>
    
    #Put those metrics somewhere easy to grab
    mutate(best_roc_auc = map_dbl(outer_full_fit,~.x |> 
                                    select(.metrics) |> 
                                    unnest(.metrics) |> 
                                    filter(.metric == "roc_auc") |> 
                                    pull(.estimate)))
      
  
  #Store the results for this fold
  d_fold_results$.results[[i]] =
    best_results_wf |>
    dplyr::select(-c(res,best_wf)) 
  
  #Store the model training performance data too
  d_fold_results$.mod_posterior[[i]] = 
    roc_mod_post
  

  # #Tabulate our metrics
  # best_results_wf |>
  #   select(wflow_id,best_roc_auc) |>
  #   arrange(-best_roc_auc) |>
  #   knitr::kable()
  # 
  # #Or plot
  # best_results_wf |>
  #   select(wflow_id,best_roc_auc) |>
  #   mutate(wflow_id = map_chr(wflow_id,str_remove,"impute_")) |>
  #   ggplot(aes(x = forcats::fct_reorder(wflow_id,best_roc_auc,max),
  #              y = best_roc_auc)) +
  #   geom_col() +
  #   theme_bw() +
  #   theme(axis.text.x.bottom = element_text(angle = 0)) +
  #   labs(x = "Model", y = "Outer ROC") +
  #   coord_cartesian(ylim = c(0.9,1))
  
  
  #lets do some cleaning up
  rm(list = c("d_impute","d_inner", "d_outer_test", "d_outer",
              "rec_impute","rec_simple",
              "wf_par","wf_seq","wf_results","wf_seq_results","wf_par_results",
              "best_results_wf",
              "roc_mod","roc_mod_post"))  
  
}

## Summarise model performance =====

#Save this (it takes an hour or so to fit!)
write_rds(d_fold_results,"nested_cv_result_all_vars.rds")

d_fold_results = read_rds("nested_cv_result_all_vars.rds")


### Plot results ===== 

#Lets plot some results

#Plot the posterior for all submodels
d_fold_results |> 
  select(-.results) |>
  unnest(.mod_posterior) |>
  ggplot(aes(x = forcats::fct_reorder(mod_type,posterior,max),
             y = posterior,
             ymin = .lower,ymax = .upper,
             colour = mod_type)) +
  geom_point(position = position_jitter(seed = 123,width = 0.1)) +
  geom_linerange(position = position_jitter(seed = 123,width = 0.1)) +
  geom_point(data  =  
               d_fold_results |> 
               select(-.results) |>
               unnest(.mod_posterior) |>
               group_by(mod_type) |>
               summarise(sd        = sd(posterior),
                         posterior = mean(posterior)) |>
               mutate(.lower = posterior-sd, .upper= posterior+sd),
             size = 4)+
  geom_linerange(data  =  
               d_fold_results |> 
               select(-.results) |>
               unnest(.mod_posterior) |>
               group_by(mod_type) |>
               summarise(sd        = sd(posterior),
                         posterior = mean(posterior)) |>
               mutate(.lower = posterior-sd, .upper= posterior+sd),
             size = 2,
             alpha = 0.6,
             lty = 1)+
  geom_hline(yintercept = 1,lty = 2) +
  theme_bw() +
  theme(axis.text.x.bottom = element_text(size = 8, angle = -45)) +
  labs(y = "Model AUC", x = "Model Type",shape = "Variable \nSet") +
  coord_cartesian(ylim = c(.85,1))


#Look at the averaged performance on the outer fold test data
d_fold_results |> 
  select(-.results) |>
  unnest(.mod_posterior) |>
  group_by(mod_type) |>
  summarise(mean_p = mean(posterior)) |>
  mutate(mod_type = case_when(mod_type == "lr.en" ~ "Penalised LR",
                              mod_type == "nnet.n" ~ "ANN",
                              mod_type == "rf" ~ "Random Forest",
                              mod_type == "SVM.linear" ~ "Linear SVM")) |>
  arrange(-mean_p)|>
  rename(Model = mod_type,
         Performance = mean_p) |>
  mutate(Performance = round(Performance,digits = 3)) |>
  knitr::kable(format = "html", booktabs = TRUE) |>
  kableExtra::kable_styling(font_size = 11)

### Bayesian GLM of model performance ##########


#Do something like the perf_mod function in tidymodels
cv_results = 
  d_fold_results |>
  select(-.mod_posterior) |>
  unnest(.results) |>
  select(id,id2,wflow_id,best_roc_auc) |> 
  mutate(id = interaction(id,id2,sep = "_")) |>
  select(-id2)



#We can directly compare with tidyposterior, which takes the best submodel from
#each workflow and then does a stan_glm on the performance metrics
cv_mod <- 
  stan_glmer(best_roc_auc ~ 0 + wflow_id + (1|id),
              data = cv_results,
              seed = 1, refresh = 0,
              chains = 4, cores = 8,
              iter = 10000, warmup = 2000,
              prior = rstanarm::normal(0,1),
              prior_intercept = rstanarm::normal(0,1.5),
              prior_aux = rstanarm::exponential(rate = 1))


#Make our own plot
cv_mod_post = 
  cv_mod |>
  broom.mixed::tidy( conf.int = T, conf.level = 0.95, conf.method = "HPDinterval") |>
  mutate(term = str_remove(term,"wflow_id")) |>
  separate(term,into = c("vars","mod_type"),sep = "_", extra = "merge") 


p_all_vars = 
  cv_mod_post |>
  mutate(mod_type = case_when(mod_type == "lr.en" ~ "Penalised LR",
                              mod_type == "nnet.n" ~ "ANN",
                              mod_type == "rf" ~ "Random Forest",
                              mod_type == "SVM.linear" ~ "Linear SVM")) |>
  ggplot(aes(x = forcats::fct_reorder(mod_type,estimate,max),
             y = estimate,ymin =  conf.low,ymax = conf.high)) +
  # geom_errorbar(position = position_dodge(width = 0.5)) +
  geom_linerange(size = 1.5,colour = "#7D3232") +
  geom_point(size = 3,colour = "#7D3232") +
  theme_bw() +
  theme(axis.text.x.bottom = element_text(size = 8, angle = -45),
        panel.grid = element_blank()) +
  labs(y = "Model AUC", x = "Model Type",shape = "Variable \nSet") +
  coord_cartesian(ylim = c(0.8,1)) +
  geom_hline(yintercept = 1,lty = 3) 


#This might become part of something
p_all_vars

#Tabulate
tab_post = 
  cv_mod_post |>
  mutate(mod_type = case_when(mod_type == "lr.en" ~ "Penalised LR",
                              mod_type == "nnet.n" ~ "ANN",
                              mod_type == "rf" ~ "Random Forest",
                              mod_type == "SVM.linear" ~ "Linear SVM")) |>
  mutate(across(where(is.double),round,digits = 3)) |>
  arrange(-estimate) |>
  transmute(Model = mod_type,
            Performance = paste(estimate, "[",conf.low, ",",conf.high,"]",sep = " "))

tab_post |>
  knitr::kable(format = "html", booktabs = TRUE) |>
  kableExtra::kable_styling(font_size = 11)


#You can do model comparison thus:
cv_mod_rel <- 
  cv_results |>
  mutate(mod_type = str_remove(wflow_id,"impute_")) |>
  mutate(mod_type = case_when(mod_type == "lr.en" ~ "Penalised LR",
                              mod_type == "nnet.n" ~ "ANN",
                              mod_type == "rf" ~ "Random Forest",
                              mod_type == "SVM.linear" ~ "Linear SVM")) |>
  mutate(mod_type = factor(mod_type,
                           levels = tab_post$Model)) %>%
  stan_glmer(best_roc_auc ~ 1 + mod_type + (1|id),
             data = .,
             seed = 1, refresh = 0,
             chains = 4, cores = 8,
             iter = 10000, warmup = 2000,
             prior = rstanarm::normal(0,1),
             prior_intercept = rstanarm::normal(0.9,0.3),
             prior_aux = rstanarm::exponential(rate = 1))


tab_mod_rel = 
  cv_mod_rel |>
  bayestestR::describe_posterior(ci_method = "HDI") |>
  as_tibble() |>
  select(Parameter,Median, CI_low,CI_high,pd) |>
  mutate(Parameter = if_else(Parameter == "(Intercept)",tab_post$Model[1],str_remove(Parameter,"mod_type"))) |>
  mutate(across(where(is.double),round,digits = 3)) |>
  transmute(Model = Parameter,
            Difference = paste(Median, "[",CI_low, ",",CI_high,"]",sep = " "),
            `Probability of Direction` = bayestestR::pd_to_p(pd))

tab_mod_rel$Difference[1] = "-"
tab_mod_rel$`Probability of Direction`[1] = 1


### Table 2 ========

#This makes a table of both performance and model comparison
left_join(tab_post,tab_mod_rel, by  ="Model")  |>
  knitr::kable(format = "html", booktabs = TRUE) |>
  kableExtra::kable_styling(font_size = 11)

#So the SVM has the best overall performance, although the LR and RF models are very similar;
#however the ANN performance is worse (although in absolute terms still very good)


# At this point we could fit the best performing models to the full training data and test on the 
#previous held out testing data. However, as we are going to do variable selection and fit models 
#with reduced sets of variables, I think we stop now and look at variable importance, and then we
# return to fit to the full data and test in the held out test data at the very end.


## Variable importance ====


### Make full permutated dataset ======

#We can use permutation based variable importance. We might need to have a full dataset for this, so 
#we need to impute it. Lets do that

#Set up a recipe with all the variables all at once, with imputation
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
d_train_impute = 
  rec_impute_prep |>
  juice()


#And make a test set using the same imputation model, but the test (assessment) 
#data
d_test_impute = 
  rec_impute_prep |>
  bake(new_data = d_test)

#And combine
#We now rebuild a split object with the imputed data
d_last_fit = 
  make_splits(d_train_impute,d_test_impute)

#Lets save this so we can use it again
write_rds(d_last_fit,"nested_cv_imputed_data.rds")
d_last_fit = read_rds("nested_cv_imputed_data.rds")


#Select only the training data for this process
d_pred = d_last_fit |> training()

### Define a function for determining variable importance =====

get_variable_importance = function(cv_fold_results,cv_data_set,model_name,np = 100,nvar = 30){
  
  
  cat(bgWhite$blue$inverse$bold('Processing ' %+% model_name %+% '... \n'))
  
  var_final_importance = 
    cv_fold_results |>
    unnest(.results) |>
    filter(str_detect(wflow_id,model_name)) |>
    slice_max(best_roc_auc) |>
    mutate(fit_engine = map(outer_full_fit,~extract_fit_engine(.x) ))
  
  
  model_explainer = 
    explain_tidymodels(model = var_final_importance$outer_full_fit[[1]]$.workflow[[1]],
                       data = cv_data_set |> select(-group), 
                       y    = cv_data_set$group == "ND-CNV")
  
  model_parts = model_parts(explainer = model_explainer,
                            type = "variable_importance",
                            B = np)
  
  best_vars = 
    model_parts |>
    as_tibble() |> 
    filter(variable %nin% c("_baseline_","_full_model_")) |>
    group_by(variable) |>
    ggdist::median_hdci(dropout_loss,.width = 0.95) |>
    filter(dropout_loss > model_parts |>
             as_tibble() |> 
             filter(variable %in% c("_full_model_")) |>
             summarise(mdl = mean(dropout_loss)) |>
             pull()) |>
    slice_max(dropout_loss,n = nvar,with_ties = F)
  
  return(list("permute_data" = model_parts,
              "variable_summary" = best_vars))
  
}


### Calculate variable importance ========

#Set a random seed for reproducibility
set.seed(05082022)

#Set up parallel computing
registerDoFuture()
plan(multisession, workers = 16)


best_vars = 
  tibble(model_name   = c("lr.en",
                          "rf",
                          "SVM",
                          "nnet"),
         title_string = c("Penalised LR",
                          "Random Forest",
                          "Linear SVM",
                          "ANN")) |>
  mutate(var_set    = map(model_name, ~get_variable_importance(cv_fold_results = d_fold_results,
                                                               model_name      = .x,
                                                               cv_data_set     = d_pred,
                                                               np = 500,nvar = 30)))


best_vars = 
  best_vars |>
  mutate(permute_data = map(var_set, ~.x$permute_data),
         variable_summary = map(var_set, ~.x$variable_summary)) |>
  select(-var_set)


#Make a plot of the variable importance from all 4 models

shift_trans = function(d = 0) {
  scales::trans_new("shift", transform = function(x) x - d, inverse = function(x) x + d)
}

plot_parts = function(part_data,title_string){
  
  part_data |>
    as_tibble() |> 
    filter(variable %nin% c("_baseline_","_full_model_")) |>
    group_by(variable) |>
    ggdist::mean_hdci(dropout_loss,.width = 0.95) |>
    slice_max(dropout_loss,n = 30,with_ties = F) |>
    mutate(variable = fct_reorder(variable,dropout_loss,max)) |>
    ggplot(aes(y = variable,x = dropout_loss,xmin = .lower,xmax = .upper)) +
    geom_col() +
    geom_point() +
    geom_linerange() +
    theme_bw() +
    theme(panel.grid = element_blank()) +
    labs(x = "Change in AUC after permutation", y = "Variable", title = title_string) +
    scale_x_continuous(trans = shift_trans(d=  part_data |>
                                             as_tibble() |> 
                                             filter(variable %in% c("_full_model_")) |>
                                             summarise(mdl = mean(dropout_loss)) |>
                                             pull(mdl)))  +
    geom_vline(xintercept = part_data |>
                 as_tibble() |> 
                 filter(variable %in% c("_full_model_")) |>
                 summarise(mdl = mean(dropout_loss)) |>
                 pull(mdl),lty = 3)
  
}


#Now we can plot
best_vars = 
  best_vars |>
  mutate(plots = map2(permute_data,title_string,plot_parts))

wrap_plots(best_vars$plots) + plot_annotation(tag_levels = "A")


### Combine Variables =====

# Compare our variables
var_set = 
  best_vars |>
  select(model_name,variable_summary) |>
  unnest(variable_summary) |>
  rename(type = model_name, 
         Variable = variable, 
         Importance = dropout_loss) |> 
  count(Variable)


#Make a plot
var_set   |>
  ggplot(aes(y = fct_reorder(Variable,n,max), x = n)) +
  geom_col() +
  theme_bw() +
  theme(axis.text.y= element_text(size = 6, angle = -0))


#Take only variables included in more than 1 method for extraction
vars_more_than_one = 
  var_set |>
  filter(n > 1)


#take the vars with the highest number of inclusions
vars_max = 
  var_set |>
  # slice_max(order_by = n)
  filter(n > 2)

#Now we have our sets of variables, lets do our ML again with these subsets


# Variable subset ML =======

#Convert the sets of variables we selected above into formulas and then into recipes
form_list = 
  bind_rows(best_vars |>
              select(model_name,variable_summary) |>
              unnest(variable_summary) |>
              rename(type = model_name, 
                     Variable = variable, 
                     Importance = dropout_loss),
            vars_max |> 
              mutate(type = "max.vars") |>
              rename(Importance = n),
            vars_more_than_one|> 
              mutate(type = "moreThanOne.vars") |>
              rename(Importance = n)) |>
  group_by(type) |>
  nest() |>
  
  #Make our formulas
  mutate(formulas = map_chr(data,~paste("group ~",.x |> 
                                          pull(Variable) |> 
                                          paste(collapse = " + ")))) |>
  
  #Make recipes which we will append our formulas into
  mutate(recipes = map(formulas, ~recipe(formula = as.formula(.x),
                                         data    = d_pred) |>
                         step_zv(all_predictors())))

#Extract the recipes and convert them into a named list
final_recipes = 
  form_list |>
  pull(recipes) 

names(final_recipes) <-  
  form_list |> pull(type)


#Save these variable sets
write_rds(final_recipes,"nested_cv_selected_vars.rds")
final_recipes = read_rds("nested_cv_selected_vars.rds")



#Now prepare for our nested cross validation

#Set a random number seed for reproducibility
set.seed(15072022)

#Use the test-train split we did above, but prepare new cv folds from the 
#training data

#Prepare for (nested) cross validation
d_folds = nested_cv(d_train, 
                    outside = vfold_cv(v = n_fold, repeats = n_outer/n_fold), 
                    inside  = bootstraps(times = n_inner)) 


#Preallocate an output tibble
d_var_select_results = 
  d_folds |>
  select(id,id2) |>
  mutate(.results       = vector(mode = "list",length = nrow(d_folds)),
         .mod_posterior = vector(mode = "list",length = nrow(d_folds)))


#Loop through each outer fold
for(i in 1:nrow(d_folds)){
  
  cat(green(
    'This is loop ' %+%
      red$underline$bold(i) %+%
      ' of ' %+%
      red$underline$bold(nrow(d_folds)) %+% '\n'
  ))
  
  
  #Get the outer data for this iteration of the loop
  d_outer = d_folds$splits[[i]]
  
  #Get the inner data for this iteration of the loop
  d_inner = d_folds$inner_resamples[[i]]
  
  
  
  #Make pre-processing recipes
  cat(fancy("Pre-processing..."), "\n")
  
  #Set up a recipe with all the variables all at once, with imputation
  rec_impute =
    d_outer |>
    analysis() |>
    recipe(group ~ .) |>
    step_zv(all_predictors()) |>
    step_impute_bag(all_predictors())   
  
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
  
  #We now rebuild a split object with the imputed data
  d_outer_split = 
    make_splits(d_impute,d_outer_test)
  
  
  ## Set models =====
  
  
  cat(fancy("Fitting Models..."), "\n")
  
  #Elastic Net Regression with glmnet
  logistic_reg_spec <- 
    logistic_reg(mode = "classification",
                 penalty = tune(), 
                 mixture = tune()) %>% 
    set_engine("glmnet")
  
  
  #Linear SVM
  svm_l_spec = 
    svm_linear(cost   = tune(),
               margin = tune()) %>%
    set_engine("kernlab") %>% 
    set_mode("classification")
  
  
  #Neural Network with nnet
  nnet_spec <- 
    mlp(hidden_units = tune(), 
        penalty      = tune(), 
        epochs       = tune()) %>% 
    set_engine("nnet", MaxNWts = 2600) %>% 
    set_mode("classification")
  
  #Random Forest
  rf_spec <- 
    rand_forest(mtry = tune(), min_n = tune(), trees = 1000) %>% 
    set_engine("ranger") %>% 
    set_mode("classification")
  
  #Make the workflow set
  
  #These are the ones that fit in parallel
  wf_par <- 
    workflow_set(
      preproc = final_recipes, 
      models  = list(lr.en        = logistic_reg_spec,
                     SVM.linear   = svm_l_spec,
                     nnet.n       = nnet_spec,
                     rf           = rf_spec)
    )
  

## Fit models =====
  
  
  #Set up controls for the finetune anova race method
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
      resamples = d_inner,
      control = race_ctrl,
      grid = 25,
      metrics = metric_set(roc_auc,mn_log_loss)
    )
  
  
  #Turn off parallel mode
  plan(sequential)
  

  #Unite the datasets
  wf_results = 
    wf_par_results
  

  #We can directly compare with tidyposterior, which takes the best submodel from
  #each workflow and then does a stan_glm on the performance metrics
  cat(fancy("Calculating Model Performance..."), "\n")
  
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
  

  # roc_mod_post |>
  #   ggplot(aes(x = forcats::fct_reorder(mod_type,posterior,max),
  #              y = posterior,ymin = .lower,ymax = .upper,group = vars,
  #              colour = mod_type)) +
  #   geom_errorbar(position = position_dodge(width = 0.5)) +
  #   geom_point(aes(shape = vars),position = position_dodge(width = 0.5)) +
  #   theme_bw() +
  #   theme(axis.text.x.bottom = element_text(size = 8, angle = -45)) +
  #   labs(y = "Model AUC", x = "Model Type",shape = "Variable \nSet")
  # 
  # #Do the ROPE Plot
  # autoplot(roc_mod, type = "ROPE", size = 0.025)
  
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
    mutate(best_wf_final  = map2(best_wf,best_params,~finalize_workflow(.x,.y))) |>
    
    #Fit the best models to the full training data for this outer fold and apply the best model
    #to the  test data for this outer fold (not previously seen by any part of the fitting process)
    #the last_fit function is the one doing the work here
    mutate(outer_full_fit = purrr::map(best_wf_final, last_fit,split = d_outer_split,metrics = metric_set(accuracy,kap,mn_log_loss,roc_auc,gain_capture))) |>
    
    #Put those metrics somewhere easy to grab
    mutate(best_roc_auc = map_dbl(outer_full_fit,~.x |> 
                                    select(.metrics) |> 
                                    unnest(.metrics) |> 
                                    filter(.metric == "roc_auc") |> 
                                    pull(.estimate)))
  
  
  #Store the results for this fold
  d_var_select_results$.results[[i]] =
    best_results_wf |>
    dplyr::select(-c(res,best_wf)) 
  
  #Store the performance model data too
  d_var_select_results$.mod_posterior[[i]] = 
    roc_mod_post
  

  # #Tabulate our metrics
  # best_results_wf |>
  #   select(wflow_id,best_roc_auc) |>
  #   arrange(-best_roc_auc) |>
  #   knitr::kable()
  # 
  # #Or plot
  # best_results_wf |>
  #   select(wflow_id,best_roc_auc) |>
  #   mutate(wflow_id = map_chr(wflow_id,str_remove,"impute_")) |>
  #   ggplot(aes(x = forcats::fct_reorder(wflow_id,best_roc_auc,max),
  #              y = best_roc_auc)) +
  #   geom_col() +
  #   theme_bw() +
  #   theme(axis.text.x.bottom = element_text(angle = -25)) +
  #   labs(x = "Model", y = "Outer ROC")

  #lets do some cleaning up
  rm(list = c("d_impute","d_inner", "d_outer_test", "d_outer",
              "rec_impute","rec_simple",
              "wf_par","wf_seq","wf_results","wf_seq_results","wf_par_results",
              "best_results_wf",
              "roc_mod","roc_mod_post"))
  
}


#Save this (it takes an hour or so to fit!)
write_rds(d_var_select_results,"nested_cv_result_selected_vars.rds")


## Summarise model performance =====


#Load all our models
d_var_select_results = read_rds("nested_cv_result_selected_vars.rds")
d_fold_results       = read_rds("nested_cv_result_all_vars.rds")

### Make plots ######

#Lets plot some results

#Plot the posterior for all submodels
d_var_select_results |> 
  select(-.results) |>
  unnest(.mod_posterior) |>
  ggplot(aes(x = forcats::fct_reorder(vars,posterior,max),
             y = posterior,
             group = vars,
             ymin = .lower,ymax = .upper,
             colour = mod_type)) +
  geom_point(aes(shape = vars),
             position = position_jitter(seed = 123,width = 0.2)) +
  geom_linerange(position = position_jitter(seed = 123,width = 0.2)) +
  geom_point(data  =  
               d_var_select_results |> 
               select(-.results) |>
               unnest(.mod_posterior) |>
               group_by(vars,mod_type) |>
               summarise(sd        = sd(posterior),
                         posterior = mean(posterior)) |>
               mutate(.lower = posterior-sd, .upper= posterior+sd),
             aes(shape = vars),
             # position = position_jitter(seed = 321,width = 0.1),
             size = 4)+
  geom_linerange(data  =  
                   d_var_select_results |> 
                   select(-.results) |>
                   unnest(.mod_posterior) |>
                   group_by(vars,mod_type) |>
                   summarise(sd        = sd(posterior),
                             posterior = mean(posterior)) |>
                   mutate(.lower = posterior-sd, .upper= posterior+sd),
                 size = 2,
                 alpha = 0.6,
                 # position = position_jitter(seed = 321,width = 0.1),
                 lty = 1) +
  theme_bw() +
  theme(axis.text.x.bottom = element_text(size = 8, angle = -45),
        legend.position = "none") +
  labs(y = "Model AUC", x = "Variable Set",shape = "Variable \nSet") +
  coord_cartesian(ylim = c(.85,1)) +
  facet_wrap(~mod_type,ncol = 4)


#Look at the averaged performance on the outer fold test data
d_var_select_results |> 
  select(-.results) |>
  unnest(.mod_posterior) |>
  group_by(vars,mod_type) |>
  summarise(mean_p = mean(posterior)) |>
  arrange(-mean_p)|>
  knitr::kable()


### Bayesian GLM of performance ==========

#Do something like the perf_mod function in tidymodels
cv_results = 
  d_var_select_results |>
  select(-.mod_posterior) |>
  unnest(.results) |>
  select(id,id2,wflow_id,best_roc_auc) |> 
  mutate(id = interaction(id,id2,sep = "_")) |>
  # separate(wflow_id,into = c("vars","mod_type"),sep = "_", extra = "merge") |>
  select(-c(id2))


#Add the data from the all variable models
cv_results = 
  cv_results |>
  bind_rows(d_fold_results |>
              select(-.mod_posterior) |>  
              unnest(.results) |>
              select(id,id2,wflow_id,best_roc_auc) |> 
              mutate(id = interaction(id,id2,sep = "_")) |>
              select(-id2) |>
              mutate(wflow_id = str_replace(wflow_id, "impute","all.vals")))

#We can directly compare with tidyposterior, which takes the best submodel from
#each workflow and then does a stan_glm on the performance metrics
cv_mod <- 
  stan_glmer(best_roc_auc ~ 0 + wflow_id + (1|id),
             data = cv_results,
             seed = 1, refresh = 0,
             chains = 4, cores = 8,
             iter = 10000, warmup = 2000,
             prior = rstanarm::normal(0,1),
             prior_intercept = rstanarm::normal(0,1.5),
             prior_aux = rstanarm::exponential(rate = 1))

# cv_mod |> summary()


### Figure 3: Alternative Top panel ======


#Lets make this into a plot


#Make our own plot
cv_mod_post = 
  cv_mod |>
  broom.mixed::tidy( conf.int = T, conf.level = 0.95, conf.method = "HPDinterval") |>
  mutate(term = str_remove(term,"wflow_id")) |>
  separate(term,into = c("vars","mod_type"),sep = "_", extra = "merge") 


p_compare_vars = 
  cv_mod_post |>
  mutate(mod_type = case_when(mod_type == "lr.en" ~ "Penalised LR",
                              mod_type == "nnet.n" ~ "ANN",
                              mod_type == "rf" ~ "Random Forest",
                              mod_type == "SVM.linear" ~ "Linear SVM"),
         vars = case_when(vars == "lr.en" ~ "Penalised LR",
                          vars == "nnet" ~ "ANN",
                          vars == "rf" ~ "Random Forest",
                          vars == "SVM" ~ "Linear SVM",
                          vars == "all.vals" ~ "All Variables",
                          vars == "moreThanOne.vars" ~ "> 1 Model",
                          vars == "max.vars" ~ "> 2 Models")) |>
  ggplot(aes(x = forcats::fct_reorder(vars,estimate,max),
             y = estimate,ymin =  conf.low,ymax = conf.high,
             colour = mod_type,
             fill   = mod_type)) +
  geom_errorbar(position = position_dodge(width = 0.5),size = 0.5) +
  geom_point(position = position_dodge(width = 0.5),size = 1) +
  theme_bw() +
  theme(axis.text.x.bottom = element_text(size = 6, angle = -45),
        axis.text.y = element_text(size = 6),
        axis.title = element_text(size = 6),
        strip.text = element_text(size = 6),
        strip.background = element_blank(),
        legend.position = "none",
        panel.grid = element_blank()) +
  labs(y = "Training Data Median AUROC", x = "Variable Set")+
  facet_wrap(~mod_type,ncol = 7) +
  geom_hline(yintercept = 1)

p_compare_vars





### Model Comparison model =====

#So it looks like the subset of variables that were included in the top 30 variables of multiple 
#ML models do best (better even that all variables)

cv_mod_post |> 
  group_by(mod_type)  |> 
  arrange(-estimate) |> 
  mutate(across(where(is.double),round,digits = 4)) |>
  ungroup() |>
  transmute(Model = mod_type,`Variable Set` = vars, Performance = paste(estimate,"[",conf.low,",",conf.high,"]", sep = " ")) |>
  print(n = 34)

#Select the best performing variable set for each model type
cv_mod_post |> 
  group_by(mod_type) |> 
  slice_max(estimate,n=1) |>
  mutate(across(where(is.double),round,digits = 3)) |>
  ungroup() |>
  arrange(-estimate) |>
  transmute(Model = mod_type,
            `Variable Set` = vars, 
            Performance = paste(estimate,"[",conf.low,",",conf.high,"]", sep = " ")) |>
  print(n = 34)




#You can do model comparison thus:
cv_mod_compare <- 
  cv_results |> 
  # mutate(wflow_id = relevel(factor(wflow_id),ref = "moreThanOne.vars_SVM.linear")) %>%
  mutate(wflow_id = relevel(factor(wflow_id),ref = "SVM_SVM.linear")) %>%
  stan_glmer(best_roc_auc ~ 1 + wflow_id + (1|id),
             data = .,
             seed = 1, refresh = 0,
             chains = 4, cores = 8,
             iter = 10000, warmup = 2000,
             prior = rstanarm::normal(0,0.5),
             prior_intercept = rstanarm::normal(0.9,0.5),
             prior_aux = rstanarm::exponential(rate = 1))


tab_mod_compare = 
  cv_mod_compare |>
  bayestestR::describe_posterior(ci_method = "HDI") |>
  as_tibble() |>
  select(Parameter,Median, CI_low,CI_high,pd) |>
  # mutate(Parameter = case_when(Parameter == "(Intercept)" ~ "moreThanOne.vars_SVM.linear",
  mutate(Parameter = case_when(Parameter == "(Intercept)" ~ "SVM_SVM.linear",
                               TRUE ~ Parameter)) |>
  mutate(Parameter = map_chr(Parameter,str_remove,"wflow_id")) |>
  separate(Parameter,into = c("var_set","mod_type"),sep = "_", extra = "merge") |>
           
  mutate(across(where(is.double),round,digits = 3)) |>
  transmute(Model = mod_type,
            `Variable Set` = var_set,
            Difference = paste(Median, "[",CI_low, ",",CI_high,"]",sep = " "),
            `Probability of Direction` = bayestestR::pd_to_p(pd))

tab_mod_compare$Difference[1] = "-"
tab_mod_compare$`Probability of Direction`[1] = 1


### Table 3 ======

#This is going to become table 3 so lets make it look nice
left_join(tab_mod_compare,
          cv_mod_post |> 
            mutate(across(where(is.double),round,digits = 3)) |>
            ungroup() |>
            transmute(Model = mod_type,
                      `Variable Set` = vars, 
                      Performance = paste(estimate,"[",conf.low,",",conf.high,"]", sep = " ")) , 
          by  = c("Model","Variable Set"))  |>
  relocate(Model, `Variable Set`, Performance,Difference,`Probability of Direction`) |>
  
  mutate(Model = case_when(Model == "lr.en" ~ "Penalised LR",
                           Model == "nnet.n" ~ "ANN",
                           Model == "rf" ~ "Random Forest",
                           Model == "SVM.linear" ~ "Linear SVM"),
         `Variable Set` = case_when(`Variable Set` == "all.vals" ~ "All Variables",
                                    `Variable Set` == "SVM" ~ "SVM",
                                    `Variable Set` == "lr.en" ~ "Penalised LR",
                                    `Variable Set` == "nnet" ~ "ANN",
                                    `Variable Set` == "rf" ~ "Random Forest",
                                    `Variable Set` == "max.vars" ~ "Variables selected by > 2 models",
                                    `Variable Set` == "moreThanOne.vars" ~ "Variables selected by > 1 model")) |>
  
  #Now make it nice and print it to the viewer
  knitr::kable(format = "html", booktabs = TRUE) |>
  kableExtra::kable_styling(font_size = 11)

#And then we do a bit of re-ordering of the table in excel/word to make the publication version

tab_mod_compare |>
  filter(Model == "SVM.linear")





# Fit final models to the test data ========


## Prepare imputed test data =====


#As we are interested in using the imputed dataset, we can load
#the set that we made earlier
d_last_fit = read_rds("nested_cv_imputed_data.rds")

#And load our models
d_var_select_results = read_rds("nested_cv_result_selected_vars.rds")

## Fit and evaluate models to full training data =====

# Now we fit our best models to the (imputed) training set

#Select the best set of parameters for each model - we will use the SVM variables as this appeared to have the best 
#performance (although differences between sets was very small)
best_mods = 
  d_var_select_results |>
  select(-.mod_posterior) |>
  unnest(.results) |> 
  separate(wflow_id,into = c("vars","mod_type"),sep = "_", extra = "merge") |>
  filter(vars == "SVM") |>
  group_by(mod_type) |>
  slice_max(best_roc_auc)

best_mods = 
  best_mods |>
  mutate(test_fit = purrr::map(best_wf_final, last_fit,split = d_last_fit,metrics = metric_set(accuracy,kap,mn_log_loss,roc_auc,gain_capture))) |>
  mutate(best_metrics = map(test_fit, ~.x$.metrics[[1]]))



## Performance measures =====

#Tabulate our final performance
tab_mods = 
  best_mods |>
  
  #We can work out the Brier score (this seems to simply be the average difference between the predicted 
  #probability and true outcome squared i.e. the mean squared error) - there is also the mn log loss measure
  #in yarstick which is the same thing on a slightly different scale (the log scale)
  mutate(brier = map_dbl(test_fit,~ .x$.predictions[[1]] |>
                           transmute(pred = `.pred_ND-CNV`,group) |>
                           mutate(group = ifelse(group == "ND-CNV",1,0)) |>
                           mutate(dist  = (group - pred)^2) |>
                           summarise(brier = mean(dist)) |>
                           pull())) |>
  select(mod_type,best_metrics,brier)|>
  unnest(best_metrics) |>
  select(mod_type,.metric,.estimate,brier) |>
  pivot_wider(names_from = .metric,values_from = .estimate) |>
  arrange(-roc_auc,mn_log_loss) |>
  mutate(across(where(is.double),round,digits = 3)) |>
  mutate(mod_type = case_when(mod_type == "lr.en" ~ "Penalised LR",
                              mod_type == "nnet.n" ~ "ANN",
                              mod_type == "rf" ~ "Random Forest",
                              mod_type == "SVM.linear" ~ "Linear SVM")) |>
  rename(gini_coefficient = gain_capture,
         model = mod_type) |>
  select(model,roc_auc,kap,mn_log_loss) 
  # relocate(model,roc_auc,gini_coefficient,mn_log_loss,brier,kap,accuracy) 

tab_mods |>
  knitr::kable(format = "html", booktabs = TRUE) |>
  kableExtra::kable_styling(font_size = 11)


## Bootstrap some confidence intervals =====

#We can estimate confidence intervals for the final performance metrics by bootstrapping
bootstrap_wf_metrics <- function(wf,data_set,n_boot = 500,c_m = metric_set(roc_auc,mn_log_loss)){
  

  #We need an input workflow and an input dataset 
  
  #We then make this into a bootstrap object from rsample
  bs = 
    data_set |>
    bootstraps(times = n_boot, strata = "group")  |>
    mutate(boot_metrics = map(splits, ~bind_cols(wf |> 
                                                   predict(.x |> analysis(),
                                                           type = "prob"),
                                                 .x |>
                                                   analysis()|>
                                                   select(group)) |>
                                mutate(group = fct_rev(group)) |>
                                c_m(group,`.pred_ND-CNV`)))
  
  return(bs)
  
}

#Now apply this function to our data
boot_mods = 
  best_mods |> 
  select(mod_type,test_fit) |> 
  unnest(test_fit)|>
  mutate(boot_metrics = map(.workflow, ~bootstrap_wf_metrics(wf = .x,
                                                            data_set = d_last_fit |>
                                                              testing(),
                                                            n_boot = 2000) ))



#Make a plot of our bootstrapped values
boot_mods |>
  select(mod_type,boot_metrics) |>
  unnest(boot_metrics)|>
  select(-splits) |>
  unnest(boot_metrics) |>
  ggplot(aes(x = .estimate, y = mod_type)) +
  stat_halfeye() +
  facet_wrap(~.metric,ncol = 1)


#We can also do model comparison using the bootstrap samples

boot_pd = left_join(boot_mods |>
                      select(mod_type,boot_metrics) |>
                      unnest(boot_metrics)|>
                      select(-splits) |>
                      unnest(boot_metrics) |>
                      filter(.metric == "roc_auc") |>
                      pivot_wider(names_from = mod_type,values_from = .estimate) |>
                      mutate(svm_lr = SVM.linear - lr.en,
                             svm_rf = SVM.linear - rf,
                             svm_nnet = SVM.linear - nnet.n) |>
                      select(svm_lr,svm_rf,svm_nnet) |>
                      pivot_longer(everything()) |>
                      group_by(name) |>
                      ggdist::mean_hdci(value),
                    boot_mods |>
                      select(mod_type,boot_metrics) |>
                      unnest(boot_metrics)|>
                      select(-splits) |>
                      unnest(boot_metrics) |>
                      filter(.metric == "roc_auc") |>
                      pivot_wider(names_from = mod_type,values_from = .estimate) |>
                      mutate(svm_lr = SVM.linear - lr.en,
                             svm_rf = SVM.linear - rf,
                             svm_nnet = SVM.linear - nnet.n) |>
                      select(svm_lr,svm_rf,svm_nnet) |>
                      pivot_longer(everything()) |>
                      group_by(name) |>
                      nest() |>
                      mutate(pd = map(data, ~p_direction(.x) |> as.numeric() |> pd_to_p(direction = "two-sided"))) |>
                      unnest(pd) |> select(-data),
                    by = "name")

boot_pd = 
  boot_pd |>
  mutate(name = factor(name,levels = c("svm_lr","svm_rf","svm_nnet"))) |>
  arrange(name) |>
  select(name,value,.lower,.upper,pd) |>
  mutate(across(where(is.double),round,digits = 3)) |>
  mutate(performance = paste(value,"[",.lower,",",.upper,"]",sep = " ")) |>
  select(name,performance,pd)


## Table 4 ======

#And we can make this into a more pleasing tabular format
tab_boot = 
  boot_mods |>
  select(mod_type,boot_metrics) |>
  unnest(boot_metrics)|>
  select(-splits) |>
  unnest(boot_metrics) |>
  group_by(mod_type,.metric) |>
  ggdist::mean_hdci(.estimate) |>
  arrange(.metric,-.estimate)


tab_boot |>
  mutate(across(where(is.double),round,digits = 3)) |>
  mutate(Performance = paste(.estimate,"[",.lower,",",.upper,"]",sep = " ")) |>
  select(mod_type,.metric,Performance) |>
  mutate(mod_type = case_when(mod_type == "lr.en" ~ "Penalised LR",
                              mod_type == "nnet.n" ~ "ANN",
                              mod_type == "rf" ~ "Random Forest",
                              mod_type == "SVM.linear" ~ "Linear SVM")) |>  
  rename(Model = mod_type) |>
  pivot_wider(names_from = .metric,values_from = Performance) |>
  mutate(Model = factor(Model, levels = c("Linear SVM","Penalised LR","Random Forest","ANN"))) |>
  arrange(Model) |>
  mutate(roc_diff = c("-",boot_pd$performance),
         roc_pd   = c("-",boot_pd$pd)) |>
  knitr::kable(format = "html", booktabs = TRUE) |>
  kableExtra::kable_styling(font_size = 11) 


## Calibration measures ======

#Performance by age

#What does performance look like with age?

## read info from our database xlsx files
DBIDs <-
  read_xlsx("MASTERDATABASE_BE_09_11_18.xlsx", sheet = "IDs")

DBIDs <-
  DBIDs |> 
  select(IDs, DOB, Gender, CAPAinterviewdateW1) |>
  mutate(DOB = ymd(DOB),
         CAPAinterviewdateW1 = ymd(CAPAinterviewdateW1),
         Age = interval(DOB,CAPAinterviewdateW1) |> 
           as.period() |>
           as.numeric("years"),
         Gender = as.numeric(Gender)) 


#We preserved the IDs of the participants in the test dataset above
tab_age = 
  bind_cols(
  DBIDs |>
    filter(IDs %in% DF$IDs) |>
    slice(d_test_id),
  best_mods$test_fit[[4]]$.predictions[[1]]) %>%
  mutate(age_binned = cut(Age,quantile(Age,0:5/5),include.lowest = T)) |>
  group_by(age_binned)|>
  mutate(group = fct_rev(group),
         .pred_class = fct_rev(.pred_class)) |>
  yardstick::metrics(truth = group,estimate=.pred_class, `.pred_ND-CNV`) |>
  select(-.estimator) |>
  pivot_wider(names_from = .metric,values_from = .estimate) |>
  mutate(across(where(is.double),round,digits = 3)) |>
  # relocate(age_binned,roc_auc,mn_log_loss,kap,accuracy)
  select(-accuracy) |>
  relocate(age_binned,roc_auc,kap,mn_log_loss)

tab_age |>
  knitr::kable(format = "html", booktabs = TRUE) |>
  kableExtra::kable_styling(font_size = 11)




#And do the same with gender
tab_gender = 
  bind_cols(
  DBIDs |>
    filter(IDs %in% DF$IDs) |>
    slice(d_test_id),
  best_mods$test_fit[[4]]$.predictions[[1]]) %>%
  mutate(Gender = case_when(
    Gender == 1 ~ "Male",
    Gender == 2 ~ "Female"  )) |>
  group_by(Gender) |>
  # yardstick::roc_auc(group,`.pred_ND-CNV`,event_level = "second")
  mutate(group = fct_rev(group),
         .pred_class = fct_rev(.pred_class)) |>
  yardstick::metrics(truth = group,estimate=.pred_class, `.pred_ND-CNV`) |>
  select(-.estimator) |>
  pivot_wider(names_from = .metric,values_from = .estimate) |>
  mutate(across(where(is.double),round,digits = 3)) |>
  relocate(Gender,roc_auc,kap,mn_log_loss,accuracy) |>
  select(-accuracy)

tab_gender |>
  knitr::kable(format = "html", booktabs = TRUE) |>
  kableExtra::kable_styling(font_size = 11)


### Table 5 ======

#Tabulate our splits
bind_rows(
  tab_age |>
    mutate(Covariate = "Age",
           age_binned = as.character(age_binned)) |>
    rename(`Covariate Value` = age_binned),
  tab_gender |>
    mutate(Covariate = "Gender") |>
    rename(`Covariate Value` = Gender)) |>
  relocate(Covariate, `Covariate Value`)  |>
  knitr::kable(format = "html", booktabs = TRUE) |>
  kableExtra::kable_styling(font_size = 11)



## Prepare plots =====


#Plot this stuff
p_mods =
  tab_mods |>
  pivot_longer(-model,names_to = ".metric",values_to = ".estimate") |>
  # mutate(.metric = factor(.metric, levels = c("roc_auc","gini_coefficient","accuracy","kap","mn_log_loss","brier"))) |>
  mutate(.metric = factor(.metric, levels = c("roc_auc","kap","mn_log_loss"))) |>
  ggplot(aes(x = model, y = .estimate)) +
  geom_point() +
  facet_wrap(~.metric,nrow = 1) +
  geom_segment(yend = 0,aes(xend = model)) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text.x.bottom = element_text(angle = -45)) +
  labs(x = "Model")


p_age = 
  tab_age |>
  pivot_longer(-age_binned,names_to = ".metric",values_to = ".estimate") |>
  mutate(age_m = map_chr(age_binned, ~as.character(.x) %>% gsub("\\)|\\]|\\(|\\[", "", .))) |>
  separate(age_m,sep = ",",into = c("l","u"),convert = T) |>
  mutate(mid = (l+u)/2) |>
  select(mid,.metric,.estimate) |>
  mutate(.metric = factor(.metric, levels = c("roc_auc","kap","mn_log_loss"))) |>
  ggplot(aes(x = mid, y = .estimate)) +
  geom_point() +
  facet_wrap(~.metric,nrow = 1) +
  geom_segment(yend = 0,aes(xend = mid)) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        strip.background = element_blank()) +
  labs(x = "Age")


p_gender = 
  tab_gender |>
  pivot_longer(-Gender,names_to = ".metric",values_to = ".estimate") |>
  mutate(.metric = factor(.metric, levels = c("roc_auc","kap","mn_log_loss"))) |>
  ggplot(aes(x = Gender, y = .estimate)) +
  geom_point() +
  facet_wrap(~.metric,nrow = 1) +
  geom_segment(yend = 0,aes(xend = Gender)) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        strip.background = element_blank()) +
  labs(x = "Gender")

p_mods/p_age/p_gender + plot_annotation(tag_levels = "A")



#Model calibration

best_mods = 
  best_mods |>
  mutate(calPlotData = map(test_fit, ~.x$.predictions[[1]] %>%
                             transmute(pred = `.pred_ND-CNV`,group) %>%
                             mutate(group = fct_rev(group)) %>%
                             caret::calibration(group ~ pred, data = .,cuts = 10) %>%
                             .$data %>%
                             as_tibble()))

p_cal = 
  best_mods |> 
  select(mod_type,calPlotData) |> 
  # filter(mod_type != "cf") |>
  filter(mod_type == "SVM.linear") |>
  mutate(mod_type = case_when(mod_type == "lr.en" ~ "Penalised LR",
                              mod_type == "nnet.n" ~ "ANN",
                              mod_type == "rf" ~ "Random Forest",
                              mod_type == "SVM.linear" ~ "Linear SVM")) |>
  unnest(calPlotData)  |>
  ggplot(aes(x = midpoint,y = Percent,ymin = Lower,ymax = Upper)) +
  geom_point(size = 0.4) +
  geom_linerange(size = .25) +
  geom_abline(intercept = 0,slope = 1, lty =  2,size = 0.25) +
  geom_smooth(formula = y ~ x, method = "lm",colour = "black",size = 0.25,lty=1,alpha = 0.25) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        legend.title = element_blank(),
        legend.position = "top",
        legend.text = element_text(size = 6),
        axis.title = element_text(size = 6),
        axis.text = element_text(size = 6),
        strip.text = element_text(size = 6)) +
  scale_x_continuous(breaks = seq(from = 0, to = 100, by = 25)) +
  labs(x = "Predicted Event %", y = "Observed Event %") +
  coord_cartesian(ylim = c(0,100), xlim = c(0,100)) +
  facet_wrap(~mod_type,ncol =1, strip.position="left")



#Full ROC curves
p_roc = 
  best_mods |>
  mutate(roc_curve_d = map(test_fit,~ .x$.predictions[[1]] |>
                           yardstick::roc_curve(group,`.pred_ND-CNV`,event_level = "second"))) |>
  select(mod_type,roc_curve_d) |>
  unnest(roc_curve_d) |>
  mutate(mod_type = case_when(mod_type == "lr.en" ~ "Penalised LR",
                              mod_type == "nnet.n" ~ "ANN",
                              mod_type == "rf" ~ "Random Forest",
                              mod_type == "SVM.linear" ~ "Linear SVM")) |>
  ggplot(aes(x = 1 - specificity, y = sensitivity,colour = mod_type)) +
  geom_path(size = 0.25) +
  geom_abline(lty = 3) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.title = element_blank(),
        legend.position = "bottom",
        legend.text = element_text(size = 6),
        axis.title = element_text(size = 6),
        axis.text = element_text(size = 6))







# Thresholds
threshold_data = 
  best_mods$test_fit[[4]]$.predictions[[1]] |>
  mutate(group = fct_rev(group)) |>
  probably::threshold_perf(group,`.pred_ND-CNV`,
                           thresholds = seq(0, 1, by = 0.0025))


threshold_data <- 
  threshold_data|>
  filter(.metric != "distance") |>
  mutate(group = case_when(
    .metric == "sens" | .metric == "spec" ~ "1",
    TRUE ~ "2"
  ))


max_j_index_threshold <- 
  threshold_data %>%
  filter(.metric == "j_index") |>
  filter(.estimate == max(.estimate)) |>
  pull(.threshold) |>
  median()

p_thres = 
  ggplot(threshold_data, aes(x = .threshold, y = .estimate, color = .metric)) +
  geom_path(size = 0.25) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.title = element_blank(),
        legend.position = "top",
        legend.text = element_text(size = 6),
        axis.title = element_text(size = 6),
        axis.text = element_text(size = 6)) +
  scale_alpha_manual(values = c(.4, 1), guide = "none") +
  geom_vline(xintercept = max_j_index_threshold, alpha = .6, color = "grey30",
             size = 0.25, lty = 2) +
  labs(
    x = "Threshold",
    y = "Metric" ) +
  theme(panel.grid = element_blank(),
        legend.position = "bottom")

p_dens =
  best_mods$test_fit[[4]] |>
  pull(.predictions) |>
  pluck(1) |>
  ggplot(aes(`.pred_ND-CNV`,colour = group,fill = group)) +
  geom_histogram(position = "dodge",alpha = 1,binwidth = 0.05,colour = "transparent",size = 0.01) +
  scale_fill_manual(values = c("#E1AF64","#32324B"))+
  geom_vline(xintercept = max_j_index_threshold, alpha = 1, color = "grey30",
             size = 0.25, lty = 2) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.title = element_blank(),
        legend.position = "top",
        legend.text = element_text(size = 6),
        axis.title = element_text(size = 6),
        axis.text = element_text(size = 6)) +
  labs(x = "Predicted Probability ND-CNV ",
       y = "Count")

(p_dens/p_thres + plot_layout(ncol = 1,heights = c(1,4)))



## Final variable Importance ======

#We could look at the importance of each of the variables in our final model ....

best_svm_importance = 
  best_mods |>
  filter(str_detect(mod_type,"SVM")) |>
  mutate(fit_engine = map(outer_full_fit,~extract_fit_engine(.x) ))


best_svm_varset = 
  best_svm_importance$best_wf_final[[1]] |>
  extract_preprocessor() %>% 
  .$var_info %>% 
  filter(role == "predictor") %>% 
  pull(variable)

best_svm_explainer = 
  explain_tidymodels(model = best_svm_importance$outer_full_fit[[1]]$.workflow[[1]],
                     data = d_last_fit |> analysis() |> select(best_svm_varset), 
                     y    = d_last_fit %>% analysis() %>% .$group == "ND-CNV")

#best_svm_performance = model_performance(best_svm_explainer)

#Permutation will involve random numbers, so set a seed for reproducibility
set.seed(22072022)

best_svm_parts = model_parts(explainer = best_svm_explainer,
                             type = "variable_importance",
                             B = 500)

# Save this given it takes a long time
write_rds(best_svm_parts, "./final_svm_variable_importance.rds")
best_svm_parts = read_rds("./final_svm_variable_importance.rds")


#Get the names for our best variables (this is a bit circular)
d_var = read_rds("nested_cv_selected_var_definitions_expanded.rds")

#Tabulate our best variables
best_svm_vars = 
  best_svm_parts |>
  as_tibble() |> 
  group_by(variable) |>
  ggdist::median_hdci(dropout_loss,.width = 0.95)|>
  arrange(-dropout_loss) 

# best_svm_vars|>
#   filter(variable %nin% c("_baseline_","_full_model_")) |>
#   left_join(d_var |> select(variable,short_name,var_def_long),by = "variable") |>
#   knitr::kable(format = "html", booktabs = TRUE) |>
#   kableExtra::kable_styling(font_size = 11)


p_best_var_imp = 
  best_svm_parts |>
  as_tibble() |> 
  group_by(variable) |>
  ggdist::median_hdci(dropout_loss,.width = 0.95)|>
  filter(variable %nin% c("_baseline_","_full_model_")) |>
  left_join(d_var |> select(variable,short_name,var_def_long),by = "variable") |>
  arrange(-dropout_loss) |>
  ggplot(aes(y = dropout_loss, x = forcats::fct_reorder(short_name,dropout_loss,min),
             ymin = .lower, ymax = .upper)) +
  geom_point(size = 0.5) +
  geom_linerange(size = 0.25) +
  geom_hline(yintercept = best_svm_vars |> filter(variable == "_full_model_") |> pull(dropout_loss),
             lty = 2,size = 0.25) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.title = element_blank(),
        legend.position = "top",
        legend.text = element_text(size = 6),
        axis.title = element_text(size = 6),
        axis.text = element_text(size = 6),
        axis.text.x.bottom = element_text(angle = 90)) +
  labs(y = "Dropout Loss (1 - AUROC)", x = "Variable") +
  scale_x_discrete(limits = rev)
       # title = "Best SVM Model",subtitle = "Vertical Line = Full Model")

p_best_var_imp


#Out of interest, how many variables have a dropout loss greater than the mean

best_svm_parts |>
  as_tibble() |> 
  group_by(variable) |>
  ggdist::median_hdci(dropout_loss,.width = 0.95) %>%
  filter(dropout_loss >  filter(., variable == "_full_model_") |> 
                          pull(dropout_loss) & 
           variable != "_baseline_") |>
  arrange(dropout_loss) |>
  print(n = 30)


#What if you scatter the two top variables?

p_12 = 
  d_last_fit |> 
  analysis() |> 
  select(group,best_svm_varset) |>
  ggplot(aes(x = P_Health_dev_15,y = P_ASQ_7, colour = group)) +
  geom_point(position = position_jitter(width = 0.3),alpha = 0.4) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  scale_color_brewer(palette = "Dark2")

p_23 = 
  d_last_fit |> 
  analysis() |> 
  select(group,best_svm_varset) |>
  ggplot(aes(x = P_ASQ_7,y = CMD_2, colour = group)) +
  geom_point(position = position_jitter(width = 0.3),alpha = 0.4) +
  theme_bw() +
  theme(panel.grid = element_blank())+
  scale_color_brewer(palette = "Dark2")


p_34 = 
  d_last_fit |> 
  analysis() |> 
  select(group,best_svm_varset) |>
  ggplot(aes(x = CMD_2,y = CMD_5, colour = group)) +
  geom_point(position = position_jitter(width = 0.3),alpha = 0.4) +
  theme_bw() +
  theme(panel.grid = element_blank())+
  scale_color_brewer(palette = "Dark2")


p_12 | p_23 | p_34

#You can see how clearly you can get good discrimination with just a few variables


## Figure 3: Performance Plot ======

#So lets make a plot with our model performance data

#What will we include? 

# ROC Curves for the different model types
# Calibration plots
# Thresholds
# Final variable importance


# Assemble
design = "123
          456
          456
          456"
p_perf = 
  plot_spacer() + p_dens + plot_spacer() + p_roc + p_thres + p_cal + plot_layout(design = design)


# p_performance = (p_perf / (p_age + p_gender) / p_best_var_imp) + plot_layout(heights = c(2,1,1))  + plot_annotation(tag_levels = "A")


ggsave("./Figures/figure_3.pdf",p_perf,width = 5, height = 3)
ggsave("./Figures/figure_3_vi.pdf",p_best_var_imp + plot_annotation(tag_levels = "A"),width = 5, height = 2)
ggsave("./Figures/figure_3_compare_vars.pdf",p_compare_vars + plot_annotation(tag_levels = "A"),width = 5, height = 2.5)


#What if we add in the training data performance



# Extract our final variables =====

#Lets get this set of variables:
best_vars = 
  final_recipes$SVM$var_info |> 
  pull(variable)

best_vars = 
  best_vars[!str_detect(best_vars,"group")]


#We can look these items up in the data dictionary....
#Read in Master participant list for confirmed genotypes
VL <- readxl::read_excel("DATA ENTRY CODEBOOK .xlsx", 
                         sheet = "VARIABLE LIST")

#Wrangle lightly
VL = 
  VL |>
  select(VARIABLE:`VARIABLE DEFINITION`) |>
  filter(VARIABLE %in% best_vars)

#Glue on the final SVM variable importance
VL =
  VL |>
  left_join(best_svm_vars |> 
              filter(variable %nin% c("_full_model_","_baseline_")) |> 
              rename(VARIABLE = variable) |>
              select(VARIABLE,dropout_loss,.lower,.upper), 
            by = "VARIABLE")

#And here we are
print(VL,n = 30)

#Save this
write_rds(VL,"nested_cv_selected_var_definitions.rds")


#Extract the smaller variable subset given this apparently offers identical performance

best_vars_small = 
  final_recipes$max.vars$var_info |> 
  pull(variable)

best_vars_small = 
  best_vars_small[!str_detect(best_vars_small,"group")]

#So just 10 items - are they all in the SVM variables?

setdiff(best_vars_small,best_vars)

#Actually no, there are 3 variables that are in the best set that are not in the 30...


#We can look these items up in the data dictionary....
#Read in Master participant list for confirmed genotypes
VL_small <- readxl::read_excel("DATA ENTRY CODEBOOK .xlsx", 
                         sheet = "VARIABLE LIST")

#Wrangle lightly
VL_small = 
  VL_small |>
  select(VARIABLE:`VARIABLE DEFINITION`) |>
  filter(VARIABLE %in% best_vars_small)

VL_small






# Other stuff =====

## Decision curve analysis =======

#We could do this, but we don't have a realistic measure of prevalence in our model because
#we have a very enriched sample for ND-CNV carriers
library(dcurves)

d_dca = 
  best_mods$test_fit[[4]] |>
  pull(.predictions) |>
  pluck(1) 

p_dca = 
  dca(group ~ `.pred_ND-CNV`, d_dca, thresholds = seq(0, .9, by = 0.01)) %>%
  plot(smooth = TRUE)


(p_dens|p_thres|p_dca) + plot_annotation(tag_levels = "A")

#' We could also think about doing ensembles here, but model performance is so good that it doesn't seem useful


