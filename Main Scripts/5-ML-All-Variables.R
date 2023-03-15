
# Introduction =========

#What we do in this script is to combine imputation, variable selection
#and model fitting into one nested cross-validation process:

# - We start with the data filtered of the non-quantitative variables
# - Then we will do our initial test/training split
# - Then with the training data we will split into nested CV
# - The inner folds will be used for the process of imputation and model training 
# - The best models will be used for variable importance calculation.
# - Next we tune machine learning models with sets of the most important variables
#   using nested CV once more
# - The best performing models are then applied to the test data and we measure
#   model calibration


# Load packages and data ========

#Lets load up
pacman::p_load(tidyverse,tidymodels,workflowsets,tidyposterior,
               gtsummary,
               finetune,probably,doFuture,furrr,themis,
               rstanarm,ggforce,doParallel,patchwork,
               vip,readxl,lubridate, crayon)

conflicted::conflicts_prefer(crayon::`%+%`)
conflicted::conflicts_prefer(dplyr::filter)

#Set our data directory
if(str_detect(Sys.info()[['nodename']],'IT088825')){
  data_dir = 'C://Users/pyxnd/OneDrive - University of Bristol/Documents/CNV Item Reduction/Data'
  
}else if(str_detect(Sys.info()[['nodename']],'AVOCADO')){
  
  data_dir = "C://Users/nadon/OneDrive - University of Bristol/Documents/CNV Item Reduction/Data/"
}

setwd(data_dir)



#Load our previously split data

d_split = read_rds("nested_cv_split_data.rds")

#Get test and training data
d_train = training(d_split)
d_test  = testing(d_split)


#Get the IDs of the participants in the testing split
d_test_id = d_split |> complement()

#Prepare for (nested) cross validation
n_outer    = 20
n_inner    = 20
n_fold     = 5


d_folds = nested_cv(d_train, 
                    outside = vfold_cv(v = n_fold, repeats = n_outer/n_fold, strata = group), 
                    inside  = bootstraps(times = n_inner, strata = group)) 


# So the state of play is that we have:

# - The full dataset with 489 individuals and 233 variables
# - The testing holdout dataset, which we will not be looking at until the very end, with 100 individuals
# - Nested CV datasets; each outer fold contains 390 individuals, split into 351 analysis, 39 assessment
# - Nested CV inner folds contain 351 individuals, 123 assessment, 351 analysis bootstraps

# Nested CV =====

#We are going to train ML models on the full variable set. To do nested CV we are going to use loops
#I know loops are note very cool in R but this is a good way to get going quickly. Plus I'm used to MATLAB

#Preallocate an output tibble
d_fold_results = 
  d_folds |>
  select(id,id2) |>
  mutate(.results       = vector(mode = "list",length = nrow(d_folds)),
         .mod_posterior = vector(mode = "list",length = nrow(d_folds)))


#Set up some preferences for the crayon package so we can print progress information to the command line
#as our ML process can take a few hours and its helpful to have updates on progress
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
  
  
  
  ### Build recipes =====
  
  #Make pre-processing recipes
  
  cat(fancy("Pre-processing..."), "\n")
  

  #Set up a recipe without any imputation
  rec_simple = 
    d_outer |>
    analysis() |>
    select(-c(IDs,Gender,GenotypeCode,multi_demo)) |>
    recipe(group ~ .) |>
    
    #Remove any variables with zero variance (shouldn't do anything)
    step_zv(all_predictors())
  
  
  
  ## Prepare imputed data ======
  
  #Set up a recipe with imputation of missing data  
  rec_impute =
    rec_simple |>
    step_impute_bag(all_predictors())   
  
  
  #Now, rather than impute on every single model training step,
  #we are going to do the imputation on this outer fold here
  
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
  
  
  
  ## Prepare up/down sampling ======
  
  
  #We also need to test means of managing the unbalanced nature of the dataset to address
  #a reviewer comment
  
  #Downsample the cases to the same number as the controls - will reduce power and make
  #everything less accurate
  rec_down =
    rec_simple |>
    step_downsample(group,
                    under_ratio = 1,
                    seed = 23022023,
                    skip = TRUE)
  
  #Upsample controls by resampling
  rec_up =
    rec_simple |>
    step_upsample(group,
                    over_ratio = 1,
                    seed = 23022023,
                    skip = TRUE)
  
  
  
  #Lets inspect the results of our up and down sampling methods

  # d_impute |>
  #   ggplot(aes(x = CMD_1,y = CMD_3, colour = group)) +
  #   geom_jitter(width = 0.1,height = 0.1)
  # 
  # 
  # rec_down |>
  #   prep(training = d_impute) |>
  #   bake(new_data = NULL) |>
  #   ggplot(aes(x = CMD_1,y = CMD_3, colour = group)) +
  #   geom_jitter(width = 0.1,height = 0.1)
  # 
  # 
  # rec_up |>
  #   prep(training = d_impute) |>
  #   bake(new_data = NULL) |>
  #   ggplot(aes(x = CMD_1,y = CMD_3, colour = group)) +
  #   geom_jitter(width = 0.1,height = 0.1)
  # 
  # d_impute |>
  #   select(group,pcb1i01,pcb3i01) |>
  #   group_by(group) |>
  #   summarise(mean_fpb = mean(pcb1i01),
  #             mean_ago = mean(pcb3i01),
  #             prop_fpb = sum(pcb1i01 == 1)/n(),
  #             prop_ago = sum(pcb3i01 == 1)/n())
  


  ## Fit Models =====
  
  #Next we set up candidate ML models. Agreed with Marianne and Valentina:

  # - Penalised Regression (i.e. elastic net regression)
  # - Radial Basis Function SVM (based on reviewer comment)
  # - Random Forests
  # - Neural Network (with nnet for speed)
  
  
  cat(fancy("Fitting Models..."), "\n")
  
  
  #Elastic Net Regression with glmnet
  en_spec <- 
    logistic_reg(mode = "classification",
                 penalty = tune(), 
                 mixture = tune()) %>% 
    set_engine("glmnet")
  
  
  #RBF SVM with the kernlab package
  svm_spec = 
    svm_rbf(cost      = tune(),
            rbf_sigma = tune(),
            margin    = tune()) %>%
    set_engine("kernlab") %>% 
    set_mode("classification")
  
  
  #Neural Network with nnet (note I have tried other neural network models like brulee/tabnet - 
  #they use a lot of cpu/gpu power without improving performance so we keep things simple)
  nnet_spec <- 
    mlp(hidden_units = tune(), 
        penalty      = tune(), 
        epochs       = tune()) %>% 
    set_engine("nnet", MaxNWts = 2600) %>% 
    set_mode("classification")
  
  
  #Random Forest (I have experimented with gradient boosted trees e.g. XGBoost and ligthGBM and
  #again they are slower and do not perform better in testing so we use random forests with ranger
  #for speed)
  rf_spec <-
    rand_forest(mtry = tune(), min_n = tune(), trees = 1000) %>%
    set_engine("ranger") %>%
    set_mode("classification")

  
  #Make workflow sets so we can run all out models together using the workflowsets package
  wf_par <- 
    workflow_set(
      preproc = list(impute  = rec_simple,
                     up      = rec_up,
                     down    = rec_down),
      models  = list(en      = en_spec,
                     svm     = svm_spec,
                     nnet    = nnet_spec,
                     rf      = rf_spec))
  
  ## Do some model fitting

  #Set up controls for the finetune anova race method
  race_ctrl <-
    control_race(
      save_pred = TRUE,
      parallel_over = "everything",
      save_workflow = TRUE
    )
  
  
  #Set up to run these in parallel (I have a 16 core CPU so we use all 16 cores)
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
  # autoplot(wf_results, id ="impute_svm")
  
  
  ## Measure Model Performance =====
  
  
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

  # #Do the ROPE Plot
  # autoplot(roc_mod, type = "ROPE", size = 0.025)

  
  #Evaluate models on the outer fold data
  
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
    
    #We finalize the workflow by combining the best parameters with the workflows
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
best_results_wf |>
  select(wflow_id,best_roc_auc) |>
  mutate(wflow_id = map_chr(wflow_id,str_remove,"impute_")) |>
  ggplot(aes(x = forcats::fct_reorder(wflow_id,best_roc_auc,max),
             y = best_roc_auc)) +
  geom_point() +
  theme_bw() +
  theme(axis.text.x.bottom = element_text(angle = 0)) +
  labs(x = "Model", y = "Outer ROC") +
  coord_cartesian(ylim = c(0.8,1))

  
  #lets do some cleaning up
  rm(list = c("d_impute","d_inner", "d_outer_test", "d_outer",
              "rec_impute","rec_simple",
              "wf_par","wf_seq","wf_results","wf_seq_results","wf_par_results",
              "best_results_wf",
              "roc_mod","roc_mod_post"))  
  
}

# Summarise model performance =====

#Save this (it takes an hour or so to fit!)
write_rds(d_fold_results,"C://Users/nadon/OneDrive - University of Bristol/Documents/CNV Item Reduction/Data/nested_cv_result_all_vars_v2.rds")

# d_fold_results = read_rds("C://Users/nadon/OneDrive - University of Bristol/Documents/CNV Item Reduction/Data/nested_cv_result_all_vars_v2.rds")


### Plot results ===== 

#Lets plot some results

#Plot the posterior for all submodels
d_fold_results |> 
  select(-.results) |>
  unnest(.mod_posterior) |>
  filter(vars == "impute") |>
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
             linewidth = 2,
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
  filter(vars == "impute") |>
  group_by(mod_type) |>
  summarise(mean_p = mean(posterior)) |>
  mutate(mod_type = case_when(mod_type == "en" ~ "Penalised LR",
                              mod_type == "nnet" ~ "ANN",
                              mod_type == "rf" ~ "Random Forest",
                              mod_type == "svm" ~ "Linear SVM")) |>
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
  filter(str_detect(wflow_id,"impute")) |>
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
  mutate(mod_type = case_when(mod_type == "en" ~ "Penalised LR",
                              mod_type == "nnet" ~ "ANN",
                              mod_type == "rf" ~ "Random Forest",
                              mod_type == "svm" ~ "RBF SVM")) |>
  ggplot(aes(x = forcats::fct_reorder(mod_type,estimate,max),
             y = estimate,ymin =  conf.low,ymax = conf.high,
         colour = vars,
         fill = vars)) +
  # geom_errorbar(position = position_dodge(width = 0.5)) +
  geom_linerange(position = position_dodge(width = 0.2),
                 linewidth = 1.5) +
  geom_point(position = position_dodge(width = 0.2),
             size = 3) +
  theme_bw() +
  theme(axis.text.x.bottom = element_text(size = 8, angle = -45),
        panel.grid = element_blank()) +
  labs(y = "Model AUC", x = "Model Type",shape = "Variable \nSet") +
  coord_cartesian(ylim = c(0.8,1)) +
  geom_hline(yintercept = 1,lty = 3) 

p_all_vars


#Tabulate
tab_post = 
  cv_mod_post |>
  mutate(mod_type = case_when(mod_type == "en" ~ "Penalised LR",
                              mod_type == "nnet" ~ "ANN",
                              mod_type == "rf" ~ "Random Forest",
                              mod_type == "svm" ~ "RBF SVM")) |>
  mutate(across(where(is.double),~round(.x,digits = 3))) |>
  arrange(-estimate) |>
  transmute(Recipe = vars,
            Model = mod_type,
            Performance = paste(estimate, "[",conf.low, ",",conf.high,"]",sep = " "))

# tab_post |>
#   knitr::kable(format = "html", booktabs = TRUE) |>
#   kableExtra::kable_styling(font_size = 11)


#You can do model comparison thus:
cv_mod_rel <- 
  cv_results |>
  mutate(wflow_id = factor(wflow_id)) |>
  mutate(wflow_id = fct_relevel(wflow_id,"impute_svm")) |>
  stan_glmer(best_roc_auc ~ 1 + wflow_id  + (1|id),
             data = _,
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
  mutate(Parameter = if_else(Parameter == "(Intercept)","impute_svm",str_remove(Parameter,"wflow_id"))) |>
  mutate(across(where(is.double),round,digits = 3)) |>
  transmute(Model = Parameter,
            Difference = paste(Median, "[",CI_low, ",",CI_high,"]",sep = " "),
            `Probability of Direction` = bayestestR::pd_to_p(pd)) |>
  separate(Model,into = c("Recipe","Model"),sep = "_", extra = "merge") |>
  mutate(Model = case_when(Model == "en" ~ "Penalised LR",
                           Model == "nnet" ~ "ANN",
                           Model == "rf" ~ "Random Forest",
                           Model == "svm" ~ "RBF SVM")) 

tab_mod_rel$Difference[1] = "-"
tab_mod_rel$`Probability of Direction`[1] = 1


### Table 3 ========

#This makes a table of both performance and model comparison
left_join(tab_post,tab_mod_rel, by  = c("Recipe","Model"))  |>
  select(-Recipe) |>
  knitr::kable(format = "html", booktabs = TRUE) |>
  kableExtra::kable_styling(font_size = 11)


#So the SVM has the best overall performance, although the LR and RF models are very similar;
#however the ANN performance is worse (although in absolute terms still very good)

# At this point we could fit the best performing models to the full training data and test on the 
# previous held out testing data. However, as we are going to do variable selection and fit models 
# with reduced sets of variables, I think we stop now and look at variable importance, and then we
# return to fit to the full data and test in the held out test data at the very end.

# Multiple Performance Measures  #######

#We could evaluate multiple performance measures at once

#Define a function to calculate the Brier score (mean squared error)
brier_score = function(data,truth,estimate){
  
  data |>
    rename(pred  = all_of(estimate),
           group = all_of(truth)) |>
    
    mutate(group = ifelse(group == "ND-GC",1,0)) |>
    mutate(dist  = (group - pred)^2) |>
    summarise(brier = round(mean(dist),digits = 3) )|>
    pull()
  
}


#We fit a model to the brier scores and ROCs together 
library(brms)
library(tidybayes)
unloadNamespace("conflicted")

d_metrics = 
  d_fold_results |>
  select(-.mod_posterior ) |>
  unnest(.results) |>
  select(id,id2,wflow_id,best_roc_auc,outer_full_fit) |>
  mutate(best_brier = map_dbl(outer_full_fit,~.x |> pluck('.predictions',1) |>
                                brier_score('group','.pred_ND-GC'))) |>
  select(-outer_full_fit) |>
  separate_wider_delim(wflow_id,names = c("sampling","model"),delim = "_") |>
  mutate(id = interaction(id,id2,sep = "_")) |>
  select(-c(id2)) |>
  filter(sampling == "impute") |>
  select(-sampling)


r_model <- bf(best_roc_auc ~ 0 + model + (1|id))
b_model <- bf(  best_brier ~ 0 + model + (1|id))


get_prior(formula = r_model + b_model + set_rescor(TRUE),
          data = d_metrics)


#Fit the multivariate BRMS model
fitM <- 
  brm(
    data = d_metrics,
    family = gaussian,
    
    formula = r_model + b_model + set_rescor(TRUE),
    
    prior = c(
      
      #ROC model
      prior(normal(0, 0.5), class = b    , resp = bestrocauc),
      prior(exponential(1), class = sd   , resp = bestrocauc),
      prior(exponential(1), class = sigma, resp = bestrocauc),
      
      
      #Brier Score model
      prior(normal(0, 0.5), class = b    , resp = bestbrier),
      prior(exponential(1), class = sd   , resp = bestbrier),
      prior(exponential(1), class = sigma, resp = bestbrier),
      
      
      #Outcome variable correlation
      prior(lkj(2), class = rescor)),
    
    chains = 4, cores = 8,
    warmup = 1000, iter = 5000, seed = 12)



summary(fitM)

pp_check(fitM,resp = "bestrocauc")
pp_check(fitM,resp = "bestbrier")


mv_p1 <- 
  fitM %>% 
  tidy_draws() %>%
  select(starts_with("b_")) %>%
  pivot_longer(everything()) %>%
  dplyr::filter(!grepl("Intercept",name)) %>%
  mutate(name = str_remove(name, "b_best")) %>%
  mutate(name = str_remove(name, "model")) %>%
  separate_wider_delim(name,names = c("outcome", "model"),delim = "_") |>
  mutate(model = case_when(model == "en" ~ "Penalised LR",
                           model == "nnet" ~ "ANN",
                           model == "rf" ~ "Random Forest",
                           model == "svm" ~ "RBF SVM")) |>
  mutate(model = factor(model,levels = c("ANN","Penalised LR","Random Forest","RBF SVM"))) |>
  mutate(outcome = case_when(outcome == "rocauc" ~ "ROC AUC",
                             outcome == "brier" ~ "Brier Score")) |>
  mutate(outcome = factor(outcome,levels = c("ROC AUC","Brier Score"))) |>
  ggplot(aes(y = value, 
             x = model,
             group = model,
             colour = model,
             shape = model)) +
  stat_pointinterval(point_interval = median_hdi, .width = c(0.66,0.95), position = position_dodge(width = 0.7)) +
  geom_vline(xintercept = 0, colour = "black", linetype = 3) +
  labs(subtitle = "BRMS Multivariate Model - All Variable Performance")+
  # coord_cartesian(xlim = c(-.5,.75)) +
  ghibli::scale_colour_ghibli_d("TotoroMedium", direction = -1) +
  theme_bw() +
  facet_wrap(~outcome,scales = "free_y")

mv_p1


#Tabulate our data
fitM %>%
  bayestestR::describe_posterior(ci = 0.95,
                                 ci_method = "hdi",
                                 centrality = "median") %>%
  as_tibble() %>%
  mutate(Parameter = str_remove(Parameter, "wflow_id"),
         Parameter = str_remove(Parameter, "b_"),) %>%
  separate_wider_delim(Parameter,names = c("outcome", "model"),delim = "_") |>
  select(-c(CI,Rhat, ESS,pd,starts_with("ROPE"))) %>%
  mutate(across(where(is.double),~round(.x,  digits = 3))) %>%
  mutate(Performance = paste(Median,"[",CI_low,",",CI_high,"]", sep = " "))  |>
  select(outcome,model,Performance) |>
  mutate(model = case_when(model == "en" ~ "Penalised LR",
                           model == "nnet" ~ "ANN",
                           model == "rf" ~ "Random Forest",
                           model == "svm" ~ "RBF SVM")) |>
  mutate(model = factor(model,levels = c("ANN","Penalised LR","Random Forest","RBF SVM"))) |>
  pivot_wider(names_from = "outcome",values_from = "Performance") |>
  knitr::kable(format = "html", booktabs = TRUE) |>
  kableExtra::kable_styling(font_size = 11)

