
# Introduction =========

#Here we take the most important variables from each of the dimensions
#Identified in our EGA analysis and we do a final round of ML with them


# Load packages and data ========

#Lets load up
pacman::p_load(tidyverse,tidymodels,workflowsets,tidyposterior,
               finetune,probably,doFuture,furrr,themis,
               rstanarm,ggforce,doParallel,patchwork,
               readxl,lubridate, crayon)

tidymodels_prefer()

`%nin%` = Negate(`%in%`)


#Set our data directory
if(str_detect(Sys.info()[['nodename']],'IT088825')){
  data_dir = 'C://Users/pyxnd/OneDrive - University of Bristol/Documents/CNV Item Reduction/Data'
  
}else if(str_detect(Sys.info()[['nodename']],'AVOCADO')){
  
  data_dir = "C://Users/nadon/OneDrive - University of Bristol/Documents/CNV Item Reduction/Data/"
}

setwd(data_dir)


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



# Load =====



#Need to load up our final recipes here
d_split = read_rds("nested_cv_split_data.rds")

#Get test and training data
d_train = training(d_split)
d_test  = testing(d_split)

final_recipes = read_rds("final_selected_vars.rds")


# Prepare dataset for ML =====

#Now prepare for our nested cross validation

#Set a random number seed for reproducibility
set.seed(08032023)

#Use the test-train split we did above, but prepare new cv folds from the 
#training data
#Prepare for (nested) cross validation
n_outer    = 20
n_inner    = 20
n_fold     = 5

#Prepare for (nested) cross validation
d_folds = nested_cv(d_train, 
                    outside = vfold_cv(v = n_fold, repeats = n_outer/n_fold), 
                    inside  = bootstraps(times = n_inner)) 


#Preallocate an output tibble
d_ega_results = 
  d_folds |>
  select(id,id2) |>
  mutate(.results       = vector(mode = "list",length = nrow(d_folds)),
         .mod_posterior = vector(mode = "list",length = nrow(d_folds)))

fancy  <- combine_styles(make_style("ivory"), 
                         make_style("grey20", bg = TRUE))


#Loop through each outer fold
for(i in 1:nrow(d_folds)){
  
  cat(blue(
    'This is loop ' ,
      red$underline$bold(i) ,
      ' of ' ,
      red$underline$bold(nrow(d_folds)) 
  ))
  
  
  #Get the outer data for this iteration of the loop
  d_outer = d_folds$splits[[i]]
  
  #Get the inner data for this iteration of the loop
  d_inner = d_folds$inner_resamples[[i]]
  
  
  
  #Make pre-processing recipes
  cat(fancy("Pre-processing..."), "\n")
  
  #Set up a recipe with all the variables all at once, with imputation
  #Set up a recipe with all the variables all at once, with imputation
  rec_impute =
    d_outer |>
    analysis() |>
    select(-c(IDs,Gender,GenotypeCode,multi_demo)) |>
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
    set_engine("ranger", importance = "impurity") %>%
    set_mode("classification")
  
  
  #Make workflow sets so we can run all out models together using the workflowsets package
  wf_par <- 
    workflow_set(
      preproc = final_recipes,
      models  = list(en      = en_spec,
                     svm     = svm_spec,
                     nnet    = nnet_spec,
                     rf      = rf_spec))
  
  
  # Fit models =====
  
  
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
      grid = 25,
      metrics = metric_set(roc_auc,mn_log_loss)
    )
  
  
  #Turn off parallel mode
  plan(sequential)
  
  
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
    mutate(outer_full_fit = purrr::map(best_wf_final, last_fit,split = d_outer_split,metrics = metric_set(kap,mn_log_loss,roc_auc))) |>
    
    #Put those metrics somewhere easy to grab
    mutate(best_roc_auc = map_dbl(outer_full_fit,~.x |> 
                                    select(.metrics) |> 
                                    unnest(.metrics) |> 
                                    filter(.metric == "roc_auc") |> 
                                    pull(.estimate)))
  
  #Store the results for this fold
  d_ega_results$.results[[i]] =
    best_results_wf |>
    dplyr::select(-c(res,best_wf)) 
  
  #Store the performance model data too
  d_ega_results$.mod_posterior[[i]] = 
    roc_mod_post
  



  #lets do some cleaning up
  rm(list = c("d_impute","d_inner", "d_outer_test", "d_outer",
              "rec_impute","rec_simple",
              "wf_par","wf_seq","wf_results","wf_seq_results","wf_par_results",
              "best_results_wf",
              "roc_mod","roc_mod_post"))
  
}


#Save this (it takes an hour or so to fit!)
write_rds(d_ega_results,"nested_cv_result_ega_vars.rds")


# Summarise model performance =====


#Load all our models
d_ega_results = read_rds("nested_cv_result_ega_vars.rds")


#Load all our models
d_var_select_results = read_rds("nested_cv_result_selected_vars_v2.rds")
d_fold_results       = read_rds("nested_cv_result_all_vars_v2.rds")


### Make plots ######

#Lets plot some results

#Plot the posterior for all submodels
d_ega_results |> 
  select(-.results) |>
  unnest(.mod_posterior) |>
  mutate(mod_type = factor(mod_type,levels = c("en","nnet","rf","svm"))) |>
  ggplot(aes(x = mod_type,
             y = posterior,
             group = mod_type,
             ymin = .lower,ymax = .upper,
             colour = mod_type)) +
  geom_point(aes(shape = mod_type),
             position = position_jitter(seed = 123,width = 0.2)) +
  geom_linerange(position = position_jitter(seed = 123,width = 0.2)) +
  geom_point(data  =  
               d_ega_results |> 
               select(-.results) |>
               unnest(.mod_posterior) |>
               group_by(mod_type) |>
               summarise(sd        = sd(posterior),
                         posterior = mean(posterior)) |>
               mutate(.lower = posterior-sd, .upper= posterior+sd),
             aes(shape = mod_type),
             size = 4) +
  geom_linerange(data  =  
                   d_ega_results |> 
                   select(-.results) |>
                   unnest(.mod_posterior) |>
                   group_by(mod_type) |>
                   summarise(sd        = sd(posterior),
                             posterior = mean(posterior)) |>
                   mutate(.lower = posterior-sd, .upper= posterior+sd),
                 linewidth = 2,
                 alpha = 0.6,
                 lty = 1) +
  theme_bw() +
  theme(axis.text.x.bottom = element_text(size = 8, angle = -45),
        legend.position = "none") +
  labs(y = "Model AUC", x = "Variable Set",shape = "Variable \nSet") +
  coord_cartesian(ylim = c(.85,1))



### Bayesian GLM of performance ==========

#We fit a model to the brier scores and ROCs together 
library(brms)
library(tidybayes)
unloadNamespace("conflicted")

d_ega_metrics = 
  d_ega_results |>
  select(-.mod_posterior ) |>
  unnest(.results) |>
  select(id,id2,wflow_id,best_roc_auc,outer_full_fit) |>
  mutate(best_brier = map_dbl(outer_full_fit,~.x |> pluck('.predictions',1) |>
                                brier_score('group','.pred_ND-GC'))) |>
  select(-outer_full_fit) |>
  mutate(wflow_id = str_remove(wflow_id,"EGA_")) |>
  mutate(id = interaction(id,id2,sep = "_")) |>
  select(-c(id2)) 


r_model <- bf(best_roc_auc ~ 0 + wflow_id + (1|id))
b_model <- bf(  best_brier ~ 0 + wflow_id + (1|id))


get_prior(formula = r_model + b_model + set_rescor(TRUE),
          data = d_ega_metrics)


#Fit the multivariate BRMS model
fitM <- 
  brm(
    data = d_ega_metrics,
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
  mutate(name = str_remove(name, "wflow_id")) %>%
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
  labs(subtitle = "BRMS Multivariate Model - Normal")+
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




# Fit final models to the test data ========

#Now we can fit our final models using the held out test data


## Prepare imputed test data =====

#As we are interested in using the use imputed dataset, we can load
#the set that we made earlier
d_last_fit    = read_rds("nested_cv_imputed_data.rds")
d_ega_results = read_rds("nested_cv_result_ega_vars.rds")


## Fit and evaluate models to full training data =====

# Now we fit our best models to the (imputed) training set

#Select the best set of parameters for each model - we will use the SVM variables as this appeared to have the best 
#performance (although differences between sets was very small)
best_mods = 
  d_ega_results |>
  unnest(.results) |> 
  mutate(model = str_remove(wflow_id,"EGA_")) |>
  select(-wflow_id) |>
  group_by(model) |>
  slice_max(best_roc_auc,n =  1)

best_mods = 
  best_mods |>
  mutate(test_fit = purrr::map(best_wf_final, 
                               last_fit,
                               split = d_last_fit,
                               metrics = metric_set(mn_log_loss,roc_auc))) |>
  mutate(best_metrics = map(test_fit, ~.x$.metrics[[1]]))



## Performance measures =====

#Tabulate our final performance
tab_mods = 
  best_mods |>
  
  #We can work out the Brier score (this seems to simply be the average difference between the predicted 
  #probability and true outcome squared i.e. the mean squared error) - there is also the mn log loss measure
  #in yarstick which is the same thing on a slightly different scale (the log scale)
  mutate(brier = map_dbl(test_fit,~ .x$.predictions[[1]] |>
                           transmute(pred = `.pred_ND-GC`,group) |>
                           mutate(group = ifelse(group == "ND-GC",1,0)) |>
                           mutate(dist  = (group - pred)^2) |>
                           summarise(brier = mean(dist)) |>
                           pull())) |>
  select(model,best_metrics,brier)|>
  unnest(best_metrics) |>
  select(model,.metric,.estimate,brier) |>
  pivot_wider(names_from = .metric,values_from = .estimate) 
 

tab_mods |>
  knitr::kable(format = "html", booktabs = TRUE) |>
  kableExtra::kable_styling(font_size = 11)



## Bootstrap some confidence intervals =====


#Estimate confidence intervals for the final performance metrics by bootstrapping
#and used downsampling to get balanced data in each bootstrap
bootstrap_wf_metrics <- function(wf,data_set,n_boot = 500){
  
  
  #We need an input workflow and an input dataset 
  
  #We then make this into a bootstrap object from rsample
  d_b = 
    data_set |>
    bootstraps(times = n_boot, strata = "group") 
  
  
  d_b = 
    d_b |> 
    mutate(d_boot = map(splits,~analysis(.x))) |>
    select(-splits) 
  
  #resample to balance groups
  d_b = 
    d_b |>   
    mutate(d_resample = map(d_boot,~.x |> 
                              recipe(group ~ .) |>
                              # step_upsample(group,
                              #               over_ratio = 1,
                              #               seed = 23022023,
                              #               skip = TRUE) |>
                              step_downsample(group,
                                              under_ratio = 1,
                                              seed = 23022023,
                                              skip = TRUE) |>
                              prep() |>
                              juice()
    )) |>
    
    #Predict the group of the resampled test data
    mutate(pred = map(d_resample, ~predict(wf,.x,type = "prob") |>
                        bind_cols(.x |>
                                    select(group)) |>
                        mutate(group = fct_rev(group)) 
    )) |>
    
    #Calculate our outcome metrics
    mutate(metrics = map(pred,~ .x |>
                           roc_auc(group,`.pred_ND-GC`) |>
                           bind_rows(
                             tibble(.metric = "brier_score",
                                    .estimator = "binary",
                                    .estimate = .x |> 
                                      brier_score("group",".pred_ND-GC")))))
  
  #Extract our final metrics
  boot_out = 
    d_b |>
    select(id,metrics) |>
    unnest(metrics)
  
  
  return(boot_out)
  
}


#Apply this function to our data
boot_mods = 
  best_mods |> 
  select(model,test_fit) |> 
  unnest(test_fit)|>
  mutate(boot_metrics = map(.workflow, ~bootstrap_wf_metrics(wf = .x,
                                                             data_set = 
                                                               d_last_fit |>
                                                               testing(),
                                                             n_boot = 999)))


#Make a plot of our bootstrapped values
boot_mods |>
  select(model,boot_metrics) |>
  unnest(boot_metrics)|>
  ggplot(aes(x = .estimate, y = model)) +
  tidybayes::stat_halfeye() +
  facet_wrap(~.metric,ncol = 1)


#And we can make this into a more pleasing tabular format
tab_boot = 
  boot_mods |>
  select(model,boot_metrics) |>
  unnest(boot_metrics)|>
  group_by(model,.metric) |>
  ggdist::mean_hdci(.estimate) |>
  arrange(.metric,-.estimate)


## Threshold -----

#Looks like the EN model has the best performance
threshold_data = 
  best_mods |>
  ungroup() |>
  filter(model == "en" ) |>
  pluck("test_fit",1) |>
  pluck(".predictions",1) |>
  mutate(group = fct_rev(group)) |>
  probably::threshold_perf(group,`.pred_ND-GC`,
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

## Confusion Matrix =====


#What are the actual numbers of test participants classified
best_mods |>
  select(model,test_fit) |>
  ungroup() |>
  filter(model == "en" ) |>
  unnest(test_fit) |> 
  select(model,.predictions) |>
  unnest(.predictions) |>
  mutate(.pred_class = if_else(`.pred_ND-GC` >= max_j_index_threshold,"ND-GC","Control")) |>
  mutate(.pred_class = factor(.pred_class)) |>
  conf_mat(group,.pred_class) 


## Table 7? ======

#We need the names of the variables
d_var = read_csv("C://Users/nadon/OneDrive - University of Bristol/Documents/CNV Item Reduction/Data/VariableDefinitionsExpanded.csv")



vars_EGA = 
  final_recipes$EGA$var_info |> 
  filter(variable != "group") %>% pull(variable) %>% as.character()


d_var |>
  filter(variable %in% vars_EGA) |>
  select(short_name,variable_definition_dict) 


tab_vars =
  tibble(vars = c("EGA"),
       var_names = d_var |>
           filter(variable %in% vars_EGA) |>
           select(short_name,variable_definition_dict) |>
           pull(short_name) |>
           paste(collapse = ", "))


#Now tabulate

tab_boot |>
  left_join(tab_vars,by = "vars") |>
  mutate(across(where(is.double),round,digits = 3)) |>
  mutate(Performance = paste(.estimate," [",.lower,", ",.upper,"]",sep = "")) |>
  select(vars,var_names,.metric,Performance) |>
  rename(`Variable Set` = vars,
         `Variable Names` = var_names) |>
  pivot_wider(names_from = .metric,values_from = Performance) |>
  rename(`Mean Log Loss` = mn_log_loss,AUROC = roc_auc) |>
  relocate(`Variable Set`,`Variable Names`,AUROC,`Mean Log Loss`) |>
  knitr::kable(format = "html", booktabs = TRUE) |>
  kableExtra::kable_styling(font_size = 11)


# Final LR -----

#Sometimes logistic regression is enough


fitted_logistic_model <- 
  
  logistic_reg() %>%
  
  # Set the engine
  set_engine("glm") %>%
  # Set the mode
  set_mode("classification") %>%
  
  # Fit the model
  fit(group~., 
      data = 
        d_last_fit |> 
        training() |> 
        dplyr::select(group, all_of(best_mods$test_fit[[1]] |> 
                                      extract_recipe() %>%
                                      .$var_info |>
                                      pull(variable)))
  )


#Model summary table 
tidy(fitted_logistic_model,exponentiate = T,conf.int = T)   |>
  filter(term != "(Intercept)") |>
  arrange(p.value) |>
  print(n = 5)


#Model coefs
fitted_logistic_model$fit$coefficients |>
  as_tibble(rownames = "term") |>
  print(n = 5)


#Predicted classes
pred_class <- 
  predict(fitted_logistic_model,
          new_data = d_last_fit |> 
            testing(),
          type = "class")


# Prediction Probabilities
pred_proba <- 
  predict(fitted_logistic_model,
          new_data = d_last_fit |> testing(),
          type = "prob")

#Join them
cnv_results <- 
  d_last_fit |> testing() %>%
  select(group) %>%
  bind_cols(pred_class, pred_proba)


#Various model metrics
conf_mat(cnv_results, truth = group,
         estimate = .pred_class)

mcc(cnv_results, truth = group,
    estimate = .pred_class)

mn_log_loss(cnv_results, truth = group,
            estimate = `.pred_ND-GC`,event_level = "second")

roc_auc(cnv_results, truth = group,
        estimate = `.pred_ND-GC`,event_level = "second")

kap(cnv_results, truth = group,
    .pred_class,event_level = "second")     

#Plot probs
cnv_results  |> 
  ggplot(aes(`.pred_ND-GC`,fill = group)) + 
  geom_histogram(binwidth = .05,
                 position = position_dodge())


#Its not terrible actually