
# Introduction =========

#Here we take the most important variables from each of the dimensions
#Identified in our EGA analysis and we do a final round of ML with them


# Load packages and data ========

#Lets load up
pacman::p_load(tidyverse,tidymodels,workflowsets,tidyposterior,
               finetune,probably,doFuture,furrr,
               rstanarm,ggforce,doParallel,patchwork,
               readxl,lubridate, crayon)

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



#Load our recipes
final_recipes = read_rds("C://Users/nadon/OneDrive - University of Bristol/Documents/CNV Item Reduction/Data/final_selected_vars.rds")


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
  

  #We only use the Linear SVM
  svm_l_spec = 
    svm_linear(cost   = tune(),
               margin = tune()) %>%
    set_engine("kernlab") %>% 
    set_mode("classification")
  
  

  #Make the workflow set
  wf_par <- 
    workflow_set(
      preproc = final_recipes, 
      models  = list(SVM.linear   = svm_l_spec)
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
  


  #lets do some cleaning up
  rm(list = c("d_impute","d_inner", "d_outer_test", "d_outer",
              "rec_impute","rec_simple",
              "wf_par","wf_seq","wf_results","wf_seq_results","wf_par_results",
              "best_results_wf",
              "roc_mod","roc_mod_post"))
  
}


#Save this (it takes an hour or so to fit!)
write_rds(d_var_select_results,"C:/Users/nadon/OneDrive - University of Bristol/Documents/CNV Item Reduction/Data/nested_cv_result_ega_vars.rds")




## Summarise model performance =====


#Load all our models
d_var_ega_fa_results = read_rds("C:/Users/nadon/OneDrive - University of Bristol/Documents/CNV Item Reduction/Data/nested_cv_result_ega_vars.rds")
d_var_select_results = read_rds("C:/Users/nadon/OneDrive - University of Bristol/Documents/CNV Item Reduction/Data/nested_cv_result_selected_vars.rds")
d_fold_results       = read_rds("C:/Users/nadon/OneDrive - University of Bristol/Documents/CNV Item Reduction/Data/nested_cv_result_all_vars.rds")

### Make plots ######

#Lets plot some results

#Plot the posterior for all submodels
d_var_ega_fa_results |> 
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
               d_var_ega_fa_results |> 
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
                   d_var_ega_fa_results |> 
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
  coord_cartesian(ylim = c(.7,1)) 

#So our 6 variable set is much better than our 4 variable set

#Look at the averaged performance on the outer fold test data
d_var_ega_fa_results |> 
  select(-.results) |>
  unnest(.mod_posterior) |>
  group_by(vars,mod_type) |>
  summarise(mean_p = mean(posterior)) |>
  arrange(-mean_p)|>
  knitr::kable(format = "html", booktabs = TRUE) |>
  kableExtra::kable_styling(font_size = 11)


### Bayesian GLM of performance ==========

#Do something like the perf_mod function in tidymodels
cv_results = 
  d_var_ega_fa_results |>
  select(-.mod_posterior) |>
  unnest(.results) |>
  select(id,id2,wflow_id,best_roc_auc) |> 
  mutate(id = interaction(id,id2,sep = "_")) |>
  # separate(wflow_id,into = c("vars","mod_type"),sep = "_", extra = "merge") |>
  select(-c(id2))


#Add the data from the all variable models and the selected variable models
cv_results = 
  cv_results |>
  bind_rows(
    #All variable results
    d_fold_results |>
              select(-.mod_posterior) |>  
              unnest(.results) |>
              select(id,id2,wflow_id,best_roc_auc) |> 
              mutate(id = interaction(id,id2,sep = "_")) |>
              select(-id2) |>
              mutate(wflow_id = str_replace(wflow_id, "impute","all.vals")),
    
    #Selected variable results
    d_var_select_results |>
              select(-.mod_posterior) |>  
              unnest(.results) |>
              select(id,id2,wflow_id,best_roc_auc) |> 
              mutate(id = interaction(id,id2,sep = "_")) |>
              select(-id2)
            )

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


#Make our own plot
cv_mod_post = 
  cv_mod |>
  broom.mixed::tidy( conf.int = T, conf.level = 0.95, conf.method = "HPDinterval") |>
  mutate(term = str_remove(term,"wflow_id")) |>
  separate(term,into = c("vars","mod_type"),sep = "_", extra = "merge") 


cv_mod_post |>
  ggplot(aes(x = forcats::fct_reorder(vars,estimate,max),
             y = estimate,ymin =  conf.low,ymax = conf.high,
             colour = mod_type,
             fill   = mod_type)) +
  geom_errorbar(position = position_dodge(width = 0.5)) +
  geom_point(position = position_dodge(width = 0.5)) +
  theme_bw() +
  theme(axis.text.x.bottom = element_text(size = 8, angle = -45),
        legend.position = "none",
        panel.grid = element_blank()) +
  labs(y = "Model AUC", x = "Model Type",shape = "Variable \nSet")+
  facet_wrap(~mod_type,ncol = 44) +
  geom_hline(yintercept = 1)


#Now lets just focus on the SVM results

cv_mod_post |> 
  filter(mod_type == "SVM.linear") |>
  group_by(mod_type)  |> 
  arrange(-estimate) |> 
  mutate(across(where(is.double),round,digits = 4)) |>
  ungroup() |>
  transmute(Model = mod_type,`Variable Set` = vars, Performance = paste(estimate,"[",conf.low,",",conf.high,"]", sep = " ")) |>
  print(n = 34)




#You can do model comparison thus:
cv_mod_compare <- 
  cv_results |> 
  filter(str_detect(wflow_id, "SVM.linear")) |>
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


### Table 6 ======

#This is going to become table 6 so lets make it look nice
left_join(tab_mod_compare,
          cv_mod_post |> 
            mutate(across(where(is.double),round,digits = 3)) |>
            ungroup() |>
            transmute(Model = mod_type,
                      `Variable Set` = vars, 
                      Performance = paste(estimate,"[",conf.low,",",conf.high,"]", sep = " ")) , 
          by  = c("Model","Variable Set"))  |>
  relocate(Model, `Variable Set`, Performance,Difference,`Probability of Direction`) |>
  
  mutate(Model = case_when(Model == "SVM.linear" ~ "Linear SVM"),
         `Variable Set` = case_when(`Variable Set` == "all.vals" ~ "All Variables",
                                    `Variable Set` == "SVM" ~ "SVM",
                                    `Variable Set` == "lr.en" ~ "Penalised LR",
                                    `Variable Set` == "nnet" ~ "ANN",
                                    `Variable Set` == "rf" ~ "Random Forest",
                                    `Variable Set` == "max.vars" ~ "Variables selected by > 2 models",
                                    `Variable Set` == "moreThanOne.vars" ~ "Variables selected by > 1 model",
                                    `Variable Set` == "EGA" ~ "Top variables from 4 EGA Dimensions",
                                    `Variable Set` == "FA" ~ "Top variables from 6 Factors")) |>
  
  #Now make it nice and print it to the viewer
  knitr::kable(format = "html", booktabs = TRUE) |>
  kableExtra::kable_styling(font_size = 11)

#So performance is worse with both of our selected sets, although not that much worse


# Fit final models to the test data ========


## Prepare imputed test data =====


#As we are interested in using the use imputed dataset, we can load
#the set that we made earlier

d_last_fit        = read_rds("C:/Users/nadon/OneDrive - University of Bristol/Documents/CNV Item Reduction/Data/nested_cv_imputed_data.rds")
d_var_ega_results = read_rds("C:/Users/nadon/OneDrive - University of Bristol/Documents/CNV Item Reduction/Data/nested_cv_result_ega_vars.rds")


## Fit and evaluate models to full training data =====

# Now we fit our best models to the (imputed) training set

#Select the best set of parameters for each model - we will use the SVM variables as this appeared to have the best 
#performance (although differences between sets was very small)
best_mods = 
  d_var_ega_results |>
  select(-.mod_posterior) |>
  unnest(.results) |> 
  separate(wflow_id,
           into = c("vars","mod_type"),
           sep = "_", 
           extra = "merge") |>
  group_by(vars) |>
  slice_max(best_roc_auc)

best_mods = 
  best_mods |>
  mutate(test_fit = purrr::map(best_wf_final, 
                               last_fit,
                               split = d_last_fit,
                               metrics = metric_set(accuracy,kap,mn_log_loss,roc_auc,gain_capture))) |>
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



#What are the actual numbers of test participants classified
best_mods |>
  filter(mod_type == "SVM.linear" & vars == "EGA") |> 
  select(test_fit) |>
  ungroup() |>
  unnest(test_fit) |> 
  select(.predictions) |>
  unnest(.predictions) |>
  conf_mat(group,.pred_class) 



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
                                                           type = "prob") ,
                                                 .x |>
                                                   analysis()|>
                                                   select(group)) |>
                                mutate(group = fct_rev(group)) |>
                                c_m(group,`.pred_ND-CNV`)))
  
  return(bs)
  
}

#Now apply this function to our data

c_m = metric_set(roc_auc,mn_log_loss)

best_mods = 
  best_mods |>
  mutate(boot_metrics = map(test_fit, ~bootstrap_wf_metrics(wf = .x$.workflow[[1]],
                                                            data_set = d_last_fit |>
                                                              testing(),
                                                            n_boot = 2000,
                                                            c_m = c_m) ))



#Make a plot of our bootstrapped values
best_mods |>
  select(mod_type,vars,boot_metrics) |>
  unnest(boot_metrics)|>
  select(-splits) |>
  unnest(boot_metrics) |>
  ungroup() |>
  ggplot(aes(x = .estimate, y = vars)) +
  tidybayes::stat_halfeye() +
  facet_wrap(~.metric,ncol = 1)

#And we can make this into a more pleasing tabular format
tab_boot = 
  best_mods |>
  select(vars,boot_metrics) |>
  unnest(boot_metrics)|>
  select(-splits) |>
  unnest(boot_metrics) |>
  group_by(vars,.metric) |>
  ggdist::mean_hdci(.estimate) |>
  arrange(.metric,-.estimate)




## Table 7? ======

#We need the names of the variables
d_var = read_csv("C://Users/nadon/OneDrive - University of Bristol/Documents/CNV Item Reduction/Data/VariableDefinitionsExpanded.csv")

d_var = 
  d_var |>
  janitor::clean_names()


vars_EFA = final_recipes$FA$var_info |> 
  filter(variable != "group") %>% pull(variable) 


vars_EGA = final_recipes$EGA$var_info |> 
  filter(variable != "group") %>% pull(variable) %>% as.character()


tab_vars =
  tibble(vars = c("FA","EGA"),
       var_names = c(
         d_var |>
           filter(variable %in% vars_EFA) |>
           select(short_name,var_def_long) |>
           pull(short_name) |>
           paste(collapse = ", "),
         d_var |>
           filter(variable %in% vars_EGA) |>
           select(short_name,var_def_long) |>
           pull(short_name) |>
           paste(collapse = ", ")))




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


