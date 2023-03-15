
# Introduction ------

# ML using the finalised selected variable datasets

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



#Need to load up our final recipes here
d_split = read_rds("nested_cv_split_data.rds")

#Get test and training data
d_train = training(d_split)
d_test  = testing(d_split)

final_recipes = read_rds("nested_cv_selected_vars.rds")

## Make training data structure =====


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
d_var_select_results = 
  d_folds |>
  select(id,id2) |>
  mutate(.results       = vector(mode = "list",length = nrow(d_folds)),
         .mod_posterior = vector(mode = "list",length = nrow(d_folds)))

fancy  <- combine_styles(make_style("ivory"), 
                         make_style("grey20", bg = TRUE))


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
  rm(list = c("d_impute","d_inner", "d_outer_test", "d_outer","d_outer_split",
              "rec_impute",
              "wf_par","wf_results","best_results_wf",
              "roc_mod","roc_mod_post"))
  
}


#Save this (it takes an hour or so to fit!)
write_rds(d_var_select_results,"nested_cv_result_selected_vars_v2.rds")


# Summarise model performance =====

#Load all our models
d_var_select_results = read_rds("nested_cv_result_selected_vars_v2.rds")
d_fold_results       = read_rds("nested_cv_result_all_vars_v2.rds")

### Make plots ######

#Lets plot some results

#Plot the posterior for all submodels
d_var_select_results |> 
  select(-.results) |>
  unnest(.mod_posterior) |>
  filter(vars == "simple") |>
  select(-vars) |>
  separate_wider_delim(mod_type,names = c("var_set","model"),delim = "_") |>
  mutate(var_set = factor(var_set,levels = c("en","nnet","rf","svm","moreThanOne","max"))) |>
  ggplot(aes(x = var_set,
             y = posterior,
             group = var_set,
             ymin = .lower,ymax = .upper,
             colour = model)) +
  geom_point(aes(shape = model),
             position = position_jitter(seed = 123,width = 0.2)) +
  geom_linerange(position = position_jitter(seed = 123,width = 0.2)) +
  geom_point(data  =  
               d_var_select_results |> 
               select(-.results) |>
               unnest(.mod_posterior) |>
               filter(vars == "simple") |>
               select(-vars) |>
               separate_wider_delim(mod_type,names = c("var_set","model"),delim = "_") |>
               group_by(var_set,model) |>
               summarise(sd        = sd(posterior),
                         posterior = mean(posterior)) |>
               mutate(.lower = posterior-sd, .upper= posterior+sd),
             aes(shape = model),
             size = 4)+
  geom_linerange(data  =  
                   d_var_select_results |> 
                   select(-.results) |>
                   unnest(.mod_posterior) |>
                   filter(vars == "simple") |>
                   select(-vars) |>
                   separate_wider_delim(mod_type,names = c("var_set","model"),delim = "_") |>
                   group_by(var_set,model) |>
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
  coord_cartesian(ylim = c(.85,1)) +
  facet_wrap(~model,ncol = 4)



### Bayesian GLM of performance ==========

#Do something like the perf_mod function in tidymodels
cv_results = 
  d_var_select_results |>
  select(-.mod_posterior) |>
  unnest(.results) |>
  select(id,id2,wflow_id,best_roc_auc) |> 
  separate_wider_delim(wflow_id,delim = "_",names = c("sampling","var_set","model")) |>
  mutate(id = interaction(id,id2,sep = "_")) |>
  mutate(id = paste(id,"v2",sep = "_")) |>
  select(-c(id2))

#Add the data from the all variable models
cv_results = 
  cv_results |>
  bind_rows(d_fold_results |>
              select(-.mod_posterior) |>  
              unnest(.results) |>
              select(id,id2,wflow_id,best_roc_auc) |> 
              mutate(wflow_id = str_replace(wflow_id,"impute","simple")) |>
              separate_wider_delim(wflow_id,delim = "_",names = c("sampling","model")) |>
              mutate(var_set = "all") |>
              mutate(id = interaction(id,id2,sep = "_")) |>
              select(-id2)) |>
  
  filter(sampling == "simple") |>
  select(-sampling) |>
  mutate(wflow_id = paste(var_set,model,sep = "_"))

#We can directly compare with tidyposterior, which takes the best submodel from
#each workflow and then does a stan_glm on the performance metrics
cv_mod <- 
  stan_glmer(best_roc_auc ~ 0 + wflow_id + (1|id),
             data = cv_results ,
             seed = 1, refresh = 0,
             chains = 4, cores = 8,
             iter = 10000, warmup = 2000,
             prior = rstanarm::normal(0,1),
             prior_intercept = rstanarm::normal(0,1.5),
             prior_aux = rstanarm::exponential(rate = 1))

# cv_mod |> summary()


### Figure 2: Top panel ======

#Lets make this into a plot
cv_mod_post = 
  cv_mod |>
  broom.mixed::tidy( conf.int = T, conf.level = 0.95, conf.method = "HPDinterval") |>
  mutate(term = str_remove(term,"wflow_id")) |>
  separate_wider_delim(term,delim = "_",names = c("var_set","model")) |>
  mutate(var_set = str_remove(var_set,".vars"))


p_compare_vars = 
  cv_mod_post |>
  mutate(model = case_when(model == "en" ~ "Penalised LR",
                           model == "nnet" ~ "ANN",
                           model == "rf" ~ "Random Forest",
                           model == "svm" ~ "RBF SVM"),
         var_set = case_when(var_set == "en" ~ "Penalised LR Variables",
                             var_set == "nnet" ~ "ANN Variables",
                             var_set == "rf" ~ "Random Forest Variables",
                             var_set == "svm" ~ "RBF SVM Variables",
                             var_set == "all" ~ "All Variables",
                             var_set == "moreThanOne" ~ "> 1 Model Variables",
                             var_set == "max" ~ "> 2 Models Variables")) |>
  mutate(model = factor(model,levels = c("ANN","Penalised LR","Random Forest","RBF SVM"))) |>
  mutate(var_set = factor(var_set,
                          levels = c("All Variables",
                                     "> 1 Model Variables",
                                     "> 2 Models Variables",
                                     "ANN Variables",
                                     "Penalised LR Variables",
                                     "Random Forest Variables",
                                     "RBF SVM Variables"))) |>
  ggplot(aes(x = var_set,
             y = estimate,
             ymin =  conf.low,ymax = conf.high,
             colour = model,
             fill   = model)) +
  geom_errorbar(position = position_dodge(width = 0.5),linewidth = 0.5) +
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
  facet_wrap(~model,ncol = 7) +
  geom_hline(yintercept = 1,lty = 2)

p_compare_vars


#Save this chunk
ggsave("C://Users/nadon/OneDrive - University of Bristol/Documents/CNV Item Reduction/Figures/Figure Parts/figure_2_top.pdf",
       p_compare_vars,width = 5, height = 2.5)


# Model Comparison =====


#Select the best performing variable set for each model type
cv_mod_post |> 
  group_by(model) |> 
  slice_max(estimate,n=1) |>
  mutate(across(where(is.double),~round(.x,digits = 3))) |>
  ungroup() |>
  arrange(-estimate) |>
  transmute(Model = model,
            `Variable Set` = var_set, 
            Performance = paste(estimate,"[",conf.low,",",conf.high,"]", sep = " ")) |>
  print(n = 30)


#So it looks like the subset of variables that gives best performance is the rf variables
#consistently across models


#You can do model comparison thus:
cv_mod_compare <- 
  cv_results |> 
  mutate(wflow_id = paste(var_set,model,sep = "_")) |>
  mutate(wflow_id = relevel(factor(wflow_id),ref = "rf_rf")) %>%
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
  mutate(Parameter = str_remove(Parameter,"wflow_id")) |>
  mutate(Parameter = str_remove(Parameter,".vars")) |>
  mutate(Parameter = case_when(Parameter == "(Intercept)" ~ "rf_rf",
                               TRUE ~ Parameter)) |>
  separate_wider_delim(Parameter,delim = "_",names = c("var_set","model")) |>
  
  mutate(across(where(is.double),~round(.x,digits = 3))) |>
  transmute(Model = model,
            `Variable Set` = var_set,
            Difference = paste(Median, "[",CI_low, ",",CI_high,"]",sep = " "),
            `Probability of Direction` = bayestestR::pd_to_p(pd))

tab_mod_compare$Difference[1] = "-"
tab_mod_compare$`Probability of Direction`[1] = 1



tab_mod_compare |>
  filter(`Probability of Direction` < 0.05) |>
  arrange(`Probability of Direction`) |>
  print(n = 30)


### Supplementary Table 5 ======

#Prepare the table
left_join(tab_mod_compare,
          cv_mod_post |> 
            mutate(across(where(is.double),~round(.x,digits = 3))) |>
            ungroup() |>
            transmute(Model = model,
                      `Variable Set` = var_set, 
                      Performance = paste(estimate,"[",conf.low,",",conf.high,"]", sep = " ")) |>
            print(n = 30), 
          by  = c("Model","Variable Set"))  |>
  relocate(Model, `Variable Set`, 
           Performance,Difference,`Probability of Direction`) |>
  mutate(Model = case_when(Model == "en" ~ "Penalised LR",
                           Model == "nnet" ~ "ANN",
                           Model == "rf" ~ "Random Forest",
                           Model == "svm" ~ "RBF SVM"),
         `Variable Set` = case_when(`Variable Set` == "all" ~ "All Variables",
                                    `Variable Set` == "svm" ~ "RBF SVM Variables",
                                    `Variable Set` == "en" ~ "Penalised LR Variables",
                                    `Variable Set` == "nnet" ~ "ANN Variables",
                                    `Variable Set` == "rf" ~ "Random Forest Variables",
                                    `Variable Set` == "max" ~ "> 2 Models Variables",
                                    `Variable Set` == "moreThanOne" ~ "> 1 Model Variables")) |>
  
  arrange(Model,`Variable Set`) |>
  
  #Now make it nice and print it to the viewer
  knitr::kable(format = "html", booktabs = TRUE) |>
  kableExtra::kable_styling(font_size = 11)

#And then we do a bit of re-ordering of the table in excel/word to make the publication version

#This time round we are seeing its the random forest model with the random forest variable
#set that has the best performance

#Given all this, which variable set do we use? I think it is keeping with the philosophy
#of the paper so far that we use the RF variable set


  
# Multiple Performance Metrics =======

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
  bind_rows(
    d_fold_results |>
      select(-.mod_posterior ) |>
      unnest(.results) |>
      separate_wider_delim(wflow_id,names = c("sampling","model"),delim = "_") |>
      mutate(wflow_id = paste(sampling,"all",model,sep = "_")) |>
      select(-c(sampling,model)),
    d_var_select_results |>
      select(-.mod_posterior) |>
      unnest(.results) ) |>
  select(id,id2,wflow_id,best_roc_auc,outer_full_fit) |>
  mutate(best_brier = map_dbl(outer_full_fit,~.x |> pluck('.predictions',1) |>
                                brier_score('group','.pred_ND-GC'))) |>
  select(-outer_full_fit) |>
  separate_wider_delim(wflow_id,names = c("sampling","var_set","model"),delim = "_") |>
  mutate(id = interaction(id,id2,sep = "_")) |>
  select(-c(id2)) 
  
  
d_metrics = 
  d_metrics |>
  filter(sampling == "impute" | sampling == "simple") |>
  select(-sampling) |>
  mutate(wflow = interaction(var_set,model,sep = "_"))
  
  
  r_model <- bf(best_roc_auc ~ 0 + wflow + (1|id))
  b_model <- bf(  best_brier ~ 0 + wflow + (1|id))
  
  
  get_prior(formula = r_model + b_model + set_rescor(TRUE),
            data = d_metrics)
  
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
      
    
      #Brier model
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
    mutate(name = str_remove(name, "wflow")) %>%
    separate_wider_delim(name,names = c("outcome", "var_set","model"),delim = "_") |>

    mutate(model = case_when(model == "en" ~ "Penalised LR",
                             model == "nnet" ~ "ANN",
                             model == "rf" ~ "Random Forest",
                             model == "svm" ~ "RBF SVM"),
           var_set = case_when(var_set == "en" ~ "Penalised LR Variables",
                               var_set == "nnet" ~ "ANN Variables",
                               var_set == "rf" ~ "Random Forest Variables",
                               var_set == "svm" ~ "RBF SVM Variables",
                               var_set == "all" ~ "All Variables",
                               var_set == "moreThanOne" ~ "> 1 Model Variables",
                               var_set == "max" ~ "> 2 Models Variables")) |>
    mutate(model = factor(model,levels = c("ANN","Penalised LR","Random Forest","RBF SVM"))) |>
    mutate(var_set = factor(var_set,
                            levels = c("All Variables",
                                       "> 1 Model Variables",
                                       "> 2 Models Variables",
                                       "ANN Variables",
                                       "Penalised LR Variables",
                                       "Random Forest Variables",
                                       "RBF SVM Variables"))) |>
    mutate(outcome = case_when(outcome == "rocauc" ~ "ROC AUC",
                               outcome == "brier" ~ "Brier Score")) |>
    mutate(outcome = factor(outcome,levels = c("ROC AUC","Brier Score"))) |>
    ggplot(aes(y = value, 
               x = model,
               group = var_set,
               colour = var_set)) +
    stat_pointinterval(point_interval = median_hdi, .width = c(0.66,0.95), position = position_dodge(width = 0.7)) +
    geom_vline(xintercept = 0, colour = "black", linetype = 3) +
    labs(subtitle = "BRMS Multivariate Model - All Variable Sets")+
    # coord_cartesian(xlim = c(-.5,.75)) +
    scale_color_brewer(palette = 1) +
    theme_bw() +
    facet_wrap(~outcome,scales = "free_y")
  
  mv_p1
  

  #Tabulate our data
  fitM %>%
    bayestestR::describe_posterior(ci = 0.95,
                       ci_method = "hdi",
                       centrality = "median",
                       test = c("p_direction","p_map")) %>%
    as_tibble() %>%
    mutate(pd_p = bayestestR::pd_to_p(pd,direction = "two-sided"),
           Parameter = str_remove(Parameter, "wflow_id"),
           Parameter = str_remove(Parameter, "b_"),) %>%
    separate_wider_delim(Parameter,names = c("outcome", "var_set","model"),delim = "_") |>
    select(-c(CI,Rhat, ESS)) %>%
    mutate(across(where(is.double),~round(.x,  digits = 3))) %>%
    mutate(Performance = paste(Median,"[",CI_low,",",CI_high,"]", sep = " "))  |>
    select(outcome,model,var_set,Performance) |>
    mutate(model = case_when(model == "en" ~ "Penalised LR",
                             model == "nnet" ~ "ANN",
                             model == "rf" ~ "Random Forest",
                             model == "svm" ~ "RBF SVM"),
           var_set = case_when(var_set == "en" ~ "Penalised LR Variables",
                               var_set == "nnet" ~ "ANN Variables",
                               var_set == "rf" ~ "Random Forest Variables",
                               var_set == "svm" ~ "RBF SVM Variables",
                               var_set == "all" ~ "All Variables",
                               var_set == "moreThanOne" ~ "> 1 Model Variables",
                               var_set == "max" ~ "> 2 Models Variables")) |>
    mutate(model = factor(model,levels = c("ANN","Penalised LR","Random Forest","RBF SVM"))) |>
    mutate(var_set = factor(var_set,
                            levels = c("All Variables",
                                       "> 1 Model Variables",
                                       "> 2 Models Variables",
                                       "ANN Variables",
                                       "Penalised LR Variables",
                                       "Random Forest Variables",
                                       "RBF SVM Variables"))) |>
    pivot_wider(names_from = "outcome",values_from = "Performance") |>
    knitr::kable(format = "html", booktabs = TRUE) |>
    kableExtra::kable_styling(font_size = 11)
  

  
  
# Up and Down Sampling ======

#One for the supplement here - within our selected model and variable set, do up and down 
#sampling during model fitting influence model performance? 

cv_results_sampling = 
    bind_rows(
      d_fold_results |>
        select(-.mod_posterior ) |>
        unnest(.results) |>
        separate_wider_delim(wflow_id,names = c("sampling","model"),delim = "_") |>
        mutate(wflow_id = paste(sampling,"all",model,sep = "_")) |>
        select(-c(sampling,model)),
      d_var_select_results |>
        select(-.mod_posterior) |>
        unnest(.results) ) |>
    select(id,id2,wflow_id,best_roc_auc,outer_full_fit) |>
    mutate(best_brier = map_dbl(outer_full_fit,~.x |> pluck('.predictions',1) |>
                                  brier_score('group','.pred_ND-GC'))) |>
    select(-outer_full_fit) |>
    separate_wider_delim(wflow_id,names = c("sampling","var_set","model"),delim = "_") |>
    mutate(id = interaction(id,id2,sep = "_")) |>
    select(-c(id2)) |>
    
    mutate(sampling = case_when(sampling == "impute" ~ "simple",
                                TRUE ~ sampling)) |>
    
    mutate(wflow_id = interaction(sampling,var_set,model,sep = "_"))
  

#Fit Model
cv_mod_sampling <- 
  stan_glmer(best_roc_auc ~ 0 + wflow_id + (1|id),
             data = cv_results_sampling,
             seed = 1, refresh = 0,
             chains = 4, cores = 8,
             iter = 10000, warmup = 2000,
             prior = rstanarm::normal(0,1),
             prior_intercept = rstanarm::normal(0,1.5),
             prior_aux = rstanarm::exponential(rate = 1))


#Make our own plot
cv_mod_post_sampling = 
  cv_mod_sampling |>
  broom.mixed::tidy( conf.int = T, conf.level = 0.95, conf.method = "HPDinterval") |>
  mutate(term = str_remove(term,"wflow_id")) |>
  separate_wider_delim(term,names = c("sampling","var_set","model"),delim = "_") 
  
  
#Tabulate
tab_post_sampling = 
  cv_mod_post_sampling |>
  mutate(model = case_when(model == "en" ~ "Penalised LR",
                           model == "nnet" ~ "ANN",
                           model == "rf" ~ "Random Forest",
                           model == "svm" ~ "RBF SVM"),
         var_set = case_when(var_set == "en" ~ "Penalised LR Variables",
                             var_set == "nnet" ~ "ANN Variables",
                             var_set == "rf" ~ "Random Forest Variables",
                             var_set == "svm" ~ "RBF SVM Variables",
                             var_set == "all" ~ "All Variables",
                             var_set == "moreThanOne" ~ "> 1 Model Variables",
                             var_set == "max" ~ "> 2 Models Variables"),
         sampling = case_when(sampling == "simple" ~ "Simple",
                              sampling == "down" ~ "Down",
                              sampling == "up" ~ "Up")) |>
  mutate(model = factor(model,levels = c("ANN","Penalised LR","Random Forest","RBF SVM"))) |>
  mutate(var_set = factor(var_set,
                          levels = c("All Variables",
                                     "> 1 Model Variables",
                                     "> 2 Models Variables",
                                     "ANN Variables",
                                     "Penalised LR Variables",
                                     "Random Forest Variables",
                                     "RBF SVM Variables"))) |>
  mutate(sampling = factor(sampling, levels = c("Simple","Down","Up"))) |>
  mutate(across(where(is.double),~round(.x,digits = 3))) |>
  arrange(model) |>
  transmute(Sampling = sampling,
            `Variable Set` = var_set,
            Model = model,
            Performance = paste(estimate, "[",conf.low, ",",conf.high,"]",sep = " "))
  
  
tab_post_sampling |>
  knitr::kable(format = "html", booktabs = TRUE) |>
  kableExtra::kable_styling(font_size = 11)