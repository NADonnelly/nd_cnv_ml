
# Introduction ------

# Here we are going to take our final models and look at how to deploy them e.g. as a shiny app
pacman::p_load(tidyverse,
               tidymodels,
               vetiver,
               crayon,
               patchwork)

tidymodels_prefer()

`%nin%` = Negate(`%in%`)


#Set our data directory
if(str_detect(Sys.info()[['nodename']],'IT088825')){
  data_dir = 'C://Users/pyxnd/OneDrive - University of Bristol/Documents/CNV Item Reduction/Data'
  
}else if(str_detect(Sys.info()[['nodename']],'AVOCADO')){
  
  data_dir = "C://Users/nadon/OneDrive - University of Bristol/Documents/CNV Item Reduction/Data/"
}

setwd(data_dir)


#We load our final models


#These are the 30 variable models
d_var_select_results = read_rds("nested_cv_result_selected_vars_v2.rds")

#This is the 5 variable model
d_var_ega_results    = read_rds("nested_cv_result_ega_vars.rds")

#This is the test data
d_last_fit           = read_rds("nested_cv_imputed_data.rds")


# Extract best models =====

best_mods = 
  bind_rows(
    d_var_select_results |>
      select(-.mod_posterior) |>
      unnest(.results) ,
    d_var_ega_results |>
      select(-.mod_posterior) |>
      unnest(.results)) |> 
  filter(str_detect(wflow_id,"simple") | str_detect(wflow_id,"EGA") ) |>
  mutate(wflow_id = str_remove(wflow_id,"simple_")) |>
  separate_wider_delim(wflow_id,
           names = c("var_set","model"),
           delim = "_") |>
  group_by(var_set,model) |>
  slice_max(best_roc_auc) |>
  
  #filter so that we have the variable sets we wish to see
  filter((var_set %in% c("rf") & model %in% c("rf")) | (var_set %in% c("EGA") & model %in% c("en")) )


best_mods = 
  best_mods |>
  mutate(test_fit = purrr::map(best_wf_final, 
                               last_fit,
                               split = d_last_fit,
                               metrics = metric_set(accuracy,kap,mn_log_loss,roc_auc,gain_capture))) |>
  mutate(best_metrics = map(test_fit, ~.x$.metrics[[1]]))



best_mods |> 
  select(var_set,best_metrics) |> 
  unnest(best_metrics) |> 
  filter(.metric %in% c("roc_auc","mn_log_loss"))




# Predicted probability distribution =====

#What I would like to do is to calculate the predicted probability from the EGA model at each possible value
#that the inputs can take

#Get all combinations of inputs for the 5 variables - this is doable with only 5 variables, only one of which 
#ordinal
vars_ega = 
  best_mods |>
  filter(var_set == "EGA") |>
  pluck("test_fit",1) |> 
  extract_recipe() %>%
  .$var_info

d_ega = 
  d_last_fit |> 
  training() |>
  select(all_of(vars_ega$variable))

d_fit_ega = 
  d_ega |>
  select(-group) |>
  crossing()


#For the 30 items we can't really do all possible combinations 
#because that would be over 2^30, which is too much. But perhaps we
#can do random sampling from the combinations available?

vars_30 = 
  best_mods |>
  filter(var_set == "rf") |>
  pluck("test_fit",1) |> 
  extract_recipe() %>%
  .$var_info

d_30 = 
  d_last_fit |> 
  training() |>
  select(all_of(vars_30$variable))

d_fit_30 = 
  d_30 |>
  select(-group) |>
  pivot_longer(everything()) |>
  group_by(name ) |>
  nest() |>
  mutate(unique_var = map(data,unique)) |>
  select(-data) |>
  mutate(r_c = map(unique_var,~slice_sample(.x, n = 1e6, replace = T) |> pull(value))) |>
  select(-unique_var) |>
  pivot_wider(names_from = "name",values_from = "r_c") |>
  unnest(everything())



#Now we need to work out how to fit these

#Note in out final ML we found a threshold of 0.7625 to be optimal for classification j-index
p_ega = 
  bind_cols(
    d_fit_ega,
    best_mods$test_fit[[1]]$.workflow[[1]] |>
      predict(d_fit_ega,type = "class"),
    best_mods$test_fit[[1]]$.workflow[[1]] |>
      predict(d_fit_ega,type = "prob")
    ) |>
  mutate(total = rowSums(across(1:5))) |>
  group_by(total) |>
  summarise(var_min  = min(`.pred_ND-GC`),
            var_prop = median(`.pred_ND-GC`),
            var_max  = max(`.pred_ND-GC`)) |>
  ggplot(aes(x = total,y = var_prop)) +
  geom_col() +
  geom_point(size = 2) +
  geom_hline(yintercept = 0.7625,lty = 2) +
  geom_hline(yintercept = 1,lty = 1) +
  geom_linerange(aes(ymin = var_min,ymax = var_max)) +
  theme_bw() +
  labs(x = "Total Score", y = "Probability ND-GC Carrier")

#Calibration of the 5 item model sucks

#For the full 30 variable RF model a threshold of 0.835 was optimal
p_30 = 
  bind_cols(
    d_fit_30,
    best_mods$test_fit[[2]]$.workflow[[1]] |>
      predict(d_fit_30,type = "class"),
    best_mods$test_fit[[2]]$.workflow[[1]] |>
      predict(d_fit_30,type = "prob")) |>
  mutate(total = rowSums(across(CMD_2:P_Health_dev_9 ))) |>
  group_by(total) |>
  summarise(var_min  = min(`.pred_ND-GC`),
            var_prop = median(`.pred_ND-GC`),
            var_max  = max(`.pred_ND-GC`)) |>
  ggplot(aes(x = total,y = var_prop)) +
  geom_col() +
  geom_point(size = 2) +
  geom_hline(yintercept = 0.835,lty = 2) +
  geom_hline(yintercept = 1,lty = 1) +
  geom_linerange(aes(ymin = var_min,ymax = var_max)) +
  theme_bw() +
  labs(x = "Total Score", y = "Probability ND-GC Carrier")


p_ega|p_30


# Simplify best model =====

#Our best model is a brisk 8 megabytes, which seems weighty
#I suspect we could do better
lobstr::obj_size(best_mods$test_fit[[2]]$.workflow[[1]])

#Most of the size seems to come from the recipe steps?
butcher::weigh(best_mods$test_fit[[2]]$.workflow[[1]])

butchered_wf = 
  butcher::butcher(best_mods$test_fit[[2]]$.workflow[[1]])

#We see a nice reduction in object size
lobstr::obj_size(butchered_wf)

#Does the butchered workflow still work?

#Make some random test data
d_test_30 = 
  d_last_fit |> 
  training() |>
  select(all_of(vars_30$variable)) |>
  select(-group) |>
  pivot_longer(everything()) |>
  group_by(name ) |>
  nest() |>
  mutate(unique_var = map(data,unique)) |>
  select(-data) |>
  mutate(r_c = map(unique_var,~slice_sample(.x, n = 10, replace = T) |> pull(value))) |>
  select(-unique_var) |>
  pivot_wider(names_from = "name",values_from = "r_c") |>
  unnest(everything())

#It does seem to work
butchered_wf |>
  predict(d_test_30 |> slice(1),type = "prob") |>
  pull(`.pred_ND-GC`) |>
  round(digits = 3) %>%
  paste("predicted probability ND-GC = ",.,sep = "")

#Now we need to load it into a shiny app

#lets save our model
write_rds(butchered_wf,"C://Users/nadon/OneDrive - University of Bristol/Documents/CNV Item Reduction/CNV_ML/cnv_ml_app/rf_model_30.rds")

#Lets do the same for the 5 item model
butchered_wf_5 = 
  butcher::butcher(best_mods$test_fit[[1]]$.workflow[[1]])


write_rds(butchered_wf_5,"C://Users/nadon/OneDrive - University of Bristol/Documents/CNV Item Reduction/CNV_ML/cnv_ml_app/en_model_5.rds")



# Plain old LR ======

## 30 variable model =====

#What if we just fit a plain old logistic regression model? It would be simpler wouldn't it


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
        dplyr::select(group, all_of(best_mods$test_fit[[2]] |> 
                                      extract_recipe() %>%
                                      .$var_info |>
                                      pull(variable)))
      )


#Model output
tidy(fitted_logistic_model,exponentiate = TRUE,conf.int = T)   |>
  filter(term != "(Intercept)") |>
  arrange(p.value) |>
  filter(p.value < 0.05)


#Model coefs
fitted_logistic_model$fit$coefficients |>
  as_tibble(rownames = "term")


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

accuracy(cnv_results, truth = group,
                      .pred_class,event_level = "second")     
        
        

## 5 variable model =====


fitted_logistic_model_ega <- 
  
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
        select(group, all_of(best_mods$test_fit[[1]] |> 
                               extract_recipe() %>%
                               .$var_info %>% 
                               pull(variable)))
  )


#Model output
tidy(fitted_logistic_model_ega,exponentiate = TRUE,conf.int = T)   |>
  filter(term != "(Intercept)") |>
  arrange(p.value) 


#Model coefs
fitted_logistic_model_ega$fit$coefficients |>
  as_tibble(rownames = "term")


#Predicted classes
pred_class_ega <- 
  predict(fitted_logistic_model_ega,
          new_data = d_last_fit |> 
            testing(),
          type = "class")


# Prediction Probabilities
pred_proba_ega <- 
  predict(fitted_logistic_model_ega,
          new_data = d_last_fit |> testing(),
          type = "prob")

#Join them
cnv_results_ega <- 
  d_last_fit |> testing() %>%
  select(group) %>%
  bind_cols(pred_class_ega, pred_proba_ega)


#Various model metrics
conf_mat(cnv_results_ega, truth = group,
         estimate = .pred_class)


mcc(cnv_results_ega, truth = group,
    estimate = .pred_class)

mn_log_loss(cnv_results_ega, truth = group,
            estimate = `.pred_ND-GC`,event_level = "second")

roc_auc(cnv_results_ega, truth = group,
        estimate = `.pred_ND-GC`,event_level = "second")

accuracy(cnv_results_ega, truth = group,
         .pred_class,event_level = "second")     


## Compare LR models ======

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
                                c_m(group,`.pred_ND-GC`)))
  
  return(bs)
  
}



boot_mods_final = 
  best_mods |>
  rename(mod_type = model) |>
  select(mod_type,test_fit) |> 
  unnest(test_fit)|>
  mutate(boot_metrics = map(.workflow, ~bootstrap_wf_metrics(wf = .x,
                                                             data_set = d_last_fit |>
                                                               testing(),
                                                             n_boot = 999) ))


#Now apply this function to our data
boot_mods_lr = 
  tibble(mod_type = c("lr_30","lr_ega"),
         model = list(fitted_logistic_model,fitted_logistic_model_ega)) |>
  mutate(boot_metrics = map(model, ~bootstrap_wf_metrics(wf = .x,
                                                             data_set = d_last_fit |>
                                                               testing(),
                                                             n_boot = 999) ))



#Make a plot of our bootstrapped values

bind_rows(
  boot_mods_lr |>
    select(mod_type,boot_metrics),
  boot_mods_final |> 
    ungroup() |>
    mutate(mod_type = paste(mod_type,"_",var_set,sep="")) |>
    select(mod_type,boot_metrics)
  ) |>
  unnest(boot_metrics)|>
  select(-splits) |>
  unnest(boot_metrics) |>
  ggplot(aes(x = .estimate, y = mod_type)) +
  ggdist::stat_halfeye() +
  facet_wrap(~.metric,ncol = 1)


#Compare models based on the bootstrapped samples
boot_pd = 
  left_join(
    bind_rows(boot_mods_final |> mutate(mod_type = paste(mod_type,var_set,sep = "_")) |> ungroup() |> select(mod_type,boot_metrics) ,
              boot_mods_lr |> select(mod_type,boot_metrics)) |>
      unnest(boot_metrics)|>
      select(-splits) |>
      unnest(boot_metrics) |>
      filter(.metric == "roc_auc") |>
      pivot_wider(names_from = mod_type,values_from = .estimate) |>
      mutate(rf30_en5  = rf_rf - en_EGA,
             rf30_lr30 = rf_rf - lr_30,
             rf30_lr5  = rf_rf - lr_ega) |>
      select(rf30_en5,rf30_lr30,rf30_lr5) |>
      pivot_longer(everything()) |>
      group_by(name) |>
      ggdist::mean_hdci(value),
    
      bind_rows(boot_mods_final |> mutate(mod_type = paste(mod_type,var_set,sep = "_")) |> ungroup() |> select(mod_type,boot_metrics) ,
                boot_mods_lr |> select(mod_type,boot_metrics)) |>
      unnest(boot_metrics)|>
      select(-splits) |>
      unnest(boot_metrics) |>
      filter(.metric == "roc_auc") |>
      pivot_wider(names_from = mod_type,values_from = .estimate) |>
      mutate(rf30_en5  = rf_rf - en_EGA,
             rf30_lr30 = rf_rf - lr_30,
             rf30_lr5  = rf_rf - lr_ega) |>
      select(rf30_en5,rf30_lr30,rf30_lr5) |>
      pivot_longer(everything()) |>
      group_by(name) |>
      nest() |>
      mutate(pd = map(data, ~bayestestR::p_direction(.x) |> as.numeric() |> bayestestR::pd_to_p(direction = "two-sided"))) |>
      unnest(pd) |> select(-data),
    by = "name")


boot_pd = 
  boot_pd |>
  arrange(name) |>
  select(name,value,.lower,.upper,pd) |>
  mutate(across(where(is.double),round,digits = 3)) |>
  mutate(performance = paste(value,"[",.lower,",",.upper,"]",sep = " ")) |>
  select(name,performance,pd)


#Lets butcher and save our 30 item LR model
butchered_lr = 
  butcher::butcher(fitted_logistic_model)

write_rds(butchered_lr,"C://Users/nadon/OneDrive - University of Bristol/Documents/CNV Item Reduction/CNV_ML/cnv_ml_app//lr_model_30.rds")