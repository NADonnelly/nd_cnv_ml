
# Introduction ------

# ML using the finalised selected variable datasets

pacman::p_load(tidyverse,tidymodels,workflowsets,tidyposterior,
               finetune,probably,doFuture,furrr,themis,
               rstanarm,ggforce,doParallel,patchwork,
               vip,readxl,lubridate, crayon,
               DALEXtra)

conflicted::conflicts_prefer(crayon::`%+%`)
conflicted::conflicts_prefer(dplyr::filter)

#Set some preferences
tidymodels_prefer()
`%nin%` = Negate(`%in%`)




#Set our data directory
if(str_detect(Sys.info()[['nodename']],'IT088825')){
  data_dir = 'C://Users/pyxnd/OneDrive - University of Bristol/Documents/CNV Item Reduction/Data'
  
}else if(str_detect(Sys.info()[['nodename']],'AVOCADO')){
  
  data_dir = "C://Users/nadon/OneDrive - University of Bristol/Documents/CNV Item Reduction/Data/"
}

setwd(data_dir)


# Load saved models and data ====== 

#We load up our testing and training data
d_last_fit = read_rds("nested_cv_imputed_data.rds")

#And load our models
d_var_select_results = read_rds("nested_cv_result_selected_vars_v2.rds")

#Load the demographic information for use in the stratification
d_demo = read_csv("DemographicTable.csv")

d_demo = 
  d_demo |>
  mutate(age_strata = rsample::make_strata(age,breaks = 4)) |>
  mutate(multi_demo = interaction(group,age_strata,gender)) 




d_split = read_rds("nested_cv_split_data.rds")
d_test_id = d_last_fit |> complement()

# Fit and evaluate models to full training data =====


# We fit our best models to the full imputed training set and evaluate on the imputed test data

#Select the best set of parameters for each model - we will use the rf
best_mods = 
  d_var_select_results |>
  select(-.mod_posterior) |>
  unnest(.results) |> 
  separate_wider_delim(wflow_id,delim = "_",names = c("sampling","var_set","model")) |>
  filter(sampling == "simple" & var_set == "rf") |>
  select(-c(sampling,var_set)) |>
  mutate(best_log_loss = map_dbl(outer_full_fit,~.x |> 
                                   select(.metrics) |>
                                   unnest(.metrics) |>
                                   filter(.metric == "mn_log_loss") |>
                                   pull(.estimate))) |>
  mutate(best_kap = map_dbl(outer_full_fit,~.x |> 
                              select(.metrics) |>
                              unnest(.metrics) |>
                              filter(.metric == "kap") |>
                              pull(.estimate))) |>
  group_by(model) |>
  slice_max(best_kap) 
  # slice_min(best_log_loss)
  

#Use last_fit to evaluate performance of the final models
best_mods = 
  best_mods |>
  mutate(test_fit = purrr::map(best_wf_final, ~last_fit(.x,
                                                        split = d_last_fit,
                                                        metrics = metric_set(kap,mn_log_loss,roc_auc,bal_accuracy)))) |>
  mutate(best_metrics = map(test_fit, ~.x$.metrics[[1]]))



# Performance measures =====


#Tabulate our final performance
tab_mods = 
  best_mods |>
  ungroup() |>
  
  #We can work out the Brier score (this seems to simply be the average difference between the predicted 
  #probability and true outcome squared i.e. the mean squared error) - there is also the mn log loss measure
  #in yardstick which is the same thing on a slightly different scale (the log scale)
  mutate(brier = map_dbl(test_fit,~ .x$.predictions[[1]] |>
                           transmute(pred = `.pred_ND-GC`,group) |>
                           mutate(group = ifelse(group == "ND-GC",1,0)) |>
                           mutate(dist  = (group - pred)^2) |>
                           summarise(brier = mean(dist)) |>
                           pull())) |>
  select(model,best_metrics,brier)|>
  unnest(best_metrics) |>
  select(model,.metric,.estimate,brier) |>
  pivot_wider(names_from = .metric,values_from = .estimate) |>
  arrange(-roc_auc,mn_log_loss) |>
  mutate(across(where(is.double),~round(.x,digits = 3))) |>
  mutate(model =  case_when(model == "en" ~ "Penalised LR",
                            model == "nnet" ~ "ANN",
                            model == "rf" ~ "Random Forest",
                            model == "svm" ~ "RBF SVM")) |>
  select(model,roc_auc,bal_accuracy,kap,brier) 


tab_mods |>
  knitr::kable(format = "html", booktabs = TRUE) |>
  kableExtra::kable_styling(font_size = 11)



tab_mods |>
  pivot_longer(roc_auc:brier) |>
  ggplot(aes(x = model,
             y = value,
             fill = name))+
  geom_col() +
  theme_bw() +
  theme(axis.text.x.bottom = element_text(angle = -90)) +
  facet_wrap(~name,nrow = 1,scales = "free_y")

#So our random forest model manages to pull off the highest ROC and lowest brier score


# Bootstrap confidence intervals =====

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
  ungroup() |>
  select(model,boot_metrics) |>
  unnest(boot_metrics) |>
  ggplot(aes(x = .estimate, y = model)) +
  ggdist::stat_halfeye() +
  facet_wrap(~.metric,ncol = 1, scale = "free_x") +
  theme_bw() +
  geom_vline(xintercept = c(0,1),lty = 1)

#The AUC measure looks very similar; however the models do differ rather on Brier score

#We can also do model comparison using the bootstrap samples
boot_pd = 
  left_join(boot_mods |>
              ungroup() |>
              select(model,boot_metrics) |>
              unnest(boot_metrics)|>
              filter(.metric == "roc_auc") |>
              pivot_wider(names_from = c(model),values_from = .estimate) |>
              mutate(rf_en   = rf - en,
                     rf_nnet = rf - nnet,
                     rf_svm   = rf -svm) |>
              select(rf_en ,rf_nnet ,rf_svm) |>
              pivot_longer(everything()) |>
              group_by(name) |>
              ggdist::mean_hdci(value),
            
            boot_mods |>
              ungroup() |>
              select(model,boot_metrics) |>
              unnest(boot_metrics)|>
              filter(.metric == "roc_auc") |>
              pivot_wider(names_from = c(model),values_from = .estimate) |>
              mutate(rf_en   = rf - en,
                     rf_nnet = rf - nnet,
                     rf_svm   = rf -svm) |>
              select(rf_en ,rf_nnet ,rf_svm) |>
              pivot_longer(everything()) |>
              group_by(name) |>
              nest() |>
              mutate(pd = map(data, ~bayestestR::p_direction(.x) |> 
                                as.numeric() |> 
                                bayestestR::pd_to_p(direction = "two-sided"))) |>
              unnest(pd) |> select(-data),
            by = "name")

boot_pd = 
  boot_pd |>
  mutate(name = factor(name, levels = c("rf_en","rf_nnet","rf_svm"))) |>
  arrange(name) |>
  select(name,value,.lower,.upper,pd) |>
  mutate(across(where(is.double),~round(.x,digits = 3))) |>
  mutate(performance = paste(value,"[",.lower,",",.upper,"]",sep = " ")) |>
  select(name,performance,pd)


## Table 4 ======

#And we can make this into a more pleasing tabular format
tab_boot = 
  boot_mods |>
  select(model, boot_metrics) |>
  unnest(boot_metrics)|>
  group_by(model,.metric) |>
  ggdist::mean_hdci(.estimate) |>
  arrange(.metric,-.estimate)


tab_boot |>
  mutate(across(where(is.double),~round(.x,digits = 3))) |>
  mutate(Performance = paste(.estimate,"[",.lower,",",.upper,"]",sep = " ")) |>
  select(model,.metric,Performance) |>
  mutate(model = case_when(model == "en" ~ "Penalised LR",
                           model == "nnet" ~ "ANN",
                           model == "rf" ~ "Random Forest",
                           model == "svm" ~ "RBF SVM")) |>  
  mutate(model = factor(model,levels = c("Random Forest","Penalised LR","ANN","RBF SVM"))) |>
  arrange(model) |>
  pivot_wider(names_from = .metric,values_from = Performance) |>
  mutate(roc_diff = c("-",boot_pd$performance),
         roc_pd   = c("-",boot_pd$pd)) |>
  knitr::kable(format = "html", booktabs = TRUE) |>
  kableExtra::kable_styling(font_size = 11) 



# Calibration measures ======

#Now we asses the calibratio nof our selected best model


## Performance by age =====

#What does performance look like with age?
tab_age = 
  bind_cols(
    d_demo |> 
      slice(d_test_id) |>
      select(-group),
    best_mods |> 
      filter(model == "rf") |> 
      pluck("test_fit",1) |> 
      pluck(".predictions",1)) |>
  select(-.config) |>
  group_by(age_strata)|>
  mutate(group = fct_rev(group),
         .pred_class = fct_rev(.pred_class)) |>
  yardstick::metrics(truth = group,
                     estimate=.pred_class, `.pred_ND-GC`) |>
  select(-.estimator) |>
  pivot_wider(names_from = .metric,values_from = .estimate) |>
  mutate(across(where(is.double),round,digits = 3)) |>
  select(-c(accuracy,kap)) |>
  relocate(age_strata,roc_auc,mn_log_loss) |>
  left_join(bind_cols(
    d_demo |> 
      slice(d_test_id) |>
      select(-group),
    best_mods |> 
      filter(model == "en") |> 
      pluck("test_fit",1) |> 
      pluck(".predictions",1)) |>
      select(-.config) |>
      mutate(group = fct_rev(group),
             .pred_class = fct_rev(.pred_class)) |>
      group_by(age_strata)|>
      nest() |>
      mutate(brier = map_dbl(data,~brier_score(.x,"group",".pred_ND-GC"))) |>
      select(-data) |>
      ungroup(),
    by = "age_strata")

tab_age |>
  knitr::kable(format = "html", booktabs = TRUE) |>
  kableExtra::kable_styling(font_size = 11)


## Performance by gender ======
tab_gender = 
  bind_cols(
    d_demo |> 
      slice(d_test_id) |>
      select(-group),
    best_mods |> 
      filter(model == "en") |> 
      pluck("test_fit",1) |> 
      pluck(".predictions",1)) |>
  select(-.config) |>
  group_by(gender)|>
  mutate(group = fct_rev(group),
         .pred_class = fct_rev(.pred_class)) |>
  yardstick::metrics(truth = group,
                     estimate=.pred_class, `.pred_ND-GC`) |>
  select(-.estimator) |>
  pivot_wider(names_from = .metric,values_from = .estimate) |>
  mutate(across(where(is.double),round,digits = 3)) |>
  select(-c(accuracy,kap)) |>
  relocate(gender,roc_auc,mn_log_loss) |>
  left_join(bind_cols(
    d_demo |> 
      slice(d_test_id) |>
      select(-group),
    best_mods |> 
      filter(model == "en") |> 
      pluck("test_fit",1) |> 
      pluck(".predictions",1)) |>
      select(-.config) |>
      mutate(group = fct_rev(group),
             .pred_class = fct_rev(.pred_class)) |>
      group_by(gender)|>
      nest() |>
      mutate(brier = map_dbl(data,~brier_score(.x,"group",".pred_ND-GC"))) |>
      select(-data) |>
      ungroup(),
    by = "gender")

tab_gender |>
  knitr::kable(format = "html", booktabs = TRUE) |>
  kableExtra::kable_styling(font_size = 11)


### Table 5 ======

#Tabulate our splits
bind_rows(
  tab_age |>
    mutate(Covariate = "Age",
           age_strata = as.character(age_strata)) |>
    rename(`Covariate Value` = age_strata),
  tab_gender |>
    mutate(Covariate = "Gender") |>
    rename(`Covariate Value` = gender)) |>
  relocate(Covariate, `Covariate Value`)  |>
  knitr::kable(format = "html", booktabs = TRUE) |>
  kableExtra::kable_styling(font_size = 11)



#Plots =====


#Plot this stuff
p_mods =
  tab_mods |>
  ungroup() |>
  pivot_longer(roc_auc:brier,names_to = ".metric",values_to = ".estimate") |>
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
  pivot_longer(-age_strata,names_to = ".metric",values_to = ".estimate") |>
  mutate(age_m = map_chr(age_strata, ~as.character(.x) %>% gsub("\\)|\\]|\\(|\\[", "", .))) |>
  separate(age_m,sep = ",",into = c("l","u"),convert = T) |>
  mutate(mid = (l+u)/2) |>
  select(mid,.metric,.estimate) |>
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
  pivot_longer(-gender,names_to = ".metric",values_to = ".estimate") |>
  ggplot(aes(x = gender, y = .estimate)) +
  geom_point() +
  facet_wrap(~.metric,nrow = 1) +
  geom_segment(yend = 0,aes(xend = gender)) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        strip.background = element_blank()) +
  labs(x = "gender")

p_mods/p_age/p_gender + plot_annotation(tag_levels = "A")



## Model calibration =====

d_cal = 
  best_mods |>
  ungroup() |>
  mutate(calPlotData = map(test_fit, ~.x$.predictions[[1]] %>%
                             transmute(pred = `.pred_ND-GC`,group) %>%
                             mutate(group = fct_rev(group)) %>%
                             caret::calibration(group ~ pred, data = .,cuts = 10) %>%
                             .$data %>%
                             as_tibble()))

p_cal = 
  d_cal |> 
  select(model,calPlotData) |> 
  unnest(calPlotData)  |>
  ggplot(aes(x = midpoint,y = Percent,ymin = Lower,ymax = Upper)) +
  geom_point(size = 0.4) +
  geom_linerange(linewidth = .25) +
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
  facet_wrap(~model,ncol =1, strip.position="left")



## Full ROC curves ------
p_roc = 
  best_mods |>
  ungroup() |>
  mutate(roc_curve_d = map(test_fit,~ .x$.predictions[[1]] |>
                             yardstick::roc_curve(group,`.pred_ND-GC`,event_level = "second"))) |>
  select(model,roc_curve_d) |>
  unnest(roc_curve_d) |>
  ggplot(aes(x = 1 - specificity, y = sensitivity,colour = model)) +
  geom_path(size = 0.25) +
  geom_abline(lty = 3) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.title = element_blank(),
        legend.position = "bottom",
        legend.text = element_text(size = 6),
        axis.title = element_text(size = 6),
        axis.text = element_text(size = 6))


## Thresholds ======
threshold_data = 
  best_mods |>
  ungroup() |>
  filter(model == "rf" ) |>
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
        legend.position = "bottom") +
  coord_cartesian(xlim = c(0,1))


## Histogram =====

p_dens =
  best_mods |>
  ungroup() |>
  select(model,test_fit)  |>
  unnest(test_fit) |>
  select(model,.predictions) |>
  unnest(.predictions) |>
  filter(model == "rf") |>
  ggplot(aes(`.pred_ND-GC`,colour = group,fill = group)) +
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
  labs(x = "Predicted Probability ND-GC",
       y = "Count") +
  coord_cartesian(xlim = c(0,1))

(p_dens/p_thres + plot_layout(ncol = 1,heights = c(1,4)))




#Now we have an optimal threshold we can do a confusion matrix

#What are the actual numbers of test participants classified
best_mods |>
  select(test_fit) |>
  ungroup() |>
  unnest(test_fit) |> 
  select(model,.predictions) |>
  unnest(.predictions) |>
  filter(model == "rf") |>
  mutate(.pred_class_new = if_else(`.pred_ND-GC` > max_j_index_threshold,"ND-GC","Control")) |>
  mutate(.pred_class_new = factor(.pred_class_new)) |>
  kap(group,.pred_class_new)
  bal_accuracy(group,.pred_class_new)
  conf_mat(group,.pred_class)

#Note here the threshold used for classification is max_j_index_threshold


  

  

# Final variable Importance ======

##Model Based ======
  
#Use model based VI for the RF model
best_mods |> 
  select(model,test_fit) |> 
  filter(model == "rf") |>  
  unnest(test_fit) |>
  select(model,.workflow) |>
  mutate(vi = map(.workflow,~.x |> 
                    extract_fit_engine() |>
                    vip::vi_model())) |>
  select(-.workflow) |>
  unnest(vi) |>
  mutate(Variable = fct_reorder(Variable, Importance,.fun = max)) |>
  ggplot(aes(x = Variable,y = Importance)) +
  geom_col(position = position_dodge()) +
  theme_bw() +
  theme(axis.text.x.bottom = element_text(angle = -90,
                                          vjust = 0)) +
    labs(y = "Model Based Variable Importance", x = "Variable") +
    scale_x_discrete(limits = rev)


## Permutation ------
  
#We could look at the importance of each of the variables in our final model
#using the same permutation based approach we did for the all variable models
#Fortunately this runs much faster with only the 30 variables, maybe about 15 
#minutes
  
best_model_importance = 
  best_mods |>
  filter(model == "rf") |>
  mutate(fit_engine = map(outer_full_fit,~extract_fit_engine(.x) ))

best_model_varset = 
  best_model_importance |>
  pluck('best_wf_final',1) |>
  extract_preprocessor() %>% 
  .$var_info %>% 
  filter(role == "predictor") %>% 
  pull(variable)

best_model_explainer = 
  explain_tidymodels(model = best_model_importance$outer_full_fit[[1]]$.workflow[[1]],
                     data = d_last_fit |> analysis() |> select(best_model_varset), 
                     y    = d_last_fit %>% analysis() %>% .$group == "ND-GC")


#best_model_performance = model_performance(best_model_explainer)

#Permutation will involve random numbers, so set a seed for reproducibility
set.seed(22072022)

best_model_parts = model_parts(explainer = best_model_explainer,
                               type = "variable_importance",
                               B = 500)


# Save this given it takes a long time
write_rds(best_model_parts, "final_variable_importance.rds")
# best_model_parts = read_rds("final_variable_importance.rds")


#Get the names for our best variables (this is a bit circular as the definitions are made elsewhere)
d_var = read_csv("nested_cv_selected_var_definitions_expanded.csv")


#Tabulate our best variables
best_model_table = 
  best_model_parts |>
  as_tibble() |> 
  group_by(variable) |>
  ggdist::median_hdci(dropout_loss,.width = 0.95)|>
  arrange(-dropout_loss) 


best_model_table |>
  filter(variable %nin% c("_baseline_","_full_model_")) |>
  left_join(d_var |> select(variable,short_name,var_def_long),by = "variable") |>
  mutate(across(where(is.double),~round(.x,digits = 3))) |>
  knitr::kable(format = "html", booktabs = TRUE) |>
  kableExtra::kable_styling(font_size = 11)


p_best_var_imp = 
  best_model_table |>
  filter(variable %nin% c("_baseline_","_full_model_")) |>
  # left_join(d_var |> select(variable,short_name,var_def_long),by = "variable") |>
  arrange(-dropout_loss) |>
  ggplot(aes(y = dropout_loss, x = forcats::fct_reorder(variable,dropout_loss,min),
             ymin = .lower, ymax = .upper)) +
  geom_point(size = 0.5) +
  geom_linerange(size = 0.25) +
  geom_hline(yintercept = best_model_table |> filter(variable == "_full_model_") |> pull(dropout_loss),
             lty = 2,size = 0.25) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.title = element_blank(),
        legend.position = "top",
        legend.text = element_text(size = 6),
        axis.title = element_text(size = 6),
        axis.text = element_text(size = 6),
        axis.text.x.bottom = element_text(angle = 90)) +
  labs(y = "Dropout Loss (1 - AUROC)", x = "Variable",
       title = "Best SVM Model",subtitle = "Vertical Line = Full Model") +
  scale_x_discrete(limits = rev)

p_best_var_imp


#Out of interest, how many variables have a dropout loss greater than the mean
best_model_parts |>
  as_tibble() |> 
  group_by(variable) |>
  ggdist::median_hdci(dropout_loss,.width = 0.95) %>%
  filter(dropout_loss >  filter(., variable == "_full_model_") |> 
           pull(dropout_loss) & 
           variable != "_baseline_") |>
  arrange(dropout_loss) |>
  print(n = 30)

#All but one



## Extract our final variables =====

#Lets get this set of variables:
best_vars = 
  best_mods$test_fit[[3]] |> 
  extract_recipe() %>%
  .$var_info |>
  pull(variable)

best_vars = 
  best_vars[!str_detect(best_vars,"group")]


#We can look these items up in the data dictionary....
VL <- readxl::read_excel("DATA ENTRY CODEBOOK .xlsx", 
                         sheet = "VARIABLE LIST")

#Wrangle lightly
VL = 
  VL |>
  select(VARIABLE:`VARIABLE DEFINITION`) |>
  filter(VARIABLE %in% best_vars)

#Glue on the final  variable importance
VL =
  VL |>
  janitor::clean_names() |>
  left_join(best_model_table |>
              select(-c(.width)), 
            by = "variable")

#And here we are
VL |>
  arrange(-dropout_loss ) |>
  print(n = 30)

#Save this
write_rds(VL,"nested_cv_selected_var_definitions.rds")





# Figure 3 Assembly: Performance Plot ======

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

ggsave("C://Users/nadon/OneDrive - University of Bristol/Documents/CNV Item Reduction/Figures/Figure Parts/figure_3.pdf",p_perf,width = 5, height = 3)
ggsave("C://Users/nadon/OneDrive - University of Bristol/Documents/CNV Item Reduction/Figures/Figure Parts/figure_3_vi.pdf",p_best_var_imp + plot_annotation(tag_levels = "A"),width = 5, height = 2)
#ggsave("C://Users/nadon/OneDrive - University of Bristol/Documents/CNV Item Reduction/Figures/Figure Parts/figure_3_compare_vars.pdf",p_compare_vars + plot_annotation(tag_levels = "A"),width = 5, height = 2.5)



# A final logistic regression model ====


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
        dplyr::select(group, all_of(best_mods$test_fit[[3]] |> 
                                      extract_recipe() %>%
                                      .$var_info |>
                                      pull(variable)))
  )



#Model summary table 
tidy(fitted_logistic_model,exponentiate = T,conf.int = T)   |>
  filter(term != "(Intercept)") |>
  arrange(p.value) |>
  print(n = 31)


#Model coefs
fitted_logistic_model$fit$coefficients |>
  as_tibble(rownames = "term") |>
  print(n = 31)


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
  geom_histogram(position = "dodge",alpha = 1,binwidth = 0.05,colour = "transparent",size = 0.01) +
  scale_fill_manual(values = c("#E1AF64","#32324B"))+
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.title = element_blank(),
        legend.position = "top",
        legend.text = element_text(size = 6),
        axis.title = element_text(size = 6),
        axis.text = element_text(size = 6)) +
  labs(x = "Predicted Probability ND-GC",
       y = "Count") +
  coord_cartesian(xlim = c(0,1))


cnv_results %>%
  caret::calibration(group ~ .pred_Control, data = .,cuts = 10) %>%
  .$data %>%
  as_tibble() |>
  ggplot(aes(x = midpoint,y = Percent,ymin = Lower,ymax = Upper)) +
  geom_point(size = 0.4) +
  geom_linerange(linewidth = .25) +
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
  coord_cartesian(ylim = c(0,100), xlim = c(0,100))
