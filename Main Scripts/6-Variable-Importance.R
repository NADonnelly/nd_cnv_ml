
# Variable importance ====

pacman::p_load(tidyverse,tidymodels,
               gtsummary,patchwork,
               vip,readxl,lubridate, crayon,
               future,doFuture,themis,
               DALEX,DALEXtra)

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


#Load the previously prepared data
d_fold_results = read_rds("nested_cv_result_all_vars_v2.rds")
d_last_fit     = read_rds("nested_cv_imputed_data.rds")


#Select only the training data for this process
d_pred = d_last_fit |> training()



## Define a function for determining variable importance =====

get_variable_importance = function(cv_fold_results,cv_data_set,model_name,np = 100,nvar = 30){
  
  #Put in some crayon to print our status to the command line
  cat(bgWhite$blue$inverse$bold('Processing ' %+% model_name %+% '... \n'))
  
  #Extract the model with the best AUC for variable importance calculation
  var_final_importance = 
    cv_fold_results |>
    unnest(.results) |>
    filter(str_detect(wflow_id,"impute") & str_detect(wflow_id,model_name)) |>
    slice_max(best_roc_auc) |>
    mutate(fit_engine = map(outer_full_fit,~extract_fit_engine(.x) ))
  
  
  cat(bgWhite$red$inverse$bold('Doing the permutation ' %+% model_name %+% '... \n'))
  
  #make a DALEX model explainer
  model_explainer = 
    explain_tidymodels(model = var_final_importance$outer_full_fit[[1]]$.workflow[[1]],
                       data  = cv_data_set |> select(-group), 
                       y     = cv_data_set$group == "ND-GC",
                       verbose = F)
  
  #Do the permutation based variable importance (this is the slow bit)
  model_parts = model_parts(explainer = model_explainer,
                            type = "variable_importance",
                            B = np)
  
  
  cat(bgWhite$green$inverse$bold('Extracting variables ' %+% model_name %+% '... \n'))
  
  #Extract the nvar most important variables from our model
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
plan(multisession, workers = 8)

#warning - this is very slow
best_vars = 
  tibble(model_name   = c("en",
                          "rf",
                          "svm",
                          "nnet"),
         title_string = c("Penalised LR",
                          "Random Forest",
                          "RBF SVM",
                          "ANN")) |>
  mutate(var_set    = map(model_name, ~get_variable_importance(cv_fold_results = d_fold_results,
                                                               model_name      = .x,
                                                               cv_data_set     = d_pred,
                                                               np = 500,nvar = 30))) |>
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
    ggdist::mean_hdci(dropout_loss,.width = 0.5) |>
    slice_max(dropout_loss,n = 30,with_ties = F) |>
    mutate(variable = fct_reorder(variable,dropout_loss,max)) |>
    ggplot(aes(y = variable,x = dropout_loss,xmin = .lower,xmax = .upper)) +
    geom_col() +
    geom_point() +
    geom_linerange() +
    theme_bw() +
    theme(panel.grid = element_blank()) +
    labs(x = "Change in AUC after permutation", y = "Variable", title = title_string) +
    scale_x_continuous(trans = shift_trans(d = part_data |>
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


#Compare the penalised LR results with the direct model coefs
left_join(
  d_fold_results |>
    unnest(.results) |>
    filter(str_detect(wflow_id,"impute") & str_detect(wflow_id,'en')) |>
    slice_max(best_roc_auc) |>
    pluck('outer_full_fit',1) |>
    pluck('.workflow',1) |>
    extract_fit_parsnip() |>
    tidy() |>
    filter(term != "(Intercept)") |>
    rename(variable = term),
  best_vars |>
    filter(model_name == "en") |>
    pluck("permute_data",1) |> 
    filter(variable %nin% c("_baseline_","_full_model_")) |>
    group_by(variable) |>
    ggdist::mean_hdci(dropout_loss,.width = 0.95) ,
  by = "variable") |>
  mutate(estimate = abs(estimate)) |>
  
  filter(estimate != 0) |>
  # arrange(-dropout_loss)
  # 
  ggplot(aes(x = estimate,y = dropout_loss)) +
  geom_point() +
  geom_smooth(method = 'lm',formula = y ~ x)

#So they do correlate, but of course the comparisons are going to be 
#inexact because the coefficients refelct the nature of the data as well as the
#importance of that variable.



### Combine Variables =====

# Compare our variables
var_set_full = 
  best_vars |>
  select(model_name,variable_summary) |>
  unnest(variable_summary) |>
  rename(type = model_name, 
         variable = variable, 
         importance = dropout_loss) |> 
  select(type,variable,importance)

var_set = 
  var_set_full |>
  count(variable)


#Notew that only 21 variables have a non-zero dropout loss in the en model
#due to penalisation

#Make a plot
var_set   |>
  ggplot(aes(y = fct_reorder(variable,n,max), x = n)) +
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
  filter(n > 2)

#Now we have our sets of variables, lets do our ML again with these subsets
#but, based on a reviewer comment, let us keep to just elastic net regression


# Make formulas ======

#Convert the sets of variables we selected above into formulas and then into recipes
form_list = 
  bind_rows(var_set_full,
            vars_max |> 
              mutate(type = "max") |>
              rename(importance = n),
            vars_more_than_one|> 
              mutate(type = "moreThanOne") |>
              rename(importance = n)) |>
  group_by(type) |>
  nest() |>
  
  #Make our formulas
  mutate(formulas = map_chr(data,~paste("group ~",.x |> 
                                          pull(variable) |> 
                                          paste(collapse = " + ")))) |>
  select(-data) |>
  
  #Make recipes which we will append our formulas into
  mutate(recipes_simple = map(formulas, ~recipe(formula = as.formula(.x),
                                                data    = d_pred) |>
                                step_zv(all_predictors())),
         recipes_up     = map(formulas, ~recipe(formula = as.formula(.x),
                                                data    = d_pred) |>
                                step_zv(all_predictors())|>
                                step_upsample(group,
                                              over_ratio = 1,
                                              seed = 23022023,
                                              skip = TRUE)),
         recipes_down   = map(formulas, ~recipe(formula = as.formula(.x),
                                                data    = d_pred) |>
                                step_zv(all_predictors())|>
                                step_downsample(group,
                                                under_ratio = 1,
                                                seed = 23022023,
                                                skip = TRUE)))


form_list = 
  form_list |>
  select(-formulas) |>
  ungroup() |>
  pivot_longer(-type) |>
  mutate(type = paste(str_remove(name,"recipes_"),type,sep = "_"))

#Extract the recipes and convert them into a named list
final_recipes = 
  form_list |>
  pull(value) 

names(final_recipes) =  
  form_list |> 
  pull(type)


#Save these variable sets
write_rds(final_recipes,"C://Users/nadon/OneDrive - University of Bristol/Documents/CNV Item Reduction/Data/nested_cv_selected_vars.rds")
# final_recipes = read_rds("C://Users/nadon/OneDrive - University of Bristol/Documents/CNV Item Reduction/Data/nested_cv_selected_vars.rds")
