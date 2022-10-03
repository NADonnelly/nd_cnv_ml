pacman::p_load(tabnet, tidymodels,modeldata)


set.seed(123)
data("lending_club", package = "modeldata")
split <- initial_split(lending_club, strata = Class)
train <- training(split)
test  <- testing(split)


rec <- recipe(Class ~ ., train) %>%
  step_normalize(all_numeric())


mod <- tabnet(epochs = 50, batch_size = 128) %>%
  set_engine("torch", verbose = TRUE) %>%
  set_mode("classification")


wf <- workflow() %>%
  add_model(mod) %>%
  add_recipe(rec)


folds <- vfold_cv(train, v = 5)


fit_rs <-
  wf %>%
  fit_resamples(folds)




#Experiment with tabnet


#We need to get imputed data ready. I think this can be done separately to make things work better
rec_0 = 
  d_outer |>
  analysis() |>
  recipe(group ~ .) |>
  step_impute_bag(all_predictors()) |>
  prep( )

d_tn = 
  rec_0 |>
  juice()


#We can make bootstraps from this dataset
d_tn_inner = 
  d_tn |>
  bootstraps(times = 25)


tn_rec <- 
  d_outer |>
  analysis() |>
  recipe(group ~ .) 

tn_spec <- 
  tabnet(epochs = 50, batch_size = 128) %>%
  set_engine("torch", verbose = TRUE) %>%
  set_mode("classification")

tn_wf <- 
  workflow() %>%
  add_model(tn_spec) %>%
  add_recipe(tn_rec)


tn_fit <-
  tn_wf %>%
  fit_resamples(d_tn_inner)

collect_metrics(tn_fit)

#Fit the model to the full training data for this outer loop
tn_model <- 
  tn_wf %>% 
  fit(d_tn)

#We can look at the variable importance
p_tv_vip = 
  tn_model |> 
  extract_fit_parsnip() |> 
  vip(num_features = 30,
      geom = "col")


tn_model |> 
  extract_fit_parsnip() |> 
  vip::vi()

#This doesn't seem to work tremendously well
tn_ex <- tabnet_explain(tn_model |> 
                            extract_fit_parsnip(), d_tn)



#apply to the test data for this loop
d_outer_test = 
  rec_0 |>
  bake(new_data = d_outer |> 
         assessment())

tn_test = 
  d_outer_test |>
  bind_cols(
    predict(tn_model, 
            d_outer_test , 
            type = "prob")
  ) 

tn_test %>% 
  roc_auc(group, .pred_Control)

#Performance is pretty good


# Lets fit with a raw tabnet model rather than a parsnip version so we can use  the pretraining ======



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

# Prepare our initial data splits
d_split = initial_split(DF1, strata = group,prop = 4/5)

#Get test and training data
d_train = training(d_split)
d_test  = testing(d_split)

#Prepare for (nested) cross validation
d_folds = vfold_cv(d_train, v = 10, repeats = 1) 



#Using the pretrain thing doesn't seem t oreally improve the fit
# mod <- 
#   tabnet_pretrain(tn_rec, d_impute, epochs = 50, valid_split = 0.2, batch_size = 128, verbose = TRUE)
# 
# 
# autoplot(mod)
# 
# model_fit <- 
#   tabnet_fit(tn_rec, d_impute , tabnet_model = mod, from_epoch=40, valid_split = 0.2, epochs = 50, verbose=TRUE)
# 
# 
# autoplot(model_fit) 
# 
# model_fit <- 
#   tabnet_fit(tn_rec, d_impute , tabnet_model = model_fit, from_epoch=50, epochs = 12, valid_split = 0.2, verbose=TRUE)
# 
# 
# d_outer_test %>% 
#   bind_cols(
#     predict(model_fit, d_outer_test, type = "prob")
#   ) %>% 
#   roc_auc(group, .pred_Control)




# Missing data?? =======



#Visualise missing data
visdat::vis_miss(d_train)

#So overall we get a global 2.2% missing

#Lets look at variable importance in the missing data

#This just errors, I don't know why
rec_missing <- 
  recipe(group ~ ., data = d_train) 

missing_pretrain <- tabnet_pretrain(rec_missing, data=d_train, epoch=50, 
                                    valid_split = 0.2, verbose=TRUE, batch=128, 
                                    pretraining_ratio=0.48)

autoplot(ames_missing_pretrain)
vip_color(ames_missing_pretrain, col_with_missings)




#Their own example doesn't work with missing data??

library(tidymodels, quietly = TRUE)
library(tabnet)
data("ames", package = "modeldata")
qplot(ames$Mas_Vnr_Area)


col_with_zero_as_na <- ames %>% 
  select(where(is.numeric)) %>% 
  select(matches("_SF|Area|Misc_Val|[Pp]orch$")) %>% 
  summarise_each(min) %>% 
  select_if(~.x==0) %>% 
  names()
ames_missing <- ames %>% mutate_at(col_with_zero_as_na, na_if, 0) %>% 
  mutate_at("Alley", na_if, "No_Alley_Access") %>% 
  mutate_at("Fence", na_if, "No_Fence") %>% 
  mutate_at(c("Garage_Cond", "Garage_Finish"), na_if, "No_Garage") %>% 
  mutate_at(c("Bsmt_Exposure", "BsmtFin_Type_1", "BsmtFin_Type_2"), na_if, "No_Basement")

visdat::vis_miss(ames_missing)


cat_emb_dim <- map_dbl(ames %>% select_if(is.factor), ~log2(nlevels(.x)) %>% round)



ames_missing_rec <- recipe(Sale_Price ~ ., data=ames_missing) %>% 
  step_normalize(all_numeric())
ames_missing_pretrain <- tabnet_pretrain(ames_missing_rec, data=ames_missing, epoch=50, 
                                         cat_emb_dim = cat_emb_dim,
                                         valid_split = 0.2, verbose=TRUE, batch=2930, 
                                         pretraining_ratio=0.37)
autoplot(ames_missing_pretrain)
vip_color(ames_missing_pretrain, col_with_missings)
