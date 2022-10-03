# Introduction ======

#As script that does variable removal based on correlation/VIF

#All of our variables are ordinal now

# Load up =====

library(tidyverse)
library(car)

library(tidymodels)
library(embed)
library(vip)
library(learntidymodels)
library(probably)
library(finetune)

library(patchwork)

library(parallel)
library(doSNOW)
library(tcltk)


#Set a value for the variable correlation cutoff we will accept

#Are we actually going to use these though?
corr_thres = 0.5
vif_thres  = 100
  

### Read in filtered data

DF <- read_csv("FilteredData2.csv")

DF0 = 
  DF |>
  mutate(group = map_chr(GenotypeCode,~ifelse(.x == 16,"Control","ND-CNV"))) |>
  mutate(group = factor(group, levels = c("Control","ND-CNV"))) |>
  select(-c(IDs,GenotypeCode)) |>
  relocate(group)

#Have as quick look at all variables
# skimr::skim(DF0)

# One imputation for development ======

#Although this isn't ideal (you should really do imputation within CV loops),
#for the purpose of development I am going to do imputation here a single time 
#on all predictor variables

#To make this work properly we have to make sure the recipe understands which
#variables are integers

#Lets work out if a variable is all integers
int_vars = 
  DF0 |> 
  pivot_longer(-c(group)) |> 
  group_by(name) |> 
  nest() |>
  mutate(all_int = map_lgl(data, ~ .x |>
                             drop_na() |> 
                             mutate(frac = map_lgl(value,~.x%%1 == 0)) |> 
                             summarise(is_int = all(frac)) |>
                             pull())) |>
  select(name,all_int) |>
  filter(all_int == T)

DF1 = 
  DF0 |>
  mutate(across(int_vars |> pull(name), as.integer)) 

#Now we have done that, we can do the imputation, using bagged trees, which seems 
#like a reasonable place to start https://recipes.tidymodels.org/reference/step_impute_bag.html 
rec_0 = 
  DF1 |>
  recipe(group ~ .) %>%
  step_impute_bag(all_predictors()) %>%
  prep(training = DF1)

DF1 = 
  rec_0 |>
  juice()

#You could save this....
write_csv(DF1, "ImputedData2.csv")


# Regularisation with LASSO ====

DF1 = read_csv("ImputedData2.csv")

d_split = initial_split(DF1, strata = group,prop = 3/4)

d_train = training(d_split)
d_test  = training(d_split)

d_val   = validation_split(d_train,strata = group,prop = 4/5)

# #This does not work tremendously well due to the very large number of predictors...
# tmwr_cols <- colorRampPalette(c("#91CBD765", "#CA225E"))
# 
# d_train %>% 
#   select(-group) %>% 
#   cor() %>% 
#   corrplot(col = tmwr_cols(200), tl.col = "black", method = "ellipse")


## Initial dimensionality reduction ======

set.seed(19052022)

#Define a function
plot_validation_results <- 
  function(recipe, dat = assessment(d_val$splits[[1]])) {
    
  recipe %>%
      
    # Estimate any additional steps
    prep() %>%
      
    # Process the data (the validation set by default)
    bake(new_data = dat) %>%
      
    # Create the scatterplot matrix
    ggplot(aes(x = .panel_x, y = .panel_y, color = group, fill = group)) +
    geom_point(alpha = 0.8, size = 1) +
    geom_autodensity(alpha = .3) +
    facet_matrix(vars(-group), layer.diag = 2) + 
    scale_color_brewer(palette = "Dark2") + 
    scale_fill_brewer(palette  = "Dark2") + 
    theme_bw()
}


analysis(d_val$splits[[1]]) |>
  recipe(group ~ .) %>%
  step_pca(all_predictors(),num_comp = 4) %>%
  
  # Estimate any additional steps
  prep() %>%
  
  # Process the data (the validation set by default)
  bake(new_data = assessment(d_val$splits[[1]]))



#Do PCA
p_1 = 
  analysis(d_val$splits[[1]]) |>
  recipe(group ~ .) %>%
  step_pca(all_predictors(),num_comp = 4) %>%
  plot_validation_results  + 
  ggtitle("Principal Components Analysis")


## Plot top loadings
p_2 = 
  analysis(d_val$splits[[1]]) |>
  recipe(group ~ .)  %>%
  step_pca(all_numeric_predictors(), num_comp = 4) %>% 
  prep() %>% 
  plot_top_loadings(component_number <= 4, n = 5) + 
  scale_fill_brewer(palette = "Paired") +
  ggtitle("Principal Component Analysis") +
  theme_bw()


# Do PLS
p_3 = 
  analysis(d_val$splits[[1]]) |>
  recipe(group ~ .) %>%
  step_pls(all_predictors(),outcome = vars(group),num_comp = 4) %>%
  plot_validation_results  + 
  ggtitle("Partial Least Squares")


p_4 = 
  analysis(d_val$splits[[1]]) |>
  recipe(group ~ .)  %>%
  step_pls(all_predictors(),outcome = vars(group),num_comp = 4) %>%
  prep() %>% 
  plot_top_loadings(component_number <= 4, n = 5, type = "pls") + 
  scale_fill_brewer(palette = "Paired") +
  ggtitle("Partial Least Squares") +
  theme_bw()

#Do ICA
p_5 = 
  analysis(d_val$splits[[1]]) |>
  recipe(group ~ .)  %>%
  step_ica(all_numeric_predictors(), num_comp = 4) %>%
  plot_validation_results() + 
  ggtitle("Independent Component Analysis")

# Supervised UMAP
p_6 =
  analysis(d_val$splits[[1]]) |>
  recipe(group ~ .)  %>%
  step_umap(all_numeric_predictors(), num_comp = 4,outcome = "group") %>%
  plot_validation_results() +
  ggtitle("UMAP")


(p_1 | p_2 | p_5) / (p_3 | p_4 | p_6 )

#So this looks a lot like you can get a good split between controls and ND CNVs 
#using these dimension reducing approachs


## Lasso =====

#Lets do a lasso with tidymodels


#Make a recipe which includes an imputation step (lets not do that for the time being)
rec_2 =
  DF1 |>
  recipe(group ~ .) %>%
  step_zv(all_numeric(), -all_outcomes())# %>%
  # step_impute_bag(all_predictors())

rec_2_prep = 
  rec_2 |>
  prep() |> 
  juice()


#Specify a basic, and not optimised model
lasso_spec = 
  logistic_reg(mode = "classification",
               engine = "glmnet",
               penalty = 0.1, mixture = 1)

#Make a workflow
lasso_wf = 
  workflow() %>%
  add_recipe(rec_2)

#Fit the model
lasso_fit = 
  lasso_wf %>%
  add_model(lasso_spec) %>%
  fit(DF0)

#Inspect the model output
# lasso_fit %>%
#   extract_fit_parsnip() %>%
#   tidy()

  
#Lets tune our lasso penalty/lambda parameter
set.seed(1234)

lasso_boot = bootstraps(DF1,strata = group)

# tune_spec = 
#   logistic_reg(mode = "classification",
#                engine = "glmnet",
#                penalty = tune(), 
#                mixture = tune())

tune_spec = 
  logistic_reg(mode = "classification",
               engine = "glmnet",
               penalty = tune(), 
               mixture = 1)

#Tuning with grid search is slow so lets try an alternative with finetune, using 
#something called simulated annealing (whatever that means)

# lambda_grid <- 
#   grid_regular(penalty(), mixture(), levels = 30) 
# 
# 
# doParallel::registerDoParallel()
# 
# lasso_grid <- 
#   tune_grid(
#     lasso_wf %>% add_model(tune_spec),
#     resamples = lasso_boot,
#     grid = lambda_grid
# )

#Can we speed things up with some finetune?
ctrl <- control_sim_anneal(verbose = TRUE)

lasso_sa = 
  lasso_wf %>% 
  add_model(tune_spec) %>%
  tune_sim_anneal(resamples = lasso_boot, iter = 20, control = ctrl, metrics = metric_set(roc_auc))



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

### Get variable importance from our LASSO model ====

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


#We can look at performance at different probability levels
library(probably)

threshold_data = 
  lasso_final %>%
  augment(rec_2_prep ) |>
  select(group,.pred_class,.pred_Control)%>%
  threshold_perf(group, .pred_Control, thresholds = seq(0, 1, by = 0.0025))


threshold_data <- 
  threshold_data %>%
  filter(.metric != "distance") %>%
  mutate(group = case_when(
    .metric == "sens" | .metric == "spec" ~ "1",
    TRUE ~ "2"
  ))

max_j_index_threshold <- 
  threshold_data %>%
  filter(.metric == "j_index") %>%
  filter(.estimate == max(.estimate)) %>%
  pull(.threshold)

ggplot(threshold_data, aes(x = .threshold, y = .estimate, color = .metric, alpha = group)) +
  geom_line() +
  theme_minimal() +
  scale_color_viridis_d(end = 0.9) +
  scale_alpha_manual(values = c(.4, 1), guide = "none") +
  geom_vline(xintercept = max_j_index_threshold, alpha = .6, color = "grey30") +
  labs(
    x = "'CNV' Threshold\n(above this value is considered 'likely CNV carrier')",
    y = "Metric Estimate",
    title = "Balancing performance by varying the threshold",
    subtitle = "Sensitivity or specificity alone might not be enough!\nVertical line = Max J-Index"
  )

#This is probably a better plot to have when we have a real cross validated model
#tested on hold out datasets, but its an interesting plot to start with



#Select variables with importance > 0
lasso_vars = 
  lasso_vi |> 
  filter(Importance > 0) |>
  mutate(type = "LASSO")

### Save selected variables ====

# we could save this
write_csv(lasso_vars, "selected_variables_lasso.csv")



# Statistical Significance =======

# I suppose one approach might be to fit a glm for every variable -
#we would need to do this to each variable separately 
DF2 = 
  DF |>
  mutate(group = map_chr(GenotypeCode,~ifelse(.x == 16,"Control","ND-CNV"))) |>
  mutate(group = factor(group, levels = c("Control","ND-CNV"))) |>
  select(-c(GenotypeCode)) |>
  relocate(IDs,group) |>
  pivot_longer(-c(IDs,group),names_to = "variable",values_to = "value") |>
  group_by(variable) |>
  nest()

DF2 = 
  DF2 |>
  mutate(pv = map_dbl(data,~glm(group ~ value,
                                data = .x, 
                                family = binomial(link = "logit")) %>%
                        broom::tidy() |>
                        filter(term == "value") |>
                        pull(p.value)))

#Lets select vraiables that are significant, p corrected by bonferroning
glm_vars = 
  DF2 |>
  filter(pv < (0.05/nrow(DF2))) |>
  select(variable,pv)

#This still produces loads of variables
write_csv(glm_vars, "selected_variables_glm.csv")

# Single variable AUC =====

#Perhaps a better approach would be to look at each variable in turn and determine the AUC

#I suppose we could do a quick tidymodels here

tidymodels_prefer()

#Make a tidymodel for a glm
log_model = 
  logistic_reg(mode = "classification",
               engine = "glm") 

log_wflow = 
  workflow() |>
  add_model(log_model) |>
  add_formula(group ~ value)

#Make 10 fold cv for each variable
DF3 = 
  DF2 |>
  mutate(folds = map(data,vfold_cv,v = 10))

#Make a control variable for model fitting
keep_pred <- 
  control_resamples(save_pred = TRUE, save_workflow = TRUE)

#Fit 10 fold CV for classification accuracy to all variables
DF3 = 
  DF3 |>
  mutate(var_metrics = map(folds,~log_wflow |> 
                             fit_resamples(resamples = .x,
                                           control = keep_pred) |> 
                             collect_metrics() |> 
                             dplyr::select(.metric,mean,std_err) |> 
                             filter(.metric == "roc_auc")
  ))

#Extract the average AUC for each variable and plot
DF3 |> 
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

# So most variables have an AUC which suggests they are at least somewhat informative

# We could take the top say 30 variables
LR_vars = 
  DF3 |> 
  dplyr::select(variable,var_metrics) |>
  unnest(var_metrics) |>
  ungroup() |>
  slice_max(order_by = mean,n = 30) |>
  mutate(type = "single_LR_auc")

write_csv(LR_vars, "selected_variables_LR_auc.csv")

#Alternatively you could leave out each variable in turn

#There is actually a neater way to do this with workflowsets


#Variance Inflation Factor ========

#Another approach, which I thought of, but which essentially doesn't work,
#is to fit a glm with all predictors at once, then try and calculate
#the variance inflation factor on this model, then remove variables with
#the highest VIF. It doesn't really work because without regularisation
#the model basically is just chaos. Also you need to do imputation to get 
#rid of enough missing data to make the model work (or maybe you could 
#use a GLMM?)

# We are going to do imputation in order to get a complete case dataset for these initial processes

#We will use the mix GB package and create 5 imputed datasets
# library(mixgb)
# 
# # Apply this to our own data
# MIXGB <- 
#   Mixgb$new(DF |>
#               select(-c(IDs,GenotypeCode)) |>
#               as.data.frame(),
#             pmm.type = "auto",
#             pmm.k = 5)
# 
# mixgb.data <- 
#   MIXGB$impute(m=5)
# 
# #Fit a model to each imputed dataset
# m_1 =
#   tibble(iter = 1:5, data = mixgb.data) |>
#   mutate(data = map(data,as_tibble)) |>
#   mutate(data = map(data,~bind_cols(DF0 |> select(group),.x))) |>
#   mutate(model = map(data,~glm(group ~ .,family = binomial(link = "logit"), data = .x)))
# 
# m_1 = 
#   m_1 |>
#   mutate(vif_data = map(model,~vif(.x) |> 
#                           as_tibble(rownames = "variable") |>
#                           arrange(value))) |>
#   mutate(plots = map(vif_data,~  ggplot(.x,aes(value)) +
#                        geom_histogram(binwidth = 50)))
# 
# patchwork::wrap_plots(m_1$plots,ncol = 5)



#This does not work tremendously well, and we can see that we get really huge variance inflation

#Lets see if we can do the correlation



# Variable Correlation =====

# Another way to (negatively) select variables might be to calculate the correlation between all pairs of variables
#and then remove variables which have the highest correlation with all other variables, on the basis that they
#probably are not going to be very informative in addition to those less correlated variables


DF1 = read_csv("ImputedData2.csv")
DF0 = DF1


## Set up dataframe for polychoric correlation calculation
DF1 = 
  DF1 |>
  dplyr::select(-group)


## Set up a matrix of combinations to iterate through
Is   <- expand_grid(x = 1:ncol(DF1),y = 1:ncol(DF1))
cors <- matrix(nrow = ncol(DF1), ncol = ncol(DF1))


## Set up parallel processing cluster

numCores <- 16
cl       <- makeSOCKcluster(numCores)

registerDoSNOW(cl)

clusterExport(cl, "Is")
clusterExport(cl, "DF1")

clusterEvalQ(cl, library(polycor))
clusterEvalQ(cl, library(dplyr))

## set up progress bar
pb       <- txtProgressBar(max = nrow(Is), style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts     <- list(progress = progress)

## calculate polychoric correlation between each variable
r <-
  
  foreach(i = 1:nrow(Is), .options.snow = opts, .combine = "c") %dopar% {
    
    vs <- (Is[i,])
    v1 <- pull(vs[1])
    v2 <- pull(vs[2])
    
    cors[v1,v2] <- polychor(DF1[,v1] |> pull(),
                            DF1[,v2] |> pull(), 
                            ML = FALSE)
    
  }

stopCluster(cl)

#Recover our correlation matrix from the parallel thing
r2 <- 
  matrix(r, 
         nrow     = sqrt(length(r)), 
         ncol     = sqrt(length(r)), 
         dimnames = list(names(DF1), names(DF1)))


#Not sure what this thing about positive definite matrixes  - does it turn the correlations in covariances?
r2 <- Matrix::nearPD(r2, corr = T)
r2 <- as.matrix(r2$mat)

#Plot this (very hard to interpret)
ggcorrplot::ggcorrplot(r2,hc.order = TRUE,tl.cex = 4)

## Identify high correlation variables =====

#We can identify high correlation variables
high_cor = caret::findCorrelation(r2, 
                                  cutoff = 0.5, 
                                  verbose = TRUE, names = TRUE, exact = FALSE)

#If we removed all variables with a > 0.5 absolute correlation, we would be left with = 
DF_lowCor = 
  tibble(variable = DF1 |>
           dplyr::select(-all_of(high_cor)) |>
           names(),
         type = "low_corr_0.5")

#Now we can save these
write_csv(DF_lowCor, "selected_variables_low_corr.csv")




## Do PCA on the correlation matrix =====

#While we're here we can do PCA on the correlation/covariance matrix

# pca <- princomp(x = DF1,fix_sign = T, cor = FALSE, scores = TRUE)
pca <- princomp(covmat = r2, fix_sign = T, cor = FALSE, scores = TRUE)

#Plot the PCA
plot(pca)


### Extract PC Scores =====

trans  <- t(pca$loadings)
inve   <- solve(as.matrix(trans))

#This isn't working for some horrible reason
scores <- as.matrix(DF1) %*% as.matrix(inve)

names(scores) <- names(pca$loadings)

Scores2 = 
  bind_cols(DF0 |> 
              dplyr::select(group) ,
            scores |>
              as_tibble()) |>
  relocate(group)


#And then what are we doing with this? We might make some plots I suppose

Scores2 |>
  dplyr::select(group,Comp.1:Comp.4) |>
  ggplot(aes(x = .panel_x, y = .panel_y, color = group, fill = group)) +
  geom_point(alpha = 0.8, size = 1) +
  geom_autodensity(alpha = .3) +
  facet_matrix(vars(-group), layer.diag = 2) + 
  scale_color_brewer(palette = "Dark2") + 
  scale_fill_brewer(palette  = "Dark2") + 
  theme_bw()




# This is hwo to do PLS with raw code
plsda_1 = mixOmics::plsda(Y = DF0 |> pull(group), X = DF0 |> dplyr::select(-group), ncomp = 4)
plotIndiv(plsda_1, ind.names = DF0 |> pull(group), ellipse = TRUE, legend =TRUE)



